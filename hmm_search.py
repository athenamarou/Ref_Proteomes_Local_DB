
"""
HMM Search Pipeline for CGLab UniProt Reference Proteome Database
===================================================================
Performs HMMER hmmscan searches using Pfam HMM profiles against all
(or a filtered subset of) proteins stored in the local UniProt database.

This pipeline handles the entire workflow:
    1. Export protein sequences from MySQL → FASTA
    2. Download and prepare Pfam HMM profiles
    3. Split FASTA into chunks for parallel processing
    4. Run hmmscan on each chunk (locally)
    5. Parse hmmscan output (domain table format)
    6. Load results into a MySQL results table

Supported modes:
    - Local (multi-process on a single machine-default)

Usage examples:
    # Full scan (all proteins, all HMMs) — local mode with 16 threads
    python hmm_search_pipeline.py \\
        --version 2026_01 \\
        --threads 16 \\
        --output-dir ./hmm_results

    # Scan only mammalian proteins (TaxID 9606, 10090)
    python hmm_search_pipeline.py \\
        --version 2026_01 \\
        --taxon-ids 9606 10090 \\
        --threads 16 \\
        --output-dir ./hmm_results


    # Only export FASTA (no search)
    python hmm_search_pipeline.py \\
        --version 2026_01 \\
        --export-only \\
        --output-dir ./hmm_results

    # Only import results (after manual hmmscan run)
    python hmm_search_pipeline.py \\
        --version 2026_01 \\
        --import-only \\
        --results-dir ./hmm_results

Requirements:
    - HMMER 3.4+ (hmmscan, hmmpress)  — http://hmmer.org/
    - Python 3.9+ with mysql-connector-python, biopython
    - MySQL access (same .env as the sync scripts)

"""

import argparse
import csv
import io
import os
import re
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import mysql.connector
from dotenv import load_dotenv

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Pfam FTP URLs (current release as of 2026-01)
PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_VERSION_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfam_version.txt"

# Default paths
DEFAULT_CHUNK_SIZE = 50000        # proteins per FASTA chunk
DEFAULT_HMMSCAN_THREADS = 4       # threads per hmmscan process
DEFAULT_LOCAL_JOBS = 4            # parallel hmmscan processes (local mode)
DEFAULT_EVALUE = 1e-5             # default significance threshold
DEFAULT_DOM_EVALUE = 1e-3         # default domain-level e-value cutoff

# Load .env from the script's directory (same as the sync scripts)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(_SCRIPT_DIR, ".env"))


# ---------------------------------------------------------------------------
# 1. DATABASE EXPORT — pull sequences from MySQL
# ---------------------------------------------------------------------------

class SequenceExporter:
    """
    Exports protein sequences from the CGLab UniProt MySQL database to FASTA.

    Memory usage is minimal: uses buffered cursors to stream rows one at
    a time from MySQL. Only the current chunk (~50K proteins, ~200 MB
    of text) is ever held in memory.

    IMPORTANT — Disk usage:
        The full export writes ~50–60 GB of FASTA files to disk.
        Use --stream mode to avoid this: it processes each chunk
        through hmmscan then immediately deletes the temp FASTA.

    This replicates the retrieval logic from get_reference_uniprot_set_lib.py
    but uses a streaming, chunk-based approach to avoid loading 133M records
    into memory.

    FASTA header format:
        >{taxon_id}.{accession}|{proteome_id}|{name}|{organism}
        MEEPQSDPSV...

    This format preserves all metadata needed for downstream analysis
    while keeping the header parseable.
    """

    def __init__(self, version, output_dir, taxon_ids=None, proteome_ids=None):
        self.version = version
        self.output_dir = Path(output_dir)
        self.taxon_ids = taxon_ids
        self.proteome_ids = proteome_ids

        self.fasta_dir = self.output_dir / "fasta_chunks"
        self.fasta_dir.mkdir(parents=True, exist_ok=True)

        self.config = {
            "host":     os.getenv("DB_HOST", "localhost"),
            "user":     os.getenv("DB_USER", "cglab_user"),
            "password": os.getenv("DB_PASSWORD", ""),
            "database": os.getenv("DB_NAME", "uniprot_db_cglab"),
        }

    def _build_query(self):
        """
        Build the SQL query with streaming (SSCursor) for large exports.

        Uses MySQL SSCursor (unbuffered) to avoid loading the entire
        result set into memory — critical for 133M proteins.
        """
        query = """
            SELECT p.accession, p.name, p.organism,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM   proteins  p
            JOIN   sequences s ON p.seq_id = s.seq_id
            WHERE  p.version = %s
        """
        params = [self.version]

        if self.taxon_ids:
            placeholders = ", ".join(["%s"] * len(self.taxon_ids))
            query += f" AND p.taxon_id IN ({placeholders})"
            params.extend(self.taxon_ids)

        if self.proteome_ids:
            placeholders = ", ".join(["%s"] * len(self.proteome_ids))
            query += f" AND p.proteome_id IN ({placeholders})"
            params.extend(self.proteome_ids)

        return query, params

    def _write_fasta_record(self, file_handle, row):
        """Write a single protein as a FASTA record."""
        accession = row["accession"]
        organism = row["organism"].replace("|", "_").replace("\x00", "")
        name = row["name"].replace("|", "_")
        header = (
            f">{row['taxon_id']}.{accession}"
            f"|{row['proteome_id']}"
            f"|{name}"
            f"|{organism}"
        )
        file_handle.write(f"{header}\n{row['sequence']}\n")

    def export(self, chunk_size=DEFAULT_CHUNK_SIZE):
        """
        Stream all matching proteins from MySQL and write FASTA chunks.

        Uses SSCursor for memory-efficient streaming of 100M+ rows.

        Returns:
            list[str]: paths to all chunk files
        """
        query, params = self._build_query()

        # Count total proteins first (regular cursor)
        count_query = query.replace(
            "SELECT p.accession, p.name, p.organism, p.taxon_id, p.proteome_id, s.sequence",
            "SELECT COUNT(*)"
        )
        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()
        cursor.execute(count_query, params)
        total = cursor.fetchone()[0]
        cursor.close()
        conn.close()

        print(f"\n{'='*60}")
        print(f"Exporting sequences from MySQL")
        print(f"  Version:  {self.version}")
        print(f"  Total:    {total:,} proteins")
        if self.taxon_ids:
            print(f"  Taxa:     {self.taxon_ids}")
        if self.proteome_ids:
            print(f"  Proteomes:{self.proteome_ids}")
        print(f"  Chunk:    {chunk_size:,} proteins/file")
        print(f"{'='*60}\n")

        if total == 0:
            print("✗ No proteins found for the given filters.")
            return []

        # Stream with SSCursor
        conn = mysql.connector.connect(**self.config, buffered=False)
        cursor = conn.cursor(dictionary=True)

        chunk_files = []
        chunk_num = 0
        count_in_chunk = 0
        total_exported = 0
        chunk_file = None
        start_time = time.time()

        try:
            cursor.execute(query, params)

            for row in cursor:
                if count_in_chunk == 0:
                    # Open new chunk file
                    chunk_num += 1
                    chunk_path = self.fasta_dir / f"chunk_{chunk_num:05d}.fasta"
                    chunk_file = open(chunk_path, "w")
                    chunk_files.append(str(chunk_path))

                self._write_fasta_record(chunk_file, row)
                count_in_chunk += 1
                total_exported += 1

                if count_in_chunk >= chunk_size:
                    chunk_file.close()
                    chunk_file = None
                    count_in_chunk = 0

                    elapsed = time.time() - start_time
                    rate = total_exported / elapsed if elapsed > 0 else 0
                    print(
                        f"  Chunk {chunk_num}: {total_exported:,}/{total:,} "
                        f"({total_exported/total*100:.1f}%) — {rate:,.0f} seq/s"
                    )

            # Close last chunk
            if chunk_file:
                chunk_file.close()

            # Write a manifest file listing all chunks
            manifest = self.fasta_dir / "manifest.txt"
            with open(manifest, "w") as f:
                f.write(f"# HMM Search Pipeline FASTA Manifest\n")
                f.write(f"# Version: {self.version}\n")
                f.write(f"# Total proteins: {total_exported:,}\n")
                f.write(f"# Chunks: {len(chunk_files)}\n")
                f.write(f"# Created: {datetime.now().isoformat()}\n")
                f.write(f"#\n")
                for fp in chunk_files:
                    f.write(f"{fp}\n")

        finally:
            cursor.close()
            conn.close()

        elapsed = time.time() - start_time
        print(f"\n✓ Exported {total_exported:,} proteins into {len(chunk_files)} chunks "
              f"in {elapsed:.0f}s ({total_exported/elapsed:,.0f} seq/s)")
        print(f"  Chunk directory: {self.fasta_dir}")
        print(f"  Manifest: {manifest}")

        return chunk_files


# ---------------------------------------------------------------------------
# 2. HMM PROFILE MANAGEMENT
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 1b. STREAMING EXPORT — MySQL → hmmscan with no persistent FASTA files
# ---------------------------------------------------------------------------

class SequenceStreamer:
    """
    Streams protein sequences directly from MySQL into temp FASTA files,
    processes each with hmmscan, then immediately deletes the FASTA.

    This is the memory-and-disk-efficient alternative to SequenceExporter.
    Only ONE chunk's worth of FASTA (~200 MB) exists on disk at any time.
    Memory usage stays at ~50–100 MB regardless of total database size.

    Uses MySQL LIMIT/OFFSET to page through the result set in chunks.
    Each chunk is:
        1. Written to a temporary FASTA file
        2. Passed to hmmscan
        3. Deleted after hmmscan completes

    This class does NOT use SSCursor because we need LIMIT/OFFSET
    pagination, which works fine with regular cursors for the small
    page sizes we use (50K rows).
    """

    def __init__(self, version, output_dir, taxon_ids=None, proteome_ids=None):
        self.version = version
        self.output_dir = Path(output_dir)
        self.taxon_ids = taxon_ids
        self.proteome_ids = proteome_ids
        self.tmp_dir = self.output_dir / "tmp_fasta"
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

        # Initialize our checkpoint variable
        self.last_accession = None

        self.config = {
            "host":     os.getenv("DB_HOST", "localhost"),
            "user":     os.getenv("DB_USER", "cglab_user"),
            "password": os.getenv("DB_PASSWORD", ""),
            "database": os.getenv("DB_NAME", "uniprot_db_cglab"),
        }

    def _build_query(self):
        """Build the base query (without LIMIT/OFFSET)."""
        query = """
            SELECT p.accession, p.name, p.organism,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM   proteins  p
            JOIN   sequences s ON p.seq_id = s.seq_id
            WHERE  p.version = %s
        """
        params = [self.version]

        if self.taxon_ids:
            placeholders = ", ".join(["%s"] * len(self.taxon_ids))
            query += f" AND p.taxon_id IN ({placeholders})"
            params.extend(self.taxon_ids)

        if self.proteome_ids:
            placeholders = ", ".join(["%s"] * len(self.proteome_ids))
            query += f" AND p.proteome_id IN ({placeholders})"
            params.extend(self.proteome_ids)

        return query, params

    def _write_fasta_record(self, file_handle, row):
        """Write a single protein as a FASTA record."""
        accession = row["accession"]
        organism = (row["organism"] or "").replace("|", "_").replace("\x00", "") #works even with null
        name = row["name"].replace("|", "_")
        header = (
            f">{row['taxon_id']}.{accession}"
            f"|{row['proteome_id']}"
            f"|{name}"
            f"|{organism}"
        )
        file_handle.write(f"{header}\n{row['sequence']}\n")

    def get_total_count(self):
        """Return total number of matching proteins."""
        # FIX: Construct the count query directly instead of using fragile string replacement
        query = "SELECT COUNT(*) FROM proteins p JOIN sequences s ON p.seq_id = s.seq_id WHERE p.version = %s"
        params = [self.version]

        if self.taxon_ids:
            placeholders = ", ".join(["%s"] * len(self.taxon_ids))
            query += f" AND p.taxon_id IN ({placeholders})"
            params.extend(self.taxon_ids)

        if self.proteome_ids:
            placeholders = ", ".join(["%s"] * len(self.proteome_ids))
            query += f" AND p.proteome_id IN ({placeholders})"
            params.extend(self.proteome_ids)

        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()
        cursor.execute(query, params)
        total = cursor.fetchone()[0]

        # Ensure we consume any remaining data so MariaDB doesn't throw an 'Unread result' error
        cursor.fetchall() 
        cursor.close()
        conn.close()
        return total

    def stream_chunk(self, chunk_num, chunk_size):
        """
        Fetch one chunk of proteins from MySQL using Keyset Pagination,
        write to a temp FASTA file, and return the file path.
        """
        query, params = self._build_query()

        # Apply the checkpoint if we have one (skip this for the very first chunk)
        if self.last_accession:
            query += " AND p.accession > %s"
            params.append(self.last_accession)

        # Only use LIMIT, not using OFFSET because its slower
        query += " ORDER BY p.accession LIMIT %s"
        page_params = params + [chunk_size]

        tmp_path = self.tmp_dir / f"stream_chunk_{chunk_num:05d}.fasta"

        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor(dictionary=True)
        cursor.execute(query, page_params)

        rows = cursor.fetchall()
        cursor.close()
        conn.close()

        if not rows:
            return None, 0

        # UPDATE THE CHECKPOINT:
        # Grabs the accession of the very last protein in this chunk
        self.last_accession = rows[-1]["accession"]

        with open(tmp_path, "w") as f:
            for row in rows:
                self._write_fasta_record(f, row)

        return str(tmp_path), len(rows)


# ---------------------------------------------------------------------------
# 2. HMM PROFILE MANAGEMENT
# ---------------------------------------------------------------------------


class HMMProfileManager:
    """
    Downloads, verifies, and prepares Pfam HMM profiles for hmmscan.

    The Pfam-A.hmm.gz file is the full library of curated HMM models
    (~20,000 families). After download, `hmmpress` indexes it for
    hmmscan searches.
    """

    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.hmm_dir = self.output_dir / "hmm_profiles"
        self.hmm_dir.mkdir(parents=True, exist_ok=True)
        self.hmm_file = self.hmm_dir / "Pfam-A.hmm"
        self.version_file = self.hmm_dir / "pfam_version.txt"

    def is_prepared(self):
        """Check if the HMM database has already been pressed."""
        if not self.hmm_file.exists():
            return False
        # hmmpress creates .h3f, .h3i, .h3m, .h3p files
        required_suffixes = [".h3f", ".h3i", ".h3m", ".h3p"]
        return all(
            (self.hmm_dir / f"Pfam-A.hmm{s}").exists()
            for s in required_suffixes
        )

    def download(self):
        """Download and decompress the Pfam-A HMM database."""
        if self.hmm_file.exists():
            print(f"  Pfam-A.hmm already exists ({self.hmm_file.stat().st_size / 1e9:.2f} GB)")
            return True

        compressed = self.hmm_dir / "Pfam-A.hmm.gz"
        print(f"  Downloading Pfam-A.hmm.gz from {PFAM_HMM_URL}...")

        try:
            import requests
            response = requests.get(PFAM_HMM_URL, stream=True)
            response.raise_for_status()
            total = int(response.headers.get("Content-Length", 0))

            downloaded = 0
            with open(compressed, "wb") as f:
                for chunk in response.iter_content(chunk_size=8 * 1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total > 0:
                            print(
                                f"\r    {downloaded/1e9:.2f}/{total/1e9:.2f} GB "
                                f"({downloaded/total*100:.0f}%)", end=""
                            )

            print(f"\n  Download complete. Decompressing...")
            import gzip
            with gzip.open(compressed, "rb") as gz_in:
                with open(self.hmm_file, "wb") as out:
                    while True:
                        block = gz_in.read(8 * 1024 * 1024)
                        if not block:
                            break
                        out.write(block)

            # Clean up compressed file
            compressed.unlink()
            print(f"  Pfam-A.hmm ready ({self.hmm_file.stat().st_size / 1e9:.2f} GB)")
            return True

        except Exception as e:
            print(f"\n  ✗ Download failed: {e}")
            return False

    def press(self):
        """Run hmmpress to create binary index files."""
        if self.is_prepared():
            print("  HMM database already indexed (hmmpress skipped).")
            return True

        print(f"  Running hmmpress on {self.hmm_file}...")
        try:
            result = subprocess.run(
                ["hmmpress", "-f", str(self.hmm_file)],
                capture_output=True, text=True, timeout=600
            )
            if result.returncode != 0:
                print(f"  ✗ hmmpress failed:\n{result.stderr}")
                return False
            print("  hmmpress completed successfully.")
            return True
        except FileNotFoundError:
            print("  ✗ hmmpress not found. Install HMMER: http://hmmer.org/")
            return False
        except subprocess.TimeoutExpired:
            print("  ✗ hmmpress timed out (600s).")
            return False

    def prepare(self):
        """Download + press in one call. Returns True on success."""
        if self.is_prepared():
            print("  HMM database ready (skipped download/press).")
            return True
        if not self.download():
            return False
        return self.press()

    def get_version(self):
        """Return Pfam version string if available."""
        if self.version_file.exists():
            return self.version_file.read_text().strip()
        return "unknown"


# ---------------------------------------------------------------------------
# 3. HMMSCAN EXECUTION
# ---------------------------------------------------------------------------

def run_single_hmmscan(args_dict):
    """
    Run hmmscan on a single FASTA chunk.

    This function is designed to be called by ProcessPoolExecutor (must
    be top-level for pickling). It runs hmmscan with the domain table
    output format (--domtblout) which is easy to parse.

    Parameters (passed as dict for pickling):
        fasta_file: str
        hmm_file: str
        output_file: str
        evalue: float
        dom_evalue: float
        threads: int
        cpus: int  (DEPRECATED: use threads)

    Returns:
        dict with stats
    """
    fasta_file = args_dict["fasta_file"]
    hmm_file = args_dict["hmm_file"]
    output_file = args_dict["output_file"]
    evalue = args_dict["evalue"]
    dom_evalue = args_dict.get("dom_evalue", DEFAULT_DOM_EVALUE)
    threads = args_dict.get("threads", DEFAULT_HMMSCAN_THREADS)

    # hmmscan command
    cmd = [
        "hmmscan",
        "--cpu", str(threads),
        "--domtblout", output_file,
        "--noali",                    # don't print alignments (faster)
        "-E", str(evalue),            # sequence-level e-value cutoff
        "--domE", str(dom_evalue),    # domain-level e-value cutoff
        hmm_file,
        fasta_file,
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,  # 2 hours per chunk
        )

        if result.returncode != 0:
            return {
                "fasta_file": fasta_file,
                "success": False,
                "error": result.stderr[:1000],
            }

        # Clean up temp FASTA after successful scan (streaming mode)
        cleanup = args_dict.get("cleanup", False)
        if cleanup and os.path.exists(fasta_file):
            os.remove(fasta_file)

        return {
            "fasta_file": fasta_file,
            "output_file": output_file,
            "success": True,
        }

    except subprocess.TimeoutExpired:
        return {
            "fasta_file": fasta_file,
            "success": False,
            "error": "Timeout (2h)",
        }
    except FileNotFoundError:
        return {
            "fasta_file": fasta_file,
            "success": False,
            "error": "hmmscan not found — install HMMER",
        }


def run_local_search(fasta_dir, hmm_file, output_dir,
                     threads=DEFAULT_HMMSCAN_THREADS,
                     parallel_jobs=DEFAULT_LOCAL_JOBS,
                     evalue=DEFAULT_EVALUE,
                     dom_evalue=DEFAULT_DOM_EVALUE,
                     cleanup=False):
    """
    Run hmmscan on all FASTA chunks locally using multiprocessing.

    Parameters:
        fasta_dir: str  — directory containing chunk_*.fasta files
        hmm_file: str   — path to pressed Pfam-A.hmm
        output_dir: str — directory for hmmscan results
        threads: int    — CPU threads per hmmscan call
        parallel_jobs: int — max concurrent hmmscan processes
        evalue: float   — sequence-level e-value cutoff
        dom_evalue: float — domain-level e-value cutoff
        cleanup: bool   — delete FASTA files after hmmscan completes

    Returns:
        dict: summary stats
    """
    fasta_dir = Path(fasta_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(fasta_dir.glob("chunk_*.fasta")) + sorted(fasta_dir.glob("stream_chunk_*.fasta"))
    if not chunks:
        print("✗ No chunk files found.")
        return {"success": False, "error": "No chunks"}

    print(f"\n{'='*60}")
    print(f"Running hmmscan — LOCAL MODE")
    print(f"  Chunks:         {len(chunks)}")
    print(f"  Threads/call:   {threads}")
    print(f"  Parallel jobs:  {parallel_jobs}")
    print(f"  E-value cutoff: {evalue}")
    print(f"  Dom E-value:    {dom_evalue}")
    print(f"  Cleanup FASTA:  {'Yes' if cleanup else 'No'}")
    print(f"{'='*60}\n")

    tasks = []
    for chunk in chunks:
        output_file = output_dir / f"{chunk.stem}_hmmscan.domtblout.txt"
        tasks.append({
            "fasta_file": str(chunk),
            "hmm_file": hmm_file,
            "output_file": str(output_file),
            "evalue": evalue,
            "dom_evalue": dom_evalue,
            "threads": threads,
            "cleanup": cleanup,
        })

    completed = 0
    failed = 0
    start_time = time.time()

    # Adjust parallel_jobs based on CPU count
    total_cpus_needed = parallel_jobs * threads
    available_cpus = os.cpu_count() or 4
    if total_cpus_needed > available_cpus:
        print(f"  ⚠ Requested {total_cpus_needed} CPUs but only {available_cpus} available.")
        print(f"  Adjusting parallel_jobs to {max(1, available_cpus // threads)}")
        parallel_jobs = max(1, available_cpus // threads)

    with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
        futures = {executor.submit(run_single_hmmscan, task): task for task in tasks}

        for future in as_completed(futures):
            task = futures[future]
            result = future.result()
            completed += 1

            if result["success"]:
                status = "✓"
            else:
                status = f"✗ ({result.get('error', 'unknown')[:50]})"
                failed += 1

            elapsed = time.time() - start_time
            rate = completed / elapsed if elapsed > 0 else 0
            eta = (len(tasks) - completed) / rate if rate > 0 else 0

            print(
                f"  [{completed}/{len(tasks)}] {status} "
                f"{os.path.basename(task['fasta_file'])} "
                f"— {rate:.2f} chunks/min, ETA: {eta/60:.0f}min"
            )

    elapsed = time.time() - start_time
    print(f"\n✓ hmmscan complete: {completed - failed}/{len(tasks)} succeeded "
          f"({failed} failed) in {elapsed/60:.1f}min")

    return {
        "success": failed == 0,
        "total": len(tasks),
        "succeeded": completed - failed,
        "failed": failed,
        "duration_s": elapsed,
    }


# ---------------------------------------------------------------------------
# 4. RESULTS PARSING & IMPORT
# ---------------------------------------------------------------------------

class HMMResultsImporter:
    """
    Parses hmmscan --domtblout files and loads results into MySQL.

    The domain table format has columns (space-delimited, fixed widths):
        target_name  target_acc  tlen  query_name  query_acc  qlen
        evalue  score  bias  #  of  HMM  from  to  hmm_cov
        #  of  dom  from  to  dom_cov  env_from  env_to  env_cov
        acc

    Results are stored in a new table: hmm_search_results
    """

    def __init__(self, version, output_dir):
        self.version = version
        self.output_dir = Path(output_dir)
        self.config = {
            "host":     os.getenv("DB_HOST", "localhost"),
            "user":     os.getenv("DB_USER", "cglab_user"),
            "password": os.getenv("DB_PASSWORD", ""),
            "database": os.getenv("DB_NAME", "uniprot_db_cglab"),
        }

    def create_results_table(self):
        """Create the hmm_search_results table."""
        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()
        try:
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS hmm_search_results (
                    id BIGINT AUTO_INCREMENT PRIMARY KEY,
                    version VARCHAR(10) NOT NULL,
                    accession VARCHAR(20) NOT NULL,
                    taxon_id INT,
                    proteome_id VARCHAR(20),
                    protein_name VARCHAR(100),

                    -- HMM profile info
                    hmm_name VARCHAR(50) NOT NULL,
                    hmm_accession VARCHAR(20),
                    hmm_type ENUM('Pfam', 'TIGRFAM', 'SUPERFAMILY', 'other') DEFAULT 'Pfam',

                    -- Alignment stats
                    full_evalue DOUBLE,
                    full_score DOUBLE,
                    full_bias DOUBLE,

                    -- Domain hit
                    domain_number INT,
                    domain_count INT,
                    domain_evalue DOUBLE,
                    domain_score DOUBLE,
                    domain_bias DOUBLE,

                    -- Coordinates (1-based, inclusive)
                    hmm_from INT,
                    hmm_to INT,
                    ali_from INT,
                    ali_to INT,
                    env_from INT,
                    env_to INT,

                    -- Coverage
                    hmm_coverage FLOAT COMMENT 'fraction of HMM matched',
                    protein_coverage FLOAT COMMENT 'fraction of protein matched',

                    -- Posterior probability
                    posterior_prob FLOAT,

                    -- Metadata
                    search_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

                    INDEX idx_version (version),
                    INDEX idx_accession (accession),
                    INDEX idx_hmm (hmm_name),
                    INDEX idx_hmm_acc (hmm_accession),
                    INDEX idx_taxon (taxon_id),
                    INDEX idx_proteome (proteome_id),
                    INDEX idx_evalue (full_evalue),
                    INDEX idx_domain_eval (domain_evalue),
                    INDEX idx_version_acc (version, accession)
                )
            """)
            conn.commit()
            print("  Results table 'hmm_search_results' ready.")
        finally:
            cursor.close()
            conn.close()

    def parse_domtblout(self, filepath):
        """
        Parse a single hmmscan --domtblout file.

        Yields:
            dict: one dict per domain hit
        """
        header_pattern = re.compile(r'^#')
        with open(filepath, 'r') as f:
            for line in f:
                if header_pattern.match(line):
                    continue
                if not line.strip():
                    continue

                fields = line.split()
                if len(fields) < 22:
                    continue

                # Parse query name to extract metadata
                # Format: {taxon_id}.{accession}|{proteome_id}|{name}|{organism}
                query_name = fields[3]
                parts = query_name.split("|")
                accession = parts[0].split(".")[-1] if "." in parts[0] else parts[0]
                taxon_id = None
                proteome_id = None
                protein_name = None

                try:
                    taxon_id = int(parts[0].split(".")[0])
                except (ValueError, IndexError):
                    pass

                if len(parts) >= 2:
                    proteome_id = parts[1]
                if len(parts) >= 3:
                    protein_name = parts[2]

                # Determine HMM type from accession format
                hmm_acc = fields[1]
                if hmm_acc.upper().startswith("PF"):
                    hmm_type = "Pfam"
                elif hmm_acc.upper().startswith("TIGR"):
                    hmm_type = "TIGRFAM"
                elif hmm_acc.upper().startswith("SSF"):
                    hmm_type = "SUPERFAMILY"
                else:
                    hmm_type = "other"

                tlen = float(fields[2])
                qlen = float(fields[5])

                yield {
                    "accession": accession,
                    "taxon_id": taxon_id,
                    "proteome_id": proteome_id,
                    "protein_name": protein_name,
                    "hmm_name": fields[0],
                    "hmm_accession": hmm_acc,
                    "hmm_type": hmm_type,
                    "full_evalue": float(fields[6]),
                    "full_score": float(fields[7]),
                    "full_bias": float(fields[8]),
                    "domain_number": int(fields[9]),
                    "domain_count": int(fields[10]),
                    "domain_evalue": float(fields[12]),
                    "domain_score": float(fields[13]),
                    "domain_bias": float(fields[14]),
                    "hmm_from": int(fields[15]),
                    "hmm_to": int(fields[16]),
                    "ali_from": int(fields[17]),
                    "ali_to": int(fields[18]),
                    "env_from": int(fields[19]),
                    "env_to": int(fields[20]),
                    "hmm_coverage": (int(fields[16]) - int(fields[15]) + 1) / tlen if tlen > 0 else 0,
                    "protein_coverage": (int(fields[18]) - int(fields[17]) + 1) / qlen if qlen > 0 else 0,
                    "posterior_prob": float(fields[21]) if len(fields) > 21 else 0.0,
                }
                

    def import_results(self, batch_size=100000):
        """
        Import all hmmscan results from the output directory into MySQL.

        Uses bulk INSERT with executemany for performance.
        """
        results_dir = self.output_dir
        domtblout_files = sorted(results_dir.glob("*.domtblout.txt"))

        if not domtblout_files:
            print("✗ No hmmscan result files found.")
            return

        print(f"\n{'='*60}")
        print(f"Importing hmmscan results into MySQL")
        print(f"  Version:  {self.version}")
        print(f"  Files:    {len(domtblout_files)}")
        print(f"{'='*60}\n")

        self.create_results_table()

        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()

        # Check for existing data
        cursor.execute(
            "SELECT COUNT(*) FROM hmm_search_results WHERE version = %s",
            (self.version,)
        )
        existing = cursor.fetchone()[0]
        if existing > 0:
            print(f"  ⚠ {existing:,} existing results for version {self.version}.")
            print(f"  Appending (use --force to clear first).")

        total_imported = 0
        total_domains = 0
        batch = []
        start_time = time.time()

        try:
            # Enable bulk mode
            cursor.execute("SET SESSION foreign_key_checks = 0")
            cursor.execute("SET SESSION unique_checks = 0")

            insert_sql = """
                INSERT INTO hmm_search_results (
                    version, accession, taxon_id, proteome_id, protein_name,
                    hmm_name, hmm_accession, hmm_type,
                    full_evalue, full_score, full_bias,
                    domain_number, domain_count, domain_evalue, domain_score, domain_bias,
                    hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,
                    hmm_coverage, protein_coverage, posterior_prob
                ) VALUES (
                    %s, %s, %s, %s, %s,
                    %s, %s, %s,
                    %s, %s, %s,
                    %s, %s, %s, %s, %s,
                    %s, %s, %s, %s, %s, %s,
                    %s, %s, %s
                )
            """

            for filepath in domtblout_files:
                file_domains = 0
                for hit in self.parse_domtblout(filepath):
                    batch.append((
                        self.version,
                        hit["accession"], hit["taxon_id"], hit["proteome_id"],
                        hit["protein_name"],
                        hit["hmm_name"], hit["hmm_accession"], hit["hmm_type"],
                        hit["full_evalue"], hit["full_score"], hit["full_bias"],
                        hit["domain_number"], hit["domain_count"],
                        hit["domain_evalue"], hit["domain_score"], hit["domain_bias"],
                        hit["hmm_from"], hit["hmm_to"], hit["ali_from"], hit["ali_to"],
                        hit["env_from"], hit["env_to"],
                        hit["hmm_coverage"], hit["protein_coverage"],
                        hit["posterior_prob"],
                    ))
                    file_domains += 1

                    if len(batch) >= batch_size:
                        cursor.executemany(insert_sql, batch)
                        conn.commit()
                        total_imported += len(batch)
                        batch = []

                total_domains += file_domains
                elapsed = time.time() - start_time
                rate = total_imported / elapsed if elapsed > 0 else 0
                print(
                    f"  ✓ {os.path.basename(filepath)}: "
                    f"{file_domains:,} domains | "
                    f"Total: {total_imported:,} | {rate:,.0f} rows/s"
                )

            # Flush remaining
            if batch:
                cursor.executemany(insert_sql, batch)
                conn.commit()
                total_imported += len(batch)

            cursor.execute("SET SESSION foreign_key_checks = 1")
            cursor.execute("SET SESSION unique_checks = 1")

            elapsed = time.time() - start_time
            print(f"\n✓ Imported {total_imported:,} domain hits "
                  f"from {len(domtblout_files)} files in {elapsed/60:.1f}min")
            print(f"  ({total_imported/elapsed:,.0f} rows/s)")

            # Summary stats
            cursor.execute("""
                SELECT
                    COUNT(*) AS total_hits,
                    COUNT(DISTINCT accession) AS unique_proteins,
                    COUNT(DISTINCT hmm_accession) AS unique_hmms,
                    COUNT(DISTINCT taxon_id) AS unique_taxa,
                    AVG(full_evalue) AS avg_evalue,
                    MIN(full_evalue) AS best_evalue
                FROM hmm_search_results
                WHERE version = %s
            """, (self.version,))
            stats = cursor.fetchone()
            if stats:
                print(f"\n  Summary:")
                print(f"    Total domain hits:  {stats[0]:,}")
                print(f"    Unique proteins:    {stats[1]:,}")
                print(f"    Unique HMMs:        {stats[2]:,}")
                print(f"    Unique taxa:        {stats[3]:,}")
                print(f"    Avg E-value:        {stats[4]:.2e}")
                print(f"    Best E-value:       {stats[5]:.2e}")

        finally:
            cursor.close()
            conn.close()


# ---------------------------------------------------------------------------
# 5. MAIN CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="HMM Search Pipeline for CGLab UniProt Reference Proteome DB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full scan — local mode with 16 threads
  python hmm_search_pipeline.py --version 2026_01 --threads 16 --output-dir ./hmm_results

  # Scan only mammals
  python hmm_search_pipeline.py --version 2026_01 --taxon-ids 9606 10090 --threads 16 --output-dir ./hmm_results

  # Generate SLURM scripts
  python hmm_search_pipeline.py --version 2026_01 --mode slurm --slurm-ncpus 8 --output-dir ./hmm_results

  # Only export FASTA (no search)
  python hmm_search_pipeline.py --version 2026_01 --export-only --output-dir ./hmm_results

  # Only import results
  python hmm_search_pipeline.py --version 2026_01 --import-only --output-dir ./hmm_results
        """
    )

    # Required
    parser.add_argument(
        "--version", required=True,
        help="UniProt DB version (e.g., 2026_01)"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Base directory for all outputs (FASTA, results, etc.)"
    )

    # Filters
    parser.add_argument(
        "--taxon-ids", nargs="+", type=int, default=None,
        help="Restrict to one or more NCBI Taxonomy IDs (e.g., 9606 10090)"
    )
    parser.add_argument(
        "--proteome-ids", nargs="+", default=None,
        help="Restrict to specific proteome IDs (e.g., UP000005640)"
    )

    # Mode
    parser.add_argument(
        "--mode", choices=["local"], default="local",
        help="Execution mode: local (multiprocess)"
    )

    # Execution parameters
    parser.add_argument(
        "--threads", type=int, default=DEFAULT_HMMSCAN_THREADS,
        help=f"CPU threads per hmmscan call (default: {DEFAULT_HMMSCAN_THREADS})"
    )
    parser.add_argument(
        "--parallel-jobs", type=int, default=DEFAULT_LOCAL_JOBS,
        help=f"Max parallel hmmscan processes in local mode (default: {DEFAULT_LOCAL_JOBS})"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
        help=f"Proteins per FASTA chunk (default: {DEFAULT_CHUNK_SIZE})"
    )
    parser.add_argument(
        "--evalue", type=float, default=DEFAULT_EVALUE,
        help=f"Sequence-level E-value cutoff (default: {DEFAULT_EVALUE})"
    )
    parser.add_argument(
        "--dom-evalue", type=float, default=DEFAULT_DOM_EVALUE,
        help=f"Domain-level E-value cutoff (default: {DEFAULT_DOM_EVALUE})"
    )


    # Streaming mode
    parser.add_argument(
        "--stream", action="store_true", default=True,
        help="(Default) Stream chunks from MySQL one at a time — minimal disk usage (~1-2 GB). Temp FASTA files are deleted after hmmscan finishes."
    )
    parser.add_argument(
        "--no-stream", action="store_true",
        help="Disable streaming: export ALL proteins to FASTA files first, then search. Uses ~50-60 GB disk."
    )
    parser.add_argument(
        "--keep-fasta", action="store_true",
        help="Keep FASTA files on disk after hmmscan completes (default: delete them to save space)"
    )

    args = parser.parse_args()

    # Resolve streaming vs batch mode
    use_stream = args.stream and not args.no_stream
    cleanup_fasta = use_stream and not args.keep_fasta

    print(f"\n{'#'*60}")
    print(f"  HMM Search Pipeline — CGLab")
    print(f"  UniProt Version: {args.version}")
    print(f"  Output Dir:      {args.output_dir}")
    print(f"  Mode:            {args.mode}")
    print(f"  Streaming:       {'Yes (minimal disk)' if use_stream else 'No (full FASTA export)'}")
    print(f"  Cleanup FASTA:   {'Yes' if cleanup_fasta else 'No (keep files)'}")
    print(f"  Started:         {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'#'*60}")

    # ── Import Only ──
    if args.import_only:
        importer = HMMResultsImporter(args.version, args.output_dir)
        importer.import_results()
        return

    # ── Step 2: Prepare HMM profiles (do this FIRST in streaming mode) ──
    hmm_mgr = HMMProfileManager(args.output_dir)
    if args.skip_hmm_download:
        print("\nSkipping HMM download (--skip-hmm-download).")
        if not hmm_mgr.is_prepared():
            print("✗ HMM database not found and download was skipped.")
            sys.exit(1)
    else:
        print("\nPreparing Pfam HMM profiles...")
        if not hmm_mgr.prepare():
            print("✗ Failed to prepare HMM profiles.")
            sys.exit(1)

    hmm_file = str(hmm_mgr.hmm_file)
    results_dir = Path(args.output_dir) / "hmmscan_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1 + 3: Export + Search (together in streaming mode) ──
    if use_stream:
        print(f"\n{'='*60}")
        print(f"STREAMING MODE — minimal disk usage")
        print(f"  Chunks are fetched from MySQL one at a time,")
        print(f"  scanned with hmmscan, then immediately deleted.")
        print(f"  Only ~1 chunk (~{args.chunk_size * 200 / 1e6:.0f} MB) on disk at any time.")
        print(f"{'='*60}\n")

        streamer = SequenceStreamer(
            version=args.version,
            output_dir=args.output_dir,
            taxon_ids=args.taxon_ids,
            proteome_ids=args.proteome_ids,
        )
        total = streamer.get_total_count()

        if total == 0:
            print("✗ No proteins found for the given filters.")
            sys.exit(1)

        n_chunks = (total + args.chunk_size - 1) // args.chunk_size
        print(f"  Total proteins: {total:,}")
        print(f"  Total chunks:   {n_chunks:,}")
        print(f"  Chunk size:     {args.chunk_size:,}\n")

        # Fetch the first chunk to seed the pipeline, then keep N in flight
        # so hmmscan doesn't starve while waiting for MySQL.
        # Strategy: maintain a small pool of fetched-but-not-yet-scanned chunks
        # to overlap MySQL I/O with hmmscan computation.

        completed = 0
        failed = 0
        start_time = time.time()
        chunk_num = 1
        chunk_files_queue = []  # (chunk_num, fasta_path) ready to scan
        prefetch = min(args.parallel_jobs * 2, 4)  # pre-fetch up to 4 chunks

        def fetch_chunks_up_to(target_num):
            """Fetch chunks from MySQL until we have target_num ready."""
            nonlocal chunk_num  # <--- THIS IS THE FIX
            while len(chunk_files_queue) < target_num and chunk_num <= n_chunks:
                tmp_path, count = streamer.stream_chunk(chunk_num, args.chunk_size)
                if tmp_path and count > 0:
                    chunk_files_queue.append((chunk_num, tmp_path))
                chunk_num += 1

        with ProcessPoolExecutor(max_workers=args.parallel_jobs) as executor:
            active_futures = {}

            while completed < n_chunks or active_futures:
                # Pre-fetch chunks to keep the queue full
                slots_available = args.parallel_jobs - len(active_futures)
                if slots_available > 0:
                    fetch_chunks_up_to(slots_available + len(chunk_files_queue))

                # Submit new work from the queue
                while chunk_files_queue and len(active_futures) < args.parallel_jobs:
                    c_num, fasta_path = chunk_files_queue.pop(0)
                    output_file = results_dir / f"chunk_{c_num:05d}_hmmscan.domtblout.txt"
                    task = {
                        "fasta_file": fasta_path,
                        "hmm_file": hmm_file,
                        "output_file": str(output_file),
                        "evalue": args.evalue,
                        "dom_evalue": args.dom_evalue,
                        "threads": args.threads,
                        "cleanup": cleanup_fasta,
                    }
                    future = executor.submit(run_single_hmmscan, task)
                    active_futures[future] = c_num

                # Wait for at least one to finish
                if active_futures:
                    done_futures = []
                    for future in as_completed(active_futures):
                        c_num = active_futures[future]
                        result = future.result()
                        done_futures.append(future)
                        completed += 1

                        if result["success"]:
                            status = "✓"
                        else:
                            status = f"✗ ({result.get('error', 'unknown')[:50]})"
                            failed += 1

                        elapsed = time.time() - start_time
                        rate = completed / elapsed if elapsed > 0 else 0
                        eta = (n_chunks - completed) / rate if rate > 0 else 0

                        # Estimate proteins done
                        prot_done = min(completed * args.chunk_size, total)

                        print(
                            f"  [{completed}/{n_chunks}] {status} "
                            f"chunk_{c_num:05d} "
                            f"({prot_done:,}/{total:,} proteins) "
                            f"— {rate:.2f} chunks/min, ETA: {eta/60:.0f}min"
                        )
                        break  # process one at a time to refill the queue

                    for f in done_futures:
                        del active_futures[f]

        elapsed = time.time() - start_time
        print(f"\n✓ hmmscan complete: {completed - failed}/{n_chunks} succeeded "
              f"({failed} failed) in {elapsed/60:.1f}min")

        if cleanup_fasta:
            # Clean up any remaining temp files
            remaining = list(Path(args.output_dir).glob("tmp_fasta/*.fasta"))
            for f in remaining:
                f.unlink()
            if remaining:
                print(f"  Cleaned up {len(remaining)} remaining temp files.")

        # Import results
        importer = HMMResultsImporter(args.version, str(results_dir))
        importer.import_results()

        print(f"\n{'#'*60}")
        print(f"  Pipeline Complete!")
        print(f"  Duration: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'#'*60}")
        return

    # ── NON-STREAMING: export all FASTA first, then search ──
    exporter = SequenceExporter(
        version=args.version,
        output_dir=args.output_dir,
        taxon_ids=args.taxon_ids,
        proteome_ids=args.proteome_ids,
    )
    chunk_files = exporter.export(chunk_size=args.chunk_size)

    if not chunk_files:
        sys.exit(1)

    if args.export_only:
        print("\n✓ Export complete (--export-only mode).")
        return

    # ── Step 3: Run hmmscan ──
    if args.mode == "local":
        stats = run_local_search(
            fasta_dir=Path(args.output_dir) / "fasta_chunks",
            hmm_file=hmm_file,
            output_dir=str(results_dir),
            threads=args.threads,
            parallel_jobs=args.parallel_jobs,
            evalue=args.evalue,
            dom_evalue=args.dom_evalue,
            cleanup=cleanup_fasta,
        )
        if not stats.get("success") and stats.get("failed", 0) > stats.get("succeeded", 0):
            print(f"\n⚠ {stats['failed']} chunks failed. Check output files.")

        if cleanup_fasta:
            fasta_dir_to_clean = Path(args.output_dir) / "fasta_chunks"
            cleaned = 0
            for f in fasta_dir_to_clean.glob("chunk_*.fasta"):
                f.unlink()
                cleaned += 1
            print(f"  Cleaned up {cleaned} FASTA chunk files.")

    # ── Step 4: Import results ──
    importer = HMMResultsImporter(args.version, str(results_dir))
    importer.import_results()

    print(f"\n{'#'*60}")
    print(f"  Pipeline Complete!")
    print(f"  Duration: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'#'*60}")


if __name__ == "__main__":
    main()
