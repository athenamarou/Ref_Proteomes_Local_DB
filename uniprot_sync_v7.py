"""
UniProt Reference Proteome Pipeline — v8
=========================================
Downloads, parses, and loads UniProt reference proteomes (.tar.gz archive)
into a versioned MySQL database with sequence deduplication.

v8 Changes from v7:
────────────────────
1. BULK sequence deduplication — replaces per-protein SELECT with per-BATCH
queries. This is the single biggest speedup (~10-50x on DB operations).

2. MySQL session tuning — temporarily disables FK checks and adjusts flush
behavior during bulk loading (2-5x faster inserts).

3. Optimised parsing — reads each .dat file once instead of building a list
then re-joining it.

4. Progress reporting — prints rate (proteins/sec), cache size, and
proteome count every 30 seconds so you can monitor long runs.

5. Larger default batch size (50 000 vs 5 000 in v6 / 20 000 in v7).

Performance tip:
────────────────
If the Synology NAS has a weak CPU (Celeron / ARM), consider running this
script on your Mac and pointing DB_HOST in .env to the NAS IP.  Python's
gzip decompression and BioPython parsing are CPU-bound and will be 5-10×
faster on a modern laptop.  Network latency for the SQL inserts is
comparable to local access on a NAS.

Usage:
    python uniprot_sync_v8.py -version 2026_01
    python uniprot_sync_v8.py -version 2026_01 --batch-size 100000 --force
"""

import requests
from Bio import SeqIO
import gzip
import os
import mysql.connector
import datetime
from dotenv import load_dotenv
import logging
import argparse
import sys
import hashlib
import tarfile
import io
import time


# ─── SYNOLOGY PATHS ────────────────────────────────────────────
BASE_PATH = "/mnt/cglab.shared/Data/DBs/Uniprot/"
LOG_FILE = os.path.join(BASE_PATH, "update_history.log")
LOCAL_DATA_FILE = os.path.join(BASE_PATH, "Reference_Proteomes_2026_01.tar.gz")

script_dir = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(script_dir, ".env"))

logging.basicConfig(
    filename=LOG_FILE, level=logging.INFO, format="%(asctime)s - %(message)s"
)


# ════════════════════════════════════════════════════════════════
#  DOWNLOADER — handles the .tar.gz archive from UniProt FTP
# ════════════════════════════════════════════════════════════════

class UniProtDownloader:
    """
    Downloads (if needed) the Reference_Proteomes tar.gz archive and
    streams individual .dat.gz proteome files from within it.
    """

    def __init__(self, url, local_filename, version):
        self.url = url
        self.local_filename = local_filename
        self.version = version

    def get_versioned_filename(self):
        """Generate filename like Reference_Proteomes_2026_01.tar.gz"""
        return os.path.join(
            os.path.dirname(self.local_filename),
            f"Reference_Proteomes_{self.version}.tar.gz",
        )

    def check_for_updates(self):
        """Return False if we already have this version on disk."""
        versioned_file = self.get_versioned_filename()

        if os.path.exists(versioned_file):
            print(f"Version {self.version} already exists locally: {versioned_file}")
            self.local_filename = versioned_file
            return False

        try:
            # HEAD request: fetch only headers (file size) without downloading
            response = requests.head(self.url)
            if response.status_code == 200:
                remote_size = int(response.headers.get("Content-Length", 0))
                print(f"Remote file size: {remote_size / (1024**3):.2f} GB")
                return True
            else:
                print(f"Warning: Could not check remote file (status {response.status_code})")
                return True
        except requests.RequestException as e:
            print(f"Warning: Could not check remote file: {e}")
            return True

    def download_file(self):
        """Stream-download the archive in 1 MB chunks with progress bar."""
        versioned_file = self.get_versioned_filename()
        temp_file = versioned_file + ".tmp"

        print(f"Downloading UniProt reference .tar.gz for version {self.version}...")
        print(f"Target: {versioned_file}")

        try:
            response = requests.get(self.url, stream=True)
            response.raise_for_status()
            total_size = int(response.headers.get("Content-Length", 0))
            downloaded = 0

            with open(temp_file, "wb") as f:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(
                                f"\rProgress: {percent:.1f}% "
                                f"({downloaded / (1024**3):.2f} GB)",
                                end="",
                            )

            print("\nDownload complete!")
            os.rename(temp_file, versioned_file)
            self.local_filename = versioned_file

        except Exception as e:
            print(f"\n✗ Download failed: {e}")
            if os.path.exists(temp_file):
                os.remove(temp_file)
            raise

    def stream_tar_contents(self):
        """
        Generator that yields (proteome_id, file_handle) for each canonical
        .dat.gz file inside the archive.

        The tar.gz layout (from the UniProt README) is:
            Archaea/UP000000242/UP000000242_2234.dat.gz
            Bacteria/UP000001234/UP000001234_83333.dat.gz
            Eukaryota/UP000005640/UP000005640_9606.dat.gz
            ...

        We skip:
        - files starting with "._"
        - Isoform/additional files (*_additional*) — only canonical sequences
        """
        versioned_file = self.get_versioned_filename()

        with tarfile.open(versioned_file, "r:gz") as tar:
            for member in tar:
                if (
                    member.isfile()
                    and member.name.endswith(".dat.gz")
                    and not member.name.startswith("._")
                    and "_additional" not in member.name
                ):
                    # Extract the UP-prefixed proteome ID from the path
                    path_parts = member.name.split("/")
                    proteome_id = None
                    for part in path_parts:
                        if part.startswith("UP"):
                            proteome_id = part
                            break

                    if not proteome_id:
                        continue

                    f = tar.extractfile(member)
                    if f:
                        # Each .dat.gz inside the tar is itself gzip-compressed
                        gz = gzip.GzipFile(fileobj=f)
                        yield proteome_id, io.TextIOWrapper(gz, encoding="utf-8")


# ════════════════════════════════════════════════════════════════
#  PARSER — extracts protein records from UniProt flat files
# ════════════════════════════════════════════════════════════════

class UniprotParser:
    """
    Parses UniProt Swiss-Prot flat-file (.dat) format.

    Each entry in a .dat file has structured lines:
        ID   P53_HUMAN    ...        ← entry name
        AC   P04637;                 ← accession number
        DE   RecName: ...            ← description
        OS   Homo sapiens            ← organism
        OX   NCBI_TaxID=9606         ← taxonomy
        DR   GO; GO:0005634; ...     ← cross-ref to Gene Ontology
        DR   Pfam; PF00870; ...      ← cross-ref to Pfam domain DB
        SQ   SEQUENCE ...            ← sequence header
             MEEPQSDPSV ...          ← amino-acid sequence
        //                           ← end of entry

    BioPython's SeqIO.parse(handle, "swiss") parses the full record
    structure.  However, it doesn't expose DR lines (GO/Pfam) in a
    convenient way, so we do a fast line-scan first to grab those,
    then let BioPython handle the rest.
    """

    def __init__(self, file_path=None):
        self.file_path = file_path
        self.annotations = {}
        self.go_map = {}
        self.pfam_map = {}

    def parse_stream(self, proteome_id, handle):
        """
        Parse a single proteome's .dat file and yield protein records.

        v8 improvement: We call handle.read() once to get the full text,
        then split it for annotation scanning and wrap in StringIO for
        BioPython.  v6/v7 called readlines() and then re-joined the
        list with "".join(lines), which created two full copies of the
        data in memory.

        Individual proteome .dat files range from a few KB (viruses with
        1-3 proteins) to ~100 MB (large eukaryotes like wheat with
        100K+ proteins), so holding one in memory is fine.
        """
        # Read entire file content once
        text = handle.read()

        # Fast line scan to extract GO and Pfam annotations
        self.annotations, self.go_map, self.pfam_map = (
            self.extract_annotations_from_text(text.splitlines(keepends=True))
        )

        # Full parse with BioPython (accession, sequence, organism, taxon, etc.)
        for record in SeqIO.parse(io.StringIO(text), "swiss"):
            taxon_ids = record.annotations.get("ncbi_taxid", [None])
            taxon_id = taxon_ids[0] if taxon_ids else None
            accession = record.id

            go_terms = self.annotations.get(accession, {}).get("go_terms", [])
            pfam_domains = self.annotations.get(accession, {}).get("pfam_domains", [])

            yield {
                "accession": accession,
                "name": record.name,
                "sequence": str(record.seq),
                "organism": record.annotations.get("organism", ""),
                "taxon_id": taxon_id,
                "proteome_id": proteome_id,
                "go_terms": go_terms,
                "pfam_domains": pfam_domains,
            }

    def extract_annotations_from_text(self, lines):
        """
        Scan raw .dat lines to extract GO and Pfam cross-references.

        DR (Database Reference) line examples:
            DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.
            DR   Pfam; PF00870; P53; 1.

        Returns three dicts:
            annotations  — {accession: {go_terms: [...], pfam_domains: [...]}}
            go_map       — {go_id: go_name}    (for the go_terms master table)
            pfam_map     — {pfam_id: pfam_name} (for the pfam_domains master table)
        """
        annotations = {}
        go_map = {}
        pfam_map = {}
        current_acc = None

        for line in lines:
            if line.startswith("AC   "):
                # AC lines can list multiple accessions separated by ";".
                # The first one is the primary accession.
                current_acc = line[5:].split(";")[0].strip()
                annotations[current_acc] = {"go_terms": [], "pfam_domains": []}

            elif line.startswith("DR   GO;") and current_acc:
                parts = line.split(";")
                if len(parts) >= 3:
                    go_id = parts[1].strip()
                    go_name = parts[2].strip()
                    if go_id.startswith("GO:"):
                        annotations[current_acc]["go_terms"].append(go_id)
                        go_map[go_id] = go_name

            elif line.startswith("DR   Pfam;") and current_acc:
                parts = line.split(";")
                if len(parts) >= 3:
                    pfam_id = parts[1].strip()
                    pfam_name = parts[2].strip()
                    if pfam_id.startswith("PF"):
                        annotations[current_acc]["pfam_domains"].append(pfam_id)
                        pfam_map[pfam_id] = pfam_name

        return annotations, go_map, pfam_map

    def get_records(self):
        """Legacy method: parse a standalone .dat.gz file (not from tar)."""
        print(f"Parsing {self.file_path}...")

        print("Extracting GO and Pfam annotations...")
        with gzip.open(self.file_path, "rt") as handle:
            self.annotations, self.go_map, self.pfam_map = (
                self.extract_annotations_from_text(handle)
            )

        print("Parsing protein records...")
        with gzip.open(self.file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "swiss"):
                taxon_ids = record.annotations.get("ncbi_taxid", [None])
                taxon_id = taxon_ids[0] if taxon_ids else None
                accession = record.id
                go_terms = self.annotations.get(accession, {}).get("go_terms", [])
                pfam_domains = self.annotations.get(accession, {}).get(
                    "pfam_domains", []
                )

                yield {
                    "accession": accession,
                    "name": record.name,
                    "sequence": str(record.seq),
                    "organism": record.annotations.get("organism", ""),
                    "taxon_id": taxon_id,
                    "go_terms": go_terms,
                    "pfam_domains": pfam_domains,
                }


# ════════════════════════════════════════════════════════════════
#  DATABASE MANAGER — MySQL operations with bulk optimisations
# ════════════════════════════════════════════════════════════════

class DataBaseManager:
    """
    Manages the 8-table MySQL schema and all insert operations.

    v8 key improvements:
    ────────────────────
    • BULK sequence deduplication: resolves all hashes in a batch with
      2–3 SQL queries instead of N individual round-trips.
    • MySQL session tuning: disables FK/unique checks during bulk load.
    • In-memory seq_cache avoids re-querying previously seen sequences.
    """

    def __init__(self, config, current_version):
        self.conn = mysql.connector.connect(**config)
        self.cursor = self.conn.cursor()
        self.version = current_version

        # ── In-memory cache: sequence_hash (str) → seq_id (int) ──
        #
        # Purpose: avoid hitting MySQL for sequences we've already resolved
        # in earlier batches.  On a first full load most sequences are new,
        # but shared orthologs across proteomes will produce cache hits.

        self.seq_cache = {}
        self.MAX_SEQ_CACHE = 100_000_000

    # ─── MySQL Bulk-Load Tuning ─────────────────────────────────

    def enable_bulk_load_mode(self):
        """
        Temporarily relax MySQL safety settings for faster inserts.

        These are SESSION-level changes (affect only this connection), so
        they won't impact other users of the database.

        foreign_key_checks = 0
            Normally MySQL verifies every FK on INSERT (e.g. checking that
            seq_id exists in the sequences table before inserting a protein).
            We control the insert order, so this check is redundant.

        unique_checks = 0
            Skips secondary unique-index verification during INSERT.
            Our INSERT IGNORE / ON DUPLICATE KEY already handles conflicts.

        innodb_flush_log_at_trx_commit = 2
            Default (=1): flush redo log to disk on EVERY COMMIT (safest).
            Value 2: flush once per second (OS crash loses ≤1 sec of data).
            For a bulk import you'd re-run anyway, so this is safe.

        Combined effect: roughly 2–5× faster INSERT throughput.
        """
        print("  Enabling bulk-load session optimizations...")
        self.cursor.execute("SET SESSION foreign_key_checks = 0")
        self.cursor.execute("SET SESSION unique_checks = 0")
        #self.cursor.execute("SET SESSION innodb_flush_log_at_trx_commit = 2")
        self.conn.commit()

    def disable_bulk_load_mode(self):
        """Restore normal MySQL safety settings after the load."""
        print("  Restoring MySQL session defaults...")
        self.cursor.execute("SET SESSION foreign_key_checks = 1")
        self.cursor.execute("SET SESSION unique_checks = 1")
        #self.cursor.execute("SET SESSION innodb_flush_log_at_trx_commit = 1")
        self.conn.commit()

    # ─── Bulk Sequence Deduplication ────────────────────────────

    def get_or_create_sequences_bulk(self, records):
        """
        Resolve sequence → seq_id for every record in a batch, using
        bulk SQL queries instead of one query per protein.


        Algorithm:
        ──────────
        1. Compute MD5 hash for every sequence in the batch.
        2. Check in-memory cache (instant, no SQL).
        3. Bulk SELECT from DB for any cache misses (1 query per chunk).
        4. Bulk INSERT truly new sequences (1 query per chunk).
        5. Bulk SELECT to retrieve auto-increment IDs for new rows.
        6. Return list of seq_ids aligned with input records.

        Parameters:
            records — list of protein dicts (must contain "sequence" key)

        Returns:
            list of int — seq_id for each record, same order as input
        """
        # Step 1: Compute MD5 for each sequence and deduplicate within batch
        # Two proteins can share the same sequence (paralogs, shared genes
        # across proteomes), so we use a dict to avoid duplicate work.
        hash_to_seq = {}      # {md5_hash: sequence_string}  — unique only
        record_hashes = []    # [md5_hash, ...]               — one per record

        for rec in records:
            h = hashlib.md5(rec["sequence"].encode()).hexdigest()
            record_hashes.append(h)
            if h not in hash_to_seq:
                hash_to_seq[h] = rec["sequence"]

        # Step 2: Identify hashes not yet in our in-memory cache
        uncached_hashes = [h for h in hash_to_seq if h not in self.seq_cache]

        # Step 3: Bulk SELECT from MySQL for uncached hashes
        # We chunk the IN (...) clause because MySQL has a packet-size limit
        # (~16 MB by default).  10 000 hashes × 32 chars ≈ 320 KB — well
        # within limits.
        CHUNK = 10_000

        if uncached_hashes:
            for i in range(0, len(uncached_hashes), CHUNK):
                chunk = uncached_hashes[i : i + CHUNK]
                placeholders = ",".join(["%s"] * len(chunk))
                self.cursor.execute(
                    f"SELECT seq_id, sequence_hash FROM sequences "
                    f"WHERE sequence_hash IN ({placeholders})",
                    chunk,
                )
                for seq_id, seq_hash in self.cursor.fetchall():
                    self.seq_cache[seq_hash] = seq_id

        # Step 4: Bulk INSERT sequences that are genuinely new
        new_hashes = [h for h in uncached_hashes if h not in self.seq_cache]

        if new_hashes:
            insert_data = [
                (hash_to_seq[h], h, self.version) for h in new_hashes
            ]
            # INSERT IGNORE: if another process inserted the same hash between
            # our SELECT (step 3) and now, we silently skip the duplicate.
            self.cursor.executemany(
                """INSERT IGNORE INTO sequences
                   (sequence, sequence_hash, first_seen_version)
                   VALUES (%s, %s, %s)""",
                insert_data,
            )

            # Step 5: Fetch back auto-increment IDs for the rows we just inserted.
            # executemany doesn't return all lastrowids, so we re-SELECT.
            for i in range(0, len(new_hashes), CHUNK):
                chunk = new_hashes[i : i + CHUNK]
                placeholders = ",".join(["%s"] * len(chunk))
                self.cursor.execute(
                    f"SELECT seq_id, sequence_hash FROM sequences "
                    f"WHERE sequence_hash IN ({placeholders})",
                    chunk,
                )
                for seq_id, seq_hash in self.cursor.fetchall():
                    self.seq_cache[seq_hash] = seq_id

        # Prevent unbounded memory growth on very large runs
        mapped_ids = [self.seq_cache[h] for h in record_hashes]

        # Prevent unbounded memory growth on very large runs
        if len(self.seq_cache) > self.MAX_SEQ_CACHE:
            print(f"\n  [Cache] Reached {self.MAX_SEQ_CACHE:,} entries, clearing...")
            self.seq_cache.clear()

        return mapped_ids

    # ─── Schema Creation ────────────────────────────────────────

    def create_table(self):
        """
        Create all 8 tables in the correct FK dependency order.
        IF NOT EXISTS means this is idempotent (safe to run repeatedly).

        Table overview:
        ───────────────
        1. sequences       — deduplicated amino-acid sequences (MD5 hash key)
        2. proteomes       — one row per reference proteome (UP ID → taxon)
        3. proteins        — versioned protein metadata (accession + version = PK)
        4. sequence_changes— tracks when a protein's sequence changes across versions
        5. go_terms        — Gene Ontology master table
        6. protein_go      — links proteins to GO terms (version-specific)
        7. pfam_domains    — Pfam domain master table
        8. protein_pfam    — links proteins to Pfam domains (version-specific)
        """
        # 1. Sequences — the deduplication backbone
        # MEDIUMTEXT supports up to 16 MB (titin is ~35 000 aa ≈ 35 KB, so fine)
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS sequences (
            seq_id INT AUTO_INCREMENT PRIMARY KEY,
            sequence MEDIUMTEXT NOT NULL,
            sequence_hash CHAR(32) UNIQUE NOT NULL,
            first_seen_version VARCHAR(10),
            INDEX idx_hash (sequence_hash)
        )
        """)

        # 2. Proteomes
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS proteomes (
            proteome_id VARCHAR(20) PRIMARY KEY,
            taxon_id INT,
            description TEXT,
            INDEX idx_taxon (taxon_id)
        )
        """)

        # 3. Proteins — composite PK (accession, version) allows tracking
        #    the same protein across different UniProt releases
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS proteins (
            accession VARCHAR(20),
            version VARCHAR(10),
            name VARCHAR(50),
            organism VARCHAR(255),
            taxon_id INT,
            seq_id INT NOT NULL,
            proteome_id VARCHAR(20),
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
            PRIMARY KEY (accession, version),
            INDEX idx_version (version),
            INDEX idx_taxon (taxon_id),
            INDEX idx_proteome (proteome_id),
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id),
            FOREIGN KEY (proteome_id) REFERENCES proteomes(proteome_id)
        )
        """)

        # 4. Sequence changes (version-to-version tracking)
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS sequence_changes (
            change_id INT AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(20),
            from_version VARCHAR(10),
            to_version VARCHAR(10),
            old_seq_id INT,
            new_seq_id INT,
            change_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (old_seq_id) REFERENCES sequences(seq_id),
            FOREIGN KEY (new_seq_id) REFERENCES sequences(seq_id)
        )
        """)

        # 5. GO terms master
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS go_terms (
            go_id VARCHAR(20) PRIMARY KEY,
            go_name VARCHAR(255),
            go_namespace VARCHAR(50),
            definition TEXT,
            INDEX idx_go_name (go_name)
        )
        """)

        # 6. Protein ↔ GO linking (version-specific)
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS protein_go (
            accession VARCHAR(20),
            version VARCHAR(10),
            go_id VARCHAR(20),
            evidence_code VARCHAR(10),
            PRIMARY KEY (accession, version, go_id),
            FOREIGN KEY (accession, version) REFERENCES proteins(accession, version),
            FOREIGN KEY (go_id) REFERENCES go_terms(go_id),
            INDEX idx_accession (accession),
            INDEX idx_go_id (go_id),
            INDEX idx_version (version)
        )
        """)

        # 7. Pfam domains master
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS pfam_domains (
            pfam_id VARCHAR(20) PRIMARY KEY,
            pfam_name VARCHAR(255),
            pfam_description TEXT,
            INDEX idx_pfam_name (pfam_name)
        )
        """)

        # 8. Protein ↔ Pfam linking (version-specific)
        self.cursor.execute("""
        CREATE TABLE IF NOT EXISTS protein_pfam (
            accession VARCHAR(20),
            version VARCHAR(10),
            pfam_id VARCHAR(20),
            start_position INT,
            end_position INT,
            e_value FLOAT,
            PRIMARY KEY (accession, version, pfam_id, start_position),
            FOREIGN KEY (accession, version) REFERENCES proteins(accession, version),
            FOREIGN KEY (pfam_id) REFERENCES pfam_domains(pfam_id),
            INDEX idx_accession (accession),
            INDEX idx_pfam_id (pfam_id),
            INDEX idx_version (version)
        )
        """)

        self.conn.commit()
        print("  Database schema initialized (8 tables).")

    # ─── Version Check ──────────────────────────────────────────

    def check_version_exists(self):
        """Check if this version already has data in the database."""
        query = "SELECT COUNT(*) FROM proteins WHERE version = %s"
        self.cursor.execute(query, (self.version,))
        count = self.cursor.fetchone()[0]
        return count > 0

    # ─── Batch Upsert ───────────────────────────────────────────

    def upsert_batch(self, batch, go_map=None, pfam_map=None):
        """
        Insert a batch of protein records into all relevant tables.

        v8 change: calls get_or_create_sequences_bulk() to resolve all
        sequence → seq_id mappings in 2–3 queries, instead of one query
        per protein (v6) or one query per cache-miss (v7).

        Insertion order respects FK dependencies:
            sequences → proteomes → proteins → protein_go / protein_pfam
        """
        # ── Bulk resolve sequences (THE key optimisation) ──
        seq_ids = self.get_or_create_sequences_bulk(batch)

        protein_data = []
        proteome_data = {}
        all_go_links = []
        all_go_master = set()
        all_pfam_links = []
        all_pfam_master = set()

        for i, record in enumerate(batch):
            seq_id = seq_ids[i]

            proteome_id = record.get("proteome_id", "Unknown")
            if proteome_id not in proteome_data:
                proteome_data[proteome_id] = record["taxon_id"]

            protein_data.append((
                record["accession"],
                self.version,
                record["name"],
                record["organism"],
                record["taxon_id"],
                seq_id,
                proteome_id,
            ))

            for go_id in record.get("go_terms", []):
                all_go_master.add(go_id)
                all_go_links.append(
                    (record["accession"], self.version, go_id, "IEA")
                )

            for pfam_id in record.get("pfam_domains", []):
                all_pfam_master.add(pfam_id)
                all_pfam_links.append(
                    (record["accession"], self.version, pfam_id, 0, 0, 0.0)
                )

        try:
            # A: Proteomes
            if proteome_data:
                self.cursor.executemany(
                    """INSERT IGNORE INTO proteomes (proteome_id, taxon_id)
                       VALUES (%s, %s)""",
                    [(pid, tid) for pid, tid in proteome_data.items()],
                )

            # B: Proteins
            if protein_data:
                self.cursor.executemany(
                    """INSERT IGNORE INTO proteins
                       (accession, version, name, organism, taxon_id, seq_id, proteome_id)
                       VALUES (%s, %s, %s, %s, %s, %s, %s)""",
                    protein_data,
                )

            # C: GO terms (master + links)
            if all_go_master and go_map:
                go_batch = [
                    (gid, go_map.get(gid, "Unknown")) for gid in all_go_master
                ]
                self.cursor.executemany(
                    """INSERT INTO go_terms (go_id, go_name)
                       VALUES (%s, %s)
                       ON DUPLICATE KEY UPDATE go_name = VALUES(go_name)""",
                    go_batch,
                )

            if all_go_links:
                self.cursor.executemany(
                    """INSERT IGNORE INTO protein_go
                       (accession, version, go_id, evidence_code)
                       VALUES (%s, %s, %s, %s)""",
                    all_go_links,
                )

            # D: Pfam domains (master + links)
            if all_pfam_master and pfam_map:
                pfam_batch = [
                    (pid, pfam_map.get(pid, "Unknown")) for pid in all_pfam_master
                ]
                self.cursor.executemany(
                    """INSERT INTO pfam_domains (pfam_id, pfam_name)
                       VALUES (%s, %s)
                       ON DUPLICATE KEY UPDATE pfam_name = VALUES(pfam_name)""",
                    pfam_batch,
                )

            if all_pfam_links:
                self.cursor.executemany(
                    """INSERT IGNORE INTO protein_pfam
                       (accession, version, pfam_id, start_position,
                        end_position, e_value)
                       VALUES (%s, %s, %s, %s, %s, %s)""",
                    all_pfam_links,
                )

            self.conn.commit()

        except mysql.connector.Error as err:
            print(f"\nError during batch upload: {err}")
            self.conn.rollback()
            raise

    # ─── Statistics ──────────────────────────────────────────────

    def get_version_stats(self):
        """Return counts for this version across all tables."""
        queries = {
            "total": "SELECT COUNT(*) FROM proteins WHERE version = %s",
            "organisms": "SELECT COUNT(DISTINCT organism) FROM proteins WHERE version = %s",
            "taxa": "SELECT COUNT(DISTINCT taxon_id) FROM proteins WHERE version = %s",
            "proteomes": "SELECT COUNT(DISTINCT proteome_id) FROM proteins WHERE version = %s",
            "unique_sequences": "SELECT COUNT(DISTINCT seq_id) FROM proteins WHERE version = %s",
            "go_terms": "SELECT COUNT(DISTINCT go_id) FROM protein_go WHERE version = %s",
            "pfam_domains": "SELECT COUNT(DISTINCT pfam_id) FROM protein_pfam WHERE version = %s",
        }

        stats = {}
        for key, query in queries.items():
            self.cursor.execute(query, (self.version,))
            result = self.cursor.fetchone()
            stats[key] = result[0] if result else 0

        if stats["total"] > 0:
            stats["dedup_ratio"] = (
                1 - stats["unique_sequences"] / stats["total"]
            ) * 100
        else:
            stats["dedup_ratio"] = 0

        return stats

    # ─── Cleanup ─────────────────────────────────────────────────

    def close(self):
        """Close the MySQL connection safely."""
        if self.conn.is_connected():
            self.cursor.close()
            self.conn.close()
            print("  MySQL connection closed.")


# ════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════

load_dotenv()


def main():
    parser = argparse.ArgumentParser(
        description="Download and sync UniProt reference proteomes from .tar.gz archive",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-version", required=True, help="UniProt version (e.g., 2026_01)"
    )
    parser.add_argument(
        "--force", action="store_true", help="Force re-sync even if version exists"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50000,
        help="Records per DB batch (default: 50000)",
    )

    args = parser.parse_args()

    DB_CONFIG = {
        "host": os.getenv("DB_HOST"),
        "user": os.getenv("DB_USER"),
        "password": os.getenv("DB_PASSWORD"),
        "database": os.getenv("DB_NAME"),
    }

    URL = f"https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Reference_Proteomes_{args.version}.tar.gz"



    start_time = datetime.datetime.now()
    print(f"\n{'='*60}")
    print(f"UniProt Reference Proteome Pipeline v8 — {args.version}")
    print(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Batch size: {args.batch_size:,}")
    print(f"{'='*60}\n")

    db = None

    try:
        # ── Database setup ──
        db = DataBaseManager(DB_CONFIG, args.version)
        db.create_table()

        if db.check_version_exists() and not args.force:
            print(f"Version {args.version} already exists in database.")
            stats = db.get_version_stats()
            print(f"  Total proteins:   {stats['total']:,}")
            print(f"  Unique proteomes: {stats['proteomes']:,}")
            return

        # ── Download if needed ──
        dl = UniProtDownloader(URL, LOCAL_DATA_FILE, args.version)
        if dl.check_for_updates():
            dl.download_file()

        # ── Enable bulk-load MySQL optimizations ──
        db.enable_bulk_load_mode()

        # ── Stream and parse ──
        uniprot_parser = UniprotParser()
        total_count = 0
        proteome_count = 0
        batch = []

        # Progress tracking
        pipeline_start = time.time()
        last_report_time = pipeline_start
        last_report_count = 0

        print(f"\nStreaming from {args.version} archive...\n")

        for proteome_id, file_handle in dl.stream_tar_contents():
            proteome_count += 1

            for record in uniprot_parser.parse_stream(proteome_id, file_handle):
                batch.append(record)

                if len(batch) >= args.batch_size:
                    db.upsert_batch(
                        batch,
                        go_map=uniprot_parser.go_map,
                        pfam_map=uniprot_parser.pfam_map,
                    )
                    total_count += len(batch)
                    batch = []

                    # ── Progress report every 30 seconds ──
                    now = time.time()
                    if now - last_report_time >= 30:
                        elapsed = now - last_report_time
                        rate = (total_count - last_report_count) / elapsed
                        total_elapsed = now - pipeline_start
                        print(
                            f"  [{datetime.timedelta(seconds=int(total_elapsed))}] "
                            f"Proteomes: {proteome_count:,} | "
                            f"Proteins: {total_count:,} | "
                            f"Rate: {rate:,.0f}/sec | "
                            f"Cache: {len(db.seq_cache):,} | "
                            f"Current: {proteome_id}"
                        )
                        last_report_time = now
                        last_report_count = total_count

        # Flush remaining records
        if batch:
            db.upsert_batch(
                batch,
                go_map=uniprot_parser.go_map,
                pfam_map=uniprot_parser.pfam_map,
            )
            total_count += len(batch)

        # ── Restore MySQL defaults ──
        db.disable_bulk_load_mode()

        # ── Final stats ──
        stats = db.get_version_stats()
        duration = datetime.datetime.now() - start_time

        print(f"\n{'='*60}")
        print("Pipeline Completed Successfully")
        print(f"{'='*60}")
        print(f"  Total proteins:       {stats['total']:,}")
        print(f"  Unique proteomes:     {stats['proteomes']:,}")
        print(f"  Unique sequences:     {stats['unique_sequences']:,}")
        print(f"  Sequence dedup ratio: {stats['dedup_ratio']:.2f}%")
        print(f"  GO terms linked:      {stats['go_terms']:,}")
        print(f"  Pfam domains linked:  {stats['pfam_domains']:,}")
        print(f"  Total duration:       {duration}")
        print(f"{'='*60}\n")

        logging.info(
            f"SUCCESS: {args.version} synced. "
            f"{total_count} proteins, {proteome_count} proteomes in {duration}"
        )

    except Exception as e:
        print(f"\nCritical Error: {e}")
        logging.error(f"FAILURE: {args.version} failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)

    finally:
        if db:
            # Ensure safety checks are restored even on error
            try:
                db.disable_bulk_load_mode()
            except Exception:
                pass
            db.close()


if __name__ == "__main__":
    main()
