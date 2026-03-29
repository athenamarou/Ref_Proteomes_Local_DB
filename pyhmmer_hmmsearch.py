# ======================================
#           NECESSARY LIBRARIES
# ======================================

import pyhmmer
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock

import argparse
import os
import time
from datetime import datetime
from pathlib import Path

import mysql.connector
from dotenv import load_dotenv

import psutil

import gc

import threading

import queue

# ======================================
#           CONFIGURATION
# ======================================

DEFAULT_CHUNK_SIZE = 50000  # proteins per MySQL fetch
PYHMMER_CPUS = 32

# Load .env from the script's directory
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(_SCRIPT_DIR, ".env"))


# ======================================
#           PYHMMER SEARCH
# ======================================


def run_pyhmmer_hmmsearch(hmm_models, protein_rows, num_cpus=32):
    # Initialize the alphabet
    alphabet = pyhmmer.easel.Alphabet.amino()

    # Create a fast lookup dictionary
    row_lookup = {row["accession"]: row for row in protein_rows}

    # Convert MYSQL row into DigitalSequence objects in RAM.
    sequences = []
    for row in protein_rows:
        seq = pyhmmer.easel.TextSequence(
            name=row["accession"].encode(), sequence=row["sequence"]
        )
        sequences.append(seq.digitize(alphabet))

    # Create a DigitalSequenceBlock (the database for the search)
    msa_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, sequences)

    # Execute the strict hmmsearch using Gathering Thresholds
    searcher = pyhmmer.hmmer.hmmsearch(
        hmm_models, msa_block, cpus=num_cpus, bit_cutoffs="gathering"
    )

    batch_results = []
    for hits in searcher:
        # Get HMM name safely (handles API changes and byte-strings)
        raw_hmm_name = hits.query.name
        safe_hmm_name = (
            raw_hmm_name.decode() if isinstance(raw_hmm_name, bytes) else raw_hmm_name
        )

        for hit in hits:
            acc = hit.name
            if isinstance(acc, bytes):
                acc = acc.decode()
            original_row = row_lookup[acc]

            # Now we loop through each individual domain found in this protein
            for i, domain in enumerate(hit.domains, 1):
                batch_results.append(
                    {
                        "accession": acc,
                        "taxon_id": original_row["taxon_id"],
                        "proteome_id": original_row["proteome_id"],
                        "protein_name": original_row["name"],
                        "hmm_name": safe_hmm_name,
                        "full_evalue": hit.evalue,
                        "full_score": hit.score,
                        "domain_number": i,
                        "domain_count": len(hit.domains),
                        "domain_evalue": domain.i_evalue,
                        "domain_score": domain.score,
                        "ali_from": domain.alignment.target_from,
                        "ali_to": domain.alignment.target_to,
                        "env_from": domain.env_from,
                        "env_to": domain.env_to,
                        "hmm_from": domain.alignment.hmm_from,
                        "hmm_to": domain.alignment.hmm_to,
                    }
                )
    return batch_results


# ======================================
#           SEQUENCE STREAMER
# ======================================


class SequenceStreamer:
    """
    Takes sequences from SQL Database
    and streams them in chunks for memory to process them.
    """

    def __init__(self, version, output_dir, taxon_ids=None, proteome_ids=None):
        self.version = version
        self.output_dir = Path(output_dir)
        self.taxon_ids = taxon_ids
        self.proteome_ids = proteome_ids

        self.config = {
            "host": os.getenv("DB_HOST", "localhost"),
            "user": os.getenv("DB_USER", "cglab_user"),
            "password": os.getenv("DB_PASSWORD", ""),
            "database": os.getenv("DB_NAME", "uniprot_db_cglab"),
        }

        # -----CHECKPOINTING-----
        self.checkpoint_file = self.output_dir / f"checkpoint_{version}.txt"
        self.last_accession = self._load_checkpoint()

    def _load_checkpoint(self):
        if self.checkpoint_file.exists():
            val = self.checkpoint_file.read_text().strip()
            if val:
                print(f"  Resuming from checkpoint: {val}")
                return val
        return None

    def _save_checkpoint(self):
        self.checkpoint_file.write_text(self.last_accession or "")

    # QUERY BUILDER
    def _build_query(self):
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

    def stream_chunk_to_memory(self, chunk_size):
        query, params = self._build_query()
        if self.last_accession:
            query += " AND p.accession > %s"
            params.append(self.last_accession)

        query += " ORDER BY p.accession LIMIT %s"
        page_params = params + [chunk_size]

        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor(dictionary=True)
        cursor.execute(query, page_params)
        rows = cursor.fetchall()
        cursor.close()
        conn.close()

        if rows:
            self.last_accession = rows[-1]["accession"]
        return rows


# ======================================
#           RESULTS IMPORTER
# ======================================


class HMMResultsImporter:
    """
    Creates results table in the DB and inserts
    results from list in my SQL DB.
    """

    def __init__(self, version, output_dir):
        self.version = version
        self.output_dir = Path(output_dir)
        self.config = {
            "host": os.getenv("DB_HOST", "localhost"),
            "user": os.getenv("DB_USER", "cglab_user"),
            "password": os.getenv("DB_PASSWORD", ""),
            "database": os.getenv("DB_NAME", "uniprot_db_cglab"),
        }

    # --------- CREATE RESULTS TABLE ---------
    def create_results_table(self):
        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()
        try:
            cursor.execute(
                """
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
            """
            )
            conn.commit()
            print("  Results table 'hmm_search_results' ready.")
        finally:
            cursor.close()
            conn.close()

    # --------- IMPORT RESULTS TO MYSQL ---------
    def import_list_to_mysql(self, batch_results, batch_size=50000):
        if not batch_results:
            return

        conn = mysql.connector.connect(**self.config)
        cursor = conn.cursor()

        insert_sql = """
            INSERT INTO hmm_search_results(
                version, accession, taxon_id, proteome_id, protein_name,
                hmm_name, full_evalue, full_score, hmm_type,
                domain_number, domain_count, domain_evalue, domain_score,
                ali_from, ali_to, hmm_from, hmm_to, env_from, env_to
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, 'Pfam', %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        try:
            data = [
                (
                    self.version,
                    r["accession"],
                    r["taxon_id"],
                    r["proteome_id"],
                    r["protein_name"],
                    r["hmm_name"],
                    r["full_evalue"],
                    r["full_score"],
                    r["domain_number"],
                    r["domain_count"],
                    r["domain_evalue"],
                    r["domain_score"],
                    r["ali_from"],
                    r["ali_to"],
                    r["hmm_from"],
                    r["hmm_to"],
                    r["env_from"],
                    r["env_to"],
                )  # Ensure these match the dictionary keys
                for r in batch_results
            ]

            for i in range(0, len(data), batch_size):
                chunk = data[i : i + batch_size]
                cursor.executemany(insert_sql, chunk)
                conn.commit()
        finally:
            cursor.close()
            conn.close()


# ======================================
#           MAIN PIPELINE
# ======================================


def main():
    parser = argparse.ArgumentParser(
        description="PyHMMER Search Pipeline for CGLab UniProt DB (In-Memory)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # -------- ARGUMENTS ----------
    parser.add_argument("--version", required=True, help="UniProt DB version")
    parser.add_argument(
        "--output-dir", required=True, help="Base directory for outputs"
    )
    parser.add_argument("--taxon-ids", nargs="+", type=int, default=None)
    parser.add_argument("--proteome-ids", nargs="+", default=None)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)

    args = parser.parse_args()

    print(f"\n{'#'*60}")
    print(f"  PyHMMER Search Pipeline — CGLab")
    print(f"  UniProt Version: {args.version}")
    print(f"  Output Dir:      {args.output_dir}")
    print(f"  CPUs Allocated:  {PYHMMER_CPUS}")
    print(f"  Started:         {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'#'*60}\n")

    # Load HMM profiles into memory using the correct path.
    HMM_DB_PATH = os.path.join(args.output_dir, "hmm_profiles", "Pfam-A.hmm")
    print(f"Loading HMM profiles from {HMM_DB_PATH} into memory..")
    with HMMFile(HMM_DB_PATH) as hmm_file:
        hmms = list(hmm_file)

    # Initialize streamer and importer
    streamer = SequenceStreamer(
        args.version, args.output_dir, args.taxon_ids, args.proteome_ids
    )
    importer = HMMResultsImporter(args.version, args.output_dir)
    importer.create_results_table()

    chunk_size = args.chunk_size
    total_processed = 0
    rate = 0
    start_time = time.time()

    # ── Three queues act as buffers between pipeline stages.
    # maxsize=2 is intentional backpressure: if inserter falls behind,
    # searcher blocks before fetching a 3rd chunk, capping memory usage
    # at roughly 3 chunks in flight simultaneously (~75k sequences at 25k chunk size).
    fetch_q = queue.Queue(maxsize=2)
    insert_q = queue.Queue(maxsize=2)

    def fetcher():
        # Runs independently of search/insert — keeps fetch_q pre-loaded
        # so pyhmmer never waits on MySQL. One lightweight I/O thread.
        while True:
            rows = streamer.stream_chunk_to_memory(chunk_size)
            fetch_q.put(rows)  # blocks here if fetch_q is full (backpressure working)
            if not rows:
                break  # empty list = DB exhausted; sentinel already in queue

    def searcher():
        # The CPU-heavy stage. pyhmmer internally uses PYHMMER_CPUS threads
        # to parallelize HMM search — this outer thread just orchestrates
        # fetching input and forwarding results; it does not add parallelism
        # to the search itself.
        while True:
            rows = fetch_q.get()  # blocks until fetcher has a chunk ready
            if not rows:
                insert_q.put(None)  # forward sentinel so inserter knows to stop
                break

            mem = psutil.virtual_memory()
            if mem.percent > 85:
                print(f"WARNING: RAM at {mem.percent:.1f}% — chunk may be too large")

            results = run_pyhmmer_hmmsearch(hmms, rows, num_cpus=PYHMMER_CPUS)

            n = len(rows)
            last_acc = rows[-1]["accession"]
            del rows  # release sequence memory before results enter insert_q;
            gc.collect()  # critical at 25k chunks — avoids two chunks overlapping in RAM

            insert_q.put((results, n, last_acc))

    def inserter():
        nonlocal total_processed  # shares the counter with main thread for progress reporting
        while True:
            item = insert_q.get()  # blocks until searcher has results ready
            if item is None:
                break  # clean exit — sentinel received, all chunks processed

            results, n, last_acc = item
            importer.import_list_to_mysql(results)

            # Checkpoint the accession belonging to this chunk specifically,

            streamer.last_accession = last_acc
            streamer._save_checkpoint()

            del results  # free result memory immediately after commit
            gc.collect()  # insert_q maxsize=2 means up to 2 result sets could sit here

            total_processed += n
            elapsed = time.time() - start_time
            rate = total_processed / elapsed if elapsed > 0 else 0
            mem = psutil.virtual_memory()
            print(
                f"Processed {total_processed:,} proteins... ({rate:,.0f} seq/s) "
                f"| RAM: {mem.percent:.1f}% ({mem.used/1e9:.1f}/{mem.total/1e9:.1f} GB)"
            )

    # fetcher and searcher are daemon threads — if main crashes they die automatically.
    # inserter is NOT daemon: even on interrupt, Python waits for it to finish
    # the current commit so you don't end up with a partial chunk in the DB.
    t1 = threading.Thread(target=fetcher, daemon=True)
    t2 = threading.Thread(target=searcher, daemon=True)
    t3 = threading.Thread(target=inserter)

    t1.start()
    t2.start()
    t3.start()

    # join() blocks main until each thread exits cleanly.
    # Order matters: t1 → t2 → t3 mirrors data flow,
    # so if you're debugging a hang you know which join() is stuck.
    t1.join()
    t2.join()
    t3.join()

    print(f"\nDone. Total processed: {total_processed:,}")


if __name__ == "__main__":
    main()
