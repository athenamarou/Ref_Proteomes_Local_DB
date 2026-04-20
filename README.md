## UniProt Reference Proteome Local Database

‼️**Status: Active Development / Under Optimization**

A pipeline for downloading, parsing, profiling, and querying versioned **UniProt reference proteomes** into a versioned local MySQL database. Designed for bioinformatics labs running self-hosted infrastructure (e.g. Synology NAS + MySQL).

---

## Overview

The pipeline consists of two independent tools:

- **`uniprot_sync_v7.py`** — downloads the UniProt reference proteome archive and streams it into a local MySQL database.
- **`get_reference_uniprot_set_lib.py`** — retrieves protein sets from that database with flexible filters and exports them as FASTA files.
- **`pyhmmer_hmmsearch.py`** - integrates MySQL queries and pyHMMER to perform hmmsearch across all the accessions for the local Reference Proteome Database built.

The database schema is versioned, meaning multiple UniProt releases can coexist in the same instance without overwriting previous data. Sequence storage uses MD5-based deduplication so identical sequences shared across proteomes or versions are stored only once.

---

## Infrastructure

| Component | Role | Path |
| --- | --- | --- |
| Synology NAS | Stores the `.tar.gz` archive | `/mnt/.../Uniprot/` |
| Server (`/var/`) | Runs MySQL; houses the database | configured via `.env` |
| `.env` file | Holds DB credentials | see setup below |

The sync script reads the archive from the NAS mount and writes parsed data to the MySQL instance on the server. Both scripts connect to MySQL using credentials from a shared `.env` file.

---

## Requirements

- Python 3.8+
- MySQL 8.0+
- BioPython (`biopython`)
- `mysql-connector-python`
- `python-dotenv`
- `requests`
- `pyhmmer`
- `psutil`

Install dependencies:

```bash
pip install biopython mysql-connector-python python-dotenv requests pyhmmer psutil`
```


---

## Configuration

Create a `.env` file in the same directory as both scripts:

```
DB_HOST=192.168.1.100       # MySQL host (IP or hostname)
DB_USER=uniprot_user        # MySQL username
DB_PASSWORD=your_password   # MySQL password
DB_NAME=uniprot_db          # Target database name`
```

For read only rights for the team , a new user can be created and then access the database like this:

```
mysql -u user -p
```

To create the new user the config block inside get_reference_uniprot_set_lib.py script can be updated like this:

```python
    # Load DB Config with Read-Only Defaults for the Lab
    DB_CONFIG = {
        "host": os.getenv("DB_HOST", "localhost"),
        "user": os.getenv("DB_USER", "new_user"),
        "password": os.getenv("DB_PASSWORD", "lab_password"), # Change to your actual lab password
        "database": os.getenv("DB_NAME", "uniprot_db"),
    }

    retriever = UniProtRetriever(DB_CONFIG)

    try:
        retriever.connect()
```

In SQL, connect as root and proceed:

```sql
CREATE USER 'cglab_user'@'localhost' IDENTIFIED BY 'lab_password';
GRANT SELECT ON uniprot_db_cglab.* TO 'cglab_user'@'localhost';
FLUSH PRIVILEGES;`
```

The read-only users can connect in the database like this, by giving the password:

```sql
mysql -u cglab_user -p`
```

The sync script also reads two hardcoded paths at the top of the file that should match your NAS mount:

```python
BASE_PATH = "/mnt/.../Uniprot/"
LOCAL_DATA_FILE = os.path.join(BASE_PATH, "Reference_Proteomes_2026_01.tar.gz")`
```

---

## Database Schema

The pipeline creates 8 tables in the correct foreign-key dependency order. All `CREATE TABLE` statements use `IF NOT EXISTS`, so the schema initialisation is idempotent. Sequence deduplication is handled automatically — identical sequences across proteomes share a single row in the `sequences` table, reducing storage significantly.

**Database Schema Overview**

```
sequences          — deduplicated amino-acid sequences (MD5 hash, auto-increment seq_id)
proteomes          — one row per reference proteome (UP-prefixed ID → taxon_id)
proteins           — versioned protein metadata; PRIMARY KEY (accession, version)
sequence_changes   — records when a protein's sequence changes between versions
go_terms           — Gene Ontology master table (go_id, go_name, namespace, definition)
protein_go         — protein ↔ GO term links, version-specific
pfam_domains       — Pfam domain master table (pfam_id, pfam_name, description)
protein_pfam       — protein ↔ Pfam domain links with positional and e-value data
hmm_search_results — optional analytical table for pyHMMER search results
```

Foreign key relationships:
```mermaid
erDiagram
  sequences {
    int seq_id PK
    char md5_hash
    text sequence
  }
  proteomes {
    varchar proteome_id PK
    int taxon_id
    varchar scientific_name
    varchar assembly
  }
  proteins {
    varchar accession PK
    int version PK
    int seq_id FK
    varchar proteome_id FK
    varchar gene_name
    varchar protein_name
    int length
    timestamp updated_at
  }
  sequence_changes {
    int change_id PK
    varchar accession FK
    int old_version
    int new_version
    int old_seq_id FK
    int new_seq_id FK
    timestamp changed_at
  }
  go_terms {
    varchar go_id PK
    varchar go_name
    varchar namespace
    text definition
  }
  protein_go {
    varchar accession FK
    int version FK
    varchar go_id FK
    varchar evidence_code
    varchar source
  }
  pfam_domains {
    varchar pfam_id PK
    varchar pfam_name
    text description
    varchar clan_id
  }
  protein_pfam {
    varchar accession FK
    int version FK
    varchar pfam_id FK
    int seq_start
    int seq_end
    float e_value
    float bit_score
  }
hmm_search_results {
    bigint id PK
    varchar version FK
    varchar accession FK
    int taxon_id
    varchar proteome_id
    varchar protein_name
    varchar hmm_name
    varchar hmm_accession
    varchar hmm_type
    double full_evalue
    double full_score
    double full_bias
    int domain_number
    int domain_count
    double domain_evalue
    double domain_score
    double domain_bias
    int hmm_from
    int hmm_to
    int ali_from
    int ali_to
    int env_from
    int env_to
    float hmm_coverage
    float protein_coverage
    float posterior_prob
    timestamp search_date
  }

  proteomes ||--o{ proteins : "contains"
  sequences ||--o{ proteins : "seq_id"
  proteins ||--o{ sequence_changes : "accession"
  sequences ||--o{ sequence_changes : "old_seq_id"
  proteins ||--o{ protein_go : "accession+version"
  go_terms ||--o{ protein_go : "go_id"
  proteins ||--o{ protein_pfam : "accession+version"
  pfam_domains ||--o{ protein_pfam : "pfam_id"
  proteins ||--o{ hmm_search_results : "accession+version"`

```
---

## Architecture

```
UniProt FTP
│
▼
UniProtDownloader          ← streams Reference_Proteomes_XXXX_XX.tar.gz
│
▼
UniprotParser              ← parses .dat.gz flat files (BioPython + line scan)
│                          extracts: proteins, GO terms, Pfam domains
▼
DataBaseManager            ← bulk upserts into MySQL schema
│                          sequence deduplication via MD5 hash
▼
MySQL Database
│
├──► pyhmmer_hmmsearch     ← [Optional] Multi-threaded HMM profiling
│
▼
UniProtRetriever           ← flexible query interface → FASTA export`
```

---

## Usage

### 1. Sync a UniProt version

```bash
# Basic sync (downloads if not present, skips if version already in DB)
python uniprot_sync.py -version 2026_01

# Force re-sync even if version already exists
python uniprot_sync.py -version 2026_01 --force

# Tune batch size for your hardware (default: 50,000)
python uniprot_sync.py -version 2026_01 --batch-size 100000`
```

If the archive already exists at the configured NAS path, the download step is skipped and parsing begins immediately. If the version already exists in the database, the script exits early unless `--force` is passed.

**Options:**

| **Flag** | **Default** | **Description** |
| --- | --- | --- |
| `-version` | required | UniProt release string, e.g. `2026_01` |
| `--batch-size` | 50000 | Number of records committed per database transaction |
| `--force` | off | Re-sync even if this version already exists in the database |

**Progress output** (printed every 30 seconds during a long run):

```
============================================================
UniProt Reference Proteome Pipeline v8 — 2026_01
Started: 2026-01-28 14:23:01
Batch size: 50,000
============================================================

Streaming from 2026_01 archive...

  [0:05:30] Proteomes: 420 | Proteins: 1,250,000 | Rate: 3,800/sec | Cache: 94,230 | Current: UP000005640
  ...

```

**Final summary on completion Example:**

```
============================================================
Pipeline Completed Successfully
============================================================
  Total proteins:       226,452,210
  Unique proteomes:     34,230
  Unique sequences:     198,340,120
  Sequence dedup ratio: 12.41%
  GO terms linked:      88,234,102
  Pfam domains linked:  71,892,445
  Total duration:       3:42:18
============================================================
```

A log entry is written to `update_history.log` in `BASE_PATH` on both success and failure.

---

### 2. Optional Downstream Analysis: HMM Domain Profiling

Once the local UniProt database is populated, you can run a high-performance profile Hidden Markov Model (HMM) search against the sequences to identify specific protein domains. This process utilizes `pyhmmer` to parallelize the search across multiple CPU cores without writing sequences to disk, ensuring maximum throughput.

Ensure your `.env` file is configured and that you have downloaded the required HMM profiles (e.g., `Pfam-A.hmm`) into your output directory structure.

```bash
python pyhmmer_hmmsearch.py \
  --version 2026_01 \
  --output-dir /path/to/your/output \
  --chunk-size 50000`
```

**Optional Arguments:**

- `-taxon-ids`: Restrict the search to specific taxonomy IDs (e.g., `10090` for Mouse, `9606` for Human).
- `-proteome-ids`: Restrict the search to specific proteomes.

**Output:** This script does not generate flat files. It automatically creates and populates the `hmm_search_results` table in the MySQL database, tracking accession IDs, HMM hits, E-values, and coordinate mappings for rapid downstream SQL querying.

---

### 3. Retrieve a protein set
List all loaded versions:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --list-versions
```

Retrieve the full human reference proteome:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606
```

Retrieve proteins from multiple taxa:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090 10116
```

Retrieve all proteins with a Homeodomain hit (all taxa):

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain
```

Retrieve human Homeodomain proteins with a strict E-value:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain -taxonomy 9606 --evalue 1e-10
```

Retrieve by Pfam accession (prefix match — PF00046 matches PF00046.36):

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name PF00046
```

Write output to a specific directory:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain --output-dir ./fastas
```

**Output:** A FASTA file in the current directory. Filename format: `uniprot_<identifier>_<version>.fasta`

FASTA header format:

```
>P04637 P53_HUMAN OX=9606 UP=UP000005640
```

---


## Performance Notes

- The sync script uses bulk MySQL session tuning (`foreign_key_checks = 0`, `unique_checks = 0`) during the load, which provides roughly 2–5x faster insert throughput. These are session-level changes and do not affect other connections.
- Sequence deduplication uses an in-memory MD5 hash cache (up to 100 million entries) to avoid redundant database round-trips. Batch lookups are chunked at 10,000 hashes per query to stay within MySQL packet limits.
- If the NAS CPU is a Celeron or ARM-based processor, running the sync script on a workstation and pointing `DB_HOST` in `.env` at the server IP will substantially improve parsing speed, as gzip decompression and BioPython parsing are CPU-bound.
- The `pyhmmer_hmmsearch.py` script implements an in-memory streaming queue with backpressure. Monitor your RAM usage carefully; if the `chunk-size` is set too high, the process may trigger an Out-Of-Memory (OOM) kill from the OS.

---

## File Structure

```
.
├── uniprot_sync_v7.py               # Sync pipeline
├── get_reference_uniprot_set_lib.py # Retrieval tool
├── pyhmmer_hmmsearch.py             # Optional HMM search pipeline
├── .env                             # DB credentials (not committed)
└── README.md`
```

The NAS archive and log file are stored outside this repository:

```
/mnt/cglab.shared/Data/DBs/Uniprot/
├── Reference_Proteomes_2026_01.tar.gz
└── update_history.log`
```

---

## Notes on the Archive Layout

The UniProt reference proteome archive is organised as:

```
Archaea/UP000000242/UP000000242_2234.dat.gz
Bacteria/UP000001234/UP000001234_83333.dat.gz
Eukaryota/UP000005640/UP000005640_9606.dat.gz
...
```

The sync script streams only canonical sequence files (`.dat.gz`), skipping `_additional` files (isoforms) and macOS metadata prefixes (`._`).
