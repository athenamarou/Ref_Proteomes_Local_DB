# 🧬 UniProt Reference Proteome Pipeline

A high-performance toolkit for downloading, parsing, and querying **UniProt Reference Proteomes** into a versioned local MySQL database. Designed for bioinformatics labs running self-hosted infrastructure (e.g. Synology NAS + MySQL).

---

## 📋 Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Database Schema](#database-schema)
- [Usage](#usage)
  - [Sync Pipeline](#1-sync-pipeline--uniprot_sync)
  - [Retrieval Tool](#2-retrieval-tool--get_reference_uniprot_set)
- [Performance Notes](#performance-notes)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

This pipeline automates the full lifecycle of UniProt reference proteome data:

1. **Downloads** the versioned `.tar.gz` archive from the UniProt FTP server
2. **Parses** all canonical `.dat.gz` proteome files (SwissProt flat-file format), extracting protein metadata, sequences, GO terms, and Pfam domain annotations
3. **Loads** everything into a normalized, versioned MySQL database with sequence deduplication
4. **Queries** allow flexible retrieval of protein sets by version, taxonomy, proteome ID, GO term, or Pfam domain — exported as FASTA files

Supports multiple UniProt release versions side-by-side (e.g. `2025_05`, `2026_01`), making it easy to track changes across releases.

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
    │                         extracts: proteins, GO terms, Pfam domains
    ▼
DataBaseManager            ← bulk upserts into 8-table MySQL schema
    │                         sequence deduplication via SHA-256 hash
    ▼
MySQL Database
    │
    ▼
UniProtRetriever           ← flexible query interface → FASTA export
```

---

## Requirements

- Python 3.9+
- MySQL 8.0+ (local or remote, e.g. Synology NAS)
- ~50 GB disk space per UniProt release (tar.gz ~25 GB + DB)

### Python Dependencies

```
biopython>=1.81
mysql-connector-python>=8.0
requests>=2.28
python-dotenv>=1.0
```

See [`requirements.txt`](requirements.txt) for pinned versions.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/uniprot-proteome-pipeline.git
cd uniprot-proteome-pipeline

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure environment variables
cp .env.example .env
# Edit .env with your database credentials
```

---

## Configuration

Create a `.env` file in the project root (never commit this file):

```dotenv
DB_HOST=192.168.1.100       # MySQL host (IP or hostname)
DB_USER=uniprot_user        # MySQL username
DB_PASSWORD=your_password   # MySQL password
DB_NAME=uniprot_db          # Target database name
```

A `.env.example` template is included in the repository.

> **Synology NAS users:** Set `DB_HOST` to your NAS's local IP. See [Performance Notes](#performance-notes) for tips on running the pipeline remotely against a NAS-hosted database.

---

## Database Schema

The pipeline uses an **8-table normalized schema**:

| Table | Description |
|---|---|
| `proteins` | Core protein metadata (accession, name, organism, taxon_id, proteome_id, version) |
| `sequences` | Deduplicated amino acid sequences keyed by SHA-256 hash (`seq_id`) |
| `protein_go` | Many-to-many: proteins ↔ Gene Ontology terms |
| `protein_pfam` | Many-to-many: proteins ↔ Pfam domain annotations (with position + e-value) |
| *(+ 4 supporting tables)* | Indexes, versioning metadata, etc. |

Sequence deduplication is handled automatically — identical sequences across proteomes share a single row in the `sequences` table, reducing storage significantly.

---

## Usage

### 1. Sync Pipeline — `uniprot_sync`

Downloads and loads a specific UniProt release into the database.

```bash
# Basic sync (downloads if not present, skips if version already in DB)
python uniprot_sync.py -version 2026_01

# Force re-sync even if version already exists
python uniprot_sync.py -version 2026_01 --force

# Tune batch size for your hardware (default: 50,000)
python uniprot_sync.py -version 2026_01 --batch-size 100000
```

**Example output:**
```
============================================================
UniProt Reference Proteome Pipeline v8 — 2026_01
Started: 2026-01-28 14:23:01
Batch size: 50,000
============================================================

Streaming from 2026_01 archive...

  [0:05:30] Proteomes: 420 | Proteins: 1,250,000 | Rate: 3,800/sec | Cache: 94,230 | Current: UP000005640
  ...

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

---

### 2. Retrieval Tool — `get_reference_uniprot_set`

Query the database and export protein sets as FASTA files.

```bash
# List all available versions in your database
python get_reference_uniprot_set.py -version 2026_01 --list-versions

# List all proteome IDs for a version
python get_reference_uniprot_set.py -version 2026_01 --list-proteomes

# Retrieve by Proteome ID (e.g. Human reference proteome)
python get_reference_uniprot_set.py -version 2026_01 --proteome-id UP000005640

# Retrieve by Taxonomy ID
python get_reference_uniprot_set.py -version 2026_01 -taxonomy 9606

# Retrieve multiple taxa
python get_reference_uniprot_set.py -version 2026_01 -taxonomy 9606 10090 10116

# Filter by GO term
python get_reference_uniprot_set.py -version 2026_01 --go-id GO:0005634

# Filter by Pfam domain
python get_reference_uniprot_set.py -version 2026_01 --pfam-id PF00870

# Combine filters (e.g. human kinases)
python get_reference_uniprot_set.py -version 2026_01 -taxonomy 9606 --pfam-id PF00069
```

Output is a `.fasta` file named `uniprot_<identifier>_<version>.fasta` in the current directory.

**FASTA header format:**
```
>P04637 Cellular tumor antigen p53 OX=9606 UP=UP000005640
```

---

## Performance Notes

| Scenario | Estimated Time | Notes |
|---|---|---|
| Full sync, local MySQL | ~2–3 hours | Modern laptop/desktop |
| Full sync, NAS-hosted MySQL | ~4–6 hours | CPU bottleneck is gzip + BioPython parsing |
| Full sync, NAS CPU (ARM/Celeron) | ~12–24 hours | Run pipeline on Mac/PC, point DB_HOST to NAS |

**Recommended:** If your NAS has a weak CPU, run `uniprot_sync.py` on your laptop with `DB_HOST` pointing to the NAS IP. Python decompression and parsing are CPU-bound; SQL inserts over LAN add minimal latency.

**Batch size tuning:**
- Low RAM machine: `--batch-size 10000`
- Default (16 GB+ RAM): `--batch-size 50000`
- High-memory server: `--batch-size 100000`

---

## Project Structure

```
uniprot-proteome-pipeline/
├── uniprot_sync.py              # Sync pipeline (download + parse + load)
├── get_reference_uniprot_set.py # Query + FASTA export tool
├── requirements.txt
├── .env.example                 # Environment variable template
├── .gitignore
└── README.md
```

---

## Contributing

Contributions are welcome! Please open an issue to discuss proposed changes before submitting a pull request.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Commit your changes (`git commit -m 'Add my feature'`)
4. Push to the branch (`git push origin feature/my-feature`)
5. Open a Pull Request

---

## License

This project is licensed under the **MIT License** — see [`LICENSE`](LICENSE) for details.

**Data attribution:** UniProt data is provided by the [UniProt Consortium](https://www.uniprot.org/) under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license. © 2002–2026 UniProt Consortium.
