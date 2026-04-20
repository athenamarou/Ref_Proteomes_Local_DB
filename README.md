## UniProt Reference Proteome Local Database

‼️**Status: Active Development / Under Optimization**

A two-script pipeline for downloading, parsing, and querying versioned **UniProt reference proteomes** into a versioned local MySQL database. Designed for bioinformatics labs
running self-hosted infrastructure (e.g. Synology NAS + MySQL).

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
pip install biopython mysql-connector-python python-dotenv requests pyhmmer psutil
