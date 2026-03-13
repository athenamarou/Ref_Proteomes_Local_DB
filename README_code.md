# UniProt Reference Proteome Pipeline

## Overview

This repository contains two scripts for building and querying a local versioned MySQL database of UniProt canonical reference proteomes.

- `uniprot_sync_v9.py` — Downloads, parses, and loads a full UniProt reference proteome release into the database.
- `get_reference_uniprot_set_v5.py` — Retrieves protein sets from the database with flexible filtering options.

The database stores canonical sequences only (no isoforms), normalized across 9 relational tables with sequence deduplication, GO term and Pfam domain annotations, and full version tracking across UniProt releases.

---

## Requirements

### Python Dependencies

```
mysql-connector-python
python-dotenv
requests
```

Install with:

```bash
pip install mysql-connector-python python-dotenv requests
```

### Environment

Create a `.env` file in the same directory as the scripts:

```
DB_HOST=localhost
DB_USER=your_db_user
DB_PASSWORD=your_db_password
DB_NAME=uniprot_db_cglab
```

### Data

The pipeline expects the UniProt reference proteome archive on local disk:

```
/mnt/cglab.shared/Data/DBs/Uniprot/Reference_Proteomes_<version>.tar.gz
```

If the file is not present, the script will attempt to download it from the UniProt FTP server. If the file is stored on a network mount (NAS), it is strongly recommended to copy it to local disk before running the pipeline to avoid network-related failures during multi-day runs.

---

## Database Schema

The schema consists of 9 tables organized around sequence deduplication and normalized metadata.

### Table Dependency Order (FK chain)

```
taxa
sequences
    |
proteomes (references taxa)
    |
proteins (references sequences, proteomes, taxa)
    |
protein_go    (references proteins, go_terms)
protein_pfam  (references proteins, pfam_domains)

go_terms      (standalone master)
pfam_domains  (standalone master)
sequence_changes (references sequences)
```

### Table Descriptions

**taxa**
Stores one row per NCBI taxon ID. Organism name and taxonomic division (Archaea, Bacteria, Eukaryota, Viruses) are stored here rather than on the proteins table to avoid redundant string storage across millions of rows.

| Column | Type | Description |
|---|---|---|
| taxon_id | INT PRIMARY KEY | NCBI Taxonomy ID |
| organism | VARCHAR(255) | Full organism name |
| division | ENUM | Archaea, Bacteria, Eukaryota, or Viruses |
| lineage | TEXT | Reserved for future lineage storage |

**sequences**
Stores unique amino acid sequences deduplicated by MD5 hash. A single sequence shared across multiple proteomes or organisms is stored only once.

| Column | Type | Description |
|---|---|---|
| seq_id | INT AUTO_INCREMENT | Internal sequence ID |
| sequence | MEDIUMTEXT | Full amino acid sequence |
| sequence_hash | CHAR(32) UNIQUE | MD5 hash for deduplication |
| first_seen_version | VARCHAR(10) | UniProt release where first observed |

**proteomes**
One row per reference proteome (UP-prefixed ID).

| Column | Type | Description |
|---|---|---|
| proteome_id | VARCHAR(20) PRIMARY KEY | UniProt Proteome ID (e.g. UP000005640) |
| taxon_id | INT | NCBI Taxonomy ID |
| description | TEXT | Reserved for future use |

**proteins**
Versioned protein records. The composite primary key (accession, version) allows the same protein to be tracked across multiple UniProt releases.

| Column | Type | Description |
|---|---|---|
| accession | VARCHAR(20) | UniProt accession (e.g. P04637) |
| version | VARCHAR(10) | UniProt release version (e.g. 2026_01) |
| name | VARCHAR(50) | Entry name (e.g. P53_HUMAN) |
| taxon_id | INT | FK to taxa |
| seq_id | INT | FK to sequences |
| proteome_id | VARCHAR(20) | FK to proteomes |
| last_updated | TIMESTAMP | Auto-updated on row change |

**go_terms / protein_go**
Master GO term vocabulary and per-protein, per-version GO annotations.

**pfam_domains / protein_pfam**
Master Pfam domain vocabulary and per-protein, per-version Pfam domain assignments.

**sequence_changes**
Tracks accessions whose sequences changed between versions. Populated externally (not by the sync pipeline itself).

---

## uniprot_sync_v9.py

### Usage

```bash
python uniprot_sync_v9.py -version 2026_01
python uniprot_sync_v9.py -version 2026_01 --batch-size 100000
python uniprot_sync_v9.py -version 2026_01 --force
```

### Arguments

| Argument | Description |
|---|---|
| `-version` | UniProt release version string, required (e.g. 2026_01) |
| `--batch-size` | Records per database batch commit, default 50000 |
| `--force` | Re-run even if this version already exists in the database |

### Pipeline Stages

**1. Schema initialization**
`create_table()` creates all 9 tables using `IF NOT EXISTS`, making it safe to call on every run. Table creation respects the FK dependency order — taxa and sequences are created before the tables that reference them.

**2. Archive download**
If the versioned `.tar.gz` file is not found on disk, the script downloads it from:
```
https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/
```
Downloads use streaming with 1 MB chunks and a temp file to prevent partial files on failure.

**3. Checkpoint loading**
On startup, the script reads `checkpoint_<version>.txt` from `BASE_PATH` if it exists. Any proteome ID listed there is skipped during the current run. This allows interrupted runs to resume without reprocessing already-loaded proteomes.

**4. MySQL bulk-load optimizations**
Before streaming begins, the following session-level settings are applied:
- `foreign_key_checks = 0` — skips FK validation on insert (order is controlled)
- `unique_checks = 0` — skips secondary index checks during insert
- `innodb_flush_log_at_trx_commit = 2` — flushes redo log once per second instead of per commit (commented out by default; safe to enable for bulk imports)

**5. Streaming and parsing**
The archive is opened with `tarfile` and streamed member by member. Each `.dat.gz` proteome file is decompressed on the fly with `gzip.GzipFile` and parsed line by line without loading the full file into memory. The parser is a single-pass manual implementation that extracts:

- Entry name (ID line)
- Primary accession (AC line)
- Organism name (OS line)
- NCBI TaxID (OX line)
- GO cross-references (DR GO lines)
- Pfam cross-references (DR Pfam lines)
- Amino acid sequence (SQ block)

Only canonical `.dat.gz` files are processed. Files matching `*_additional*` (isoforms) are skipped.

**6. Batch upsert**
Records accumulate in memory until `batch_size` is reached. Each batch is committed in a single transaction covering all 9 tables in FK dependency order:

```
taxa -> sequences -> proteomes -> proteins -> go_terms -> protein_go -> pfam_domains -> protein_pfam
```

Sequence deduplication is handled by bulk MD5 hash resolution — all sequences in a batch are resolved in 2-3 SQL queries rather than one per protein.

**7. Checkpointing**
After each proteome's records are fully committed, its ID is appended to the checkpoint file. The checkpoint write happens after the commit, never before, guaranteeing consistency between the file and the database state.

**8. Progress reporting**
Every 30 seconds the pipeline prints:
```
[HH:MM:SS] Proteomes: N | Proteins: N | Rate: N/sec | Cache: N | Current: UP000XXXXXX
```

### Performance Notes

The in-memory sequence cache (`seq_cache`) stores resolved hash-to-seq_id mappings across batches. The cache ceiling is set to 100,000,000 entries. On servers with limited RAM, reducing this to 5,000,000 has negligible effect on deduplication quality while significantly reducing memory usage.

The parsing step is CPU-bound. On the same server that hosts MySQL, CPU contention is not a concern. If running the script remotely against a network-hosted MySQL instance, network latency on `executemany` batch inserts is negligible relative to parsing time.

---

## get_reference_uniprot_set_v5.py

### Usage

```bash
python get_reference_uniprot_set_v5.py -version 2026_01
python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606
python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606 10090
python get_reference_uniprot_set_v5.py -version 2026_01 --proteome-id UP000005640
python get_reference_uniprot_set_v5.py -version 2026_01 --go-id GO:0005634
python get_reference_uniprot_set_v5.py -version 2026_01 --pfam-id PF00870
python get_reference_uniprot_set_v5.py -version 2026_01 --division Eukaryota
python get_reference_uniprot_set_v5.py -version 2026_01 --list-versions
python get_reference_uniprot_set_v5.py -version 2026_01 --list-proteomes
```

### Arguments

| Argument | Description |
|---|---|
| `-version` | UniProt release version, required |
| `-taxonomy` | One or more NCBI Taxonomy IDs (space-separated) |
| `--proteome-id` | Filter by proteome ID (e.g. UP000005640) |
| `--go-id` | Filter by GO term (e.g. GO:0005634) |
| `--pfam-id` | Filter by Pfam domain (e.g. PF00870) |
| `--division` | Filter by division: Archaea, Bacteria, Eukaryota, or Viruses |
| `--list-versions` | Print all versions present in the database |
| `--list-proteomes` | Print all proteome IDs for the specified version |

### Output

Results are exported as a FASTA file in the working directory:

```
uniprot_<identifier>_<version>.fasta
```

FASTA headers follow standard UniProt format:

```
>ACCESSION ENTRY_NAME OS=Organism name OX=TaxonID UP=ProteomeID DIV=Division
```

Filters can be combined. For example, taxonomy + GO ID returns proteins from the specified taxon that are annotated with the specified GO term.

---

## Checkpoint File

The checkpoint file is a plain text file, one proteome ID per line, written to `BASE_PATH`:

```
checkpoint_2026_01.txt
```

To force a complete re-run from scratch, delete this file before running with `--force`. To resume an interrupted run, simply re-run the script without deleting the checkpoint — already-completed proteomes will be skipped automatically.

---

## Log File

All successful and failed runs are recorded to:

```
/mnt/cglab.shared/Data/DBs/Uniprot/update_history.log
```

Format:
```
YYYY-MM-DD HH:MM:SS - SUCCESS: 2026_01 synced. N proteins, N proteomes in HH:MM:SS
YYYY-MM-DD HH:MM:SS - FAILURE: 2026_01 failed: <error message>
```
