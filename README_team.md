# UniProt Local Database — Team Documentation

## Purpose

This document describes the local UniProt reference proteome database maintained by the Computational Genomics Laboratory. It is intended for lab members who need to retrieve protein sequences or annotations from the database, and for anyone responsible for maintaining or updating it.

---

## What the Database Contains

The database stores canonical reference proteomes from UniProt for one or more versioned releases. A reference proteome is a proteome selected by UniProt as representative for a given species — covering well-studied model organisms and organisms of biomedical or phylogenetic interest.

The current loaded release is **2026_01** (published January 28, 2026), which covers approximately 34,000 species across all domains of life and contains over 200 million canonical protein sequences.

The following information is stored per protein:

- Accession number and entry name
- Amino acid sequence (deduplicated — identical sequences across species share a single stored copy)
- Organism and NCBI Taxonomy ID
- Taxonomic division (Archaea, Bacteria, Eukaryota, Viruses)
- Proteome ID (the UniProt UP-prefixed identifier for the source proteome)
- GO term annotations
- Pfam domain annotations
- Version of the UniProt release in which the protein appears

Isoforms and additional sequences are not stored. Only canonical sequences are loaded.

---

## Database Access

### Connection Details

The database runs on the lab server. Connection credentials are stored in the `.env` file in the scripts directory. Do not share or commit this file.

```
Host:     [lab server hostname or IP]
Database: uniprot_db_cglab
User:     uniprot_user
```

Contact the database administrator if you need access credentials.

### Querying the Database Directly

Standard MySQL client:

```bash
mysql -h <host> -u uniprot_user -p uniprot_db_cglab
```

---

## Retrieving Protein Sets

Use `get_reference_uniprot_set_v5.py` to export protein sequences as FASTA files. This script queries the database and writes output to the current working directory.

### Basic Usage

All queries require a version argument:

```bash
python get_reference_uniprot_set_v5.py -version 2026_01
```

### Filtering Options

**By species (NCBI Taxonomy ID)**

Single species:
```bash
python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606
```

Multiple species in one query:
```bash
python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606 10090 7955
```

Common taxonomy IDs:

| Organism | Taxonomy ID |
|---|---|
| Homo sapiens | 9606 |
| Mus musculus | 10090 |
| Rattus norvegicus | 10116 |
| Danio rerio | 7955 |
| Drosophila melanogaster | 7227 |
| Caenorhabditis elegans | 6239 |
| Saccharomyces cerevisiae | 559292 |
| Arabidopsis thaliana | 3702 |
| Escherichia coli K-12 | 83333 |

**By proteome ID**

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --proteome-id UP000005640
```

UP000005640 is the human reference proteome.

**By GO term**

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --go-id GO:0005634
```

GO:0005634 is the nucleus. This returns all proteins annotated to that GO term in the specified version.

**By Pfam domain**

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --pfam-id PF00870
```

**By taxonomic division**

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --division Bacteria
```

Valid values: Archaea, Bacteria, Eukaryota, Viruses.

**Combining filters**

Filters can be combined. For example, all human proteins annotated to a specific Pfam domain:

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606 --pfam-id PF00870
```

### Utility Commands

List all versions currently loaded in the database:

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --list-versions
```

List all proteome IDs available for a version:

```bash
python get_reference_uniprot_set_v5.py -version 2026_01 --list-proteomes
```

### Output Format

The output is a FASTA file written to the current working directory:

```
uniprot_<identifier>_<version>.fasta
```

FASTA header format:

```
>ACCESSION ENTRY_NAME OS=Organism name OX=TaxonID UP=ProteomeID DIV=Division
```

Example:

```
>P04637 P53_HUMAN OS=Homo sapiens OX=9606 UP=UP000005640 DIV=Eukaryota
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNLNKEYILQKLIDISQTPAIQKLQVNQTLKKKDDTLNIFSEDIFGLRNYNSNPDQKPIVHSQEMIDLQTKEDSSDSHYNMTCPQLKFLNREKDSVSQSRQDLQRGEADLLNQLKNQALQKPSVTVNLYQDQLSEDDVAHYKEVHREQKKEENLEAKLNRMKENLTKVLENLEQSPQKPTCVNQFHQKPLHNLWEDLMSQMRETPVKQKFHNSQLQKMDLQDLSYERQ
```

---

## Database Updates

UniProt releases a new version approximately every 8 weeks. When a new release is available, the database is updated by running the sync pipeline. The pipeline is additive — new version data is inserted alongside existing versions without removing previous data.

### Update Procedure

1. Download or confirm the new `.tar.gz` archive is present on the server:
   ```
   /mnt/cglab.shared/Data/DBs/Uniprot/Reference_Proteomes_<version>.tar.gz
   ```
   Copy it to a local disk path before running the pipeline if it is stored on the NAS mount.

2. Run the sync pipeline on the server:
   ```bash
   python uniprot_sync_v9.py -version 2026_03
   ```

3. Monitor progress. The pipeline reports status every 30 seconds:
   ```
   [02:14:33] Proteomes: 8420 | Proteins: 48,200,000 | Rate: 3,200/sec | Cache: 12,400,000 | Current: UP000005640
   ```

4. A full load of a new release takes approximately 3-7 days depending on server load. Do not interrupt the process. If the process is interrupted for any reason, restart it with the same command — the checkpoint system will resume from the last completed proteome.

5. Verify the load completed successfully:
   ```bash
   python get_reference_uniprot_set_v5.py -version 2026_03 --list-versions
   ```

### If the Pipeline Fails

If the pipeline crashes mid-run:

- Do not drop the database.
- Do not delete the checkpoint file.
- Identify the error from the terminal output or the log file at:
  ```
  /mnt/cglab.shared/Data/DBs/Uniprot/update_history.log
  ```
- Fix the underlying issue if applicable.
- Restart the pipeline with the same command. The checkpoint system will skip all already-loaded proteomes and resume from where it stopped.

Using `--force` is only necessary if you want to reload a version that is already marked as complete in the database. Do not use `--force` during a recovery restart.

---

## Storage

### Database

The database is hosted on the lab server. Expected storage for a full 2026_01 load:

| Table | Approximate size |
|---|---|
| sequences | 200-250 GB |
| proteins | 30-40 GB |
| protein_go | 15-20 GB |
| protein_pfam | 10-15 GB |
| taxa, proteomes, go_terms, pfam_domains | < 1 GB combined |

Total: approximately 260-320 GB per loaded version.

### Archive Files

The raw `.tar.gz` archive for each release is approximately 311 GB compressed. Archives are stored at:

```
/mnt/cglab.shared/Data/DBs/Uniprot/
```

Archives from previous releases can be removed once the corresponding database version is verified, unless version comparison studies require them.

---

## Versioning Policy

Multiple UniProt versions can coexist in the database simultaneously. All tables use a `version` column as part of their primary key. Queries always require a version argument to ensure reproducibility — the same query against version 2025_05 and 2026_01 may return different results as proteins are added, removed, or updated between releases.

It is recommended to record the UniProt version used for any analysis in lab notebooks and publications.

---

## Contact

For access issues, database errors, or questions about the pipeline, contact the person responsible for database maintenance in the lab.

For questions about UniProt data content, annotations, or proteome selection criteria, refer to the official UniProt documentation at https://www.uniprot.org/help.
