# get_reference_uniprot_set_lib.py

## Overview

A retrieval script-interface for the **Local UniProt Reference Proteome database**.
Connects to the local MySQL database and allows querying of protein sequences from any
loaded UniProt reference proteome release.

The script serves two purposes simultaneously:

- Run directly from the **command line** to export a FASTA file.
- Import as a **Python library** into other scripts for programmatic access to sequences,
  enabling downstream analyses such as multiple sequence alignment and phylogenetic tree construction.

The database contains **133,566,280 canonical proteins** from **34,230 reference proteomes**
across Bacteria, Eukaryota, Viruses, and Archaea as of release `2026_01`.
Domain annotations from a complete **Pfam-A hmmsearch, using pyHMMER** are available for retrieval,
enabling HMM-based protein set construction directly from the command line or library API.

---

## Requirements

**📚 Libraries**

```text
- mysql-connector-python
- python-dotenv
- biopython  (optional — required only for to_biopython())
```

Install dependencies:

```bash
pip install mysql-connector-python python-dotenv biopython
```

**Database**

- MySQL
- Tables expected: `proteins`, `sequences` (joined on `seq_id`), `protein_go`, `protein_pfam`, `hmm_search_results`
- User with `SELECT` privileges

**Environment Variables**

A `.env` file must be present in the same directory as the script:

```
DB_HOST=localhost
DB_USER=cglab_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

---

## Command-Line Usage

```bash
python get_reference_uniprot_set_lib.py -version VERSION [filters]
```

### Required argument

```
-version VERSION        UniProt release string, e.g. 2026_01
```

### Optional filters

```
-taxonomy ID [ID ...]   One or more NCBI Taxonomy IDs
--proteome-id ID        UniProt Proteome ID, e.g. UP000005640
--go-id ID              Gene Ontology term, e.g. GO:0005634
--pfam-id ID            Pfam domain accession via protein_pfam table, e.g. PF00870
--hmm-name NAME         Pfam HMM name or accession via local hmmsearch results,
                        e.g. Homeodomain or PF00046
--evalue FLOAT          Maximum full-sequence E-value for HMM hits (default: 1e-5).
                        Only applies when --hmm-name is used.
--output-dir DIR        Directory where the FASTA file will be written
--list-versions         Print all versions present in the database and exit
--list-proteomes        List all proteome IDs for the given version and exit
```

All filters are optional and combinable (AND logic).

### Examples

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

### Output format

```
uniprot_<identifier>_<version>.fasta
```

FASTA headers:

```
>{taxon_id}.{accession}
```

---

## Library Usage

### Import

```python
from get_reference_uniprot_set_lib import (
    UniProtRetriever,
    get_db_config,
    fetch_sequences,
    fetch_sequences_by_hmm_hit,
    fetch_fasta_string,
    fetch_fasta_string_by_hmm_hit,
)
```

### Standard retrieval (one-liner)

```python
from get_reference_uniprot_set_lib import fetch_sequences

records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
for r in records:
    print(r["accession"], r["organism"])
```

### HMM-based retrieval (one-liner)

The simplest way to get sequences for a domain family — directly from the local hmmsearch results:

```python
from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit

# All Homeodomain proteins across all taxa
records = fetch_sequences_by_hmm_hit(version="2026_01", hmm_query="Homeodomain")

# Human + mouse kinases with strict E-value
records = fetch_sequences_by_hmm_hit(
    version       = "2026_01",
    hmm_query     = "PF00069",
    evalue_cutoff = 1e-10,
    taxon_ids     = [9606, 10090],
)
```

### HMM-based retrieval → FASTA string (one-liner)

Pipe directly into BioPython without writing a file:

```python
import io
from Bio import SeqIO
from get_reference_uniprot_set_lib import fetch_fasta_string_by_hmm_hit

fasta = fetch_fasta_string_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
seqs  = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
```

### Context manager (recommended for multiple queries or tree building pipelines)

```python
from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

with UniProtRetriever(get_db_config()) as db:

    # Retrieve by HMM hit
    records = db.get_proteins_by_hmm_hit(
        version       = "2026_01",
        hmm_query     = "Homeodomain",
        evalue_cutoff = 1e-5,
        taxon_ids     = [9606, 10090],
    )

    # Convert to BioPython SeqRecord objects → pipe into MAFFT, ETE3, DendroPy etc.
    seqrecords = db.to_biopython(records)

    # Or export directly to FASTA
    db.export_fasta(records, "2026_01", "Homeodomain_human_mouse", output_dir="./fastas")
```

### Full tree-building example

```python
from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

with UniProtRetriever(get_db_config()) as db:

    # 1. Retrieve all bacterial proteins with a specific domain
    records = db.get_proteins_by_hmm_hit(
        version       = "2026_01",
        hmm_query     = "PF00069",   # kinase domain
        evalue_cutoff = 1e-5,
        taxon_ids     = [2],         # NCBI taxid for Bacteria
    )
    print(f"Retrieved {len(records):,} sequences")

    # 2. Convert to BioPython SeqRecord
    seqrecords = db.to_biopython(records)

    # 3. Write to FASTA for alignment
    db.export_fasta(records, "2026_01", "kinase_bacteria", output_dir="./fastas")

    # 4. Pass seqrecords directly into alignment/tree pipeline
    # e.g. with MAFFT via subprocess, or ETE3, or DendroPy
```

### Overriding database credentials at runtime

```python
config = get_db_config(host="10.0.0.5", user="myuser", password="mypass")
with UniProtRetriever(config) as db:
    records = db.get_proteins(version="2026_01", taxon_ids=[6087])
```

---

## 💻 Architecture

```
INPUT ARGUMENTS (-version, -taxonomy, --hmm-name, etc.)
                           |
                           V
CONNECT to MySQL (via UniProtRetriever or convenience function)
                           |
                           V
BUILD QUERY
  ├── Standard: proteins JOIN sequences
  │   Optional joins: protein_go, protein_pfam
  │   WHERE clauses from provided filters (AND logic)
  │
  └── HMM-based: hmm_search_results JOIN proteins JOIN sequences
      WHERE hmm_name or hmm_accession LIKE query
      AND full_evalue <= cutoff
      SELECT DISTINCT (deduplicates multi-domain proteins)
                           |
                           V
EXECUTE QUERY → fetchall() → list of dicts
                           |
                           V
OUTPUT (choose one or combine)
  ├── to_fasta_string()         → in-memory FASTA string
  ├── to_biopython()            → list of BioPython SeqRecord objects
  └── export_fasta()            → write .fasta file to disk
```

---

## API Reference

### Convenience functions (no class management needed)

| Function | Description |
|---|---|
| `fetch_sequences(version, ...)` | Standard filter retrieval → list of dicts |
| `fetch_sequences_by_hmm_hit(version, hmm_query, ...)` | HMM hit retrieval → list of dicts |
| `fetch_fasta_string(version, ...)` | Standard filter retrieval → FASTA string |
| `fetch_fasta_string_by_hmm_hit(version, hmm_query, ...)` | HMM hit retrieval → FASTA string |

### `get_db_config(host, user, password, database)`

Builds the database configuration dictionary from environment variables or explicit arguments.

### `class UniProtRetriever`

#### `connect()` / `close()`
Manual connection management. Not needed when using the context manager.

#### `list_available_versions()`
Prints and returns a summary of all UniProt versions in the database.

#### `get_proteome_ids(version)`
Returns a list of all UniProt Proteome IDs for the given version.

#### `get_proteins(version, taxon_ids, proteome_id, go_id, pfam_id)`
Standard retrieval with optional filters. Returns list of dicts.

#### `get_proteins_by_hmm_hit(version, hmm_query, evalue_cutoff, taxon_ids)`
Retrieves sequences for proteins with a hit to a given Pfam HMM profile in the local
`hmm_search_results` table. Results are deduplicated — one sequence per protein regardless
of domain copy count. Accepts either a Pfam name or accession as `hmm_query`.

#### `to_fasta_string(records)`
Converts records to a FASTA-formatted string. Does not write any file.

#### `to_biopython(records)`
Converts records to a list of BioPython `SeqRecord` objects. Requires BioPython.
Each SeqRecord has `id` = `{taxon_id}.{accession}` and `description` = `{protein_name} [{organism}] UP={proteome_id}`.

#### `export_fasta(records, version, identifier, output_dir, filename)`
Writes records to a FASTA file. Creates `output_dir` automatically if needed.
Returns the absolute path of the written file.

---

## Output Example

```
============================================================
UniProt Reference Set Retrieval
Version : 2026_01
HMM    : Homeodomain (E-value ≤ 1e-05)
============================================================

✓ Exported 142,301 sequences → /home/user/fastas/uniprot_hmm_Homeodomain_2026_01.fasta
Successfully retrieved 142,301 sequences.
```

---

## Notes

- All records returned are plain Python dicts — no special dependencies to access or iterate.
- `to_biopython()` is the recommended bridge to alignment tools (MUSCLE, MAFFT) and tree-building libraries (ETE3).
- When issuing multiple queries in the same script, always prefer the context manager to avoid connection overhead.
- The database stores only **canonical sequences**. Isoforms are not present.
- `get_proteins_by_hmm_hit()` queries `hmm_search_results` which was populated by a full Pfam-A hmmsearch using **gathering thresholds**. Results are therefore already filtered at the profile-specific trusted cutoff — `evalue_cutoff` provides an additional filter on top of that.
- For large taxon sets, consider filtering by HMM hit or GO term first to reduce the result set before loading sequences into memory.

--
