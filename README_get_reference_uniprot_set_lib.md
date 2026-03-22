# get_reference_uniprot_set_lib.py

## Overview

This script is the retrieval interface for the CGLab local UniProt Reference Proteome database.
It connects to a MySQL instance populated by `uniprot_sync.py` and allows querying of protein
sequences from any loaded UniProt reference proteome release.

It serves two purposes simultaneously. It can be run directly from the command line to export
a FASTA file, and it can be imported as a Python library into other scripts, allowing
programmatic access to sequences for downstream analyses such as multiple sequence alignment
and phylogenetic tree construction.

The database it connects to contains 133,566,280 canonical proteins from 34,230 reference
proteomes across Bacteria, Eukaryota, Viruses, and Archaea as of release 2026_01.

---

## Requirements

- Python 3.9 or higher
- mysql-connector-python
- python-dotenv
- BioPython (optional, required only for the `to_biopython()` method)

Install dependencies:

    pip install mysql-connector-python python-dotenv biopython

---

## Configuration

The script reads database credentials from a `.env` file located in the same directory as the
script. The file should contain:

    DB_HOST=192.168.1.100
    DB_USER=cglab_user
    DB_PASSWORD=your_password
    DB_NAME=uniprot_db_cglab

If no `.env` file is present, the script falls back to the following defaults:

    host:     localhost
    user:     cglab_user
    password: cglab2026
    database: uniprot_db_cglab

These defaults can also be overridden at runtime by passing explicit values to `get_db_config()`.

---

## Command-Line Usage

    python get_reference_uniprot_set_lib.py -version VERSION [filters]

### Required argument

    -version VERSION        UniProt release string, e.g. 2026_01

### Optional filters

    -taxonomy ID [ID ...]   One or more NCBI Taxonomy IDs
    --proteome-id ID        UniProt Proteome ID, e.g. UP000005640
    --go-id ID              Gene Ontology term, e.g. GO:0005634
    --pfam-id ID            Pfam domain accession, e.g. PF00870
    --output-dir DIR        Directory where the FASTA file will be written
    --list-versions         Print all versions present in the database and exit
    --list-proteomes        List all proteome IDs for the given version and exit

All filters are optional and combinable. When combined, they apply with AND logic, meaning
a protein must satisfy all specified filters to be included in the output.

### Examples

List all loaded versions in the database:

    python get_reference_uniprot_set_lib.py -version 2026_01 --list-versions

List all proteome IDs for a given version:

    python get_reference_uniprot_set_lib.py -version 2026_01 --list-proteomes

Retrieve the full human reference proteome:

    python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606

Retrieve proteins from multiple taxa:

    python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090 10116

Retrieve by Proteome ID:

    python get_reference_uniprot_set_lib.py -version 2026_01 --proteome-id UP000005640

Filter by GO term:

    python get_reference_uniprot_set_lib.py -version 2026_01 --go-id GO:0005634

Filter by Pfam domain:

    python get_reference_uniprot_set_lib.py -version 2026_01 --pfam-id PF00870

Combine filters (e.g. human kinases):

    python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 --pfam-id PF00069

Write output to a specific directory:

    python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 --output-dir ./fastas

### Output format

The command-line mode writes a FASTA file to the current directory or to `--output-dir` if
specified. The filename follows this format:

    uniprot_<identifier>_<version>.fasta

FASTA headers use the format:

    >{taxon_id}.{accession}

---

## Library Usage

The script is designed to be imported into other Python scripts for programmatic use.
This is the intended mode for tree-building pipelines and other downstream analyses in the lab.

### Import

    from get_reference_uniprot_set_lib import (
        UniProtRetriever,
        get_db_config,
        fetch_sequences,
        fetch_fasta_string,
    )

### Quick one-liner

    from get_reference_uniprot_set_lib import fetch_sequences

    records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
    for r in records:
        print(r["accession"], r["organism"])

### Context manager (recommended for multiple queries)

The context manager handles connection opening and closing automatically and is the recommended
pattern when issuing more than one query in the same script. It avoids the overhead of opening
a new database connection for every call.

    from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

    with UniProtRetriever(get_db_config()) as db:
        records = db.get_proteins(version="2026_01", taxon_ids=[9606])
        fasta_str   = db.to_fasta_string(records)
        seqrecords  = db.to_biopython(records)
        db.export_fasta(records, "2026_01", "human", output_dir="./fastas")

### Overriding database credentials at runtime

    config = get_db_config(host="10.0.0.5", user="myuser", password="mypass")
    with UniProtRetriever(config) as db:
        records = db.get_proteins(version="2026_01", taxon_ids=[6087])

### Piping directly into Bio.SeqIO without writing a file

    import io
    from Bio import SeqIO
    from get_reference_uniprot_set_lib import fetch_fasta_string

    fasta = fetch_fasta_string("2026_01", taxon_ids=[9606])
    seqs  = list(SeqIO.parse(io.StringIO(fasta), "fasta"))

---

## API Reference

### get_db_config(host, user, password, database)

Builds and returns the database configuration dictionary. Reads from environment variables
if set, otherwise uses lab defaults. Any argument passed explicitly takes priority over
environment variables and defaults.

Returns a dict suitable for passing directly into `UniProtRetriever(config)`.

---

### class UniProtRetriever

The main class. Manages the database connection and exposes all query and export methods.

#### connect() / close()

Open and close the database connection manually. These are not needed when using the
context manager, which calls them automatically.

#### list_available_versions()

Prints a summary table of all UniProt versions loaded in the database and returns a list
of dicts. Each dict contains the keys: `version`, `protein_count`, `taxon_count`,
and `proteome_count`.

    with UniProtRetriever(get_db_config()) as db:
        versions = db.list_available_versions()

#### get_proteome_ids(version)

Returns a list of all UniProt Proteome IDs (UP-prefixed strings) present in the database
for the given version.

    with UniProtRetriever(get_db_config()) as db:
        ids = db.get_proteome_ids("2026_01")

#### get_proteins(version, taxon_ids, proteome_id, go_id, pfam_id)

The main retrieval method. All filter arguments are optional. Returns a list of dicts.

Each dict in the returned list contains the following keys:

    accession    str    UniProt accession
    name         str    Protein name
    organism     str    Organism name string
    taxon_id     int    NCBI Taxonomy ID
    proteome_id  str    UniProt Proteome ID (UP-prefixed)
    sequence     str    Canonical amino acid sequence

When `go_id` is supplied, the query joins against the `protein_go` table.
When `pfam_id` is supplied, the query joins against the `protein_pfam` table.
Both can be combined with taxon or proteome filters simultaneously.

    records = db.get_proteins(
        version    = "2026_01",
        taxon_ids  = [9606, 10090],
        pfam_id    = "PF00069"
    )

#### to_fasta_string(records)

Converts a list of records (output of `get_proteins()`) to a FASTA-formatted string.
Does not write any file. Useful for piping directly into parsers.

    fasta = db.to_fasta_string(records)

#### to_biopython(records)

Converts records to a list of BioPython `SeqRecord` objects. Requires BioPython.
This is the recommended bridge to alignment tools (MUSCLE, MAFFT) and tree-building
libraries (ETE3, DendroPy, dendropy).

    seqrecords = db.to_biopython(records)

Each SeqRecord has:
- `id` set to `{taxon_id}.{accession}`
- `name` set to `{accession}`
- `description` set to `{protein_name} [{organism}] UP={proteome_id}`

#### export_fasta(records, version, identifier, output_dir, filename)

Writes records to a FASTA file on disk. The output directory is created automatically
if it does not exist. Returns the absolute path of the written file.

    path = db.export_fasta(records, "2026_01", "hydra", output_dir="./fastas")

If `filename` is not specified, the file is named:

    uniprot_<identifier>_<version>.fasta

---

### fetch_sequences(version, taxon_ids, proteome_id, go_id, pfam_id, db_config)

Module-level convenience function. Opens a connection, runs the query, closes the
connection, and returns the list of record dicts. Equivalent to using `UniProtRetriever`
as a context manager with a single `get_proteins()` call. Accepts an optional `db_config`
dict to override connection parameters.

    records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])

---

### fetch_fasta_string(version, taxon_ids, proteome_id, go_id, pfam_id, db_config)

Same as `fetch_sequences()` but returns a multi-sequence FASTA string instead of a list
of dicts. Accepts the same arguments.

---

## Notes for Downstream Use

- All records returned by `get_proteins()` are plain Python dicts and require no special
  dependencies to access or iterate.
- The `to_biopython()` method is the recommended bridge to alignment and tree-building tools.
- When issuing multiple queries in the same script, always prefer the context manager over
  repeated calls to `fetch_sequences()` to avoid opening and closing a new database connection
  for each query.
- The database stores only canonical sequences. Isoforms are not present, consistent with the
  UniProt reference proteome specification.
- For queries spanning very large taxon sets, consider whether filtering by Pfam domain or
  GO term first reduces the result set to a manageable size before loading sequences into memory.

---

## Known Issues

The docstring example in line 21 of the source file still references the old module name
`get_reference_uniprot_set_v5`. The correct module name for import is
`get_reference_uniprot_set_lib`. This will be corrected in the next version.
