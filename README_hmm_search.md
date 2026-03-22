# hmm_search.py

## Overview

This script runs HMMER hmmscan searches against the CGLab local UniProt Reference Proteome
database. It handles the complete workflow from sequence export to result storage in a single
pipeline, and is designed to operate at the scale of the full 2026_01 UniProt release
(133.5 million canonical proteins).

The pipeline proceeds in the following order:

1. Export protein sequences from MySQL to FASTA (chunked or streaming)
2. Download and prepare the Pfam-A HMM profile library
3. Run hmmscan on each FASTA chunk, in parallel
4. Parse the hmmscan domain table output
5. Load results into a dedicated MySQL table (`hmm_search_results`)

Searches can be scoped to specific taxa or proteomes rather than the full database, which
is the typical usage for comparative genomics analyses targeting selected species.

---

## Requirements

- Python 3.9 or higher
- HMMER 3.4 or higher (hmmscan and hmmpress must be on PATH)
- mysql-connector-python
- python-dotenv
- requests (for Pfam download)
- BioPython (optional)

Install Python dependencies:

    pip install mysql-connector-python python-dotenv requests

HMMER installation: http://hmmer.org/

---

## Configuration

The script reads database credentials from a `.env` file in the same directory as the script.
This is the same `.env` file used by `uniprot_sync.py` and `get_reference_uniprot_set_lib.py`.

    DB_HOST=192.168.1.100
    DB_USER=cglab_user
    DB_PASSWORD=your_password
    DB_NAME=uniprot_db_cglab

Default values if no `.env` is present:

    host:     localhost
    user:     cglab_user
    password: (empty)
    database: uniprot_db_cglab

---

## Usage

    python hmm_search.py --version VERSION --output-dir DIR [options]

### Required arguments

    --version VERSION       UniProt release string, e.g. 2026_01
    --output-dir DIR        Base output directory for all generated files

### Taxonomy and proteome filters

    --taxon-ids ID [ID ...]         Restrict to one or more NCBI Taxonomy IDs
    --proteome-ids ID [ID ...]      Restrict to one or more UniProt Proteome IDs

If neither filter is specified, all proteins in the given version are processed.

### Execution parameters

    --threads N             CPU threads per hmmscan process (default: 4)
    --parallel-jobs N       Number of concurrent hmmscan processes (default: 4)
    --chunk-size N          Proteins per FASTA chunk (default: 50000)
    --evalue FLOAT          Sequence-level E-value cutoff (default: 1e-5)
    --dom-evalue FLOAT      Domain-level E-value cutoff (default: 1e-3)

### Streaming and disk usage

    --stream                Stream chunks from MySQL one at a time (default).
                            Only one FASTA chunk exists on disk at a time (~200 MB).
                            Total disk usage is minimal (~1-2 GB for results only).
    --no-stream             Export ALL proteins to FASTA files first, then run hmmscan.
                            Requires approximately 50-60 GB of disk space for the full DB.
    --keep-fasta            Keep FASTA files after hmmscan completes (default: delete).

### Pipeline control

    --export-only           Export sequences to FASTA and stop. Does not run hmmscan.
    --import-only           Parse and import existing hmmscan result files. Does not
                            export sequences or run hmmscan. Useful when hmmscan was
                            run manually outside the pipeline.
    --skip-hmm-download     Skip the Pfam download step. Requires Pfam-A.hmm to already
                            exist in the output directory.

---

## Examples

Run a full scan against Hydra vulgaris proteins (taxon 6087) with 8 threads:

    python hmm_search.py \
        --version 2026_01 \
        --taxon-ids 6087 \
        --threads 8 \
        --output-dir ./hmm_results_hydra

Run a scan against multiple metazoan species with 4 parallel jobs, 4 threads each:

    python hmm_search.py \
        --version 2026_01 \
        --taxon-ids 6087 9606 10090 7227 6239 \
        --threads 4 \
        --parallel-jobs 4 \
        --output-dir ./hmm_results_metazoa

Export FASTA only (no search), for manual inspection or external use:

    python hmm_search.py \
        --version 2026_01 \
        --taxon-ids 9606 \
        --export-only \
        --output-dir ./fasta_export

Import previously generated hmmscan results without re-running the search:

    python hmm_search.py \
        --version 2026_01 \
        --import-only \
        --output-dir ./hmm_results_hydra

Run with a stricter E-value threshold and keep all FASTA files:

    python hmm_search.py \
        --version 2026_01 \
        --taxon-ids 9606 \
        --evalue 1e-10 \
        --dom-evalue 1e-5 \
        --keep-fasta \
        --output-dir ./hmm_strict

---

## Output Directory Structure

After a complete run, the output directory contains the following:

    <output-dir>/
        hmm_profiles/
            Pfam-A.hmm           Full Pfam HMM text file (~20,000 families, ~4 GB)
            Pfam-A.hmm.h3f       Binary index files created by hmmpress
            Pfam-A.hmm.h3i
            Pfam-A.hmm.h3m
            Pfam-A.hmm.h3p
            pfam_version.txt     Pfam release version string
        fasta_chunks/            (non-streaming mode only, deleted if --keep-fasta not set)
            chunk_00001.fasta
            chunk_00002.fasta
            ...
            manifest.txt         List of all chunk file paths with metadata
        tmp_fasta/               (streaming mode only, files deleted after each scan)
            stream_chunk_NNNNN.fasta
        hmmscan_results/
            chunk_00001_hmmscan.domtblout.txt
            chunk_00002_hmmscan.domtblout.txt
            ...

---

## FASTA Header Format

All FASTA files produced by the pipeline use this header format:

    >{taxon_id}.{accession}|{proteome_id}|{protein_name}|{organism}

Example:

    >6087.A0A0A0MP48|UP000007141|Wnt-3|Hydra vulgaris

This format preserves all metadata needed for downstream analyses and is parseable by
splitting on the pipe character.

---

## Database Table: hmm_search_results

Results are stored in the `hmm_search_results` table, which is created automatically on
the first import. The table schema is:

    id               BIGINT         Auto-increment primary key
    version          VARCHAR(10)    UniProt version (e.g. 2026_01)
    accession        VARCHAR(20)    UniProt protein accession
    taxon_id         INT            NCBI Taxonomy ID
    proteome_id      VARCHAR(20)    UniProt Proteome ID
    protein_name     VARCHAR(100)   Protein name (from FASTA header)
    hmm_name         VARCHAR(50)    HMM profile name (e.g. Wnt)
    hmm_accession    VARCHAR(20)    HMM accession (e.g. PF00858)
    hmm_type         ENUM           Pfam / TIGRFAM / SUPERFAMILY / other
    full_evalue      DOUBLE         Sequence-level E-value
    full_score       DOUBLE         Sequence-level bit score
    full_bias        DOUBLE         Sequence-level composition bias
    domain_number    INT            Index of this domain hit within the protein
    domain_count     INT            Total number of domain hits in the protein
    domain_evalue    DOUBLE         Domain-level E-value (conditional)
    domain_score     DOUBLE         Domain-level bit score
    domain_bias      DOUBLE         Domain-level composition bias
    hmm_from         INT            Start position in HMM (1-based)
    hmm_to           INT            End position in HMM (1-based)
    ali_from         INT            Start position in protein alignment (1-based)
    ali_to           INT            End position in protein alignment (1-based)
    env_from         INT            Start of domain envelope in protein
    env_to           INT            End of domain envelope in protein
    hmm_coverage     FLOAT          Fraction of the HMM length matched
    protein_coverage FLOAT          Fraction of the protein length matched
    posterior_prob   FLOAT          Mean posterior probability of aligned residues
    search_date      TIMESTAMP      When the row was inserted

Indexes are created on: `version`, `accession`, `hmm_name`, `hmm_accession`, `taxon_id`,
`proteome_id`, `full_evalue`, `domain_evalue`, and the composite `(version, accession)`.

### Useful queries against hmm_search_results

Retrieve all proteins in Hydra vulgaris with a hit to PF00858 (Wnt domain):

    SELECT accession, protein_name, full_evalue, domain_evalue,
           hmm_from, hmm_to, ali_from, ali_to, hmm_coverage
    FROM   hmm_search_results
    WHERE  version = '2026_01'
      AND  taxon_id = 6087
      AND  hmm_accession = 'PF00858'
    ORDER  BY full_evalue;

Count hits per HMM family across a set of taxa:

    SELECT hmm_accession, hmm_name,
           COUNT(DISTINCT accession) AS protein_count,
           COUNT(DISTINCT taxon_id)  AS taxon_count
    FROM   hmm_search_results
    WHERE  version  = '2026_01'
      AND  taxon_id IN (6087, 9606, 10090)
    GROUP  BY hmm_accession, hmm_name
    ORDER  BY protein_count DESC
    LIMIT  20;

Find proteins with multiple domain hits (domain repeats):

    SELECT accession, protein_name, taxon_id, hmm_accession,
           COUNT(*) AS domain_copies
    FROM   hmm_search_results
    WHERE  version = '2026_01'
      AND  taxon_id = 6087
    GROUP  BY accession, protein_name, taxon_id, hmm_accession
    HAVING domain_copies > 1
    ORDER  BY domain_copies DESC;

---

## Performance Notes

### Thread and parallelism settings

Each hmmscan call runs on `--threads` CPU threads. The pipeline launches up to
`--parallel-jobs` such calls simultaneously. The total CPU usage at peak is:

    --threads x --parallel-jobs

Set this to match your server's available CPUs. For example, on a 32-core machine:

    --threads 8 --parallel-jobs 4    (uses 32 cores)

The script auto-adjusts `--parallel-jobs` downward if the requested total exceeds
`os.cpu_count()`.

### Streaming mode vs full export

Streaming mode (the default) keeps disk usage minimal by writing one FASTA chunk,
scanning it, then deleting it before writing the next. The downside is that MySQL
must serve each chunk sequentially, so the pipeline cannot prefetch the next chunk
while hmmscan is running on the current one. In practice this is rarely a bottleneck
for scoped queries (single taxon or small species sets).

Full export mode (`--no-stream`) writes all FASTA chunks first, then runs hmmscan
across all of them in parallel. This can be faster for large multi-taxon runs but
requires approximately 50-60 GB of disk space for the full database, or proportionally
less for filtered subsets.

### Chunk size

The default chunk size of 50,000 proteins is a balance between memory usage and the
overhead of launching hmmscan processes. For taxon-scoped runs where the total protein
count is small (e.g. a single cnidarian proteome of ~20,000 proteins), a single chunk
will cover the entire dataset. For very large runs, increasing the chunk size to
100,000 reduces the number of hmmscan launches and can improve throughput.

### Re-running and appending results

If results for a given version already exist in the database, the importer appends to
them rather than replacing them. This is intentional and allows partial imports to be
resumed. If you need a clean re-import, truncate the relevant rows first:

    DELETE FROM hmm_search_results WHERE version = '2026_01' AND taxon_id = 6087;

---

## Notes for Tree-Building Workflows

The intended downstream use of this pipeline in the lab is to identify protein families
of interest across multiple reference proteomes and feed the results into alignment and
tree construction. The recommended workflow is:

1. Run `hmm_search.py` for the species set of interest.
2. Query `hmm_search_results` to retrieve accessions for the target HMM family.
3. Use `get_reference_uniprot_set_lib.py` (imported as a library) to retrieve full
   sequences for those accessions.
4. Pass the sequences to an alignment tool (MUSCLE, MAFFT) via BioPython.
5. Build the tree with IQ-TREE, FastTree, or RAxML.

A method `get_proteins_by_hmm_hit()` to support step 2 and 3 in a single library call
is planned for a future version of `get_reference_uniprot_set_lib.py`.
