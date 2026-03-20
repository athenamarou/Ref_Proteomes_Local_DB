# HMM Search Pipeline for CGLab UniProt Reference Proteome Database

## Overview

This pipeline performs **HMMER `hmmscan`** searches using all Pfam HMM profiles against every protein stored in your local UniProt reference proteome MySQL database. It handles the entire workflow from sequence export through search execution to results storage — designed to work with the **133.5 million proteins** in the 2026_01 release.

---

## Table of Contents

1. [Background: What is HMM Search?](#1-background-what-is-hmm-search)
2. [Why Run Your Own HMM Search?](#2-why-run-your-own-hmm-search)
3. [Pipeline Architecture](#3-pipeline-architecture)
4. [Prerequisites](#4-prerequisites)
5. [Installation](#5-installation)
6. [Quick Start](#6-quick-start)
7. [Configuration](#7-configuration)
8. [Usage — Full Examples](#8-usage--full-examples)
9. [Scaling Considerations](#9-scaling-considerations)
10. [Output & Results](#10-output--results)
11. [MySQL Results Table Schema](#11-mysql-results-table-schema)
12. [Querying Results](#12-querying-results)
13. [Troubleshooting](#13-troubleshooting)
14. [How It Works — Detailed Explanation](#14-how-it-works--detailed-explanation)

---

## 1. Background: What is HMM Search?

### HMMs (Hidden Markov Models)

A **Hidden Markov Model** is a statistical model that captures the patterns of conservation and variation within a protein family. Unlike a simple sequence alignment or BLAST, an HMM encodes:

- **Match states**: Which positions are conserved (which amino acids are tolerated)
- **Insert states**: Where gaps/insertions occur
- **Delete states**: Where positions are skipped
- **Transition probabilities**: How likely each state is to follow another

### Pfam

[Pfam](https://pfam.xfam.org/) is the world's largest database of protein families, each represented by a curated HMM. As of 2026, Pfam contains **~20,000 protein families** (Pfam-A) covering the vast majority of known protein domains.

### HMMER

[HMMER](http://hmmer.org/) is the standard bioinformatics toolkit for sequence similarity searches using profile HMMs. The key commands are:

| Command | Purpose |
|---------|---------|
| `hmmscan` | Search a **single sequence** against a **database of HMMs** (what we use) |
| `hmmsearch` | Search a **single HMM** against a **database of sequences** |
| `hmmpress` | Index an HMM database for fast searching |
| `hmmfetch` | Extract specific HMMs from a database |

**Why `hmmscan` and not `hmmsearch`?**  
Your task is to find *which Pfam domains* exist in *every protein*. That means scanning each protein sequence against the entire Pfam library — exactly what `hmmscan` does. `hmmsearch` is for the opposite question: "which proteins contain this specific domain?"

---

## 2. Why Run Your Own HMM Search?

You might ask: *UniProt already annotates Pfam domains — why re-run this?*

Several important reasons:

1. **Different thresholds**: UniProt uses specific Pfam-defined thresholds (gathering thresholds / GA cutoffs). You can apply custom E-value or score cutoffs for your research question.

2. **Unannotated proteins**: Some reference proteome proteins (especially from less-studied organisms) may lack complete Pfam annotation in UniProt. Running your own scan catches these.

3. **Domain-level detail**: `hmmscan` provides per-domain coordinates, E-values, and scores for multi-domain proteins. The raw domain architecture is more detailed than what UniProt stores.

4. **Custom databases**: You can combine Pfam with other HMM libraries (TIGRFAMs, SUPERFAMILY, custom models) in a single pass.

5. **Reproducibility**: Having the full HMM search results in your own database means you can re-analyze with different cutoffs without re-running the computationally expensive search.

6. **Cross-version comparison**: Track how domain annotations change between UniProt releases using your own results table.

---

## 3. Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     CGLab HMM Search Pipeline                   │
└─────────────────────────────────────────────────────────────────┘

  Step 1: EXPORT                          Step 2: PREPARE HMMs
  ┌──────────────┐                         ┌──────────────────┐
  │   MySQL DB   │──── sequences ──────▶   │  Pfam FTP Server │
  │  133.5M prot │    as FASTA chunks      │  Pfam-A.hmm.gz   │
  └──────────────┘                         └────────┬─────────┘
                                                   │ download + decompress
                                                   ▼
  ┌────────────────────────┐              ┌──────────────────┐
  │   fasta_chunks/        │              │  Pfam-A.hmm      │
  │   chunk_00001.fasta    │              │  Pfam-A.hmm.h3f  │
  │   chunk_00002.fasta    │              │  Pfam-A.hmm.h3i  │
  │   ... (2,671 files)    │              │  Pfam-A.hmm.h3m  │
  │   chunk_02671.fasta    │              │  Pfam-A.hmm.h3p  │
  └───────────┬────────────┘              └────────┬─────────┘
              │                                    │
              │  Step 3: HMMSCAN                   │
              │  ┌─────────────────────────────────┐│
              ▼  │  hmmscan --domtblout            ││
  ┌───────────────────────────────────────────────┐││
  │            LOCAL MODE (multiprocessing)        │││
  │                                               │││
  │  Process 1: chunk_00001 → results_00001       │││
  │  Process 2: chunk_00002 → results_00002       │││
  │  Process 3: chunk_00003 → results_00003       │││
  │  Process 4: chunk_00004 → results_00004       │││
  │  ...                                          │││
  └───────────────────────────────────────────────┘││
                                                    ││
              OR                                    ││
              │                                     ││
              ▼                                     ││
  ┌───────────────────────────────────────────────┐││
  │            SLURM MODE (HPC cluster)           │││
  │                                               │││
  │  Job[1-100]:    chunk_00001 → results_00001   │││
  │  Job[101-200]:  chunk_00002 → results_00002   │││
  │  ...                                          │││
  └───────────────────────────────────────────────┘││
              │                                     ││
              ▼  Step 4: IMPORT                     ││
  ┌────────────────────────┐                       ││
  │  hmmscan_results/      │                       ││
  │  *_domtblout.txt files │                       ││
  └───────────┬────────────┘                       ││
              │ parse + bulk insert                 ││
              ▼                                     ││
  ┌────────────────────────┐                       ││
  │   MySQL DB             │◀──────────────────────┘│
  │   hmm_search_results   │                        │
  │   (new table)          │                        │
  └────────────────────────┘                        │
                                                   │
  ┌────────────────────────────────────────────────┘
  Output directory structure:
  
  hmm_results/
  ├── fasta_chunks/           # Step 1 output
  │   ├── chunk_00001.fasta
  │   ├── chunk_00002.fasta
  │   ├── ...
  │   └── manifest.txt
  ├── hmm_profiles/           # Step 2 output
  │   ├── Pfam-A.hmm
  │   ├── Pfam-A.hmm.h3f
  │   ├── Pfam-A.hmm.h3i
  │   ├── Pfam-A.hmm.h3m
  │   └── Pfam-A.hmm.h3p
  ├── hmmscan_results/        # Step 3 output
  │   ├── chunk_00001_hmmscan.domtblout.txt
  │   ├── chunk_00002_hmmscan.domtblout.txt
  │   └── ...
  ├── submit_all.sh           # (SLURM mode only)
  ├── hmmscan_job.sh          # (SLURM mode only)
  └── slurm_logs/             # (SLURM mode only)
```

---

## 4. Prerequisites

### Software

| Software | Version | Purpose | Install |
|----------|---------|---------|---------|
| **HMMER** | ≥ 3.3 | `hmmscan`, `hmmpress` | `conda install -c bioconda hmmer` or [hmmer.org](http://hmmer.org/) |
| **Python** | ≥ 3.9 | Pipeline script | System package manager |
| **MySQL Connector** | latest | Database access | `pip install mysql-connector-python` |
| **BioPython** | latest | Sequence handling | `pip install biopython` |
| **python-dotenv** | latest | .env file loading | `pip install python-dotenv` |
| **Requests** | latest | Pfam download | `pip install requests` |
| **SLURM** (optional) | — | HPC cluster scheduling | Pre-installed on HPC |

### Disk Space

| Resource | Estimated Size |
|----------|---------------|
| Pfam-A.hmm (pressed) | ~3 GB |
| hmmscan results (domtblout) | ~100–200 GB |
| **Streaming mode** | **~3–6 GB total** (only 1 chunk at a time) |
| **Batch mode** (full FASTA export) | **~50–60 GB** for FASTA + above |

> ⚡ **Streaming mode is the default.** It fetches chunks directly from MySQL, runs hmmscan, then deletes the temp FASTA. Only ~1 chunk (~200 MB) exists on disk at any time. Use `--no-stream` only if you need to keep the FASTA files for other tools.

### Database Access

The pipeline uses the same `.env` configuration as your existing scripts:

```env
DB_HOST=192.168.1.100
DB_USER=cglab_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

---

## 5. Installation

### Option A: Conda (Recommended)

```bash
conda create -n hmmer_pipeline python=3.11 -y
conda activate hmmer_pipeline
conda install -c bioconda hmmer -y
pip install mysql-connector-python biopython python-dotenv requests
```

### Option B: System Packages

```bash
# Ubuntu/Debian
sudo apt-get install hmmer python3-pip
pip3 install mysql-connector-python biopython python-dotenv requests

# CentOS/RHEL
sudo yum install hmmer python3-pip
pip3 install mysql-connector-python biopython python-dotenv requests

# macOS (Homebrew)
brew install hmmer
pip3 install mysql-connector-python biopython python-dotenv requests
```

### Verify Installation

```bash
hmmscan --version     # Should show HMMER 3.3.x or later
hmmpress -h           # Should show usage
python3 -c "import mysql.connector; print('MySQL OK')"
```

---

## 6. Quick Start

### Step 0: Place the pipeline files

```bash
# Copy hmm_search_pipeline.py and .env to your working directory
# The .env file should be the same one used by uniprot_sync_v7.py
ls -la .env hmm_search_pipeline.py
```

### Step 1: Export only (verify setup)

```bash
# Start with a small subset to test
python hmm_search_pipeline.py \
    --version 2026_01 \
    --taxon-ids 9606 \
    --chunk-size 10000 \
    --export-only \
    --output-dir ./hmm_test
```

This exports only human proteins (~20K) as a sanity check.

### Step 2: Run a test search

```bash
python hmm_search_pipeline.py \
    --version 2026_01 \
    --taxon-ids 9606 \
    --chunk-size 10000 \
    --threads 4 \
    --output-dir ./hmm_test
```

> **Memory & disk usage**: In streaming mode (default), this only uses ~50–100 MB RAM and writes ~200 MB of temp FASTA at a time. No large persistent files.

### Step 3: Full production run

```bash
python hmm_search_pipeline.py \
    --version 2026_01 \
    --threads 4 \
    --parallel-jobs 8 \
    --chunk-size 50000 \
    --output-dir /mnt/cglab.shared/Data/DBs/hmm_results_2026_01
```

---

## 7. Configuration

### Command-Line Parameters

#### Required

| Parameter | Description |
|-----------|-------------|
| `--version` | UniProt DB version (e.g., `2026_01`) |
| `--output-dir` | Base output directory |

#### Filters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--taxon-ids` | All taxa | One or more NCBI Taxonomy IDs |
| `--proteome-ids` | All proteomes | Specific proteome IDs |

#### Execution

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mode` | `local` | `local` or `slurm` |
| `--threads` | 4 | CPU threads per hmmscan call |
| `--parallel-jobs` | 4 | Concurrent hmmscan processes (local mode) |
| `--chunk-size` | 50,000 | Proteins per FASTA chunk |
| `--evalue` | 1e-5 | Sequence-level E-value cutoff |
| `--dom-evalue` | 1e-3 | Domain-level E-value cutoff |

#### SLURM-Specific

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--slurm-ncpus` | 8 | CPUs per SLURM task |
| `--slurm-partition` | compute | SLURM partition name |
| `--slurm-time` | 04:00:00 | Time limit per task |

#### Partial Execution

| Parameter | Description |
|-----------|-------------|
| `--export-only` | Only export FASTA, skip search (implies `--no-stream`) |
| `--import-only` | Only import results, skip everything else |
| `--skip-hmm-download` | Use already-downloaded HMM profiles |
| `--stream` | (Default) Stream chunks from MySQL — minimal disk (~1-2 GB) |
| `--no-stream` | Export ALL FASTA to disk first (~50-60 GB), then search |
| `--keep-fasta` | Don't delete FASTA files after hmmscan completes |

---

## 8. Usage — Full Examples

### Example 1: Full database scan (local mode)

Run hmmscan against all 133.5M proteins on a machine with 32 cores:

```bash
python hmm_search_pipeline.py \
    --version 2026_01 \
    --output-dir /mnt/cglab.shared/Data/DBs/hmm_full_2026_01 \
    --threads 4 \
    --parallel-jobs 8 \
    --chunk-size 50000
```

**Estimated time**: ~3–5 days on a 32-core server.

### Example 2: Eukaryota only

Filter to Eukaryota (you'll need to look up the major taxon IDs or use proteome-level filtering):

```bash
# Common model organism TaxIDs
python hmm_search_pipeline.py \
    --version 2026_01 \
    --taxon-ids 9606 10090 10116 7227 6239 3702 \
    --output-dir ./hmm_eukaryota \
    --threads 8 \
    --parallel-jobs 4
```

### Example 3: Single proteome

```bash
python hmm_search_pipeline.py \
    --version 2026_01 \
    --proteome-ids UP000005640 \
    --output-dir ./hmm_human \
    --threads 8
```

### Example 4: SLURM (HPC cluster)

```bash
# Step 1: Export + prepare (fast)
python hmm_search_pipeline.py \
    --version 2026_01 \
    --export-only \
    --output-dir /scratch/$USER/hmm_2026_01

# Step 2: Prepare HMMs (fast)
python hmm_search_pipeline.py \
    --version 2026_01 \
    --skip-hmm-download \
    --mode slurm \
    --slurm-ncpus 8 \
    --slurm-partition compute \
    --slurm-time 02:00:00 \
    --output-dir /scratch/$USER/hmm_2026_01

# Step 3: Submit jobs
cd /scratch/$USER/hmm_2026_01
bash submit_all.sh

# Step 4: Monitor
squeue -u $USER
watch -n 60 'ls hmmscan_results/*_domtblout.txt 2>/dev/null | wc -l'

# Step 5: Import results after all jobs complete
python hmm_search_pipeline.py \
    --version 2026_01 \
    --import-only \
    --output-dir /scratch/$USER/hmm_2026_01
```

### Example 5: Tighter E-value cutoff

For high-confidence domain assignments only:

```bash
python hmm_search_pipeline.py \
    --version 2026_01 \
    --taxon-ids 9606 \
    --evalue 1e-10 \
    --dom-evalue 1e-5 \
    --output-dir ./hmm_strict
```

### Example 6: Export, then manually run hmmscan

```bash
# Export FASTA
python hmm_search_pipeline.py \
    --version 2026_01 \
    --export-only \
    --output-dir ./my_hmm

# Manually run hmmscan (e.g., with custom flags)
hmmscan --cpu 16 --domtblout results.txt /path/to/Pfam-A.hmm ./my_hmm/fasta_chunks/chunk_00001.fasta

# Import results
python hmm_search_pipeline.py \
    --version 2026_01 \
    --import-only \
    --output-dir ./my_hmm
```

---

## 9. Scaling Considerations

### Streaming vs Batch Mode

The pipeline runs in **streaming mode by default**. This is the right choice for almost all use cases:

| Aspect | Streaming (default) | Batch (--no-stream) |
|--------|-------------------|---------------------|
| **Disk usage** | ~3–6 GB total | ~50–60 GB for FASTA |
| **Memory usage** | ~50–100 MB | ~50–100 MB |
| **FASTA files kept?** | No (auto-deleted) | Yes |
| **Resume after crash** | Restart from beginning | Can restart from failed chunks |
| **SLURM support** | Not applicable | Yes (needs files on shared storage) |
| **When to use** | Local runs, limited disk | HPC clusters, need FASTA for other tools |

### How Streaming Mode Works

1. MySQL query fetches 50,000 proteins using `LIMIT/OFFSET`
2. They're written to a single temp FASTA file (~200 MB)
3. `hmmscan` processes that file
4. The temp FASTA is **immediately deleted**
5. The next 50,000 proteins are fetched

At no point does the pipeline hold more than one chunk in memory or on disk.

### Full Database (133.5M proteins)

Running hmmscan against all 133.5M proteins with ~20,000 Pfam HMMs is a **massive** computational task. Here's what to expect:

| Resource | Single Machine (32 cores) | HPC (100 nodes × 16 cores) |
|----------|--------------------------|---------------------------|
| **Export** | ~2–4 hours | ~2–4 hours |
| **HMM download + press** | ~30 min | ~30 min |
| **hmmscan** | ~3–7 days | ~6–12 hours |
| **Import results** | ~2–4 hours | ~2–4 hours |
| **Disk I/O** | ~50 GB FASTA + ~200 GB results | Same |

### Strategy Recommendations

1. **Start with a test subset** — always verify with `--taxon-ids 9606` before committing to a full run.

2. **Use HPC if available** — SLURM mode is strongly recommended for the full database. A single machine will take days.

3. **Adjust chunk size** — smaller chunks (10K–25K) allow finer-grained parallelism and better load balancing on HPC. Larger chunks (100K) reduce overhead on local machines with many cores.

4. **Monitor disk space** — the FASTA export alone is ~50 GB. Ensure you have 250+ GB free.

5. **Use `nohup` or `screen`/`tmux`** for long-running local jobs:
   ```bash
   nohup python hmm_search_pipeline.py --version 2026_01 --threads 4 --output-dir ./hmm_full > pipeline.log 2>&1 &
   ```

6. **Resume capability** — if the pipeline crashes during hmmscan, you can identify which chunks succeeded (by looking at `hmmscan_results/`) and re-run only those. The import step appends, so it's safe to re-run.

### Resource Requirements per hmmscan Process

| Metric | Value |
|--------|-------|
| RAM per process | ~2–4 GB |
| RAM with `--cpu 8` | ~8–16 GB |
| Disk read (FASTA) | ~50–200 MB/s |
| Disk write (results) | ~100–300 MB/s |

---

## 10. Output & Results

### hmmscan Domain Table Format (--domtblout)

The pipeline uses `hmmscan --domtblout` which produces a tab-delimited output with one line per domain hit. Each line contains:

```
target_name  target_acc  tlen  query_name  query_acc  qlen
e-value  score  bias  #  of  HMM  from  to  hmm_cov
#  of  dom  from  to  dom_cov  env_from  env_to  env_cov  acc
```

**Key columns explained:**

| Column | Meaning |
|--------|---------|
| `target_name` | Pfam family name (e.g., `P53`) |
| `target_acc` | Pfam accession (e.g., `PF00870`) |
| `tlen` | Total length of the HMM model |
| `query_name` | Your protein (e.g., `9606.P04637\|UP000005640\|P53_HUMAN\|Homo sapiens`) |
| `e-value` | Sequence-level E-value (how significant is this HMM matching this protein) |
| `score` | Bit score |
| `domain_evalue` | Domain-level E-value (significance of this specific domain hit) |
| `hmm_from` / `hmm_to` | Coordinates within the HMM model |
| `ali_from` / `ali_to` | Coordinates within the protein sequence |
| `hmm_coverage` | Fraction of the HMM model matched |
| `protein_coverage` | Fraction of the protein covered |

### Understanding E-values

- **E-value < 1e-5**: Likely a true domain hit (default threshold)
- **E-value < 1e-10**: High-confidence domain assignment
- **E-value < 1e-50**: Nearly certain domain assignment
- **E-value > 1e-3**: Likely a false positive

The pipeline uses two E-value cutoffs:
- `--evalue` (default 1e-5): Controls which proteins are reported at all
- `--dom-evalue` (default 1e-3): Controls which individual domain hits within a protein are reported

---

## 11. MySQL Results Table Schema

The pipeline creates a new table called `hmm_search_results` in your existing `uniprot_db_cglab` database:

```sql
CREATE TABLE IF NOT EXISTS hmm_search_results (
    id BIGINT AUTO_INCREMENT PRIMARY KEY,
    version VARCHAR(10) NOT NULL,          -- UniProt version (e.g., '2026_01')

    -- Protein identification
    accession VARCHAR(20) NOT NULL,        -- UniProt accession (e.g., 'P04637')
    taxon_id INT,                          -- NCBI Taxonomy ID
    proteome_id VARCHAR(20),               -- Proteome ID (e.g., 'UP000005640')
    protein_name VARCHAR(100),             -- Protein name (e.g., 'P53_HUMAN')

    -- HMM profile information
    hmm_name VARCHAR(50) NOT NULL,         -- Pfam family name (e.g., 'P53')
    hmm_accession VARCHAR(20),             -- Pfam accession (e.g., 'PF00870')
    hmm_type ENUM('Pfam', 'TIGRFAM', 'SUPERFAMILY', 'other'),

    -- Full sequence match statistics
    full_evalue DOUBLE,                    -- Best E-value for this protein+HMM pair
    full_score DOUBLE,                     -- Corresponding bit score
    full_bias DOUBLE,                      -- Composition bias correction

    -- Individual domain hit
    domain_number INT,                     -- Which domain hit (1st, 2nd, etc.)
    domain_count INT,                      -- Total domains of this type in the protein
    domain_evalue DOUBLE,                  -- Domain-level E-value
    domain_score DOUBLE,                   -- Domain-level bit score
    domain_bias DOUBLE,

    -- Alignment coordinates (1-based, inclusive)
    hmm_from INT,                          -- Start position in HMM
    hmm_to INT,                            -- End position in HMM
    ali_from INT,                          -- Start position in protein
    ali_to INT,                            -- End position in protein
    env_from INT,                          -- Envelope start (wider than alignment)
    env_to INT,                            -- Envelope end

    -- Coverage metrics
    hmm_coverage FLOAT,                    -- Fraction of HMM model matched
    protein_coverage FLOAT,                -- Fraction of protein sequence matched
    posterior_prob FLOAT,                  -- Posterior probability of the domain

    -- Metadata
    search_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Indexes for fast querying
    INDEX idx_version (version),
    INDEX idx_accession (accession),
    INDEX idx_hmm (hmm_name),
    INDEX idx_hmm_acc (hmm_accession),
    INDEX idx_taxon (taxon_id),
    INDEX idx_proteome (proteome_id),
    INDEX idx_evalue (full_evalue),
    INDEX idx_domain_eval (domain_evalue),
    INDEX idx_version_acc (version, accession)
);
```

This table is **separate** from your existing `protein_pfam` table, which stores UniProt's curated annotations. Your HMM search results are independent and can be compared against UniProt's annotations.

---

## 12. Querying Results

Here are some useful SQL queries to explore your results:

### Count of total hits per version
```sql
SELECT version,
       COUNT(*) AS total_hits,
       COUNT(DISTINCT accession) AS unique_proteins,
       COUNT(DISTINCT hmm_accession) AS unique_domains
FROM   hmm_search_results
GROUP BY version;
```

### Top 20 most common Pfam domains
```sql
SELECT hmm_accession, hmm_name, COUNT(*) AS hit_count,
       COUNT(DISTINCT accession) AS protein_count,
       COUNT(DISTINCT taxon_id) AS taxon_count
FROM   hmm_search_results
WHERE  version = '2026_01'
GROUP BY hmm_accession, hmm_name
ORDER BY hit_count DESC
LIMIT 20;
```

### Domain architecture of a specific protein
```sql
SELECT hmm_name, hmm_accession,
       domain_number, domain_count,
       ali_from, ali_to,
       domain_evalue, domain_score,
       hmm_coverage, protein_coverage
FROM   hmm_search_results
WHERE  version = '2026_01'
  AND  accession = 'P04637'
ORDER BY ali_from;
```

### Find proteins with a specific domain in specific organisms
```sql
SELECT r.accession, r.protein_name, r.organism,
       r.full_evalue, r.ali_from, r.ali_to
FROM   hmm_search_results r
WHERE  r.version = '2026_01'
  AND  r.hmm_accession = 'PF00069'   -- Protein kinase domain
  AND  r.taxon_id IN (9606, 10090)   -- Human + Mouse
ORDER BY r.full_evalue;
```

### Compare your results vs UniProt's annotations
```sql
-- Domains found by your search but NOT in UniProt annotation
SELECT DISTINCT r.hmm_accession, r.hmm_name,
       r.accession, r.taxon_id
FROM   hmm_search_results r
LEFT JOIN protein_pfam u
  ON r.accession = u.accession
 AND r.version = u.version
 AND r.hmm_accession = u.pfam_id
WHERE  r.version = '2026_01'
  AND  u.pfam_id IS NULL;
```

### Domain statistics per taxonomic group
```sql
SELECT 
    CASE 
        WHEN taxon_id <= 2      THEN 'Bacteria'
        WHEN taxon_id BETWEEN 2759 AND 2157 THEN 'Eukaryota'
        WHEN taxon_id BETWEEN 2157 AND 6191 THEN 'Archaea'
        ELSE 'Other'
    END AS kingdom,
    COUNT(DISTINCT accession) AS proteins,
    COUNT(DISTINCT hmm_accession) AS domains,
    AVG(protein_coverage) AS avg_coverage
FROM   hmm_search_results
WHERE  version = '2026_01'
GROUP BY kingdom;
```

---

## 13. Troubleshooting

### "hmmscan: not found"
HMMER is not installed or not in your PATH:
```bash
# Check if installed
which hmmscan

# Install with conda
conda install -c bioconda hmmer

# Or download from http://hmmer.org/ and add to PATH
export PATH=/path/to/hmmer/bin:$PATH
```

### MySQL connection refused
- Verify the MySQL server is running: `mysql -u cglab_user -p -h 192.168.1.100`
- Check your `.env` file is in the same directory as the script
- If running on a different machine, ensure the server allows remote connections

### Out of disk space
The full pipeline needs ~200 GB. Check available space:
```bash
df -h /mnt/cglab.shared/
```

Reduce space by:
- Using `--taxon-ids` to filter
- Deleting FASTA chunks after search: `rm -rf hmm_results/fasta_chunks/`

### hmmscan is very slow
- Increase `--threads` (more CPU per call)
- Increase `--parallel-jobs` (more concurrent calls)
- Use `--mode slurm` on HPC
- Consider reducing `--chunk-size` for better load balancing

### Process killed (OOM)
hmmscan uses ~2–4 GB RAM per process. With `--parallel-jobs 8`, you need ~16–32 GB free RAM. Reduce parallel jobs if you hit memory limits.

### SLURM jobs failing
- Check logs: `cat hmm_results/slurm_logs/hmmscan_*.err`
- Common issue: HMMER not available in the SLURM node's PATH. Add `module load hmmer` to the job script
- Time limit too short: increase `--slurm-time`

---

## 14. How It Works — Detailed Explanation

### Step 1: Sequence Export (MySQL → FASTA)

The `SequenceExporter` class connects to your MySQL database and streams protein sequences using an unbuffered cursor (SSCursor). This is critical because loading 133.5 million rows into memory would require ~50+ GB of RAM.

**What happens:**
1. A `COUNT(*)` query determines the total number of matching proteins
2. The main query streams rows one at a time using `dictionary=True, buffered=True` cursor
3. Each protein is written as a FASTA record with a metadata-rich header:
   ```
   >9606.P04637|UP000005640|P53_HUMAN|Homo sapiens
   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
   ...
   ```
4. Every `chunk_size` (default 50,000) proteins, the current file is closed and a new one opens
5. A `manifest.txt` is written listing all chunk files

**Why chunk?** HMMER is single-threaded per invocation. By splitting into chunks, we can run multiple `hmmscan` processes in parallel, fully utilizing multi-core hardware.

### Step 2: HMM Profile Preparation

The `HMMProfileManager` handles downloading and indexing the Pfam HMM library:

1. **Download**: Fetches `Pfam-A.hmm.gz` (~1.5 GB compressed, ~3 GB uncompressed) from the EBI FTP server
2. **Decompress**: Gunzips to `Pfam-A.hmm`
3. **Press**: Runs `hmmpress` which creates 4 binary index files (`.h3f`, `.h3i`, `.h3m`, `.h3p`). These enable HMMER to search the HMM library efficiently using a pre-computed profile index (similar to how BLAST uses a sequence database index)

The pressed files reduce search time by ~10–100× compared to searching uncompressed HMM files.

### Step 3: hmmscan Execution

The pipeline runs `hmmscan` with these key flags:

```
hmmscan --cpu 4 \           # Multi-threaded within each call
        --domtblout out.txt \ # Domain table output (one line per hit)
        --noali \              # Don't print alignments (saves I/O time)
        -E 1e-5 \              # Sequence-level E-value cutoff
        --domE 1e-3 \          # Domain-level E-value cutoff
        Pfam-A.hmm \           # The HMM database
        chunk.fasta            # Input sequences
```

**`--domtblout` vs standard output**: The domain table format is much more compact than the default human-readable output and contains all the information needed for downstream analysis. Each line is a single domain hit, making parsing trivial.

**`--noali`**: Suppresses alignment printing. Since we're interested in *which* domains are present (not the specific alignments), this saves significant disk I/O.

**Two E-value cutoffs**: `hmmscan` applies E-value filtering at two levels:
- **Sequence level** (`-E`): "Is this protein likely to contain this domain family?" — controls which protein–HMM pairs appear in the output at all
- **Domain level** (`--domE`): "Is this specific domain instance real?" — filters out spurious domain hits within a protein that passes the sequence-level filter

**Local mode**: Uses Python's `ProcessPoolExecutor` to run multiple `hmmscan` processes simultaneously. The number of parallel processes is controlled by `--parallel-jobs`, and each process uses `--threads` CPU cores.

**SLURM mode**: Generates a SLURM array job script where each array task processes one FASTA chunk. This distributes work across cluster nodes and handles automatic retries on node failures.

### Step 4: Results Import (domtblout → MySQL)

The `HMMResultsImporter`:

1. **Creates the table**: `hmm_search_results` with indexes on all commonly queried columns
2. **Parses each file**: Reads domtblout format, extracting protein metadata from the FASTA header (taxon ID, proteome ID, protein name) and HMM metadata from the Pfam columns
3. **Bulk inserts**: Accumulates rows in batches of 100,000 and uses `executemany` for efficient MySQL insertion
4. **Computes coverage**: Calculates both HMM coverage and protein coverage for each domain hit
5. **Detects HMM type**: Classifies each HMM as Pfam, TIGRFAM, or SUPERFAMILY based on the accession prefix

**FASTA header parsing**: The pipeline encodes all protein metadata in the FASTA header during export, so the import step doesn't need to query MySQL for each protein — it extracts everything from the header string. This makes the import step self-contained and fast.

---

## File Manifest

```
hmm_search_pipeline/
├── README.md                    # This file
├── hmm_search_pipeline.py       # Main pipeline script
└── .env                         # Database credentials (reuse from sync scripts)
```

---

## Relationship to Existing Scripts

```
Your existing workflow:
  uniprot_sync_v7.py  →  MySQL DB  →  get_reference_uniprot_set_lib.py  →  FASTA

This pipeline:
  MySQL DB  →  hmm_search_pipeline.py  →  MySQL DB (hmm_search_results table)
```

The pipeline **reads from** the same database built by `uniprot_sync_v7.py` and **writes to** a new table in that same database. It can optionally use `get_reference_uniprot_set_lib.py`'s retrieval functions for programmatic access to results.
