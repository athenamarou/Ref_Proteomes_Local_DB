## README PYHMMER HMMSEARCH ACROSS LOCAL UNIPROT DB

A script that integrates **MySQL** queries and **pyHMMER** to perform **hmmsearch** across all the accessions for the local **Reference Proteome** Database built. 

**Requirements**

**📚 Libraries**

```text
- pyhmmer
- mysql-connector-python
- python-dotenv
- psutil
```

Install dependencies:

```bash
pip install pyhmmer mysql-connector-python python-dotenv psutil
```

**Database**

- MySQL
- Τables expected : proteins, sequences (joined on seq_id)
- User with INSERT and SELECT privilleges

**Environment Variables**

A **.env** file must be present in the same directory as the script, containing:

```
DB_HOST=localhost
DB_USER=cglab_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

**📁 Input files**

Pfam-A.hmm (all hmm profiles) placed at `{output_dir}/hmm_profiles/Pfam-A.hmm`

(The output_dir has to be created before running, since the checkpoint.txt file and the HMM profiles path are both derived from it.)

**💾 Hardware**

- **CPU**: The script is built to use 32 cores (PYHMMER_CPUS=32), but this constant can be adjusted if runninf in a machine with different capabilities.
- **RAM**: Recommended at least 16-32GB of RAM. The script is written to print a warning message if there is a memory consuming over 85%.

**💻 Architechture**

```
INPUT ARGUMENTS (--version, --chunk-size, --output-dir, etc.)
                           |
                           V
LOAD HMM PROFILES into memory (from HMM_DB_PATH)
                           |
                           V
INITIALIZE STREAMER & IMPORTER
  - Streamer: Prepares to read sequence chunks from MySQL.
  - Importer: Creates 'hmm_search_results' table in MySQL.
                            |
                            V
INITIALIZE THREAD QUEUES (Backpressure Control)
  - fetch_q (maxsize=2)
  - insert_q (maxsize=2)
                            |
                            V
START MULTI-THREADED PIPELINE (Concurrent Execution)
  |
  +--> [THREAD 1: FETCHER] (Daemon, Lightweight I/O)
  |      1. Streamer fetches chunk of sequences from MySQL DB.
  |      2. Push chunk into fetch_q (Blocks if queue is full).
  |      3. When DB is empty, push Sentinel (None) to signal end.
  |
  +--> [THREAD 2: SEARCHER] (Daemon, CPU Heavy)
  |      1. Pull sequence chunk from fetch_q.
  |      2. Run run_pyhmmer_hmmsearch (Parallelized internally).
  |      3. Delete sequence variables and force Garbage Collection.
  |      4. Push results and last accession ID to insert_q.
  |      5. If Sentinel received, forward it to insert_q and exit.
  |
  +--> [THREAD 3: INSERTER] (Non-Daemon, Safe I/O)
         1. Pull search results from insert_q.
         2. Importer bulk-inserts results into MySQL DB.
         3. Update and save accession checkpoint to disk.
         4. Delete result variables and force Garbage Collection.
         5. Calculate and print real-time metrics (RAM, seq/s).
         6. If Sentinel received, exit gracefully.
                            |
                            V
JOIN THREADS (Synchronized Teardown)
  - t1.join() -> t2.join() -> t3.join()
  - Main thread waits here. Mirrors data flow to ensure no partial DB commits.
                            |
                            V
PIPELINE COMPLETE (Print total processed summary)
```

✔️ The script makes sure that no double writes are made, because of **IF NOT EXISTS** for the **hmmsearch_results_table** and the **INSERT IGNORE** for the entries into the table. If thre script stops for whatever reason, the checkpoint file is used to retrieve the last accession and start again from there, without processing all the previous ones again.

**EXAMPLE USE**

```bash
# Example: Run hmmsearch across the whole DB, processing chunks of 25,000 sequences to limit RAM usage.
python3 -u /path/pyhmmer_hmmer_search.py --version 2026_01 \
	--output-dir/path/hmm_results \
	--chunk-size 25000 \
	2>&1 | tee -a /path/hmm_results/hmmsearch_2026_01.log
```

***`python3 -u`:** Forces unbuffered output. This ensures  real-time progress metrics (RAM, seq/s) print to the screen immediately instead of getting caught in Python's internal text buffer (if not prefered, can be removed).

***`2>&1 | tee -a`:** Captures both standard output (progress) and standard errors (crashes/warnings), displays them on the user's screen in real-time, and simultaneously appends them to a log file so    the data are not lost when the  terminal is closed.

‼️ Recommended to run the process inside a screen or tmux session to ensure it keeps running even if the computer is turned off. 

**OUTPUT EXAPLE**

```
############################################################
PyHMMER Search Pipeline — CGLab
UniProt Version: 2026_01
Output Dir:      /home/amarougka/uniprot_scripts/hmm_results
CPUs Allocated:  32
Started:         2026-03-29 21:12:47
############################################################
Loading HMM profiles from /home/amarougka/uniprot_scripts/hmm_results/hmm_profiles/Pfam-A.hmm into memory..
Resuming from checkpoint: A0A9Q9SG71
Results table 'hmm_search_results' ready.
```

```
Processed 28,650,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,675,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,700,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,725,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,750,000 proteins... (460 seq/s) | RAM: 14.3% (19.3/134.2 GB)
Processed 28,775,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,800,000 proteins... (460 seq/s) | RAM: 14.2% (19.1/134.2 GB)
Processed 28,825,000 proteins... (460 seq/s) | RAM: 14.4% (19.4/134.2 GB)
Processed 28,850,000 proteins... (460 seq/s) | RAM: 14.4% (19.4/134.2 GB)
Processed 28,875,000 proteins... (460 seq/s) | RAM: 14.5% (19.4/134.2 GB)
Processed 28,900,000 proteins... (460 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 28,925,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 28,950,000 proteins... (460 seq/s) | RAM: 14.2% (19.1/134.2 GB)
Processed 28,975,000 proteins... (460 seq/s) | RAM: 14.3% (19.1/134.2 GB)
Processed 29,000,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 29,025,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 29,050,000 proteins... (460 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 29,075,000 proteins... (460 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,100,000 proteins... (459 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,125,000 proteins... (459 seq/s) | RAM: 14.3% (19.3/134.2 GB)
Processed 29,150,000 proteins... (459 seq/s) | RAM: 14.2% (19.1/134.2 GB)
Processed 29,175,000 proteins... (459 seq/s) | RAM: 14.5% (19.4/134.2 GB)
Processed 29,200,000 proteins... (459 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 29,225,000 proteins... (459 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,250,000 proteins... (459 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,275,000 proteins... (459 seq/s) | RAM: 14.3% (19.2/134.2 GB)
Processed 29,300,000 proteins... (459 seq/s) | RAM: 14.3% (19.3/134.2 GB)
Processed 29,325,000 proteins... (459 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,350,000 proteins... (458 seq/s) | RAM: 14.4% (19.3/134.2 GB)
Processed 29,375,000 proteins... (458 seq/s) | RAM: 14.4% (19.4/134.2 GB)
```
