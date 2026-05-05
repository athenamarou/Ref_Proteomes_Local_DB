"""
UniProt Reference Set Retrieval Library
========================================
Retrieves protein sequences from the local UniProt database.

Can be used in two ways:

1. As a command-line tool (for quick retrieval without writing custom scripts):
   python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090

2. As an importable library in your scripts:

   # --- Quickstart (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences

   records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
   for r in records:
       print(r["accession"], r["organism"])

   # --- HMM hit retrieval (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit

   records = fetch_sequences_by_hmm_hit(
       version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
   )

   # --- Context manager (recommended for multiple queries) ---
   from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

   with UniProtRetriever(get_db_config()) as db:
       records    = db.get_proteins(version="2026_01", taxon_ids=[9606, 10090])
       hmm_recs   = db.get_proteins_by_hmm_hit(version="2026_01", hmm_query="Homeodomain")
       fasta_str  = db.to_fasta_string(records)       # in-memory FASTA string
       seqrecords = db.to_biopython(records)           # list of BioPython SeqRecord
       db.export_fasta(records, "2026_01", "mammals", output_dir="./fastas")
"""

import mysql.connector
import argparse
import os
import sys
from dotenv import load_dotenv

load_dotenv()

# Public API — what users get when they do `from ... import *`
__all__ = [
    "UniProtRetriever",
    "get_db_config",
    "fetch_sequences",
    "fetch_sequences_by_hmm_hit",
    "fetch_fasta_string",
    "fetch_fasta_string_by_hmm_hit",
    "fetch_domains_by_accession",
]


# ---------------------------------------------------------------------------
# Configuration helper
# ---------------------------------------------------------------------------

def get_db_config(host=None, user=None, password=None, database=None):
    """
    Build the database config dict, falling back to environment variables
    and then to the lab defaults.

    users can override any value when calling this function, or by
    setting environment variables in a .env file:
        DB_HOST, DB_USER, DB_PASSWORD, DB_NAME

    Parameters
    ----------
    host, user, password, database : str, optional
        Explicit overrides — take priority over everything else.

    Returns
    -------
    dict  ready to pass into UniProtRetriever(config)
    """
    return {
        "host":     host     or os.getenv("DB_HOST",     "localhost"),
        "user":     user     or os.getenv("DB_USER",     "user"),
        "password": password or os.getenv("DB_PASSWORD", "your_password"),
        "database": database or os.getenv("DB_NAME",     "uniprot_db"),
    }


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

class UniProtRetriever:
    """
    Retrieves UniProt reference sets from the local CGLab database.

    Supports filtering by TaxID, Proteome ID, Pfam domain, GO term, and HMM search
    results from the local Pfam-A hmmsearch conducted by pyhmmer hmmsearch.

    Results can be returned as raw dicts, a FASTA string, BioPython
    SeqRecord objects, or written directly to a .fasta file.

    Recommended usage — context manager (handles connect/close):

        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins(version="2026_01", taxon_ids=[9606])
            seqs = db.to_biopython(records)

        # HMM-based retrieval for tree building:
        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins_by_hmm_hit(
                version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
            )
            seqs = db.to_biopython(records)
    """

    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict
            mysql.connector connection parameters.
            Use get_db_config() to build this from environment variables.
        """
        self.config = config
        self.conn   = None
        self.cursor = None

    # ------------------------------------------------------------------
    # Connection management
    # ------------------------------------------------------------------

    def connect(self):
        """
        Open the database connection.

        """
        self.conn   = mysql.connector.connect(**self.config)
        self.cursor = self.conn.cursor(dictionary=True)

    def close(self):
        """Close the database connection (safe to call even if not connected)."""
        if self.conn and self.conn.is_connected():
            self.cursor.close()
            self.conn.close()

    # Context-manager support — in order to be able to use `with UniProtRetriever(...) as db:`
    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    # ------------------------------------------------------------------
    # Utility queries
    # ------------------------------------------------------------------

    def list_available_versions(self):
        """
        Print a summary table of all versions present in the database.

        Returns
        -------
        list[dict]  one dict per version with keys:
                    version, protein_count, taxon_count, proteome_count
        """
        query = """
            SELECT version,
                   COUNT(*)                  AS protein_count,
                   COUNT(DISTINCT taxon_id)  AS taxon_count,
                   COUNT(DISTINCT proteome_id) AS proteome_count
            FROM   proteins
            GROUP  BY version
            ORDER  BY version DESC
        """
        self.cursor.execute(query)
        versions = self.cursor.fetchall()
        if versions:
            print("\nAvailable versions in database:")
            print("-" * 70)
            for v in versions:
                print(
                    f"Version: {v['version']} | "
                    f"Proteins: {v['protein_count']:,} | "
                    f"Taxa: {v['taxon_count']:,} | "
                    f"Proteomes: {v['proteome_count']:,}"
                )
            print("-" * 70)
        else:
            print("\nNo versions found in database.")
        return versions

    def get_proteome_ids(self, version):
        """
        Return all unique Proteome IDs present for a given version.

        Parameters
        ----------
        version : str  e.g. "2026_01"

        Returns
        -------
        list[str]
        """
        query = "SELECT DISTINCT proteome_id FROM proteins WHERE version = %s"
        self.cursor.execute(query, (version,))
        return [row["proteome_id"] for row in self.cursor.fetchall()]

    # ------------------------------------------------------------------
    # Main retrieval
    # ------------------------------------------------------------------

    def get_proteins(
        self,
        version,
        taxon_ids=None,
        proteome_id=None,
        go_id=None,
        pfam_id=None,
    ):
        """
        Retrieve protein records from the database.

        All filters are optional and combinable.

        Parameters
        ----------
        version : str
            UniProt release version, e.g. "2026_01".
        taxon_ids : int or list[int], optional
            One or more NCBI Taxonomy IDs.
        proteome_id : str, optional
            UniProt Proteome ID, e.g. "UP000005640".
        go_id : str, optional
            Gene Ontology term, e.g. "GO:0005634".
        pfam_id : str, optional
            Pfam domain accession, e.g. "PF00870".

        Returns
        -------
        list[dict]
            Each dict has keys:
            accession, name, organism, taxon_id, proteome_id, sequence
        """
        query = """
            SELECT p.accession, p.name, p.organism,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM   proteins  p
            JOIN   sequences s ON p.seq_id = s.seq_id
        """
        if go_id:
            query += (
                " JOIN protein_go pg"
                " ON p.accession = pg.accession AND p.version = pg.version"
            )
        if pfam_id:
            query += (
                " JOIN protein_pfam pp"
                " ON p.accession = pp.accession AND p.version = pp.version"
            )

        where_clauses = ["p.version = %s"]
        params        = [version]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                where_clauses.append(f"p.taxon_id IN ({placeholders})")
                params.extend(taxon_ids)
            else:
                where_clauses.append("p.taxon_id = %s")
                params.append(taxon_ids)

        if proteome_id:
            where_clauses.append("p.proteome_id = %s")
            params.append(proteome_id)

        if go_id:
            where_clauses.append("pg.go_id = %s")
            params.append(go_id)

        if pfam_id:
            where_clauses.append("pp.pfam_id = %s")
            params.append(pfam_id)

        query += " WHERE " + " AND ".join(where_clauses)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ------------------------------------------------------------------
    # HMM hit based retrieval
    # ------------------------------------------------------------------
    def get_proteins_by_hmm_hit(
        self,
        version,
        hmm_query,
        evalue_cutoff=None,
        taxon_ids=None,
    ):
        """
        Fetch sequences for all proteins with a hit to a given HMM profile.
        Suitable for phylogenetic tree building pipelines.

        Queries the local hmm_search_results table that contains results from a full
        Pfam-A hmmsearch using gathering thresholds. Results are already
        filtered at the profile-specific trusted cutoff. So, evalue_cutoff
        provides an additional filter on top of that.

        Results are deduplicated: a protein with five Homeodomain repeats
        appears only once in the output, since only its sequence is needed.

        Suitable for phylogenetic tree building pipelines.

        Parameters
        ----------
        version : str
            UniProt release version, e.g. "2026_01".
        hmm_query : str
            Pfam name (e.g., "Homeodomain") OR Accession (e.g., "PF00046").
        evalue_cutoff : float, optional
            Maximum full sequence E-value to consider a hit valid.
            Defaults to None.
        taxon_ids : int or list[int], optional
            Filter by specific NCBI Taxonomy IDs.

        Returns
        -------
        list[dict]
            Each dict has keys: accession, name, organism, taxon_id,
            proteome_id, sequence.

        Examples
        --------
        # All Homeodomain proteins across all taxa:
        records = db.get_proteins_by_hmm_hit("2026_01", "Homeodomain")

        # Human + mouse kinases with strict E-value:
        records = db.get_proteins_by_hmm_hit(
            "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
        )
        seqrecords = db.to_biopython(records)
        db.export_fasta(records, "2026_01", "kinase_human_mouse", output_dir="./fastas")
        """

        # Uses DISTINCT because a protein might have 5 repeating Homeodomains, but we only need to fetch its sequence once.
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""
        query = f"""
            SELECT DISTINCT p.accession, p.name, p.organism,
                            p.taxon_id, p.proteome_id, s.sequence
            FROM   hmm_search_results h
            JOIN   proteins p  ON h.accession = p.accession AND h.version = p.version
            JOIN   sequences s ON p.seq_id = s.seq_id
            WHERE  h.version = %s
              {evalue_clause}
              AND  (h.hmm_name = %s OR h.hmm_accession LIKE %s)
        """


        # We use LIKE %s with a wildcard (%) so if the user searches for "PF00046",
        # it successfully matches "PF00046.36" in your database.
        params = [version]
        if evalue_cutoff is not None:
            params.append(evalue_cutoff)
        params += [hmm_query, f"{hmm_query}%"]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                query += f" AND h.taxon_id IN ({placeholders})"
                params.extend(taxon_ids)
            else:
                query += " AND h.taxon_id = %s"
                params.append(taxon_ids)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err


    # ------------------------------------------------------------------
    # Domain lookup by accession
    # ------------------------------------------------------------------

    def get_domains_by_accession(
        self,
        version,
        accessions,
        evalue_cutoff=None,
        chunk_size=5000,
    ):
        """
        Return all Pfam domain hits from hmm_search_results for one or
        more protein accessions.

        One row is returned per domain occurrence — a protein with five
        Ankyrin repeats produces five rows, one per repeat. This is
        intentional: use this method when you need domain coordinates and
        scores, not just the sequence.

        Large accession lists are chunked automatically (default 5,000
        per query) so the IN (...) clause never becomes a bottleneck.

        Parameters
        ----------
        version : str
            UniProt release version, e.g. "2026_01".
        accessions : str or list[str]
            One or more UniProt accessions, e.g. ["P04637", "P10275"].
        evalue_cutoff : float, optional
            Maximum full-sequence E-value. Default: 1e-5.
            Note: all results already passed gathering thresholds during
            the hmmsearch — this is an additional filter on top of that.
        chunk_size : int, optional
            Maximum number of accessions per SQL query. Default: 5000.

        Returns
        -------
        list[dict]
            One dict per domain occurrence. Keys:
            accession, protein_name, organism, taxon_id, proteome_id,
            hmm_name, hmm_accession,
            domain_number, domain_count,
            full_evalue, full_score,
            domain_evalue, domain_score,
            ali_from, ali_to, hmm_from, hmm_to, env_from, env_to.
            Sorted by accession, then ali_from.

        Examples
        --------
        # Single protein — all its Pfam domains
        domains = db.get_domains_by_accession("2026_01", "P04637")

        # Multiple proteins
        domains = db.get_domains_by_accession(
            "2026_01", ["P04637", "P10275", "Q8BFR5"]
        )
        for d in domains:
            print(f"{d['accession']}  {d['hmm_name']}  "
                  f"{d['ali_from']}–{d['ali_to']}  E={d['full_evalue']:.2e}")
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        query_template = """
            SELECT
                h.accession,
                p.name          AS protein_name,
                p.organism,
                p.taxon_id,
                p.proteome_id,
                h.hmm_name,
                h.hmm_accession,
                h.domain_number,
                h.domain_count,
                h.full_evalue,
                h.full_score,
                h.domain_evalue,
                h.domain_score,
                h.ali_from,
                h.ali_to,
                h.hmm_from,
                h.hmm_to,
                h.env_from,
                h.env_to
            FROM  hmm_search_results h
            JOIN  proteins p
              ON  h.accession = p.accession
             AND  h.version   = p.version
            WHERE h.version      = %s
              {evalue_clause}
              AND h.accession  IN ({{placeholders}})
            ORDER BY h.accession, h.ali_from
        """

        all_results = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))
                query  = query_template.format(
                    evalue_clause=evalue_clause,
                    placeholders=placeholders,
                )
                params = [version]
                if evalue_cutoff is not None:
                    params.append(evalue_cutoff)
                params += chunk
                self.cursor.execute(query, tuple(params))
                all_results.extend(self.cursor.fetchall())
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        return all_results


    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def to_fasta_string(self, records):
        """
        Convert a list of protein records to a FASTA-formatted string.

        This does NOT write any file — the string is returned so callers
        can pass it directly to parsers (e.g. Bio.SeqIO.parse via StringIO).

        Parameters
        ----------
        records : list[dict]  output of get_proteins()

        Returns
        -------
        str  multi-sequence FASTA string
        """
        lines = []
        for rec in records:
            header = f">{rec['taxon_id']}.{rec['accession']}"
            lines.append(f"{header}\n{rec['sequence']}")
        return "\n".join(lines) + "\n"

    def to_biopython(self, records):
        """
        Convert records to a list of BioPython SeqRecord objects.

        Requires BioPython (pip install biopython).
        Useful for alignment tools, tree building (ETE3, DendroPy, etc.),
        and any workflow that expects SeqRecord objects.

        Parameters
        ----------
        records : list[dict]  output of get_proteins()

        Returns
        -------
        list[Bio.SeqRecord.SeqRecord]
            Each SeqRecord has:
            - id          = "{taxon_id}.{accession}"
            - name        = "{accession}"
            - description = "{protein_name} [{organism}] UP={proteome_id}"
        """
        try:
            from Bio.Seq      import Seq
            from Bio.SeqRecord import SeqRecord
        except ImportError as e:
            raise ImportError(
                "BioPython is required for to_biopython(). "
                "Install it with:  pip install biopython"
            ) from e

        seqrecords = []
        for rec in records:
            sr = SeqRecord(
                Seq(rec["sequence"]),
                id          = f"{rec['taxon_id']}.{rec['accession']}",
                name        = rec["accession"],
                description = f"{rec['name']} [{rec['organism']}] "
                              f"UP={rec['proteome_id']}",
            )
            seqrecords.append(sr)
        return seqrecords

    def export_fasta(self, records, version, identifier, output_dir=None, filename=None):
        """
        Write protein records to a FASTA file.

        Parameters
        ----------
        records : list[dict]  output of get_proteins()
        version : str         e.g. "2026_01"  (used in auto-generated filename)
        identifier : str      label for the file  (e.g. proteome ID, taxon tag)
        output_dir : str, optional
            Directory to write the file into.  Defaults to the current
            working directory.  Created automatically if it does not exist.
        filename : str, optional
            Override the auto-generated filename entirely.
            If not given, the file is named:
            uniprot_<identifier>_<version>.fasta

        Returns
        -------
        str  absolute path of the written file
        """
        if filename is None:
            filename = f"uniprot_{identifier}_{version}.fasta"

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            filepath = os.path.join(output_dir, filename)
        else:
            filepath = filename

        with open(filepath, "w") as f:
            f.write(self.to_fasta_string(records))

        print(f" Exported {len(records):,} sequences → {os.path.abspath(filepath)}")
        return os.path.abspath(filepath)


# ---------------------------------------------------------------------------
# Module-level convenience functions
# (users can call these without managing a class instance)
# ---------------------------------------------------------------------------

def fetch_sequences(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → disconnect → return records.

    Parameters
    ----------
    version : str
        UniProt release version, e.g. "2026_01".
    taxon_ids : int or list[int], optional
    proteome_id : str, optional
    go_id : str, optional
    pfam_id : str, optional
    db_config : dict, optional
        Override connection parameters. Defaults to get_db_config()
        (reads from environment variables / .env file).

    Returns
    -------
    list[dict]  keys: accession, name, organism, taxon_id, proteome_id, sequence

    Example
    -------
    >>> from get_reference_uniprot_set_lib import fetch_sequences
    >>> records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])
    >>> print(len(records), "proteins retrieved")
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins(
            version     = version,
            taxon_ids   = taxon_ids,
            proteome_id = proteome_id,
            go_id       = go_id,
            pfam_id     = pfam_id,
        )


def fetch_sequences_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=None,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → disconnect → return records.

    Retrieves sequences for all proteins with a hit to the given Pfam HMM profile.
    Results are deduplicated — one sequence per protein regardless of domain copy count.
    Suitable as a direct input to alignment and tree-building pipelines.

    Parameters
    ----------
    version : str
        UniProt release version, e.g. "2026_01".
    hmm_query : str
        Pfam name (e.g. "Homeodomain") or accession (e.g. "PF00046").
        Accession matching is prefix-based so "PF00046" matches "PF00046.36".
    evalue_cutoff : float, optional
        Maximum full-sequence E-value. Default: 1e-5.
    taxon_ids : int or list[int], optional
        Filter by one or more NCBI Taxonomy IDs.
    db_config : dict, optional
        Override connection parameters. Defaults to get_db_config().

    Returns
    -------
    list[dict]  keys: accession, name, organism, taxon_id, proteome_id, sequence

    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit
    >>> # All Homeodomain proteins across all taxa
    >>> records = fetch_sequences_by_hmm_hit("2026_01", "Homeodomain")
    >>> print(len(records), "proteins retrieved")

    >>> # Human + mouse kinases, strict E-value, pipe into BioPython
    >>> from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config
    >>> records = fetch_sequences_by_hmm_hit(
    ...     "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
    ... )
    >>> with UniProtRetriever(get_db_config()) as db:
    ...     seqrecords = db.to_biopython(records)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins_by_hmm_hit(
            version       = version,
            hmm_query     = hmm_query,
            evalue_cutoff = evalue_cutoff,
            taxon_ids     = taxon_ids,
        )


def fetch_fasta_string(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → return a FASTA string.

    Useful for piping directly into Bio.SeqIO.parse() without writing a file.

    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string
    >>> fasta = fetch_fasta_string("2026_01", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins(
            version     = version,
            taxon_ids   = taxon_ids,
            proteome_id = proteome_id,
            go_id       = go_id,
            pfam_id     = pfam_id,
        )
        return db.to_fasta_string(records)

def fetch_fasta_string_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=None,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → return a FASTA string.

    Useful for piping HMM-filtered sequences directly into Bio.SeqIO.parse()
    without writing a file.

    Parameters
    ----------
    version : str
    hmm_query : str   Pfam name or accession
    evalue_cutoff : float, optional  Default: 1e-5
    taxon_ids : int or list[int], optional
    db_config : dict, optional

    Returns
    -------
    str  multi-sequence FASTA string

    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string_by_hmm_hit
    >>> fasta = fetch_fasta_string_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins_by_hmm_hit(
            version       = version,
            hmm_query     = hmm_query,
            evalue_cutoff = evalue_cutoff,
            taxon_ids     = taxon_ids,
        )
        return db.to_fasta_string(records)


def fetch_domains_by_accession(
    version,
    accessions,
    evalue_cutoff=None,
    chunk_size=5000,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → disconnect →
    return all domain hits for the given accession(s).

    Returns one dict per domain occurrence (not one per protein). Use
    this when you need domain coordinates, scores, and HMM accessions
    rather than sequences.

    Parameters
    ----------
    version : str
        UniProt release version, e.g. "2026_01".
    accessions : str or list[str]
        One or more UniProt protein accessions.
    evalue_cutoff : float, optional
        Maximum full-sequence E-value. Default: 1e-5.
    chunk_size : int, optional
        Accessions per SQL query — chunked automatically. Default: 5000.
    db_config : dict, optional
        Override connection parameters. Defaults to get_db_config().

    Returns
    -------
    list[dict]
        Keys: accession, protein_name, organism, taxon_id, proteome_id,
        hmm_name, hmm_accession, domain_number, domain_count,
        full_evalue, full_score, domain_evalue, domain_score,
        ali_from, ali_to, hmm_from, hmm_to, env_from, env_to.

    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_domains_by_accession
    >>> domains = fetch_domains_by_accession("2026_01", "P04637")
    >>> for d in domains:
    ...     print(d["hmm_name"], d["hmm_accession"], d["ali_from"], d["ali_to"])

    >>> # Works with any size list — chunked automatically
    >>> domains = fetch_domains_by_accession("2026_01", my_list_of_1000_accessions)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_domains_by_accession(
            version       = version,
            accessions    = accessions,
            evalue_cutoff = evalue_cutoff,
            chunk_size    = chunk_size,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description="UniProt Reference Set Retrieval Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-version",        required=True, help="UniProt version (e.g., 2026_01)")
    parser.add_argument("-taxonomy",       nargs="+", type=int, help="One or more Taxonomy IDs")
    parser.add_argument("--proteome-id",   help="Filter by Proteome ID (e.g., UP000005640)")
    parser.add_argument("--go-id",         help="Filter by GO ID (e.g., GO:0005634)")
    parser.add_argument("--pfam-id",       help="Filter by Pfam ID (e.g., PF00870)")
    parser.add_argument("--list-versions", action="store_true", help="List all available versions")
    parser.add_argument("--list-proteomes",action="store_true", help="List proteome IDs for this version")
    parser.add_argument("--output-dir",    default=None, help="Directory for the output FASTA file")
    parser.add_argument("--hmm-name",  default=None, help="Pfam HMM name or accession, e.g. Homeodomain or PF00046")
    parser.add_argument("--evalue",    type=float, default=1e-5, help="Max full-sequence E-value for HMM hits (default: 1e-5)")
    return parser


def main():
    args = _build_parser().parse_args()

    retriever = UniProtRetriever(get_db_config())

    try:
        retriever.connect()

        if args.list_versions:
            retriever.list_available_versions()
            return

        if args.list_proteomes:
            proteomes = retriever.get_proteome_ids(args.version)
            print(f"\nProteomes in version {args.version}:")
            for p_id in proteomes:
                print(f"  - {p_id}")
            return

        print(f"\n{'='*60}")
        print(f"UniProt Reference Set Retrieval")
        print(f"Version : {args.version}")
        if args.taxonomy:    print(f"Taxonomy: {args.taxonomy}")
        if args.proteome_id: print(f"Proteome: {args.proteome_id}")
        if args.go_id:       print(f"GO Term : {args.go_id}")
        if args.pfam_id:     print(f"Pfam ID : {args.pfam_id}")

        if args.hmm_name:
            print(f"HMM    : {args.hmm_name} (E-value ≤ {args.evalue})")
            print(f"{'='*60}\n")
            records = retriever.get_proteins_by_hmm_hit(
                version       = args.version,
                hmm_query     = args.hmm_name,
                evalue_cutoff = args.evalue,
                taxon_ids     = args.taxonomy,
            )
            identifier = f"hmm_{args.hmm_name}"
        else:
            print(f"{'='*60}\n")
            records = retriever.get_proteins(
                version     = args.version,
                taxon_ids   = args.taxonomy,
                proteome_id = args.proteome_id,
                go_id       = args.go_id,
                pfam_id     = args.pfam_id,
            )
            identifier = "filtered_set"
            if args.proteome_id:
                identifier = args.proteome_id
            elif args.taxonomy:
                identifier = f"tax_{args.taxonomy[0]}"

        if not records:
            print("\n✗ No matching data found.")
            return

        retriever.export_fasta(
            records, args.version, identifier,
            output_dir=args.output_dir,
        )
        print(f"\nSuccessfully retrieved {len(records):,} sequences.")

    except mysql.connector.Error as err:
        print(f"\n✗ Database error: {err}")
        sys.exit(1)
    except Exception as e:
        print(f"\n Error: {e}")
        sys.exit(1)
    finally:
        retriever.close()


if __name__ == "__main__":
    main()
"""
UniProt Reference Set Retrieval Library
========================================
Retrieves protein sequences from the local UniProt database.

Can be used in two ways:

1. As a command-line tool (for quick retrieval without writing custom scripts):
   python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090

2. As an importable library in your scripts:

   # --- Quickstart (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences

   records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
   for r in records:
       print(r["accession"], r["organism"])

   # --- HMM hit retrieval (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit

   records = fetch_sequences_by_hmm_hit(
       version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
   )

   # --- Context manager (recommended for multiple queries) ---
   from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

   with UniProtRetriever(get_db_config()) as db:
       records    = db.get_proteins(version="2026_01", taxon_ids=[9606, 10090])
       hmm_recs   = db.get_proteins_by_hmm_hit(version="2026_01", hmm_query="Homeodomain")
       fasta_str  = db.to_fasta_string(records)       # in-memory FASTA string
       seqrecords = db.to_biopython(records)           # list of BioPython SeqRecord
       db.export_fasta(records, "2026_01", "mammals", output_dir="./fastas")
"""

import mysql.connector
import argparse
import os
import sys
from dotenv import load_dotenv

load_dotenv()

# Public API — what users get when they do `from ... import *`
__all__ = [
    "UniProtRetriever",
    "get_db_config",
    "fetch_sequences",
    "fetch_sequences_by_hmm_hit",
    "fetch_fasta_string",
    "fetch_fasta_string_by_hmm_hit",
]


# ---------------------------------------------------------------------------
# Configuration helper
# ---------------------------------------------------------------------------

def get_db_config(host=None, user=None, password=None, database=None):
    """
    Build the database config dict, falling back to environment variables
    and then to the lab defaults.

    users can override any value when calling this function, or by
    setting environment variables in a .env file:
        DB_HOST, DB_USER, DB_PASSWORD, DB_NAME

    Parameters
    ----------
    host, user, password, database : str, optional
        Explicit overrides — take priority over everything else.

    Returns
    -------
    dict  ready to pass into UniProtRetriever(config)
    """
    return {
        "host":     host     or os.getenv("DB_HOST",     "localhost"),
        "user":     user     or os.getenv("DB_USER",     "user"),
        "password": password or os.getenv("DB_PASSWORD", "your_password"),
        "database": database or os.getenv("DB_NAME",     "uniprot_db"),
    }


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

class UniProtRetriever:
    """
    Retrieves UniProt reference sets from the local CGLab database.

    Supports filtering by TaxID, Proteome ID, Pfam domain, GO term, and HMM search
    results from the local Pfam-A hmmsearch conducted by pyhmmer hmmsearch.

    Results can be returned as raw dicts, a FASTA string, BioPython
    SeqRecord objects, or written directly to a .fasta file.

    Recommended usage — context manager (handles connect/close):

        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins(version="2026_01", taxon_ids=[9606])
            seqs = db.to_biopython(records)

        # HMM-based retrieval for tree building:
        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins_by_hmm_hit(
                version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
            )
            seqs = db.to_biopython(records)
    """

    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict
            mysql.connector connection parameters.
            Use get_db_config() to build this from environment variables.
        """
        self.config = config
        self.conn   = None
        self.cursor = None

    # ------------------------------------------------------------------
    # Connection management
    # ------------------------------------------------------------------

    def connect(self):
        """
        Open the database connection.

        """
        self.conn   = mysql.connector.connect(**self.config)
        self.cursor = self.conn.cursor(dictionary=True)

    def close(self):
        """Close the database connection (safe to call even if not connected)."""
        if self.conn and self.conn.is_connected():
            self.cursor.close()
            self.conn.close()

    # Context-manager support — in order to be able to use `with UniProtRetriever(...) as db:`
    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    # ------------------------------------------------------------------
    # Utility queries
    # ------------------------------------------------------------------

    def list_available_versions(self):
        """
        Print a summary table of all versions present in the database.

        Returns
        -------
        list[dict]  one dict per version with keys:
                    version, protein_count, taxon_count, proteome_count
        """
        query = """
            SELECT version,
                   COUNT(*)                  AS protein_count,
                   COUNT(DISTINCT taxon_id)  AS taxon_count,
                   COUNT(DISTINCT proteome_id) AS proteome_count
            FROM   proteins
            GROUP  BY version
            ORDER  BY version DESC
        """
        self.cursor.execute(query)
        versions = self.cursor.fetchall()
        if versions:
            print("\nAvailable versions in database:")
            print("-" * 70)
            for v in versions:
                print(
                    f"Version: {v['version']} | "
                    f"Proteins: {v['protein_count']:,} | "
                    f"Taxa: {v['taxon_count']:,} | "
                    f"Proteomes: {v['proteome_count']:,}"
                )
            print("-" * 70)
        else:
            print("\nNo versions found in database.")
        return versions

    def get_proteome_ids(self, version):
        """
        Return all unique Proteome IDs present for a given version.

        Parameters
        ----------
        version : str  e.g. "2026_01"

        Returns
        -------
        list[str]
        """
        query = "SELECT DISTINCT proteome_id FROM proteins WHERE version = %s"
        self.cursor.execute(query, (version,))
        return [row["proteome_id"] for row in self.cursor.fetchall()]

    # ------------------------------------------------------------------
    # Main retrieval
    # ------------------------------------------------------------------

    def get_proteins(
        self,
        version,
        taxon_ids=None,
        proteome_id=None,
        go_id=None,
        pfam_id=None,
    ):
        """
        Retrieve protein records from the database.

        All filters are optional and combinable.

        Parameters
        ----------
        version : str
            UniProt release version, e.g. "2026_01".
        taxon_ids : int or list[int], optional
            One or more NCBI Taxonomy IDs.
        proteome_id : str, optional
            UniProt Proteome ID, e.g. "UP000005640".
        go_id : str, optional
            Gene Ontology term, e.g. "GO:0005634".
        pfam_id : str, optional
            Pfam domain accession, e.g. "PF00870".

        Returns
        -------
        list[dict]
            Each dict has keys:
            accession, name, organism, taxon_id, proteome_id, sequence
        """
        query = """
            SELECT p.accession, p.name, p.organism,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM   proteins  p
            JOIN   sequences s ON p.seq_id = s.seq_id
        """
        if go_id:
            query += (
                " JOIN protein_go pg"
                " ON p.accession = pg.accession AND p.version = pg.version"
            )
        if pfam_id:
            query += (
                " JOIN protein_pfam pp"
                " ON p.accession = pp.accession AND p.version = pp.version"
            )

        where_clauses = ["p.version = %s"]
        params        = [version]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                where_clauses.append(f"p.taxon_id IN ({placeholders})")
                params.extend(taxon_ids)
            else:
                where_clauses.append("p.taxon_id = %s")
                params.append(taxon_ids)

        if proteome_id:
            where_clauses.append("p.proteome_id = %s")
            params.append(proteome_id)

        if go_id:
            where_clauses.append("pg.go_id = %s")
            params.append(go_id)

        if pfam_id:
            where_clauses.append("pp.pfam_id = %s")
            params.append(pfam_id)

        query += " WHERE " + " AND ".join(where_clauses)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ------------------------------------------------------------------
    # HMM hit based retrieval
    # ------------------------------------------------------------------
    def get_proteins_by_hmm_hit(
        self,
        version,
        hmm_query,
        evalue_cutoff=1e-5,
        taxon_ids=None,
    ):
        """
        Fetch sequences for all proteins with a hit to a given HMM profile.
        Suitable for phylogenetic tree building pipelines.

        Queries the local hmm_search_results table that contains results from a full
        Pfam-A hmmsearch using gathering thresholds. Results are already
        filtered at the profile-specific trusted cutoff. So, evalue_cutoff
        provides an additional filter on top of that.

        Results are deduplicated: a protein with five Homeodomain repeats
        appears only once in the output, since only its sequence is needed.

        Suitable for phylogenetic tree building pipelines.

        Parameters
        ----------
        version : str
            UniProt release version, e.g. "2026_01".
        hmm_query : str
            Pfam name (e.g., "Homeodomain") OR Accession (e.g., "PF00046").
        evalue_cutoff : float
            Maximum full sequence E-value to consider a hit valid.
            Defaults to 1e-5.
        taxon_ids : int or list[int], optional
            Filter by specific NCBI Taxonomy IDs.

        Returns
        -------
        list[dict]
            Each dict has keys: accession, name, organism, taxon_id,
            proteome_id, sequence.

        Examples
        --------
        # All Homeodomain proteins across all taxa:
        records = db.get_proteins_by_hmm_hit("2026_01", "Homeodomain")

        # Human + mouse kinases with strict E-value:
        records = db.get_proteins_by_hmm_hit(
            "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
        )
        seqrecords = db.to_biopython(records)
        db.export_fasta(records, "2026_01", "kinase_human_mouse", output_dir="./fastas")
        """

        # Uses DISTINCT because a protein might have 5 repeating Homeodomains, but we only need to fetch its sequence once.

        query = """
            SELECT DISTINCT p.accession, p.name, p.organism,
                            p.taxon_id, p.proteome_id, s.sequence
            FROM   hmm_search_results h
            JOIN   proteins p  ON h.accession = p.accession AND h.version = p.version
            JOIN   sequences s ON p.seq_id = s.seq_id
            WHERE  h.version = %s
              AND  h.full_evalue <= %s
              AND  (h.hmm_name = %s OR h.hmm_accession LIKE %s)
    """

        # We use LIKE %s with a wildcard (%) so if the user searches for "PF00046",
        # it successfully matches "PF00046.36" in your database.
        params = [version, evalue_cutoff, hmm_query, f"{hmm_query}%"]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                query += f" AND h.taxon_id IN ({placeholders})"
                params.extend(taxon_ids)
            else:
                query += " AND h.taxon_id = %s"
                params.append(taxon_ids)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def to_fasta_string(self, records):
        """
        Convert a list of protein records to a FASTA-formatted string.

        This does NOT write any file — the string is returned so callers
        can pass it directly to parsers (e.g. Bio.SeqIO.parse via StringIO).

        Parameters
        ----------
        records : list[dict]  output of get_proteins()

        Returns
        -------
        str  multi-sequence FASTA string
        """
        lines = []
        for rec in records:
            header = f">{rec['taxon_id']}.{rec['accession']}"
            lines.append(f"{header}\n{rec['sequence']}")
        return "\n".join(lines) + "\n"

    def to_biopython(self, records):
        """
        Convert records to a list of BioPython SeqRecord objects.

        Requires BioPython (pip install biopython).
        Useful for alignment tools, tree building (ETE3, DendroPy, etc.),
        and any workflow that expects SeqRecord objects.

        Parameters
        ----------
        records : list[dict]  output of get_proteins()

        Returns
        -------
        list[Bio.SeqRecord.SeqRecord]
            Each SeqRecord has:
            - id          = "{taxon_id}.{accession}"
            - name        = "{accession}"
            - description = "{protein_name} [{organism}] UP={proteome_id}"
        """
        try:
            from Bio.Seq      import Seq
            from Bio.SeqRecord import SeqRecord
        except ImportError as e:
            raise ImportError(
                "BioPython is required for to_biopython(). "
                "Install it with:  pip install biopython"
            ) from e

        seqrecords = []
        for rec in records:
            sr = SeqRecord(
                Seq(rec["sequence"]),
                id          = f"{rec['taxon_id']}.{rec['accession']}",
                name        = rec["accession"],
                description = f"{rec['name']} [{rec['organism']}] "
                              f"UP={rec['proteome_id']}",
            )
            seqrecords.append(sr)
        return seqrecords

    def export_fasta(self, records, version, identifier, output_dir=None, filename=None):
        """
        Write protein records to a FASTA file.

        Parameters
        ----------
        records : list[dict]  output of get_proteins()
        version : str         e.g. "2026_01"  (used in auto-generated filename)
        identifier : str      label for the file  (e.g. proteome ID, taxon tag)
        output_dir : str, optional
            Directory to write the file into.  Defaults to the current
            working directory.  Created automatically if it does not exist.
        filename : str, optional
            Override the auto-generated filename entirely.
            If not given, the file is named:
            uniprot_<identifier>_<version>.fasta

        Returns
        -------
        str  absolute path of the written file
        """
        if filename is None:
            filename = f"uniprot_{identifier}_{version}.fasta"

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            filepath = os.path.join(output_dir, filename)
        else:
            filepath = filename

        with open(filepath, "w") as f:
            f.write(self.to_fasta_string(records))

        print(f" Exported {len(records):,} sequences → {os.path.abspath(filepath)}")
        return os.path.abspath(filepath)


# ---------------------------------------------------------------------------
# Module-level convenience functions
# (users can call these without managing a class instance)
# ---------------------------------------------------------------------------

def fetch_sequences(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → disconnect → return records.

    Parameters
    ----------
    version : str
        UniProt release version, e.g. "2026_01".
    taxon_ids : int or list[int], optional
    proteome_id : str, optional
    go_id : str, optional
    pfam_id : str, optional
    db_config : dict, optional
        Override connection parameters. Defaults to get_db_config()
        (reads from environment variables / .env file).

    Returns
    -------
    list[dict]  keys: accession, name, organism, taxon_id, proteome_id, sequence

    Example
    -------
    >>> from get_reference_uniprot_set_lib import fetch_sequences
    >>> records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])
    >>> print(len(records), "proteins retrieved")
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins(
            version     = version,
            taxon_ids   = taxon_ids,
            proteome_id = proteome_id,
            go_id       = go_id,
            pfam_id     = pfam_id,
        )


def fetch_sequences_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=1e-5,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → disconnect → return records.

    Retrieves sequences for all proteins with a hit to the given Pfam HMM profile.
    Results are deduplicated — one sequence per protein regardless of domain copy count.
    Suitable as a direct input to alignment and tree-building pipelines.

    Parameters
    ----------
    version : str
        UniProt release version, e.g. "2026_01".
    hmm_query : str
        Pfam name (e.g. "Homeodomain") or accession (e.g. "PF00046").
        Accession matching is prefix-based so "PF00046" matches "PF00046.36".
    evalue_cutoff : float, optional
        Maximum full-sequence E-value. Default: 1e-5.
    taxon_ids : int or list[int], optional
        Filter by one or more NCBI Taxonomy IDs.
    db_config : dict, optional
        Override connection parameters. Defaults to get_db_config().

    Returns
    -------
    list[dict]  keys: accession, name, organism, taxon_id, proteome_id, sequence

    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit
    >>> # All Homeodomain proteins across all taxa
    >>> records = fetch_sequences_by_hmm_hit("2026_01", "Homeodomain")
    >>> print(len(records), "proteins retrieved")

    >>> # Human + mouse kinases, strict E-value, pipe into BioPython
    >>> from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config
    >>> records = fetch_sequences_by_hmm_hit(
    ...     "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
    ... )
    >>> with UniProtRetriever(get_db_config()) as db:
    ...     seqrecords = db.to_biopython(records)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins_by_hmm_hit(
            version       = version,
            hmm_query     = hmm_query,
            evalue_cutoff = evalue_cutoff,
            taxon_ids     = taxon_ids,
        )


def fetch_fasta_string(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → return a FASTA string.

    Useful for piping directly into Bio.SeqIO.parse() without writing a file.

    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string
    >>> fasta = fetch_fasta_string("2026_01", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins(
            version     = version,
            taxon_ids   = taxon_ids,
            proteome_id = proteome_id,
            go_id       = go_id,
            pfam_id     = pfam_id,
        )
        return db.to_fasta_string(records)

def fetch_fasta_string_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=1e-5,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → return a FASTA string.

    Useful for piping HMM-filtered sequences directly into Bio.SeqIO.parse()
    without writing a file.

    Parameters
    ----------
    version : str
    hmm_query : str   Pfam name or accession
    evalue_cutoff : float, optional  Default: 1e-5
    taxon_ids : int or list[int], optional
    db_config : dict, optional

    Returns
    -------
    str  multi-sequence FASTA string

    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string_by_hmm_hit
    >>> fasta = fetch_fasta_string_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins_by_hmm_hit(
            version       = version,
            hmm_query     = hmm_query,
            evalue_cutoff = evalue_cutoff,
            taxon_ids     = taxon_ids,
        )
        return db.to_fasta_string(records)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description="UniProt Reference Set Retrieval Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-version",        required=True, help="UniProt version (e.g., 2026_01)")
    parser.add_argument("-taxonomy",       nargs="+", type=int, help="One or more Taxonomy IDs")
    parser.add_argument("--proteome-id",   help="Filter by Proteome ID (e.g., UP000005640)")
    parser.add_argument("--go-id",         help="Filter by GO ID (e.g., GO:0005634)")
    parser.add_argument("--pfam-id",       help="Filter by Pfam ID (e.g., PF00870)")
    parser.add_argument("--list-versions", action="store_true", help="List all available versions")
    parser.add_argument("--list-proteomes",action="store_true", help="List proteome IDs for this version")
    parser.add_argument("--output-dir",    default=None, help="Directory for the output FASTA file")
    parser.add_argument("--hmm-name",  default=None, help="Pfam HMM name or accession, e.g. Homeodomain or PF00046")
    parser.add_argument("--evalue",    type=float, default=1e-5, help="Max full-sequence E-value for HMM hits (default: 1e-5)")
    return parser


def main():
    args = _build_parser().parse_args()

    retriever = UniProtRetriever(get_db_config())

    try:
        retriever.connect()

        if args.list_versions:
            retriever.list_available_versions()
            return

        if args.list_proteomes:
            proteomes = retriever.get_proteome_ids(args.version)
            print(f"\nProteomes in version {args.version}:")
            for p_id in proteomes:
                print(f"  - {p_id}")
            return

        print(f"\n{'='*60}")
        print(f"UniProt Reference Set Retrieval")
        print(f"Version : {args.version}")
        if args.taxonomy:    print(f"Taxonomy: {args.taxonomy}")
        if args.proteome_id: print(f"Proteome: {args.proteome_id}")
        if args.go_id:       print(f"GO Term : {args.go_id}")
        if args.pfam_id:     print(f"Pfam ID : {args.pfam_id}")

        if args.hmm_name:
            print(f"HMM    : {args.hmm_name} (E-value ≤ {args.evalue})")
            print(f"{'='*60}\n")
            records = retriever.get_proteins_by_hmm_hit(
                version       = args.version,
                hmm_query     = args.hmm_name,
                evalue_cutoff = args.evalue,
                taxon_ids     = args.taxonomy,
            )
            identifier = f"hmm_{args.hmm_name}"
        else:
            print(f"{'='*60}\n")
            records = retriever.get_proteins(
                version     = args.version,
                taxon_ids   = args.taxonomy,
                proteome_id = args.proteome_id,
                go_id       = args.go_id,
                pfam_id     = args.pfam_id,
            )
            identifier = "filtered_set"
            if args.proteome_id:
                identifier = args.proteome_id
            elif args.taxonomy:
                identifier = f"tax_{args.taxonomy[0]}"

        if not records:
            print("\n✗ No matching data found.")
            return

        retriever.export_fasta(
            records, args.version, identifier,
            output_dir=args.output_dir,
        )
        print(f"\nSuccessfully retrieved {len(records):,} sequences.")

    except mysql.connector.Error as err:
        print(f"\n✗ Database error: {err}")
        sys.exit(1)
    except Exception as e:
        print(f"\n Error: {e}")
        sys.exit(1)
    finally:
        retriever.close()


if __name__ == "__main__":
    main()
