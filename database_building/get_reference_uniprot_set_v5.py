"""
UniProt Reference Set Retrieval Tool — v5
==========================================
Retrieves specific versions of UniProt reference proteomes from local database.

Changes for v5:
────────────────────
1. Updated base query to JOIN the new taxa table — organism and division
   are now fetched from taxa instead of the proteins table.

2. Added --division filter (Archaea, Bacteria, Eukaryota, Viruses).

3. Updated FASTA headers to standard UniProt format:
   >ACCESSION ENTRY_NAME OS=Organism OX=TaxonID UP=ProteomeID DIV=Division

4. Removed unused --filtered-headers argument.

5. Removed unused import_datetime helper function.

Usage:
    python get_reference_uniprot_set_v5.py -version 2026_01
    python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606
    python get_reference_uniprot_set_v5.py -version 2026_01 -taxonomy 9606 10090
    python get_reference_uniprot_set_v5.py -version 2026_01 --proteome-id UP000005640
    python get_reference_uniprot_set_v5.py -version 2026_01 --go-id GO:0005634
    python get_reference_uniprot_set_v5.py -version 2026_01 --pfam-id PF00870
    python get_reference_uniprot_set_v5.py -version 2026_01 --division Eukaryota
    python get_reference_uniprot_set_v5.py -version 2026_01 --list-versions
    python get_reference_uniprot_set_v5.py -version 2026_01 --list-proteomes
"""

import mysql.connector
import argparse
import os
from dotenv import load_dotenv
import sys

load_dotenv()


class UniProtRetriever:
    """
    Retrieves UniProt reference sets from the local database.
    Supports filtering by TaxID, Proteome ID, GO term, Pfam domain,
    and taxonomic division.
    """

    def __init__(self, config):
        self.config = config
        self.conn = None
        self.cursor = None

    def connect(self):
        """Establish database connection."""
        try:
            self.conn = mysql.connector.connect(**self.config)
            self.cursor = self.conn.cursor(dictionary=True)
            print("Connected to database successfully")
        except mysql.connector.Error as err:
            print(f"Database connection error: {err}")
            sys.exit(1)

    def list_available_versions(self):
        """List all available versions in the database."""
        query = """
            SELECT version,
                   COUNT(*) as protein_count,
                   COUNT(DISTINCT taxon_id) as taxon_count,
                   COUNT(DISTINCT proteome_id) as proteome_count
            FROM proteins
            GROUP BY version
            ORDER BY version DESC
        """
        try:
            self.cursor.execute(query)
            versions = self.cursor.fetchall()
            if versions:
                print("\nAvailable versions in database:")
                print("-" * 70)
                for v in versions:
                    print(
                        f"Version: {v['version']} | Proteins: {v['protein_count']:,} | "
                        f"Taxa: {v['taxon_count']:,} | Proteomes: {v['proteome_count']:,}"
                    )
                print("-" * 70)
            else:
                print("\nNo versions found in database.")
        except mysql.connector.Error as err:
            print(f"Error listing versions: {err}")

    def get_proteome_ids(self, version):
        """Returns all unique Proteome IDs for a specific version."""
        query = "SELECT DISTINCT proteome_id FROM proteins WHERE version = %s"
        self.cursor.execute(query, (version,))
        return [row["proteome_id"] for row in self.cursor.fetchall()]

    def get_proteins(
        self,
        version,
        taxon_ids=None,
        proteome_id=None,
        go_id=None,
        pfam_id=None,
        division=None,
    ):
        """
        Flexible retrieval joining proteins, sequences, and taxa.
        Supports filtering by TaxID, Proteome ID, GO term, Pfam domain,
        and taxonomic division.

        organism and division are fetched from the taxa table via JOIN.
        """
        # Base query — taxa JOIN is always present since organism and
        # division now live there rather than on the proteins table.
        query = """
            SELECT p.accession, p.name, t.organism, t.division,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM proteins p
            JOIN sequences s ON p.seq_id = s.seq_id
            JOIN taxa t ON p.taxon_id = t.taxon_id
        """

        # Conditional JOINs for functional annotation filters
        if go_id:
            query += " JOIN protein_go pg ON p.accession = pg.accession AND p.version = pg.version"
        if pfam_id:
            query += " JOIN protein_pfam pp ON p.accession = pp.accession AND p.version = pp.version"

        # Build WHERE clause
        where_clauses = ["p.version = %s"]
        params = [version]

        # TaxID filter — accepts single int or list of ints
        if taxon_ids:
            if isinstance(taxon_ids, list):
                placeholder = ", ".join(["%s"] * len(taxon_ids))
                where_clauses.append(f"p.taxon_id IN ({placeholder})")
                params.extend(taxon_ids)
            else:
                where_clauses.append("p.taxon_id = %s")
                params.append(taxon_ids)

        # Proteome ID filter
        if proteome_id:
            where_clauses.append("p.proteome_id = %s")
            params.append(proteome_id)

        # GO term filter
        if go_id:
            where_clauses.append("pg.go_id = %s")
            params.append(go_id)

        # Pfam domain filter
        if pfam_id:
            where_clauses.append("pp.pfam_id = %s")
            params.append(pfam_id)

        # Division filter — uses the taxa table column
        if division:
            where_clauses.append("t.division = %s")
            params.append(division)

        query += " WHERE " + " AND ".join(where_clauses)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            print(f"Query Error: {err}")
            return []

    def export_fasta(self, records, version, identifier):
        """
        Exports protein records to a FASTA file.

        Header format follows standard UniProt convention:
            >ACCESSION ENTRY_NAME OS=Organism OX=TaxonID UP=ProteomeID DIV=Division
        """
        filename = f"uniprot_{identifier}_{version}.fasta"
        try:
            with open(filename, "w") as f:
                for rec in records:
                    header = (
                        f">{rec['accession']} {rec['name']} "
                        f"OS={rec['organism']} "
                        f"OX={rec['taxon_id']} "
                        f"UP={rec['proteome_id']} "
                        f"DIV={rec['division']}"
                    )
                    f.write(f"{header}\n{rec['sequence']}\n")
            print(f"Exported {len(records):,} sequences to {filename}")
        except Exception as e:
            print(f"Error exporting FASTA: {e}")

    def close(self):
        """Close database connection."""
        if self.conn and self.conn.is_connected():
            self.cursor.close()
            self.conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="UniProt Reference Set Retrieval Tool (9-Table Schema)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    parser.add_argument(
        "-version", required=True, help="UniProt version (e.g., 2026_01)"
    )

    # Retrieval filters
    parser.add_argument(
        "-taxonomy", nargs="+", type=int, help="One or more NCBI Taxonomy IDs"
    )
    parser.add_argument(
        "--proteome-id", help="Filter by Proteome ID (e.g., UP000005640)"
    )
    parser.add_argument(
        "--go-id", help="Filter by GO term (e.g., GO:0005634)"
    )
    parser.add_argument(
        "--pfam-id", help="Filter by Pfam domain (e.g., PF00870)"
    )
    parser.add_argument(
        "--division",
        choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
        help="Filter by taxonomic division",
    )

    # Utility commands
    parser.add_argument(
        "--list-versions", action="store_true", help="List all versions in the database"
    )
    parser.add_argument(
        "--list-proteomes",
        action="store_true",
        help="List all proteome IDs for this version",
    )

    args = parser.parse_args()

    DB_CONFIG = {
        "host": os.getenv("DB_HOST"),
        "user": os.getenv("DB_USER"),
        "password": os.getenv("DB_PASSWORD"),
        "database": os.getenv("DB_NAME"),
    }

    retriever = UniProtRetriever(DB_CONFIG)

    try:
        retriever.connect()

        if args.list_versions:
            retriever.list_available_versions()
            return

        if args.list_proteomes:
            proteomes = retriever.get_proteome_ids(args.version)
            print(f"\nProteomes in version {args.version}:")
            for p_id in proteomes:
                print(f"  {p_id}")
            return

        # Print query summary
        print(f"\n{'='*60}")
        print(f"UniProt Reference Set Retrieval — v5")
        print(f"Version:  {args.version}")
        if args.taxonomy:
            print(f"Taxonomy: {args.taxonomy}")
        if args.proteome_id:
            print(f"Proteome: {args.proteome_id}")
        if args.go_id:
            print(f"GO Term:  {args.go_id}")
        if args.pfam_id:
            print(f"Pfam ID:  {args.pfam_id}")
        if args.division:
            print(f"Division: {args.division}")
        print(f"{'='*60}\n")

        records = retriever.get_proteins(
            version=args.version,
            taxon_ids=args.taxonomy,
            proteome_id=args.proteome_id,
            go_id=args.go_id,
            pfam_id=args.pfam_id,
            division=args.division,
        )

        if not records:
            print("No matching records found.")
            return

        # Build output filename identifier from the most specific filter used
        identifier = "filtered_set"
        if args.proteome_id:
            identifier = args.proteome_id
        elif args.taxonomy:
            identifier = f"tax_{args.taxonomy[0]}"
        elif args.division:
            identifier = args.division.lower()

        retriever.export_fasta(records, args.version, identifier)
        print(f"Retrieved {len(records):,} sequences.")

    except Exception as e:
        print(f"\nError: {e}")
    finally:
        retriever.close()


if __name__ == "__main__":
    main()
