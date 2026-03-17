"""
UniProt Reference Set Retrieval Tool
Retrieves specific versions of UniProt reference proteomes from local database
"""

import mysql.connector
import argparse
import os
from dotenv import load_dotenv
import sys

load_dotenv()


class UniProtRetriever:
    """
    Retrieves UniProt reference sets from the local database
    supporting TaxID, Proteome ID, Pfam, and GO filters.
    """

    def __init__(self, config):
        self.config = config
        self.conn = None
        self.cursor = None

    def connect(self):
        """Establish database connection"""
        try:
            self.conn = mysql.connector.connect(**self.config)
            self.cursor = self.conn.cursor(dictionary=True)
            print("Connected to database successfully")
        except mysql.connector.Error as err:
            print(f"Database connection error: {err}")
            sys.exit(1)

    def list_available_versions(self):
        """List all available versions in the database"""
        query = """
            SELECT version, COUNT(*) as protein_count,
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
        """Returns all unique Proteome IDs for a specific version"""
        query = "SELECT DISTINCT proteome_id FROM proteins WHERE version = %s"
        self.cursor.execute(query, (version,))
        return [row["proteome_id"] for row in self.cursor.fetchall()]

    def get_proteins(
        self, version, taxon_ids=None, proteome_id=None, go_id=None, pfam_id=None
    ):
        """
        Flexible retrieval joining metadata and sequences.
        Supports querying by TaxID, Proteome, GO ID, or Pfam ID.
        """
        # Base query joining proteins and sequences to get the strings
        query = """
            SELECT p.accession, p.name, p.organism, p.taxon_id, p.proteome_id, s.sequence
            FROM proteins p
            JOIN sequences s ON p.seq_id = s.seq_id
        """

        # Join functional tables if filters are provided
        if go_id:
            query += " JOIN protein_go pg ON p.accession = pg.accession AND p.version = pg.version"
        if pfam_id:
            query += " JOIN protein_pfam pp ON p.accession = pp.accession AND p.version = pp.version"

        # Initialize filter clauses and parameters
        where_clauses = ["p.version = %s"]
        params = [version]

        # Apply TaxID filter
        if taxon_ids:
            if isinstance(taxon_ids, list):
                placeholder = ", ".join(["%s"] * len(taxon_ids))
                where_clauses.append(f"p.taxon_id IN ({placeholder})")
                params.extend(taxon_ids)
            else:
                where_clauses.append("p.taxon_id = %s")
                params.append(taxon_ids)

        # Apply Proteome ID filter
        if proteome_id:
            where_clauses.append("p.proteome_id = %s")
            params.append(proteome_id)

        # Apply Functional filters
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
            print(f"Query Error: {err}")
            return []

    def export_fasta(self, records, version, identifier, filtered=False):
        """Exports protein records to a FASTA file"""
        filename = f"uniprot_{identifier}_{version}.fasta"
        try:
            with open(filename, "w") as f:
                for rec in records:
                    #header = f">{rec['accession']} {rec['name']} OX={rec['taxon_id']} UP={rec['proteome_id']}"
                    header = f">{rec['taxon_id']}.{rec['accession']}"
                    f.write(f"{header}\n{rec['sequence']}\n")
            print(f"✓ Exported {len(records):,} sequences to {filename}")
        except Exception as e:
            print(f"Error exporting FASTA: {e}")

    def close(self):
        """Close database connection"""
        if self.conn and self.conn.is_connected():
            self.cursor.close()
            self.conn.close()


def import_datetime():
    """Helper to get current timestamp"""
    import datetime

    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main():
    parser = argparse.ArgumentParser(
        description="UniProt Reference Set Retrieval Tool (8-Table Schema Version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required Version
    parser.add_argument(
        "-version", required=True, help="UniProt version (e.g., 2026_01)"
    )

    # Retrieval Filters (Your Wishes)
    parser.add_argument(
        "-taxonomy", nargs="+", type=int, help="One or more Taxonomy IDs"
    )
    parser.add_argument(
        "--proteome-id", help="Filter by specific Proteome ID (e.g., UP000005640)"
    )
    parser.add_argument("--go-id", help="Filter by GO ID (e.g., GO:0005634)")
    parser.add_argument("--pfam-id", help="Filter by Pfam ID (e.g., PF00870)")

    # Utility Commands
    parser.add_argument(
        "--list-versions", action="store_true", help="List all available versions"
    )
    parser.add_argument(
        "--list-proteomes",
        action="store_true",
        help="List all proteome IDs in this version",
    )
    parser.add_argument(
        "--filtered-headers", action="store_true", help="Use simplified FASTA headers"
    )

    args = parser.parse_args()

    # Load DB Config
    DB_CONFIG = {
        "host": os.getenv("DB_HOST"),
        "user": os.getenv("DB_USER"),
        "password": os.getenv("DB_PASSWORD"),
        "database": os.getenv("DB_NAME"),
    }

    retriever = UniProtRetriever(DB_CONFIG)

    try:
        retriever.connect()

        # Handle Utility Commands
        if args.list_versions:
            retriever.list_available_versions()
            return

        if args.list_proteomes:
            proteomes = retriever.get_proteome_ids(args.version)
            print(f"\nProteomes in version {args.version}:")
            for p_id in proteomes:
                print(f" - {p_id}")
            return

        # Execute Protein Retrieval
        print(f"\n{'='*60}")
        print(f"UniProt Reference Set Retrieval")
        print(f"Version: {args.version}")
        if args.taxonomy:
            print(f"Taxonomy: {args.taxonomy}")
        if args.proteome_id:
            print(f"Proteome: {args.proteome_id}")
        if args.go_id:
            print(f"GO Term:  {args.go_id}")
        if args.pfam_id:
            print(f"Pfam ID:  {args.pfam_id}")
        print(f"{'='*60}\n")

        records = retriever.get_proteins(
            version=args.version,
            taxon_ids=args.taxonomy,
            proteome_id=args.proteome_id,
            go_id=args.go_id,
            pfam_id=args.pfam_id,
        )

        if not records:
            print("\n✗ No matching data found.")
            return

        # Define an identifier for the output filename
        identifier = "filtered_set"
        if args.proteome_id:
            identifier = args.proteome_id
        elif args.taxonomy:
            identifier = f"tax_{args.taxonomy[0]}"

        # Export to FASTA
        retriever.export_fasta(records, args.version, identifier)

        print(f"\nSuccessfully retrieved {len(records):,} sequences.")

    except Exception as e:
        print(f"\nError: {e}")
    finally:
        retriever.close()


if __name__ == "__main__":
    main()
