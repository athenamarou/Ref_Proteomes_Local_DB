import subprocess
import os
from ete3 import Tree
from get_reference_set_lib import UniProtRetriever, get_db_config

TARGET_VERSION = "2026_01"
TARGET_TAXON = 6087
TARGET_DOMAIN = "PF00858"

print(f"1. Fetching {TARGET_DOMAIN} in Taxon {TARGET_TAXON} from database...")

with UniProtRetriever(get_db_config()) as db:
    records = db.get_proteins_by_hmm_hit(
        version=TARGET_VERSION,
        hmm_query=TARGET_DOMAIN,
        taxon_ids=TARGET_TAXON,
        evalue_cutoff=1e-5
    )

    if not records:
        print("No proteins found.")
        exit()

    # Save the raw unaligned FASTA
    raw_fasta = db.export_fasta(
        records=records, 
        version=TARGET_VERSION, 
        identifier=f"hmm_{TARGET_DOMAIN}_tax_{TARGET_TAXON}",
        output_dir="./tree_data"
    )

# Set up file paths for the next steps
aligned_fasta = raw_fasta.replace(".fasta", "_aligned.fasta")
tree_file = raw_fasta.replace(".fasta", ".nwk") # .nwk is the standard Newick tree format

print(f"\n2. Aligning sequences with MAFFT...")
# Runs: mafft --auto raw_fasta > aligned_fasta
with open(aligned_fasta, "w") as out_f:
    subprocess.run(["mafft", "--auto", raw_fasta], stdout=out_f, stderr=subprocess.DEVNULL)

print(f"3. Building phylogenetic tree with FastTree...")
# Runs: FastTree aligned_fasta > tree_file
with open(tree_file, "w") as out_f:
    subprocess.run(["FastTree", aligned_fasta], stdout=out_f, stderr=subprocess.DEVNULL)

print(f"\n✓ Pipeline Complete!")
print(f"  Raw FASTA: {raw_fasta}")
print(f"  Alignment: {aligned_fasta}")
print(f"  Tree File: {tree_file}")


print("\n4. Drawing tree with ETE3...")

os.environ['QT_QPA_PLATFORM'] = 'offscreen'
if os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
    with open(tree_file, 'r') as f:
        newick_text = f.read().strip()
    
    # Try loading the string directly
    try:
        t = Tree(newick_text, format=1)
        output_pdf = raw_fasta.replace(".fasta", ".pdf")
        t.render(output_pdf)
        print(f"✓ Tree rendered to {output_pdf}")
    except Exception as e:
        print(f"✗ ETE3 still struggling: {e}")
else:
    print("✗ Tree file is empty or missing. Check your alignment!")