import sys
import os
import subprocess
import pandas as pd
import logging
from dfu_resistance_analyzer import load_gene_mappings, parse_blast_results

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("outputs/debug.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def debug_blast(fasta_file, db_path, output_file, card_tsv):
    """Run BLAST and parse ARGs."""
    possible_paths = [
        fasta_file,
        os.path.join(os.getcwd(), fasta_file),
        os.path.join(os.getcwd(), "data", fasta_file)
    ]
    fasta_path = None
    for path in possible_paths:
        if os.path.exists(path):
            fasta_path = path
            break
    if not fasta_path:
        logger.error(f"FASTA file {fasta_file} not found in {possible_paths}")
        return
    if not all(os.path.exists(f"{db_path}.{ext}") for ext in ["nhr", "nin", "nsq"]):
        logger.error(f"BLAST database {db_path} incomplete")
        return
    cmd = [
        "blastn",
        "-query", fasta_path,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-num_threads", "4",
        "-max_target_seqs", "100"
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if os.path.exists(output_file):
            hit_count = sum(1 for _ in open(output_file))
            logger.info(f"BLAST completed: {hit_count} hits in {output_file}")
            with open(output_file, 'r') as f:
                logger.info(f"First 5 hits:\n{f.read(500)}...")
            # Parse ARGs
            gene_mappings = load_gene_mappings(card_tsv)
            df = parse_blast_results(output_file, gene_mappings, min_identity=20.0, min_coverage=0.0)
            logger.info(f"Detected {len(df)} ARGs")
            if not df.empty:
                logger.info(f"Sample ARGs:\n{df[['Gene', 'Antibiotic', 'Percent_Identity']].head().to_string()}")
            df.to_csv(output_file.replace(".txt", "_args.csv"), index=False)
        else:
            logger.error("BLAST output file not created")
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST failed: {e.stderr}")

if __name__ == "__main__":
    fasta_file = sys.argv[1] if len(sys.argv) > 1 else "GCA_048961785.1_ASM4896178v1_genomic.fna"
    db_path = "card_database/card_db"
    card_tsv = "card_database/aro_index.tsv"
    output_file = "outputs/debug_blast_results.txt"
    os.makedirs("outputs", exist_ok=True)
    debug_blast(fasta_file, db_path, output_file, card_tsv)