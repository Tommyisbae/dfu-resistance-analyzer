import sys
import os
import subprocess
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("outputs/analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def validate_card_database(card_tsv, db_path):
    """Validate CARD database integrity."""
    if not os.path.exists(card_tsv):
        raise FileNotFoundError(f"CARD index file {card_tsv} not found")
    if not all(os.path.exists(f"{db_path}.{ext}") for ext in ["nhr", "nin", "nsq"]):
        raise FileNotFoundError(f"BLAST database {db_path} incomplete")
    tsv_size = os.path.getsize(card_tsv) / (1024 * 1024)  # MB
    if tsv_size < 0.5:
        raise ValueError(f"CARD index file {card_tsv} is too small ({tsv_size:.2f} MB)")
    df = pd.read_csv(card_tsv, sep='\t')
    if len(df) < 6000:
        raise ValueError(f"CARD index contains only {len(df)} AROs; expected ~6439")
    logger.info(f"CARD database validated: {len(df)} AROs, {tsv_size:.2f} MB")

def load_gene_mappings(card_tsv):
    """Load CARD gene mappings from aro_index.tsv, handling duplicate AROs."""
    df = pd.read_csv(card_tsv, sep='\t')
    if 'ARO Accession' not in df.columns or 'Model Name' not in df.columns or 'Drug Class' not in df.columns:
        raise ValueError("Invalid CARD index format: Missing required columns")
    df = df.drop_duplicates(subset='ARO Accession', keep='first')
    logger.info(f"Loaded {len(df)} unique ARO mappings from {card_tsv}")
    return df.set_index('ARO Accession')[['Model Name', 'Drug Class']].to_dict('index')

def run_blast(fasta_file, db_path, output_file):
    """Run BLASTn with error handling and timeout."""
    cmd = ["blastn", "-version"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if "2.12.0+" not in result.stdout:
            raise RuntimeError(f"Unsupported BLAST version: {result.stdout}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"BLAST version check failed: {e.stderr}")
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file {fasta_file} not found")
    cmd = [
        "blastn",
        "-query", fasta_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-num_threads", "4",
        "-max_target_seqs", "100"
    ]
    for attempt in range(3):
        try:
            logger.info(f"Running BLAST attempt {attempt + 1} for {fasta_file}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
            if os.path.exists(output_file):
                hit_count = sum(1 for _ in open(output_file))
                logger.info(f"BLAST completed: {hit_count} hits in {output_file}")
                return
            else:
                logger.error("BLAST output file not created")
                raise RuntimeError("BLAST output file not created")
        except subprocess.TimeoutExpired:
            logger.warning(f"BLAST timed out on attempt {attempt + 1}")
            if attempt == 2:
                raise TimeoutError("BLAST execution timed out after 3 attempts")
        except subprocess.CalledProcessError as e:
            logger.warning(f"BLAST attempt {attempt + 1} failed: {e.stderr}")
            if attempt == 2:
                raise RuntimeError(f"BLAST failed: {e.stderr}")

def parse_blast_results(blast_file, gene_mappings, min_identity=20.0, min_coverage=0.0):
    """Parse BLAST results and map to antibiotics."""
    if not os.path.exists(blast_file):
        raise FileNotFoundError(f"BLAST output {blast_file} not found")
    if os.path.getsize(blast_file) == 0:
        logger.warning("BLAST output is empty; no hits found")
        return pd.DataFrame()
    results = []
    total_hits = 0
    with open(blast_file, 'r') as f:
        for line in f:
            total_hits += 1
            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue
            aro = fields[1].split('|')[-2] if 'ARO' in fields[1] else fields[1]
            pident = float(fields[2])
            qlen = int(fields[12])
            slen = int(fields[13])
            coverage = (int(fields[3]) / min(qlen, slen)) * 100
            if pident >= min_identity and coverage >= min_coverage:
                if aro in gene_mappings:
                    results.append({
                        'Gene': gene_mappings[aro]['Model Name'],
                        'Antibiotic': gene_mappings[aro]['Drug Class'],
                        'Percent_Identity': pident,
                        'Coverage': coverage,
                        'Evalue': float(fields[10])
                    })
    df = pd.DataFrame(results)
    logger.info(f"Parsed {total_hits} BLAST hits, found {len(df)} ARGs with {min_identity}% identity, {min_coverage}% coverage")
    return df

def plot_results(df, output_file):
    """Generate Plotly bar plot for top 10 ARGs."""
    if df.empty:
        raise ValueError("No resistance genes to plot")
    top_df = df.nlargest(10, 'Percent_Identity')
    fig = px.bar(
        top_df,
        x="Gene",
        y="Percent_Identity",
        color="Antibiotic",
        title="Top 10 Resistance Genes",
        height=500,
        text="Percent_Identity",
        hover_data=["Coverage", "Evalue"]
    )
    fig.update_traces(texttemplate="%{text:.1f}%", textposition="auto")
    fig.update_layout(xaxis_tickangle=45)
    html_file = output_file.replace(".png", ".html")
    fig.write_html(html_file)
    logger.info(f"Generated HTML plot: {html_file}")
    try:
        fig.write_image(output_file, format="png")
        logger.info(f"Generated PNG plot: {output_file}")
    except Exception as e:
        logger.warning(f"Failed to generate PNG plot: {str(e)}. Using HTML plot only.")
        if not os.path.exists(html_file):
            raise RuntimeError("Failed to generate both PNG and HTML plots")

def save_results(df, output_file):
    """Save results to CSV and summary table."""
    if df.empty:
        raise ValueError("No resistance genes to save")
    df.to_csv(output_file, index=False)
    summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
    summary.to_csv(output_file.replace(".csv", "_summary.csv"), index=False)
    logger.info(f"Saved results: {output_file}, summary: {output_file.replace('.csv', '_summary.csv')}")

def validate_fasta(fasta_file):
    """Validate FASTA file format and content."""
    logger.info(f"Checking FASTA file: {fasta_file}, path: {os.path.abspath(fasta_file)}")
    try:
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA file {fasta_file} not found")
        if not os.access(fasta_file, os.R_OK):
            raise PermissionError(f"FASTA file {fasta_file} not readable")
        with open(fasta_file, 'r') as f:
            records = list(SeqIO.parse(f, "fasta"))
        if not records:
            raise ValueError("FASTA file is empty or invalid")
        total_length = sum(len(record.seq) for record in records)
        if total_length < 1000:
            raise ValueError(f"FASTA file too small: {total_length} bp")
        base_name = os.path.basename(fasta_file)
        if base_name.startswith("GCF_") and not os.path.exists(fasta_file.replace("GCF_", "GCA_")):
            logger.warning(f"Checking for GCA_ prefix: {fasta_file.replace('GCF_', 'GCA_')}")
        logger.info(f"Validated FASTA: {fasta_file}, {len(records)} sequences, {total_length} bp")
        return True
    except Exception as e:
        logger.error(f"Invalid FASTA: {str(e)}")
        raise ValueError(f"Invalid FASTA: {str(e)}")

def main(fasta_file=None):
    """Main function to run ARG detection."""
    if fasta_file is None:
        if len(sys.argv) < 2:
            raise ValueError("No FASTA file provided. Usage: python dfu_resistance_analyzer.py <fasta_file>")
        fasta_file = sys.argv[1]
    
    # Validate inputs
    validate_fasta(fasta_file)
    card_tsv = "card_database/aro_index.tsv"
    db_path = "card_database/card_db"
    validate_card_database(card_tsv, db_path)
    
    # Paths
    output_dir = "outputs"
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    blast_output = f"{output_dir}/{base_name}_blast_results.txt"
    csv_output = f"{output_dir}/{base_name}_report.csv"
    plot_output = f"{output_dir}/{base_name}_plot.png"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load CARD mappings
    gene_mappings = load_gene_mappings(card_tsv)
    
    # Run BLAST
    run_blast(fasta_file, db_path, blast_output)
    
    # Parse results
    df = parse_blast_results(blast_output, gene_mappings, min_identity=20.0, min_coverage=0.0)
    
    # Save results
    if not df.empty:
        save_results(df, csv_output)
        plot_results(df, plot_output)
    else:
        logger.warning("No resistance genes detected with thresholds (20% identity, 0% coverage).")
        raise ValueError("No resistance genes detected with thresholds (20% identity, 0% coverage). Check FASTA file, CARD database, or BLAST output in outputs/analysis.log.")

if __name__ == "__main__":
    main()