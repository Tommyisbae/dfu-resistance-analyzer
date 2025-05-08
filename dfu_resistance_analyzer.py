import sys
import os
import subprocess
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import logging
import re
import shutil
import tempfile
import psutil

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

def sanitize_filename(filename):
    """Sanitize file name to remove invalid characters."""
    return re.sub(r'[^a-zA-Z0-9._-]', '_', os.path.basename(filename))

def validate_card_database(card_tsv, db_path):
    """Validate CARD database integrity."""
    card_tsv = os.path.abspath(card_tsv)
    db_path = os.path.abspath(db_path)
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
    logger.info(f"CARD database validated: {len(df)} AROs, {tsv_size:.2f} MB at {card_tsv}")

def load_gene_mappings(card_tsv):
    """Load CARD gene mappings from aro_index.tsv, handling duplicate AROs."""
    card_tsv = os.path.abspath(card_tsv)
    df = pd.read_csv(card_tsv, sep='\t')
    if 'ARO Accession' not in df.columns or 'Model Name' not in df.columns or 'Drug Class' not in df.columns:
        raise ValueError("Invalid CARD index format: Missing required columns")
    df = df.drop_duplicates(subset='ARO Accession', keep='first')
    logger.info(f"Loaded {len(df)} unique ARO mappings from {card_tsv}")
    return df.set_index('ARO Accession')[['Model Name', 'Drug Class']].to_dict('index')

def check_environment():
    """Check system environment for BLAST execution."""
    blast_path = shutil.which("blastn")
    if not blast_path:
        logger.error("blastn executable not found in PATH")
        raise RuntimeError("blastn executable not found in PATH")
    
    mem = psutil.virtual_memory()
    if mem.available < 2 * 1024 * 1024 * 1024:  # Less than 2GB
        logger.error(f"Insufficient memory: {mem.available / (1024 * 1024)} MB available")
        raise RuntimeError("Insufficient memory. Free up at least 2GB.")
    
    logger.info(f"Environment: blastn at {blast_path}, {mem.available / (1024 * 1024)} MB free, {os.cpu_count()} CPUs")
    return blast_path

def run_blast(fasta_file, db_path, output_file):
    """Run BLASTn with error handling and fallback."""
    blast_path = check_environment()
    fasta_file = os.path.abspath(fasta_file)
    db_path = os.path.abspath(db_path)
    output_file = os.path.abspath(output_file)
    
    cmd = [blast_path, "-version"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if "2.12.0+" not in result.stdout:
            raise RuntimeError(f"Unsupported BLAST version: {result.stdout}")
        logger.info(f"BLAST version: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST version check failed: {e.stderr}")
        raise RuntimeError(f"BLAST version check failed: {e.stderr}")
    
    if not os.path.exists(fasta_file):
        logger.error(f"FASTA file {fasta_file} not found")
        raise FileNotFoundError(f"FASTA file {fasta_file} not found")
    
    if not os.access(fasta_file, os.R_OK):
        logger.error(f"No read permission for FASTA file {fasta_file}")
        raise PermissionError(f"No read permission for FASTA file {fasta_file}")
    
    if not all(os.path.exists(f"{db_path}.{ext}") for ext in ["nhr", "nin", "nsq"]):
        logger.error(f"BLAST database {db_path} not found")
        raise FileNotFoundError(f"BLAST database {db_path} not found")
    
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    if not os.access(output_dir, os.W_OK):
        logger.error(f"No write permission for output directory {output_dir}")
        raise PermissionError(f"No write permission for output directory {output_dir}")
    
    # Primary BLAST command
    cmd = [
        blast_path,
        "-query", fasta_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-num_threads", "4",
        "-max_target_seqs", "100"
    ]
    
    for attempt in range(3):
        try:
            logger.info(f"Running BLAST attempt {attempt + 1} for {fasta_file}: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
            logger.info(f"BLAST stdout: {result.stdout[:1000] if result.stdout else 'No stdout'}...")
            logger.info(f"BLAST stderr: {result.stderr[:1000] if result.stderr else 'No stderr'}...")
            if os.path.exists(output_file):
                hit_count = sum(1 for _ in open(output_file))
                logger.info(f"BLAST completed: {hit_count} hits in {output_file}")
                # Log first few lines of output
                with open(output_file, "r") as f:
                    logger.info(f"BLAST output (first 500 chars): {f.read(500)}...")
                return
            else:
                logger.error(f"BLAST output file {output_file} not created")
                raise RuntimeError(f"BLAST output file {output_file} not created")
        except subprocess.TimeoutExpired:
            logger.warning(f"BLAST timed out on attempt {attempt + 1}")
            if attempt == 2:
                raise TimeoutError("BLAST execution timed out after 3 attempts")
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST attempt {attempt + 1} failed: {e.stderr}")
            if attempt == 2:
                # Try fallback command
                logger.info("Trying fallback BLAST command")
                cmd = [
                    blast_path,
                    "-query", fasta_file,
                    "-db", db_path,
                    "-out", output_file,
                    "-outfmt", "6"
                ]
                try:
                    logger.info(f"Running fallback BLAST: {' '.join(cmd)}")
                    result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
                    logger.info(f"Fallback BLAST stdout: {result.stdout[:1000] if result.stdout else 'No stdout'}...")
                    logger.info(f"Fallback BLAST stderr: {result.stderr[:1000] if result.stderr else 'No stderr'}...")
                    if os.path.exists(output_file):
                        hit_count = sum(1 for _ in open(output_file))
                        logger.info(f"Fallback BLAST completed: {hit_count} hits in {output_file}")
                        with open(output_file, "r") as f:
                            logger.info(f"Fallback BLAST output (first 500 chars): {f.read(500)}...")
                        return
                    else:
                        raise RuntimeError(f"Fallback BLAST output file {output_file} not created")
                except Exception as e:
                    logger.error(f"Fallback BLAST failed: {str(e)}")
                    raise RuntimeError(f"Fallback BLAST failed: {str(e)}")
        except Exception as e:
            logger.error(f"BLAST attempt {attempt + 1} failed with unexpected error: {str(e)}")
            if attempt == 2:
                raise RuntimeError(f"BLAST failed: {str(e)}")

def parse_blast_results(blast_file, gene_mappings, min_identity=20.0, min_coverage=0.0):
    """Parse BLAST results and map to antibiotics."""
    blast_file = os.path.abspath(blast_file)
    if not os.path.exists(blast_file):
        # Check for temp_ prefixed file as fallback
        temp_blast_file = os.path.join(os.path.dirname(blast_file), f"temp_{os.path.basename(blast_file)}")
        if os.path.exists(temp_blast_file):
            logger.warning(f"BLAST output {blast_file} not found, using {temp_blast_file}")
            blast_file = temp_blast_file
        else:
            logger.error(f"BLAST output {blast_file} not found")
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
    output_file = os.path.abspath(output_file)
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
    try:
        fig.write_html(html_file)
        logger.info(f"Generated HTML plot: {html_file}")
    except Exception as e:
        logger.error(f"Failed to generate HTML plot: {str(e)}")
        raise RuntimeError(f"Failed to generate HTML plot: {str(e)}")
    try:
        fig.write_image(output_file, format="png")
        logger.info(f"Generated PNG plot: {output_file}")
    except Exception as e:
        logger.warning(f"Failed to generate PNG plot: {str(e)}. HTML plot available at {html_file}")

def save_results(df, output_file):
    """Save results to CSV and summary table."""
    output_file = os.path.abspath(output_file)
    if df.empty:
        raise ValueError("No resistance genes to save")
    df.to_csv(output_file, index=False)
    summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
    summary.to_csv(output_file.replace(".csv", "_summary.csv"), index=False)
    logger.info(f"Saved results: {output_file}, summary: {output_file.replace('.csv', '_summary.csv')}")

def validate_fasta(fasta_file):
    """Validate FASTA file format and content."""
    fasta_file = os.path.abspath(fasta_file)
    logger.info(f"Checking FASTA file: {fasta_file}")
    try:
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA file {fasta_file} not found")
        if not os.access(fasta_file, os.R_OK):
            raise PermissionError(f"FASTA file {fasta_file} not readable")
        file_size = os.path.getsize(fasta_file) / (1024 * 1024)  # MB
        logger.info(f"FASTA file size: {file_size:.2f} MB")
        with open(fasta_file, 'r') as f:
            records = list(SeqIO.parse(f, "fasta"))
        if not records:
            raise ValueError("FASTA file is empty or invalid")
        total_length = sum(len(record.seq) for record in records)
        if total_length < 1000:
            raise ValueError(f"FASTA file too small: {total_length} bp")
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
    
    # Check environment
    check_environment()
    
    # Validate inputs
    validate_fasta(fasta_file)
    card_tsv = os.path.abspath("card_database/aro_index.tsv")
    db_path = os.path.abspath("card_database/card_db")
    validate_card_database(card_tsv, db_path)
    
    # Sanitize file name
    base_name = sanitize_filename(os.path.splitext(os.path.basename(fasta_file))[0])
    
    # Paths
    output_dir = os.path.abspath("outputs")
    blast_output = os.path.join(output_dir, f"{base_name}_blast_results.txt")
    csv_output = os.path.join(output_dir, f"{base_name}_report.csv")
    plot_output = os.path.join(output_dir, f"{base_name}_plot.png")
    
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