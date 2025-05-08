import sys
import os
import subprocess
import pandas as pd
import plotly.express as px
from Bio import SeqIO

def load_gene_mappings(card_tsv):
    """Load CARD gene mappings from aro_index.tsv, handling duplicate AROs."""
    if not os.path.exists(card_tsv):
        raise FileNotFoundError(f"CARD index file {card_tsv} not found")
    df = pd.read_csv(card_tsv, sep='\t')
    if 'ARO Accession' not in df.columns or 'Model Name' not in df.columns or 'Drug Class' not in df.columns:
        raise ValueError("Invalid CARD index format: Missing required columns")
    # Deduplicate by keeping first occurrence of each ARO
    df = df.drop_duplicates(subset='ARO Accession', keep='first')
    return df.set_index('ARO Accession')[['Model Name', 'Drug Class']].to_dict('index')

def run_blast(fasta_file, db_path, output_file):
    """Run BLASTn with error handling."""
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file {fasta_file} not found")
    if not os.path.exists(db_path + ".nhr"):
        raise FileNotFoundError(f"BLAST database {db_path} not found")
    cmd = [
        "blastn",
        "-query", fasta_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-num_threads", "4",
        "-max_target_seqs", "20"
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"BLAST failed: {e.stderr}")

def parse_blast_results(blast_file, gene_mappings):
    """Parse BLAST results and map to antibiotics."""
    if not os.path.exists(blast_file):
        raise FileNotFoundError(f"BLAST output {blast_file} not found")
    results = []
    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue
            aro = fields[1]  # sseqid (ARO accession)
            pident = float(fields[2])
            qlen = int(fields[12])
            slen = int(fields[13])
            coverage = (int(fields[3]) / min(qlen, slen)) * 100
            if pident >= 50 and coverage >= 10:
                if aro in gene_mappings:
                    results.append({
                        'Gene': gene_mappings[aro]['Model Name'],
                        'Antibiotic': gene_mappings[aro]['Drug Class'],
                        'Percent_Identity': pident,
                        'Coverage': coverage,
                        'Evalue': float(fields[10])
                    })
    return pd.DataFrame(results)

def plot_results(df, output_file):
    """Generate Plotly bar plot for top 10 ARGs."""
    if df.empty:
        raise ValueError("No resistance genes to plot")
    top_df = df.nlargest(10, 'Percent_Identity')
    fig = px.bar(top_df, x="Gene", y="Percent_Identity", color="Antibiotic",
                 title="Top 10 Resistance Genes", height=500)
    fig.update_layout(xaxis_tickangle=45)
    fig.write_html(output_file.replace(".png", ".html"))
    fig.write_image(output_file, format="png")

def save_results(df, output_file):
    """Save results to CSV and summary table."""
    if df.empty:
        raise ValueError("No resistance genes to save")
    df.to_csv(output_file, index=False)
    summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
    summary.to_csv(output_file.replace(".csv", "_summary.csv"), index=False)

def validate_fasta(fasta_file):
    """Validate FASTA file format."""
    try:
        with open(fasta_file, 'r') as f:
            records = list(SeqIO.parse(f, "fasta"))
        if not records:
            raise ValueError("FASTA file is empty or invalid")
        return True
    except Exception as e:
        raise ValueError(f"Invalid FASTA format: {str(e)}")

def main(fasta_file=None):
    """Main function to run ARG detection."""
    if fasta_file is None:
        if len(sys.argv) < 2:
            raise ValueError("No FASTA file provided. Usage: python dfu_resistance_analyzer.py <fasta_file>")
        fasta_file = sys.argv[1]
    
    # Validate inputs
    validate_fasta(fasta_file)
    
    # Paths
    card_tsv = "card_database/aro_index.tsv"
    db_path = "card_database/card_db"
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
    df = parse_blast_results(blast_output, gene_mappings)
    
    # Save results
    if not df.empty:
        save_results(df, csv_output)
        plot_results(df, plot_output)
    else:
        raise ValueError("No resistance genes detected with given thresholds (50% identity, 10% coverage).")

if __name__ == "__main__":
    main()