import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
import sys
import re

def load_aro_mappings(aro_file):
    """Load ARO mappings from aro_index.tsv."""
    mappings = {}
    df = pd.read_csv(aro_file, sep='\t')
    for _, row in df.iterrows():
        aro = row['ARO Accession']
        gene = row['Model Name']
        antibiotic = row.get('Drug Class', 'Unknown')
        mappings[aro] = {'gene': gene, 'antibiotic': antibiotic}
    return mappings

def run_blast(fasta_file, db_path, output_file):
    """Run BLASTn on input FASTA against CARD database."""
    cmd = [
        'blastn',
        '-query', fasta_file,
        '-db', db_path,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-num_threads', '4',
        '-max_target_seqs', '20'
    ]
    subprocess.run(cmd, check=True)

def parse_blast_results(blast_file, mappings, identity_threshold=50.0, coverage_threshold=10.0):
    """Parse BLAST results and map to resistance genes."""
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    df = pd.read_csv(blast_file, sep='\t', names=columns)
    
    print(f"Raw BLAST hits: {len(df)}")
    if not df.empty:
        print("Sample BLAST hits:")
        print(df[['sseqid', 'pident', 'length', 'evalue', 'slen']].head().to_string(index=False))
    
    # Calculate coverage based on subject length (resistance genes are short)
    df['coverage'] = (df['length'] / df['slen']) * 100
    
    # Filter by identity and coverage
    df = df[(df['pident'] >= identity_threshold) & (df['coverage'] >= coverage_threshold)]
    
    print(f"Filtered BLAST hits (>{identity_threshold}% identity, >{coverage_threshold}% coverage): {len(df)}")
    
    # Map to ARO and antibiotics
    results = []
    unmatched_aros = set()
    for _, row in df.iterrows():
        sseqid = row['sseqid']
        aro_match = re.search(r'ARO:(\d+)', sseqid)
        if aro_match:
            aro = f"ARO:{aro_match.group(1)}"
        else:
            print(f"Warning: No ARO ID found in {sseqid}")
            continue
        
        if aro in mappings:
            gene = mappings[aro]['gene']
            antibiotic = mappings[aro]['antibiotic']
            results.append({
                'Query': row['qseqid'],
                'Gene': gene,
                'Antibiotic': antibiotic,
                'Percent_Identity': row['pident'],
                'Alignment_Length': row['length'],
                'E-value': row['evalue'],
                'Bitscore': row['bitscore'],
                'Coverage': row['coverage']
            })
        else:
            unmatched_aros.add(aro)
    
    if unmatched_aros:
        print(f"Warning: {len(unmatched_aros)} AROs not found in mappings: {', '.join(sorted(unmatched_aros))}")
    
    results_df = pd.DataFrame(results)
    print(f"Mapped resistance genes: {len(results_df)}")
    return results_df

def plot_results(df, output_file):
    """Generate a bar plot of resistance genes by percent identity."""
    if df.empty:
        print("No resistance genes detected, creating empty plot")
        plt.figure(figsize=(12, 6))
        plt.title('Detected Resistance Genes')
        plt.xlabel('Gene')
        plt.ylabel('Percent Identity')
        plt.text(0.5, 0.5, 'No resistance genes detected', horizontalalignment='center', verticalalignment='center')
        plt.savefig(output_file)
        plt.close()
        return
    
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df, x='Gene', y='Percent_Identity', hue='Antibiotic')
    plt.title('Detected Resistance Genes')
    plt.xlabel('Gene')
    plt.ylabel('Percent Identity')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python dfu_resistance_analyzer.py <input_fasta>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    aro_file = 'card_database/aro_index.tsv'
    db_path = 'card_database/card_db'
    blast_output = 'outputs/blast_results.txt'
    report_output = 'outputs/resistance_report.csv'
    plot_output = 'outputs/resistance_plot.png'
    
    # Create outputs directory
    os.makedirs('outputs', exist_ok=True)
    
    # Load ARO mappings
    print(f"Loading gene mappings from {aro_file}")
    mappings = load_aro_mappings(aro_file)
    print(f"Loaded {len(mappings)} gene mappings")
    
    # Run BLAST
    print(f"Running BLAST against {db_path}")
    run_blast(fasta_file, db_path, blast_output)
    
    # Parse results
    print("Parsing BLAST results")
    results_df = parse_blast_results(blast_output, mappings)
    
    if results_df.empty:
        print("No resistance genes detected")
    else:
        print(f"Detected {len(results_df)} resistance gene hit(s)")
        print(results_df[['Gene', 'Antibiotic', 'Percent_Identity', 'Alignment_Length']].to_string(index=False))
    
    # Save report
    results_df.to_csv(report_output, index=False)
    print(f"Report saved to {report_output}")
    
    # Generate plot
    plot_results(results_df, plot_output)
    print(f"Plot saved to {plot_output}")

if __name__ == '__main__':
    main()