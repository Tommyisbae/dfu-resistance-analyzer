import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time
import urllib.error
import urllib.request
import os
import ssl

# Map resistance genes to antibiotics
GENE_TO_ANTIBIOTIC = {
    "mecA": "Methicillin",
    "blaNDM-1": "Carbapenems",
    "tetA": "Tetracycline",
    "blaKPC": "Carbapenems"
}

def parse_fasta(fasta_file):
    """
    Read FASTA file (DNA from DFU bacteria, e.g., S. aureus, P. aeruginosa).
    Why: Automates reading raw data, often manually processed in labs.
    Returns: List of sequences or empty list on error.
    """
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file '{fasta_file}' not found.")
        return []
    try:
        if os.path.getsize(fasta_file) > 1_000_000:
            print("Error: FASTA file too large (>1MB). Use smaller file for testing.")
            return []
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if not sequences:
            print(f"Error: FASTA file '{fasta_file}' is empty or invalid.")
            return []
        print(f"Loaded {len(sequences)} sequences from {fasta_file}")
        return sequences
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []

def run_blast(sequence, use_mock=True, retries=3, delay=5):
    """
    Run BLAST to find resistance genes (online or mock).
    Mock: Returns fake results for testing.
    Online: NCBI server with SSL and retry logic.
    Returns: XML file or mock list.
    """
    if use_mock:
        print("Using mock BLAST results (for testing).")
        return [
            {"Gene": "mecA", "Antibiotic": "Methicillin", "E-value": 1e-20},
            {"Gene": "blaNDM-1", "Antibiotic": "Carbapenems", "E-value": 1e-15}
        ]

    for attempt in range(retries):
        try:
            print(f"Running online BLAST query (attempt {attempt+1}/{retries})... (1-2 minutes)")
            context = ssl._create_unverified_context()
            result_handle = NCBIWWW.qblast("blastn", "nr", sequence, hitlist_size=5, ssl_context=context)
            blast_file = "blast_results.xml"
            with open(blast_file, "w") as out_handle:
                content = result_handle.read()
                out_handle.write(content)
            result_handle.close()
            print(f"BLAST results saved to {blast_file}")
            # Debug: Check if results contain expected genes
            with open(blast_file) as f:
                content = f.read()
                for gene in GENE_TO_ANTIBIOTIC:
                    if gene.lower() in content.lower():
                        print(f"Found {gene} in BLAST results")
            return blast_file
        except urllib.error.URLError as e:
            print(f"Online BLAST failed: {e}")
            if attempt < retries - 1:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print("All online retries failed.")
        except Exception as e:
            print(f"Online BLAST query failed: {e}")
            return None
    return None

def parse_blast_results(blast_input, e_value_threshold=1e-10):
    """
    Parse BLAST results (XML or mock list) to detect resistance genes.
    Returns: List of detected genes and antibiotics.
    """
    resistance_hits = []
    if isinstance(blast_input, list):
        return blast_input
    if not blast_input or not os.path.exists(blast_input):
        print("Error: No BLAST results to parse.")
        return resistance_hits
    try:
        with open(blast_input) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    print(f"Alignment title: {alignment.title}")  # Debug
                    for hsp in alignment.hsps:
                        if hsp.expect < e_value_threshold:
                            gene = None
                            for known_gene in GENE_TO_ANTIBIOTIC:
                                if known_gene.lower() in alignment.title.lower():
                                    gene = known_gene
                                    break
                            if gene:
                                resistance_hits.append({
                                    "Gene": gene,
                                    "Antibiotic": GENE_TO_ANTIBIOTIC[gene],
                                    "E-value": hsp.expect
                                })
        if not resistance_hits:
            print("No resistance genes matched in BLAST results.")
        return resistance_hits
    except Exception as e:
        print(f"Error parsing BLAST results: {e}")
        return []

def generate_report(hits, output_csv="resistance_report.csv"):
    """
    Create CSV, plot, and summary for clinicians.
    Outputs:
    - CSV: Gene, Prevalence, Antibiotic.
    - Plot: Gene frequencies.
    - Summary: Clinical advice.
    """
    if not hits:
        summary = "No resistance genes detected.\nRecommendation: Standard antibiotics may be effective, but confirm with lab tests."
        print(summary)
        with open("resistance_summary.txt", "w") as f:
            f.write(summary)
        return

    df = pd.DataFrame(hits)
    prevalence = df["Gene"].value_counts(normalize=True).reset_index()
    prevalence.columns = ["Gene", "Prevalence"]
    prevalence["Antibiotic"] = prevalence["Gene"].map(GENE_TO_ANTIBIOTIC)
    prevalence.to_csv(output_csv, index=False)
    print(f"Report saved to {output_csv}")
    print(prevalence)

    plt.figure(figsize=(8, 6))
    sns.barplot(x="Gene", y="Prevalence", data=prevalence)
    plt.title("Resistance Gene Prevalence in DFU Samples")
    plt.xlabel("Resistance Gene")
    plt.ylabel("Prevalence")
    plt.savefig("prevalence_plot.png")
    plt.close()

    summary = "Resistance Analysis for Diabetic Foot Ulcer Samples:\n"
    for _, row in prevalence.iterrows():
        summary += f"- {row['Gene']} found in {row['Prevalence']*100:.1f}% of samples, " \
                   f"indicating resistance to {row['Antibiotic']}.\n"
    summary += "Recommendation: Avoid listed antibiotics. Consider vancomycin or linezolid, and consult a microbiologist."
    print(summary)
    with open("resistance_summary.txt", "w") as f:
        f.write(summary)

def main():
    """
    Analyze DFU bacterial DNA for resistance genes.
    Why: Automates detection of genes (mecA, blaNDM-1) in local samples,
         addressing Nigeria's high resistance (92.9% multidrug-resistant DFUs).
    """
    parser = argparse.ArgumentParser(description="Analyze DFU bacterial DNA for resistance genes.")
    parser.add_argument("fasta_file", help="Path to FASTA file (e.g., test_sequences.fasta)")
    args = parser.parse_args()

    sequences = parse_fasta(args.fasta_file)
    if not sequences:
        print("No sequences to process. Exiting.")
        return

    all_hits = []
    for i, seq in enumerate(sequences[:3]):  # Process all 3 sequences
        print(f"Processing sequence {i+1}/{len(sequences)}: {seq.id}")
        blast_result = run_blast(seq.seq, use_mock=True)  # Mock results for testing
        if blast_result:
            hits = parse_blast_results(blast_result)
            all_hits.extend(hits)
        time.sleep(2)

    generate_report(all_hits)

if __name__ == "__main__":
    main()