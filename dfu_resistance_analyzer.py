import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time
import urllib.error
import os

# Map resistance genes to antibiotics (known genes, but need local confirmation)
GENE_TO_ANTIBIOTIC = {
    "mecA": "Methicillin",
    "blaNDM-1": "Carbapenems",
    "tetA": "Tetracycline",
    "blaKPC": "Carbapenems"
}

def parse_fasta(fasta_file):
    """
    Read FASTA file (DNA from DFU bacteria, e.g., S. aureus, P. aeruginosa).
    FASTA: '>ID\nATGAAATATG...' (DNA sequences from sequencing).
    Why: Automates reading raw data, which is often manually processed in labs.
    Source: Public (NCBI SRA, e.g., PRJNA644678) or professor's data.
    Returns: List of sequences or empty list on error.
    """
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file '{fasta_file}' not found.")
        return []
    try:
        # Check file size (limit to 1MB for MVP to avoid slowdown)
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

def run_blast(sequence):
    """
    Run BLAST to find resistance genes (e.g., mecA) in NCBI database.
    Why: Automates gene detection, critical since known genes vary by region.
    BLAST: Compares DNA to known genes, like searching a library.
    Returns: XML file with results or None on error.
    """
    try:
        print("Running BLAST query... (1-2 minutes per sequence)")
        result_handle = NCBIWWW.qblast("blastn", "nr", sequence, hitlist_size=5)
        blast_file = "blast_results.xml"
        with open(blast_file, "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        return blast_file
    except urllib.error.URLError:
        print("BLAST failed: No internet connection.")
        return None
    except Exception as e:
        print(f"BLAST query failed: {e}")
        return None

def parse_blast_results(blast_file, e_value_threshold=1e-10):
    """
    Parse BLAST results to detect resistance genes.
    Why: Confirms known genes (e.g., mecA) in local samples, addressing regional gaps.
    Filters: E-value < 1e-10 for accuracy.
    Returns: List of detected genes and antibiotics.
    """
    resistance_hits = []
    if not blast_file or not os.path.exists(blast_file):
        print("Error: No BLAST results to parse.")
        return resistance_hits
    try:
        with open(blast_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
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
        return resistance_hits
    except Exception as e:
        print(f"Error parsing BLAST results: {e}")
        return []

def generate_report(hits, output_csv="resistance_report.csv"):
    """
    Create CSV, plot, and summary for clinicians.
    Why: Translates known genes into local, actionable advice (e.g., avoid Methicillin).
    Outputs:
    - CSV: Gene, Prevalence, Antibiotic.
    - Plot: Gene frequencies.
    - Summary: Clinical advice.
    Example Run Log:
        Loaded 3 sequences from test_sequences.fasta
        Processing sequence 1/3: Sequence1_S_aureus_mecA
        Running BLAST query... (1-2 minutes)
        Processing sequence 2/3: Sequence2_P_aeruginosa_blaNDM-1
        ...
        Report saved to resistance_report.csv
           Gene  Prevalence Antibiotic
        0  mecA        0.67 Methicillin
        1  blaNDM-1    0.33 Carbapenems
        Resistance Analysis for Diabetic Foot Ulcer Samples:
        - mecA found in 66.7% of samples, indicating resistance to Methicillin.
        - blaNDM-1 found in 33.3% of samples, indicating resistance to Carbapenems.
        Recommendation: Avoid listed antibiotics. Consider vancomycin or linezolid.
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
    Why: Automates detection of known genes (mecA, blaNDM-1) in local samples,
         addressing Nigeria's high resistance (92.9% multidrug-resistant DFUs).
    Steps:
    1. Read FASTA (e.g., NCBI SRA PRJNA644678 or test_sequences.fasta).
    2. Run BLAST to find genes.
    3. Generate clinical reports.
    """
    parser = argparse.ArgumentParser(description="Analyze DFU bacterial DNA for resistance genes.")
    parser.add_argument("fasta_file", help="Path to FASTA file (e.g., test_sequences.fasta)")
    args = parser.parse_args()

    sequences = parse_fasta(args.fasta_file)
    if not sequences:
        print("No sequences to process. Exiting.")
        return

    all_hits = []
    for i, seq in enumerate(sequences[:3]):  # Limit to 3 for testing
        print(f"Processing sequence {i+1}/{min(3, len(sequences))}: {seq.id}")
        blast_file = run_blast(seq.seq)
        if blast_file:
            hits = parse_blast_results(blast_file)
            all_hits.extend(hits)
        time.sleep(2)  # Avoid overwhelming NCBI server

    generate_report(all_hits)

if __name__ == "__main__":
    main()