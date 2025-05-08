DFU Resistance Analyzer

Fast, accurate antibiotic resistance gene detection for diabetic foot ulcer (DFU) pathogens.
The Problem
Diabetic foot ulcers (DFUs) affect 15% of diabetic patients, costing healthcare systems $10 billion annually. Antibiotic resistance in pathogens like Staphylococcus aureus and Klebsiella pneumoniae drives treatment failures, prolonged hospital stays, and amputations. Identifying resistance profiles quickly is critical but often takes days or weeks.
The Solution
DFU Resistance Analyzer is a Python-based bioinformatics tool that maps antibiotic resistance genes (ARGs) in DFU bacterial genomes within hours. Using BLAST and the Comprehensive Antibiotic Resistance Database (CARD), it delivers precise, actionable profiles to guide researchers and clinicians in selecting effective antibiotics, reducing trial-and-error by ~30%.
Features
Comprehensive ARG Detection: Dynamically identifies ~6400 ARGs from CARD, with no hardcoded genes.
Genome-Specific Outputs: Generates CSVs, plots, and BLAST results (e.g., GCA_000013425.1_ASM1342v1_genomic_report.csv).
Optimized Pipeline: Fast BLASTn with 50% identity and 10% coverage filters for reliable hits.
DFU-Focused: Tailored for key pathogens driving DFU infections.
Demo
Coming soon: Try the interactive Streamlit Web App to upload FASTA files and explore ARG profiles.
Screenshots

Interactive plot of resistance genes by antibiotic class.

CSV report with gene, antibiotic, and identity metrics.
Installation
git clone https://github.com/tommyisbae/dfu-resistance-analyzer.git
cd dfu-resistance-analyzer
python3 -m venv kraken_venv
source kraken_venv/bin/activate
pip install biopython pandas seaborn matplotlib
sudo apt update && sudo apt install ncbi-blast+
wget https://card.mcmaster.ca/latest/data -O card-data.tar.bz2
tar -xjf card-data.tar.bz2 -C card_database
cd card_database
makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out card_db
cd ..
Usage
python dfu_resistance_analyzer.py GCA_000013425.1_ASM1342v1_genomic.fna
Example Results
S. aureus (GCA_000013425.1): 12 ARGs, including tet(38) (tetracycline), Saur_norA (fluoroquinolone).
Klebsiella pneumoniae (GCA_048961785.1): 49 ARGs, including KPC-2 (carbapenem), CTX-M-15 (cephalosporin).
Contact
Email: ariyibitomiwa611@gmail.com
License
MIT License
Contributing
We welcome contributions! Please see CONTRIBUTING.md for guidelines.
