# DFU Resistance Analyzer

A tool for detecting antibiotic resistance genes using the CARD database, designed to rival ResFinder.

## Setup

### Windows Setup
1. **Install Conda**:
   - Download and install Miniconda or Anaconda for Windows: https://docs.conda.io/en/latest/miniconda.html
   - Open Anaconda Prompt and create an environment:
     ```
     conda create -n dfu_analyzer python=3.9
     conda activate dfu_analyzer
     ```
2. **Install BLAST+**:
   - Download BLAST+ from NCBI: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
   - Install and add the BLAST+ bin directory (e.g., `C:\Program Files\NCBI\blast-<version>\bin`) to your PATH:
     - Right-click 'This PC' → Properties → Advanced system settings → Environment Variables → Edit 'Path' → Add the bin directory.
3. **Clone the Repository**:
   ```
   git clone <repo-url>
   cd dfu-resistance-analyzer
   ```
4. **Install Dependencies**:
   ```
   pip install -r requirements.txt
   ```
   - Optional: For PDF export, install reportlab:
     ```
     pip install reportlab
     ```
5. **Set Up CARD Database**:
   - Ensure the CARD database (`card_database/`) is present with `aro_index.tsv` and BLAST database files (`card_db.*`).

### Linux/WSL Setup
1. **Install Conda**:
   - Follow the Linux instructions for Miniconda: https://docs.conda.io/en/latest/miniconda.html
   - Create an environment:
     ```
     conda create -n dfu_analyzer python=3.9
     conda activate dfu_analyzer
     ```
2. **Install BLAST+**:
   - Install BLAST+ (e.g., on Ubuntu: `sudo apt-get install ncbi-blast+`).
3. **Clone and Install**:
   ```
   git clone <repo-url>
   cd dfu-resistance-analyzer
   pip install -r requirements.txt
   ```
4. **Set Up CARD Database**:
   - Ensure the CARD database is in `card_database/`.

## Usage
1. Activate environment:
   - Windows: `conda activate dfu_analyzer` (Anaconda Prompt)
   - Linux/WSL: `conda activate dfu_analyzer`
2. Run the app:
   ```
   streamlit run app.py
   ```
3. Upload a FASTA file (e.g., `inputs/SRR25645458.fasta`).
4. Adjust thresholds and plot settings as needed.

## Outputs
- Results: `~/dfu_outputs/<sample>_report.csv` (Windows: `C:\Users\<USER>\dfu_outputs`)
- Plot: `~/dfu_outputs/<sample>_plot.html/png`

## Features
- Uses CARD database (6,439 AROs) for comprehensive ARG detection.
- Adjustable thresholds (default: 90% identity, 80% coverage).
- User-friendly plot with customizable dimensions.

