# DFU Resistance Analyzer

A web-based tool for detecting antibiotic resistance genes using the CARD database.

## Web Usage
- Visit the deployed app at [https://dfu-resistance-analyzer.onrender.com](https://dfu-resistance-analyzer.onrender.com).
- Upload a FASTA file (≤100MB) to detect resistance genes.
- Select a threshold preset (e.g., "Strict (90%, 80%)") and adjust plot settings.
- Download results (CSV, HTML/PNG plots, optional PDF).

## For Files >100MB
- Convert .tgz to FASTA (e.g., `tar -xzf file.tgz && seqtk seq -a file.fastq > file.fasta`).
- Split into 100MB chunks (e.g., `seqkit split -s 100000000 file.fasta`).
- Upload chunks separately and combine results manually.

## Secure Handling
- Uploaded files are sensitive and handled securely. They are stored temporarily in memory, processed, and deleted immediately after results are generated.
- No sensitive data is logged.

## Deployment Notes (Render)
- Hosted on Render using Docker.
- Requirements: Internet access and a modern web browser.
- Source code: Available at [GitHub Repository URL](https://github.com/yourusername/dfu-resistance-analyzer) (replace with your repo URL).
- Built with Python 3.9, BLAST+, and the CARD database.

## Local Development (Optional)
- Clone the repo: `git clone <repo-url>`
- Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
- Create environment: `conda create -n dfu_analyzer python=3.9 && conda activate dfu_analyzer`
- Install BLAST+ (e.g., `sudo apt-get install ncbi-blast+` on Ubuntu/WSL).
- Install dependencies: `pip install -r requirements.txt`
- Run locally: `streamlit run app.py`
- Ensure `card_database/` is present with `aro_index.tsv` and `card_db.*`.

## Features
- Uses CARD database (6,439 AROs) for comprehensive ARG detection.
- Adjustable thresholds (default: 90% identity, 80% coverage).
- Interactive plot with customizable dimensions.
