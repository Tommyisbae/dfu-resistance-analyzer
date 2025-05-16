import streamlit as st
import pandas as pd
import plotly.express as px
import os
import logging
import tempfile
import re
import shutil
import psutil
from dfu_resistance_analyzer import main as run_analysis, parse_blast_results, load_gene_mappings
import io

# Get environment variables for paths
CARD_DATABASE_PATH = os.getenv("CARD_DATABASE_PATH", "/app/card_database")
OUTPUT_DIR = os.getenv("OUTPUT_DIR", "/tmp/dfu_outputs")

# Attempt to import reportlab for PDF export
try:
    from reportlab.lib.pagesizes import letter
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
    from reportlab.lib import colors
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    logging.warning("reportlab not installed; PDF export disabled")

# Ensure output directory exists before setting up logging
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Setup logging (avoid sensitive data)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(OUTPUT_DIR, "streamlit.log")),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def sanitize_filename(filename):
    """Sanitize file name to remove invalid characters."""
    return re.sub(r'[^a-zA-Z0-9._-]', '_', filename)

def check_environment():
    """Check system environment for analysis."""
    blast_path = shutil.which("blastn")
    if not blast_path:
        logger.error("blastn executable not found in PATH")
        raise RuntimeError("blastn executable not found in PATH")
    
    mem = psutil.virtual_memory()
    if mem.available < 2 * 1024 * 1024 * 1024:  # Less than 2GB
        logger.error(f"Insufficient memory: {mem.available / (1024 * 1024)} MB available")
        raise RuntimeError("Insufficient memory. Free up at least 2GB.")
    
    disk = psutil.disk_usage('/')
    if disk.free < 1 * 1024 * 1024 * 1024:  # Less than 1GB
        logger.error(f"Insufficient disk space: {disk.free / (1024 * 1024)} MB free")
        raise RuntimeError("Insufficient disk space. Free up at least 1GB.")
    
    card_tsv = os.path.abspath(os.path.join(CARD_DATABASE_PATH, "aro_index.tsv"))
    db_path = os.path.abspath(os.path.join(CARD_DATABASE_PATH, "card_db"))
    if not os.path.exists(card_tsv):
        logger.error(f"CARD index file {card_tsv} not found")
        raise FileNotFoundError(f"CARD index file {card_tsv} not found")
    if not all(os.path.exists(f"{db_path}.{ext}") for ext in ["nhr", "nin", "nsq"]):
        logger.error(f"BLAST database {db_path} incomplete")
        raise FileNotFoundError(f"BLAST database {db_path} incomplete")
    
    logger.info(f"Environment: blastn at {blast_path}, {mem.available / (1024 * 1024)} MB free, {disk.free / (1024 * 1024)} MB disk, {os.cpu_count()} CPUs")
    return True

st.title("DFU Resistance Analyzer")
st.markdown("Upload a bacterial genome FASTA file to detect antibiotic resistance genes. Files are processed securely and deleted after use.")

# Disclaimer about file handling
st.markdown("**Note**: Uploaded files are sensitive and handled securely. They are stored temporarily in memory, processed, and deleted immediately after results are generated.")

# File splitting guide
st.markdown("**For files >100MB**: Convert .tgz to FASTA (e.g., `tar -xzf file.tgz && seqtk seq -a file.fastq > file.fasta`) and split into 100MB chunks (e.g., `seqkit split -s 100000000 file.fasta`). Upload chunks separately.")

# Threshold settings
st.sidebar.header("Analysis Settings")
preset = st.sidebar.selectbox("Threshold Preset", ["Strict (90%, 80%)", "Moderate (70%, 60%)", "Lenient (50%, 40%)"])
if preset == "Strict (90%, 80%)":
    identity_threshold, coverage_threshold = 90.0, 80.0
elif preset == "Moderate (70%, 60%)":
    identity_threshold, coverage_threshold = 70.0, 60.0
else:
    identity_threshold, coverage_threshold = 50.0, 40.0
plot_height = st.sidebar.slider("Plot Height", 400, 1000, 600)
plot_width = st.sidebar.slider("Plot Width", 600, 1200, 800)

# Environment check
try:
    check_environment()
except Exception as e:
    st.error(f"Environment setup failed: {str(e)}. Contact support.")
    logger.error(f"Environment check failed: {str(e)}")
    st.stop()

# File upload with 100MB limit
fasta_file = st.file_uploader("Choose FASTA file", type=["fna", "fasta"], accept_multiple_files=False)
MAX_FILE_SIZE = 100 * 1024 * 1024  # 100MB in bytes
if fasta_file:
    if fasta_file.size > MAX_FILE_SIZE:
        st.error(f"File size exceeds 100MB limit. Please split into smaller files (e.g., using `seqkit split -s 100000000 file.fasta`).")
    elif not fasta_file.name.lower().endswith((".fna", ".fasta")):
        st.warning("File does not end with .fna or .fasta. Ensure it’s a valid FASTA file.")
    else:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_fasta = os.path.join(temp_dir, f"temp_{sanitize_filename(fasta_file.name)}")
            logger.info(f"Saving uploaded file: {fasta_file.name} to {temp_fasta}")
            try:
                with open(temp_fasta, "wb") as f:
                    content = fasta_file.read()
                    f.write(content)
                    logger.info(f"Wrote {len(content)} bytes to {temp_fasta}")
                if not os.path.exists(temp_fasta):
                    st.error(f"Failed to save {temp_fasta}")
                    logger.error(f"Failed to save {temp_fasta}")
                    raise FileNotFoundError(f"Failed to save {temp_fasta}")
                
                file_size = os.path.getsize(temp_fasta) / (1024 * 1024)  # MB
                logger.info(f"Saved {temp_fasta}, size: {file_size:.2f} MB")
                if file_size < 1:
                    st.warning(f"Uploaded file is small ({file_size:.2f} MB). Verify the FASTA file.")
                
                with st.spinner("Running BLAST analysis..."):
                    try:
                        run_analysis(temp_fasta)
                        base_name = sanitize_filename(os.path.splitext(fasta_file.name)[0])
                        output_csv = os.path.join(OUTPUT_DIR, f"{base_name}_report.csv")
                        output_summary = output_csv.replace(".csv", "_summary.csv")
                        blast_output = os.path.join(OUTPUT_DIR, f"{base_name}_blast_results.txt")
                        plot_html = os.path.join(OUTPUT_DIR, f"{base_name}_plot.html")
                        plot_png = os.path.join(OUTPUT_DIR, f"{base_name}_plot.png")
                        
                        logger.info(f"Expected output files: {blast_output}, {output_csv}, {plot_html}, {plot_png}")
                        
                        raw_hits = []
                        if os.path.exists(blast_output):
                            with open(blast_output, 'r') as f:
                                raw_hits = [line.strip().split('\t') for line in f if line.strip()]
                            hit_count = len(raw_hits)
                            logger.info(f"Loaded {hit_count} raw BLAST hits from {blast_output}")
                        else:
                            hit_count = 0
                            logger.warning(f"BLAST output {blast_output} not found")
                        
                        card_tsv = os.path.abspath(os.path.join(CARD_DATABASE_PATH, "aro_index.tsv"))
                        gene_mappings = load_gene_mappings(card_tsv)
                        df = parse_blast_results(blast_output, gene_mappings, min_identity=identity_threshold, min_coverage=coverage_threshold)
                        
                        if not df.empty:
                            antibiotics = sorted(df["Antibiotic"].unique())
                            selected_antibiotic = st.multiselect("Filter by Antibiotic Class", antibiotics, default=antibiotics)
                            df = df[df["Antibiotic"].isin(selected_antibiotic)]
                            
                            df.to_csv(output_csv, index=False)
                            summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
                            summary.to_csv(output_summary, index=False)
                            
                            st.subheader(f"Resistance Genes Detected ({len(df)})")
                            st.dataframe(
                                df[["Gene", "Antibiotic", "Percent_Identity", "Coverage", "Evalue"]]
                                .sort_values("Percent_Identity", ascending=False),
                                use_container_width=True
                            )
                            
                            st.subheader("Summary by Antibiotic")
                            st.dataframe(summary, use_container_width=True)
                            
                            top_df = df.nlargest(10, "Percent_Identity")
                            fig = px.bar(
                                top_df,
                                x="Gene",
                                y="Percent_Identity",
                                color="Antibiotic",
                                title="Top 10 Resistance Genes",
                                height=plot_height,
                                width=plot_width,
                                text="Percent_Identity",
                                hover_data=["Coverage", "Evalue"],
                                labels={"Percent_Identity": "% Identity"}
                            )
                            fig.update_traces(texttemplate="%{text:.1f}%", textposition="auto")
                            fig.update_layout(xaxis_tickangle=45)
                            st.plotly_chart(fig, use_container_width=True)
                            
                            if not os.path.exists(plot_png):
                                st.warning("PNG plot generation failed.")
                                logger.warning(f"PNG plot missing: {plot_png}")
                            if not os.path.exists(plot_html):
                                st.warning("HTML plot generation failed. Displaying interactive plot above.")
                                logger.warning(f"HTML plot missing: {plot_html}")
                            else:
                                st.download_button(
                                    label="Download Plot (HTML)",
                                    data=open(plot_html, "rb").read(),
                                    file_name=f"{base_name}_plot.html",
                                    mime="text/html"
                                )
                                st.download_button(
                                    label="Download Plot (PNG)",
                                    data=open(plot_png, "rb").read(),
                                    file_name=f"{base_name}_plot.png",
                                    mime="image/png"
                                )
                            
                            st.download_button(
                                label="Download Report (CSV)",
                                data=df.to_csv(index=False),
                                file_name=f"{base_name}_report.csv",
                                mime="text/csv"
                            )
                            if REPORTLAB_AVAILABLE:
                                pdf_buffer = io.BytesIO()
                                doc = SimpleDocTemplate(pdf_buffer, pagesize=letter)
                                data = [df.columns.tolist()] + df[["Gene", "Antibiotic", "Percent_Identity", "Coverage", "Evalue"]].values.tolist()
                                table = Table(data)
                                table.setStyle(TableStyle([
                                    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                                    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                                    ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                                    ('GRID', (0, 0), (-1, -1), 1, colors.black)
                                ]))
                                doc.build([table])
                                pdf_buffer.seek(0)
                                st.download_button(
                                    label="Download Report (PDF)",
                                    data=pdf_buffer,
                                    file_name=f"{base_name}_report.pdf",
                                    mime="application/pdf"
                                )
                            else:
                                st.warning("PDF export unavailable. Install reportlab for PDF functionality.")
                                logger.warning("PDF export skipped due to missing reportlab")
                        else:
                            st.warning(f"No resistance genes detected with thresholds ({identity_threshold}% identity, {coverage_threshold}% coverage).")
                            st.write(f"Raw BLAST hits: {hit_count}")
                            if raw_hits:
                                st.subheader("Sample of Raw BLAST Hits")
                                raw_df = pd.DataFrame(
                                    raw_hits[:5],
                                    columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
                                )
                                st.dataframe(raw_df[["sseqid", "pident", "length", "evalue"]], use_container_width=True)
                            st.write("Check logs for details or try lowering thresholds.")
                        
                        # Auto-delete temporary files
                        for file in os.listdir(temp_dir):
                            file_path = os.path.join(temp_dir, file)
                            if os.path.isfile(file_path):
                                os.remove(file_path)
                        for file in os.listdir(OUTPUT_DIR):
                            file_path = os.path.join(OUTPUT_DIR, file)
                            if os.path.isfile(file_path):
                                os.remove(file_path)
                        logger.info(f"Deleted temporary files from {temp_dir} and {OUTPUT_DIR}")
                    
                    except FileNotFoundError as e:
                        st.error(f"Error: {str(e)}. Check logs.")
                        logger.error(f"Analysis error: {str(e)}")
                    except TimeoutError as e:
                        st.error(f"Error: BLAST execution timed out after 30 minutes. Try a smaller FASTA file.")
                        logger.error(f"Analysis timeout: {str(e)}")
                    except ValueError as e:
                        st.warning(f"Error: {str(e)}")
                        logger.warning(f"Analysis warning: {str(e)}")
                    except Exception as e:
                        st.error(f"Unexpected error: {str(e)}. Check logs.")
                        logger.error(f"Unexpected error: {str(e)}")
            
            except Exception as e:
                st.error(f"Failed to process uploaded file: {str(e)}")
                logger.error(f"File processing error: {str(e)}")