import streamlit as st
import pandas as pd
import plotly.express as px
import os
import logging
from dfu_resistance_analyzer import main as run_analysis, parse_blast_results, load_gene_mappings
import io

# Attempt to import reportlab for PDF export
try:
    from reportlab.lib.pagesizes import letter
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
    from reportlab.lib import colors
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    logging.warning("reportlab not installed; PDF export disabled")

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("outputs/streamlit.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

st.title("DFU Resistance Analyzer")
st.markdown("Upload a bacterial genome FASTA file to detect antibiotic resistance genes.")

# Threshold settings
st.sidebar.header("Analysis Settings")
preset = st.sidebar.selectbox("Threshold Preset", ["Lenient (20%, 0%)", "Moderate (30%, 5%)", "Strict (50%, 10%)"])
if preset == "Lenient (20%, 0%)":
    identity_threshold, coverage_threshold = 20.0, 0.0
elif preset == "Moderate (30%, 5%)":
    identity_threshold, coverage_threshold = 30.0, 5.0
else:
    identity_threshold, coverage_threshold = 50.0, 10.0
plot_height = st.sidebar.slider("Plot Height", 400, 800, 500)
plot_title = st.sidebar.text_input("Plot Title", "Top 10 Resistance Genes")

# File upload
fasta_file = st.file_uploader("Choose FASTA file", type=["fna"], accept_multiple_files=False)
if fasta_file:
    # Save uploaded file
    temp_fasta = f"temp_{fasta_file.name}"
    logger.info(f"Saving uploaded file: {fasta_file.name} to {os.path.abspath(temp_fasta)}")
    try:
        with open(temp_fasta, "wb") as f:
            f.write(fasta_file.read())
        if not os.path.exists(temp_fasta):
            st.error(f"Failed to save {temp_fasta}")
            logger.error(f"Failed to save {temp_fasta}")
        else:
            file_size = os.path.getsize(temp_fasta) / (1024 * 1024)  # MB
            logger.info(f"Saved {temp_fasta}, size: {file_size:.2f} MB")
            if file_size < 1:
                st.warning(f"Uploaded file is small ({file_size:.2f} MB). Verify the FASTA file.")
            
            # Run analysis
            try:
                run_analysis(temp_fasta)
                output_csv = f"outputs/{os.path.splitext(fasta_file.name)[0]}_report.csv"
                output_summary = output_csv.replace(".csv", "_summary.csv")
                blast_output = f"outputs/{os.path.splitext(fasta_file.name)[0]}_blast_results.txt"
                plot_html = f"outputs/{os.path.splitext(fasta_file.name)[0]}_plot.html"
                plot_png = f"outputs/{os.path.splitext(fasta_file.name)[0]}_plot.png"
                
                # Load raw BLAST hits
                raw_hits = []
                if os.path.exists(blast_output):
                    with open(blast_output, 'r') as f:
                        raw_hits = [line.strip().split('\t') for line in f if line.strip()]
                    hit_count = len(raw_hits)
                    logger.info(f"Loaded {hit_count} raw BLAST hits from {blast_output}")
                else:
                    hit_count = 0
                    logger.warning(f"BLAST output {blast_output} not found")
                
                card_tsv = "card_database/aro_index.tsv"
                gene_mappings = load_gene_mappings(card_tsv)
                df = parse_blast_results(blast_output, gene_mappings, min_identity=identity_threshold, min_coverage=coverage_threshold)
                
                if not df.empty:
                    # Antibiotic filter
                    antibiotics = sorted(df["Antibiotic"].unique())
                    selected_antibiotic = st.multiselect("Filter by Antibiotic Class", antibiotics, default=antibiotics)
                    df = df[df["Antibiotic"].isin(selected_antibiotic)]
                    
                    # Save filtered results
                    df.to_csv(output_csv, index=False)
                    summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
                    summary.to_csv(output_summary, index=False)
                    
                    # Display sortable table
                    st.subheader(f"Resistance Genes Detected ({len(df)})")
                    st.dataframe(
                        df[["Gene", "Antibiotic", "Percent_Identity", "Coverage", "Evalue"]]
                        .sort_values("Percent_Identity", ascending=False),
                        use_container_width=True
                    )
                    
                    # Display summary table
                    st.subheader("Summary by Antibiotic")
                    st.dataframe(summary, use_container_width=True)
                    
                    # Plot top 10 ARGs
                    top_df = df.nlargest(10, "Percent_Identity")
                    fig = px.bar(
                        top_df,
                        x="Gene",
                        y="Percent_Identity",
                        color="Antibiotic",
                        title=plot_title,
                        height=plot_height,
                        text="Percent_Identity",
                        hover_data=["Coverage", "Evalue"],
                        labels={"Percent_Identity": "% Identity"}
                    )
                    fig.update_traces(texttemplate="%{text:.1f}%", textposition="auto")
                    fig.update_layout(xaxis_tickangle=45)
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Check plot files
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
                            file_name=f"{os.path.splitext(fasta_file.name)[0]}_plot.html",
                            mime="text/html"
                        )
                    
                    # Download buttons
                    st.download_button(
                        label="Download Report (CSV)",
                        data=df.to_csv(index=False),
                        file_name=f"{os.path.splitext(fasta_file.name)[0]}_report.csv",
                        mime="text/csv"
                    )
                    # PDF export
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
                            file_name=f"{os.path.splitext(fasta_file.name)[0]}_report.pdf",
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
                    st.write("Check outputs/analysis.log or outputs/streamlit.log for details or try lowering thresholds.")
            except FileNotFoundError as e:
                st.error(f"Error: {str(e)}. Check logs at outputs/analysis.log or outputs/streamlit.log.")
                logger.error(f"Analysis error: {str(e)}")
            except TimeoutError as e:
                st.error(f"Error: BLAST execution timed out after 1 hour. Try a smaller FASTA file or contact support.")
                logger.error(f"Analysis timeout: {str(e)}")
            except ValueError as e:
                st.warning(f"Error: {str(e)}")
                logger.warning(f"Analysis warning: {str(e)}")
            except Exception as e:
                st.error(f"Unexpected error: {str(e)}. Check outputs/streamlit.log.")
                logger.error(f"Unexpected error: {str(e)}")
            
            # Clean up
            if os.path.exists(temp_fasta):
                logger.info(f"Cleaning up {temp_fasta}")
                os.remove(temp_fasta)
    except Exception as e:
        st.error(f"Failed to process uploaded file: {str(e)}")
        logger.error(f"File processing error: {str(e)}")