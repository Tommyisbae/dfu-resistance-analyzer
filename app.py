import streamlit as st
import pandas as pd
import plotly.express as px
import os
import logging
from dfu_resistance_analyzer import main as run_analysis, parse_blast_results, load_gene_mappings

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
                        title="Top 10 Resistance Genes",
                        height=plot_height,
                        text="Percent_Identity",
                        hover_data=["Coverage", "Evalue"],
                        labels={"Percent_Identity": "% Identity"}
                    )
                    fig.update_traces(texttemplate="%{text:.1f}%", textposition="auto")
                    fig.update_layout(xaxis_tickangle=45)
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Download buttons
                    st.download_button(
                        label="Download Report",
                        data=df.to_csv(index=False),
                        file_name=f"{os.path.splitext(fasta_file.name)[0]}_report.csv",
                        mime="text/csv"
                    )
                    st.download_button(
                        label="Download Plot (HTML)",
                        data=open(output_csv.replace(".csv", "_plot.html"), "rb").read(),
                        file_name=f"{os.path.splitext(fasta_file.name)[0]}_plot.html",
                        mime="text/html"
                    )
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
                st.error(f"Error: {str(e)}. Ensure the FASTA file name matches the uploaded file (e.g., GCA_ vs. GCF_).")
                logger.error(f"Analysis error: {str(e)}")
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