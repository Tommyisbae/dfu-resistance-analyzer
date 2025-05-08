import streamlit as st
import pandas as pd
import plotly.express as px
import os
from dfu_resistance_analyzer import main as run_analysis, parse_blast_results, load_gene_mappings

st.title("DFU Resistance Analyzer")
st.markdown("Upload a bacterial genome FASTA file to detect antibiotic resistance genes.")

# Threshold settings
st.sidebar.header("Analysis Settings")
identity_threshold = st.sidebar.slider("Minimum % Identity", 20.0, 100.0, 30.0, step=5.0)
coverage_threshold = st.sidebar.slider("Minimum % Coverage", 0.0, 50.0, 5.0, step=1.0)

# File upload
fasta_file = st.file_uploader("Choose FASTA file", type=["fna"], accept_multiple_files=False)
if fasta_file:
    # Save uploaded file
    temp_fasta = f"temp_{fasta_file.name}"
    with open(temp_fasta, "wb") as f:
        f.write(fasta_file.read())
    
    # Run analysis
    try:
        run_analysis(temp_fasta)
        output_csv = f"outputs/{os.path.splitext(fasta_file.name)[0]}_report.csv"
        output_summary = output_csv.replace(".csv", "_summary.csv")
        blast_output = f"outputs/{os.path.splitext(fasta_file.name)[0]}_blast_results.txt"
        
        if os.path.exists(output_csv):
            df = pd.read_csv(output_csv)
            
            # Re-filter with user thresholds
            card_tsv = "card_database/aro_index.tsv"
            gene_mappings = load_gene_mappings(card_tsv)
            df = parse_blast_results(blast_output, gene_mappings)
            df = df[(df["Percent_Identity"] >= identity_threshold) & (df["Coverage"] >= coverage_threshold)]
            
            if not df.empty:
                # Save filtered results
                df.to_csv(output_csv, index=False)
                summary = df.groupby("Antibiotic").size().reset_index(name="ARG_Count")
                summary.to_csv(output_summary, index=False)
                
                # Display sortable table
                st.subheader("Resistance Genes Detected")
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
                    height=500,
                    text="Percent_Identity",
                    hover_data=["Coverage", "Evalue"],
                    labels={"Percent_Identity": "% Identity"}
                )
                fig.update_traces(texttemplate="%{text:.1f}%", textposition="auto")
                fig.update_layout(xaxis_tickangle=45)
                st.plotly_chart(fig, use_container_width=True)
                
                # Download button
                st.download_button(
                    label="Download Report",
                    data=df.to_csv(index=False),
                    file_name=f"{os.path.splitext(fasta_file.name)[0]}_report.csv",
                    mime="text/csv"
                )
            else:
                st.warning(f"No resistance genes detected with thresholds ({identity_threshold}% identity, {coverage_threshold}% coverage). Try lowering thresholds in the sidebar.")
        else:
            st.error("Analysis failed: No resistance genes detected or invalid FASTA format.")
    except FileNotFoundError as e:
        st.error(f"Error: {str(e)}")
    except ValueError as e:
        st.warning(f"Error: {str(e)}")
    except Exception as e:
        st.error(f"Unexpected error: {str(e)}. Please check the FASTA file and CARD database, then try again.")
    
    # Clean up
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)