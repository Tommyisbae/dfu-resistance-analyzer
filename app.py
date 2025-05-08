import streamlit as st
import pandas as pd
import plotly.express as px
import os
from dfu_resistance_analyzer import main as run_analysis

st.title("DFU Resistance Analyzer")
st.markdown("Upload a bacterial genome FASTA file to detect antibiotic resistance genes.")

# File upload
fasta_file = st.file_uploader("Choose FASTA file", type=["fna"], accept_multiple_files=False)
if fasta_file:
    # Save uploaded file
    temp_fasta = "temp.fna"
    with open(temp_fasta, "wb") as f:
        f.write(fasta_file.read())
    
    # Run analysis
    try:
        run_analysis(temp_fasta)
        output_csv = f"outputs/{os.path.splitext(fasta_file.name)[0]}_report.csv"
        if os.path.exists(output_csv):
            df = pd.read_csv(output_csv)
            
            # Display sortable table
            st.subheader("Resistance Genes Detected")
            st.dataframe(
                df[["Gene", "Antibiotic", "Percent_Identity", "Coverage", "Evalue"]]
                .sort_values("Percent_Identity", ascending=False),
                use_container_width=True
            )
            
            # Plot top 10 ARGs
            top_df = df.nlargest(10, "Percent_Identity")
            fig = px.bar(
                top_df,
                x="Gene",
                y="Percent_Identity",
                color="Antibiotic",
                title="Top 10 Resistance Genes",
                height=500,
                labels={"Percent_Identity": "% Identity"}
            )
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
            st.error("Analysis failed: No resistance genes detected or invalid FASTA format.")
    except FileNotFoundError as e:
        st.error(f"Error: {str(e)}")
    except ValueError as e:
        st.error(f"Error: {str(e)}")
    except Exception as e:
        st.error(f"Unexpected error: {str(e)}. Please check the FASTA file and try again.")
    
    # Clean up
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)