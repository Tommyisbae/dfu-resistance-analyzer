# Use official Python 3.9 image as base
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies including BLAST+
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file (removed redundant copy)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files and CARD database
COPY app.py .
COPY dfu_resistance_analyzer.py .
COPY card_database/ /app/card_database/

# Expose port for Streamlit
EXPOSE 8501

# Set environment variables for configurable paths
ENV CARD_DATABASE_PATH=/app/card_database
ENV OUTPUT_DIR=/tmp/dfu_outputs

# Run Streamlit
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.enableCORS=false"]