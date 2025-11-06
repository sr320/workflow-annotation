FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    ncbi-blast+ \
    curl \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy application files
COPY requirements.txt .
COPY app.py .
COPY blast2slim.sh .
COPY templates/ templates/
COPY start_webapp.sh .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt gunicorn

# Create necessary directories with proper permissions
RUN mkdir -p uploads output blastdb && \
    chmod 755 blast2slim.sh start_webapp.sh

# Expose port
EXPOSE 5000

# Set environment variables
ENV PORT=5000
ENV PYTHONUNBUFFERED=1

# Run the application with gunicorn for production
CMD ["gunicorn", "--workers", "4", "--bind", "0.0.0.0:5000", "--timeout", "600", "app:app"]
