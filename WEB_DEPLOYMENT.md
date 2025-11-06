# Web Application Deployment Guide

## Overview

This web application provides a user-friendly interface for the blast2slim annotation workflow. Users can submit FASTA files (via URL or upload), configure parameters, and download results.

## Features

- **Multiple Input Methods**: URL or file upload
- **Parameter Configuration**: Sequence type, DIAMOND/BLAST+, thread count
- **Job Tracking**: Monitor job status in real-time
- **Result Downloads**: Download TSV annotations, summaries, or all files as ZIP
- **Job History**: View all submitted jobs

## Quick Start (Local Development)

### Prerequisites

1. **Python 3.8+**
2. **BLAST+ or DIAMOND**
3. **System dependencies**: curl, gzip, bash

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/sr320/workflow-annotation.git
cd workflow-annotation

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Install BLAST+ or DIAMOND (via conda recommended)
conda install -c bioconda blast
# OR for faster processing:
conda install -c bioconda diamond

# 4. Start the web application
python app.py
```

The application will be available at `http://localhost:5000`

### Configuration

Environment variables:
- `PORT`: Server port (default: 5000)
- `DEBUG`: Enable debug mode (default: False)
- `SECRET_KEY`: Flask secret key for sessions (set in production)

```bash
# Example with custom settings
export PORT=8080
export DEBUG=True
python app.py
```

## Production Deployment

### Option 1: Traditional Server (Ubuntu/Debian)

#### 1. Setup Server

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install system dependencies
sudo apt install -y python3 python3-pip nginx supervisor curl gzip

# Install BLAST+ or DIAMOND
sudo apt install -y ncbi-blast+
# OR
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo mv diamond /usr/local/bin/
```

#### 2. Setup Application

```bash
# Create application user
sudo useradd -m -s /bin/bash annotationapp

# Clone repository
sudo -u annotationapp git clone https://github.com/sr320/workflow-annotation.git /home/annotationapp/workflow-annotation
cd /home/annotationapp/workflow-annotation

# Install Python dependencies
sudo -u annotationapp pip3 install -r requirements.txt

# Install gunicorn for production
sudo -u annotationapp pip3 install gunicorn
```

#### 3. Configure Supervisor

Create `/etc/supervisor/conf.d/annotation-app.conf`:

```ini
[program:annotation-app]
command=/usr/local/bin/gunicorn --workers 4 --bind 127.0.0.1:8000 app:app
directory=/home/annotationapp/workflow-annotation
user=annotationapp
autostart=true
autorestart=true
stderr_logfile=/var/log/annotation-app/err.log
stdout_logfile=/var/log/annotation-app/out.log
environment=PATH="/usr/local/bin:/usr/bin:/bin",SECRET_KEY="change-this-secret-key"
```

```bash
# Create log directory
sudo mkdir -p /var/log/annotation-app
sudo chown annotationapp:annotationapp /var/log/annotation-app

# Start the application
sudo supervisorctl reread
sudo supervisorctl update
sudo supervisorctl start annotation-app
```

#### 4. Configure Nginx

Create `/etc/nginx/sites-available/annotation-app`:

```nginx
server {
    listen 80;
    server_name your-domain.com;  # Change this
    
    client_max_body_size 50M;
    
    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 600s;
        proxy_connect_timeout 600s;
    }
}
```

```bash
# Enable site
sudo ln -s /etc/nginx/sites-available/annotation-app /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

#### 5. Setup SSL (Optional but Recommended)

```bash
# Install certbot
sudo apt install -y certbot python3-certbot-nginx

# Get SSL certificate
sudo certbot --nginx -d your-domain.com
```

### Option 2: Docker Deployment

#### 1. Create Dockerfile

Create `Dockerfile` in the repository root:

```dockerfile
FROM ubuntu:22.04

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
COPY . .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt gunicorn

# Create necessary directories
RUN mkdir -p uploads output blastdb

# Expose port
EXPOSE 5000

# Run the application
CMD ["gunicorn", "--workers", "4", "--bind", "0.0.0.0:5000", "--timeout", "600", "app:app"]
```

#### 2. Build and Run

```bash
# Build image
docker build -t annotation-webapp .

# Run container
docker run -d \
  --name annotation-app \
  -p 5000:5000 \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/blastdb:/app/blastdb \
  -e SECRET_KEY="your-secret-key" \
  annotation-webapp

# View logs
docker logs -f annotation-app
```

#### 3. Docker Compose (Recommended)

Create `docker-compose.yml`:

```yaml
version: '3.8'

services:
  web:
    build: .
    ports:
      - "5000:5000"
    volumes:
      - ./output:/app/output
      - ./blastdb:/app/blastdb
      - ./uploads:/app/uploads
    environment:
      - SECRET_KEY=${SECRET_KEY:-change-me-in-production}
      - PORT=5000
    restart: unless-stopped
```

```bash
# Start services
docker-compose up -d

# View logs
docker-compose logs -f

# Stop services
docker-compose down
```

### Option 3: Cloud Platforms

#### Heroku

```bash
# Install Heroku CLI
curl https://cli-assets.heroku.com/install.sh | sh

# Login and create app
heroku login
heroku create your-annotation-app

# Add buildpacks
heroku buildpacks:add --index 1 heroku-community/apt
heroku buildpacks:add --index 2 heroku/python

# Create Aptfile for system dependencies
echo "ncbi-blast+" > Aptfile

# Create Procfile
echo "web: gunicorn --workers 4 --timeout 600 app:app" > Procfile

# Set config
heroku config:set SECRET_KEY=$(python -c 'import secrets; print(secrets.token_hex(32))')

# Deploy
git add .
git commit -m "Deploy to Heroku"
git push heroku main

# Open app
heroku open
```

#### AWS EC2

Similar to Traditional Server deployment above, but:
1. Launch Ubuntu EC2 instance (t3.medium or larger recommended)
2. Configure security group (allow HTTP/HTTPS)
3. Follow Traditional Server setup steps
4. Consider using AWS RDS for persistent job storage (future enhancement)

#### Google Cloud Platform

```bash
# Install gcloud CLI
# Follow: https://cloud.google.com/sdk/docs/install

# Deploy to Cloud Run (containerized)
gcloud run deploy annotation-app \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --timeout 600 \
  --memory 4Gi
```

#### DigitalOcean App Platform

1. Connect your GitHub repository
2. Select "Python" as the app type
3. Add environment variables:
   - `SECRET_KEY`: Random secret key
4. Add run command: `gunicorn --workers 4 --timeout 600 app:app`
5. Deploy!

## Usage

### Web Interface

1. **Navigate to home page**
2. **Choose input method**:
   - URL: Paste FASTA URL
   - Upload: Select local file
3. **Configure parameters**:
   - Sequence type: Nucleotide or Protein
   - DIAMOND: Enable for faster processing
   - Threads: Adjust based on server capacity
4. **Submit job**
5. **Monitor progress**: Page auto-refreshes
6. **Download results**: TSV, summary, or ZIP of all files

### API Endpoints

#### Submit Job (Form)
```
POST /submit
Content-Type: multipart/form-data

Parameters:
- input_type: "url" or "file"
- fasta_url: URL to FASTA file (if input_type=url)
- fasta_file: File upload (if input_type=file)
- sequence_type: "nucleotide" or "protein"
- use_diamond: "on" or omit
- threads: Integer (1-40)
```

#### Check Job Status (API)
```
GET /api/job/<job_id>
Returns: JSON with job status and details
```

#### Download Results
```
GET /download/<job_id>/annotation  # TSV file
GET /download/<job_id>/summary     # Summary markdown
GET /download/<job_id>/all         # ZIP of all files
```

## Monitoring and Maintenance

### Check Application Status

```bash
# Supervisor
sudo supervisorctl status annotation-app

# Docker
docker ps
docker logs annotation-app

# Check disk usage
df -h
du -sh output/ blastdb/
```

### Cleanup Old Jobs

Jobs accumulate in the `output/` directory. Create a cleanup script:

```bash
#!/bin/bash
# cleanup-old-jobs.sh

# Remove job outputs older than 7 days
find /home/annotationapp/workflow-annotation/output/webapp_* -type d -mtime +7 -exec rm -rf {} +
find /home/annotationapp/workflow-annotation/uploads/ -type f -mtime +7 -delete

echo "Cleanup completed: $(date)"
```

Add to crontab:
```bash
sudo crontab -e -u annotationapp
# Add: 0 2 * * * /home/annotationapp/workflow-annotation/cleanup-old-jobs.sh
```

### View Logs

```bash
# Application logs (Supervisor)
sudo tail -f /var/log/annotation-app/out.log
sudo tail -f /var/log/annotation-app/err.log

# Nginx logs
sudo tail -f /var/log/nginx/access.log
sudo tail -f /var/log/nginx/error.log

# Docker logs
docker logs -f annotation-app
```

## Performance Considerations

### Resource Requirements

**Minimum:**
- CPU: 2 cores
- RAM: 4 GB
- Disk: 20 GB

**Recommended:**
- CPU: 4+ cores
- RAM: 8+ GB
- Disk: 50+ GB (depends on job volume)

### Optimization Tips

1. **Use DIAMOND**: 10-20x faster than BLAST+
2. **Limit concurrent jobs**: Adjust gunicorn workers
3. **Database caching**: Share `blastdb/` across deployments
4. **Job cleanup**: Regularly remove old results
5. **Reverse proxy**: Use nginx for static files and load balancing

## Troubleshooting

### Issue: Jobs stay in "running" state

**Solution**: Check that BLAST+/DIAMOND is installed and in PATH
```bash
which makeblastdb
which blastx
which diamond
```

### Issue: File upload fails

**Solution**: Check nginx client_max_body_size and Flask MAX_CONTENT_LENGTH

### Issue: Out of disk space

**Solution**: Cleanup old jobs and database files
```bash
df -h
du -sh output/* | sort -h
rm -rf output/webapp_*/
```

### Issue: Permission errors

**Solution**: Ensure application user owns directories
```bash
sudo chown -R annotationapp:annotationapp /home/annotationapp/workflow-annotation
sudo chmod -R 755 /home/annotationapp/workflow-annotation
```

## Security Considerations

1. **Change SECRET_KEY**: Use a random secret in production
2. **Input validation**: The app validates file types and URLs
3. **File size limits**: Maximum 50MB per upload
4. **HTTPS**: Always use SSL certificates in production
5. **Firewall**: Restrict access to necessary ports only
6. **Updates**: Keep dependencies updated
7. **Isolation**: Consider containerization (Docker)

## Support

For issues or questions:
- GitHub Issues: https://github.com/sr320/workflow-annotation/issues
- Documentation: See README.md for workflow details

## License

Same license as the main workflow-annotation project.
