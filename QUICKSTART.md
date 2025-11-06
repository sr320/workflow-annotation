# Quick Start Guide - Web Application

This guide will get you up and running with the annotation workflow web application in minutes.

## 🚀 Fastest Start (3 steps)

### 1. Install Dependencies

```bash
pip install Flask requests pandas goatools matplotlib
```

### 2. Start the Server

```bash
python app.py
```

### 3. Open Browser

Navigate to: **http://localhost:5000**

That's it! 🎉

---

## 📋 What You Can Do

### Submit a Job

1. **Via URL**: Paste a FASTA file URL (like from GitHub raw)
2. **Via Upload**: Upload a local FASTA file (max 50MB)

### Configure Parameters

- **Sequence Type**: Nucleotide (BLASTX) or Protein (BLASTP)
- **Performance**: Use DIAMOND for 10-20x faster processing
- **Threads**: Adjust CPU usage (1-40 threads)

### Monitor & Download

- Real-time job status updates (auto-refresh every 5 seconds)
- Download annotations (TSV), summary (Markdown), or all files (ZIP)
- View job history

---

## 🐳 Docker Quick Start

If you have Docker installed:

```bash
# Build and start
docker-compose up -d

# View logs
docker-compose logs -f

# Stop
docker-compose down
```

Access at: **http://localhost:5000**

---

## 🧪 Try an Example

1. Go to http://localhost:5000
2. Use this example protein URL:
   ```
   https://raw.githubusercontent.com/urol-e5/timeseries_molecular/refs/heads/main/D-Apul/data/Machinery.fasta
   ```
3. Select "Protein (BLASTP)" as sequence type
4. Click "Start Annotation"
5. Wait for processing (usually 2-5 minutes for small datasets)
6. Download your results!

---

## 📊 What You Get

After processing completes:

- **annotation_with_goslim.tsv**: Main results with GO-Slim terms
- **summary.md**: Human-readable summary with statistics
- **goslim_chart.png**: Visualization of GO categories
- All files as ZIP archive

---

## 🛠️ Requirements

**Minimum:**
- Python 3.8+
- BLAST+ or DIAMOND (for annotation engine)
- 2 CPU cores, 4GB RAM

**Recommended:**
- Python 3.10+
- DIAMOND (faster processing)
- 4+ CPU cores, 8GB+ RAM

---

## 💡 Tips

- **Use DIAMOND** for datasets with >1000 sequences
- **Adjust threads** based on your system (use `nproc` to see available cores)
- **Check job history** at /jobs to see past submissions
- **Use URLs** when possible - faster than uploading large files

---

## 🔧 Troubleshooting

### "Connection refused"
- Make sure the server is running: `python app.py`
- Check the port isn't blocked: try accessing `http://127.0.0.1:5000`

### "BLAST not found"
- Install BLAST+: `conda install -c bioconda blast`
- Or install DIAMOND: `conda install -c bioconda diamond`

### "Job stays in running state"
- Normal for large datasets - be patient!
- Check logs for errors: Look at terminal output

### "File upload fails"
- Check file size (max 50MB)
- Verify file extension (.fasta, .fa, .fas, .cds, .faa)

---

## 📚 More Information

- Full deployment guide: [WEB_DEPLOYMENT.md](WEB_DEPLOYMENT.md)
- Workflow details: [README.md](README.md)
- Issues: [GitHub Issues](https://github.com/sr320/workflow-annotation/issues)

---

## 🎯 Next Steps

Once comfortable with the web interface:

1. **Deploy to production** - See WEB_DEPLOYMENT.md
2. **Integrate with pipelines** - Use the API endpoints
3. **Customize** - Modify templates and parameters
4. **Scale up** - Use Docker and load balancers for multiple users

Happy annotating! 🧬
