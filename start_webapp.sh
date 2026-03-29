#!/usr/bin/env bash
# Simple startup script for the web application

set -e

echo "🧬 Starting Annotation Workflow Web Application..."
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: Python 3 is not installed"
    exit 1
fi

# Check for BLAST or DIAMOND
if ! command -v makeblastdb &> /dev/null && ! command -v diamond &> /dev/null; then
    echo "⚠️  Warning: Neither BLAST+ nor DIAMOND found"
    echo "   Please install one of them:"
    echo "   - BLAST+: conda install -c bioconda blast"
    echo "   - DIAMOND: conda install -c bioconda diamond"
    echo ""
fi

# Check dependencies
echo "Checking Python dependencies..."
if ! python3 -c "import flask, requests, pandas, goatools, matplotlib" 2>/dev/null; then
    echo "Installing Python dependencies..."
    pip install -r requirements.txt
fi

# Create necessary directories
mkdir -p uploads output blastdb

# Set default environment if not set
export PORT="${PORT:-5000}"
export DEBUG="${DEBUG:-False}"

echo ""
echo "✅ Starting web server on http://0.0.0.0:${PORT}"
echo ""
echo "   Access the application at:"
echo "   → http://localhost:${PORT}"
echo ""
echo "Press Ctrl+C to stop"
echo ""

# Start the application
python3 app.py
