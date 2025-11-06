#!/usr/bin/env python3
"""
Flask web application for the blast2slim annotation workflow.
Allows users to submit FASTA URLs and download results.
"""
import os
import subprocess
import uuid
import time
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
import threading
import json

app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max file size

# Job tracking
jobs = {}
jobs_lock = threading.Lock()

ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fas', 'cds', 'faa'}

def allowed_file(filename):
    """Check if file has allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def run_annotation(job_id, input_source, is_protein, use_diamond, threads, is_url=False):
    """Run the annotation workflow in background."""
    try:
        with jobs_lock:
            jobs[job_id]['status'] = 'running'
            jobs[job_id]['start_time'] = datetime.now().isoformat()
        
        # Build command
        cmd = ['bash', 'blast2slim.sh', '-i', input_source, '-o', f'webapp_{job_id}', '--threads', str(threads)]
        if is_protein:
            cmd.append('--protein')
        if use_diamond:
            cmd.append('--diamond')
        
        # Run the annotation
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent
        )
        
        # Find the output directory
        output_base = Path(__file__).parent / 'output' / f'webapp_{job_id}'
        if output_base.exists():
            # Find the most recent run directory
            run_dirs = sorted(output_base.glob('run_*'))
            if run_dirs:
                output_dir = run_dirs[-1]
                
                # Check for output files
                annotation_file = output_dir / 'annotation_with_goslim.tsv'
                summary_file = output_dir / 'summary.md'
                
                with jobs_lock:
                    jobs[job_id]['status'] = 'completed' if result.returncode == 0 else 'failed'
                    jobs[job_id]['end_time'] = datetime.now().isoformat()
                    jobs[job_id]['output_dir'] = str(output_dir)
                    jobs[job_id]['stdout'] = result.stdout
                    jobs[job_id]['stderr'] = result.stderr
                    jobs[job_id]['returncode'] = result.returncode
                    
                    if annotation_file.exists():
                        jobs[job_id]['files'] = {
                            'annotation': str(annotation_file),
                            'summary': str(summary_file) if summary_file.exists() else None
                        }
                        # Count results
                        with open(annotation_file) as f:
                            line_count = sum(1 for _ in f) - 1  # exclude header
                        jobs[job_id]['result_count'] = line_count
            else:
                with jobs_lock:
                    jobs[job_id]['status'] = 'failed'
                    jobs[job_id]['error'] = 'No run directory created'
        else:
            with jobs_lock:
                jobs[job_id]['status'] = 'failed'
                jobs[job_id]['error'] = 'No output directory created'
                jobs[job_id]['stdout'] = result.stdout
                jobs[job_id]['stderr'] = result.stderr
    
    except Exception as e:
        with jobs_lock:
            jobs[job_id]['status'] = 'failed'
            jobs[job_id]['error'] = str(e)

@app.route('/')
def index():
    """Main page with submission form."""
    return render_template('index.html')

@app.route('/submit', methods=['POST'])
def submit():
    """Handle job submission."""
    # Generate unique job ID
    job_id = str(uuid.uuid4())[:8]
    
    # Get parameters
    input_type = request.form.get('input_type', 'url')
    is_protein = request.form.get('sequence_type') == 'protein'
    use_diamond = request.form.get('use_diamond') == 'on'
    threads = int(request.form.get('threads', 4))
    
    # Validate threads
    if threads < 1:
        threads = 1
    elif threads > 40:
        threads = 40
    
    # Handle input
    if input_type == 'url':
        fasta_url = request.form.get('fasta_url', '').strip()
        if not fasta_url:
            return render_template('error.html', error='Please provide a FASTA URL')
        input_source = fasta_url
        is_url = True
    else:
        # Handle file upload
        if 'fasta_file' not in request.files:
            return render_template('error.html', error='No file provided')
        file = request.files['fasta_file']
        if file.filename == '':
            return render_template('error.html', error='No file selected')
        if not allowed_file(file.filename):
            return render_template('error.html', error=f'Invalid file type. Allowed: {", ".join(ALLOWED_EXTENSIONS)}')
        
        # Save uploaded file
        upload_dir = Path(__file__).parent / 'uploads'
        upload_dir.mkdir(exist_ok=True)
        filename = secure_filename(f'{job_id}_{file.filename}')
        filepath = upload_dir / filename
        file.save(str(filepath))
        input_source = str(filepath)
        is_url = False
    
    # Create job entry
    with jobs_lock:
        jobs[job_id] = {
            'id': job_id,
            'status': 'queued',
            'input_source': input_source if is_url else Path(input_source).name,
            'is_protein': is_protein,
            'use_diamond': use_diamond,
            'threads': threads,
            'created_time': datetime.now().isoformat()
        }
    
    # Start background thread
    thread = threading.Thread(
        target=run_annotation,
        args=(job_id, input_source, is_protein, use_diamond, threads, is_url)
    )
    thread.daemon = True
    thread.start()
    
    return redirect(url_for('job_status', job_id=job_id))

@app.route('/job/<job_id>')
def job_status(job_id):
    """Show job status page."""
    with jobs_lock:
        job = jobs.get(job_id)
    
    if not job:
        return render_template('error.html', error='Job not found')
    
    return render_template('job_status.html', job=job)

@app.route('/api/job/<job_id>')
def api_job_status(job_id):
    """API endpoint for job status."""
    with jobs_lock:
        job = jobs.get(job_id)
    
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    
    return jsonify(job)

@app.route('/download/<job_id>/<file_type>')
def download(job_id, file_type):
    """Download result files."""
    with jobs_lock:
        job = jobs.get(job_id)
    
    if not job or job['status'] != 'completed':
        return render_template('error.html', error='Job not found or not completed')
    
    files = job.get('files', {})
    
    if file_type == 'annotation' and files.get('annotation'):
        filepath = Path(files['annotation'])
        if filepath.exists():
            return send_file(filepath, as_attachment=True, download_name=f'annotation_{job_id}.tsv')
    elif file_type == 'summary' and files.get('summary'):
        filepath = Path(files['summary'])
        if filepath.exists():
            return send_file(filepath, as_attachment=True, download_name=f'summary_{job_id}.md')
    elif file_type == 'all' and job.get('output_dir'):
        # Create a zip of all outputs
        import zipfile
        from io import BytesIO
        
        output_dir = Path(job['output_dir'])
        if output_dir.exists():
            memory_file = BytesIO()
            with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
                for file in output_dir.iterdir():
                    if file.is_file():
                        zf.write(file, file.name)
            memory_file.seek(0)
            return send_file(
                memory_file,
                mimetype='application/zip',
                as_attachment=True,
                download_name=f'annotation_results_{job_id}.zip'
            )
    
    return render_template('error.html', error='File not found')

@app.route('/jobs')
def list_jobs():
    """List all jobs."""
    with jobs_lock:
        job_list = list(jobs.values())
    # Sort by creation time, newest first
    job_list.sort(key=lambda x: x.get('created_time', ''), reverse=True)
    return render_template('jobs_list.html', jobs=job_list)

if __name__ == '__main__':
    # Create necessary directories
    Path('uploads').mkdir(exist_ok=True)
    
    # Run the app
    port = int(os.environ.get('PORT', 5000))
    debug = os.environ.get('DEBUG', 'False').lower() == 'true'
    
    print(f"Starting annotation web application on port {port}")
    print(f"Debug mode: {debug}")
    app.run(host='0.0.0.0', port=port, debug=debug)
