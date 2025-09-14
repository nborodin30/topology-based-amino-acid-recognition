import os
import subprocess
from flask import Flask, render_template, request, redirect, url_for, flash, send_from_directory
from werkzeug.utils import secure_filename

# Import your processing scripts
import Topology1
import Topology2

# Initialize Flask app
app = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'

# Define folders for uploads and processed files
UPLOAD_FOLDER = 'uploads'
PROCESSED_FOLDER = 'processed'
ALLOWED_EXTENSIONS = {'sdf'}

# Configure app to use the folders
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['PROCESSED_FOLDER'] = PROCESSED_FOLDER

# Create the folders if they don't exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(PROCESSED_FOLDER, exist_ok=True)

def allowed_file(filename):
    """
    Checks if the uploaded file has an allowed extension.
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def home():
    """
    Renders the home page with a file upload form.
    """
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    """
    Handles file upload and redirects to the version selection page.
    """
    if 'file' not in request.files:
        flash('No file part')
        return redirect(url_for('home'))

    file = request.files['file']

    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('home'))

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        upload_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(upload_path)

        flash(f'File {filename} uploaded successfully. Please select a version to process it.')
        return redirect(url_for('select_version', uploaded_filename=filename))
    else:
        flash('Invalid file type. Please upload a .sdf file.')
        return redirect(url_for('home'))

@app.route('/select-version/<uploaded_filename>')
def select_version(uploaded_filename):
    """
    Renders the page for selecting the processing version.
    """
    return render_template('select_version.html', uploaded_filename=uploaded_filename)

@app.route('/process-file/<uploaded_filename>', methods=['POST'])
def process_file(uploaded_filename):
    """
    Handles version selection, processes the file, and displays download link.
    """
    version = request.form.get('version')
    upload_path = os.path.join(app.config['UPLOAD_FOLDER'], uploaded_filename)
    processed_output_path = os.path.join(app.config['PROCESSED_FOLDER'], f'annotated_{version}_{uploaded_filename}')

    try:
        if version == 'v1':
            Topology1.process_sdf(upload_path, processed_output_path)
        elif version == 'v2':
            Topology2.process_sdf(upload_path, processed_output_path)
        else:
            flash('Please select a valid version.', 'error')
            return redirect(url_for('select_version', uploaded_filename=uploaded_filename))

        flash(f'File processed successfully with {version.upper()}!')
        with open(processed_output_path, 'r', encoding='utf-8') as f:
            file_content = f.read()
        return render_template('index.html', processed_filename=os.path.basename(processed_output_path), processed_content=file_content)
    
    except Exception as e:
        flash(f'An error occurred during processing: {e}', 'error')
        return redirect(url_for('select_version', uploaded_filename=uploaded_filename))

@app.route('/download/<filename>')
def download_file(filename):
    """
    Allows the user to download the processed file.
    """
    return send_from_directory(app.config['PROCESSED_FOLDER'], filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)