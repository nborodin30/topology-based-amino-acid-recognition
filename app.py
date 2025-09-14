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
    Handles file upload, processing, and download.
    """
    if 'file' not in request.files:
        flash('No file part')
        return redirect(url_for('home'))

    file = request.files['file']
    version = request.form.get('version')

    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('home'))

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        upload_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(upload_path)
        
        # Determine which script to run based on the user's selection
        processed_output_path = os.path.join(app.config['PROCESSED_FOLDER'], f'annotated_v1_{filename}')
        
        try:
            if version == 'v1':
                # Call the main function of Topology1.py
                Topology1.process_sdf(upload_path, processed_output_path)
            elif version == 'v2':
               # Call the main function of Topology2.py
                Topology2.process_sdf(upload_path, processed_output_path)
            else:
                flash('Please select a version.', 'error')
                return redirect(url_for('home'))

            flash(f'File processed successfully with Version {version.upper()}! You can download it below.')
            
            # Pass the processed file name to the template for download link
            return render_template('index.html', processed_filename=os.path.basename(processed_output_path))

        except Exception as e:
            flash(f'An error occurred during processing: {e}', 'error')
            return redirect(url_for('home'))

    else:
        flash('Invalid file type. Please upload a .sdf file.')
        return redirect(url_for('home'))

@app.route('/download/<filename>')
def download_file(filename):
    """
    Allows the user to download the processed file.
    """
    return send_from_directory(app.config['PROCESSED_FOLDER'], filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
