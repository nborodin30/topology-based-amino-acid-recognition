# Topology-based Amino Acid Recognition

The **Topology-based Amino Acid Recognition** project is a web-based tool designed for the **automated recognition of amino acid residues** in molecular structures. It addresses the critical problem of reconstructing amino acid sequences when original residue identifiers are lost during data processing or file format conversions. The tool uses only the topology of atoms and their bonds to accurately identify amino acids.

---

## Usage (Online)

The simplest way to use the tool is to access the live version online. You don't need to install anything.

1.  Open the website: [https://aarecognition.pythonanywhere.com/](https://aarecognition.pythonanywhere.com/)
2.  Upload your `.sdf` file, and the system will process it automatically.

---

## Installation

If you would like to run the project locally, follow these instructions.

### Prerequisites

* **Python 3.13**
* **pip**

### Step-by-Step Guide

1.  Clone the project repository:
    ```bash
    git clone [https://github.com/nborodin30/topology-based-amino-acid-recognition.git](https://github.com/nborodin30/topology-based-amino-acid-recognition.git)
    cd topology-based-amino-acid-recognition
    ```

2.  Install all necessary dependencies from the `requirements.txt` file:
    ```bash
    pip install -r requirements.txt
    ```
    Key dependencies include **Flask** (for the web interface), **networkx** (for graph operations), and **GraKeL** (for structural comparisons).

3.  Run the application:
    ```bash
    python app.py
    ```

4.  Open your web browser and navigate to the address shown in the console (usually `http://127.0.0.1:5000`).

---

## About the Project

Two distinct versions of the algorithm were developed, each with unique characteristics tailored for different use cases:

* **Version 1 (Pattern Matching)**: This version uses a streamlined pattern matching approach. It relies on sequential atom IDs and simple heuristics to detect common amino acids. It is optimized for speed and is well-suited for rapid analysis of small, well-ordered peptides.
* **Version 2 (Graph-Based Recognition)**: This approach leverages graph theory, representing the molecule as a graph and using algorithms like graph isomorphism and the Weisfeiler-Lehman kernel to compare it against known amino acid templates. This version is more robust to "noisy" data and can accurately identify all 20 canonical amino acids.
