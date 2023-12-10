from flask import Flask, request, render_template, redirect, url_for
import nglview as nv
from Bio import PDB
from Bio.PDB import PDBParser, is_aa
from Bio.PDB import DSSP
from bs4 import BeautifulSoup
import os
app = Flask(__name__)
def amino_acid_to_codon():
    codon_table = {
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'N': ['AAU', 'AAC'],
        'D': ['GAU', 'GAC'],
        'C': ['UGU', 'UGC'],
        'Q': ['CAA', 'CAG'],
        'E': ['GAA', 'GAG'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG'],
        'H': ['CAU', 'CAC'],
        'I': ['AUU', 'AUC', 'AUA'],
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'K': ['AAA', 'AAG'],
        'M': ['AUG'],
        'F': ['UUU', 'UUC'],
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'W': ['UGG'],
        'Y': ['UAU', 'UAC'],
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        '*': ['UAA', 'UAG', 'UGA']  # Stop codons
    }
    return codon_table

def count_codons(sequence):
    codon_table = amino_acid_to_codon()  # Use the previously defined amino_acid_to_codon function

    codon_count = {}
    for amino_acid in sequence:
        possible_codons = codon_table.get(amino_acid, [])
        if possible_codons:
            # Counting the codons
            for codon in possible_codons:
                if codon in codon_count:
                    codon_count[codon] += 1
                else:
                    codon_count[codon] = 1
    
    return codon_count
def amino_acid_sequence_to_codons(sequence):
    codon_table = amino_acid_to_codon()
    codon_sequence = []

    for amino_acid in sequence:
        possible_codons = codon_table.get(amino_acid, [])
        if possible_codons:
            codon_sequence.append(possible_codons[0])  # Choosing the first codon for simplicity
    
    return codon_sequence

def calculate_molecular_mass(structure):
    atomic_weights = {
        'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'H': 1.0078  
    }
    molecular_mass = 0.0

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_name = atom.element
                    if atom_name in atomic_weights:
                        molecular_mass += atomic_weights[atom_name]

    return molecular_mass

def extract_amino_acid_sequence(structure):
    amino_acids = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    amino_acids.append(residue.get_resname())

    return ''.join(amino_acids)

def search_codon(codon_counts, codon):
    return codon_counts.get(codon, 0)

@app.route('/', methods=['GET', 'POST'])
def process_data():
    pdb_id = ''
    viz = ''
    template_folder = 'templates'
    
    if request.method == 'POST':
        pdb_id = request.form['input_name']  # Access form data for POST request
        # Process the submitted data or perform actions
        file_path = f"PDB_Files/{pdb_id.lower()}.pdb"
        parser = PDBParser()
        structure = parser.get_structure(pdb_id.upper(), file_path)
        molecular_weight = calculate_molecular_mass(structure)
        amino_acid = extract_amino_acid_sequence(structure)
        codons = amino_acid_sequence_to_codons(amino_acid)
        codon_count = count_codons(amino_acid)
        # Handle the form submission here
        view = nv.show_structure_file(file_path)
        parser = PDBParser()
        structure = parser.get_structure(pdb_id.upper(), file_path)
        viz = structure
        view.add_component(file_path)  # Add the component again to ensure it's loaded
        filepath = os.path.join(template_folder, "new_html_file.html")
        new_filepath = os.path.join(template_folder, "protein-data.html")
        # open(filepath, 'w')
        nv.write_html(filepath, [view])
        print(f"Visualization saved to {filepath}")
        with open(filepath, 'r') as file:
            html_content = file.read()
        
        soup = BeautifulSoup(html_content, 'html.parser')

        # Create css tag in html file and append it to head tag
        css_link = soup.new_tag('link', rel='stylesheet', href='styles.css')
        head_tag = soup.find('head')
        if head_tag:
            head_tag.append(css_link)

        # Create new div component to store information of protein strand
        new_div = soup.new_tag('div')
        new_div['class'] = 'new-component'

        # Create paragraph tag for displaying data of protein strand
        data = soup.new_tag('p')
        data['class'] = 'data'
        data.string = f'Molecular Weight: {molecular_weight}\nAmino Acid Structure: {amino_acid}\nCodon Count: {codon_count}'

        # Create the form tag
        form_tag = soup.new_tag('form', action='/data', method='post')

        # Create the input tag
        input_tag = soup.new_tag('input')
        input_tag['type'] = 'text'
        input_tag['name'] = 'input_tag'

        # Create the button tag
        button_tag = soup.new_tag('button')
        button_tag['type'] = 'submit'
        button_tag.string = 'Submit'

        # Append the input and button tags to the form tag
        form_tag.append(input_tag)
        form_tag.append(button_tag)

        # Append the form, data, and new_div tags to the appropriate elements in the document
        if new_div:
            new_div.append(data)
            new_div.append(form_tag)

        existing_element = soup.find('body')
        if existing_element:
            existing_element.append(new_div)

        # Save the modified content to a new HTML file
        with open(new_filepath, 'w') as file:
            file.write(str(soup))

        return redirect(url_for('protein_data'))

    return render_template("index.html")

@app.route('/protein-data', methods=['GET'])
def protein_data():
    # Use the 'data' parameter to pass data to the new HTML file
    return render_template('protein-data.html')
if __name__ == '__main__':
    app.run(debug=True)




