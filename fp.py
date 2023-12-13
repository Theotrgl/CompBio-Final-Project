from flask import Flask, request, render_template, redirect, url_for, session
import nglview as nv
from Bio import PDB
from Bio.PDB import PDBParser, is_aa
from Bio.PDB import DSSP
from bs4 import BeautifulSoup
import os
app = Flask(__name__, template_folder="templates", static_folder="static_files")
app.secret_key = 'aksjfgaisfg917f19vf197vbsbfbf1s9VF(!&U(SB1u9bsf9B!F(bs)))'
def amino_acid_to_codon():
    codon_table = {
        'ALA': ['GCU', 'GCC', 'GCA', 'GCG'],
        'ARG': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'ASN': ['AAU', 'AAC'],
        'ASP': ['GAU', 'GAC'],
        'CYS': ['UGU', 'UGC'],
        'GLN': ['CAA', 'CAG'],
        'GLU': ['GAA', 'GAG'],
        'GLY': ['GGU', 'GGC', 'GGA', 'GGG'],
        'HIS': ['CAU', 'CAC'],
        'ILE': ['AUU', 'AUC', 'AUA'],
        'LEU': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'LYS': ['AAA', 'AAG'],
        'MET': ['AUG'],
        'PHE': ['UUU', 'UUC'],
        'PRO': ['CCU', 'CCC', 'CCA', 'CCG'],
        'SER': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'THR': ['ACU', 'ACC', 'ACA', 'ACG'],
        'TRP': ['UGG'],
        'TYR': ['UAU', 'UAC'],
        'VAL': ['GUU', 'GUC', 'GUA', 'GUG'],
        'STOP': ['UAA', 'UAG', 'UGA']  # Stop codons
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
    codon_table = amino_acid_to_codon()  # Use the previously defined amino_acid_to_codon function

    codons = []
    for amino_acid in sequence:
        possible_codons = codon_table.get(amino_acid, [])
        if possible_codons:
            for codon in possible_codons:
                if codon not in codons:
                    codons.append(codon)
    
    return codons  # Converting the set back to a list before returning





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

    # Combine all concatenated codes into a single string
    concatenated_codes = ''.join(amino_acids)

    # Split the concatenated string into individual three-letter codes
    separated_amino_acids = [concatenated_codes[i:i + 3] for i in range(0, len(concatenated_codes), 3)]

    return separated_amino_acids



def search_codon(codon_counts, codon):
    return codon_counts.get(codon, 0)

@app.route('/', methods=['GET', 'POST'])
def process_data():
    pdb_id = ''
    template_folder = 'templates'
    
    if request.method == 'POST':
        pdb_id = request.form['input_name']  # Access form data for POST request
        # Process the submitted data or perform actions
        file_path = f"PDB_Files/{pdb_id.lower()}.pdb"
        if not os.path.exists(file_path):
            # Append an error message to the rendered template
            error_message = "File not found. Please enter a valid PDB ID."
            return render_template("index.html", error_message=error_message)
        else:
            parser = PDBParser()
            structure = parser.get_structure(pdb_id.upper(), file_path)
            mw = calculate_molecular_mass(structure)
            molecular_weight = format(mw, '3f')
            amino_acid = extract_amino_acid_sequence(structure)
            codons = amino_acid_sequence_to_codons(amino_acid)
            codon_count = count_codons(amino_acid)
            session['codon_count'] = codon_count
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
            css_link = soup.new_tag('link', rel='stylesheet', href="{{ url_for('static', filename='style.css') }}")
            head_tag = soup.find('head')
            if head_tag:
                head_tag.append(css_link)

            container_div = soup.new_tag('div', **{'class':'container'})

            data_types = ['Molecular Weight', 'Amino Acid Structure', 'Existing Codons']

            for idx, data_type in enumerate(data_types, start=1):
                accordion_div = soup.new_tag('div', **{'class':'accordion'})

                input_checkbox = soup.new_tag('input', type='checkbox', id=f'Acc{idx}')
                label_for_checkbox = soup.new_tag('label', **{'for':f'Acc{idx}'})
                label_for_checkbox.string = data_type

                chevron_down_div = soup.new_tag('div', **{'class':'fas fa-chevron-down rotate'})

                content_div = soup.new_tag('div', **{'class':'content'})
                # Replace the static content with your dynamic variables
                if data_type == 'Molecular Weight':
                    content_div.string = f'{molecular_weight}'
                elif data_type == 'Amino Acid Structure':
                    content_div.string = f'{amino_acid}'
                elif data_type == 'Existing Codons':
                    content_div.string = f'{codons}'

                accordion_div.append(input_checkbox)
                accordion_div.append(label_for_checkbox)
                accordion_div.append(chevron_down_div)
                accordion_div.append(content_div)

                container_div.append(accordion_div)
            
            # Additional HTML structure
            container2_div = soup.new_tag('div', **{'class': 'container2'})
            box_div = soup.new_tag('div', **{'class': 'box'})
            search_bar_div = soup.new_tag('div', **{'class': 'search-bar'})
            form_tag = soup.new_tag('form', action='/protein-data', method='post')
            input_tag = soup.new_tag('input', type='text', placeholder='Search Codons')
            input_tag['name'] = 'codon'
            button_tag = soup.new_tag('button', type='submit')
            i_tag = soup.new_tag('i', **{'class': 'fas fa-search'})

            # Append the tags to form the structure
            button_tag.append(i_tag)
            form_tag.append(input_tag)
            form_tag.append(button_tag)
            search_bar_div.append(form_tag)
            box_div.append(search_bar_div)
            container2_div.append(box_div)

            # Append the new HTML structure to the existing container_div
            container_div.append(container2_div)

            existing_element = soup.find('body')
            if existing_element:
                existing_element.append(container_div)

            # Save the modified content to a new HTML file
            with open(new_filepath, 'w') as file:
                file.write(str(soup))

            return redirect(url_for('protein_data'))
    return render_template("index.html")

@app.route('/protein-data', methods=['GET', 'POST'])
def protein_data():
    if request.method == 'POST':
        codon = request.form.get('codon', '')
        codon_count = session.get('codon_count')
        searchedCodon = search_codon(codon_count, codon)
        
        new_filepath = os.path.join('templates', "protein-data.html")
        with open(new_filepath, 'r') as file:
            html_content = file.read()
        
        soup = BeautifulSoup(html_content, 'html.parser')

        existing_results = soup.find_all(class_='codon-result')
        for result in existing_results:
            result.extract()
        
        codon_search_result = soup.new_tag('p')
        codon_search_result['class'] = "codon-result"
        codon_search_result.string = f"there are {searchedCodon} {codon} sequence in this protein strand!"

        existing_element = soup.find(class_='search-bar')
        if existing_element:
            existing_element.append(codon_search_result)

        updated_content = str(soup)
        with open(new_filepath, 'w') as file:
            file.write(updated_content)
        
    return render_template('protein-data.html')

if __name__ == '__main__':
    app.run(debug=True)




