# Import Necessary Modules
from flask import Flask, request, render_template, redirect, url_for, session
import nglview as nv
from Bio import PDB, SeqIO
from Bio.PDB import PDBParser, is_aa
from bs4 import BeautifulSoup
import os
import time
import tracemalloc
import psutil
# Initialize flask app
app = Flask(__name__, template_folder="templates", static_folder="static_files")
app.secret_key = 'aksjfgaisfg917f19vf197vbsbfbf1s9VF(!&U(SB1u9bsf9B!F(bs)))'

#Functions for extracting and manipulating protein data

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
        'STOP': ['UAA', 'UAG', 'UGA'] 
    }
    return codon_table
#Convert Codon triplet to single alphabet amino acid notation
def change_amino_acid_format(codon_sequence):
    #Initializes dictionary for matching with codon sequence
    three_letter_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

    #Reads every triplet (codon triplet) and matches against dictionary
    amino_acid_sequence = ''
    for i in range(0, len(codon_sequence), 3):
        codon = codon_sequence[i:i + 3]
        amino_acid_sequence += three_letter_to_one_letter.get(codon, '')

    return amino_acid_sequence

#Counts the unique codon/amino acid in the sequence
def count_codons(sequence):
    codon_table = amino_acid_to_codon()  
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


#Calculate molecular mass by adding the existing atoms found in the PDB file
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

#Extracting amino acid using PDB.is_aa() function
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

#Extract Title of protein strand
def extract_protein_title(file_path):
    titles = []
    with open(file_path, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if line.startswith("TITLE"):
                title = line.strip().split(' ', 1)[1]
                titles.append(title)

    return titles if titles else None

#Extract polypeptide chain from PDB file 
def extract_chain_from_pdb(structure, chain_id):

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                # Collect residues to form the sequence
                residues = [residue.get_resname() for residue in chain if residue.get_id()[0] == ' ']
                sequence = ''.join(residues)
                sequence = change_amino_acid_format(sequence)
                # Convert the sequence to FASTA format to match clustal omega template
                fasta_seq = f">{chain.id}\n{sequence}"
                return fasta_seq

    return None
#Create index page
@app.route('/', methods=['GET', 'POST'])
def process_data():
    pdb_id = ''
    template_folder = 'templates'
    #Pre-renders the /protein-data page based on user input
    if request.method == 'POST':
        #Fetching user input data
        pdb_id = request.form['input_name']  
        #Loading /protein-data page based on user input (PDB ID)
        file_path = f"PDB_Files/{pdb_id.lower()}.pdb"
        #Error checks for valid PDB ID
        if not os.path.exists(file_path):
            error_message = "File not found. Please enter a valid PDB ID."
            return render_template("index.html", error_message=error_message)
        else:
            #Initialize function values into variables
            parser = PDBParser()
            structure = parser.get_structure(pdb_id.upper(), file_path)
            molecular_weight = calculate_molecular_mass(structure)
            amino_acid = extract_amino_acid_sequence(structure)
            protein_A = extract_chain_from_pdb(structure, 'A')
            protein_B = extract_chain_from_pdb(structure, 'B')
            protein_C = extract_chain_from_pdb(structure, 'C')
            pdb_title = extract_protein_title(file_path)
            mw = "{:.2f}".format(molecular_weight)
            codons = list(set(amino_acid))
            codon_count = count_codons(amino_acid)

            #Add session data to carry over to other pages
            session['codon_count'] = codon_count
            session['aminoAcid'] = amino_acid

            #Use nglview to display PDB structure
            view = nv.show_structure_file(file_path)
            #Initialize PDBParser()
            parser = PDBParser()
            #Fetches pdb data from PDB_Files folder
            structure = parser.get_structure(pdb_id.upper(), file_path)
            viz = structure
            view.add_component(file_path)
            filepath = os.path.join(template_folder, "new_html_file.html")
            new_filepath = os.path.join(template_folder, "protein-data.html")
            nv.write_html(filepath, [view])
            print(f"Visualization saved to {filepath}")
            #Open new html file for /protein-data page and modify using soup
            with open(filepath, 'r') as file:
                html_content = file.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            title_tag = soup.title

            # Update the title content
            title_tag.string = 'Protein-data'
            # Create css tag in html file and append it to head tag
            css_link = soup.new_tag('link', rel='stylesheet', href="{{ url_for('static', filename='style.css') }}")
            head_tag = soup.find('head')
            if head_tag:
                head_tag.append(css_link)
            
            container_div = soup.new_tag('div', **{'class':'container'})

            data_types = ['ID','Title', 'Molecular Weight', 'Polypeptide Chain A', 'Polypeptide Chain B', 'Polypeptide Chain C','Amino Acid Structure', 'Existing Codons']
            session['dataType'] = data_types
            for idx, data_type in enumerate(data_types, start=1):
                accordion_div = soup.new_tag('div', **{'class':'accordion'})

                input_checkbox = soup.new_tag('input', type='checkbox', id=f'Acc{idx}')
                label_for_checkbox = soup.new_tag('label', **{'for':f'Acc{idx}'})
                label_for_checkbox.string = data_type

                chevron_down_div = soup.new_tag('div', **{'class':'fas fa-chevron-down rotate'})

                content_div = soup.new_tag('div', **{'class':'content'})
                # Replace the static content with your dynamic variables
                if data_type == 'ID':
                    content_div.string = f'{pdb_id.upper()}'
                elif data_type == 'Title':
                    content_div.string = f'{pdb_title}'
                elif data_type == 'Polypeptide Chain A':
                    content_div.string = f'{protein_A}'
                elif data_type == 'Polypeptide Chain B':
                    content_div.string = f'{protein_B}'
                elif data_type == 'Polypeptide Chain C':
                    content_div.string = f'{protein_C}'
                elif data_type == 'Molecular Weight':
                    content_div.string = f'{mw, "Da"}'
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
            SubHeading = soup.new_tag('h3')
            SubHeading['class'] = 'sub-heading'
            SubHeading.string = 'Input Codon To Search'
            # Append the tags to form the structure
            button_tag.append(i_tag)
            form_tag.append(input_tag)
            form_tag.append(button_tag)
            search_bar_div.append(form_tag)
            box_div.append(search_bar_div)
            container2_div.append(box_div)

            # Append the new HTML structure to the existing container_div
            container_div.append(SubHeading)
            container_div.append(container2_div)
        
            link_tag = soup.new_tag('a', href='/')
            link_tag.string = 'Go back to Homepage'
            container_div.append(link_tag)
            existing_element = soup.find('body')
            if existing_element:
                existing_element.append(container_div)

            # Save the modified content to a new HTML file
            with open(new_filepath, 'w') as file:
                file.write(str(soup))

            return redirect(url_for('protein_data'))
    return render_template("index.html")

@app.route('/protein-data/', methods=['GET', 'POST'])
def protein_data():
    if request.method == 'POST':
        codon = request.form.get('codon', '')
        codon = codon.upper()
        codon_count = session.get('codon_count')
        aminoAcid = session.get('aminoAcid')
        data_types = session.get('dataType')
        codon_idx = [i for i, x in enumerate(aminoAcid) if x == codon]
        searchedCodon = aminoAcid.count(codon.upper())
        aminoAcid_copy = aminoAcid
        new_filepath = os.path.join('templates', "protein-data.html")
        with open(new_filepath, 'r') as file:
            html_content = file.read()

        for index in codon_idx:
            if index < len(aminoAcid_copy):
                aminoAcid_copy[index] = f'<span style="color: red;">{aminoAcid_copy[index]}</span>'
        highlighted_str = ''.join(aminoAcid_copy)
        
        soup = BeautifulSoup(html_content, 'html.parser')
        content_div_separated_aa = soup.find('label', string='Amino Acid Structure').find_next('div', {'class': 'content'})
        # Clear previous content in the div
        content_div_separated_aa.clear()

        # Create a new tag for the modified content
        modified_content_tag = soup.new_tag('span')
        modified_content_tag.append(BeautifulSoup(highlighted_str, 'html.parser'))

        # Append the modified content to the div
        content_div_separated_aa.append(modified_content_tag)
        existing_results = soup.find_all(class_='codon-result')
        for result in existing_results:
            result.extract()
            
        if (searchedCodon == 0):
            codon_search_result = soup.new_tag('p')
            codon_search_result['class'] = "codon-result"
            codon_search_result.string = f"{codon} sequence cannot be found in this protein stand"
        else:
            codon_search_result = soup.new_tag('p')
            codon_search_result['class'] = "codon-result"
            codon_search_result.string = f"There are {searchedCodon} {codon} sequence in this protein strand! found in indexes {codon_idx}"

        existing_element = soup.find(class_='search-bar')
        if existing_element:
            existing_element.append(codon_search_result)

        updated_content = str(soup)
        with open(new_filepath, 'w') as file:
            file.write(updated_content)
        
    return render_template('protein-data.html')

if __name__ == '__main__':
    app.run(debug=True)




