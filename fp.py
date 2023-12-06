import nglview as nv
from Bio import PDB
from Bio.PDB import PDBParser, is_aa
from Bio.PDB import DSSP
from bs4 import BeautifulSoup
import os

exit = False
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
def main():
    global exit
    folder_path = "./PDB_Files"
    while exit == False:
        pdb_id = input("Please Input a valid PDB ID to query: ")
        print(pdb_id.lower())
        file_path = f"/{pdb_id.lower()}.pdb"
        if pdb_id.lower() == "exit":
            print("Successfully exitted!!")
            exit = True
            break

        if os.path.exists(folder_path) and os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, f"{pdb_id.lower()}.pdb")
            if os.path.exists(file_path):  # Check if the file exists
                # Create an NGLview widget instance
                view = nv.show_structure_file(file_path)
                parser = PDBParser()
                structure = parser.get_structure(pdb_id.upper(), file_path)
                molecular_weight = calculate_molecular_mass(structure)
                amino_acid = extract_amino_acid_sequence(structure)
                codons = amino_acid_sequence_to_codons(amino_acid)
                codon_count = count_codons(amino_acid)
                view.add_component(file_path)  # Add the component again to ensure it's loaded
                # Export the widget as an HTML file
                filename = f"{pdb_id.lower()}_protein_structure.html"
                nv.write_html(filename, [view])
                print(f"Visualization saved to {filename}")
                # Save the visualization to an HTML file
                with open(filename, 'r') as file:
                    html_content = file.read()

                soup = BeautifulSoup(html_content, 'html.parser')
                new_div = soup.new_tag('div')
                new_div['class'] = 'new-component'
                form_tag = soup.new_tag('form', action='/', method='post')
                input_tag = soup.new_tag('input')
                button_tag = soup.new_tag('button', type='submit')
                button_tag.string = 'Submit'
                # new_p = soup.new_tag('p')
                # new_p['class'] = 'test'
                # new_p.string = '{{amino_acid}}'
                # new_div.string = f'Molecular Weight: {molecular_weight}\nAmino Acid Structure: {amino_acid}\nCodon Count: {codon_count}'
                title_element = soup.find(class_='title')
                if title_element:
                    title_element.string = 'New Title'
                existing_element = soup.find('body')
                if existing_element:
                    existing_element.append(new_div)
                    existing_element.append(form_tag)
                elem = soup.find(class_ = 'new-component')
                if elem:
                    elem.append(input_tag)
                    elem.append(button_tag)
                # find_new_div = soup.find(class_ = 'new-component')
                # if find_new_div:
                #     find_new_div.append(new_p)
                # Save the modified content back to the file
                with open(f'new-{pdb_id.lower()}_html_file.html', 'w') as file:
                    file.write(str(soup))
                
                if os.path.exists(filename):
                    os.remove(filename)
                    print(f"{filename} has been deleted.")
                else:
                    print(f"{filename} does not exist.")
            else:
                print(f"{pdb_id} does not exist in our database or is an invalid PDB ID!!")
        else:
            print(f"{folder_path} does not exist or is not a directory.")

if __name__ == "__main__":
    main()