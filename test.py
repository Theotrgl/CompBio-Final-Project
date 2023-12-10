import nglview as nv
from Bio import PDB
from Bio.PDB import PDBParser, is_aa
from Bio.PDB import DSSP
from bs4 import BeautifulSoup
import os

exit = False

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
                print("Molecular Weight: ", molecular_weight)
                print("Amino Acid: ", amino_acid)
                # Save the visualization to an HTML file
                view.add_component(file_path)  # Add the component again to ensure it's loaded
                # view._set_size("600px", "400px")  # Set size (optional)
                # Export the widget as an HTML file
                filename = "protein_structure.html"
                nv.write_html(filename, [view])
                print(f"Visualization saved to {filename}")

                with open(filename, 'r') as file:
                    html_content = file.read()

                soup = BeautifulSoup(html_content, 'html.parser')

                title_element = soup.find(class_='title')
                if title_element:
                    title_element.string = 'New Title'

                # Save the modified content back to the file
                with open('new_html_file.html', 'w') as file:
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