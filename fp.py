import os
import gzip
from Bio import PDB


def load_specific_file(folder_path):
    # List all files in the folder
    files = os.listdir(folder_path)
    pdb_files = [f for f in files if f.endswith('.pdb.gz')]

    # Display available files to the user
    print("Available files:")
    for i, file_name in enumerate(pdb_files):
        print(f"{i + 1}. {file_name}")

    # Prompt user to select a file
    file_index = int(input("Enter the index of the file to load: ")) - 1

    # Check if the index is valid
    if 0 <= file_index < len(pdb_files):
        selected_file = pdb_files[file_index]
        file_path = os.path.join(folder_path, selected_file)

        # Decompress the selected file
        with gzip.open(file_path, 'rb') as f_in:
            # Remove the '.gz' extension to get the file name without compression
            decompressed_file_path = file_path[:-3]
            with open(decompressed_file_path, 'wb') as f_out:
                f_out.write(f_in.read())

        # Load the decompressed file using Biopython's PDB module
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('selected_structure', decompressed_file_path)

        # Process the structure (example: print chain information)
        for model in structure:
            for chain in model:
                print(f"Chain {chain.get_id()} in the selected structure")

        # Clean up: Remove the decompressed file
        os.remove(decompressed_file_path)

        return structure
    else:
        print("Invalid index. Please select a valid file index.")

# Provide the path to your PDB_Files folder containing .pdb.gz files
folder_path = 'PDB_Files'

# Call the function with the folder path to load a specific file based on user input
load_specific_file(folder_path)
