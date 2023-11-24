import nglview as nv
from bs4 import BeautifulSoup
import os

exit = False

folder_path = "./PDB_Files"
while exit == False:
    pdb_id = input("Please Input a valid PDB ID to query: ")
    
    if pdb_id.lower() == "exit":
        print("Successfully exitted!!")
        exit_loop = True
        break

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        # Check if the given_filename exists in the folder
        file_exists = pdb_id in os.listdir(folder_path)
        if file_exists:
            # Create an NGLview widget instance
            view = nv.show_structure_file(f"./PDB_Files/{pdb_id.lower()}.pdb")

            # Save the visualization to an HTML file
            view.add_component(f"./PDB_Files/{pdb_id.lower()}.pdb")  # Add the component again to ensure it's loaded
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
