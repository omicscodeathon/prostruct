import subprocess
import os

def model(alignment_directory,structure_directory):
    aligned_files_list = []
    structure_files_list = []
    host_directory = os.getcwd()

    for path,directories,alignment_files in os.walk(alignment_directory):
        for aligned_file in alignment_files:
            aligned_files_list.append(aligned_file)

    for path,directories,structure_files in os.walk(structure_directory):
        for structure_file in structure_files:
            structure_files_list.append(structure_file)

    for i in range(len(aligned_files_list)):
        docker_command = f"docker run -v {host_directory}:/prostruc registry.scicore.unibas.ch/schwede/promod3 build-model -f /prostruc/alignment_files/{aligned_files_list[i]} -p /prostruc/promod3_input_pdb_files/{structure_files_list[i]} -o /prostruc/modeled_structures/{structure_files_list[i]}_output_model.pdb"
        result = subprocess.run(docker_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print(f"{aligned_files_list[i]} built successfully")
        else:
            print(f"ERROR: {aligned_files_list[i]} build failed")
