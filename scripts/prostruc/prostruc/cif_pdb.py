from Bio.PDB import MMCIFParser, PDBIO
import os

def cif_to_pdb(filename,structure_id,directory):
    # Create a parser object to read the CIF file
    parser = MMCIFParser()
    structure = parser.get_structure(structure_id=structure_id, filename=filename)

    # Create an output object to write the structure to a PDB file
    io = PDBIO()
    io.set_structure(structure)
    path = os.path.join(f"{directory}",f"{structure_id}_pdb.pdb")
    io.save(path)


def convert(source_directory,output_directory):
    for path,directories,files in os.walk(source_directory):
        for pdb_file in files:
            path = os.path.join(source_directory,pdb_file)
            cif_to_pdb(filename=path,structure_id=pdb_file,directory=output_directory)


# convert(source_directory="cap5_pdb_files",output_directory="promod3_pdb_files")


