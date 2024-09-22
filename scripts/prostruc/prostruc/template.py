from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
from Bio.PDB import parse_pdb_header
from prostruc.cif_pdb import convert
import os

class Template:
    def __init__(self,file_system_info):
        self.blast_result_path = file_system_info['blast search output path']
        self.templates_directory = file_system_info['templates directory']
        self.blast_search_result = file_system_info['blast search result']


    def check_sequence_similarity(self,e_value_threshold,min_identity=30):
        try:
            match_list = []
            with open(self.blast_result_path,'r') as file_handler:
                blast_records = NCBIXML.parse(file_handler)

                for blast_record in blast_records:
                    if blast_record.alignments:
                        for i, alignment in enumerate(blast_record.alignments):
                            for hsp in alignment.hsps:
                                percentage_identity = (hsp.identities / hsp.align_length) * 100
                                if hsp.expect < e_value_threshold: #and percentage_identity > min_identity:
                                    param = {
                                        "id": alignment.hit_id.split("|")[1].split(".")[0],
                                        "evalue": hsp.expect,
                                        "identity":percentage_identity
                                    }
                                    match_list.append(param)
            return match_list

        except Exception as error:
            print(error)
            return None


    # Arguement for this method is a list of pdb file IDS
    # The return value is the list of filenames that have been downloaded
    # If return value is False, then file download failed
    def download_pdb_files(self,match_list):
        try:
            downloader = PDBList()
            id_list = []
            for item in match_list:
                item_id = item['id']
                id_list.append(item_id)
            download_result = downloader.download_pdb_files(pdb_codes=id_list, pdir= self.templates_directory)
            convert(source_directory=self.templates_directory,output_directory="promod3_input_pdb_files")
            return download_result

        except Exception as error:
            return None


    def check_for_missing_residues(self):
        try:
            missing_residue = True
            if os.path.exists(self.templates_directory):
                print(f"[*] Checking pdb file header for missing residues")
                for path,directories,files in os.walk(self.templates_directory):
                    for pdb_file in files:
                        pdb_file_path = os.path.join(self.templates_directory,pdb_file)
                        header_details = parse_pdb_header(pdb_file_path)
                        if header_details['has_missing_residues'] == True:
                            # print(f"[*] {pdb_file} has missing residues")
                            missing_residue = True
                            self.delete_file(pdb_file_path, pdb_file)
                        else:
                            missing_residue = False

                if missing_residue == False:
                    return False

                else:
                    return True

            # else:
            #     print(f"{self.hits_directory} does not exist")
            #     return f"{self.hits_directory} does not exist"

        except Exception as error:
            return error



    # Deletes all files that have missing residues
    def delete_file(self,file_path,pdb_file):
        try:
            print(f"[*] Deleting file {pdb_file}")
            os.remove(file_path)

        except Exception as error:
            return False


