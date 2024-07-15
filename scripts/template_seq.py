from Bio import SearchIO
from Bio.PDB import PDBList
from Bio.PDB import parse_pdb_header
import os

class Template:
    def __init__(self,file_system_info):
        self.blast_result_path = file_system_info['blast search output path']
        self.hits_directory = file_system_info['hits directory']
        self.blast_search_result = file_system_info['blast search result']


    def check_sequence_similarity(self, e_value_threshold):
        try:
            records = SearchIO.read(self.blast_result_path, 'blast-xml')
            counter = 0
            hit_list = []
            for record in records:
                # print(record)
                counter = counter +1
                for data in record:
                    if data.evalue < e_value_threshold:
                        param = {
                            "id": record.id,
                            "evalue": data.evalue,
                            "bit score": data.bitscore
                        }
                        hit_list.append(param)
            return hit_list

        except Exception as error:
            return None


    # Arguement for this method is a list of pdb file IDS
    # The return value is the list of filenames that have been downloaded
    # If return value is False, then file download failed
    def download_pdb_files(self,hit_list):
        try:
            downloader = PDBList()
            id_list = []
            for item in hit_list:
                item_id = (item['id'].split("|"))[1]
                id_list.append(item_id)
            download_result = downloader.download_pdb_files(pdb_codes=id_list, pdir= self.hits_directory)
            return download_result

        except Exception as error:
            return None


    def check_for_missing_residues(self):
        try:
            counter = 0
            missing_residue = True
            if os.path.exists(self.hits_directory):
                print(f"[*] Checking pdb file header for missing residues")
                for path,directories,files in os.walk(self.hits_directory):
                    for pdb_file in files:
                        counter += 1
                        pdb_file_path = os.path.join(self.hits_directory,pdb_file)
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


