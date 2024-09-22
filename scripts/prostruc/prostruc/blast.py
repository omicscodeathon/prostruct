from Bio.Blast import NCBIWWW,NCBIXML
import urllib.error
import time
from datetime import date
import os
# from Bio.PDB import MMCIFParser


amino_acid_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
                    'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
                    'TRP': 'W', 'TYR': 'Y'}


def blast_search_xml(sequence,job_name):
    try:
        blast_search_result = f"{job_name}_{date.today()}.xml"
        blast_search_output_path = os.path.join('blast_search_output',f"{job_name}_{date.today()}.xml")
        path = f"{job_name}_{date.today()}.xml"
        start_time = time.time()
        search_result = NCBIWWW.qblast(program='blastp', database='pdb', sequence=sequence)
        end_time = time.time()
        search_duration = end_time - start_time

        # print(f"[*] Saving search result as {self.job_name}_{date.today()}.xml\n")
        with open(blast_search_output_path,"w") as file_handler:
            result = search_result.read()
            file_handler.write(result)
            file_handler.close()

        result = {
            'result': search_result,
            'duration': search_duration
        }

        return result


    except urllib.error.URLError as error:
        return None

