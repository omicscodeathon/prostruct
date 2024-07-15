from Bio.Blast import NCBIWWW
from Bio import Align, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo
import urllib.error
import pdbreader
import time
from datetime import date
import os
from Bio.PDB import MMCIFParser

amino_acid_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
                    'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
                    'TRP': 'W', 'TYR': 'Y'}

def blast_search(sequence,job_name):
    try:
        blast_search_result = f"{job_name}_{date.today()}.xml"
        blast_search_output_path = os.path.join('blast_search_output',f"{job_name}_{date.today()}.xml")
        path = f"{job_name}_{date.today()}.xml"
        # print("[*] Searching database for templates")
        start_time = time.time()
        search_result = NCBIWWW.qblast(program='blastp', database='pdb', sequence=sequence)
        end_time = time.time()
        search_duration = end_time - start_time
        # print(f"[*] Search done, operation took {search_duration} seconds")

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


def get_sequence_from_cif(file):
    parser = MMCIFParser()
    structure = parser.get_structure('structure',file)
    sequence = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    sequence.append(residue.resname)

    sequence = ''.join([amino_acid_dict.get(aa, 'X') for aa in sequence])
    return sequence



def align(target, directory):
    alignment_records = []
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    # aligner.substitution_matrix = MatrixInfo.blosum65
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5


    for path, directories, files in os.walk(directory):
        for pdb_file in files:
            file_path = os.path.join(path, pdb_file)
            sequence = get_sequence_from_cif(file_path)
            alignments = aligner.align(target,sequence)
            for alignment in alignments:
                target_seq = (''.join([target[i:j] if i != j else '-' * (j - i) for i, j in alignment.aligned[0]]))
                template_seq = (''.join([sequence[i:j] if i != j else '-' * (j - i) for i, j in alignment.aligned[1]]))

                target_record = SeqRecord(Seq(target_seq),id='target')
                template_record = SeqRecord(Seq(template_seq), id=f"{pdb_file}")

                alignment_records.append(target_record)
                alignment_records.append(template_record)

                alignments_obj = Align.MultipleSeqAlignment(alignment_records)

                with open(f"{pdb_file}.aln",'w') as handler:
                    AlignIO.write(alignments_obj,handler,'clustal')





# SEQ = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
# align(target=SEQ,directory='cap1_pdb_files')



