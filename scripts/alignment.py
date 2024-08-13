# from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from blast import amino_acid_dict


class Alignment:
    def __init__(self,target, directory):
        self.target = target
        self.directory = directory

        # Setting up the alignment tool
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # Change to 'local' for local alignment
        self.aligner.match_score = 2.0
        self.aligner.mismatch_score = -1.0
        self.aligner.open_gap_score = -0.5
        self.aligner.extend_gap_score = -0.1


    def get_sequence_from_cif(self,file):
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


    def pairwise(self,template,id):
        try:
            seq_alignments = self.aligner.align(self.target,template)
            best_aligned = seq_alignments[0]

            target_aligned, template_aligned = best_aligned[0],best_aligned[1]
            target_record = SeqRecord(seq=Seq(target_aligned),id='target')
            template_record = SeqRecord(seq=Seq(template_aligned),id="template")

            path = os.path.join("alignment_files", f"{id}_aln.fasta")
            with open(path,'w') as file_handler:
                SeqIO.write([target_record,template_record],handle=file_handler,format='fasta')


        except Exception as error:
            print(error)
            return "Failed"


    # def pairwise(self,template,id):
    #     try:
    #         seq_alignments = pairwise2.align.globalxx(self.target,template)
    #         best_aligned = seq_alignments[0]
    #
    #         target_aligned, template_aligned = best_aligned[0],best_aligned[1]
    #         target_record = SeqRecord(seq=Seq(target_aligned),id='target')
    #         template_record = SeqRecord(seq=Seq(template_aligned),id="template")
    #
    #         path = os.path.join("alignment_files", f"{id}_aln.fasta")
    #         with open(path,'w') as file_handler:
    #             SeqIO.write([target_record,template_record],handle=file_handler,format='fasta')
    #
    #
    #     except Exception as error:
    #         print(error)
    #         return "Failed"



    def align_pairwise(self):
        for path, directories, files in os.walk(self.directory):
            for pdb_file in files:
                pdb_id = pdb_file.split('.')
                file_path = os.path.join(path, pdb_file)
                sequence = self.get_sequence_from_cif(file_path)
                alignment_state = self.pairwise(template=sequence, id=pdb_id[0])
                # print(f"alignment state: {alignment_state}")
                if alignment_state == "Failed":
                    print(f"!!! alignment for target and {pdb_file} failed")

