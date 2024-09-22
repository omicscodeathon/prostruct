from Bio import SeqIO
from Bio.Seq import Seq
import os
from datetime import date

class Target:
    # Creates filesystem structure
    # Extracts sequence as text input or from a fasta file
    def __init__(self,job_name,seq=None,file_path=None):
        self.job_name = job_name
        self.file_mamanger_info = self.filemanager()

        if seq == None:
            self.sequence = self.sequence_from_fasta(file_path)

        if file_path == None:
            self.sequence = seq



    # Returns information about the file management system
    # If return value is None, then operation failed
    def filemanager(self):
        try:
            self.create_directory(directory_name='blast_search_output')
            self.create_directory(directory_name='input_fasta_files')
            self.create_directory(directory_name='alignment_files')
            self.create_directory(directory_name="promod3_input_pdb_files")
            self.create_directory(directory_name="modeled_structures")
            # self.create_directory(directory_name=self.job_name)

            self.templates_directory = f"{self.job_name}_pdb_files"
            self.alignment_files_directory = "alignment_files"
            self.blast_search_result = f"{self.job_name}_{date.today()}.xml"
            self.blast_search_output_path = os.path.join('blast_search_output',f"{self.job_name}_{date.today()}.xml")



            file_manager_info = {
                'templates directory': self.templates_directory,
                'blast search result': self.blast_search_result,
                'blast search output path': self.blast_search_output_path,
                'alignment files directory': self.alignment_files_directory
            }
            return file_manager_info

        except Exception as error:
            return None


    # Checks for and creates a directory
    def create_directory(self,directory_name):
        if os.path.exists(directory_name) == True:
            pass
        else:
            os.mkdir(directory_name)
            # print(f"{directory_name} created")



    # Extracts the sequence from a fasta file
    # If return value is None, then the operation failed
    def sequence_from_fasta(self,fasta_filepath):
        global seq
        try:
            if self.is_file_fasta(filepath=fasta_filepath) == False:
                fasta_filename = self.convert_to_fasta(filepath=fasta_filepath)
                sequence_obj = SeqIO.parse(fasta_filename,format='fasta')
                for record in sequence_obj:
                    seq = record.seq
                sequence = str(seq)


            else:
                sequence_obj = SeqIO.parse(fasta_filepath,format='fasta')
                sequence = sequence_obj.seq

            return sequence

        except Exception as error:
            print(f"!!! {error}:: Invalid file path")
            return None



    # Validates whether a given sequence file is in the fasta format
    # If return value is True, then the file is in fasta format
    # If return value is False, then the file is not in fasta format
    def is_file_fasta(self,filepath):
        try:
            with open(filepath,'r') as file_handler:
                sequence = file_handler.read()
                seq = sequence.seq
                return True

        except AttributeError:
            return False


    # Converts a given file into fasta file format
    # Stores the conversion in fasta_files directory
    # If return value is None, the operation failed
    def convert_to_fasta(self,filepath):
        try:
            new_file = f"{self.job_name}.fasta"
            new_file_path = os.path.join('input_fasta_files',f"{new_file}")
            with open(filepath,'r') as file_handler:
                seq = file_handler.read()

            with open(new_file_path,'w') as new_file_handler:
                new_file_handler.write(seq)
                # return SeqIO.SeqRecord(sequence).seq
            return new_file_path

        except Exception as error:
            return None


    # Returns the target sequences
    def get_sequence(self):
        return self.sequence

    def get_job_name(self):
        return self.job_name

    # Returns information about the file management structure
    def get_file_manager_info(self):
        return self.file_mamanger_info


    # Returns teh path to the blast search result
    def get_blast_search_result_path(self):
        return self.blast_search_output_path

#
# SEQ = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
# target = Target(job_name='new job',file_path='C:\\Users\\Wilson\\Downloads\\fasta.txt')
# result = target.convert_to_fasta(filepath='C:\\Users\\Wilson\\Downloads\\fasta.txt')
# # print(target.sequence_from_fasta(fasta_filepath='C:\\Users\\Wilson\\Downloads\\fasta.txt'))
# print(target.get_sequence())
# # print(target.analyse_target())
