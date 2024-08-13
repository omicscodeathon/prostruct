import urllib3
from target import Target
from blast import blast_search_xml
from template import Template
from analyser import basic_analysis
from alignment import Alignment
from model import model
from deeplearning_prediction import predict_with_deep_learning
from validation import validate
from cleanup import clean
import time
import warnings
from Bio import BiopythonWarning, BiopythonDeprecationWarning
import argparse

# Supressing warnings
warnings.simplefilter('ignore', urllib3.exceptions.InsecureRequestWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore',BiopythonDeprecationWarning)


# Setting up the command-line arguement parser instance
parser = argparse.ArgumentParser(description="ProStruc terminal version")

# Defining arguments
parser.add_argument("--fasta_file", metavar="Path to the target fasta file",type=str, required=False)
parser.add_argument("--sequence", metavar="Amino acid Sequence", type=str, required=False)
parser.add_argument("--job_name", metavar="Use an appropriate job name", type=str)
parser.add_argument("--email", metavar="Use a valid email address", type=str)

args = parser.parse_args()

# Ignore Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)


# SEQ = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
# FILE = "C:\\Users\\Wilson\\Downloads\\fasta3.txt"

if args.fasta_file:
    target_protein = Target(job_name=args.job_name,file_path=args.fasta_file)
elif args.sequence:
    target_protein = Target(job_name=args.job_name,seq=args.sequence)
else:
    print("Either fasta file or sequence is required")
    exit(1)

# target_protein = Target(job_name=args.job_name,file_path=args.fasta_file)
sequence = target_protein.get_sequence()
job_name = target_protein.get_job_name()
print('[*] Initialization Done!!')

path = target_protein.get_blast_search_result_path()
file_system_info = target_protein.get_file_manager_info()
if file_system_info == None:
    print("!!! Operation failed")

print(f'[*] Analysing the sequence')
analysis_result = basic_analysis(sequence=sequence)
if analysis_result == None:
    print('!!! Analysis failed')
# else:
#     print(analysis_result)

print("[*] Generating validation model using ESMFold")
validation_model = predict_with_deep_learning(sequence=sequence)
if validation_model == "ERROR":
    print("!!! Internet access is required")
    exit(code=1)

print('[*] Performing blast search')
blast_result = blast_search_xml(sequence=sequence,job_name=job_name)
if blast_result == None:
    print("!!! Blast search operation requires internet")
else:
    print(f"[*] Blast search took {blast_result['duration']} seconds")

    template_obj = Template(file_system_info=file_system_info)

    print("[*] Checking for similar sequences")
    similar_sequences = template_obj.check_sequence_similarity(e_value_threshold=1e-50,min_identity=30)
    if similar_sequences == None:
        print("!!! Similarity check failed")

    else:
        print('[*] Preparing to download files')
        download_list = template_obj.download_pdb_files(match_list=similar_sequences)
        template_obj.check_for_missing_residues()

        print("[*] Performing sequence alignment")
        sequence_aligner = Alignment(target=sequence,directory=f"{job_name}_pdb_files")
        sequence_aligner.align_pairwise()

        print("[*] Beginning modelling")
        print("[*] Modelling Engine = docker image of ProMod3")
        start_time = time.time()
        model(alignment_directory="alignment_files", structure_directory="promod3_input_pdb_files")
        end_time = time.time()
        print(f"[*] Modelling took {end_time - start_time} seconds")


print("[*] Validating Structures")
validation_state = validate(modeled_structures_directory="modeled_structures", validation_model_structure=validation_model,email=args.email,job_name=job_name)
if validation_state == True:
    print("[*] Validation successful!!")

    # Deleting all temporary files and directories
    print("[*] Cleaning up")
    clean(directory="blast_search_output")
    clean(directory="input_fasta_files")
    clean(directory="obsolete")
    clean(directory="alignment_files")
    clean(directory=f"{job_name}_pdb_files")
    clean(directory="promod3_input_pdb_files")

    print("[***] Check modeled_structures directory for predicted structures")
