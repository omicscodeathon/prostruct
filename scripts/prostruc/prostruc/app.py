import os
import urllib3
from prostruc.target import Target
from prostruc.blast import blast_search_xml
from prostruc.template import Template
from prostruc.analyser import basic_analysis
from prostruc.alignment import Alignment
from prostruc.model import model
from prostruc.deeplearning_prediction import predict_with_deep_learning
from prostruc.initial_checks import check_for_internet, check_for_docker, docker_check_image
from prostruc.validation import validate
from prostruc.cleanup import clean
import time
import warnings
from Bio import BiopythonWarning, BiopythonDeprecationWarning
import argparse


# Ignore Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', urllib3.exceptions.InsecureRequestWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore',BiopythonDeprecationWarning)
warnings.filterwarnings("ignore")


def main():
    # Setting up the command-line arguement parser instance
    parser = argparse.ArgumentParser(description="ProStruc: A Comprehensive Command-line Tool for Protein Structure Prediction and Validation")

    # Defining arguments
    parser.add_argument("--fasta_file", metavar="FILE",type=str, required=False, help="Path to the target fasta file")
    parser.add_argument("--sequence", metavar= "SEQUENCE", type=str, required=False, help="Amino acid Sequence")
    parser.add_argument("--job_name", metavar= "JOB NAME", type=str,help="Use an appropriate job name")
    parser.add_argument("--email", metavar="example@domain.com", type=str,help="Use a valid email address")

    args = parser.parse_args()


    ################################### Fasta File ###############################################
    if args.fasta_file and args.job_name and args.email:
        target_protein = Target(job_name=args.job_name,file_path=args.fasta_file)
        path = target_protein.get_blast_search_result_path()
        file_system_info = target_protein.get_file_manager_info()

        print("[*] Checking for internet connectivity and docker engine")
        internet_state = check_for_internet()
        if internet_state == False:
            print("Quitting, Internet connectivity is required")
            exit(1)

        docker_state = check_for_docker()
        if docker_state == False:
            print("Quitting, Docker is either not running or not installed")
            exit(1)

        else:
            print("[*] Docker engine is installed")

            print(f"[*] Checking ProMod3 image :: registry.scicore.unibas.ch/schwede/promod3:latest")
            promod3_image_state = docker_check_image(image_tag="registry.scicore.unibas.ch/schwede/promod3:latest")

            print(f"[*] Checking TMAlign image :: edraizen/tmalign:latest")
            tmalign_image_state = docker_check_image(image_tag="edraizen/tmalign:latest")

        print(f'[*] Analysing the sequence')
        analysis_result = basic_analysis(sequence=args.sequence)
        if analysis_result == None:
            print('ERROR: Analysis failed')
            clean(directory="alignment_files")
            clean(directory="blast_search_output")
            clean(directory="input_fasta_files")
            clean(directory="obsolete")
            clean(directory="promod3_input_pdb_files")
            clean(directory="modeled_structures")
            exit(1)

        elif analysis_result["Length"][0] > 400:
            print("ERROR: Maximum sequence length is 400")
            clean(directory="alignment_files")
            clean(directory="blast_search_output")
            clean(directory="input_fasta_files")
            clean(directory="obsolete")
            clean(directory="promod3_input_pdb_files")
            clean(directory="modeled_structures")
            exit(1)

        else:
            print("[*] Generating validation model using ESMFold")
            validation_model = predict_with_deep_learning(sequence=args.sequence,job_name=args.job_name)
            if validation_model == "ERROR":
                print("ERROR: Internet access is required")
                clean(directory="alignment_files")
                clean(directory="blast_search_output")
                clean(directory="input_fasta_files")
                clean(directory="obsolete")
                clean(directory="promod3_input_pdb_files")
                exit(code=1)

            else:
                print('[*] Performing blast search')
                blast_result = blast_search_xml(sequence=args.sequence,job_name=args.job_name)
                print(f"[*] Blast search took {blast_result['duration']} seconds")
                if blast_result == None:
                    print("ERROR: Blast search operation requires internet")
                    clean(directory="alignment_files")
                    clean(directory="blast_search_output")
                    clean(directory="input_fasta_files")
                    clean(directory="obsolete")
                    clean(directory="promod3_input_pdb_files")
                    clean(directory="modeled_structures")
                    os.remove(f"{args.job_name}_validation_model")
                    exit(code=1)

                else:
                    print("[*] Performing similarity check")
                    template_obj = Template(file_system_info=file_system_info)
                    similar_sequences = template_obj.check_sequence_similarity(e_value_threshold=0.01,min_identity=30)
                    if similar_sequences is None or len(similar_sequences) == 0:
                        print("INFO: No similar sequences found. Exiting pipeline.")
                        clean(directory="alignment_files")
                        clean(directory="blast_search_output")
                        clean(directory="input_fasta_files")
                        clean(directory="obsolete")
                        clean(directory="promod3_input_pdb_files")
                        clean(directory="modeled_structures")
                        os.remove(f"{args.job_name}_validation_model")
                        exit(1)

                    else:
                        print('[*] Preparing to download files')
                        download_list = template_obj.download_pdb_files(match_list=similar_sequences)
                        template_obj.check_for_missing_residues()

                        print("[*] Performing sequence alignment")
                        sequence_aligner = Alignment(target=args.sequence,directory=f"{args.job_name}_pdb_files")
                        sequence_aligner.align_pairwise()

                        print("[*] Beginning modelling")
                        print("[*] Modelling Engine = docker image of ProMod3")
                        start_time = time.time()
                        model(alignment_directory="alignment_files", structure_directory="promod3_input_pdb_files")
                        end_time = time.time()
                        print(f"[*] Modelling took {end_time - start_time} seconds")

                        print("[*] Validating predictions")
                        validation_state = validate(modeled_structures_directory="modeled_structures",email=args.email,job_name=args.job_name)
                        if validation_state['state'] == True:
                            print("[*] Validation successful!!")
                            clean(directory="alignment_files")
                            clean(directory="blast_search_output")
                            clean(directory="input_fasta_files")
                            clean(directory="obsolete")
                            clean(directory="promod3_input_pdb_files")
                            os.remove(f"{args.job_name}_validation_model")


    ############################# Inline Sequence  ##################################################

    elif args.sequence and args.job_name and args.email:
        target_protein = Target(job_name=args.job_name,seq=args.sequence)
        path = target_protein.get_blast_search_result_path()
        file_system_info = target_protein.get_file_manager_info()

        print("[*] Checking for internet connectivity and docker engine")
        internet_state = check_for_internet()
        if internet_state == False:
            print("Quitting, Internet connectivity is required")
            exit(1)

        docker_state = check_for_docker()
        if docker_state == False:
            print("Quitting, Docker is either not running or not installed")
            exit(1)

        else:
            print("[*] Docker engine is installed")

            print(f"[*] Checking ProMod3 image :: registry.scicore.unibas.ch/schwede/promod3:latest")
            promod3_image_state = docker_check_image(image_tag="registry.scicore.unibas.ch/schwede/promod3:latest")

            print(f"[*] Checking TMAlign image :: edraizen/tmalign:latest")
            tmalign_image_state = docker_check_image(image_tag="edraizen/tmalign:latest")

        print(f'[*] Analysing the sequence')
        analysis_result = basic_analysis(sequence=args.sequence)
        if analysis_result == None:
            print('ERROR: Analysis failed')
            clean(directory="alignment_files")
            clean(directory="blast_search_output")
            clean(directory="input_fasta_files")
            clean(directory="obsolete")
            clean(directory="promod3_input_pdb_files")
            clean(directory="modeled_structures")
            exit(1)

        elif analysis_result["Length"][0] > 400:
            print("ERROR: Maximum sequence length is 400")
            clean(directory="alignment_files")
            clean(directory="blast_search_output")
            clean(directory="input_fasta_files")
            clean(directory="obsolete")
            clean(directory="promod3_input_pdb_files")
            clean(directory="modeled_structures")
            exit(1)

        else:
            print("[*] Generating validation model using ESMFold")
            validation_model = predict_with_deep_learning(sequence=args.sequence,job_name=args.job_name)
            if validation_model == "ERROR":
                print("ERROR: Internet access is required")
                clean(directory="alignment_files")
                clean(directory="blast_search_output")
                clean(directory="input_fasta_files")
                clean(directory="obsolete")
                clean(directory="promod3_input_pdb_files")
                exit(code=1)

            else:
                print('[*] Performing blast search')
                blast_result = blast_search_xml(sequence=args.sequence,job_name=args.job_name)
                print(f"[*] Blast search took {blast_result['duration']} seconds")
                if blast_result == None:
                    print("ERROR: Blast search operation requires internet")
                    clean(directory="alignment_files")
                    clean(directory="blast_search_output")
                    clean(directory="input_fasta_files")
                    clean(directory="obsolete")
                    clean(directory="promod3_input_pdb_files")
                    clean(directory="modeled_structures")
                    os.remove(f"{args.job_name}_validation_model")
                    exit(code=1)

                else:
                    print("[*] Performing similarity check")
                    template_obj = Template(file_system_info=file_system_info)
                    similar_sequences = template_obj.check_sequence_similarity(e_value_threshold=0.01,min_identity=30)
                    if similar_sequences is None or len(similar_sequences) == 0:
                        print("INFO: No similar sequences found. Exiting pipeline.")
                        clean(directory="alignment_files")
                        clean(directory="blast_search_output")
                        clean(directory="input_fasta_files")
                        clean(directory="obsolete")
                        clean(directory="promod3_input_pdb_files")
                        clean(directory="modeled_structures")
                        os.remove(f"{args.job_name}_validation_model")
                        exit(1)

                    else:
                        print('[*] Preparing to download files')
                        download_list = template_obj.download_pdb_files(match_list=similar_sequences)
                        template_obj.check_for_missing_residues()

                        print("[*] Performing sequence alignment")
                        sequence_aligner = Alignment(target=args.sequence,directory=f"{args.job_name}_pdb_files")
                        sequence_aligner.align_pairwise()

                        print("[*] Beginning modelling")
                        print("[*] Modelling Engine = docker image of ProMod3")
                        start_time = time.time()
                        model(alignment_directory="alignment_files", structure_directory="promod3_input_pdb_files")
                        end_time = time.time()
                        print(f"[*] Modelling took {end_time - start_time} seconds")

                        print("[*] Validating predictions")
                        validation_state = validate(modeled_structures_directory="modeled_structures",email=args.email,job_name=args.job_name)
                        if validation_state['state'] == True:
                            print("[*] Validation successful!!")
                            clean(directory="alignment_files")
                            clean(directory="blast_search_output")
                            clean(directory="input_fasta_files")
                            clean(directory="obsolete")
                            clean(directory="promod3_input_pdb_files")
                            os.remove(f"{args.job_name}_validation_model")

    else:
        print("ERROR: Either fasta file or sequence is required. If one of these are provided, then you are missing either job name or email")
        exit(1)


# if __name__ =="__main__":
#     main()
