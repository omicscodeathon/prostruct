from target import Target
from alignment import blast_search
from template_seq import Template
from analyser import basic_analysis

SEQ = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"


target_protein = Target(job_name='cap3',file_path='fasta file path')
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
else:
    print(analysis_result)


print('[*] Performing blast search')
blast_result = blast_search(sequence=sequence,job_name=job_name)
if blast_result == None:
    print("!!! Blast search operation requires internet")
else:
    print(f"[*] Blast search took {blast_result['duration']} seconds")

    template_obj = Template(file_system_info=file_system_info)

    print("[*] Checking for similar sequences")
    similar_sequences = template_obj.check_sequence_similarity(e_value_threshold=1e-50)
    if similar_sequences == None:
        print("!!! Similarity check failed")
        
    else:
        print('[*] Preparing to download files')
        download_list = template_obj.download_pdb_files(hit_list=similar_sequences)
        template_obj.check_for_missing_residues()
