import os
import re
import time
from Bio.PDB import PDBParser, Superimposer
import requests
import subprocess


def calculate_rmsd_old(target_model, validation_model_structure):
    parser = PDBParser(QUIET=True)
    target_structure = parser.get_structure('target', target_model)

    target_atoms = []
    validation_atoms = []

    # Extract only the alpha carbon atoms (CA)
    for model in target_structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    target_atoms.append(residue['CA'])

    for model in validation_model_structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    validation_atoms.append(residue['CA'])

    # if len(target_atoms) != len(validation_atoms):
    #     raise PDBException("The number of alpha carbon atoms in the target and validation models do not match.")

    calculator = Superimposer()
    calculator.set_atoms(target_atoms, validation_atoms)
    rmsd = calculator.rms

    return rmsd


def calculate_qmean(target_model,email,job_name,method='qmeandisco'):
    SERVER = "https://swissmodel.expasy.org/qmean/submit"
    data = {
        'email': email,
        'project_name': job_name,
        'method': method
    }
    response = requests.post(SERVER, files={'structure': open(target_model,'rb')},data=data)

    if response.status_code == 200:
       job_info = response.json()
       return job_info['results_json']
    # else:
    #     raise Exception("Error in fetching QMEAN score")


def get_qmean_result(results_url):
    # global qmean
    try:
        response = requests.get(results_url)
        result = dict(response.json())
        if result['status'] == "RUNNING":
            print("waiting for result")
            time.sleep(60)
            response = requests.get(results_url)
            result = dict(response.json())
            qmean = result['models']['model_001']['scores']['global_scores']['avg_local_score']
        else:
            qmean = result['models']['model_001']['scores']['global_scores']['avg_local_score']  #["qmean6_norm_score"]

        return qmean
        # else:
        #     return qmean

    except Exception as error:
        print(error)
        return None


def calculate_tm_score(target_model, validation_model):
    host_directory = os.getcwd()

    tm_command = f"docker run -v {host_directory}:/validation edraizen/tmalign /validation/{target_model} /validation/{validation_model}"

    result = subprocess.run(tm_command,capture_output=True, text=True)
    # print(result)
    output = result.stdout

    # Extract RMSD
    rmsd_match = re.search(r"RMSD=\s+(\d+\.\d+)", output)
    rmsd = float(rmsd_match.group(1)) if rmsd_match else None

    # Extract TM-scores
    tm_scores = re.findall(r"TM-score=\s+(\d+\.\d+)", output)
    tm_score_chain1 = float(tm_scores[0]) if len(tm_scores) > 0 else None
    tm_score_chain2 = float(tm_scores[1]) if len(tm_scores) > 1 else None

    #return rmsd, tm_score_chain1, tm_score_chain2
    return {
        "rmsd":rmsd,
        "tm_score_chain1": tm_score_chain1,
        "tm_score_chain2": tm_score_chain2
    }

# print(calculate_tm_score(target_model="1efe.cif_pdb.pdb_output_model.pdb",validation_model="pro4_validation_model"))

# validation_model = predict_with_deep_learning(sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
# result = calculate_rmsd(target_model="modeled_structures/2kqp.cif_pdb.pdb_output_model.pdb",validation_model_structure=validation_model)
# print(result)
