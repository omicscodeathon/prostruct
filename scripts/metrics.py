import time
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.PDBExceptions import PDBException
from bs4 import BeautifulSoup
import requests



def calculate_rmsd(target_model, validation_model_structure):
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

    if len(target_atoms) != len(validation_atoms):
        raise PDBException("The number of alpha carbon atoms in the target and validation models do not match.")

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
    else:
        raise Exception("Error in fetching QMEAN score")


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


