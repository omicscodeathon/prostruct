import os
import requests
from Bio import PDB
import numpy as np
import time

def validate(modeled_structures_directory,validation_model_structure,email,job_name):
    print("[*] Validating models using RMSD")
    if os.path.exists(modeled_structures_directory):
        for path,directories,files in os.walk(modeled_structures_directory):
            for target_model in files:
                target_model_filepath = os.path.join(modeled_structures_directory,target_model)
                rmsd = calculate_rsmd(target_model=target_model_filepath, validation_model_structure=validation_model_structure)
                print(f"{target_model}: {rmsd}")
                if rmsd > 2:
                    os.remove(path=target_model_filepath)
                    print(f"{target_model} deleted")

        print("[*] Validating models using QMEAN")
        for path,directories,files in os.walk(modeled_structures_directory):
            for target_model in files:
                target_model_filepath = os.path.join(modeled_structures_directory,target_model)
                result_url = calculate_qmean(target_model=target_model_filepath,email=email,job_name=job_name)
                time.sleep(200)
                qmean = get_qmean_result(results_url=result_url)
                print(f"{target_model}: {qmean} :: {result_url}")
                if qmean < float(0.6):
                    os.remove(path=target_model_filepath)
                    print(f"{target_model} deleted")

        return True


def calculate_rsmd(target_model,validation_model_structure):
    global diff
    parser = PDB.PDBParser()
    target_structure = parser.get_structure(file=target_model,id='target')
    # print(type(target_structure))
    # validation_structure = parser.get_structure(file=validation_model, id='validation')

    target_atoms = []
    validation_atoms = []

    for model in target_structure:

        for chain in model:
            for residue in chain:
                for atom in residue:
                    target_atoms.append(atom)

    for model in validation_model_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    validation_atoms.append(atom)

    for target_atom,validation_atom in zip(target_atoms,validation_atoms):
        diff = np.array([target_atom.coord - validation_atom.coord])

    rmsd = np.sqrt((diff ** 2).sum() / len(target_atoms))
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


# validation_model = predict_with_deep_learning(sequence="FEWRWADIAAECERFLGPNGFGGVQISPPNDHIVLNNPWRPWWQRYQPIGYNLCSRSGSENELRDMITRCNNVGVNIYVDAVINHMCGAGGGEGTHSSCGSWFSAGRRDFPTVPYSHLDFNDNKCRTGSGDIENYGDSNQVRDCRLVGLLDLALEKEYVRGKVVDFMNKLIDMGVAGFRVDACKHMWP")
# validate(modeled_structures_directory="modeled_structures",validation_model_structure=validation_model,email="senawilson123@gmail.com",job_name='cap10')

# print(get_qmean_result(results_url="https://swissmodel.expasy.org/qmean/t48vbG.json"))
