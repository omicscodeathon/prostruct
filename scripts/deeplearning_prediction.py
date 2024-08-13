import requests
from Bio.PDB import PDBParser
from io import StringIO

# Basic header definition
header = {
    'Content-Type': 'application/x-www-form-urlencoded',
}

# URL to the deep learning engine
url = 'https://api.esmatlas.com/foldSequence/v1/pdb'


def predict_with_deep_learning(sequence):
    try:
        parser = PDBParser()
        response = requests.post(url=url,headers=header, data=sequence, verify=False)
        prediction = StringIO(response.content.decode('utf-8'))
        predictions = parser.get_structure("validation",prediction)
        return predictions
    except Exception as error:
        print(error)
        return "ERROR"

