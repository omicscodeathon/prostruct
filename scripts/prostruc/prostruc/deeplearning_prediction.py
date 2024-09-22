import requests
from Bio.PDB import PDBParser,PDBIO
from io import StringIO
# Basic header definition
header = {
    'Content-Type': 'application/x-www-form-urlencoded',
}

# URL to the deep learning model
url = 'https://api.esmatlas.com/foldSequence/v1/pdb'



def predict_with_deep_learning(job_name,sequence):
    try:
        parser = PDBParser()
        response = requests.post(url=url,headers=header, data=sequence, verify=False)
        prediction = StringIO(response.content.decode('utf-8'))
        predictions = parser.get_structure("validation",prediction)
        io=PDBIO()
        io.set_structure(predictions)
        io.save(f"{job_name}_validation_model")

    except Exception as error:
        print(error)
        return "ERROR"

# SEQ = """
# MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
# """
# prediction = predict_with_deep_learning(sequence=SEQ)
# # print(type(prediction))
# print(prediction)
# # print(type(prediction))
