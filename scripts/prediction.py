import requests

# Basic header definition
header = {
    'Content-Type': 'application/x-www-form-urlencoded',
}

# URL to the deep learning model
url = 'https://api.esmatlas.com/foldSequence/v1/pdb'


def predict_with_deep_learning(sequence):
    try:
         response = requests.post(url=url,headers=header, data=sequence, verify=False)
         prediction = response.content.decode('utf-8')
         return prediction

    except Exception as error:
        print(error)
        return "ERROR"
