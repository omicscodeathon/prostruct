import os
import time
from metrics import calculate_rmsd, calculate_qmean, get_qmean_result


def validate(modeled_structures_directory,validation_model_structure,email,job_name):
    print("[*] Validating models using RMSD")
    if os.path.exists(modeled_structures_directory):
        for path,directories,files in os.walk(modeled_structures_directory):
            for target_model in files:
                target_model_filepath = os.path.join(modeled_structures_directory,target_model)
                rmsd = calculate_rmsd(target_model=target_model_filepath, validation_model_structure=validation_model_structure)
                print(f"{target_model}: {rmsd}")
                if rmsd > 2:
                    os.remove(path=target_model_filepath)
                    print(f"{target_model} deleted")


        print("[*] Validating models using QMEAN")
        for path,directories,files in os.walk(modeled_structures_directory):
            for target_model in files:
                target_model_filepath = os.path.join(modeled_structures_directory,target_model)
                result_url = calculate_qmean(target_model=target_model_filepath,email=email,job_name=job_name)
                time.sleep(150)
                qmean = get_qmean_result(results_url=result_url)
                print(f"{target_model}: {qmean} :: {result_url}")
                if qmean < float(0.6):
                    os.remove(path=target_model_filepath)
                    print(f"{target_model} deleted")

        return True


    except Exception as error:
        print(error)
        return None


