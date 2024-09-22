import os
import shutil
def clean(directory):
    try:
        if os.path.exists(directory):
            print(f"Removing {directory}")
            shutil.rmtree(path=directory)

        else:
            print(f"{directory} doesn't exist")

    except Exception as error:
        print(error)
        return "ERROR"

