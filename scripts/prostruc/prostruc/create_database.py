import sqlite3

def create_database():
    try:
        connection = sqlite3.connect("prostruc_job_3.db")
        cursor = connection.cursor()

        cursor.execute(
            '''CREATE TABLE IF NOT EXISTS job (
                        job_id TEXT PRIMARY KEY,
                        job_name TEXT NOT NULL,
                        sequence TEXT,
                        email TEXT NOT NULL,
                        status TEXT NOT NULL,
                        current_operation TEXT NOT NULL
                    )'''
        )


        connection.commit()
        connection.close()

        return True

    except Exception as error:
        return False

create_database()
