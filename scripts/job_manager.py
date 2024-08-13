import uuid
import time
import sqlite3

def get_database():
    connection = sqlite3.connect("prostruc_jobs.db")
    connection.row_factory = sqlite3.Row
    return connection

def submit_job(job_id,job_name,sequence,email,status):
    connection = get_database()
    connection.execute(
        'INSERT INTO jobs (job_id, job_name, sequence, email, status) VALUES (?, ?, ?, ?, ?)',
        (job_id, job_name, sequence, email, status)
    )
    connection.commit()
    connection.close()


def update_job_status(job_id,status):
    connection = get_database()
    connection.execute('UPDATE jobs SET status = ? WHERE job_id = ?',(status, job_id))
    connection.commit()
    connection.close()


def retrieve_job(job_id):
    connection = get_database()
    job_details = connection.execute('SELECT * FROM jobs WHERE job_id = ?', (job_id,)).fetchone()
    connection.close()
    return job_details
