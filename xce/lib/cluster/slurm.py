import os
import ssl
import json
import httplib
import paramiko
import time
import traceback
from PyQt4 import QtGui
from datetime import datetime
from xce.lib.XChemLog import updateLog
from uuid import uuid4

CLUSTER_BASTION = "wilson.diamond.ac.uk"
CLUSTER_USER = os.environ.get("CLUSTER_USER", os.getlogin())
CLUSTER_HOST = "slurm-rest.diamond.ac.uk"
CLUSTER_PORT = 8443
CLUSTER_PARTITION = "cs05r"

TOKEN = None
TOKEN_EXPIRY = None


def get_token(error=None):
    global TOKEN
    global TOKEN_EXPIRY

    if TOKEN is None or TOKEN_EXPIRY is None or TOKEN_EXPIRY < time.clock() + 60:
        password_prompt = error + "\n" + "Password:" if error else "Password:"
        password, ok = QtGui.QInputDialog.getText(
            None, "SLURM Authentication", password_prompt
        )
        if not ok:
            return None

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.load_system_host_keys()
        try:
            ssh.connect(CLUSTER_BASTION, username=CLUSTER_USER, password=str(password))
        except paramiko.AuthenticationException:
            print(traceback.format_exc())
            return get_token(error="SSH Authentication Failed")
        stdin, stdout, stderr = ssh.exec_command("scontrol token lifespan=3600")
        if stdout.channel.recv_exit_status() != 0:
            return get_token(error="Token Acquisition Failed")
        TOKEN = stdout.next().split("=")[1].strip()
        TOKEN_EXPIRY = time.clock() + 3600
    return TOKEN


def construct_headers(token):
    return {
        "Content-Type": "application/json",
        "X-SLURM-USER-NAME": CLUSTER_USER,
        "X-SLURM-USER-TOKEN": token,
    }


def submit_cluster_job(
    name, file, xce_logfile, token, array=None, exclusive=False, memory=None, tasks=None
):
    with open(file) as script_file:
        script = "\n".join(script_file.readlines())
    payload = dict(
        script=script,
        job=dict(
            partition=CLUSTER_PARTITION,
            name=str(name),
            account=CLUSTER_USER,
            environment=dict(PLACE="HOLDER"),
            standard_output=os.path.join(os.getcwd(), "{}.stdout".format(name)),
            standard_error=os.path.join(os.getcwd(), "{}.stderr".format(name)),
        ),
    )
    if array is not None:
        payload["job"]["array"] = array
    if exclusive is True:
        payload["job"]["exclusive"] = "mcs"
        payload["job"]["mcs_label"] = str(uuid4())
    if memory is not None:
        payload["job"]["memory_per_node"]["set"] = True
        payload["job"]["memory_per_node"]["number"] = memory
    if tasks is not None:
        payload["job"]["tasks_per_node"]["set"] = True
        payload["job"]["tasks_per_node"]["number"] = tasks
    body = json.dumps(payload)
    logfile = updateLog(xce_logfile)
    logfile.insert("Submitting job, '{}', to Slurm with body: {}".format(name, body))
    connection = httplib.HTTPSConnection(
        CLUSTER_HOST, CLUSTER_PORT, context=ssl._create_unverified_context()
    )
    connection.request(
        "POST", "/slurm/v0.0.38/job/submit", body=body, headers=construct_headers(token)
    )
    response = connection.getresponse().read()
    logfile.insert("Got response: {}".format(response))


def query_running_jobs(xce_logfile, token):
    connection = httplib.HTTPSConnection(
        CLUSTER_HOST,
        CLUSTER_PORT,
        context=ssl._create_unverified_context(),
    )
    connection.request("GET", "/slurm/v0.0.38/jobs", headers=construct_headers(token))
    response = connection.getresponse()
    response_body = response.read()

    if response.status != 200:
        logifle = updateLog(xce_logfile)
        logifle.insert("Got response: {}".format(response_body))

    jobs = []
    for job in json.loads(response_body)["jobs"]:
        if job["user_name"] != CLUSTER_USER:
            continue

        job_id = job["job_id"]
        job_name = job["name"]
        job_status = job["job_state"]

        start_time = datetime.utcfromtimestamp(job["start_time"])
        run_time = datetime.now() - start_time

        jobs.append((job_id, job_name, job_status, run_time))

    return jobs
