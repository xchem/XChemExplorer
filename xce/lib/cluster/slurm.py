import os
import ssl
import json
import httplib
from datetime import datetime
from xce.lib.XChemLog import updateLog
from uuid import uuid4

CLUSTER_USER = os.getlogin()
CLUSTER_TOKEN = os.environ.get("SLURM_JWT", "")
CLUSTER_HOST = "slurm-rest.diamond.ac.uk"
CLUSTER_PORT = 8443
CLUSTER_PARTITION = "cs04r"
CLUSTER_ACCOUNT = "labxchem"

HEADERS = {
    "Content-Type": "application/json",
    "X-SLURM-USER-NAME": CLUSTER_USER,
    "X-SLURM-USER-TOKEN": CLUSTER_TOKEN,
}


def submit_cluster_job(
    name, file, xce_logfile, array=None, exclusive=False, memory=None, tasks=None
):
    with open(file) as script_file:
        script = "\n".join(script_file.readlines())
    payload = dict(
        script=script,
        job=dict(
            partition=CLUSTER_PARTITION,
            name=str(name),
            account=CLUSTER_ACCOUNT,
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
    connection.request("POST", "/slurm/v0.0.38/job/submit", body=body, headers=HEADERS)
    response = connection.getresponse().read()
    logfile.insert("Got response: {}".format(response))


def query_running_jobs(xce_logfile):
    connection = httplib.HTTPSConnection(
        CLUSTER_HOST,
        CLUSTER_PORT,
        context=ssl._create_unverified_context(),
    )
    connection.request("GET", "/slurm/v0.0.38/jobs", headers=HEADERS)
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
