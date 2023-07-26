import os
import json
import httplib
from datetime import datetime
from xce.lib.XChemLog import updateLog

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


def submit_cluster_job(name, file, xce_logfile, array=None):
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
    body = json.dumps(payload)
    logfile = updateLog(xce_logfile)
    logfile.insert(
        "Submitting job, '{}', to Slurm with body: {} and headers: {}".format(
            name, body, HEADERS
        )
    )
    connection = httplib.HTTPSConnection(CLUSTER_HOST, CLUSTER_PORT)
    connection.request("POST", "/slurm/v0.0.38/job/submit", body=body, headers=HEADERS)
    response = connection.getresponse().read()
    logfile.insert("Got response: {}".format(response))


def query_running_jobs():
    connection = httplib.HTTPSConnection(CLUSTER_HOST, CLUSTER_PORT)
    connection.request("GET", "/slurm/v0.0.38/jobs", headers=HEADERS)
    response = connection.getresponse()
    response_body = response.read()

    if response.status != 200:
        print("Got response: {}".format(response_body))

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
