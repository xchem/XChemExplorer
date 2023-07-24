import os

LOGDIR = os.environ["XCE_LOGDIR"]
PROJECT_NAME = "labxchem"
QUEUE = "medium.q"


def submit_cluster_job(name, file, resources=""):
    os.system(
        "qsub -P {} -q {} -N {} -l {} {}".format(
            PROJECT_NAME, QUEUE, name, resources, file
        )
    )
