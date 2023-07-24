import os

LOGDIR = os.environ["XCE_LOGDIR"]
PROJECT_NAME = "labxchem"
QUEUE = "medium.q"


def submit_cluster_job(
    name, file, resources=None, parallel_environment=None, logfile=None, errfile=None
):
    base_command = "qsub -P {} -q {} -N {}".format(PROJECT_NAME, QUEUE, name)
    resource_arg = "-l {}".format(resources) if resources is not None else ""
    parallel_environment_arg = (
        "-pe {}".format(parallel_environment)
        if parallel_environment is not None
        else ""
    )
    logfile_arg = "-o {}".format(logfile) if logfile is not None else ""
    errfile_arg = "-e {}".format(errfile) if errfile is not None else ""
    os.system(
        " ".join(
            [
                base_command,
                resource_arg,
                parallel_environment_arg,
                logfile_arg,
                errfile_arg,
                file,
            ]
        )
    )
