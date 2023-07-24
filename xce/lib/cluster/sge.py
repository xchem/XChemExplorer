import os

LOGDIR = os.environ["XCE_LOGDIR"]
PROJECT_NAME = "labxchem"
QUEUE = "medium.q"


def submit_cluster_job(
    name,
    file,
    resources=None,
    parallel_environment=None,
    outfile=None,
    errfile=None,
    tasks=None,
    concurrent=None,
):
    base_command = "qsub -P {} -q {} -N {}".format(PROJECT_NAME, QUEUE, name)
    resource_arg = "-l {}".format(resources) if resources is not None else ""
    parallel_environment_arg = (
        "-pe {}".format(parallel_environment)
        if parallel_environment is not None
        else ""
    )
    outfile_arg = "-o {}".format(outfile) if outfile is not None else ""
    errfile_arg = "-e {}".format(errfile) if errfile is not None else ""
    tasks_arg = "-t {}".format(tasks) if tasks is not None else ""
    concurrent_arg = "-tc {}".format(concurrent) if concurrent is not None else ""
    os.system(
        " ".join(
            [
                base_command,
                resource_arg,
                parallel_environment_arg,
                outfile_arg,
                errfile_arg,
                tasks_arg,
                concurrent_arg,
                file,
            ]
        )
    )
