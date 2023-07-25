import os
from xce.lib.XChemLog import updateLog

LOGDIR = os.environ["XCE_LOGDIR"]
PROJECT_NAME = "labxchem"
QUEUE = "medium.q"


def submit_cluster_job(
    name,
    file,
    xce_logfile,
    resources=None,
    parallel_environment=None,
    outfile=None,
    errfile=None,
    tasks=None,
    concurrent=None,
):
    base_command = "qsub -P {} -q {} -N {}".format(PROJECT_NAME, QUEUE, name)
    resource_arg = "-l {}".format(resources) if resources is not None else None
    parallel_environment_arg = (
        "-pe {}".format(parallel_environment)
        if parallel_environment is not None
        else None
    )
    outfile_arg = "-o {}".format(outfile) if outfile is not None else None
    errfile_arg = "-e {}".format(errfile) if errfile is not None else None
    tasks_arg = "-t {}".format(tasks) if tasks is not None else None
    concurrent_arg = "-tc {}".format(concurrent) if concurrent is not None else None
    command = " ".join(
        part
        for part in [
            base_command,
            resource_arg,
            parallel_environment_arg,
            outfile_arg,
            errfile_arg,
            tasks_arg,
            concurrent_arg,
            file,
        ]
        if part is not None
    )
    logfile = updateLog(xce_logfile)
    logfile.insert(
        "Submitting job, '{}', to SGE with command: {}".format(name, command)
    )
    os.system(command)
