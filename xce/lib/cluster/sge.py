import os
import subprocess
from datetime import datetime
from xce.lib.XChemLog import updateLog

PROJECT_NAME = "labxchem"
QUEUE = "medium.q"


def submit_cluster_job(
    name,
    file,
    xce_logfile,
    exclusive=False,
    memory=None,
    parallel_environment=None,
    outfile=None,
    errfile=None,
    tasks=None,
    concurrent=None,
):
    base_command = "qsub -P {} -q {} -N {}".format(PROJECT_NAME, QUEUE, name)
    resource_params = ",".join(
        resource
        for resource in [
            "exclusive" if exclusive else None,
            "m_mem_free={}".format(memory) if memory is not None else None,
        ]
        if resource is not None
    )
    resource_arg = "-l {}".format(resource_params) if resource_params else None
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


def query_running_jobs():
    qstat = subprocess.Popen(["qstat -u $user"], stdout=subprocess.PIPE)

    jobs = []
    for line in qstat.stdout.readlines():
        fields = line.split()

        if len(fields) < 7:
            continue

        job_id = fields[0]
        job_name = fields[2]
        job_status = fields[4]

        start_time = datetime.strptime(
            "{} {}".format(fields[5], fields[6]), "%m/%d/%Y %H:%M:%S"
        )
        run_time = datetime.now() - start_time
        jobs.append((job_id, job_name, job_status, run_time))

    return jobs
