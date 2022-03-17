"""

Tools for running jobs via Slurm.

"""

import subprocess, sys

SLURM_TEMPLATE="""#!/bin/sh
#SBATCH --nodes=$NODES
#SBATCH --ntasks-per-node=$TASKS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --output=$JOBNAME.log
#SBATCH --error=$JOBNAME.err

$CMD
"""

def submitJob(cmd, jobName, dependentJobID = None, nodes = 1, tasks = 1, mem = 8000, time = "12:00:00"):
    """Submit a command to run on a node via Slurm.

    Args:
        cmd (str): The command to run.
        jobName (str): Name of the job (used for e.g. log files).
        dependentJobID (int): If this job depends on a previous job completing sucessfully, give the
            job ID number here.
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run per node.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.

    Returns:
        The ID number of the submitted Slurm job.

    Note:
        The job is submitted from the current working directory. So you *might*
        want to use absolute paths here.

    """

    script=SLURM_TEMPLATE
    script=script.replace("$CMD", cmd)
    script=script.replace("$JOBNAME", jobName)
    script=script.replace("$NODES", str(nodes))
    script=script.replace("$TASKS", str(tasks))
    script=script.replace("$MEM", str(mem))
    script=script.replace("$TIME", str(time))

    fileName=jobName+".slurm"
    with open(fileName, "w") as outFile:
        outFile.write(script)

    args=['sbatch']
    if dependentJobID is not None:
        if type(dependentJobID) != int:
            raise Exception("dependentJobID must be an int")
        args=args+['-d', "afterok:%s" % (str(dependentJobID))]
    args=args+[fileName]

    process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise Exception("Non-zero return code when submitting job %s" % (jobName))
    assert(process.stdout[:19] == "Submitted batch job")
    jobID=int(process.stdout.split("Submitted batch job")[-1])

    return jobID

