"""

Tools for running jobs via Slurm.

"""

import subprocess, sys

SLURM_TEMPLATE="""#!/bin/sh
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$TASKS
#SBATCH --cpus-per-task=$CPUSPERTASK
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --output=$JOBNAME.log
#SBATCH --error=$JOBNAME.err

$CMD
"""

def writeJobScript(cmd, jobName, nodes = 1, tasks = 1, cpusPerTask = 1, mem = 8000, time = "12:00:00"):
    """Write a job script to run a command via slurm.

    Args:
        cmd (str): The command to run.
        jobName (str): Name of the job (used for e.g. log files). Must not contain any spaces.
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks.
        cpusPerTask (int): Number of CPUs to use per task.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.

    Returns:
        The file name for the slurm batch file.

    Note:
        This doesn't submit jobs... see ``submitJob`` for a routine that does (it uses this routine).

    """

    script=SLURM_TEMPLATE
    script=script.replace("$CMD", cmd)
    script=script.replace("$JOBNAME", jobName)
    script=script.replace("$NODES", str(nodes))
    script=script.replace("$TASKS", str(tasks))
    script=script.replace("$CPUSPERTASK", str(cpusPerTask))
    script=script.replace("$MEM", str(mem))
    script=script.replace("$TIME", str(time))

    fileName=jobName+".slurm"
    with open(fileName, "w") as outFile:
        outFile.write(script)

    return fileName


def submitJob(cmd, jobName, dependentJobIDs = None, nodes = 1, tasks = 1, cpusPerTask = 1, mem = 8000, time = "12:00:00",
              cmdIsBatchScript = False):
    """Submit a command to run on a node via Slurm.

    Args:
        cmd (str): The command to run, OR the path to the batch script to submit (if cmdIsBatchScript = True).
        jobName (str): Name of the job (used for e.g. log files). Must not contain any spaces.
        dependentJobIDs (list): If this job depends on a previous job completing sucessfully, give the
            job ID numbers here, as a list.
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run.
        cpusPerTask (int): Number of CPUs to use per task.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.
        cmdIsBatchScript (bool): If False, a batch script is created for the command. If True, the cmd must be the path
            to a pre-existing slurm batch script (in which case most parameters e.g. nodes, tasks etc. are redundant).

    Returns:
        The ID number of the submitted Slurm job.

    Note:
        The job is submitted from the current working directory. So you *might*
        want to use absolute paths here.

    """

    #script=SLURM_TEMPLATE
    #script=script.replace("$CMD", cmd)
    #script=script.replace("$JOBNAME", jobName)
    #script=script.replace("$NODES", str(nodes))
    #script=script.replace("$TASKS", str(tasks))
    #script=script.replace("$MEM", str(mem))
    #script=script.replace("$TIME", str(time))

    #fileName=jobName+".slurm"
    #with open(fileName, "w") as outFile:
        #outFile.write(script)

    if cmdIsBatchScript == True:
        fileName=cmd
    else:
        fileName=writeJobScript(cmd, jobName, nodes = nodes, tasks = tasks, mem = mem, time = time)

    args=['sbatch']
    if dependentJobIDs is not None:
        if type(dependentJobIDs) != list:
            raise Exception("dependentJobIDs must be given as a list")
        dependStr="afterok"
        for dependentJobID in dependentJobIDs:
            dependStr=dependStr+":"+str(dependentJobID)
        args=args+['--dependency=%s' % (dependStr)]
    args=args+[fileName]

    process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise Exception("Non-zero return code when submitting job %s" % (jobName))
    assert(process.stdout[:19] == "Submitted batch job")
    jobID=int(process.stdout.split("Submitted batch job")[-1])

    return jobID

