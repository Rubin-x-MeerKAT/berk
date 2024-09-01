"""

Tools for running jobs via a workload manager (Slurm or PBS).

"""

import subprocess, os, sys

SLURM_TEMPLATE="""#!/bin/sh
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$TASKS
#SBATCH --cpus-per-task=$CPUSPERTASK
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --output=$JOBNAME.log
#SBATCH --error=$JOBNAME.err

$CMD
sleep 10
"""

PBS_TEMPLATE="""#!/bin/sh
#PBS -l nodes=$NODES:ppn=$TASKS,mem=$MEM
#PBS -P $PBS_PROJECT
#PBS -q $PBS_QUEUE
#PBS -l walltime=$TIME
#PBS -o $JOBNAME.log
#PBS -e $JOBNAME.err
#PBS -m abe
#PBS -M $PBS_EMAIL
ulimit -s unlimited
cd $CWD

$CMD
sleep 10
"""

def writeJobScript(cmd, jobName, nodes = 1, tasks = 1, mem = 8000, time = "12:00:00",
                   workloadManager = 'slurm'):
    """Write a job script to run a command via slurm.

    Args:
        cmd (str): The command to run.
        jobName (str): Name of the job (used for e.g. log files).
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run per node.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.
        workloadManager (str): Either 'slurm' or 'pbs'. Note that for PBS, several
            environment variables also need to be defined (`PBS_NODETYPE`,
            `PBS_PROJECT`, `PBS_QUEUE`, `PBS_EMAIL`).

    Returns:
        The file name for the batch file.

    Note:
        This doesn't submit jobs... see ``submitJob`` for a routine that does (it uses this routine).

    """

    if workloadManager == 'slurm':
        script=SLURM_TEMPLATE
    elif workloadManager == 'pbs':
        script=PBS_TEMPLATE
    else:
        raise Exception("workloadManager should be either 'slurm' or 'pbs'")
    script=script.replace("$CMD", cmd)
    script=script.replace("$JOBNAME", jobName)
    script=script.replace("$NODES", str(nodes))
    script=script.replace("$TASKS", str(tasks))
    if workloadManager == 'slurm':
        script=script.replace("$MEM", str(mem)) # MB
    elif workloadManager == 'pbs':
        script=script.replace("$MEM", "%.dgb" % (mem/1000))
    script=script.replace("$TIME", str(time))
    if workloadManager == 'pbs':
        envVars=['PBS_NODETYPE', 'PBS_PROJECT', 'PBS_QUEUE', 'PBS_EMAIL']
        for e in envVars:
            script=script.replace("$%s" % (e), os.environ[e])
        script=script.replace("$CWD", os.path.abspath(os.path.curdir))

    fileName=jobName+"."+workloadManager
    with open(fileName, "w") as outFile:
        outFile.write(script)

    return fileName


def submitJob(cmd, jobName, dependentJobIDs = None, nodes = 1, tasks = 1, mem = 8000,
              time = "12:00:00", workloadManager = 'slurm', cmdIsBatchScript = False):
    """Submit a command to run on a node via Slurm or PBS.

    Args:
        cmd (str): The command to run.
        jobName (str): Name of the job (used for e.g. log files).
        dependentJobIDs (list): If this job depends on a previous job completing sucessfully, give the
            job ID numbers here, as a list.
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run per node.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.
        workloadManager (str): Either 'slurm' or 'pbs'. Note that for PBS, several
            environment variables also need to be defined (`PBS_NODETYPE`,
            `PBS_PROJECT`, `PBS_QUEUE`, `PBS_EMAIL`).

    Returns:
        The ID number of the submitted Slurm job.

    Note:
        The job is submitted from the current working directory. So you *might*
        want to use absolute paths here.

    """

    if cmdIsBatchScript == True:
        fileName=cmd
    else:
        fileName=writeJobScript(cmd, jobName, nodes = nodes, tasks = tasks, mem = mem, time = time,
                                workloadManager = workloadManager)

    if workloadManager == 'slurm':
        args=['sbatch']
        if dependentJobIDs is not None:
            if type(dependentJobIDs) != list:
                raise Exception("dependentJobIDs must be given as a list")
            dependStr="afterok"
            for dependentJobID in dependentJobIDs:
                dependStr=dependStr+":"+str(dependentJobID)
            args=args+['--dependency=%s' % (dependStr)]
        args=args+[fileName]
    elif workloadManager == 'pbs':
        args=['qsub']
        if dependentJobIDs is not None:
            if type(dependentJobIDs) != list:
                raise Exception("dependentJobIDs must be given as a list")
            dependStr="afterok"
            for dependentJobID in dependentJobIDs:
                dependStr=dependStr+":"+str(dependentJobID)
            args=args+['-W depend=%s' % (dependStr)]
        args=args+[fileName]
    else:
        raise Exception("workloadManager should be either 'slurm' or 'pbs'")

    process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise Exception("Non-zero return code when submitting job %s" % (jobName))
    if workloadManager == 'slurm':
        assert(process.stdout[:19] == "Submitted batch job")
        jobID=int(process.stdout.split("Submitted batch job")[-1])
    elif workloadManager == 'pbs':
        jobID=int(process.stdout.split(".")[0])

    return jobID

