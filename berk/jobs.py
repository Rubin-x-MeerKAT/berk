"""

Tools for running jobs via a workload manager (Slurm or PBS).

"""

import os
import subprocess

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

def write_job_script(cmd, job_name, nodes = 1, tasks = 1, mem = 8000, time = "12:00:00",
                   workload_manager = 'slurm'):
    """Write a job script to run a command via slurm.

    Args:
        cmd (str): The command to run.
        job_name (str): Name of the job (used for e.g. log files).
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run per node.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.
        workload_manager (str): Either 'slurm' or 'pbs'. Note that for PBS, several
            environment variables also need to be defined (`PBS_NODETYPE`,
            `PBS_PROJECT`, `PBS_QUEUE`, `PBS_EMAIL`).

    Returns:
        The file name for the batch file.

    Note:
        This doesn't submit jobs... see ``submit_job`` for a routine that does
        (it uses this routine).

    """

    if workload_manager == 'slurm':
        script=SLURM_TEMPLATE
    elif workload_manager == 'pbs':
        script=PBS_TEMPLATE
    else:
        raise Exception("workload_manager should be either 'slurm' or 'pbs'")
    script=script.replace("$CMD", cmd)
    script=script.replace("$JOBNAME", job_name)
    script=script.replace("$NODES", str(nodes))
    script=script.replace("$TASKS", str(tasks))
    if workload_manager == 'slurm':
        script=script.replace("$MEM", str(mem)) # MB
    elif workload_manager == 'pbs':
        script=script.replace("$MEM", "%.dgb" % (mem/1000))
    script=script.replace("$TIME", str(time))
    if workload_manager == 'pbs':
        env_vars=['PBS_NODETYPE', 'PBS_PROJECT', 'PBS_QUEUE', 'PBS_EMAIL']
        for e in env_vars:
            script=script.replace("$%s" % (e), os.environ[e])
        script=script.replace("$CWD", os.path.abspath(os.path.curdir))

    filename=job_name+"."+workload_manager
    with open(filename, "w", encoding = 'utf8') as out_file:
        out_file.write(script)

    return filename


def submit_job(cmd, job_name, dependent_job_ids = None, nodes = 1, tasks = 1, mem = 8000,
              time = "12:00:00", workload_manager = 'slurm', cmd_is_batch_script = False):
    """Submit a command to run on a node via Slurm or PBS.

    Args:
        cmd (str): The command to run.
        job_name (str): Name of the job (used for e.g. log files).
        dependent_job_ids (list): If this job depends on a previous job completing sucessfully,
            give the job ID numbers here, as a list.
        nodes (int): Number of nodes on which the job will run.
        tasks (int): Number of tasks to run per node.
        mem (int): Requested memory (KB) for the job.
        time (str): Wall time limit for the job.
        workload_manager (str): Either 'slurm' or 'pbs'. Note that for PBS, several
            environment variables also need to be defined (`PBS_NODETYPE`,
            `PBS_PROJECT`, `PBS_QUEUE`, `PBS_EMAIL`).

    Returns:
        The ID number of the submitted Slurm job.

    Note:
        The job is submitted from the current working directory. So you *might*
        want to use absolute paths here.

    """

    if cmd_is_batch_script is True:
        filename=cmd
    else:
        filename=write_job_script(cmd, job_name, nodes = nodes, tasks = tasks, mem = mem,
                                  time = time, workload_manager = workload_manager)

    if workload_manager == 'slurm':
        args=['sbatch']
        if dependent_job_ids is not None:
            if isinstance(dependent_job_ids, list) is False:
                raise Exception("dependent_job_ids must be given as a list")
            depend_str="afterok"
            for dependent_job_id in dependent_job_ids:
                depend_str=depend_str+":"+str(dependent_job_id)
            args=args+['--dependency=%s' % (depend_str)]
        args=args+[filename]
    elif workload_manager == 'pbs':
        args=['qsub']
        if dependent_job_ids is not None:
            if isinstance(dependent_job_ids, list) is False:
                raise Exception("dependent_job_ids must be given as a list")
            depend_str="afterok"
            for dependent_job_id in dependent_job_ids:
                depend_str=depend_str+":"+str(dependent_job_id)
            args=args+['-W depend=%s' % (depend_str)]
        args=args+[filename]
    else:
        raise Exception("workload_manager should be either 'slurm' or 'pbs'")

    process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise Exception("Non-zero return code when submitting job %s" % (job_name))
    if workload_manager == 'slurm':
        assert process.stdout[:19] == "Submitted batch job"
        job_id=int(process.stdout.split("Submitted batch job")[-1])
    elif workload_manager == 'pbs':
        job_id=int(process.stdout.split(".")[0])
    else:
        raise Exception("workload_manager should be either 'slurm' or 'pbs'")

    return job_id
