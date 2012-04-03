#!/bin/sh
###
### PBS script for submitting a multithreaded job to a single (SMP) 
### node of the PSI batch cluster.
###
### You MUST edit this skeleton and insert the parameters and commands
### for your particular program.
###
### NOTE:  "#PBS" at the beginning of a line is NOT a comment.  To 
### comment out one of those lines, add a "### " to the beginning of
### the line.
###

### Set the job name.
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
### YOU MUST CHANGE THE JOB NAME
#PBS -N hw1-sequential-test

### Declare myprogram non-rerunable
#PBS -r n

### Set the queue to "batch", the only available queue. 
### This is only an issue for systems with multiple queues.
#PBS -q batch

### Specify the number of nodes and processors per node for your job.  
### This gives you one node (you only want one node, unless you are 
### using MPI) with 8 processors (a.k.a. "cores") on that node.  If 
### you get an error message here, it probably means that you're 
### running on Citris (which doesn't have 8-core nodes) instead of PSI.
#PBS -l nodes=1:ppn=8

### You can tell PBS how much memory you expect your job will use per
### node.  You can give it in gigabytes (for example, 1g) or megabytes
### (for example, 1024m).  PBS may ignore this value, but it's still a
### good idea to pick a reasonable upper bound.  Make sure that it's
### not more memory than the node actually has!
#PBS -l mem=1g

### Real-world time limit (here it's ten minutes).
#PBS -l walltime=00:10:00

### Change to the working directory, from which processes are
### launched.  By default Torque launches processes from your home
### directory.  Jobs should only be run from /home, /project, or
### /work, since Torque returns results via NFS.  /work is much
### preferable as it avoids overloading EECS home directories.
###
### Note: do NOT comment out this line!
cd $PBS_O_WORKDIR 

### Run some informational commands.  You can comment these out for
### production runs, but they are handy for debugging.
echo Working directory is $PBS_O_WORKDIR
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`

### Define number of processors in this batch script, which can be
### handy for invoking mpirun.  You can comment out the "echo" line.
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

### Set up the process environment.  YOU NEED THIS.  It only works
### if you run hw1.x out of the same directory in which you built it.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PBS_O_WORKDIR/dcmt0.4/lib

### Launch the process.
### YOU MUST CHANGE THE SECOND COMMAND-LINE ARGUMENT TO BE THE NUMBER
### OF THREADS.
$PBS_O_WORKDIR/hw1.x params.txt 1

