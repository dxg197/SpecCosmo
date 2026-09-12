#! /bin/bash
#BSUB -P @ALLOCATION@
#BSUB -q @QUEUE@
#BSUB -W @WALLTIME_HH@:@WALLTIME_MM@
#BSUB -n @PROCS@
#BSUB -R 'span[ptile=@PPN@]'
#BSUB @('@CHAINED_JOB_ID@' ne '' ? '-w ended(@CHAINED_JOB_ID@)' : '')@
#BSUB -J @SIMULATION_NAME@
#BSUB -N
#BSUB -o @RUNDIR@/@SIMULATION_NAME@.out
#BSUB -e @RUNDIR@/@SIMULATION_NAME@.err

echo "Preparing:"
set -x                          # Output commands
set -e                          # Abort on errors

if test -n '@CHAINED_JOB_ID@'; then
    # This script was pre-submitted.  Clean up the simulation, then
    # submit it again without queuing a new job.
    echo "This is a chained script"
    echo "This is restart id \"@RESTART_ID@\" on host \"@HOSTNAME@\""
    echo "This is job \"$PBS_JOBID\", depending on job \"@CHAINED_JOB_ID@\""
    echo "Preparing simulation for restart:"
    cd @SOURCEDIR@
    # Specify --recover to ensure that the new restart either recovers
    # or fails
    ./simfactory/sim cleanup-submit @SIMULATION_NAME@ --hostname=@HOSTNAME@ --postsubmit=@RESTART_ID@
fi

cd @RUNDIR@-active
mkdir -p @SCRATCHDIR@
ln -s @SCRATCHDIR@ scratch
export GMON_OUT_PREFIX=gmon.out

echo "Checking:"
pwd
hostname
date

echo "LSB nodes:"
cat ${LSB_NODEFILE}
cat ${LSB_NODEFILE} > SIMFACTORY/NODES

echo "Environment:"
env > SIMFACTORY/ENVIRONMENT

echo "Starting:"
export CACTUS_STARTTIME=$(date +%s)
ibrun ./@EXECUTABLE@ -L 3 @PARFILE@

echo "Stopping:"
date

echo "Done."
