#! /bin/bash
#$ -A @ALLOCATION@
#$ -q @QUEUE@
#$ -r n
#$ -l h_rt=@WALLTIME@
#$ -pe @(@PPN_USED@/@NUM_THREADS@)@way @PROCS_REQUESTED@
#$ @('@CHAINED_JOB_ID@' ne '' ? '-hold_jid @CHAINED_JOB_ID@' : '')@
#$ -V
#$ -N @SHORT_SIMULATION_NAME@
#$ -M @EMAIL@
#$ -m abe
#$ -o @RUNDIR@/@SIMULATION_NAME@.out
#$ -e @RUNDIR@/@SIMULATION_NAME@.err

echo "Preparing:"
set -x                          # Output commands
set -e                          # Abort on errors

if test -n '@CHAINED_JOB_ID@'; then
    # This script was pre-submitted.  Clean up the simulation, then
    # submit it again without queuing a new job.
    echo "This is a chained script"
    echo "This is restart id \"@RESTART_ID@\" on host \"@HOSTNAME@\""
    echo "This is job \"$JOB_ID\", depending on job \"@CHAINED_JOB_ID@\""
    echo "Preparing simulation for restart:"
    cd @SOURCEDIR@
    # Specify --recover to ensure that the new restart either recovers
    # or fails
    ./simfactory/sim cleanup-submit @SIMULATION_NAME@ --hostname=@HOSTNAME@ --postsubmit=@RESTART_ID@ --recover
fi

cd @RUNDIR@-active
mkdir -p @SCRATCHDIR@
ln -s @SCRATCHDIR@ scratch
export GMON_OUT_PREFIX=gmon.out

# NOTE: The call to ibrun below depends on environment settings; make
# sure these module choices correspond to the MPI selection in the
# options file
module unload openmpi
module unload mvapich
module unload mvapich-devel
module unload mvapich-ud
module unload mvapich2
module load mvapich/0.9.9

echo "Checking:"
pwd
hostname
date

echo "PAPI:"
module load papi/3.6.0
export PAPIDIR=/opt/apps/papi/papi-3.6.0/bin
echo "papi_avail:"        && $PAPIDIR/papi_avail        || true
echo "papi_clockres:"     && $PAPIDIR/papi_clockres     || true
echo "papi_command_line:" && $PAPIDIR/papi_command_line || true
echo "papi_cost:"         && $PAPIDIR/papi_cost         || true
echo "papi_decode:"       && $PAPIDIR/papi_decode       || true
echo "papi_mem_info:"     && $PAPIDIR/papi_mem_info     || true
echo "papi_native_avail:" && $PAPIDIR/papi_native_avail || true

echo "Environment:"
export OMP_NUM_THREADS=@NUM_THREADS@
env > SIMFACTORY/ENVIRONMENT

cat $PE_HOSTFILE > SIMFACTORY/NODES || true

echo "Starting:"
export CACTUS_STARTTIME=$(date +%s)
time /share/sge6.2/default/pe_scripts/ibrun /share/sge6.2/default/pe_scripts/tacc_affinity ./@EXECUTABLE@ -L 3 @PARFILE@
# This should work as well
#time mpiexec -machinefile NODES -np @NUM_PROCS@ tacc_affinity ./@EXECUTABLE@ -L 3 @PARFILE@

echo "Stopping:"
date

echo "Done."
