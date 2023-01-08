#!/bin/bash
# PBS main script

#PBS -S /bin/bash
#PBS -l select=2:ncpus=20:mem=4gb
#PBS -l place=free
#PBS -l walltime=01:00:00
#PBS -q charon_2h
#PBS -N endorse_01
#PBS -j oe

set -x

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
WORKDIR="${SCRIPTPATH}/../sandbox/mlmc_run"

which python3



# singularity SIF image path (preferably create in advance)
# SING_FLOW=docker://flow123d:endorse_ci:91e460
SING_IMG="${SCRIPTPATH}/../endorse_ci.sif"
SING_IMG=${1:-${SING_IMG}}


# possibly set container mpiexec path
# IMG_MPIEXEC="/usr/local/mpich_3.4.2/bin/mpiexec"

# program and its arguments
PROG="endorse_mlmc run -c -np=4 --dim=2 'edz noedz' '2 5 10'"

mkdir -p ${WORKDIR}
cp ${SCRIPTPATH}/../test_data/*.yaml ${WORKDIR}
cd "${WORKDIR}"

# main call
# TODO: simplify swrap to work without installation
export PYTHONPATH=${SCRIPTPATH}/../../submodules/swrap/src/swrap
ls -l ${PYTHONPATH}

singularity exec $SING_IMG python3 -m pip install --upgrade --user -e ${SCRIPTPATH}/../..
#singularity exec docker://flow123d/geomop-gnu:2.0.0 venv/bin/python3 -m mlmc.tool.pbs_job /auto/liberec3-tul/home/jan_brezina/workspace/Endorse/tests/sandbox/mlmc_run/edz-002/output 0000 >/auto/liberec3-tul/home/jan_brezina/workspace/Endorse/tests/sandbox/mlmc_run/edz-002/output/jobs/0000_STDOUT 2>&1
python3 -m sexec -i $SING_IMG $PROG
# python3 $SING_SCRIPT -i $SING_FLOW -m $IMG_MPIEXEC -- $PROG
