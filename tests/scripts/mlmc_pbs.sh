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

# collect arguments:

# singularity_exec_mpi script path
SING_SCRIPT="${SCRIPTPATH}/../../submodules/swrap/src/swrap/sexec.py"

# singularity SIF image path (preferably create in advance)
# SING_FLOW=docker://flow123d:endorse_ci:91e460

SING_IMG="${SCRIPTPATH}/../endorse_ci.sif"


# possibly set container mpiexec path
# IMG_MPIEXEC="/usr/local/mpich_3.4.2/bin/mpiexec"

# program and its arguments
PROG="python3 endorse_mlmc run -c -np=4  'edz noedz' '2 5 10'"

cd "${WORKDIR}"

# main call
python3 $SING_SCRIPT -i $SING_IMG -- $PROG
# python3 $SING_SCRIPT -i $SING_FLOW -m $IMG_MPIEXEC -- $PROG
