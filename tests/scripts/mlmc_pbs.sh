#!/bin/bash
# PBS main script

#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l place=free
#PBS -l walltime=10:00:00
#PBS -q charon
#PBS -N endorse-main
#PBS -j oe

set -x

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
WORKDIR="${SCRIPTPATH}/../sandbox/mlmc_run"

which python3



# singularity SIF image path (preferably create in advance)
# SING_FLOW=docker://flow123d:endorse_ci:91e460
tag=7d9354
source_image="docker://flow123d/endorse_ci:${tag}"
image_file=${SCRIPTPATH}/../"endorse_ci_${tag}.sif"
if [ ! -f ${image_file} ]
then
    singularity build  ${image_file} ${source_image}  
fi


command=$1

# possibly set container mpiexec path
# IMG_MPIEXEC="/usr/local/mpich_3.4.2/bin/mpiexec"


# program and its arguments
endorse_cmd=${SCRIPTPATH}/../../venv/bin/endorse_mlmc
PARAM_SET="edz,noedz 2,5,10"

if [ "${command}" == "sample" ]
then
    SAMPLE_CMD="${endorse_cmd} run -c -nt=6 --dim=3 ${PARAM_SET}"

    mkdir -p ${WORKDIR}
    cd ${SCRIPTPATH}/../test_data
    cp *.yaml large_model_local.msh2 "${WORKDIR}"
    

    # main call
    # TODO: simplify swrap to work without installation
    export PYTHONPATH=${SCRIPTPATH}/../../submodules/swrap/src/swrap
    ls -l ${PYTHONPATH}

    #singularity exec $SING_IMG python3 -m pip install --upgrade --user -e ${SCRIPTPATH}/../..
    #singularity exec docker://flow123d/geomop-gnu:2.0.0 venv/bin/python3 -m mlmc.tool.pbs_job /auto/liberec3-tul/home/jan_brezina/workspace/Endorse/tests/sandbox/mlmc_run/edz-002/output 0000 >/auto/liberec3-tul/home/jan_brezina/workspace/Endorse/tests/sandbox/mlmc_run/edz-002/output/jobs/0000_STDOUT 2>&1
    cd "${WORKDIR}"
    python3 -m sexec -i $image_file $SAMPLE_CMD
    # python3 $SING_SCRIPT -i $SING_FLOW -m $IMG_MPIEXEC -- $PROG
elif [ "${command}" == "plot" ]
then

    PLOT_CMD="${endorse_cmd} plot cases ${PARAM_SET}"
    cd "${WORKDIR}"
    singularity exec $image_file $PLOT_CMD
fi
