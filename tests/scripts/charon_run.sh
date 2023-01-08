#!/bin/bash

set -x

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
WORKDIR="${SCRIPTPATH}/../sandbox/mlmc_run"

tag=7d9354
source_image="docker://flow123d/endorse_ci:${tag}"
image_file=${SCRIPTPATH}/../"endorse_ci_${tag}.sif"
if [ ! -f ${image_file} ]
then
    singularity build  ${image_file} ${source_image}  
fi

QSUB=""
#QSUB=qsub --
${QSUB} ${SCRIPTPATH}/mlmc_pbs.sh ${image_file}

