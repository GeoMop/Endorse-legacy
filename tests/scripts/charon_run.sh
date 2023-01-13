#!/bin/bash

set -x

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

qsub ${SCRIPTPATH}/mlmc_pbs.sh
#QSUB=""
#QSUB="qsub -o ${SCRIPTPATH}/endorse-main.OUT --"
#QSUB="qsub --"
#echo ${QSUB}
#${QSUB} ${SCRIPTPATH}/mlmc_pbs.sh  sample



