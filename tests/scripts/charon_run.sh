#!/bin/bash

set -x

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
WORKDIR="${SCRIPTPATH}/../sandbox/mlmc_run"

qsub mlmc_pbs.sh

