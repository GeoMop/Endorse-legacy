#!/bin/bash

export SINGULARITY_TMPDIR=$SCRATCHDIR

#singularity build docker://flow123d/endorse:latest

#singularity build endorse.sif docker://flow123d/endorse:latest


cd /storage/liberec3-tul/home/martin_spetlik/Endorse_full_transport

singularity exec docker://flow123d/endorse:latest ./setup.sh

#singularity exec docker://flow123d/endorse:latest . venv/bin/activate

cd /storage/liberec3-tul/home/martin_spetlik/Endorse_MS_full_transport/tests/mlmc

singularity exec docker://flow123d/endorse:latest python3 fullscale_transport.py run ../ --clean


#./storage/liberec3-tul/home/martin_spetlik/Endorse_full_transport/setup.sh