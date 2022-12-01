#!/bin/bash

set -x

#py_script=`pwd`/$1
pbs_script="run_endorse.pbs"#`pwd`/$1.pbs
#script_path=${py_script%/*}

#work_dir=$2

cat >$pbs_script <<EOF
#!/bin/bash
#PBS -S /bin/bash
#PBS -l select=1:ncpus=16:cgroups=cpuacct:mem=16Gb
#PBS -l walltime=2:00:00
#PBS -q charon_2h
#PBS -N endorse
#PBS -j oe

#export TMPDIR=$SCRATCHDIR

export SINGULARITY_TMPDIR=$SCRATCHDIR


cd /storage/liberec3-tul/home/martin_spetlik/Endorse_full_transport

singularity exec docker://flow123d/endorse:latest ./setup_sing.sh

singularity exec docker://flow123d/endorse:latest python3 /storage/liberec3-tul/home/martin_spetlik/Endorse_MS_full_transport/tests/test_transport.py


EOF

qsub $pbs_script
