#!/bin/sh

UPSTREAMURL="https://github.com/SKSNSim/SKSNSim/releases/download/v1.2.0-data/supernova_data_260508.tar.gz"

if [ -z "${SKSNSIMDATADIR}" ]; then
  echo "Environmental variable \"SKSNSIMDATADIR\" is not defined."
  echo "Please configure the variable."
  exit 1
fi

extractdir=$(dirname "${SKSNSIMDATADIR}")

if [ ! -d ${extractdir} ]; then
  echo "The data directory does not exist."
  echo "Making directory "${extractdir}
  mkdir -p ${extractdir}
fi

target=${extractdir}/supernova_data_260508.tar.gz

which wget || (echo "no wget on your system. Please download manually from \"${UPSTREAMURL}\"." && exit 1)
wget -O $target $UPSTREAMURL

origdir=$(pwd)
cd ${extractdir} && tar -xvzf $target
cd $origdir
