#!/bin/sh

UPSTREAMURL="https://github.com/SKSNSim/SKSNSim/releases/download/v1.2.0-data/supernova_data.tar.gz"

if [ -z "${SKSNSIMDATADIR}" ]; then
  echo "Environmental variable \"SKSNSIMDATADIR\" is not defined."
  echo "Please configure the variable."
  exit 1
fi

if [ ! -d ${SKSNSIMDATADIR} ]; then
  echo "The data directory does not exit."
  echo "Making directory "$SKSNSIMDATADIR
  mkdir -p $SKSNSIMDATADIR
fi

target=${SKSNSIMDATADIR}/supernova_data.tar.gz

which wget || (echo "no wget on your system. Please download manually from \"${UPSTREAMURL}\"." && exit 1)
wget -O $target $UPSTREAMURL

origdir=$(pwd)
cd $SKSNSIMDATADIR && tar -xvzf $target
cd $origdir

