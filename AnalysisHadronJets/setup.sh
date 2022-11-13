#!/bin/sh

if [ -z $LD_LIBRARY_PATH ]; then
export LD_LIBRARY_PATH
fi

if [ -z $DYLD_LIBRARY_PATH ]; then 
export DYLD_LIBRARY_PATH
fi

export BASEDIR=${HOME}
export PYTHIAINSTALLDIR=/wsu/home/fy/fy41/fy4125/Software

export JetScape=${PWD}/lib
export LD_LIBRARY_PATH=${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
export DYLD_LIBRARY_PATH=${JetScape}:${DYLD_LIBRARY_PATH}

# PYTHIA8 directory
export PYTHIA8_DIR=${PYTHIAINSTALLDIR}/pythia8230
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8230
#export PYTHIA8_INCLUDE_DIR=`${PYTHIA8DIR}/bin/pythia8-config --includedir`/Pythia8
#export PYTHIA8_LIBRARIES=`${PYTHIA8DIR}/bin/pythia8-config --libdir`
export LD_LIBRARY_PATH=${PYTHIA8_DIR}/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${PYTHIA8_DIR}/lib:${DYLD_LIBRARY_PATH}

#ROOT setup
##export ROOTSYS=/home/amit/root
##export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

#BOOST
##export BOOST_ROOT=/home/amit/JET/boost_1_64_0

if [ -z ${TERM} -o -z ${SHELL} ]; then exit 0
fi

echo ''
echo 'Setup JetScape Library'
echo '======================'
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "JetScape: " ${JetScape}
echo "Pythia8: "${PYTHIA8_DIR}/lib    
echo "<I>---------------Info--------------------<I>"
echo ""
