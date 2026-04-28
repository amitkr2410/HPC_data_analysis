#!/bin/bash

JETSCAPEDIR=$1
OutputXML=$2
TestOutName=$3
HadronFile=$4
PartonFile=$5
LogFile=$6

echo $1 $2  $3 $4  $5 $6 

cd ${JETSCAPEDIR}/build

echo ${PWD}
InitialTime=${SECONDS}

./runJetscape  ${OutputXML}  >> ${LogFile}

./FinalStateHadrons  ${TestOutName}.dat  ${HadronFile}   >>  ${LogFile}

./FinalStatePartons  ${TestOutName}.dat   ${PartonFile}  >>  ${LogFile}


FinalTime=${SECONDS}
Duration=$((FinalTime - InitialTime))   #seconds

echo ${Duration} >>  ${LogFile}
HH=$((Duration/3600))
MM=$((Duration/60 - HH*60))
SS=$((Duration - HH*3600 - MM*60))
echo "HH::MM::SS=${HH}::${MM}::${SS}"  >> ${LogFile}
