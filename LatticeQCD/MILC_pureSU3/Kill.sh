#!/bin/bash

First=2907135
 Last=2907141


LastPlusOne=$(( ${Last} - ${First} +1))
for(( i=0; i< ${LastPlusOne}; i++ ))
do
JOBID=$((${First}+${i}))
echo "JOBID is " ${JOBID}
qdel ${JOBID}
done
