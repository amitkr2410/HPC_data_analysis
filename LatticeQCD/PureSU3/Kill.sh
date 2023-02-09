#!/bin/bash

First=455732
 Last=455743


LastPlusOne=$(( ${Last} - ${First} +1))
for(( i=0; i< ${LastPlusOne}; i++ ))
do
JOBID=$((${First}+${i}))
echo "JOBID is " ${JOBID}
qdel ${JOBID}
done