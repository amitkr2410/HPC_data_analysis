#!/bin/bash

First=45302
 Last=45334


LastPlusOne=$(( ${Last} - ${First} +1))
for(( i=0; i< ${LastPlusOne}; i++ ))
do
JOBID=$((${First}+${i}))
echo "JOBID is " ${JOBID}
qdel ${JOBID}
done