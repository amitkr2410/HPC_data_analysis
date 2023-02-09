#!/bin/bash

Nt=10
Ns=40
BetaArray=(4.0 4.5 5.0 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.3 6.5 7.0 7.5 8.0 8.5 9.0)
Configuration=5000

## 18 values of beta

for (( Trial=0; Trial<18; Trial++ ))
do 
    cd /wsu/home/fy/fy41/fy4125/Lattice/PolyakovLoopSU3/
    echo "Working on i=" $Trial
    LogFile=/dev/null
    ErrFile=/dev/null
    #LogFile=/wsu/home/fy/fy41/fy4125/Lattice/SU3/Log/Output${Trial}_${BetaArray[${Trial}]}.dat
    #ErrFile=/wsu/home/fy/fy41/fy4125/Lattice/SU3/Log/Err${Trial}_${BetaArray[${Trial}]}.dat
    Exec=/wsu/home/fy/fy41/fy4125/Lattice/PolyakovLoopSU3/AveragePolyakovLoopSU3
    N=$((202+((${Trial}/12))))
    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
    ## qsub -V -q eamxq -l mem=2gb   -N AmitQueue  -o $LogFile -e $ErrFile --  $Exec $Trial
    qsub -V -q wsuq  -l mem=2gb  -N APNt${Nt}_i$Trial  -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaArray[${Trial}]} ${Configuration}
    #qsub -V -q wsuq -l vnode=wsu${N} -l mem=2gb  -N ${Nt}_${Ns}_$Trial  -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaArray[${Trial}]} ${Configuration}

done



# wsuq mwsuq eamxq mwsuq zflq 
# mwsuq ezfhq zfh1-17 