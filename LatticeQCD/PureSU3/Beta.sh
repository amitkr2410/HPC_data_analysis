#!/bin/bash

Nt=4
Ns=16
BetaArray=(5.81 5.83 6.27 6.34)
Configuration=5000

## 33 values of beta

for (( Trial=3; Trial<4; Trial++ ))
do 
    cd /wsu/home/fy/fy41/fy4125/Lattice/SU3/
    echo "Working on i=" $Trial
    LogFile=/dev/null
    #ErrFile=/dev/null
    #LogFile=/wsu/home/fy/fy41/fy4125/Lattice/SU3/Log/Output${Trial}_${BetaArray[${Trial}]}.dat
    ErrFile=/wsu/home/fy/fy41/fy4125/Lattice/SU3/Log/Err${Trial}_${BetaArray[${Trial}]}.dat
    Exec=/wsu/home/fy/fy41/fy4125/Lattice/SU3/AverageCorrelatorSU3
    #Exec=/wsu/home/fy/fy41/fy4125/Lattice/SU3/AveragePlaquetteSU3
    #N=$((202+((${Trial}/12))))
    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
    #qsub -V -q eamxq -l mem=2gb   -N ${Nt}_${Ns}_${Trial}  -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaArray[${Trial}]} ${Configuration}
qsub -V -q eamxq  -l mem=2gb  -N PG4_$Trial  -o $LogFile -e $ErrFile -- $Exec $Nt $Ns ${Configuration} ${BetaArray[${Trial}]}
    #qsub -V -q wsuq -l vnode=wsu${N} -l mem=2gb  -N ${Nt}_${Ns}_$Trial  -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaArray[${Trial}]} ${Configuration}

done



# wsuq mwsuq eamxq mwsuq zflq 
# mwsuq ezfhq zfh1-17 
# qselect -u fy4125 | xargs qdel
