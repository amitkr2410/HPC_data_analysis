#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/RUN/FastJet
cd ${MyDIR}
for((Trial=0; Trial<1; Trial++))
do    
i=0           #i=0-5 pp,   0-1 pbpb/pp
M=1           #0=pp, 1-2=pp-pbpb
K=0           #eloss type K=0a, 1b, 2c, 3d, 4e
#PPorPbPb=(pp pp PbPb)     # pp or pp or PbPb
#BGS=(-1  0   1)           # -1 or 0 or  1
#BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
#ElossType=("0p20_Q0_1p0"  "0p20_Q0_1p4"  "0p20_Q0_2p0" "0p25_Q0_1p0"  "0p25_Q0_1p4" )
    
    echo Submitting Job \# $Trial
    #LogFile=${MyDIR}/Log/ScreenOutput${Trial}.out
    ErrFile=${MyDIR}/Log/Commodor_$i_$M_$k.err
    LogFile=/dev/null
    #  runPbPb.sh  PPrunSingleHadron.sh  PbPbrunSingleHadron.sh
    Exec=${MyDIR}/runPbPb.sh
    
    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
    ## qsub -V -q eamxq -l mem=2gb   -N AmitQueue  -o $LogFile -e $ErrFile --  $Exec $Trial
    qsub -V -q eamxq  -l mem=4gb  -l ncpus=1 -N $i${M}Jet$K   -o $LogFile -e $ErrFile --  $Exec $i $M $K 
done

##qsub -V -I -q eamxq -- /wsu/home/fy/fy41/fy4125/RUN/PP1/simple_job 0




##      Script is submitted to this Queue: 
## To delete all the jobs given by user say fy4125
## then type in terminal " qselect -u fy4125 | xargs qdel

## To delete a process with job ID say 458393
## then type in terminal "qdel 458393" 
##   qstat -Q    "to know jobs running on different queue"
##   wsuq mwsuq eamxq mwsuq zflq ezfhq
##PBS  -q eamxq                 
##PBS -q mwsuq
                                                                             
##PBS -q zfhq 
                                                                             
##      One core and 1GB of RAM selected:
                                                  

##PBS -l select=1:ncpus=1:mem=10GB
                                                         
##PBS -l select=1:ncpus=4:mem=60GB:cpu_speed=3.3 
                                          
##PBS -l select=1:ncpus=1:mem=50GB:cpu_speed=3.0
                                           
## vmem=25GB vnode=amx3 
                                                                   

##export DISPLAY=localhost:0.0 
                                                            
##ulimit -s unlimited 
