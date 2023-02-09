#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Measurement/FullQCDCloverTraceless/MyFFApp
cd ${MyDIR}
Nt=$1
Ns=$2
a=$3  #BetaIndex=$3
b=$4 
ncpus=$5   #Total
quename=$6
model=$7
Select=$8
Mpiprocs=$9
Mem=${10}
for (( Trial=$a; Trial<$b; Trial++ )) #69
do
    BetaIndex=$Trial
    echo Submitting Job \# $Trial "," $Nt  $Ns ${BetaIndex} $ncpus  $quename  $model "," $Select $Mpiprocs  $Mem
    echo $model 
    #LogFile=/wsu/home/fy/fy41/fy4125/LogLattice/CommodoreFullQCD_Nt${Nt}Output${BetaIndex}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/LogLattice/CommodoreFullQCDClover_Nt${Nt}Err${BetaIndex}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job

  qsub -V  -q ${quename} -l select=${Select}:ncpus=${Mpiprocs}:mem=${Mem}gb:mpiprocs=${Mpiprocs}:cpu_model=${model}  -N Nt${Nt}FQC${BetaIndex} -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaIndex} ${ncpus}
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
