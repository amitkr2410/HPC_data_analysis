#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Measurement/FullQCDCloverTraceless/MyFFApp
cd ${MyDIR}
Nt=$1
Ns=$2
a=$3
b=$4
#BetaIndex=$3
ncpus=$5
mem=$6
QueName=$7
CPUModel=$8

for (( Trial=a; Trial<b; Trial++ )) #69
do
    BetaIndex=$Trial
    echo Submitting Job \# $Trial "," $Nt  $Ns ${BetaIndex} $ncpus

    #LogFile=/wsu/home/fy/fy41/fy4125/LogLattice/ScreenNewNt4Output${BetaIndex}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/LogLattice/CommodoreFullQCD_Nt${Nt}${BetaIndex}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job
# ./simple_job  $Nt  $Ns   $Trial  4
    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
   ##qsub -V -q eamxq -l mem=32gb -l ncpus=8 -l mpiprocs=8  -N Nt4i$Trial  -o $LogFile -e $ErrFile --  $Exec $Trial
    ## wsu201-209 (40cpus each node, 1536 GBRAM)  cpu_model=E5-4627v4 -l cpu_type=Intel
   #qsub -V -q wsuq  -l select=4:ncpus=8:mem=16gb:mpiprocs=8:cpu_type=Intel -N Nt16i$Trial -o $LogFile -e $ErrFile --  $Exec $Trial 
   qsub -V  -q ${QueName} -l select=${ncpus}:ncpus=1:mem=${mem}gb:mpiprocs=1:cpu_model=${CPUModel}  -N Nt${Nt}FQ${BetaIndex} -o $LogFile -e $ErrFile --  $Exec $Nt $Ns ${BetaIndex} ${ncpus}
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
