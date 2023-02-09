#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Measurement/PureGaugeCloverTraceless/AmitPureGauge
cd ${MyDIR}
NT=(4  6  8  16 24 32)
NS=(16 24 32 16 24 32)

i=1
BetaIndex=$1
queue=$2
Select=$3
Nodes=$4
CPUperNode=$5
MEM=$6
    echo Submitting Job \# Nt${NT[$i]} BetaIndex= ${BetaIndex}
    #LogFile=/wsu/home/fy/fy41/fy4125/Log/ScreenNewNt4Output${Trial}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/Log/ScreenPureGaugeNt${NT[$i]}New${BetaIndex}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job_generate

   ## wsu161-wsu184 (28cpus each node, 128GB)
   ##   qsub -V -q mwsuq  -l cpu_type=Intel -l cpu_model=E5-2697v3 -l mem=16gb -l ncpus=8 -l mpiprocs=8 -N Nt4i$Trial -o $LogFile -e $ErrFile --        $Exec $Trial
    sbatch -q ${queue}  -N ${Nodes} -n ${CPUperNode} --mem=${MEM}G -t 7-0   --job-name  ${NT[${i}]}   -o ${LogFile} -e ${ErrFile} --  ${Exec} ${NT[$i]} ${NS[$i]} ${BetaIndex} ${Select}

   #qsub -V   -l select=1:ncpus=16:mem=16gb:mpiprocs=16:cpu_type=Intel  -N Nt16i$Trial -o $LogFile -e $ErrFile --  $Exec $Trial





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
