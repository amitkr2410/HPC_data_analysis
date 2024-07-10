#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/RUN/Analysis
cd ${MyDIR}
for((Trial=0; Trial<1; Trial++))
do    
i=0       ##i=0-5 pp,   0-1 pbpb/pp
M=2       #0=pp, 1-2=pp-pbpb
K=$1        #eloss type pp, Type5_alphas_0p3_Q0_2GeV 
C=$2         #Centrality=("10-20" "20-30"   "30-40" "40-50" "50-60"  )
HadronDirectory=$3
SoftV2Flag=$4
quename=$5
#PPorPbPb=(pp pp PbPb)     # pp or pp or PbPb
#BGS=(-1  0   1)           # -1 or 0 or  1
#BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
#ElossType=("0p20_Q0_1p0"  "0p20_Q0_1p4"  "0p20_Q0_2p0" "0p25_Q0_1p0"  "0p25_Q0_1p4" )
    
    echo Submitting Job \# $Trial "," $i  $M   $K $C ${SoftV2Flag} ${quename}
    #LogFile=${MyDIR}/Log/ScreenOutput${Trial}.out
    ErrFile=${MyDIR}/Log/CommodorHadronV2_${i}_${M}_${K}_${C}_${SoftV2Flag}_${quename}.err
    LogFile=/dev/null
    #  runPbPb.sh  PPrunSingleHadron.sh  PbPbrunSingleHadron.sh
    Exec=${MyDIR}/PbPbrunHadronV2.sh
    
    sbatch -q ${quename}  -N 1 -n 1 --mem=8G -t 2-0   --job-name ${i}${M}HV2${C}_${SoftV2Flag}  -o ${LogFile} -e ${ErrFile} --  ${Exec} $i $M  $K $C ${HadronDirectory} ${SoftV2Flag}
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
