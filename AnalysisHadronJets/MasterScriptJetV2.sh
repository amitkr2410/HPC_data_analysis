#!/bin/bash 
pTMIN=(0.0    0.2  0.5  0.8 1.0 )
pTMINS=(0p0   0p2  0p5  0p8 1p0 )
MyDIR=/wsu/home/fy/fy41/fy4125/RUN/Analysis
cd ${MyDIR}
for((Trial=0; Trial<1; Trial++))
do    
i=${1}           #i=0-5 pp,   0-1 pbpb/pp                                                       
M=$2           #0=pp, 1-2=pp-pbpb                                                                                           
K=$3           #eloss type:  pp, Type5_alphas_0p3_Q0_2GeV 
C=${4}          #(0="0-10"  1="0-5" 2="5-10" 3="10-20" 4="20-30"  5="30-40" 6="40-50" 7="50-60" 8="20-40" )              
HadronDirectoryFullPath=$5
EtaJetMIN=$6    #0.0                                                                                                        
EtaJetMINS=$7  #0p0                                                                                                         
EtaJetMAX=$8   #2.8                                                                                                         
EtaJetMAXS=$9  #2p8                                                                                                         
#pTMIN=$9 #pTTrackMin =0.2 0.5 1.0
#pTMINS=$10 #0p0 0p2 0p5 0p8 1p0                                                                                            
quename=${10}
JetRADIUS=${11}
JetRADIUSS=${12}
EXP=${13}          #Atlas Alice CMS
EtaCutFlAG=${14}   ## for |y| 0, for eta 1
PPorPbPb=${15} 
ECM=${16}       #200 2760 5020
SoftV2Flag=${17}
#PPorPbPb=(pp ppPbPb PbPb)     # pp or ppPbPb or pppPb or PbPb
#BGS=(-1  0   1)           # -1 or 0 or  1
#BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
#ElossType=("0p20_Q0_1p0"  "0p20_Q0_1p4"  "0p20_Q0_2p0" "0p25_Q0_1p0"  "0p25_Q0_1p4" )    
    echo Submitting Job \# ${Trial} "," ${EtaJetMINS}_${EtaJetMAXS} $i $M  $K $C ${SoftV2Flag} ${quename}
    #LogFile=${MyDIR}/Log/ScreenOutput${Trial}.out
    ErrFile=${MyDIR}/Log/CommodorJet_${EtaJetMINS}_${EtaJetMAXS}_${i}_${M}_${K}_${C}_${SoftV2Flag}.err
    LogFile=/dev/null
    #  runPbPb.sh  PPrunSingleHadron.sh  PbPbrunSingleHadron.sh
    Exec=${MyDIR}/runPbPbJetV2.sh
    
    sbatch -q ${quename}  -N 1 -n 1 --mem=8G -t 4-0   --job-name ${i}${M}J${K}_${C}_${EtaJetMINS}_${EtaJetMAXS}_${Trial}   -o ${LogFile} -e ${ErrFile} --  ${Exec} $i  $M $K $C ${EtaJetMIN} ${EtaJetMINS} ${EtaJetMAX} ${EtaJetMAXS} ${pTMIN[${Trial}]} ${pTMINS[${Trial}]} ${HadronDirectoryFullPath}  ${JetRADIUS}  ${JetRADIUSS}  ${EXP}   ${EtaCutFlAG}  ${PPorPbPb} ${ECM} ${SoftV2Flag}
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
