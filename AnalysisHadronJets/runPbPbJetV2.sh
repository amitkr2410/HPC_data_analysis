#!/bin/bash
# [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb) [15] Centrality 30-40
#module load ${HOME}/.bashrcROOT
module load root
#.  ${ROOTSYS}/bin/thisroot.sh
i=$1           #i=0-5 pp,   0-1 pbpb/pp
M=$2           #0=pp, 1-2=pp-pbpb
K=$3           # pp, Type5_alphas_0p3_Q0_2GeV
C=$4           #("0-10"  "0-5" "5-10" "10-20" "20-30"  "30-40" "40-50" "50-60" "20-40" )
EtaJetMIN=$5    #0.0
EtaJetMINS=$6  #0p0
EtaJetMAX=$7   #2.8
EtaJetMAXS=$8  #2p8
pTMIN=$9 #pTTrackMin =0.2 0.5 1.0
pTMINS=${10} #0p0 0p2 0p5 0p8 1p0
HadronDirectoryFullPath=${11}
JetRADIUS=${12}
JetRADIUSS=${13}
EXP=${14}
EtaCutFlAG=${15}
PPorPbPb=${16}             # pp or pPb or ppPbPb or PbPb or dAu 
ECM=${17}
SoftV2Flag=${18}  # 0= without soft v2, 1=with soft v2 

#PPorPbPb=(pp pp PbPb)     # pp or pp or PbPb
BGS=(-1  0   1)           # -1 or 0 or  1
BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
Centrality=("0-10"  "0-5" "5-10" "10-20" "20-30"  "30-40" "40-50" "50-60" "20-40" )
#ElossType=("alphas_0p25_Q0_2p0")
ElossType=${K}
ElossName=(""   "_Cent_${Centrality[$C]}_Matter_${ElossType}"   "_Cent_${Centrality[$C]}_MatterLBT_withRecoil_${ElossType}" ) 
JetRadius=${JetRADIUS}
JetRadiusS=${JetRADIUSS}
           MC=("Jetscape"  "Jetscape"  "Jetscape"  "Jetscape"  "Pythia"    "Pythia" )   ## Jetscape or Pythia
       Matter=("Hadron"    "Parton"    "Hadron"    "Parton"    "Hadron"    "Parton")      ## Hadron or Parton
  ChargedFlag=(    0          0           0         0           0          0)     ## for charged particle 1, for all particle 0         
ChargedFlagName=("all"     "all"       "all"       "all"       "all"       "all") ## either all or Charged
Hadronization=("Colorless" "Colorless" "Colored"   "Colored"   "Pythia"    "Pythia")
Exp=${EXP}       ## CMS or Atlas or Alice or Star
Ecm=${ECM}        ## 200 or 2760 or 5020 or 7000 GeV

EtaJetMin=${EtaJetMIN}
EtaJetMinS=${EtaJetMINS}
EtaJetMax=${EtaJetMAX}
EtaJetMaxS=${EtaJetMAXS}
EtaCutFlag=${EtaCutFlAG}      ## for |y| 0, for eta 1

FileNameOutPut=${Hadronization[$i]}JetV2_R_${JetRadiusS}_AntikT_${PPorPbPb}_${MC[$i]}_${Exp}_${Ecm}GeV_EtaJetMin_${EtaJetMinS}_EtaJetMax_${EtaJetMaxS}_FS_${ChargedFlagName[$i]}${Matter[$i]}${ElossName[$M]}${BGSstring[$M]}_pTTrackMin_${pTMINS}GeVSoftV2Flag${SoftV2Flag} 

EtaTrackMax=4.0         # Default is 4.0
pTTrackMin=${pTMIN}     # Default is 0.01 GeV
pTTrackMinS=${pTMINS}     # Default is 0.01 GeV 
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV  Pythia7000GeV  # PPColorlessData/Jetscape7000GeV                      
# /wsu/home/fy/fy41/fy4125   RUN/PP5TeVColored/OutPutFile RUN/PP5TeVColorless/OutPutFile   RUN/PP5TeV${Hadronization[i]}/OutPutFile    
# HadronDirectoryRun=( /wsu/home/fy/fy41/fy4125/RUN/PP5TeV${Hadronization[$i]}/OutPutFile /wsu/home/fy/fy41/fy4125/RUN/PP5TeV${Hadronization[$i]}/OutPutFile  /wsu/tmp/AmitJETSCAPE/PbPb5TeV${ElossName[$M]} ) 
HadronDirectoryRun=(${HadronDirectoryFullPath}   ${HadronDirectoryFullPath}   ${HadronDirectoryFullPath})

RunDIR=/wsu/home/fy/fy41/fy4125/RUN/Analysis
cd ${RunDIR}
File=${RunDIR}/Log/ErrNodeJet${FileNameOutPut}.dat
Exec=JetpTV2 
#echo "JOB INFO: Queue=" $PBS_O_QUEUE ",  NodeName="$HOSTNAME   ",  JobName="$PBS_JOBNAME  ",  JobID="$PBS_JOBID >> ${File}
echo "Command: " $Exec  ${HadronDirectoryRun[$M]}  ${FileNameOutPut} ${Exp}  ${Ecm}  ${JetRadius}  ${EtaJetMin}  ${EtaJetMax} ${EtaTrackMax}  ${pTTrackMin} ${EtaCutFlag}  ${ChargedFlag[$i]}  ${MC[$i]}  ${Matter[$i]} ${BGS[$M]} ${Centrality[$C]} ${PPorPbPb} ${SoftV2Flag}
echo "Command: " $Exec  ${HadronDirectoryRun[$M]}  ${FileNameOutPut} ${Exp}  ${Ecm}  ${JetRadius}  ${EtaJetMin}  ${EtaJetMax} ${EtaTrackMax}  ${pTTrackMin} ${EtaCutFlag}  ${ChargedFlag[$i]}  ${MC[$i]}  ${Matter[$i]} ${BGS[$M]} ${Centrality[$C]} ${PPorPbPb} ${SoftV2Flag}  >> ${File}
echo ":upto here" >> ${File}
# make JetpTCrossSectionPartonHadronPbPb >> ${File}
 
./$Exec  ${HadronDirectoryRun[$M]}  ${FileNameOutPut} ${Exp}  ${Ecm}  ${JetRadius}  ${EtaJetMin}  ${EtaJetMax} ${EtaTrackMax}  ${pTTrackMin} ${EtaCutFlag}  ${ChargedFlag[$i]}  ${MC[$i]}  ${Matter[$i]} ${BGS[$M]} ${Centrality[$C]}  ${PPorPbPb} ${SoftV2Flag}   > /dev/null

# [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb) [15] Centrality [16] Species= pp or PbPb pPb   [17] SoftV2Flag (0= without soft v2, 1=with soft v2)
