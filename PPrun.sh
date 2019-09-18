#!/bin/bash

## [1]=JetRadius,  [2]=Jetscape or Pythia, [3]=CMS or Atlas or Alice or Star, [4]= 200 or 2760 or 7000, [5]=EtaJetMin, [6]=EtaJetMax, [7]EtaCutFlag, [8]=ChargedFlag=0 or 1 (all or charged), [9]Hadron or Parton, [10]FileNameOutPut  [11] EtaTrackMax [12] pTTrackMin [13]HadronDirectory

for (( i=3;i<4;i++ )) #i=0, 1, 2, 3
do
Trial=$i
MC=("Jetscape" "Jetscape" "Pythia" "Pythia")   ## Jetscape or Pythia
Exp=Atlas       ## CMS or Atlas or Alice or Star
Ecm=5020       ## 200 or 2760 or 5020 or 7000 GeV

JetRadius=0.4
JetRadiusS=0p4

EtaJetMin=0.0
EtaJetMinS=0p0
EtaJetMax=2.8
EtaJetMaxS=2p8
EtaCutFlag=0      ## for |y| 0, for eta 1
ChargedFlag=0     ## for charged particle 1, for all particle 0

Matter=("Hadron"  "Parton"  "Hadron"  "Parton")      ## Hadron or Parton

ChargedFlagName=all  ## either all or Charged

FileNameOutPut=JetCrossSection_R_${JetRadiusS}_AntikT_pp_${MC[$i]}_${Exp}_${Ecm}GeV_EtaJetMin_${EtaJetMinS}_EtaJetMax_${EtaJetMaxS}_FS_${ChargedFlagName}${Matter[$i]}

EtaTrackMax=3.0         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV  Pythia7000GeV
# PPColorlessData/Jetscape7000GeV
# /wsu/home/fy/fy41/fy4125   RUN/PP5TeVColored/OutPutFile RUN/PP5TeVColorless/OutPutFile   RUN/PythiaRun/OutPutFileHadron
HadronDirectoryRun=/wsu/home/fy/fy41/fy4125/RUN/PythiaRun/OutPutFileParton
    echo "Working on i=" $i
    #LogFile=/wsu/home/fy/fy41/fy4125/RUN/FastJet/Log/Out_${FileNameOutPut}_${i}.dat
    ErrFile=/wsu/home/fy/fy41/fy4125/RUN/FastJet/Log/Err_${FileNameOutPut}_${i}.dat  
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetpTCrossSectionPartonHadronPP

qsub -V -q accq   -l mem=2gb  -N ${EtaJetMaxS}J${Trial}  -o $LogFile -e $ErrFile --  $Exec  ${JetRadius}  ${MC[$i]}   ${Exp}  ${Ecm}         ${EtaJetMin}  ${EtaJetMax}   ${EtaCutFlag}  ${ChargedFlag}      ${Matter[$i]}  ${FileNameOutPut}    ${EtaTrackMax}  ${pTTrackMin}   ${HadronDirectoryRun}

done
