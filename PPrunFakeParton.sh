#!/bin/bash

## [1]=JetRadius,  [2]=Jetscape or Pythia, [3]=CMS or Atlas or Alice or Star, [4]= 200 or 2760 or 7000, [5]=EtaJetMin, [6]=EtaJetMax, [7]EtaCutFlag, [8]=ChargedFlag=0 or 1 (all or charged), [9]Hadron or Parton, [10]FileNameOutPut  [11] EtaTrackMax [12] pTTrackMin 

FakePartonPz=(20 40 60 80 100 160 200)
for (( i=0;i<7;i++ )) #i=0, 1, 2, 3
do
MC=Jetscape
Exp=Star       ## CMS or Atlas or Alice or Star
Ecm=200       ## 200 or 2760 or 7000 GeV

JetRadius=0.6
JetRadiusS=0p6

EtaJetMin=0.0
EtaJetMinS=0p0
EtaJetMax=0.5
EtaJetMaxS=0p5
EtaCutFlag=1      ## for |y| 0, for eta 1
ChargedFlag=0     ## for charged particle 1, for all particle 0

Matter=Hadron      ## Hadron or Parton

ChargedFlagName=all  ## either all or Charged
Hadronization=Colorless ## Colored or Colorless

FileNameOutPut=JetCrossSection_R_${JetRadiusS}_AntikT_pp_${MC}_${Exp}_${Ecm}GeV_EtaJetMin_${EtaJetMinS}_EtaMax_${EtaJetMaxS}_FS_${ChargedFlagName}${Matter}${Hadronization}FakePartonPz${FakePartonPz[$i]}

EtaTrackMax=2.5         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV
# PPColorlessData/Jetscape7000GeV
HadronDirectoryRun=FakeParton${Hadronization}200GeV/JetscapeHadronListFakePartonPz${FakePartonPz[$i]}.out

    echo "Working on i=" $i
    LogFile=/wsu/home/fy/fy41/fy4125/RUN/FastJet/Log/Out${i}.dat
    ErrFile=/wsu/home/fy/fy41/fy4125/RUN/FastJet/Log/Err${i}.dat 

    #LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetpTCrossSectionPartonHadronPP

qsub -V -q wsuq   -l mem=2gb  -N NoC${i}  -o $LogFile -e $ErrFile --  $Exec  ${JetRadius}  ${MC}   ${Exp}  ${Ecm}         ${EtaJetMin}  ${EtaJetMax}   ${EtaCutFlag}  ${ChargedFlag}      ${Matter}  ${FileNameOutPut}    ${EtaTrackMax}  ${pTTrackMin}   ${HadronDirectoryRun}

done
