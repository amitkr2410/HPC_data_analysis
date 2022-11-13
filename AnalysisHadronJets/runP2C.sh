#!/bin/bash

#[1]=Jetscape or Pythia, [2]=CMS or Atlas or Alice or Star, [3]= 200 or 2760 or 7000, [4]EtaCutFlag, [5]=ChargedFlag=0 or 1 (all or charged), [6]=Hadron or Parton, [7]=EtaTrackMax, [8]=pTTrackMin, [9]=pTTrigMin [10]=pTTrigMax [11]=pTAssMin  [12]=pTAssMax [13]=FileNameOutPut  [14]=HadronDirectoryRun [15] NumberOfFiles=10

MyDIR=/wsu/home/fy/fy41/fy4125/RUN/FastJet
cd ${MyDIR}
for (( i=1;i<2;i++ ))
do
MC=("Jetscape" "Jetscape" "Pythia" "Pythia")   ## Jetscape or Pythia
Exp=CMS       ## CMS or Atlas or Alice or Star
Ecm=7000        ## 200 or 2760 or 7000 GeV

pTTrigMin=10.0
pTTrigMinS=10p0
pTTrigMax=100.0 
pTTrigMaxS=100p0 
pTAssMin=${pTTrigMin}
pTAssMinS=${pTTrigMinS}
pTAssMax=${pTTrigMax}
pTAssMaxS=${pTTrigMaxS}

EtaCutFlag=1      ## for |y| 0, for eta 1
ChargedFlag=0     ## for charged particle 1, for all particle 0
ChargedFlagName=all  ## either all or Charged

Matter=("Hadron"  "Parton"  "Hadron"  "Parton")      ## Hadron or Parton

Cent=0_5      ## 0_5  or 20_30  or 30_40
ElossModule=MatterLBT0p2_Q0_2     #MatterinvacLBT0p25  MatterinvacLBT0p3 MatterinvacLBT0p4  MatterLBT_alphas_0p25_Q0_2  MatterQhat0p5_Q0_1 MatterQhat1p0_Q0_1
#FileNameOutPut=YieldDeltaEtaDeltaPhi_PbPb_${MC[$i]}_${Exp}_${Ecm}GeV_pTTrigMin_${pTTrigMinS}_pTTrigMax_${pTTrigMaxS}_FS_${ChargedFlagName}${Matter[$i]}_Cent${Cent}_${ElossModule}
FileNameOutPut=OutputP2CorrelationMinpT10p0Eta2p0
EtaTrackMax=2.0         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
HadronDirectoryRun=TwoStageHydro/Data/JetscapeHadronList
NumberOfFiles=104

#qsub -V -q wsuq   -l mem=2gb  -N 2pC_${Trial}  -o $LogFile -e $ErrFile --  $Exec  ${MC[$i]}   ${Exp}  ${Ecm}   ${EtaCutFlag}  ${ChargedFlag}      ${Matter[$i]}     ${EtaTrackMax}  ${pTTrackMin}  ${pTTrigMin} ${pTTrigMax}   ${pTAssMin}  ${pTAssMax} ${FileNameOutPut}  ${HadronDirectoryRun}
echo "./P2Correlation  ${MC[$i]}   ${Exp}  ${Ecm}   ${EtaCutFlag}  ${ChargedFlag}      ${Matter[$i]}     ${EtaTrackMax}  ${pTTrackMin}  ${pTTrigMin} ${pTTrigMax}   ${pTAssMin}  ${pTAssMax} ${FileNameOutPut}  ${HadronDirectoryRun}" 

./P2Correlation  ${MC[$i]}   ${Exp}  ${Ecm}   ${EtaCutFlag}  ${ChargedFlag}      ${Matter[$i]}     ${EtaTrackMax}  ${pTTrackMin}  ${pTTrigMin} ${pTTrigMax}   ${pTAssMin}  ${pTAssMax}    ${FileNameOutPut}  ${HadronDirectoryRun}   ${NumberOfFiles}
done


