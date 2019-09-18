#!/bin/bash

#[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,      [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10] BGS

#for (( i=0; i<6; i++ )) #i=0, 1,     2, 3       4,5

i=$1            #i=0, 1,     2, 3       4,5
PPorPbPb=pp     # pp or PbPb
BGS=0             # 0 or  1

MC=("Jetscape" "Jetscape" "Jetscape"  "Jetscape" "Pythia" "Pythia" )   ## Jetscape or Pythia
Exp=CMS        ## CMS or ATLAS or ALICE or STAR
Ecm=5020       ## 200 or 2760 or 5020, 7000 GeV

SingleHadronEtaCut=1.0
SingleHadronEtaCutS=1p0
EtaCutFlag=1      ## for |y| 0, for eta 1
         Matter=("Hadron"  "Parton"  "Hadron"  "Parton"  "Hadron"   "Parton")      ## Hadron or Parton 
    ChargedFlag=(1          0         1         0         1         0)     ## for charged particle 1, for all particle 0
ChargedFlagName=("Charged"  "all"     "Charged" "all"     "Charged" "all") ## either all or Charged

Hadronization=("Colorless" "Colorless" "Colored" "Colored" "Pythia" "Pythia")
FileNameOutPut=${Hadronization[$i]}SingleHadronYield_${PPorPbPb}_${MC[$i]}_${Exp}_${Ecm}GeV_EtaMax_${SingleHadronEtaCutS}_FS_${ChargedFlagName[$i]}${Matter[$i]}
EtaTrackMax=3.0         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV  Pythia7000GeV  # PPColorlessData/Jetscape7000GeV
# /wsu/home/fy/fy41/fy4125   RUN/PP5TeVColored/OutPutFile RUN/PP5TeVColorless/OutPutFile   RUN/PP5TeV${Hadronization[i]}/OutPutFile
HadronDirectoryRun=/wsu/home/fy/fy41/fy4125/RUN/PP5TeV${Hadronization[$i]}/OutPutFile
RunDIR=/wsu/home/fy/fy41/fy4125/RUN/FastJet
cd ${RunDIR}
File=${RunDIR}/Log/ErrNode${FileNameOutPut}.dat

echo "JOB INFO: Queue=" $PBS_O_QUEUE ",  NodeName="$HOSTNAME   ",  JobName="$PBS_JOBNAME  ",  JobID="$PBS_JOBID >> ${File}

./AnalysisSpectraSingleHadron  $Exec  ${HadronDirectoryRun} ${FileNameOutPut} ${Exp}      ${Ecm}        ${SingleHadronEtaCut}    ${EtaCutFlag}  ${ChargedFlag[$i]}     ${MC[$i]}  ${Matter[$i]}  ${BGS} > /dev/null

#done

#[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,      [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10]  BGS
