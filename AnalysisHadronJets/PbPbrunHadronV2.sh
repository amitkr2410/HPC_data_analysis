#!/bin/bash
#[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,      [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10] BGS
module load root
i=$1           #i=0-5 pp,   0-1 pbpb/pp                                                              
M=$2           #0=pp, 1-2=pp-pbpb          
K=$3           #eloss type pp, Type5_alphas_0p3_Q0_2GeV     
C=$4           #Centrality=("10-20" "20-30"   "30-40" "40-50" "50-60"  ) 
HadronDirectory=$5
SoftV2Flag=$6  #0 or 1, off or on

ElossType=${K}
PPorPbPb=(pp pp PbPb)     # pp or pp or PbPb                                                          
BGS=(-1  0   1)           # -1 or 0 or  1                               
BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
Centrality=("10-20" "20-30"   "30-40" "40-50" "50-60"  )
#ElossType=("0p15_Q0_1p0"  "0p20_Q0_1p0"  "0p20_Q0_1p4"  "0p20_Q0_2p0" "0p25_Q0_1p0"    "0p25_Q0_1p4"   "0p25_Q0_2p0"  "0p25_Q0_3p0"  "0p25_Q0_10p0")
ElossName=(""   ""   "_Cent_${Centrality[$C]}_MatterLBT_withRecoil_alphas_${ElossType}" )
Exp=Atlas        ## CMS or ATLAS or ALICE or STAR
Ecm=5020       ## 200 or 2760 or 5020, 7000 GeV

SingleHadronEtaCut=2.5
SingleHadronEtaCutS=2p5
EtaCutFlag=1      ## for |y| 0, for eta 1
             MC=("Jetscape" "Jetscape" "Jetscape"  "Jetscape" "Pythia" "Pythia" )   ## Jetscape or Pythia
         Matter=("Hadron"  "Parton"  "Hadron"  "Parton"  "Hadron"   "Parton")      ## Hadron or Parton 
    ChargedFlag=(1          0         1         0         1         0)     ## for charged particle 1, for all particle 0
ChargedFlagName=("Charged"  "all"     "Charged" "all"     "Charged" "all") ## either all or Charged
Hadronization=("Colorless" "Colorless" "Colored" "Colored" "Pythia" "Pythia")
FileNameOutPut=${Hadronization[$i]}ChargedParticleV2_${PPorPbPb[$M]}_${MC[$i]}_${Exp}_${Ecm}GeV_EtaMax_${SingleHadronEtaCutS}_FS_${ChargedFlagName[$i]}${Matter[$i]}${ElossName[$M]}_SoftV2Flag${SoftV2Flag}${BGSstring[$M]}
EtaTrackMax=4.0         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV  Pythia7000GeV  # PPColorlessData/Jetscape7000GeV
# /wsu/home/fy/fy41/fy4125   RUN/PP5TeVColored/OutPutFile RUN/PP5TeVColorless/OutPutFile   RUN/PP5TeV${Hadronization[i]}/OutPutFile
#HadronDirectoryRun=(/wsu/home/groups/maj-shen/JETSCAPEDataFile/PPColorlessData/Jetscape2760GeV/JetscapePPColorless2760GeV  /wsu/home/groups/maj-shen/JETSCAPEDataFile/PPColorlessData/Jetscape2760GeV/JetscapePPColorless2760GeV   /wsu/tmp/AmitJETSCAPE2TeV/PbPb2TeV_Cent_0-5_MatterLBT_withRecoil_alphas_0p25_Q0_2p0)
HadronDirectoryRun=(${HadronDirectory}  ${HadronDirectory}  ${HadronDirectory}  )
RunDIR=/wsu/home/fy/fy41/fy4125/RUN/Analysis
cd ${RunDIR}
File=${RunDIR}/Log/ErrNodeHadronV2${FileNameOutPut}.dat

echo "JOB INFO: Queue=" $PBS_O_QUEUE ",  NodeName="$HOSTNAME   ",  JobName="$PBS_JOBNAME  ",  JobID="$PBS_JOBID >> ${File}
echo "AnalysisSpectraV2"    ${HadronDirectoryRun[$M]} ${FileNameOutPut} ${Exp}      ${Ecm}        ${SingleHadronEtaCut}    ${EtaCutFlag}  ${ChargedFlag[$i]}     ${MC[$i]}   ${Matter[$i]}     ${BGS[$M]} .${Centrality[$C]}  ${SoftV2Flag}

./AnalysisSpectraV2   ${HadronDirectoryRun[$M]} ${FileNameOutPut} ${Exp}      ${Ecm}        ${SingleHadronEtaCut}    ${EtaCutFlag}  ${ChargedFlag[$i]}     ${MC[$i]}  ${Matter[$i]}  ${BGS[$M]} ${Centrality[$C]}  ${SoftV2Flag}   > /dev/null

#done

#[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,      [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10]  BGS    [11] Centrality =0-10  [12] SoftV2Flag=0 or 1 (off, on)       
