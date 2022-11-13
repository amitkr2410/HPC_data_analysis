#!/bin/bash
# [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb)
i=$1           #i=0-5 pp,   0-1 pbpb/pp
M=$2           #0=pp, 1-2=pp-pbpb
K=$3           #eloss type K=0a, 1b, 2c, 3d, 4e
PPorPbPb=(pp pp PbPb)     # pp or pp or PbPb
BGS=(-1  0   1)           # -1 or 0 or  1
BGSstring=(""  "_BGS${BGS[$M]}"   "_BGS${BGS[$M]}")
ElossType=("0p20_Q0_1p0"  "0p20_Q0_1p4"  "0p20_Q0_2p0" "0p25_Q0_1p0"  "0p25_Q0_1p4" )
ElossName=(""   ""   "_Cent_0-10_MatterLBT_withRecoil_alphas_${ElossType[$K]}" ) 
JetRadius=0.4
JetRadiusS=0p4
           MC=("Jetscape"  "Jetscape"  "Jetscape"  "Jetscape"  "Pythia"    "Pythia" )   ## Jetscape or Pythia
       Matter=("Hadron"    "Parton"    "Hadron"    "Parton"    "Hadron"    "Parton")      ## Hadron or Parton
  ChargedFlag=(    0          0           0         0           0          0)     ## for charged particle 1, for all particle 0         
ChargedFlagName=("all"     "all"       "all"       "all"       "all"       "all") ## either all or Charged
Hadronization=("Colorless" "Colorless" "Colored"   "Colored"   "Pythia"    "Pythia")
Exp=Atlas       ## CMS or Atlas or Alice or Star
Ecm=5020        ## 200 or 2760 or 5020 or 7000 GeV

EtaJetMin=0.0
EtaJetMinS=0p0
EtaJetMax=2.8 
EtaJetMaxS=2p8 
EtaCutFlag=0      ## for |y| 0, for eta 1

FileNameOutPut=${Hadronization[$i]}JetCrossSection_R_${JetRadiusS}_AntikT_${PPorPbPb[$M]}_${MC[$i]}_${Exp}_${Ecm}GeV_EtaJetMin_${EtaJetMinS}_EtaJetMax_${EtaJetMaxS}_FS_${ChargedFlagName[$i]}${Matter[$i]}${ElossName[$M]}${BGSstring[$M]} 

EtaTrackMax=3.0         # Default is 2.5
pTTrackMin=0.01          # Default is 0.01
# PPPaperData/JETSCAPE_200GeV  Pythia_200GeV  Pythia7000GeV  # PPColorlessData/Jetscape7000GeV                      
# /wsu/home/fy/fy41/fy4125   RUN/PP5TeVColored/OutPutFile RUN/PP5TeVColorless/OutPutFile   RUN/PP5TeV${Hadronization[i]}/OutPutFile    
HadronDirectoryRun=( /wsu/home/fy/fy41/fy4125/RUN/PP5TeV${Hadronization[$i]}/OutPutFile /wsu/home/fy/fy41/fy4125/RUN/PP5TeV${Hadronization[$i]}/OutPutFile  /wsu/tmp/AmitJETSCAPE/PbPb5TeV${ElossName[$M]} ) 
RunDIR=/wsu/home/fy/fy41/fy4125/RUN/FastJet
cd ${RunDIR}
File=${RunDIR}/Log/ErrNode${FileNameOutPut}.dat
Exec=JetpTCrossSectionPartonHadronPbPb
echo "JOB INFO: Queue=" $PBS_O_QUEUE ",  NodeName="$HOSTNAME   ",  JobName="$PBS_JOBNAME  ",  JobID="$PBS_JOBID >> ${File}
echo "Command: " $Exec  ${HadronDirectoryRun[$M]}  ${FileNameOutPut} ${Exp}  ${Ecm}  ${JetRadius}  ${EtaJetMin}  ${EtaJetMax} ${EtaTrackMax}  ${pTTrackMin} ${EtaCutFlag}  ${ChargedFlag[$i]}  ${MC[$i]}  ${Matter[$i]} ${BGS[$M]} >> ${File}
echo ":upto here" >> ${File}

./$Exec  ${HadronDirectoryRun[$M]}  ${FileNameOutPut} ${Exp}  ${Ecm}  ${JetRadius}  ${EtaJetMin}  ${EtaJetMax} ${EtaTrackMax}  ${pTTrackMin} ${EtaCutFlag}  ${ChargedFlag[$i]}  ${MC[$i]}  ${Matter[$i]} ${BGS[$M]} > /dev/null

# [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb)
