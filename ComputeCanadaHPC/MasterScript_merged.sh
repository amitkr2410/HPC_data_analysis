#!/bin/bash 
    EcmMaxIndex=(66 51 23)
    Ecm=(5020 2760 200)

    DirInput=$1
    DirOutput=$2
    HadronORParton=$3
    EcmIndex=$4
    queue=$5
    Exec=/scratch/amitkr/JSWork/simple_job_merged${Ecm[${EcmIndex}]}.sh

    for(( Bin=0; Bin<${EcmMaxIndex[${EcmIndex}]}; Bin++ ))
    do        
	    ErrFile=/scratch/amitkr/JSWork/Log/Merged${Ecm[${EcmIndex}]}${HadronORParton}${Bin}.err
	    LogFile=/scratch/amitkr/JSWork/Log/Merged${Ecm[${EcmIndex}]}${HadronORParton}${Bin}.out
	    sbatch --account=rrg-jeon-ac   -N 1 -n 1 --mem=4G --time=00:30:00   --job-name  Merg${HadronORParton}${Bin}   -o ${LogFile} -e ${ErrFile} --  ${Exec} ${DirInput}  ${DirOutput} ${HadronORParton} ${Bin} 
    done

