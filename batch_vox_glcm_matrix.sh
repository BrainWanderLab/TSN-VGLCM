#! /bin/bash
source /etc/profile
#$ -S /bin/bash # SGE shell type


vglcm_datadir=$1
mT1_datadir=$2
MNI_datafile=$3
outdir=$4
subj=$5
toolpath=$6
logpath=${outdir}/log_vglcm_matrix/${subj}/pipeline_vglcm.log
mkdir ${outdir}/log_vglcm_matrix/${subj}/ -p
#######################
tooldir=${toolpath}
est1="addpath('${tooldir}');Batch_radi_simW('${outdir}','${vglcm_datadir}','${mT1_datadir}','${subj}','${MNI_datafile}');exit"
cmd="matlab -nosplash -nodesktop -r \"$est1\""
script=/tmp/VGLCM_pipe_`date +%Y-%m-%d-%h-%m-%s-%N`
echo ${cmd} >$script
chmod +x $script
$script >> $logpath 2>&1
