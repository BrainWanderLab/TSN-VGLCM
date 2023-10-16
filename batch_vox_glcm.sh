#! /bin/bash
source /etc/profile
#$ -S /bin/bash # SGE shell type


mT1datafile=$1
GMdatafile=$2
WMdatafile=$3
outdir=$4
subj=$5
toolpath=$6
logpath=${outdir}/log_vglcm/${subj}/pipeline_vglcm.log
mkdir ${outdir}/log_vglcm/${subj}/ -p
resultdir=${outdir}/${subj}
if [ ! -d ${resultdir} ];then mkdir ${resultdir} -p; fi
#######################
tooldir=${toolpath}
est1="addpath('${tooldir}');Batch_VoxelTexture('${mT1datafile}','${GMdatafile}','${WMdatafile}','${resultdir}');exit"
cmd="matlab -nosplash -nodesktop -r \"$est1\""
script=/tmp/VGLCM_pipe_`date +%Y-%m-%d-%h-%m-%s-%N`
echo ${cmd} >$script
chmod +x $script
$script >> $logpath 2>&1
