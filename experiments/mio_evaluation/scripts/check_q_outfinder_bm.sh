#!/bin/bash
#shell for qsub to use
#$ -S /bin/bash
# name for the job on qstat
#$ -N check_q_outfinder
# tell SGE that it's an array and job numbers
#$ -t 1-70:1
# tell SGE to run at most 2 jobs at once
#$ -tc 2

FOLNAME=bm # folder where data is
EXP=mio_evaluation # experiment
JOBNAME=check_q_outfinder # name of job on SGE
TIMELIMIT=3600
DEP_VAR=TRUE

SRCLOC=$HOME/HPFit/experiments/src
DATALOC=$HOME/HPFit/experiments/$EXP/data/$FOLNAME
RESLOC=$HOME/HPFit/experiments/$EXP/results/$JOBNAME/$FOLNAME
SEEDFILE=$DATALOC/$FOLNAME.in
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SGE_TASK_ID

Q=0
echo $Q


SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/usr/local/MATLAB/R2018b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_mio.R\")" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbmio3\", $TIMELIMIT, \"$RESLOC\", FALSE)" >> $RESLOC/hyper.$ID.in

/usr/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/log/$ID.Rout
rm $HOME/$JOBNAME.e*
rm $HOME/$JOBNAME.o*
