#!/bin/bash
#shell for qsub to use
#$ -S /bin/bash
# name for the job on qstat
#$ -N evaluation
# tell SGE that it's an array and job numbers
#$ -t 1-10:1
# tell SGE to run at most 3 jobs at once
#$ -tc 3

FOLNAME=bm-nox # folder where data is
EXP=rbm_evaluation # experiment
JOBNAME=evaluation # name of job on SGE and location of results
TIMELIMIT=300 # used for RBM-MIO3
Q=0.50
DEP_VAR=TRUE

SRCLOC=$HOME/hyperplane_fitting/src
DATALOC=$HOME/hyperplane_fitting/$EXP/data/$FOLNAME
RESLOC=$HOME/hyperplane_fitting/$EXP/results/$JOBNAME/$FOLNAME
MOSEKLOC=$HOME/src/mosek/9.3/toolbox/r2015a
SEEDFILE=$DATALOC/$FOLNAME.in
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SGE_TASK_ID
  


SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/usr/local/MATLAB/R2018b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/mpack.txt\")" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_competitors.R\")" >> $RESLOC/hyper.$ID.in
echo "get_dists(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in

/usr/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

#rm $RESLOC/hyper.$ID.in
#rm $RESLOC/$ID.Rout
rm $HOME/$JOBNAME.e*
rm $HOME/$JOBNAME.o*
