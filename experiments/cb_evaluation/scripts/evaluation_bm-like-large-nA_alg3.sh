#!/bin/bash
#shell for qsub to use
#$ -S /bin/bash
# name for the job on qstat
#$ -N evaluation
# tell SGE that it's an array and job numbers
#$ -t 1-40:1
# tell SGE to run at most 3 jobs at once
#$ -tc 3
# hard wall clock time of 1 hour
#$ -l h_rt=1:00:00

FOLNAME=bm-like-large-nA # folder where data is
EXP=cb_evaluation # experiment
JOBNAME=evaluation # name of job on SGE
TIMELIMIT=60 # used for CB-MIO3
Q=0.50
DEP_VAR=TRUE

SRCLOC=$HOME/HPFit/experiments/src
DATALOC=$HOME/HPFit/experiments/$EXP/data/$FOLNAME
RESLOC=$HOME/HPFit/experiments/$EXP/results/$JOBNAME/$FOLNAME
MOSEKLOC=$HOME/src/mosek/9.3/toolbox/r2015a
SEEDFILE=$DATALOC/$FOLNAME.in
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SGE_TASK_ID
  


SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "library(robust)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/usr/local/MATLAB/R2018b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/mpack.txt\")" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_competitors.R\")" >> $RESLOC/hyper.$ID.in

#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_cb(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "get_alg3(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_mh(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "hbreg_only(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_dists(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_lts(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_lqs(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_mm(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_rewlse(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in


/usr/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/$ID.Rout
rm $HOME/$JOBNAME.e*
rm $HOME/$JOBNAME.o*
