#!/bin/bash
#SBATCH --job-name ltsbmlikelargenB
#SBATCH --cpus-per-task=1
#SBATCH --mem 4G
#SBATCH --partition cpu-small
#SBATCH -t 0-1:00 # time limit of one hour
##SBATCH --output $HOME/HPFit/experiments/cb_evaluation/results/evaluation/olive/log/slurm-%A_%a.out
#SBATCH --array=1-40

module load matlab/R2023b
module load R/4.3.1

FOLNAME=bm-like-large-nB # folder where data is
EXP=cb_evaluation # experiment
JOBNAME=evaluation # name of job on SGE
TIMELIMIT=60 # used for CB-MIO3
Q=0.50
DEP_VAR=TRUE

SRCLOC=$HOME/HPFit/experiments/src
DATALOC=$HOME/HPFit/experiments/$EXP/data/$FOLNAME
RESLOC=$HOME/HPFit/experiments/$EXP/results/$JOBNAME/$FOLNAME
MOSEKLOC=$HOME/src/mosek/9.3/toolbox/r2015a
GUROBILOC=$HOME/opt/gurobi1100/linux64/matlab
SEEDFILE=$DATALOC/$FOLNAME.in
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SLURM_ARRAY_TASK_ID
  

SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "library(robust)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/opt/matlab2023b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/mpack.txt\")" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_competitors.R\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_alg3(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_mh(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_mm(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_rewlse(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_cb(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "hbreg_only(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\", \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_lm(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "get_lts(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in
#echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
#echo "get_lqs(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, $TIMELIMIT, \"$RESLOC\")" >> $RESLOC/hyper.$ID.in

/opt/R-4.3.1/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID_lts.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/$ID.Rout
rm $HOME/HPFit/experiments/$EXP/scripts/slurm-${SLURM_ARRAY_JOB_ID}_$ID.out
