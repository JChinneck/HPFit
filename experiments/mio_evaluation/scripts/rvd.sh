#!/bin/bash
#SBATCH --job-name rvd
#SBATCH --cpus-per-task=14
#SBATCH --mem 16G
#SBATCH --partition cpu-small
##SBATCH --output $HOME/HPFit/experiments/cb_evaluation/results/evaluation/olive/log/slurm-%A_%a.out
#SBATCH --array=1-60

module load matlab/R2023b
module load R/4.3.1

GRB_LICENSE_FILE=$HOME/gurobi.lic

FOLNAME=rvd # folder where data is
EXP=mio_evaluation # experiment
JOBNAME=miostarts # name of job on SGE
TIMELIMIT=60
Q=0.50
DEP_VAR=TRUE

SRCLOC=$HOME/HPFit/experiments/src
DATALOC=$HOME/HPFit/experiments/$EXP/data/$FOLNAME
RESLOC=$HOME/HPFit/experiments/$EXP/results/$JOBNAME/$FOLNAME
SEEDFILE=$DATALOC/$FOLNAME.in
MOSEKLOC=$HOME/src/mosek/9.3/toolbox/r2015a
GUROBILOC=$HOME/opt/gurobi1100/linux64/matlab
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SLURM_ARRAY_TASK_ID

SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/opt/matlab2023b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_mio.R\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio1\",        $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio-bm\",      $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio1\",   $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio1\",    $TIMELIMIT, \"$RESLOC\", TRUE,  \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\",  $TIMELIMIT, \"$RESLOC\", TRUE,  \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio1\",    $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in
echo "set.seed(12345)" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio-bm\",  $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\", \"$GUROBILOC\")" >> $RESLOC/hyper.$ID.in

/opt/R-4.3.1/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/log/$ID.Rout
rm $HOME/HPFit/experiments/$EXP/scripts/slurm-${SLURM_ARRAY_JOB_ID}_$ID.out
