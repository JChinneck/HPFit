#!/bin/bash
#shell for qsub to use
#$ -S /bin/bash
# name for the job on qstat
#$ -N comparison
# tell SGE that it's an array and job numbers

FOLNAME=bm # folder where data is
EXP=mio_evaluation # experiment
JOBNAME=comparison # name of job on SGE
TIMELIMIT=3600
Q=0.50
DEP_VAR=TRUE

SRCLOC=$HOME/HPFit/experiments/src
DATALOC=$HOME/HPFit/experiments/$EXP/data/$FOLNAME
RESLOC=$HOME/HPFit/experiments/$EXP/results/$JOBNAME/$FOLNAME
SEEDFILE=$DATALOC/$FOLNAME.in
MOSEKLOC=$HOME/src/mosek/9.3/toolbox/r2015a
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

#ID=$SGE_TASK_ID
ID=1


#SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/usr/local/MATLAB/R2018b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_mio.R\")" >> $RESLOC/hyper.$ID.in
SEED=m3001n11m_outliers2000i33.csv
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
SEED=m3001n11m_outliers2000i34.csv
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
SEED=m3001n11m_outliers2000i37.csv
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
SEED=m3001n11m_outliers2000i39.csv
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in

#SEED=m6001n21m_outliers4000i63.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i64.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i65.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i66.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio1\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio1\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i67.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"alg3-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"cbq-mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i68.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio-bm\", $TIMELIMIT, \"$RESLOC\", FALSE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in
#SEED=m6001n21m_outliers4000i69.csv
#echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"lqs-mio-bm\", $TIMELIMIT, \"$RESLOC\", TRUE, \"$MOSEKLOC\")" >> $RESLOC/hyper.$ID.in

/usr/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/log/$ID.Rout
rm $HOME/$JOBNAME.e*
rm $HOME/$JOBNAME.o*
