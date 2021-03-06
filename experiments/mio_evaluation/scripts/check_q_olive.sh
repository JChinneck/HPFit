#!/bin/bash
#shell for qsub to use
#$ -S /bin/bash
# name for the job on qstat
#$ -N check_q
# tell SGE that it's an array and job numbers
#$ -t 1-6:1
# tell SGE to run at most 3 jobs at once
#$ -tc 3

FOLNAME=olive # folder where data is
EXP=mio_evaluation # experiment
JOBNAME=check_q # name of job on SGE
TIMELIMIT=3600
DEP_VAR=TRUE

SRCLOC=$HOME/hyperplane_fitting/src
DATALOC=$HOME/hyperplane_fitting/$EXP/data/$FOLNAME
RESLOC=$HOME/hyperplane_fitting/$EXP/results/$JOBNAME/$FOLNAME
SEEDFILE=$DATALOC/$FOLNAME.in
mkdir -p $RESLOC # make folder for results
mkdir -p $RESLOC/log # make folder for logs

ID=$SGE_TASK_ID
  
case $ID in

    1)
        Q=0.99
        ;;
    2)
        Q=0.974
        ;;
    3)
        Q=0.814
        ;;
    4)
        Q=0.943
        ;;
    5)
        Q=0.75
        ;;
    6)
        Q=0.8
        ;;
esac

echo $Q


SEED=$(sed -n -e "$ID p" $SEEDFILE)
echo "library(MASS)" > $RESLOC/hyper.$ID.in
echo "library(matlabr)" >> $RESLOC/hyper.$ID.in
echo "options(matlab.path='/usr/local/MATLAB/R2018b/bin')" >> $RESLOC/hyper.$ID.in
echo "source(\"$SRCLOC/run_mio.R\")" >> $RESLOC/hyper.$ID.in
echo "run_mio(\"$DATALOC\", \"$SRCLOC\", \"$SEED\", $Q, $DEP_VAR, \"mio1\", $TIMELIMIT, \"$RESLOC\", TRUE)" >> $RESLOC/hyper.$ID.in

/usr/bin/R CMD BATCH $RESLOC/hyper.$ID.in $RESLOC/log/$ID.Rout

rm $RESLOC/hyper.$ID.in
#rm $RESLOC/log/$ID.Rout
rm $HOME/$JOBNAME.e*
rm $HOME/$JOBNAME.o*
