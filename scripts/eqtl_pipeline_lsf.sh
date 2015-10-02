#!/bin/bash
################################################################
# simple wrapper to run a single instance of eqtl_pipeline in a
# cluster with the lsf job scheduler
debug=1

function stop_job {
    if [ "$DEBUG-" == "1-" ]; then
	echo "stop/suspend job $1"
    else
	bstop -J $1
    fi
}

function resume_job {
    if [ "$DEBUG-" == "1-" ]; then
	echo "resume job $1"
    else
	bresume -J $1
    fi
}

function submit_job {
    local jobname=$1
    local waitforids=$2
    shift
    shift
    local cmd2e=$*
    local ECHO=
    if [ "$DEBUG-" = "1-" ]; then
        ECHO=echo
    fi

    #########################################################
    # limit the number of parallel jobs by using lsf's groups
    
    # default group: irap
    # TODO: define groups by stage: 
    #    irap_qc, ...
    local GROUP=
    if [ "$LSF_GROUP-" != "-" ]; then
	GROUP="-g $LSF_GROUP"
    fi
    #########################################################
    #-R  "span[ptile=$THREADS]"
    local MAX_MEM=`get_maxmem $MEM`
    if [ "$WAIT_FOR_IDS-" != "-no" ]; then
	$ECHO bsub $LSF_PARAMS -q $QUEUE -n $THREADS -R "span[hosts=1]"  -M $MAX_MEM -R "select[mem>=$MEM] rusage[mem=$MEM]"  -w "done($WAIT_FOR_IDS)"  -cwd `pwd` -o "`get_path2logfile`/$jobname-%J.out" -e "`get_path2logfile`/$jobname-%J.err" -J $jobname  $cmd2e max_threads=$THREADS  max_mem=$MEM
    else
	$ECHO bsub $LSF_PARAMS -q $QUEUE  $GROUP -n $THREADS -R "span[hosts=1]"  -M $MAX_MEM -R "select[mem>=$MEM]  rusage[mem=$MEM]"   -cwd `pwd` -o "`get_path2logfile`/$jobname-%J.out" -e "`get_path2logfile`/$jobname-%J.err" -J $jobname  $cmd2e max_threads=$THREADS   max_mem=$MEM	
    fi
}


function submit_jobs {
    # uses ARGS and targets
    local jobname_prefix=$1
    local wait_for_id=$2
    shift
    shift
    local cmd=$*
    let j=1
    for t in $targets; do
	submit_job "$jobname_prefix[$j]" $wait_for_id $cmd $t
	let j=$j+1
    done
}
################################################################
if [ "$*-" == "- " ]; then
    echo "ERROR: missing arguments"
    echo "usage: eqtl_pipeline_lsf.sh conf=file [extra options]"
    exit 1
fi

eqtl_pipeline $* -n -q

if [ $? -eq 0 ]; then
    echo "All done - no need to submit jobs"
    echo 0
fi
ARGS=$*

RAND=`perl -e "print int(rand()*10);"`
DATE=`date "+%w%H%M%S"`
JOBNAME_SUF="$DATE$RAND"

# submit the jobs
submit_job "eqtl0_$JOBNAME_SUF" $wait_for_id $cmd $t
stop_job "eqtl0_$JOBNAME_SUF"

targets=`eqtl_pipeline $* targets1 | tail -n 1`
submit_jobs eqtl1_$JOBNAME_SUF  "eqtl0_$JOBNAME_SUF" eqtl_pipeline $ARGS

targets=`eqtl_pipeline $* targets2 | tail -n 1`
submit_jobs eqtl2_$JOBNAME_SUF "eqtl_$JOBNAME_SUF1[*]" eqtl_pipeline $ARGS

targets=`eqtl_pipeline $* targets3 | tail -n 1`
submit_jobs eqtl3_$JOBNAME_SUF "eqtl2_$JOBNAME_SUF[*]" eqtl_pipeline $ARGS

targets=`eqtl_pipeline $* targets4 | tail -n 1`
submit_jobs eqtl4_$JOBNAME_SUF "eqtl3_$JOBNAME_SUF[*]" eqtl_pipeline $ARGS


resume_job "eqtl0_$JOBNAME_SUF"
exit 0