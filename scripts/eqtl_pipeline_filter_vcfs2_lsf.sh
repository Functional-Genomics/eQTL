#!/bin/bash
################################################################
# simple wrapper to run a single instance of eqtl_pipeline in a
# cluster with the lsf job scheduler
#DEBUG=1

if [ "$MEM-" =  "-" ]; then
    MEM=12000
fi

if [ "$THREADS-" = "-" ]; then
    THREADS=1
fi

if [ "$QUEUE-" = "-" ]; then
    echo "QUEUE not defined. please run: export QUEUE=<lsf_queue_to_use>" 
    exit 1
fi


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
    # TODO: define groups by stage: 
    local GROUP=
    if [ "$LSF_GROUP-" != "-" ]; then
	GROUP="-g $LSF_GROUP"
    fi
    #########################################################
    #-R  "span[ptile=$THREADS]"
    local MAX_MEM=$MEM
    if [ "$waitforids-" != "-" ]; then
	$ECHO bsub $LSF_PARAMS -q $QUEUE -n $THREADS -R "span[hosts=1]"  -M $MAX_MEM -R "select[mem>=$MEM] rusage[mem=$MEM]"  -w "ended($waitforids)"  -cwd `pwd` -o "${LOGS_FOLDER}$jobname-%J.out" -e "${LOGS_FOLDER}$jobname-%J.err" -J $jobname  $cmd2e 
    else
	$ECHO bsub $LSF_PARAMS -q $QUEUE  $GROUP -n $THREADS -R "span[hosts=1]"  -M $MAX_MEM -R "select[mem>=$MEM]  rusage[mem=$MEM]"    -cwd `pwd` -o "${LOGS_FOLDER}$jobname-%J.out" -e "${LOGS_FOLDER}$jobname-%J.err" -J $jobname  $cmd2e 
    fi
}

function submit_job_get_email {
    local jobname=$1
    local waitforids=$2
    shift
    shift
    local cmd2e=$*
    local ECHO=
    if [ "$DEBUG-" = "1-" ]; then
        ECHO=echo
    fi
    local GROUP=
    if [ "$LSF_GROUP-" != "-" ]; then
	GROUP="-g $LSF_GROUP"
    fi
    $ECHO bsub $LSF_PARAMS -q $QUEUE  $GROUP -n $THREADS -R "span[hosts=1]"  -M $MEM -R "select[mem>=$MEM]  rusage[mem=$MEM]" -w "ended($waitforids)"  -cwd `pwd` -J $jobname  $cmd2e
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
echo "info: env variables passed to LSF - LSF_PARAMS QUEUE THREADS MEM LSF_GROUP"

ARGS=$*

RAND=`perl -e "print int(rand()*10);"`
DATE=`date "+%w%H%M%S"`
JOBNAME_SUF="$DATE$RAND"

LOGS_FOLDER=logs/$JOBNAME_SUF/
mkdir -p $LOGS_FOLDER

echo "set of jobs 2 job sthat should be completed..."
targets=`eqtl_pipeline $* targets2 | tail -n 1`
for t in $targets; do
    eqtl_pipeline $ARGS $t -n -q
    if [ $? != 0 ]; then
	echo "ERROR: $t is still not done!"
	exit 1
    fi
done
echo "set of jobs 1 completed."
# submit the jobs
echo "step0 jobs..."
submit_job "eqtl0_$JOBNAME_SUF" ""  eqtl_pipeline $ARGS step0
stop_job "eqtl0_$JOBNAME_SUF"

echo "set of jobs 3..."
targets=`eqtl_pipeline $* targets3 | tail -n 1`
submit_jobs eqtl3_$JOBNAME_SUF "eqtl0_$JOBNAME_SUF*" eqtl_pipeline $ARGS


# final job
targets=filter_vcfs
submit_job_get_email  eqtl3f_$JOBNAME_SUF "eqtl3_$JOBNAME_SUF*" eqtl_pipeline $ARGS

resume_job "eqtl0_$JOBNAME_SUF"
exit 0
