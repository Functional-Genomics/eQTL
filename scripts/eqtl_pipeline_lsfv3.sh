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
	$ECHO bsub $LSF_PARAMS -q $QUEUE $GROUP -n $THREADS -R "span[hosts=1]"  -M $MAX_MEM -R "select[mem>=$MEM] rusage[mem=$MEM]"  -w "ended($waitforids)"  -cwd `pwd` -o "${LOGS_FOLDER}$jobname-%J.out" -e "${LOGS_FOLDER}$jobname-%J.err" -J $jobname  $cmd2e 
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
    # uses EQTL_ARGS and targets
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
function stop_now {
    cur_status=$1
    if [ "$wave-" == "$cur_status-" ]; then
	targets=`eqtl_pipeline $EQTL_ARGS targets$wave | tail -n 1`
	submit_job_get_email  eqtl8_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS $targets
	resume_job "eqtl0_$JOBNAME_SUF"
	exit
    fi
}

################################################################
if [ "$*-" == "- " ]; then
    echo "ERROR: missing arguments"
    echo "usage: eqtl_pipeline_lsf_wave.sh wave conf=file [extra options]"
    exit 1
fi

wave=$1
shift 1

eqtl_pipeline $* -n -q

if [ $? -eq 0 ]; then
    echo "All done - no need to submit jobs"
    echo 0
fi
echo "info: env variables passed to LSF - LSF_PARAMS QUEUE THREADS MEM LSF_GROUP"

EQTL_ARGS="$* lsf_mode=1"

RAND=`perl -e "print int(rand()*10);"`
DATE=`date "+%w%H%M%S"`
JOBNAME_SUF="$DATE$RAND"

LOGS_FOLDER=logs/$JOBNAME_SUF/
mkdir -p $LOGS_FOLDER
# submit the jobs
echo "step0 jobs..."
submit_job "step0_$JOBNAME_SUF" ""  eqtl_pipeline $EQTL_ARGS step0
stop_job "step0_$JOBNAME_SUF"

PREV=step0_$JOBNAME_SUF
echo "set of jobs 0..."
targets=`eqtl_pipeline $EQTL_ARGS targets0 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping submission of jobs...
else
    submit_jobs eqtl0_$JOBNAME_SUF  "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl0_$JOBNAME_SUF*"
fi

echo "set of jobs 1..."
targets=`eqtl_pipeline $EQTL_ARGS targets1 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping submission of jobs...
else
    submit_jobs eqtl1_$JOBNAME_SUF  "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl1_$JOBNAME_SUF*"
fi

stop_now 1

echo "set of jobs 2..."
targets=`eqtl_pipeline $EQTL_ARGS targets2 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl2_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl2_$JOBNAME_SUF*"
fi

stop_now 2

echo "set of jobs 3..."
targets=`eqtl_pipeline $EQTL_ARGS targets3 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl3_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl3_$JOBNAME_SUF*"
fi

stop_now 3

echo "set of jobs 4..."
targets=`eqtl_pipeline $EQTL_ARGS targets4 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl4_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl4_$JOBNAME_SUF*"
fi

stop_now 4

echo "set of jobs 5..."
targets=`eqtl_pipeline $EQTL_ARGS targets5 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
if [ $? -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl5_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl5_$JOBNAME_SUF*"
fi

stop_now 5

echo "set of jobs 6..."
targets=`eqtl_pipeline $EQTL_ARGS targets6 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl6_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl6_$JOBNAME_SUF*"
fi

stop_now 6

echo "set of jobs 7..."
targets=`eqtl_pipeline $EQTL_ARGS targets7 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl7_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl7_$JOBNAME_SUF*"
fi

stop_now 7

echo "set of jobs 8..."
targets=`eqtl_pipeline $EQTL_ARGS targets8 | tail -n 1`
eqtl_pipeline $EQTL_ARGS $targets -n -q
status=$?
if [ ${#targets} -lt 2 ] || [ $status -eq 0 ]; then
    echo skipping  wave of jobs...
else
    submit_jobs eqtl8_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS
    PREV="eqtl8_$JOBNAME_SUF*"
fi

stop_now 8

echo "set of jobs 9..."
targets=`eqtl_pipeline $EQTL_ARGS targets9 | tail -n 1`
submit_jobs eqtl9_$JOBNAME_SUF "$PREV" eqtl_pipeline $EQTL_ARGS


    
# final job
targets=step4
submit_job_get_email  eqtl8_$JOBNAME_SUF "eqtl9_$JOBNAME_SUF*" eqtl_pipeline $EQTL_ARGS

resume_job "eqtl0_$JOBNAME_SUF"
exit 0
