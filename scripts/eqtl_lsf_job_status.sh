#!/bin/bash
#######################################################################
#
# Copyright 2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
jobname=$1
JOB_ID=$2
MEM=$3
LOG_DIR=$4

shift 4
if [  $? -gt 0 ]; then
    echo "ERROR! Usage:  jobname jobid MEM log_dir job_cmd" >&2
    exit 1
fi
JOB_CMD=$*

if [ "$JOB_CMD-" == "-" ]; then
    echo "ERROR! Usage:  jobname jobid MEM log_dir job_cmd" >&2
    exit 1
fi

echo "CMD=$JOB_CMD" 
# sleep for a while otherwise you get the dependent job as being in the RUN state
sleep 60
echo "Checking jobs status...."
STATUS=`bjobs  -a -J $jobname| cut -c 17-22|tail -n 1| sed -E "s/([a-zA-Z]*).*/\1/"`
echo "Checking jobs status....done."
if [ "$STATUS-" == "-" ]; then
    echo "ERROR: Job $job information not found" > /dev/stderr
    exit 1
fi

# print some info about the job
bjobs  -a -l -J $jobname

if [ "$STATUS-" == "DONE-" ]; then
 EXIT_CODE=0
 echo "Job finished successfully ($STATUS)"
 exit $EXIT_CODE
else
 EXIT_CODE=1
 echo "Job failed ($STATUS)"
 EXIT_CODE=$?
fi

INSUF_MEM=0
# Provide a bit more information
function check_jobs_status {
  dir=$1
  let success_jobs=0
  let failed_jobs=0

  function check_jobs_in_dir {    
    local f
    local out_files
    out_files=`ls -1 $1/*.out 2>/dev/null`
    for f in $out_files; do
       OK=`grep -c -E  "^Successfully completed." $f`
       #echo "OK=$OK" > /dev/stderr
       if [ "$OK" == "1" ]; then
         let success_jobs=$success_jobs+1
       else 
         # error
	 err_f=`echo $f|sed "s/.out$/.err/"`
	 if [ -e $err_f ]; then
	     ERR_LAST_LINES=`tail -n 10 $err_f`
	 else
	     ERR_LAST_LINES="$err_f not found"
	 fi
         LAST_LINES=`tail -n 50 $f|head -n 9`
         let failed_jobs=$failed_jobs+1
         EXIT_STATUS=`grep "Exited with exit code" $f|tail -n 1 |cut -f 5 -d\ |sed "s/\.$//"`
	 MEM_INSUF=`grep -c "^TERM_MEMLIMIT: job killed after reaching LSF memory usage limit." $f`
         CMD=`grep "LSBATCH: User input" -A 1 $f|tail -n 1`
	 echo "------------------------------------------" 
	 echo "------------------------------------------" 
	 echo "Output file: $f" 
	 echo "Error file: $err_f" 
         echo "Exit status: $EXIT_STATUS" 
	 if [ "$MEM_INSUF" == "1" ]; then 
	     INSUF_MEM=1
	     echo "Insufficient memory"
	 fi
	 echo "$CMD" 
	 echo "*********** $f"
	 echo "$LAST_LINES" 
	 echo "*********** $err_f"
	 echo "$ERR_LAST_LINES"
    fi
    done
    # checking subfolders
    local DIRS=`ls -d $1/*/ 2> /dev/null`
    local d
    for d in $DIRS; do
      check_jobs_in_dir $d
    done 
   }
   check_jobs_in_dir $dir
   
   echo JOBS OK=$success_jobs
   if [ $failed_jobs \> 0 ]; then
       echo JOBS FAILED=$failed_jobs
       if [ $INSUF_MEM == 1 ]; then 
	   EXIT_STATUS=1
       else
	   EXIT_STATUS=2
       fi
   else
       EXIT_STATUS=0
   fi
   
}
EXIT_STATUS=

check_jobs_status $LOG_DIR

exit $EXIT_STATUS
