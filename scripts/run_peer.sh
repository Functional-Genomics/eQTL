#!/bin/bash
set -e
OUT=`mktemp`
runpeer.py $* | tee $OUT
STATUS=$?
if [ $STATUS -eq 0 ]; then
    # grep $OUT for convergence
    N=`grep -i -c "Converged (.* after" $OUT`
    rm -f $OUT
    if [ "$N-" == "0-" ]; then
	exit 1
    fi
else
    rm -f $OUT
    exit $STATUS
fi

exit 0
