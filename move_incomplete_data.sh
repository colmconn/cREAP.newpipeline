#!/bin/bash

## set -x 

programName=`basename $0`

trap exit SIGHUP SIGINT SIGTERM

studyName=cREAP.newpipeline

GETOPT=$( which getopt )
ROOT=/data/sanDiego/$studyName
DATA=$ROOT/data
## this variable is used to search for a preexisting non-linear warped
## anatomy to avoid having to compute the warp again since it's
## so time consuming to do so
MACHLEARN_DATA=/data/sanDiego/cEMM.machlearn/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

subjects=$( cd $DATA  ; ls -d [0-9][0-9][0-9]_{A,B,C,C,D,E} [0-9][0-9][0-9]_{A,B,C,C,D,E}2 2> /dev/null | grep -v 999 | sort -u )

for subject in $subjects ; do
    regressorFile=$DATA/regressors/${subject}_reap.wav.1D
    
    if  [[ ! -f ${DATA}/$subject/${subject}reap+orig.HEAD ]] || \
	    [[ ! -f ${DATA}/$subject/${subject}+orig.HEAD ]] || \
	    [[ ! -f ${regressorFile} ]]	; then 

	mv -f $DATA/${subject} ${DATA}/incomplete
    fi
    
done

rm -fr $DATA/incomplete
