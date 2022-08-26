#!/bin/bash

# clean up previously set env
if [[ -z $FORCE_SUPERA_DIR ]]; then
    where="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    SUPERA_DIR=${where}
else
    SUPERA_DIR=$FORCE_SUPERA_DIR
fi

experiment=$1
if [ -z $experiment ]; then
    echo "You must specify an experiment in the argument of setup.sh"
    echo "Options: `ls $SUPERA_DIR/experiments`"
    unset experiment;
    unset SUPERA_DIR;
else
    if [ ! -d "$SUPERA_DIR/experiments/$experiment" ]; then
	echo "Could not locate $experiment directory... make sure that exists under experiments dir"
    else
	export SUPERADIR=$SUPERA_DIR;
	ln -sf $SUPERA_DIR/experiments/$experiment/FMWKInterface.* $SUPERA_DIR;
	ln -sf $SUPERA_DIR/experiments/$experiment/CMakeLists.txt $SUPERA_DIR;
	ln -sf $SUPERA_DIR/experiments/$experiment/ExperimentTypes.h $SUPERA_DIR;
	echo "set up for $experiment"
    fi
fi


