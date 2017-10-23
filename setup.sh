#!/bin/bash

# clean up previously set env
if [[ -z $FORCE_SUPERA_DIR ]]; then
    where="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    SUPERA_DIR=${where}
else
    SUPERA_DIR=$FORCE_SUPERA_DIR
fi

experiment=$1
if [ ! -d "$SUPERA_DIR/experiments/$experiment" ]; then
    echo "Could not locate $experiment directory... make sure that exists under experiments dir"
else
    export SUPERADIR=$SUPERA_DIR;
    ln -sf $SUPERA_DIR/experiments/$experiment/FMWKInterface.* $SUPERA_DIR;
    ln -sf $SUPERA_DIR/experiments/$experiment/CMakeLists.txt $SUPERA_DIR;
fi


