#!/bin/bash

CURRENT_DIR=$PWD
export LD_LIBRARY_PATH=$CURRENT_DIR/libs/gsl-1.15/:$CURRENT_DIR/libs/gsl-1.15/blas/:$CURRENT_DIR/libs/gsl-1.15/.libs/:$CURRENT_DIR/libs/gsl-1.15/blas/.libs/:$CURRENT_DIR/libs/gsl-1.15/cblas/.libs/:$LD_LIBRARY_PATH

