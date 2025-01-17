#!/bin/bash
export LD_LIBRARY_PATH=$LIBP4EST_PATH:$LD_LIBRARY_PATH
mpiexecjl --project=. -n $1 julia $2