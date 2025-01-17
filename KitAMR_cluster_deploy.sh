#!/bin/bash
./bootstrap
./configure CC=mpicc --enable-mpi --enable-shared --prefix=../libp4est
make
lib_path="$pwd/libp4est/lib/"
BASHRC_FILE="$HOME/.bashrc"
if grep -q "^export LIBP4EST_PATH=" "$BASHRC_FILE"; then
    sed -i "s#^export LIBP4EST_PATH=.*#export LIBP4EST_PATH=$lib_path#" "$BASHRC_FILE"
else
    echo "export LIBP4EST_PATH=$lib_path" >> "$BASHRC_FILE"
fi
source $HOME/.bashrc
