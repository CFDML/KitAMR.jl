#!/bin/bash
git clone https://github.com/cburstedde/p4est.git
cd ./p4est
git checkout -b develop
git submodule init
git submodule update
