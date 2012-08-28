#!/bin/bash

EXE_FILE=main.out
cd src
make clean
make
./$EXE_FILE >/dev/null
cd ..
gcov src/Energy.cc 
lcov --directory $PWD/src --capture  --output-file stat.info
mkdir -p Stats
genhtml -o Stats/ stat.info
