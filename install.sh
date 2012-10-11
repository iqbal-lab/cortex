#!/bin/bash

if [ ! -d libs ]; then
    exit 0
fi


cd libs

cd samtools-0.1.18
make clean
make
cd ..

cd gsl-1.15
make clean
./configure
make
cd ..

cd string_buffer
make clean
make
cd ..

cd seq_file
make clean
make STRING_BUF_PATH=../string_buffer SAMTOOLS_PATH=../samtools-0.1.18
cd ..

cd ..