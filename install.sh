#!/bin/bash

if [ ! -d libs ]
then
  echo "Directory libs does not exist" 1>&2
  exit 1
fi

make -C libs/htslib clean
make -C libs/htslib

make -C libs/string_buffer clean
make -C libs/string_buffer

make -C libs/seq_file clean
make -C libs/seq_file STRING_BUF_PATH=../string_buffer HTS_PATH=../htslib

make -C scripts/analyse_variants/seq-align clean
make -C scripts/analyse_variants/seq-align
