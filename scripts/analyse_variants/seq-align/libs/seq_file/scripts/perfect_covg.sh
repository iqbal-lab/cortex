#!/bin/bash

# Simulate reads starting at every base, with a given read length
# Requires dnacat (https://github.com/noporpoise/seq_file)

set -euo pipefail

if [[ $# -ne 2 || !( $1 =~ ^[0-9]+$ ) ]]
then
  echo "usage ./perfect_covg.sh <readlen> <ref.fa>"
  exit -1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DNACAT="$DIR/../bin/dnacat"

if [[ !(-x $DNACAT) ]]
then
  echo "Please compile dnacat with 'make'"
  exit -1
fi

$DNACAT -P $2 | \
  awk 'BEGIN{N='$1'} {l=length($0); for(i=1; i+N-1<=l; i++){ print substr($0,i,N) }}'

exit $?
