#!/bin/bash

# Now only checks for uneven number of entries
# Script should also take in account if each entry is an alloc or dealloc
# to protect against for instance two allocs and zero deallocs

if [[ $# -eq  0 ]] ; then
  echo 'First argument must specify directory to meminfo file'
  exit 0
fi

TFILE=$(mktemp)
awk 'NR>1 {print $2}' $1 > $TFILE
echo 'Recorded allocs that misses a recorded dealloc:'
while read allocated; do
  if [ "$(($(grep -wc $allocated $TFILE) % 2))" -eq "1" ] ; then
    echo $allocated
  fi
done < $TFILE
rm $TFILE
