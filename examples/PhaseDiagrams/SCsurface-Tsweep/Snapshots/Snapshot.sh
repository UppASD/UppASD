#!/bin/bash
iname=`ls -l inp.* | tail -1 | awk '{ print $NF }'`
if [[ "$iname" == inp.*.out.gz ]]; then
   nat=`zcat $iname | grep Natom | awk '{ print $2 }'`
else
   nat=`grep Natom $iname | awk '{ print $2 }'`
fi
sed "s/NMAX/$nat/" Snapshot.py > mySnapshot.py
zcat coord.* > atomsparsed.out
mname=`ls -l moment.* | tail -1 | awk '{ print $NF }'`
if [[ "$mname" == moment.*.out.gz ]]; then
   zcat $mname | tail -$nat > momentsparsed.out
else
   cat $mname | tail -$nat > momentsparsed.out
fi
vtkpython mySnapshot.py
rm -f momentsparsed.out mySnapshot.py
