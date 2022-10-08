#!/bin/bash
# Creates q-file for UppASD bcc (110) surface
qlen=$1
rm -f tmpfile0 tmpfile 
# neg x 
for q in `seq -$qlen 0`
do 
   echo $q | awk -v qlen=$qlen '{ printf "%f %f %f\n", 0,0,$1/qlen/2 }' >> tmpfile 
done
# pos x 
for q in `seq 1 $qlen`
do 
   echo $q | awk -v qlen=$qlen '{ printf "%f %f %f\n", 0,0,$1/qlen/2 }' >> tmpfile 
done
# number of lines
nq=`wc tmpfile | awk '{ print $1 }'` 
echo $nq | awk '{ printf "            %i\n",$1 }' > tmpfile0 
cat tmpfile0 tmpfile > qfile
rm -f tmpfile0 tmpfile 

