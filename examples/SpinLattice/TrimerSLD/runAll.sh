# Script to run all the examples
#!/bin/bash

for n in \
TrimerauSLD2site \
TrimerauSLD2siteSym \
TrimerauSLD3site \
TrimerauSLD3siteSym \
TrimerauSLD3siteTilted
do
    cd $n
    echo Starts example $n
    ../../../../source/sd > out.log
    echo Ends example $n
    cd ../
done
