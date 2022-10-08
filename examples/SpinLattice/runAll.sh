# Script to run all the examples
#!/bin/bash

for n in \
DimerLD \
DimerSLD \
TetramerSLD \
Triang2DSLD \
TrimerSLD \
bccFeSLD
do
    cd $n
    echo Starts example $n
    ./runAll.sh > out.log
    echo Ends example $n
    cd ../
done
