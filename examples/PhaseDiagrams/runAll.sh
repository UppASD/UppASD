# Script to run all the examples
#!/bin/bash

for n in \
Ising-Tsweep \
SCsurface-Tsweep \
bccFe-Tsweep \
bccFe-qhb
do
    cd $n
    echo Starts example $n
    time -p ./runme.sh > out.log
    echo Ends example $n
    cd ../
done
