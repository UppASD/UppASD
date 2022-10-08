# Script to run all the examples
#!/bin/bash

for n in \
GNEB_HeisChain \
SkyrmionLattice \
SkyrmionTransport
do
    cd $n
    echo Starts example $n
    time ../../../source/sd > out.log
    echo Ends example $n
    cd ../
done
