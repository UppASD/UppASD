# Script to run all the examples
#!/bin/bash

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

for n in \
GNEB_HeisChain \
SkyrmionLattice \
SkyrmionTransport
do
    cd $n
    echo Starts example $n
    time ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
