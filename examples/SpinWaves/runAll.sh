# Script to run all the examples
#!/bin/bash

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

for n in \
ChiralSpiral \
Kagome_ncAMS \
Triangular_ncAMS \
bccFe
do
    cd $n
    echo Starts example $n
    time ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
