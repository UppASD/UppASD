# Script to run all the examples
#!/bin/bash

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

for n in \
Triang2DLDN20T001
do
    cd $n
    echo Starts example $n
    ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
