# Script to run all the examples
#!/bin/bash

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

for n in \
TetramerpuSLD \
TetramerpuSLDAijk \
TetramerpuSLDAijkfp \
TetramerpuSLDfp/
do
    cd $n
    echo Starts example $n
    ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
