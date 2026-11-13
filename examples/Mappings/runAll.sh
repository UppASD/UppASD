# Script to run all the examples
#!/bin/bash

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

for n in \
FeCo
do
    cd $n
    for m in \
    B2 \
    random
    do
        cd $m
        echo Starts example $n''/$m
        time ${SD_BINARY} > out.log
        echo Ends example $n $m
        cd ../
    done
    cd ../
done

for n in \
GaMnAs-onecell \
RSPtElkAFM \
RandomAlloy \
bccFe-tensor \
bccFe-variants
do
    cd $n
    echo Starts example $n
    time ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
