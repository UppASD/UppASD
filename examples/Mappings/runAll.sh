# Script to run all the examples
#!/bin/bash

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
        time ../../../../source/sd > out.log
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
    time ../../../source/sd > out.log
    echo Ends example $n
    cd ../
done
