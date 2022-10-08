# Script to run all the examples
#!/bin/bash

for n in \
bccFe3TM \
bccFeN10LDT300a000 \
bccFeN10SLDT300a000
do
    cd $n
    echo Starts example $n
    ../../../../source/sd > out.log
    echo Ends example $n
    cd ../
done
