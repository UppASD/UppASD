# Script to run all the examples
#!/bin/bash

for n in \
ChiralSpiral \
Kagome_ncAMS \
Triangular_ncAMS \
bccFe
do
    cd $n
    echo Starts example $n
    time -p ${SD_BINARY} > out.log
    echo Ends example $n
    cd ../
done
