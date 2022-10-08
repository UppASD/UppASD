# Script to run all the examples
#!/bin/bash

cd SingleSpin
./SWEEP_DAMP.sh*
./SWEEP_FIELD.sh*
./SWEEP_TEMP.sh*
cd ../

for n in \
HeisChain \
HeisChainAF \
HeisChainDM \
bccFe \
fcc001 \
fccCo
do
    cd $n
    echo Starts example $n
    time ../../../source/sd > out.log
    echo Ends example $n
    cd ../
done
