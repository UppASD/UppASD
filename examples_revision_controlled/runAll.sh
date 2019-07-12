# Script to run all the examples
#!/bin/bash

# GNEB_HeisChain broken
# xctensortest broken
# fccCo long execution time
for n in \
Fe \
GaMnAs \
HeisChain \
HeisChainAF \
HeisChainDM \
RSPtElkAFM \
Ralloy \
bccFe
do
    cd $n
    echo Starts example $n
    time ../../source/sd > out.log
    echo Ends example $n
    cd ../
done

for n in \
2d-systems/fcc001 \
FeCo/B2 \
FeCo/random
do
    cd $n
    echo Starts example $n
    time ../../../source/sd > out.log
    echo Ends example $n
    cd ../../
done

cd bccFe-Tsweep
./runme.sh
cd ..

cd 2d-systems/sc
echo Starts example 2d-systems/sc
./runme.sh
echo Ends example 2d-systems/sc
cd ../../

cd SingleSpin
echo Starts example SingleSpin
./SWEEP_DAMP.sh
./SWEEP_FIELD.sh
./SWEEP_TEMP.sh
echo Ends example SingleSpin
cd ..
