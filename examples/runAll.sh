# Script to run all the examples
#!/bin/bash

for n in \
Mappings \
PhaseDiagrams \
SimpleSystems \
SpecialFeatures \
SpinLattice \
SpinWaves
do
    cd $n
    echo Starts example $n
    ./runAll.sh > out.log
    echo Ends example $n
    cd ../
done
