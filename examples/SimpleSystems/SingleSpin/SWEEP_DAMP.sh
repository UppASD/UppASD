#! /bin/bash

# Perform a loop over various values of the damping
for damp in 0 1 10 100 1000 10000 100000
do

   # Transform the damping to acceptable values for ASD
   damping=$(echo "scale=4; $damp/1000" | bc | awk '{printf "%4.4f\n",$0}')

   # Create a folder for the calculation
   mkdir DAMP$damping
 
   # Copy input from the Base folder to the calculation folder
   cp BASE/* DAMP$damping

   # Enter the calculation folder
   cd DAMP$damping
 
      # Modify the files for the actual calculation
      # Setting the damping
      sed -i "s/GILBERT/$damping/g" inpsd.dat

      # The temperature is set to zero for simplicity
      sed -i "s/TEMPE/0.0000/g" inpsd.dat

      echo -e "Performing simulation for damping $damping"
      # Running the ASD simulation
      ../../../../source/sd > out.out
      echo -e "done"

   cd ../

done

exit
