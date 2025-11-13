#! /bin/bash

# Perform a loop over various values of the field
for Field in 0 20 25 30 35 40 45
do

   damping=$(echo "scale=4; 100/1000" | bc | awk '{printf "%4.4f\n",$0}')
   bfield=$(echo "scale=4; $Field/100" | bc | awk '{printf "%4.2f\n",$0}')

   # Create a folder for the calculation
   mkdir FIELD$bfield
 
   # Copy input from the Base folder to the calculation folder
   cp BASE/* FIELD$bfield

   # Enter the calculation folder
   cd FIELD$bfield
 
      # Modify the files for the actual calculation
      # Setting the damping
      sed -i "s/GILBERT/$damping/g" inpsd.dat

      # The temperature is set to zero for simplicity
      sed -i "s/TEMPE/0.0000/g" inpsd.dat

      # Add the external magnetic field
      echo -e "hfield 0.00  0.00  -$bfield" >> inpsd.dat

      echo -e "Performing simulation for bfield $bfield"
      # Running the ASD simulation
      ${SD_BINARY} > out.out
      echo -e "done"

   cd ../

done

exit
