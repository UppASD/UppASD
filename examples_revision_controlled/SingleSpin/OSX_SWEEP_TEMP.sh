#! /bin/bash

# Perform a loop over various values of the temperature
for TEMP in 1 2 3 4 5 10
do

   damping=$(echo "scale=4; 100/1000" | bc | awk '{printf "%4.4f\n",$0}')

   # Create a folder for the calculation
   mkdir T$TEMP
 
   # Copy input from the Base folder to the calculation folder
   cp BASE/* T$TEMP

   # Enter the calculation folder
   cd T$TEMP
 
      # Modify the files for the actual calculation
      # Setting the damping
      sed -i '' "s/GILBERT/$damping/g" inpsd.dat

      # Setting the temperature
      sed -i '' "s/TEMPE/$TEMP/g" inpsd.dat
      # Set the spin along the z-axis
      sed -i '' "s/roteul    1/roteul    0/g" inpsd.dat

      # Changing the number of replicas or ensembles that are performed
      for ENS in 200
      do

         mkdir EN$ENS

         cp inpsd.dat jfile momfile posfile kfile EN$ENS

         cd EN$ENS
  
            sed -i '' "s/Mensemble 1/Mensemble $ENS/g" inpsd.dat
            echo -e "Performing simulation for temperature $TEMP and $ENS emsembles"
            # Running the ASD simulation
            ../../../../source/sd > out.out
            echo -e "done"

         cd ../
      done
     
     # For T=4 perform the simulation for a single ensemble
     if [ $TEMP -eq 4 ] 
     then
 
        mkdir EN1

        cp inpsd.dat jfile momfile kfile posfile EN1

        cd EN1


            sed -i '' "s/Nstep     100000/Nstep     1000000/g" inpsd.dat
            sed -i '' "s/timestep  1.000e-14/timestep  1.000e-16/g" inpsd.dat
            sed -i '' "s/Mensemble 1/Mensemble 1/g" inpsd.dat
            echo -e "Performing simulation for temperature $TEMP and 1 emsembles"
            # Running the ASD simulation
            ../../../../source/sd > out.out
            echo -e "done"


        cd ../   
 
     fi

   cd ../

done

exit
