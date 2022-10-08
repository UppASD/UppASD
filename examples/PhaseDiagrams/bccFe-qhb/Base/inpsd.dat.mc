simid bccFe100                                    
ncell     10        10        10                System size            
BC        P         P         P                 Boundary conditions (0=vacuum, P=periodic)
cell      1.00000   0.00000   0.00000         
          0.00000   1.00000   0.00000
          0.00000   0.00000   1.00000
Sym       1                                     Symmetry of lattice (0 for no, 1 for cubic, 2 for 2d cubic, 3 for hexagonal)

posfile   ./posfile
momfile   ./momfile
exchange  ./jfile

SDEalgh   1                                     SDE solver: 1=midpoint, 2=heun, 3=heun3, 4=Heun_proper, 5=Depondt
Initmag   1                                     Initial config of moments (1=random, 2=cone, 3=spec., 4=file)
#restartfile ./restart.bccFe100.out
Mensemble 5

ip_mode   H                                     Initial phase parameters
ip_temp   TEMP                                  --
ip_mcNstep     10000                            --

mode      M                                     S=SD, M=MC
temp      TEMP                                  Measurement phase parameters
mcNstep   20000 

do_avrg   Y                                     Measure averages

do_cumu Y

do_tottraj N                                    Measure moments

plotenergy 1

# Quantum statistics (R is best choice)
do_qhb R
# Algorithm for rescaling (TC=from given T_curie, TM=from Mean-field T_curie, TR=from RPA T_curie, MT=from measured magnetization)
qhb_mode TC
# Given T_curie if qhb_mode=TC
tcurie 1000.0

# Calculate adiabatic magnon spectrum
do_ams Y
# q-points from full cell
qpoints C
# Read magnon DOS from file
do_magdos F
# Precalculated magnon DOS (with higher resolution)
magdosfile ./magdos.bccFe.dat
