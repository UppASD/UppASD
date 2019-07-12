##########################################
## Makefile Settings for ifort profile  ##
##########################################

#------------------------------------------------#
# The different compilers
#------------------------------------------------#

current_dir = $(shell pwd)

INCLUDEDIR=${current_dir}/Third_party/ovf/include

#------------------------------------------------#
# Flags for FORTRAN compilation
#------------------------------------------------#
# Basic optimization settings explained
# -ip         Inline function, substantioal speed up
# -O3         Optimization, faster execution, slow make times
# -ipo        Inline between files
# -xP         Intel processor specific optimizations
# -fast       Uses -ipo -O3 -xP  -static
OVFFCFLAGS = -O3 -ip  -xHost

# Basic debug settings explained
# -g           Debug with gdb
# -CB          Array bound checks (same as check)
# -traceback   Gdb traceback
# -p           Function profiling with gprof
# -warn all    Show all warning messages
# -check all   Check all
# -xT          Optimization for intel(R) core(TM)2 Duo
#OVFFCDEBUG = -g -traceback
OVFFCDEBUG =

OVFCOMPFLAGS = -fPIC

OVFCCOMPLAGS = -DFMT_HEADER_ONLY -I $(INCLUDEDIR) -fPIC -std=c++11
