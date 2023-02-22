#!/usr/bin/env python3
##############UppASD#profile#tool##################
# Written by Anders Bergman
#####################################################

from __future__ import print_function
from os import listdir, walk,system,environ, path
#import os.path as path
import subprocess as sub
from glob import glob
import re
import shutil
import sys

# Start of Program #
def main():
    print("-------------UppASD-profile-helper---------------")
    print("    __  __          ___   _______    ____  ___   ")
    print("   / / / /__  ___  / _ | / __/ _ \  / __/ / _ \  ")
    print("  / /_/ / _ \/ _ \/ __ |_\ \/ // / /__ \_/ // /  ")
    print("  \____/ .__/ .__/_/ |_/___/____/ /____(_)___/   ")
    print("      /_/  /_/                                   ")
    print("------------------------------------------------")
    guessEnvironmentSettings()          # Output files for graphical overview of dependencies.
    print("------------------------------------------------")

#Probe the system for library paths
#and available compilers
#Also tries to find the level of OpenMP support
#according to http://www.openmp.org/resources/openmp-compilers/
#
def guessEnvironmentSettings():
    print("    Probing system for compilers/libraries      ")
    print("                                                ")
    p = sub.Popen(['which', 'gfortran'],stdout=sub.PIPE,stderr=sub.PIPE)
    gfortran, errors = p.communicate()
    if (len(gfortran)>0): 
        have_gfortran=True
        p = sub.Popen(['gfortran', '-v'],stdout=sub.PIPE,stderr=sub.PIPE)
        gfortranver, errors = p.communicate()
        errors=errors.decode()
        sta=errors.find('gcc version ')+12
        sto=errors[sta:].find(' ')+sta
        gv=errors[sta:sto]
        print("   GNU Fortran compiler found, version",gv)
        print("    FC = gfortran")
        gvmaj=int(gv[:gv.find('.')])
        gvmin=float(gv[gv.find('.')+1:])
        gsimd=(gvmaj>=4)and(gvmin>=9.1)
        #if(gsimd):
        #    print("    FCOMPFLAGS = -fopenmp -fopenmp-simd")
        #else:
        #    print("    FCOMPFLAGS = -fopenmp ")
    else:
        have_gfortran=False
        print(" GNU Fortran compiler not found")

    p = sub.Popen(['which', 'ftn'],stdout=sub.PIPE,stderr=sub.PIPE)
    ftn, errors = p.communicate()
    if (len(ftn)>0):
        have_ftn=True
    else:
        have_ftn=False

    p = sub.Popen(['which', 'ifort'],stdout=sub.PIPE,stderr=sub.PIPE)
    ifort, errors = p.communicate()
    if (len(ifort)>0): 
        have_ifort=True
        p = sub.Popen(['ifort', '-v'],stdout=sub.PIPE,stderr=sub.PIPE)
        ifortver, errors = p.communicate()
        errors=errors.decode()
        sta=errors.find('ifort version ')+14
        iv=errors[sta:-1]
        #print("                                  ")
        print("   Intel Fortran compiler found, version",iv)
        #print("    FC = ifort")
        #ivmaj=int(iv[:iv.find('.')])
        #ivmin=float(iv[iv.find('.')+1:])
        #if(ivmaj>=16):
        #    print("    FCOMPFLAGS = -qopenmp ")
        #elif(ivmaj==13 and ivmin==3):
        #    print("    FCOMPFLAGS = -openmp -openmp-simd")
        #else:
        #    print("    FCOMPFLAGS = -openmp")
    else:
        have_ifort=False
        #print(" Intel Fortran compiler not found")

    p = sub.Popen(['which', 'pgf90'],stdout=sub.PIPE,stderr=sub.PIPE)
    pgf90, errors = p.communicate()
    if (len(pgf90)>0): 
        have_pgf90=True
        p = sub.Popen(['pgf90', '--version'],stdout=sub.PIPE,stderr=sub.PIPE)
        pgf90ver, errors = p.communicate()
        pgf90ver=pgf90ver.decode().split(' ')
        iv=pgf90ver[1]
        #print("                                  ")
        print("   PGI Fortran compiler found, version",iv)
        #print("    FC = pgf90")
        #ivmaj=int(iv[:iv.find('.')])
        #ivmin=float(iv[iv.find('.')+1:])
        #if(ivmaj>=16):
        #    print("    FCOMPFLAGS = -qopenmp ")
        #elif(ivmaj==13 and ivmin==3):
        #    print("    FCOMPFLAGS = -openmp -openmp-simd")
        #else:
        #    print("    FCOMPFLAGS = -openmp")
    else:
        have_pgf90=False
        #print(" Intel Fortran compiler not found")

    mkl=environ.get('MKLROOT')
    if (mkl!=None): 
        have_mkl=True
        #print("                                  ")
        print("   Intel Math Kernel Library (MKL) found")
        #print("    FLIBFLAGS = -mkl=sequential")
        #print("    USE_VSL = YES ")
    else:
        have_mkl=False
        #print(" Intel Math Kernel Library (MKL) not found")
        #print("    USE_VSL = NO")

    osx=path.isfile('/System/Library/Frameworks/Accelerate.framework/Accelerate')
    if (osx):
        have_osx=True
        #print("                                  ")
        print("   Apple Accelerate framework found ")
        #print("    FLIBFLAGS = -framework=Accelerate ")
        #print("    USE_VSL = NO  ")
    else:
        have_osx=False
        #print(" Apple Accelerate framework not found ")

    got_locate=system("which locate > /dev/null 2>&1 ")
    #got_locate=system("which locate &> /dev/null")
    got_ldconf=system("which ldconfig > /dev/null 2>&1")
    have_blas=False
    if(got_locate==0):
        p = sub.Popen(['locate', 'libblas.so'],stdout=sub.PIPE,stderr=sub.PIPE)
        blas, errors = p.communicate()
        p = sub.Popen(['locate', 'liblapack.so'],stdout=sub.PIPE,stderr=sub.PIPE)
        lapack, errors = p.communicate()
        if (len(lapack)>0 and len(blas)>0): 
            have_blas=True
    elif(got_ldconf==0):
        blas=system("ldconfig -p | grep blas &> /dev/null")
        lapack=system("ldconfig -p | grep lapack &> /dev/null")
        if (blas==0 and lapack==0):
            have_blas=True

    if (have_blas==True):
        #print("                                  ")
        print("   Standard BLAS/LAPACK library found")
        #print("    FLIBFLAGS = -lblas -llapack")
    #else:
    #    #print(" Generic BLAS/LAPACK library not found")

    p = sub.Popen(['which', 'nvcc'],stdout=sub.PIPE,stderr=sub.PIPE)
    nvcc, errors = p.communicate()
    if (len(nvcc)>0): 
        have_cuda=True
        #print("                                  ")
        print("   Nvidia CUDA compiler found")
        #print("    USE_CUDA = YES")
        #print("                                  ")
        #print("    CUDA_INSTALL_PATH = ",nvcc[:-10] )
        #print("    CUDA_INCLUDE_PATH = ",nvcc[:-10]+"/include")
        #print("    CUDA_LIBRARY_PATH = ",nvcc[:-10]+"/lib64")
    else:
        have_cuda=False
        #print(" Nvidia CUDA compiler not found")
        #print("    USE_CUDA = NO ")

    print("                                  ")
    print("        Suggested compilation profiles  ")
    print("           (may need manual tuning) ")
    prof_list=[]
    if(have_ifort):
        if(have_mkl):
            prof_list.append('ifort')
            if(have_cuda):
                prof_list.append('ifort-cuda')
        else:
            prof_list.append('ifort-nomkl')
            if(have_cuda):
                prof_list.append('ifort-cuda-nomkl')
    if(have_pgf90):
        if(have_mkl):
            prof_list.append('pgf90')
        else:
            prof_list.append('pgf90-nomkl')
    if(have_gfortran):
        if(have_osx):
            prof_list.append('gfortran-osx')
        elif(have_blas):
            prof_list.append('gfortran')
            if(have_cuda):
                prof_list.append('gfortran-cuda')
    if(have_ftn):
        prof_list.append('crayftn-ftn')
        prof_list.append('gfortran-ftn')
        prof_list.append('aocc-ftn')
        print("   Cray Fortran compiler wrapper found")
        print("    FC = ftn")
    print('  '+'; '.join('{}'.format(k) for k in prof_list))

    if (have_cuda):
        nvcc=nvcc.decode()
        print("                                              ")
        print("  Suggested CUDA paths: (please edit profile)" )
        print("    CUDA_INSTALL_PATH = ",nvcc[:-10] )
        print("    CUDA_INCLUDE_PATH = ",nvcc[:-10]+"/include")
        print("    CUDA_LIBRARY_PATH = ",nvcc[:-10]+"/lib64")

    print("------------------------------------------------")
    if (have_ifort or have_gfortran):
        print("     Compile the code by: make <profile>        ")
        print("     or: make PROFILE=<profile> PROG=target      ")
    else:
        print("     No compiler was auto-detected. ")
        print("     Please modify an own profile. ")




main()
