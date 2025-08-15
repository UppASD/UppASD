#!/usr/bin/env python
##############UppASD#configure#tool##################
# Written by Jonathan Werpers
# Modified, debugged and extended by Thomas Nystrand
# Extended further by Anders Bergman
#####################################################

from __future__ import print_function
from os import listdir, walk,system,environ, path, sep
#import os.path as path
import subprocess as sub
from glob import glob
import re
import shutil
import sys

# TODO these should be found and exluded by search
# Dangerous to use for both fort and C versions
# Hard to detect name conflicts can appear
external_modules = [
    'omp_lib',
    'iso_fortran_env',
]

external_CCCC = [
    '<cstdio>',
    '<cstring>',
    '<cstlib>',
    '<cuda.h>',
    '<vector>',
    '<sys/time.h>',
    '<sys/types.h>',
    '<string>',
    '<string.h>',
    '<stdio.h>',
    '<stdlib.h>',
    '<time.h>',
    '<stddef.h>',
    '<real_type.h>',
    '<sched.h>',
    '<pthread.h>',
    '<cmath>',
    '<queue>',
    '<cstdlib>',
    '<cuda.h>',
    '<curand.h>',
    '<cstddef>',
    '<cuda_runtime.h>',
    '<nvToolsExtCuda.h>',
    '<signal.h>',
    '<csignal>',
    '<list>',
    '<ctime>',
    '<iostream>',
    '<detail/Filter_File_Handle.hpp>',
    '<detail/Helpers.hpp>',
    '<fmt/format.h>',
]

third_party_CCCC=['"ovf.h"']
third_party_FORT=['ovf']

scriptpath     = path.realpath(__file__)
scriptlocation = path.dirname(scriptpath)
source_folder  = path.dirname(scriptlocation)
gpu_folder     = path.join(source_folder,'gpu_files')
third_party_folder = path.join(source_folder,'Third_party')

################################################################################
# Start of Program                                                             #
################################################################################
def main():
    print ("-----------UppASD-dependencies-tool--------------")
    print ("    __  __          ___   _______    ____  ___   ")
    print ("   / / / /__  ___  / _ | / __/ _ \  / __/ / _ \  ")
    print ("  / /_/ / _ \/ _ \/ __ |_\ \/ // / /__ \_/ // /  ")
    print ("  \____/ .__/ .__/_/ |_/___/____/ /____(_)___/   ")
    print ("      /_/  /_/                                   ")
    print ("------------------------------------------------")
    createMakefileBackup()              # Makefiles are defined in templates which overwrites the old.
    createProgramObjsDepsFile()         # List all objs. Needed for linking.
    createDependencyFile()              # Gnerate all rules needed for Makefile.
    createDependencyDotFile()           # Output files for graphical overview of dependencies.
    print ("------------------------------------------------")
    return
################################################################################
# Print dependency file which shows all rules that is needed and included in Makefile
################################################################################
def createDependencyFile():
    sourcefilesFORT = findFORTSourcefiles(source_folder)
    sourcefilesCCCC = findCUDAandCPPSourcefiles(source_folder)
    ThirdPartysourcefilesFORT =findThirdPartyFORTSourcefiles(third_party_folder,source_folder)
    ThirdPartysourcefilesCCCC =findThirdPartyCPPSourcefiles(third_party_folder,source_folder)
    sourcefilesFORT=ThirdPartysourcefilesFORT+sourcefilesFORT
    ############################################################################
    # Create the dependency file
    ############################################################################
    dependenciesFORT = readDependenciesFORT(source_folder, sourcefilesFORT,"use ")
    dependenciesCCCC = readDependenciesCCCC(source_folder, sourcefilesCCCC)
    #dependenciesCCCC_ThirdParty = readDependenciesCCCC(source_folder, ThirdPartysourcefilesCCCC)
    f = path.join(scriptlocation, "dependencies.make")
    writeDependencyFile(f, dependenciesFORT)

    fg = path.join(scriptlocation, "dependencies_c.make")
    writeDependencyFileCCCC(fg, dependenciesCCCC)

    #fg = path.join(scriptlocation, "dependencies_tp_c.make")
    #writeDependencyFileCCCC(fg, dependenciesCCCC_ThirdParty)

    #print "------------------------------------------------"

    #if len(external_modules) > 0:
    #    print "The following modules were removed from dependencies \n" + \
    #          "because they are specified as external:"
    #    for mod in external_modules:
    #        print mod

    #    for mod in external_CCCC:
    #        print mod
    print("------------------------------------------------")
    print("Module dependencies written to dependencies.make")

################################################################################
# Backup old makefiles. Copy templates to makefiles.
# The backup is used as a protection against accidental changes and commits.
# TODO: Do this for all .template endings in a loop
################################################################################
def createMakefileBackup():
    createBackup(scriptlocation,source_folder,"makefile.template","Makefile")
    createBackup(scriptlocation,scriptlocation,"makefileCUDA.template","makefileCUDA")

def createBackup(scriptlocation,source_folder,scriptname,destname):
    template_path = path.join(scriptlocation,scriptname)
    dest = path.join(source_folder,destname)
    if path.isfile(dest):
        shutil.copy(dest,dest+".backup")

    shutil.copy(template_path,dest)

################################################################################
# Creates a file which contains all obj files that Makefile will need before linking
# TODO: gpu_folder should be located automatically? All subfolders for instance.
################################################################################
def createProgramObjsDepsFile():
    mainFile = findMainFile(source_folder)

    sourcefilesFORT = findFORTSourcefiles(source_folder)
    #sourcefilesCCCC = findCUDAandCPPSourcefiles(gpu_folder)
    sourcefilesCCCC = findCUDAandCPPSourcefiles(source_folder)
    # Third party software dependencies
    ThirdPartysourcefilesFORT =findThirdPartyFORTSourcefiles(third_party_folder,source_folder)
    ThirdPartysourcefilesCCCC =findThirdPartyCPPSourcefiles(third_party_folder,source_folder)
    sourcefilesFORT=ThirdPartysourcefilesFORT+sourcefilesFORT
    #sourcefilesCCCC=ThirdPartysourcefilesCCCC+sourcefilesCCCC
    dependenciesFORT = readDependenciesFORT(source_folder, sourcefilesFORT, "use ")
    implicitFORT = calculateImplicitDependencies(dependenciesFORT, mainFile)

    f = path.join(scriptlocation, "objs.make")
    ofp = open(f, "w")
    ofp.write("OBJS = \\\n")

    fg = path.join(scriptlocation, "objs_c.make")
    ofgp = open(fg, "w")
    ofgp.write("COBJS = \\\n")

    fp = path.join(scriptlocation, "objs_tp_c.make")
    ofgptp = open(fp, "w")
    ofgptp.write("COBJSTP = \\\n")

    #maxLen = maxStringLengthInSet(implicitFORT.union(sourcefilesCCCC)) + 2
    maxLen = maxStringLengthInSet(implicitFORT) + 2
    maxLen_g = maxStringLengthInSet(sourcefilesCCCC) + 2
    maxLen_tp = maxStringLengthInSet(ThirdPartysourcefilesCCCC) + 2

    for dep in sorted(implicitFORT):
        s = "       " + (dep+".o").ljust(maxLen)+"  \\\n"
        ofp.write(s)

    for src in sorted(sourcefilesCCCC):
        s = "       " + (path.splitext(src)[0]+".o").ljust(maxLen)+"  \\\n"
        ofgp.write(s)

    for src in sorted(ThirdPartysourcefilesCCCC):
        s = "       " + (path.splitext(src)[0]+".o").ljust(maxLen)+"  \\\n"
        ofgptp.write(s)

    print("Objects written to objs.make, objs_c.make and objs_tp_c.make")

# Create the file iwsed
def createDependencyDotFile():
    sourcefiles = findFORTSourcefiles(source_folder)
    dependencies = readDependenciesFORT(source_folder, sourcefiles,"use ")
    dot_file=path.join(scriptlocation,"dependencies.dot")
    writeDotGraphFile(dot_file, dependencies)
    print("Direct dependencies written to dependencies.dot")

def maxStringLengthInSet(setOfStrings):
    length = 0
    for s in setOfStrings:
        if len(s) > length:
            length = len(s)
    return length

# Make sure we only include files that are used in program
def calculateImplicitDependencies(dependencies, filename):
    implict = set()
    moduleName = path.splitext(filename)[0]
    implict.add(moduleName)
    recureseDependencies(dependencies, implict, moduleName)
    return implict

def recureseDependencies(deps, implicit, toAdd):
    for dep in deps[toAdd+".f90"]:
        if dep not in implicit:
            implicit.add(dep)
            recureseDependencies(deps, implicit, dep)

def findFORTSourcefiles(folder):
    path_list=[]
    for root, dirs, files in walk(folder):
        dirs[:] = [d for d in dirs if d not in ['Third_party']]
        for name in files:
            if name.endswith((".F90", ".f90")):
                file_path=path.join(root, name)
                path_list.append(path.relpath(file_path,folder))
    return path_list

# Returns a the set of files with .c .cpp and .cu ending in folder
# TODO merge into other f90 function
def findCUDAandCPPSourcefiles(folder):
    path_list=[]
    for root, dirs, files in walk(folder):
        dirs[:] = [d for d in dirs if d not in ['Third_party']]
        for name in files:
            if name.endswith((".c", ".cpp", ".cu")):
                file_path=path.join(root, name)
                path_list.append(path.relpath(file_path,folder))
    return path_list

def findThirdPartyFORTSourcefiles(folder,main_file):
    path_list=[]
    for root, dirs, files in walk(folder):
        dirs[:] = [d for d in dirs if d not in ['test']]
        for name in files:
            if name.endswith((".F90", ".f90")):
                file_path=path.join(root, name)
                path_list.append(path.relpath(file_path,main_file))
    return path_list

def findThirdPartyCPPSourcefiles(folder,main_file):
    path_list=[]
    for root, dirs, files in walk(folder):
        dirs[:] = [d for d in dirs if d not in ['test']]
        for name in files:
            if name.endswith((".cpp", ".c")):
                file_path=path.join(root, name)
                path_list.append(path.relpath(file_path,main_file))
    return path_list

def readDependenciesFORT(folder, files, keyword):
    dependencies = {}
    for f in files:
        fil = path.join(folder, f)
        fp = open(fil)
        deps = set()
        for line in fp:
            if keyword in line and line.strip().startswith(keyword):
                package = line.split()[1].split(',')[0]
        	# deps.add(package)
		# Obtaining the Correct case version
		# E.g. file is called uppASD.f90
		# and fortran uses UpPaSd => uppASD is saved
                for recur_files in files:
                    fff = recur_files.split('/')[len(recur_files.split('/'))-1]
                    if package.lower()==fff[:-4].lower():
                        deps.add(recur_files[:-4])
                        break

        dependencies[f] = deps
    removeExternalDependencies(dependencies, external_modules)
    return dependencies

################################################################################
# Goes through each file and add to a set which other files it depends on
# via keywords such as #include A.cpp or A.hpp would add A.cpp and A.hpp
# Removes all libs which are known to be external
################################################################################
def readDependenciesCCCC(folder,files):
    dependencies = {}
    for f in files:
        fil = path.join(folder, f)
        deps = set()
        recursiveAddDepsCCCC(fil,deps)
        dependencies[f] = deps
    #removeExternalDependencies(dependencies, external_modules)
    return dependencies

################################################################################
# Traverse each line of file looking for includes
# Traverse each found include
################################################################################
def recursiveAddDepsCCCC(fil,deps):
    fp = open(fil)
    for line in fp:
        if "#include" in line and line.strip().startswith("#include"):
            package = line.split()[1]
            if not (package in external_CCCC):
                if "<" in package and ">" in package:
                   print("ERROR: C/CUDA/C++ depedencies not written")
                   print("       Add "+package+" to list of excluded libs")
                   sys.exit()
                if not (package in deps):
                    if package not in third_party_CCCC:
                        depsstring = path.join(gpu_folder,package.replace('"',''))
                        deps.add(depsstring)
                        nfile = path.join(sep,depsstring)
                        recursiveAddDepsCCCC(nfile,deps)
                    else:
                        for root, dirs, files in walk(third_party_folder):
                            dirs[:] = [d for d in dirs if d not in ['test']]
                            for name in files:
                                if name.endswith((package.replace('"',''))):
                                    file_path=path.join(root, name)
                                    depsstring = file_path
                                    deps.add(depsstring)
                                    nfile = path.join(sep,depsstring)
                                    recursiveAddDepsCCCC(nfile,deps)
    fp.close()
    return deps



# Removes all external used libraries that are listed in ext_libs
# from the current set of dependencies
def removeExternalDependencies(dependencies, ext_libs):
    for f, deps in dependencies.items():
        for lib in ext_libs:
            deps.discard(lib)


# Locates the fortran file which starts the execution
# Used as a starting point for implicit dependencies
def findMainFile(folder):
    sourcefiles = findFORTSourcefiles(folder)
    programs = []
    for f in sourcefiles:
        fil = path.join(folder, f)
        fp = open(fil)
        for line in fp:
            if "program" in line and line.strip().startswith("program"):
                programs.append(f)

    if len(programs) != 1:
        print("WARNING: Found more than one program statement.", programs)
    return programs[0]


# Write info which can be used to graphically view the
# dependencies between the different modules
def writeDotGraphFile(f, dependencies):
    fp = open(f, "w")
    fp.write('digraph G {\n')
    fp.write('  overlap=scale\n')
    fp.write('  concentrate=true\n')
    fp.write('  splines=true\n')
    fp.write('   ratio="fill"\n')
    fp.write('    size="8.3,11.7!"\n')
    fp.write('     margin=0\n')


    for f, deps in sorted(dependencies.items()):
        f = f[:-4].lower()
        for dep in sorted(deps):
            # New order
            fp.write("  \"" + f +"\" -> \"" + dep.lower() + "\" ;\n")
            # Old order
            #fp.write("  " + dep.lower()+" -> " + f + ";\n")
        fp.write("\n")
    fp.write("}\n")

# Generate the rules needed for Makefile
# All must be done in lower
# Toherwise make problems will appear
# Fortran compilers always outputs lower case mod files for some reason...
def writeDependencyFile(f, dependencies):
    ofp = open(f, "w")
    for f, deps in sorted(dependencies.items()):
        if len(deps) > 0:
            target = f[:-4]+".o " + f[:-4]+".mod"
            #target = f[:-4]+".o " + f[:-4]+".o"
            ofp.write(target+":")
            for dep in sorted(deps):
                #ofp.write(" "+dep+".mod")
                ofp.write(" "+dep+".o")
            ofp.write("\n")


# Generate the rules needed for Makefile
def writeDependencyFileCCCC(f, dependencies):
    ofp = open(f, "w")
    for f, deps in sorted(dependencies.items()):
        if len(deps) > 0:
            target = f.split('.')[0]+".o "
            ofp.write(target+":")
            for dep in sorted(deps):
                ofp.write(" "+dep)
            ofp.write("\n")



# Prints the file dependencies to screen
def printDependencies(srcFiles, dependencies):
    for f in srcFiles:
        deps = dependencies[f]
        if len(deps) > 0:
            print(f+":")
            for dep in sorted(deps):
                print("    " + dep)
        else:
            print(f + " has no dependendcies.")
        print


def removeLeadingDigits(s):
    return re.sub()

main()
