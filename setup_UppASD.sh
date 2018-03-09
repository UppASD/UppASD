#!/bin/bash

echo  "---------------UppASD-setup-script---------------"
echo  "    __  __          ___   _______    ____  ___   "
echo  "   / / / /__  ___  / _ | / __/ _ \  / __/ / _ \  "
echo  "  / /_/ / _ \/ _ \/ __ |_\ \/ // / /__ \_/ // /  "
echo  "  \____/ .__/ .__/_/ |_/___/____/ /____(_)___/   "
echo  "      /_/  /_/                                   "
echo  "-------------------------------------------------"
echo  "  This scripts sets up the UppASD build system   "
echo  " No files will be installed outside this folder  "
#read -p "       Press any key to continue... "
echo  "-------------------------------------------------"
echo  "  *Creating makefile and object dependencies    "
python source/make/generateDependencies.py > /dev/null
echo  "  *Scanning the system for available compilers   "
python source/make/suggestProfiles.py
