#!/bin/bash

# Deletes all .out files.
# Add more cleanup fileformates as needed
find . -name '*.out' -type f -delete
find . -name 'out*' -type f -delete
find . -name 'meminfo' -type f -delete
find . -name 'snap*.p??' -type f -delete
find . -name 'fort.*' -type f -delete
find . -name 'averages.*' -type f -delete
find . -name 'cumulants.*' -type f -delete
find . -name 'totenergy.*' -type f -delete
