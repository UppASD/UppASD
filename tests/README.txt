The tests are normally run with ctest from the build directory:
ctest                      (all enabled tests)
ctest -L asd-tests         (one group: asd-tests, regression-test, sld-tests, gneb-tests, cuda-tests, hip-tests)
ctest -R asd-tests.bccFe   (a single test directory)
ctest -j 8                 (in parallel)
ctest runs the tests in a copy of this directory, <build>/tests.

For help
./bergtest.py --help

For a standard simulation run as
./bergtest.py --binary /home/user/.../build/bin/sd.gfortran
The binary can also be given with the SD_BINARY environment variable.

Before the first case in a directory is run, old output files in that directory are removed.
If one still wished to append out files from previous simulation, use the dirty flag:
./bergtest.py --dirty --binary /home/user/.../build/bin/sd.gfortran

To only remove the output files in all directories, use:
./bergtest.py --cleanup

To run a specific .yaml or other test file use:
--file NAME_OF_TEST_FILE_IN_SAME_FOLDER

Running specific tests in the .yaml or other test file use:
--case # [# ...]
where # are integers, starting from 1, signifying which tests are to be run.

bergtest.py exits with an error if a test file is missing or if no tests were performed.


--------------------------------------
nondeallocs.sh
Used to analyze a meminfo file.
Does not cover all possible leaks.
