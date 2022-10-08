For help
./bergtest.py --help

For a standard simulation run as
./bergtest.py --cleanup
./bergtest.py --binary /home/user/.../source/sd

The first cleanup step is not necessary since cleaning is done by default.
If one still wished to append out files from previous simulation, use the dirty flag:
./bergtest.py --dirty --binary /home/user/.../source/sd

To run a specific .yaml or other test file use:
--file NAME_OF_TEST_FILE_IN_SAME_FOLDER

Running a specific test in the .yaml or other test file use:
--case #
where # is an integer signifying which test is to be run.


--------------------------------------
nondeallocs.sh
Used to analyze a meminfo file.
Does not cover all possible leaks.
