#
#                                  Maintain SD
#
SHELL = /bin/sh

# Program name
PROG = sd

# List of available building configs (found in make/default-profiles/systemname.make)
DEFAULT_SYSTEMS := gfortran gfortran-cuda gfortran-osx gfortran-cuda-osx ifort ifort-cuda ifort-nomkl ifort-cuda-nomkl pathscale pgf90 gfortran-win64
LOCAL_SYSTEMS := $(filter-out $(DEFAULT_SYSTEMS),$(shell ls ./source/make/user_profiles/*.make | sed 's/..source.make.user_profiles.//' | sed 's/.make//'))
SYSTEMS := $(DEFAULT_SYSTEMS) $(LOCAL_SYSTEMS)

.PHONY: deps PRINT nocopyprofile copyprofile help clean probe docs tests dist dist_minimal $(SYSTEMS)

# Including the help files
include ./source/make/makefileHELP
include ./source/make/makefileHELPDEFAULTS

PRINT:
	@echo $@ 

deps:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@python -B ./source/make/generateDependencies.py

probe:
	@python -B ./source/make/suggestProfiles.py

docs:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@cd ./docs; doxygen Doxyfile; pdflatex UppASDmanual.tex ; pdflatex UppASDmanual.tex

regression-test:
	@cd ./codeTester; python -B -u ./bergtest.py --file regression799.yaml | tee regression-tests.log ; \
	./cleanAll.sh

tests:
	@cd ./codeTester; python -B -u ./bergtest.py --file regulartests.yaml | tee tests.log ; ./cleanAll.sh

# Clean all .mod and .o files as well as mod and obj folders
clean:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	rm -f ./source/*.o ./source/*/*.o ./source/mod/*.mod ./source/$(PROG)

# Run same make file but with proper profile syntax and parallel make and print to log file
$(SYSTEMS):
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@cd ./source; $(MAKE) PROFILE=$@ 

# Generate an compressed archive
dist:
	@echo "Packaging source, examples, documentation, and tests to ./UppASD_dist.tar.gz"
	@cd codeTester ; ./cleanAll.sh ; cd ..
	@tar cf ./UppASD_dist.tar Makefile setup_UppASD.sh LICENSE AUTHORS README.md \
	./source/*.f90 ./source/*/*.f90 ./source/make/ ./source/gpu_files/ \
	./examples_revision_controlled ./docs/Doxyfile ./docs/*.pdf \
	./codeTester/ ; \
	gzip --best -f ./UppASD_dist.tar

dist_minimal:
	@echo "Packaging source to ./UppASD_src.tar.gz"
	@tar cf ./UppASD_src.tar Makefile setup_UppASD.sh LICENSE AUTHORS README.md \
	 ./source/*.f90 ./source/*/*.f90 ./source/make/ ./source/gpu_files/  ;\
	gzip --best -f ./UppASD_src.tar

