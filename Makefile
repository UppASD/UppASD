#
#                                  Maintain SD
#
SHELL = /bin/sh


# Program name
PROG = sd

# List of available building configs (found in make/default-profiles/systemname.make)
DEFAULT_SYSTEMS := gfortran gfortran-cuda gfortran-osx gfortran-cuda-osx ifort ifort-cuda ifort-nomkl ifort-cuda-nomkl pathscale pgf90 pgf90-nomkl jureca gfortran-win64 gfortran-FC
LOCAL_SYSTEMS := $(filter-out $(DEFAULT_SYSTEMS),$(shell ls ./source/make/user_profiles/*.make | sed 's/..source.make.user_profiles.//' | sed 's/.make//'))
SYSTEMS := $(DEFAULT_SYSTEMS) $(LOCAL_SYSTEMS)

.PHONY: deps PRINT nocopyprofile copyprofile help clean probe docs tests asd-tests sld-tests dist dist_minimal gneb-tests $(SYSTEMS)

# Including the help files
include ./source/make/makefileHELP
include ./source/make/makefileHELPDEFAULTS

PRINT:
	@echo $@

deps:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@python3 ./source/make/generateDependencies.py

probe:
	@python3 ./source/make/suggestProfiles.py

docs:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@cd ./docs; doxygen Doxyfile; cd Manual; pdflatex UppASDmanual.tex ; pdflatex UppASDmanual.tex

tests:
	@echo ''
	@echo 'To run tests for selected functionalies, run:' 
	@echo '`make asd-tests`, `make sld-tests`, and/or `make gneb-tests`'
	@echo ''
	@echo 'For a quick regression test, run `make regression-test`'
	@echo ''

asd-tests:
	@cd ./tests; python3 -u ./bergtest.py --file regulartests.yaml | tee tests.log
	@cd ./tests; python3 -u ./bergtest.py --clean

sld-tests:
	@cd ./tests; python3 -u ./bergtest.py --file sldtests.yaml | tee tests.log
	@cd ./tests; python3 -u ./bergtest.py --clean

gneb-tests:
	@cd ./tests; python3 -u ./bergtest.py --file regressionGNEB_15d7169b.yaml | tee tests.log
	@cd ./tests; python3 -u ./bergtest.py --clean

regression-test:
	@cd ./tests; python3 -u ./bergtest.py --file regressionHisingen.yaml | tee regression-tests.log
	@cd ./tests; python3 -u ./bergtest.py --clean

# Clean all .mod and .o files as well as mod and obj folders
clean:
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	rm -f ./source/*.o ./source/*/*.o ./source/mod/*.mod ./source/$(PROG) ./source/*/*/*.o ./source/*/*/*/*.o 

# Run same make file but with proper profile syntax and parallel make and print to log file
$(SYSTEMS):
	@if [ ! -d source/make/user_profiles ] ; then mkdir source/make/user_profiles ; fi
	@cd ./source; $(MAKE) PROFILE=$@
	@if [ ! -d ./bin ] ; then mkdir ./bin ; fi
	@cp ./source/sd ./bin/sd.$@

# Generate an compressed archive
dist:
	@echo "Packaging source, examples, documentation, and tests to ./UppASD_dist.tar.gz"
	@cd tests ; ./cleanAll.sh ; cd ..
	@tar cf ./UppASD_dist.tar Makefile setup_UppASD.sh \
	./source/*.f90 ./source/*/*.f90 ./source/make/ ./source/gpu_files/ ./source/README/ \
	./source/Third_party/ \
	./examples ./docs/Doxyfile ./docs/*.pdf ASD_GUI/*.py \
	./docs/Manual/*.tex ./docs/Manual/*.ist ./docs/Manual/Pictures/*.png \
	./tests/ ; \
	gzip --best -f ./UppASD_dist.tar

dist_minimal:
	@echo "Packaging source to ./UppASD_src.tar.gz"
	@cd tests ; ./cleanAll.sh ; cd ..
	@tar cf ./UppASD_src.tar Makefile setup_UppASD.sh \
	./source/*.f90 ./source/*/*.f90 ./source/make/ ./source/gpu_files/ ./source/README/ \
	./source/Third_party/ ; \
	gzip --best -f ./UppASD_src.tar

sources:
	@echo "Packaging source to ./UppASD_src.tar.gz"
	@cd tests ; ./cleanAll.sh ; cd ..
	@tar cf ./UppASD_src.tar ./CMakeLists.txt \
	./source/*.f90 ./source/CMakeLists.txt ./source/*/CMakeLists.txt ./source/*/*.f90 ./source/gpu_files/ \
	./source/Third_party/  ./cmake
	@gzip --best -f ./UppASD_src.tar
