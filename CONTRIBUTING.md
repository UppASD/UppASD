# Contributing to UppASD

UppASD is an open software simulation package for the study of the magnetization at the atomic scale. 
The code is written mostly in `FORTRAN 90`, so please adhere to that standard. The code has been optimized mostly with the `Intel Fortran Compiler` in mind, making use of several of the functions and capabilities coming from `Intel MKL`. The code is of course compatible with other compilers, such as `gfortran`, but many of the optimizations done were not made with this compiler in mind.

## Git conventions and standard
- If you find a bug or are planning to add a feature [*write an issue about it*](https://guides.github.com/features/issues/). Of this way development can be traced, double work can be avoided and releases can be more easily managed.
- Please write [informative commit messages](https://www.freecodecamp.org/news/writing-good-commit-messages-a-practical-guide/). If you are addressing a given issue in the git, remember to [refer to it in your comment](https://about.gitlab.com/blog/2016/03/08/gitlab-tutorial-its-all-connected/).
- [Be mindful to tests](https://docs.gitlab.com/ee/development/testing_guide/best_practices.html):
   - If you introduce a new feature write a test about it. `UppASD` uses the script `tests/bergtest.py` to handle its tests with its reference values are in `yaml` files found in the `tests` folder.
   - If you fix a bug or change a behaviour that alters the values from any test, make sure that this is correct before updating the defaults in the test suite.
   - Make sure that your tests are added to the `.gitlab-ci.yml` to ensure that the CI/CD automatically runs your tests in the future.
- Please try to use [feature branches](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) as much as possible during development.
- The `develop` branch will keep the most up to date version of the code, with the latest changes pushed to it.
   - If you are implementing into develop remember to first `pull` the remote to the local, resolve conflicts, perform tests and then `push` your changes.
   - The agreed `git` practices for working in `develop` are:
      * Write your changes.
      * Perform `make asd-tests`, `make sld-tests`, `make gneb-tests` and `make regression-test`.
      * Perform `git pull`.
      * Solve any conflicts.
      * Perform `make asd-tests`, `make sld-tests`, `make gneb-tests` and `make regression-test` again.
      * If everything is okay do `git commit` and leave an explanatory comment of the changes.
      * After that `git push`.
- `UppASD` is a code under constant development, to avoid divergence between branches, remember to update your local implementation as often as possible.
- If you have a question the `UppASD` you can join the UppASD slack channel (contact the [administrators](ander.bergman@physics.uu.se) about how to join) where you can always ask questions and discuss your ideas.

## Coding conventions and formats

UppASD has the following conventions:

- `real` numbers are all of `kind=dblprec`, where `dblprec` is defined as `selected_real_kind(15, 307)`. Write **all** float numbers as `real(dblprec)`. Float constants should thus be written like `1.0_dblprec` or `1.0e-12_dblprec`. Also avoid intrinsics with fixed kind, like `dcos()`.
- The indentation level is set to **3** spaces.
- `UppASD` is not fixed format, so please adhere to that convention.
- General mathematical functions that can be used by many subroutines should be added to the `math_functions.f90` module found in `source/Tools/`.
   * Any function and or subroutine defined here **must** start as `f_dummy_name`.
- To define an array the following convention is used `real(dblprec), dimension(:,:,:), allocatable, optional, intent(in) :: dummy_1`
   * Always indicate the `intent` of the variable.
   * Pass the `dimension` option explicitly whenever possible.
- `UppASD` makes use of [Doxygen](http://www.doxygen.nl/) for automatic documentation, please document your routines.
   * Indicate a short description of the purpose of the piece of code with `@brief`, a more detailed explanation when needed with `@details` and the authors of the code with `@author`.
   * Please use the `!<` command to document **all** input/output variables for a routine.
- Whenever writing a routine for `UppASD` try for it to be **modular** i.e. as self contained as possible.
- Python code should be written in python3 following the [PEP8](https://www.python.org/dev/peps/pep-0008/) coding standard.
- **Always** test the code with the `make tests` command before pushing a revision to the repository.
