# Contributing to UppASD

UppASD is an open software simulation package for the study of the magnetization at the atomic scale. 
The code is written mostly in `FORTRAN 90`, so please adhere to that standard. The code has been optimized mostly with the `Intel Fortran Compiler` in mind, making use of several of the functions and capabilities coming from `Intel MKL`. The code is of course compatible with other compilers, such as `gfortran`, but many of the optimizations done were not made with this compiler in mind.

## Conventions and formats

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
   * Indicate a short description of the purpose of the piece of code with `@brief`, a more detailed explanation when needed with `@details` and the authrs of the code with `@author`.
   * Please use the `!<` command to document **all** input/output variables for a routine.
- Whenever writing a routine for `UppASD` try for it to be **modular** i.e. as self contained as possible.
- **Always** test the code with the `make tests` command before pushing a revision to the repository.
- The agreed `git` practices are:
   * Write your changes.
   * Perform `make tests` and `make regression-test`.
   * Perform `git pull`.
   * Solve any conflicts.
   * Perform `make tests` and `make regression-test` again.
   * If everything is okay do `git commit` and leave an explanatory comment of the changes.
   * After that `git push`.
  
