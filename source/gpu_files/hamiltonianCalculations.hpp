/* Written by Niklas Fejes based on UppASD FORTRAN version
 * 
 * Date: 2014/09/11 - Thomas Nystrand
 * Added dzyalonshinskii effect and moved heisenberg to separate function
 */

#ifndef __HAMILTONIAN_CALCULATIONS_HPP__
#define __HAMILTONIAN_CALCULATIONS_HPP__

#include "real_type.h"

#include "matrix.hpp"


class HamiltonianCalculations {
private:
	real                         beff_s[3];
	const matrix<real,2>         &ncoup;
	const matrix<unsigned int,2> &nlist;
	const unsigned int           *nlistsize;

	const unsigned int           do_dm;
	const matrix<real,3,3>       &dm_vect;
	const matrix<unsigned int,2> &dmlist;
	const unsigned int           *dmlistsize;
	

	// Contributions to the hamiltonian
	void heisenberg_field             (const size_t i, const size_t k, const matrix<real,3,3> &emomM, real *beff_s);
	void dzyalonshinskii_moriya_field (const size_t i, const size_t k, const matrix<real,3,3> &emomM, real *beff_s);

public:
	// Constructor
	HamiltonianCalculations(
			const matrix<real,2> 		&p1,
			const matrix<unsigned int,2> 	&p2,
			const unsigned int 		*p3,
			const unsigned int 		v4,
			const matrix<real,3,3>		&p5,
			const matrix<unsigned int,2> 	&p6,
			const unsigned int           	*p7)
			: ncoup(p1), nlist(p2), nlistsize(p3), do_dm(v4), dm_vect(p5), dmlist(p6), dmlistsize(p7) {}

	// Destructor
	~HamiltonianCalculations() {}

	// Variants
	void heisge_jij(
		matrix<real,3,3> &beff, 
		const matrix<real,3,3> &emomM,
		const matrix<real,3,3> &emom, 
		const matrix<real,3,3> &external_field);
};


#endif

