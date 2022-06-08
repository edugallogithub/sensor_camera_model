#ifndef MATH_TEST_INTERP
#define MATH_TEST_INTERP

#include "math_test.h"
#include <math/classifiers.h>
#include <jail/unit_test.h>

/*
This file is dedicated to the testing of cubic and spline interpolations, which 
are too complex as to be checked with the rest of the predicates.

The "test_search" test verifies the proper behavior of the interpolation functions
search_binary, search_equispaced, and find_index, which together obtain those 
members of a vector that should be used for the interpolation.

The "test_spline" test verifies the proper computation of splines by comparing
the results with those obtained in MatLab.

The "test_hermite_1d" verifies the proper behavior of the Hermite interpolation in
one dimension (both 1st and 2nd order differentials), checking the values and differentials
at the nodes, and the symmetry in the results with a symmetric table outside the nodes.

The "test_hermite_2d" verifies the proper behavior of the Hermite interpolation in
two dimensions (both 1st and 2nd order differentials), checking the values and differentials
at the nodes, and the symmetry in the results with a symmetric table outside the nodes.

The "test_hermite_3d" verifies the proper behavior of the Hermite interpolation in
three dimensions (both 1st and 2nd order differentials), checking the values and differentials
at the nodes, and the symmetry in the results with a symmetric table outside the nodes.

The "report_speed" test evaluates the time required to construct and evaluate tables
of different dimensions, sizes, and interpolation methods, creating a report with 
the results as a text file.

 */

namespace math {
namespace test {
	
class Tinterp: public ::jail::unit_test {
public:
	Tinterp(jail::counter&);
	/**< constructor based on counter */
	void run();	
	/**< execute tests and write results on console */

	void test_search();	
	void test_spline();
	void test_hermite_1d(math::logic::INTERP_MODE);
	void test_hermite_2d(math::logic::INTERP_MODE);
	void test_hermite_3d(math::logic::INTERP_MODE);

};

}; // closes namespace test
}; // closes namespace math
#endif

