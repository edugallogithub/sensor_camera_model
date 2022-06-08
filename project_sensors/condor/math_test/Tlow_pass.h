#ifndef MATH_TEST_LOW_PASS
#define MATH_TEST_LOW_PASS

#include "math_test.h"
#include <math/classifiers.h>
#include <jail/unit_test.h>

/*
This file is dedicated to the testing of the low pass filter implementation.
 */

namespace math {
namespace test {
	
class Tlow_pass: public ::jail::unit_test {
public:
	Tlow_pass(jail::counter&);
	/**< constructor based on counter */
	void run();	
	/**< execute tests and write results on console */

	void test1();

};

}; // closes namespace test
}; // closes namespace math
#endif

