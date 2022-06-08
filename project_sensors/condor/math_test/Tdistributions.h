#ifndef MATH_TEST_DISTRIBUTIONS
#define MATH_TEST_DISTRIBUTIONS

#include "math_test.h"
#include "jail/unit_test.h"

namespace math {
namespace test {
	
class Tdistributions: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tdistributions(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;
	/**< specific tests */

	void test_normal_distribution();
    void test_generator();
    void test_normal_one_distribution_two_generators();
    void test_random();
    void test_compare_generators();

};

} // closes namespace test

} // closes namespace math
#endif

