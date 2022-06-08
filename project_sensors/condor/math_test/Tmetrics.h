#ifndef MATH_TEST_METRICS
#define MATH_TEST_METRICS

#include "math_test.h"
#include "math/templates/metrics_.h"
#include "jail/unit_test.h"

/*
This file test the different metrics computations.
*/

namespace math {
namespace test {
	
class Tmetrics : public ::jail::unit_test {
public:
	explicit Tmetrics(jail::counter&);
	/**< constructor based on counter */
	void run() override;
	/**< execute tests and write results on console */

	void test_classic();
    void test_robust();
    void test_comparison_std_mad();
    void test_plot();
};

}; // closes namespace test
}; // closes namespace math
#endif

