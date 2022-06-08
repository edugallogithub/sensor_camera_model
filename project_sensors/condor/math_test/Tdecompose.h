#ifndef MATH_TEST_DECOMPOSE
#define MATH_TEST_DECOMPOSE

#include "math_test.h"
#include "jail/unit_test.h"


namespace math {
namespace test {
	
class Tdecompose: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tdecompose(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;
	/**< specific tests */

    void test_PartialPivLU();
    void test_FullPivLU();
    void test_LDLT();
    void test_BDCSVD();
    void test_matrix_exp();
    void test_linear_fit_lsq();
};

} // closes namespace test

} // closes namespace math
#endif

