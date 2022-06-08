#ifndef ATT_TEST_SE3
#define ATT_TEST_SE3

#include "ang_test.h"
#include <jail/unit_test.h>

/*
 * This file contains tests to verify the proper behavior of the different methods
 * of representing a transformation (speu_rodrigues, speu_dcm, homogeneous, twist, screw, dual).
 */

namespace ang {
namespace test {
	
class Tse3: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tse3(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
	void test_se3();            // rotation plus translation
    void test_exp_log_maps();
    void test_power();
    void test_sclerp();
    void test_plus_minus();
    void test_adjoint();
    void test_velocity();
    void test_dual_set();
};

} // closes namespace test
} // closes namespace ang

#endif

