#ifndef ATT_TEST_SE3_JACOBIAN
#define ATT_TEST_SE3_JACOBIAN

#include "ang_test.h"
#include <jail/unit_test.h>

/*
 * This file contains tests to verify the jacobians of the different methods
 * of representing a transformation (speu_rodrigues, speu_dcm, homogeneous, twist, screw, dual).
 */

namespace ang {
namespace test {
	
class Tse3_jacobian: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tse3_jacobian(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
    void test_jac_right_left();

    void test_jac_right_inverse();
    void test_jac_right_composition_wrt_first();
    void test_jac_right_composition_wrt_second();
    void test_jac_right_forward_motion_wrt_motion();
    void test_jac_right_backward_motion_wrt_motion();
    void test_jac_right_log();
    void test_jac_right_plus_wrt_first();
    void test_jac_right_plus_wrt_second();
    void test_jac_right_minus_wrt_first();
    void test_jac_right_minus_wrt_second();
    void test_jac_right_forward_adjoint_wrt_motion();
    void test_jac_right_backward_adjoint_wrt_motion();

    void test_jac_left_inverse();
    void test_jac_left_composition_wrt_first();
    void test_jac_left_composition_wrt_second();
    void test_jac_left_forward_motion_wrt_motion();
    void test_jac_left_backward_motion_wrt_motion();
    void test_jac_left_log();
    void test_jac_left_plus_wrt_first();
    void test_jac_left_plus_wrt_second();
    void test_jac_left_minus_wrt_first();
    void test_jac_left_minus_wrt_second();
    void test_jac_left_forward_adjoint_wrt_motion();
    void test_jac_left_backward_adjoint_wrt_motion();

    void test_jac_euclidean_forward_motion_wrt_point();
    void test_jac_euclidean_backward_motion_wrt_point();
    void test_jac_euclidean_forward_adjoint_wrt_tangent();
    void test_jac_euclidean_backward_adjoint_wrt_tangent();

    void test_jac_euclidean_wrt_speu_dcm();
    void test_jac_euclidean_wrt_speu_dcm_linear();
    void test_jac_euclidean_wrt_speu_rodrigues();
    void test_jac_euclidean_wrt_homogeneous();
    void test_jac_euclidean_wrt_homogeneous_linear();
    void test_jac_euclidean_wrt_trfv();
    void test_jac_euclidean_wrt_trfv_bis();
};

} // closes namespace test
} // closes namespace ang

#endif

