#ifndef ATT_TEST_SO3_JACOBIAN
#define ATT_TEST_SO3_JACOBIAN

#include "ang_test.h"
#include "jail/unit_test.h"

/*
 * This file contains tests to verify the jacobians of the different rotation
 * representations (euler, dcm, rodrigues, rotv).
 */

namespace ang {
namespace test {
	
class Tso3_jacobian: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tso3_jacobian(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
    void test_jac_right_left();

    void test_jac_right_inverse();
    void test_jac_right_composition_wrt_first();
    void test_jac_right_composition_wrt_second();
    void test_jac_right_forward_rotation_wrt_rotation();
    void test_jac_right_backward_rotation_wrt_rotation();
    void test_jac_right_log();
    void test_jac_right_plus_wrt_first();
    void test_jac_right_plus_wrt_second();
    void test_jac_right_minus_wrt_first();
    void test_jac_right_minus_wrt_second();
    void test_jac_right_forward_adjoint_wrt_rotation();
    void test_jac_right_backward_adjoint_wrt_rotation();

    void test_jac_left_inverse();
    void test_jac_left_composition_wrt_first();
    void test_jac_left_composition_wrt_second();
    void test_jac_left_forward_rotation_wrt_rotation();
    void test_jac_left_backward_rotation_wrt_rotation();
    void test_jac_left_log();
    void test_jac_left_plus_wrt_first();
    void test_jac_left_plus_wrt_second();
    void test_jac_left_minus_wrt_first();
    void test_jac_left_minus_wrt_second();
    void test_jac_left_forward_adjoint_wrt_rotation();
    void test_jac_left_backward_adjoint_wrt_rotation();

    void test_jac_euclidean_forward_rotation_wrt_vector();
    void test_jac_euclidean_backward_rotation_wrt_vector();
    void test_jac_euclidean_forward_adjoint_wrt_tangent();
    void test_jac_euclidean_backward_adjoint_wrt_tangent();

    void test_jac_euclidean_wrt_dcm();
    void test_jac_euclidean_wrt_dcm_linear();
    void test_jac_euclidean_wrt_rodrigues();
    void test_jac_euclidean_wrt_rodrigues_bilinear();
    void test_jac_euclidean_wrt_rotv();

};

} // closes namespace test

} // closes namespace ang

#endif

