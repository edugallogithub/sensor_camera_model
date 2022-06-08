#include "Tso3_jacobian.h"

#include "../ang/rotate/dcm.h"
#include "../ang/rotate/rodrigues.h"
#include "../ang/rotate/euler.h"
#include "../ang/rotate/rotv.h"
#include "../ang/rotate/so3_tangent.h"
#include "../ang/tools.h"
#include <iostream>

using namespace std;

ang::test::Tso3_jacobian::Tso3_jacobian(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tso3_jacobian::run() {
	::jail::unit_test::run();

    test_jac_right_left();                              cout << endl << endl;

    test_jac_right_inverse();                           cout << endl << endl;
    test_jac_right_composition_wrt_first();             cout << endl << endl;
    test_jac_right_composition_wrt_second();            cout << endl << endl;
    test_jac_right_forward_rotation_wrt_rotation();     cout << endl << endl;
    test_jac_right_backward_rotation_wrt_rotation();    cout << endl << endl;
    test_jac_right_log() ;                              cout << endl << endl;
    test_jac_right_plus_wrt_first();                    cout << endl << endl;
    test_jac_right_plus_wrt_second();                   cout << endl << endl;
    test_jac_right_minus_wrt_first();                   cout << endl << endl;
    test_jac_right_minus_wrt_second();                  cout << endl << endl;
    test_jac_right_forward_adjoint_wrt_rotation();      cout << endl << endl;
    test_jac_right_backward_adjoint_wrt_rotation();     cout << endl << endl;

    test_jac_left_inverse();                            cout << endl << endl;
    test_jac_left_composition_wrt_first();              cout << endl << endl;
    test_jac_left_composition_wrt_second();             cout << endl << endl;
    test_jac_left_forward_rotation_wrt_rotation();      cout << endl << endl;
    test_jac_left_backward_rotation_wrt_rotation();     cout << endl << endl;
    test_jac_left_log() ;                               cout << endl << endl;
    test_jac_left_plus_wrt_first();                     cout << endl << endl;
    test_jac_left_plus_wrt_second();                    cout << endl << endl;
    test_jac_left_minus_wrt_first();                    cout << endl << endl;
    test_jac_left_minus_wrt_second();                   cout << endl << endl;
    test_jac_left_forward_adjoint_wrt_rotation();       cout << endl << endl;
    test_jac_left_backward_adjoint_wrt_rotation();      cout << endl << endl;

    test_jac_euclidean_forward_rotation_wrt_vector();   cout << endl << endl;
    test_jac_euclidean_backward_rotation_wrt_vector();  cout << endl << endl;
    test_jac_euclidean_forward_adjoint_wrt_tangent();   cout << endl << endl;
    test_jac_euclidean_backward_adjoint_wrt_tangent();  cout << endl << endl;

    test_jac_euclidean_wrt_dcm();                       cout << endl << endl;
    test_jac_euclidean_wrt_dcm_linear();                cout << endl << endl;
    test_jac_euclidean_wrt_rodrigues();                 cout << endl << endl;
    test_jac_euclidean_wrt_rodrigues_bilinear();        cout << endl << endl;
    test_jac_euclidean_wrt_rotv();                      cout << endl << endl;

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_left() {

    // It checks the proper behavior of the "rotv" right and left jacobian methods,
    // validating the following equations:
    // exp(rv + Deltarv) ~= exp(rv) plus (JR * Deltarv)
    // exp(rv + Deltarv) ~= (JL * Deltarv) plus_left exp(rv)
    // rv + JRinv * Deltarv ~= log[exp(rv) plus Deltarv]
    // rv + JLinv * Deltarv ~= log[Deltarv plus_left exp(rv)]

    rotv rv(0.4, -0.7, 0.3);
    dcm R = rv.exp_map_dcm();

    rotv Delta_rv(0.0002, 0.00025, -0.00015); // needs to be small
    dcm Delta_R = Delta_rv.exp_map_dcm();

    // I do an algebraic sum of the rotation vectors (not a concatenation - composition).
    // Then I compute the perturbation as the product of the jacobian by the rotation vector
    rotv Srv(rv() + Delta_rv());
    dcm SR = Srv.exp_map_dcm();

    rotv Delta_rv_right(rv.jac_right() * Delta_rv());
    dcm SR_right = R.plus_right(Delta_rv_right);

    check("jac right 0 - 0  ", SR()(0,0), SR_right()(0,0), 1e-8);
    check("jac right 0 - 1  ", SR()(0,1), SR_right()(0,1), 1e-8);
    check("jac right 0 - 2  ", SR()(0,2), SR_right()(0,2), 1e-8);
    check("jac right 1 - 0  ", SR()(1,0), SR_right()(1,0), 1e-8);
    check("jac right 1 - 1  ", SR()(1,1), SR_right()(1,1), 1e-8);
    check("jac right 1 - 2  ", SR()(1,2), SR_right()(1,2), 1e-8);
    check("jac right 2 - 0  ", SR()(2,0), SR_right()(2,0), 1e-8);
    check("jac right 2 - 1  ", SR()(2,1), SR_right()(2,1), 1e-8);
    check("jac right 2 - 2  ", SR()(2,2), SR_right()(2,2), 1e-8);

    rotv Delta_rv_left(rv.jac_left() * Delta_rv());
    dcm SR_left = R.plus_left(Delta_rv_left);

    check("jac left 0 - 0   ", SR()(0,0), SR_left()(0,0), 1e-8);
    check("jac left 0 - 1   ", SR()(0,1), SR_left()(0,1), 1e-8);
    check("jac left 0 - 2   ", SR()(0,2), SR_left()(0,2), 1e-8);
    check("jac left 1 - 0   ", SR()(1,0), SR_left()(1,0), 1e-8);
    check("jac left 1 - 1   ", SR()(1,1), SR_left()(1,1), 1e-8);
    check("jac left 1 - 2   ", SR()(1,2), SR_left()(1,2), 1e-8);
    check("jac left 2 - 0   ", SR()(2,0), SR_left()(2,0), 1e-8);
    check("jac left 2 - 1   ", SR()(2,1), SR_left()(2,1), 1e-8);
    check("jac left 2 - 2   ", SR()(2,2), SR_left()(2,2), 1e-8);

    Eigen::Vector3d rv_right_add = rv() + rv.jac_right_inv() * Delta_rv();
    Eigen::Vector3d rv_right_add2 = (rv.exp_map_dcm().plus_right(Delta_rv)).log_map()();

    check("jac right inv 0  ", rv_right_add(0), rv_right_add2(0), 1e-8);
    check("jac right inv 1  ", rv_right_add(1), rv_right_add2(1), 1e-8);
    check("jac right inv 2  ", rv_right_add(2), rv_right_add2(2), 1e-8);

    Eigen::Vector3d rv_left_add = rv() + rv.jac_left_inv() * Delta_rv();
    Eigen::Vector3d rv_left_add2 = (rv.exp_map_dcm().plus_left(Delta_rv)).log_map()();

    check("jac left inv 0   ", rv_left_add(0), rv_left_add2(0), 1e-8);
    check("jac left inv 1   ", rv_left_add(1), rv_left_add2(1), 1e-8);
    check("jac left inv 2   ", rv_left_add(2), rv_left_add2(2), 1e-8);

    // check that what is called inverse is in fact the inverse
    Eigen::Matrix3d J_right      = rv.jac_right();
    Eigen::Matrix3d J_right_inv  = rv.jac_right_inv();
    Eigen::Matrix3d J_right_inv2 = J_right.inverse();

    check("jac right inverse 0 - 0  ", J_right_inv(0,0), J_right_inv2(0,0), 1e-12);
    check("jac right inverse 0 - 1  ", J_right_inv(0,1), J_right_inv2(0,1), 1e-12);
    check("jac right inverse 0 - 2  ", J_right_inv(0,2), J_right_inv2(0,2), 1e-12);
    check("jac right inverse 1 - 0  ", J_right_inv(1,0), J_right_inv2(1,0), 1e-12);
    check("jac right inverse 1 - 1  ", J_right_inv(1,1), J_right_inv2(1,1), 1e-12);
    check("jac right inverse 1 - 2  ", J_right_inv(1,2), J_right_inv2(1,2), 1e-12);
    check("jac right inverse 2 - 0  ", J_right_inv(2,0), J_right_inv2(2,0), 1e-12);
    check("jac right inverse 2 - 1  ", J_right_inv(2,1), J_right_inv2(2,1), 1e-12);
    check("jac right inverse 2 - 2  ", J_right_inv(2,2), J_right_inv2(2,2), 1e-12);

    Eigen::Matrix3d J_left      = rv.jac_left();
    Eigen::Matrix3d J_left_inv  = rv.jac_left_inv();
    Eigen::Matrix3d J_left_inv2 = J_left.inverse();

    check("jac left inverse 0 - 0  ", J_left_inv(0,0), J_left_inv2(0,0), 1e-12);
    check("jac left inverse 0 - 1  ", J_left_inv(0,1), J_left_inv2(0,1), 1e-12);
    check("jac left inverse 0 - 2  ", J_left_inv(0,2), J_left_inv2(0,2), 1e-12);
    check("jac left inverse 1 - 0  ", J_left_inv(1,0), J_left_inv2(1,0), 1e-12);
    check("jac left inverse 1 - 1  ", J_left_inv(1,1), J_left_inv2(1,1), 1e-12);
    check("jac left inverse 1 - 2  ", J_left_inv(1,2), J_left_inv2(1,2), 1e-12);
    check("jac left inverse 2 - 0  ", J_left_inv(2,0), J_left_inv2(2,0), 1e-12);
    check("jac left inverse 2 - 1  ", J_left_inv(2,1), J_left_inv2(2,1), 1e-12);
    check("jac left inverse 2 - 2  ", J_left_inv(2,2), J_left_inv2(2,2), 1e-12);

} // closes test_jacobian_rotv_right_left

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_inverse() {

    // It verifies the behavior of the rotation classes jac_right_inverse method, validating the following:
    // (R plus DeltarB)^-1 = R^-1 plus J * DeltarB

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);
    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rodrigues Sq_nb  = q_nb.plus_right(DeltarB);
    rodrigues q_bn   = q_nb.inverse();
    rodrigues Sq_bn  = Sq_nb.inverse();
    Eigen::Matrix3d Jac_q = q_nb.jac_right_inverse();
    rotv Deltar_bn_q = ang::rotv(Jac_q * DeltarB());
    rodrigues SSq_bn = q_bn.plus_right(Deltar_bn_q);

    check("q jac right inverse 1 taylor    ", Sq_bn()(0), SSq_bn()(0), 1e-12);
    check("q jac right inverse 2 taylor    ", Sq_bn()(1), SSq_bn()(1), 1e-12);
    check("q jac right inverse 3 taylor    ", Sq_bn()(2), SSq_bn()(2), 1e-12);
    check("q jac right inverse 4 taylor    ", Sq_bn()(3), SSq_bn()(3), 1e-12);

    dcm R_bn         = R_nb.inverse();
    dcm SR_nb        = R_nb.plus_right(DeltarB);
    dcm SR_bn        = SR_nb.inverse();
    Eigen::Matrix3d Jac_R = R_nb.jac_right_inverse();
    rotv Deltar_bn_R = ang::rotv(Jac_R * DeltarB());
    dcm SSR_bn       = R_bn.plus_right(Deltar_bn_R);

    check("R jac right inverse 11 taylor   ", SR_bn()(0,0), SSR_bn()(0,0), 1e-12);
    check("R jac right inverse 12 taylor   ", SR_bn()(0,1), SSR_bn()(0,1), 1e-12);
    check("R jac right inverse 13 taylor   ", SR_bn()(0,2), SSR_bn()(0,2), 1e-12);
    check("R jac right inverse 21 taylor   ", SR_bn()(1,0), SSR_bn()(1,0), 1e-12);
    check("R jac right inverse 22 taylor   ", SR_bn()(1,1), SSR_bn()(1,1), 1e-12);
    check("R jac right inverse 23 taylor   ", SR_bn()(1,2), SSR_bn()(1,2), 1e-12);
    check("R jac right inverse 31 taylor   ", SR_bn()(2,0), SSR_bn()(2,0), 1e-12);
    check("R jac right inverse 32 taylor   ", SR_bn()(2,1), SSR_bn()(2,1), 1e-12);
    check("R jac right inverse 33 taylor   ", SR_bn()(2,2), SSR_bn()(2,2), 1e-12);

    rotv r_bn        = r_nb.inverse();
    rotv Sr_nb       = r_nb.plus_right(DeltarB);
    rotv Sr_bn       = Sr_nb.inverse();
    Eigen::Matrix3d Jac_r = r_nb.jac_right_inverse();
    rotv Deltar_bn_r = ang::rotv(Jac_r * DeltarB());
    rotv SSr_bn      = r_bn.plus_right(Deltar_bn_r);

    check("r jac right inverse 1 taylor    ", Sr_bn()(0), SSr_bn()(0), 1e-12);
    check("r jac right inverse 2 taylor    ", Sr_bn()(1), SSr_bn()(1), 1e-12);
    check("r jac right inverse 3 taylor    ", Sr_bn()(2), SSr_bn()(2), 1e-12);

} // closes test_jac_right_inverse

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_composition_wrt_first() {

    // It verifies the behavior of the rotation classes jac_right_composition_wrt_first method, validating the following:
    // (R1 plus DeltarB1) * R2 = R1 * R2 plus J * DeltarB1

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rotv r_bh(0.4, -1.2, 0.72);
    rodrigues q_bh(r_bh);
    dcm R_bh(r_bh);

    rodrigues Sq_nb       = q_nb.plus_right(DeltarB);
    rodrigues q_nh        = q_nb * q_bh;
    rodrigues Sq_nh       = Sq_nb * q_bh;
    Eigen::Matrix3d Jac_q = q_nb.jac_right_composition_wrt_first(q_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarB());
    rodrigues SSq_nh      = q_nh.plus_right(Deltar_bh_q);

    check("q jac right composition first 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-12);
    check("q jac right composition first 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-12);
    check("q jac right composition first 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-12);
    check("q jac right composition first 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-12);

    dcm SR_nb             = R_nb.plus_right(DeltarB);
    dcm R_nh              = R_nb * R_bh;
    dcm SR_nh             = SR_nb * R_bh;
    Eigen::Matrix3d Jac_R = R_nb.jac_right_composition_wrt_first(R_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarB());
    dcm SSR_nh            = R_nh.plus_right(Deltar_bh_R);

    check("R jac right composition first 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac right composition first 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac right composition first 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac right composition first 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac right composition first 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac right composition first 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac right composition first 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac right composition first 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac right composition first 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-12);

    rotv Sr_nb            = r_nb.plus_right(DeltarB);
    rotv r_nh             = r_nb * r_bh;
    rotv Sr_nh            = Sr_nb * r_bh;
    Eigen::Matrix3d Jac_r = r_nb.jac_right_composition_wrt_first(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarB());
    rotv SSr_nh           = r_nh.plus_right(Deltar_bh_r);

    check("r jac right composition first 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-12);
    check("r jac right composition first 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-12);
    check("r jac right composition first 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_right_composition_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_composition_wrt_second() {

    // It verifies the behavior of the rotation classes jac_right_composition_wrt_second method, validating the following:
    // R1 * (R2 plus DeltarB2) = R1 * R2 plus J * DeltarB2

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv r_bh(0.4, -1.2, 0.72);
    rodrigues q_bh(r_bh);
    dcm R_bh(r_bh);

    rotv DeltarH(1e-3, -5e-4, 2e-3);

    rodrigues Sq_bh       = q_bh.plus_right(DeltarH);
    rodrigues q_nh        = q_nb * q_bh;
    rodrigues Sq_nh       = q_nb * Sq_bh;
    Eigen::Matrix3d Jac_q = q_nb.jac_right_composition_wrt_second(q_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarH());
    rodrigues SSq_nh      = q_nh.plus_right(Deltar_bh_q);

    check("q jac right composition second 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-12);
    check("q jac right composition second 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-12);
    check("q jac right composition second 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-12);
    check("q jac right composition second 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-12);

    dcm SR_bh             = R_bh.plus_right(DeltarH);
    dcm R_nh              = R_nb * R_bh;
    dcm SR_nh             = R_nb * SR_bh;
    Eigen::Matrix3d Jac_R = R_nb.jac_right_composition_wrt_second(R_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarH());
    dcm SSR_nh            = R_nh.plus_right(Deltar_bh_R);

    check("R jac right composition second 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac right composition second 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac right composition second 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac right composition second 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac right composition second 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac right composition second 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac right composition second 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac right composition second 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac right composition second 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-12);

    rotv Sr_bh            = r_bh.plus_right(DeltarH);
    rotv r_nh             = r_nb * r_bh;
    rotv Sr_nh            = r_nb * Sr_bh;
    Eigen::Matrix3d Jac_r = r_nb.jac_right_composition_wrt_second(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarH());
    rotv SSr_nh           = r_nh.plus_right(Deltar_bh_r);

    check("r jac right composition second 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-12);
    check("r jac right composition second 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-12);
    check("r jac right composition second 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_right_composition_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_forward_rotation_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_right_forward_rotation_wrt_rotation method, validating the following:
    // (R plus DeltarB) * v = R * v + J * DeltarB

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vB(0.33, -0.65, 0.17);

    Eigen::Vector3d vN_q      = q_nb * vB;
    rodrigues Sq_nb           = q_nb.plus_right(DeltarB);
    Eigen::Vector3d SvN_q     = Sq_nb * vB;
    Eigen::Matrix3d Jac_q     = q_nb.jac_right_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_q = Jac_q * DeltarB();
    Eigen::Vector3d SSvN_q    = vN_q + DeltavN_q;

    check("q jac right forward rotation rotation 1 taylor    ", SvN_q(0), SSvN_q(0), 1e-5);
    check("q jac right forward rotation rotation 2 taylor    ", SvN_q(1), SSvN_q(1), 1e-5);
    check("q jac right forward rotation rotation 3 taylor    ", SvN_q(2), SSvN_q(2), 1e-5);

    Eigen::Vector3d vN_R      = R_nb * vB;
    dcm SR_nb                 = R_nb.plus_right(DeltarB);
    Eigen::Vector3d SvN_R     = SR_nb * vB;
    Eigen::Matrix3d Jac_R     = R_nb.jac_right_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_R = Jac_R * DeltarB();
    Eigen::Vector3d SSvN_R    = vN_R + DeltavN_R;

    check("R jac right forward rotation rotation 1 taylor    ", SvN_R(0), SSvN_R(0), 1e-5);
    check("R jac right forward rotation rotation 2 taylor    ", SvN_R(1), SSvN_R(1), 1e-5);
    check("R jac right forward rotation rotation 3 taylor    ", SvN_R(2), SSvN_R(2), 1e-5);

    Eigen::Vector3d vN_r      = r_nb * vB;
    rotv Sr_nb                = r_nb.plus_right(DeltarB);
    Eigen::Vector3d SvN_r     = Sr_nb * vB;
    Eigen::Matrix3d Jac_r     = r_nb.jac_right_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_r = Jac_r * DeltarB();
    Eigen::Vector3d SSvN_r    = vN_r + DeltavN_r;

    check("r jac right forward rotation rotation 1 taylor    ", SvN_r(0), SSvN_r(0), 1e-5);
    check("r jac right forward rotation rotation 2 taylor    ", SvN_r(1), SSvN_r(1), 1e-5);
    check("r jac right forward rotation rotation 3 taylor    ", SvN_r(2), SSvN_r(2), 1e-5);
} // closes test_jac_right_forward_rotation_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_backward_rotation_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_right_backward_rotation_wrt_rotation method, validating the following:
    // (R plus DeltarB) / v = R / v + J * DeltarB

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vN(0.33, -0.65, 0.17);

    Eigen::Vector3d vB_q      = q_nb / vN;
    rodrigues Sq_nb           = q_nb.plus_right(DeltarB);
    Eigen::Vector3d SvB_q     = Sq_nb / vN;
    Eigen::Matrix3d Jac_q     = q_nb.jac_right_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_q = Jac_q * DeltarB();
    Eigen::Vector3d SSvB_q    = vB_q + DeltavB_q;

    check("q jac right backward rotation rotation 1 taylor    ", SvB_q(0), SSvB_q(0), 1e-5);
    check("q jac right backward rotation rotation 2 taylor    ", SvB_q(1), SSvB_q(1), 1e-5);
    check("q jac right backward rotation rotation 3 taylor    ", SvB_q(2), SSvB_q(2), 1e-5);

    Eigen::Vector3d vB_R      = R_nb / vN;
    dcm SR_nb                 = R_nb.plus_right(DeltarB);
    Eigen::Vector3d SvB_R     = SR_nb / vN;
    Eigen::Matrix3d Jac_R     = R_nb.jac_right_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_R = Jac_R * DeltarB();
    Eigen::Vector3d SSvB_R    = vB_R + DeltavB_R;

    check("R jac right backward rotation rotation 1 taylor    ", SvB_R(0), SSvB_R(0), 1e-5);
    check("R jac right backward rotation rotation 2 taylor    ", SvB_R(1), SSvB_R(1), 1e-5);
    check("R jac right backward rotation rotation 3 taylor    ", SvB_R(2), SSvB_R(2), 1e-5);

    Eigen::Vector3d vB_r      = r_nb / vN;
    rotv Sr_nb                = r_nb.plus_right(DeltarB);
    Eigen::Vector3d SvB_r     = Sr_nb / vN;
    Eigen::Matrix3d Jac_r     = r_nb.jac_right_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_r = Jac_r * DeltarB();
    Eigen::Vector3d SSvB_r    = vB_r + DeltavB_r;

    check("r jac right backward rotation rotation 1 taylor    ", SvB_r(0), SSvB_r(0), 1e-5);
    check("r jac right backward rotation rotation 2 taylor    ", SvB_r(1), SSvB_r(1), 1e-5);
    check("r jac right backward rotation rotation 3 taylor    ", SvB_r(2), SSvB_r(2), 1e-5);
} // closes test_jac_right_backward_rotation_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_log() {

    // It checks the proper behavior of the "rotv" right logarithmic jacobian method,
    // validating the following equations:
    // Log(R plus Deltarv) ~= Log(R) + (J * Deltarv)
    // Log(q plus Deltarv) ~= Log(q) + (J * Deltarv)
    // Log(rv plus Deltarv) ~= Log(rv) + (J * Deltarv)

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);
    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rodrigues Sq_nb       = q_nb.plus_right(DeltarB);
    rotv Sr_nb_q          = Sq_nb.log_map();
    rotv r_nb_q           = q_nb.log_map();
    Eigen::Matrix3d Jac_q = q_nb.jac_right_log();
    rotv Deltar_nb_q      = rotv(Jac_q * DeltarB());
    rotv SSr_nb_q(r_nb_q() + Deltar_nb_q());

    check("q jac right log 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-6);
    check("q jac right log 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-6);
    check("q jac right log 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-6);

    dcm SR_nb             = R_nb.plus_right(DeltarB);
    rotv Sr_nb_R          = SR_nb.log_map();
    rotv r_nb_R           = R_nb.log_map();
    Eigen::Matrix3d Jac_R = R_nb.jac_right_log();
    rotv Deltar_nb_R      = rotv(Jac_R * DeltarB());
    rotv SSr_nb_R(r_nb_R() + Deltar_nb_R());

    check("R jac right log 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-6);
    check("R jac right log 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-6);
    check("R jac right log 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-6);

    rotv Sr_nb_r          = r_nb.plus_right(DeltarB);
    Eigen::Matrix3d Jac_r = r_nb.jac_right_log();
    rotv Deltar_nb_r      = rotv(Jac_r * DeltarB());
    rotv SSr_nb_r(r_nb() + Deltar_nb_r());

    check("r jac right log 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-6);
    check("r jac right log 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-6);
    check("r jac right log 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-6);

} // closes test_jacobian_right_log

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_plus_wrt_first() {

    // It verifies the behavior of the rotation classes jac_right_plus_wrt_first method, validating the following:
    // (R1 plus Deltar1B) plus rv2 = (R1 plus rv2) plus J * Deltar1B

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rotv r_bh(0.4, -1.2, 0.72);

    rodrigues q_nh        = q_nb.plus_right(r_bh);
    rodrigues Sq_nb       = q_nb.plus_right(DeltarB);
    rodrigues Sq_nh       = Sq_nb.plus_right(r_bh);
    Eigen::Matrix3d Jac_q = q_nb.jac_right_plus_wrt_first(r_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarB());
    rodrigues SSq_nh      = q_nh.plus_right(Deltar_bh_q);

    check("q jac right plus first 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-12);
    check("q jac right plus first 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-12);
    check("q jac right plus first 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-12);
    check("q jac right plus first 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-12);

    dcm R_nh              = R_nb.plus_right(r_bh);
    dcm SR_nb             = R_nb.plus_right(DeltarB);
    dcm SR_nh             = SR_nb.plus_right(r_bh);
    Eigen::Matrix3d Jac_R = R_nb.jac_right_plus_wrt_first(r_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarB());
    dcm SSR_nh            = R_nh.plus_right(Deltar_bh_R);

    check("R jac right plus first 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac right plus first 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac right plus first 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac right plus first 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac right plus first 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac right plus first 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac right plus first 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac right plus first 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac right plus first 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-12);

    rotv r_nh             = r_nb.plus_right(r_bh);
    rotv Sr_nb            = r_nb.plus_right(DeltarB);
    rotv Sr_nh            = Sr_nb.plus_right(r_bh);
    Eigen::Matrix3d Jac_r = r_nb.jac_right_plus_wrt_first(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarB());
    rotv SSr_nh           = r_nh.plus_right(Deltar_bh_r);

    check("r jac right plus first 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-12);
    check("r jac right plus first 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-12);
    check("r jac right plus first 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_right_plus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_plus_wrt_second() {

    // It verifies the behavior of the rotation classes jac_right_plus_wrt_second method, validating the following:
    // R1 plus (rv2 + Deltar2B) = (R1 plus rv2) plus J * Deltar2B

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarH(1e-3, -5e-4, 2e-3);

    rotv r_bh(0.4, -1.2, 0.72);
    rotv rr_bh(r_bh() + DeltarH());

    rodrigues q_nh        = q_nb.plus_right(r_bh);
    rodrigues Sq_nh       = q_nb.plus_right(rr_bh);
    Eigen::Matrix3d Jac_q = q_nb.jac_right_plus_wrt_second(r_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarH());
    rodrigues SSq_nh      = q_nh.plus_right(Deltar_bh_q);

    check("q jac right plus second 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-6);
    check("q jac right plus second 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-6);
    check("q jac right plus second 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-6);
    check("q jac right plus second 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-6);

    dcm R_nh              = R_nb.plus_right(r_bh);
    dcm SR_nh             = R_nb.plus_right(rr_bh);
    Eigen::Matrix3d Jac_R = R_nb.jac_right_plus_wrt_second(r_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarH());
    dcm SSR_nh            = R_nh.plus_right(Deltar_bh_R);

    check("R jac right plus second 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-6);
    check("R jac right plus second 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-6);
    check("R jac right plus second 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-6);
    check("R jac right plus second 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-6);
    check("R jac right plus second 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-6);
    check("R jac right plus second 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-6);
    check("R jac right plus second 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-6);
    check("R jac right plus second 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-6);
    check("R jac right plus second 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-6);

    rotv r_nh             = r_nb.plus_right(r_bh);
    rotv Sr_nh            = r_nb.plus_right(rr_bh);
    Eigen::Matrix3d Jac_r = r_nb.jac_right_plus_wrt_second(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarH());
    rotv SSr_nh           = r_nh.plus_right(Deltar_bh_r);

    check("r jac right plus second 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-6);
    check("r jac right plus second 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-6);
    check("r jac right plus second 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-6);

} // closes test_jac_right_plus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_minus_wrt_first() {

    // It verifies the behavior of the rotation classes jac_right_minus_wrt_first method, validating the following:
    //* (R2 plus Deltar2B) minus R1 = (R2 minus R1) + J * Deltar2B

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rotv Pr_nb      = r_nb.plus_right(DeltarB);
    rodrigues Pq_nb = q_nb.plus_right(DeltarB);
    dcm PR_nb       = R_nb.plus_right(DeltarB);

    rotv DeltarB_q        = Pq_nb.minus_right(q_nb);
    rodrigues Sq_nb       = Pq_nb.plus_right(DeltarB);
    rotv Sr_nb_q          = Sq_nb.minus_right(q_nb);
    Eigen::Matrix3d Jac_q = Pq_nb.jac_right_minus_wrt_first(q_nb);
    rotv Deltar_nb_q      = ang::rotv(Jac_q * DeltarB());
    rotv SSr_nb_q(DeltarB_q() + Deltar_nb_q());

    check("q jac right minus first 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-12);
    check("q jac right minus first 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-12);
    check("q jac right minus first 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-12);

    rotv DeltarB_R        = PR_nb.minus_right(R_nb);
    dcm SR_nb             = PR_nb.plus_right(DeltarB);
    rotv Sr_nb_R          = SR_nb.minus_right(R_nb);
    Eigen::Matrix3d Jac_R = PR_nb.jac_right_minus_wrt_first(R_nb);
    rotv Deltar_nb_R      = ang::rotv(Jac_R * DeltarB());
    rotv SSr_nb_R(DeltarB_R() + Deltar_nb_R());

    check("R jac right minus first 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-12);
    check("R jac right minus first 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-12);
    check("R jac right minus first 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-12);

    rotv DeltarB_r        = Pr_nb.minus_right(r_nb);
    rotv Sr_nb            = Pr_nb.plus_right(DeltarB);
    rotv Sr_nb_r          = Sr_nb.minus_right(r_nb);
    Eigen::Matrix3d Jac_r = Pr_nb.jac_right_minus_wrt_first(r_nb);
    rotv Deltar_nb_r      = ang::rotv(Jac_r * DeltarB());
    rotv SSr_nb_r(DeltarB_r() + Deltar_nb_r());

    check("r jac right minus first 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-12);
    check("r jac right minus first 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-12);
    check("r jac right minus first 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-12);

} // closes test_jac_right_minus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_minus_wrt_second() {

    // It verifies the behavior of the rotation classes jac_right_minus_wrt_second method, validating the following:
    // R2 minus(R1 plus Deltar1B) = (R2 minus R1) + J * Deltar1B

    rotv r1_nb(-0.7, 0.45, 1.0);
    rodrigues q1_nb(r1_nb);
    dcm R1_nb(r1_nb);

    rotv DeltarBA(1e-3, -5e-4, 2e-3);

    rotv r2_nb      = r1_nb.plus_right(DeltarBA);
    rodrigues q2_nb = q1_nb.plus_right(DeltarBA);
    dcm R2_nb       = R1_nb.plus_right(DeltarBA);

    rotv DeltarB(-3.2e-4, -1.7e-4, 2.2e-3);

    rotv DeltarB_q        = q2_nb.minus_right(q1_nb);
    rodrigues Sq_nb       = q1_nb.plus_right(DeltarB);
    rotv Sr_nb_q          = q2_nb.minus_right(Sq_nb);
    Eigen::Matrix3d Jac_q = q2_nb.jac_right_minus_wrt_second(q1_nb);
    rotv Deltar_nb_q      = ang::rotv(Jac_q * DeltarB());
    rotv SSr_nb_q(DeltarB_q() + Deltar_nb_q());

    check("q jac right minus second 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-7);
    check("q jac right minus second 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-7);
    check("q jac right minus second 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-7);

    rotv DeltarB_R        = R2_nb.minus_right(R1_nb);
    dcm SR_nb             = R1_nb.plus_right(DeltarB);
    rotv Sr_nb_R          = R2_nb.minus_right(SR_nb);
    Eigen::Matrix3d Jac_R = R2_nb.jac_right_minus_wrt_second(R1_nb);
    rotv Deltar_nb_R      = ang::rotv(Jac_R * DeltarB());
    rotv SSr_nb_R(DeltarB_R() + Deltar_nb_R());

    check("R jac right minus second 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-7);
    check("R jac right minus second 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-7);
    check("R jac right minus second 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-7);

    rotv DeltarB_r        = r2_nb.minus_right(r1_nb);
    rotv Sr_nb            = r1_nb.plus_right(DeltarB);
    rotv Sr_nb_r          = r2_nb.minus_right(Sr_nb);
    Eigen::Matrix3d Jac_r = r2_nb.jac_right_minus_wrt_second(r1_nb);
    rotv Deltar_nb_r      = ang::rotv(Jac_r * DeltarB());
    rotv SSr_nb_r(DeltarB_r() + Deltar_nb_r());

    check("r jac right minus second 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-7);
    check("r jac right minus second 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-7);
    check("r jac right minus second 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-7);

} // closes test_jac_right_minus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_forward_adjoint_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_right_forward_adjoint_wrt_rotation method, validating the following:
    // Ad(R plus DeltarB) | w = AdR | w + J * DeltarB

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    so3_tangent wB(0.33, -0.65, 0.17);

    so3_tangent wN_q          = q_nb | wB;
    rodrigues Sq_nb           = q_nb.plus_right(DeltarB);
    so3_tangent SwN_q         = Sq_nb | wB;
    Eigen::Matrix3d Jac_q     = q_nb.jac_right_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_q = Jac_q * DeltarB();
    so3_tangent SSwN_q        = wN_q + DeltawN_q;

    check("q jac right forward adjoint rotation 1 taylor    ", SwN_q()(0), SSwN_q()(0), 1e-5);
    check("q jac right forward adjoint rotation 2 taylor    ", SwN_q()(1), SSwN_q()(1), 1e-5);
    check("q jac right forward adjoint rotation 3 taylor    ", SwN_q()(2), SSwN_q()(2), 1e-5);

    so3_tangent wN_R          = R_nb | wB;
    dcm SR_nb                 = R_nb.plus_right(DeltarB);
    so3_tangent SwN_R         = SR_nb | wB;
    Eigen::Matrix3d Jac_R     = R_nb.jac_right_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_R = Jac_R * DeltarB();
    so3_tangent SSwN_R        = wN_R + DeltawN_R;

    check("R jac right forward adjoint rotation 1 taylor    ", SwN_R()(0), SSwN_R()(0), 1e-5);
    check("R jac right forward adjoint rotation 2 taylor    ", SwN_R()(1), SSwN_R()(1), 1e-5);
    check("R jac right forward adjoint rotation 3 taylor    ", SwN_R()(2), SSwN_R()(2), 1e-5);

    so3_tangent wN_r          = r_nb | wB;
    rotv Sr_nb                = r_nb.plus_right(DeltarB);
    so3_tangent SwN_r         = Sr_nb | wB;
    Eigen::Matrix3d Jac_r     = r_nb.jac_right_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_r = Jac_r * DeltarB();
    so3_tangent SSwN_r        = wN_r + DeltawN_r;

    check("r jac right forward adjoint rotation 1 taylor    ", SwN_r()(0), SSwN_r()(0), 1e-5);
    check("r jac right forward adjoint rotation 2 taylor    ", SwN_r()(1), SSwN_r()(1), 1e-5);
    check("r jac right forward adjoint rotation 3 taylor    ", SwN_r()(2), SSwN_r()(2), 1e-5);
} // closes test_jac_right_forward_adjoint_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_right_backward_adjoint_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_right_backward_adjoint_wrt_rotation method, validating the following:
    // Ad(R plus DeltarB) % w = AdR % w + J * DeltarB

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    so3_tangent wN(0.33, -0.65, 0.17);

    so3_tangent wB_q          = q_nb % wN;
    rodrigues Sq_nb           = q_nb.plus_right(DeltarB);
    so3_tangent SwB_q         = Sq_nb % wN;
    Eigen::Matrix3d Jac_q     = q_nb.jac_right_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_q = Jac_q * DeltarB();
    so3_tangent SSwB_q        = wB_q + DeltawB_q;

    check("q jac right backward adjoint rotation 1 taylor    ", SwB_q()(0), SSwB_q()(0), 1e-5);
    check("q jac right backward adjoint rotation 2 taylor    ", SwB_q()(1), SSwB_q()(1), 1e-5);
    check("q jac right backward adjoint rotation 3 taylor    ", SwB_q()(2), SSwB_q()(2), 1e-5);

    so3_tangent wB_R          = R_nb % wN;
    dcm SR_nb                 = R_nb.plus_right(DeltarB);
    so3_tangent SwB_R         = SR_nb % wN;
    Eigen::Matrix3d Jac_R     = R_nb.jac_right_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_R = Jac_R * DeltarB();
    so3_tangent SSwB_R        = wB_R + DeltawB_R;

    check("R jac right backward adjoint rotation 1 taylor    ", SwB_R()(0), SSwB_R()(0), 1e-5);
    check("R jac right backward adjoint rotation 2 taylor    ", SwB_R()(1), SSwB_R()(1), 1e-5);
    check("R jac right backward adjoint rotation 3 taylor    ", SwB_R()(2), SSwB_R()(2), 1e-5);

    so3_tangent wB_r          = r_nb % wN;
    rotv Sr_nb                = r_nb.plus_right(DeltarB);
    so3_tangent SwB_r         = Sr_nb % wN;
    Eigen::Matrix3d Jac_r     = r_nb.jac_right_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_r = Jac_r * DeltarB();
    so3_tangent SSwB_r        = wB_r + DeltawB_r;

    check("r jac right backward adjoint rotation 1 taylor    ", SwB_r()(0), SSwB_r()(0), 1e-5);
    check("r jac right backward adjoint rotation 2 taylor    ", SwB_r()(1), SSwB_r()(1), 1e-5);
    check("r jac right backward adjoint rotation 3 taylor    ", SwB_r()(2), SSwB_r()(2), 1e-5);
} // closes test_jac_right_backward_adjoint_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_inverse() {

    // It verifies the behavior of the rotation classes jac_left_inverse method, validating the following:
    // (DeltarN plus R)^-1 = J * DeltarN plus R^-1

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);
    rotv DeltarN(1e-3, -5e-4, 2e-3);

    rodrigues Sq_nb  = q_nb.plus_left(DeltarN);
    rodrigues q_bn   = q_nb.inverse();
    rodrigues Sq_bn  = Sq_nb.inverse();
    Eigen::Matrix3d Jac_q = q_nb.jac_left_inverse();
    rotv Deltar_bn_q = ang::rotv(Jac_q * DeltarN());
    rodrigues SSq_bn = q_bn.plus_left(Deltar_bn_q);

    check("q jac left inverse 1 taylor    ", Sq_bn()(0), SSq_bn()(0), 1e-12);
    check("q jac left inverse 2 taylor    ", Sq_bn()(1), SSq_bn()(1), 1e-12);
    check("q jac left inverse 3 taylor    ", Sq_bn()(2), SSq_bn()(2), 1e-12);
    check("q jac left inverse 4 taylor    ", Sq_bn()(3), SSq_bn()(3), 1e-12);

    dcm R_bn         = R_nb.inverse();
    dcm SR_nb        = R_nb.plus_left(DeltarN);
    dcm SR_bn        = SR_nb.inverse();
    Eigen::Matrix3d Jac_R = R_nb.jac_left_inverse();
    rotv Deltar_bn_R = ang::rotv(Jac_R * DeltarN());
    dcm SSR_bn       = R_bn.plus_left(Deltar_bn_R);

    check("R jac left inverse 11 taylor   ", SR_bn()(0,0), SSR_bn()(0,0), 1e-12);
    check("R jac left inverse 12 taylor   ", SR_bn()(0,1), SSR_bn()(0,1), 1e-12);
    check("R jac left inverse 13 taylor   ", SR_bn()(0,2), SSR_bn()(0,2), 1e-12);
    check("R jac left inverse 21 taylor   ", SR_bn()(1,0), SSR_bn()(1,0), 1e-12);
    check("R jac left inverse 22 taylor   ", SR_bn()(1,1), SSR_bn()(1,1), 1e-12);
    check("R jac left inverse 23 taylor   ", SR_bn()(1,2), SSR_bn()(1,2), 1e-12);
    check("R jac left inverse 31 taylor   ", SR_bn()(2,0), SSR_bn()(2,0), 1e-12);
    check("R jac left inverse 32 taylor   ", SR_bn()(2,1), SSR_bn()(2,1), 1e-12);
    check("R jac left inverse 33 taylor   ", SR_bn()(2,2), SSR_bn()(2,2), 1e-12);

    rotv r_bn        = r_nb.inverse();
    rotv Sr_nb       = r_nb.plus_left(DeltarN);
    rotv Sr_bn       = Sr_nb.inverse();
    Eigen::Matrix3d Jac_r = r_nb.jac_left_inverse();
    rotv Deltar_bn_r = ang::rotv(Jac_r * DeltarN());
    rotv SSr_bn      = r_bn.plus_left(Deltar_bn_r);

    check("r jac left inverse 1 taylor    ", Sr_bn()(0), SSr_bn()(0), 1e-12);
    check("r jac left inverse 2 taylor    ", Sr_bn()(1), SSr_bn()(1), 1e-12);
    check("r jac left inverse 3 taylor    ", Sr_bn()(2), SSr_bn()(2), 1e-12);

} // closes test_jac_left_inverse

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_composition_wrt_first() {

    // It verifies the behavior of the rotation classes jac_left_composition_wrt_first method, validating the following:
    // (DeltarN1 plus R1) * R2 = J * DeltarB1 plus R1 * R2

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    rotv r_bh(0.4, -1.2, 0.72);
    rodrigues q_bh(r_bh);
    dcm R_bh(r_bh);

    rodrigues Sq_nb       = q_nb.plus_left(DeltarN);
    rodrigues q_nh        = q_nb * q_bh;
    rodrigues Sq_nh       = Sq_nb * q_bh;
    Eigen::Matrix3d Jac_q = q_nb.jac_left_composition_wrt_first(q_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarN());
    rodrigues SSq_nh      = q_nh.plus_left(Deltar_bh_q);

    check("q jac left composition first 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-12);
    check("q jac left composition first 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-12);
    check("q jac left composition first 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-12);
    check("q jac left composition first 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-12);

    dcm SR_nb             = R_nb.plus_left(DeltarN);
    dcm R_nh              = R_nb * R_bh;
    dcm SR_nh             = SR_nb * R_bh;
    Eigen::Matrix3d Jac_R = R_nb.jac_left_composition_wrt_first(R_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarN());
    dcm SSR_nh            = R_nh.plus_left(Deltar_bh_R);

    check("R jac left composition first 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac left composition first 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac left composition first 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac left composition first 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac left composition first 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac left composition first 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac left composition first 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac left composition first 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac left composition first 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-12);

    rotv Sr_nb            = r_nb.plus_left(DeltarN);
    rotv r_nh             = r_nb * r_bh;
    rotv Sr_nh            = Sr_nb * r_bh;
    Eigen::Matrix3d Jac_r = r_nb.jac_left_composition_wrt_first(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarN());
    rotv SSr_nh           = r_nh.plus_left(Deltar_bh_r);

    check("r jac left composition first 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-12);
    check("r jac left composition first 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-12);
    check("r jac left composition first 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_left_composition_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_composition_wrt_second() {

    // It verifies the behavior of the rotation classes jac_left_composition_wrt_second method, validating the following:
    // R1 * (DeltarN2 plus R2) = J * DeltarB2 plus R1 * R2

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv r_bh(0.4, -1.2, 0.72);
    rodrigues q_bh(r_bh);
    dcm R_bh(r_bh);

    rotv DeltarB(1e-3, -5e-4, 2e-3);

    rodrigues Sq_bh       = q_bh.plus_left(DeltarB);
    rodrigues q_nh        = q_nb * q_bh;
    rodrigues Sq_nh       = q_nb * Sq_bh;
    Eigen::Matrix3d Jac_q = q_nb.jac_left_composition_wrt_second(q_bh);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarB());
    rodrigues SSq_nh      = q_nh.plus_left(Deltar_bh_q);

    check("q jac left composition second 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-12);
    check("q jac left composition second 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-12);
    check("q jac left composition second 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-12);
    check("q jac left composition second 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-12);

    dcm SR_bh             = R_bh.plus_left(DeltarB);
    dcm R_nh              = R_nb * R_bh;
    dcm SR_nh             = R_nb * SR_bh;
    Eigen::Matrix3d Jac_R = R_nb.jac_left_composition_wrt_second(R_bh);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarB());
    dcm SSR_nh            = R_nh.plus_left(Deltar_bh_R);

    check("R jac left composition second 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac left composition second 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac left composition second 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac left composition second 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac left composition second 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac left composition second 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac left composition second 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac left composition second 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac left composition second 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-12);

    rotv Sr_bh            = r_bh.plus_left(DeltarB);
    rotv r_nh             = r_nb * r_bh;
    rotv Sr_nh            = r_nb * Sr_bh;
    Eigen::Matrix3d Jac_r = r_nb.jac_left_composition_wrt_second(r_bh);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarB());
    rotv SSr_nh           = r_nh.plus_left(Deltar_bh_r);

    check("r jac left composition second 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-12);
    check("r jac left composition second 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-12);
    check("r jac left composition second 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_left_composition_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_forward_rotation_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_left_forward_rotation_wrt_rotation method, validating the following:
    // (DeltarN plus R) * v = R * v + J * DeltarN

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vB(0.33, -0.65, 0.17);

    Eigen::Vector3d vN_q      = q_nb * vB;
    rodrigues Sq_nb           = q_nb.plus_left(DeltarN);
    Eigen::Vector3d SvN_q     = Sq_nb * vB;
    Eigen::Matrix3d Jac_q     = q_nb.jac_left_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_q = Jac_q * DeltarN();
    Eigen::Vector3d SSvN_q    = vN_q + DeltavN_q;

    check("q jac left forward rotation rotation 1 taylor    ", SvN_q(0), SSvN_q(0), 1e-5);
    check("q jac left forward rotation rotation 2 taylor    ", SvN_q(1), SSvN_q(1), 1e-5);
    check("q jac left forward rotation rotation 3 taylor    ", SvN_q(2), SSvN_q(2), 1e-5);

    Eigen::Vector3d vN_R      = R_nb * vB;
    dcm SR_nb                 = R_nb.plus_left(DeltarN);
    Eigen::Vector3d SvN_R     = SR_nb * vB;
    Eigen::Matrix3d Jac_R     = R_nb.jac_left_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_R = Jac_R * DeltarN();
    Eigen::Vector3d SSvN_R    = vN_R + DeltavN_R;

    check("R jac left forward rotation rotation 1 taylor    ", SvN_R(0), SSvN_R(0), 1e-5);
    check("R jac left forward rotation rotation 2 taylor    ", SvN_R(1), SSvN_R(1), 1e-5);
    check("R jac left forward rotation rotation 3 taylor    ", SvN_R(2), SSvN_R(2), 1e-5);

    Eigen::Vector3d vN_r      = r_nb * vB;
    rotv Sr_nb                = r_nb.plus_left(DeltarN);
    Eigen::Vector3d SvN_r     = Sr_nb * vB;
    Eigen::Matrix3d Jac_r     = r_nb.jac_left_forward_rotation_wrt_rotation(vB);
    Eigen::Vector3d DeltavN_r = Jac_r * DeltarN();
    Eigen::Vector3d SSvN_r    = vN_r + DeltavN_r;

    check("r jac left forward rotation rotation 1 taylor    ", SvN_r(0), SSvN_r(0), 1e-5);
    check("r jac left forward rotation rotation 2 taylor    ", SvN_r(1), SSvN_r(1), 1e-5);
    check("r jac left forward rotation rotation 3 taylor    ", SvN_r(2), SSvN_r(2), 1e-5);
} // closes test_jac_left_forward_rotation_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_backward_rotation_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_left_backward_rotation_wrt_rotation method, validating the following:
    // (DeltarN plus R) / v = R / v + J * DeltarN

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vN(0.33, -0.65, 0.17);

    Eigen::Vector3d vB_q      = q_nb / vN;
    rodrigues Sq_nb           = q_nb.plus_left(DeltarN);
    Eigen::Vector3d SvB_q     = Sq_nb / vN;
    Eigen::Matrix3d Jac_q     = q_nb.jac_left_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_q = Jac_q * DeltarN();
    Eigen::Vector3d SSvB_q    = vB_q + DeltavB_q;

    check("q jac left backward rotation rotation 1 taylor    ", SvB_q(0), SSvB_q(0), 1e-5);
    check("q jac left backward rotation rotation 2 taylor    ", SvB_q(1), SSvB_q(1), 1e-5);
    check("q jac left backward rotation rotation 3 taylor    ", SvB_q(2), SSvB_q(2), 1e-5);

    Eigen::Vector3d vB_R      = R_nb / vN;
    dcm SR_nb                 = R_nb.plus_left(DeltarN);
    Eigen::Vector3d SvB_R     = SR_nb / vN;
    Eigen::Matrix3d Jac_R     = R_nb.jac_left_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_R = Jac_R * DeltarN();
    Eigen::Vector3d SSvB_R    = vB_R + DeltavB_R;

    check("R jac left backward rotation rotation 1 taylor    ", SvB_R(0), SSvB_R(0), 1e-5);
    check("R jac left backward rotation rotation 2 taylor    ", SvB_R(1), SSvB_R(1), 1e-5);
    check("R jac left backward rotation rotation 3 taylor    ", SvB_R(2), SSvB_R(2), 1e-5);

    Eigen::Vector3d vB_r      = r_nb / vN;
    rotv Sr_nb                = r_nb.plus_left(DeltarN);
    Eigen::Vector3d SvB_r     = Sr_nb / vN;
    Eigen::Matrix3d Jac_r     = r_nb.jac_left_backward_rotation_wrt_rotation(vN);
    Eigen::Vector3d DeltavB_r = Jac_r * DeltarN();
    Eigen::Vector3d SSvB_r    = vB_r + DeltavB_r;

    check("r jac left backward rotation rotation 1 taylor    ", SvB_r(0), SSvB_r(0), 1e-5);
    check("r jac left backward rotation rotation 2 taylor    ", SvB_r(1), SSvB_r(1), 1e-5);
    check("r jac left backward rotation rotation 3 taylor    ", SvB_r(2), SSvB_r(2), 1e-5);
} // closes test_jac_left_backward_rotation_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_log() {

    // It checks the proper behavior of the "rotv" left logarithmic jacobian method,
    // validating the following equations:
    // Log(Deltarv plus R) ~= Log(R) + (J * Deltarv)

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);
    rotv DeltarN(1e-3, -5e-4, 2e-3);

    rodrigues Sq_nb       = q_nb.plus_left(DeltarN);
    rotv Sr_nb_q          = Sq_nb.log_map();
    rotv r_nb_q           = q_nb.log_map();
    Eigen::Matrix3d Jac_q = q_nb.jac_left_log();
    rotv Deltar_nb_q      = rotv(Jac_q * DeltarN());
    rotv SSr_nb_q(r_nb_q() + Deltar_nb_q());

    check("q jac left log 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-6);
    check("q jac left log 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-6);
    check("q jac left log 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-6);

    dcm SR_nb             = R_nb.plus_left(DeltarN);
    rotv Sr_nb_R          = SR_nb.log_map();
    rotv r_nb_R           = R_nb.log_map();
    Eigen::Matrix3d Jac_R = R_nb.jac_left_log();
    rotv Deltar_nb_R      = rotv(Jac_R * DeltarN());
    rotv SSr_nb_R(r_nb_R() + Deltar_nb_R());

    check("R jac left log 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-6);
    check("R jac left log 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-6);
    check("R jac left log 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-6);

    rotv Sr_nb_r          = r_nb.plus_left(DeltarN);
    Eigen::Matrix3d Jac_r = r_nb.jac_left_log();
    rotv Deltar_nb_r      = rotv(Jac_r * DeltarN());
    rotv SSr_nb_r(r_nb() + Deltar_nb_r());

    check("r jac left log 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-6);
    check("r jac left log 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-6);
    check("r jac left log 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-6);

} // closes test_jacobian_left_log

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_plus_wrt_first() {

    // It verifies the behavior of the rotation classes jac_left_plus_wrt_first method, validating the following:
    // (rv1 + Deltar1N) plus R2 = J * Deltar1N plus (rv1 plus R2) */

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    rotv r_en(0.4, -1.2, 0.72);
    rotv rr_bh(r_en() + DeltarN());

    rodrigues q_eb        = q_nb.plus_left(r_en);
    rodrigues Sq_nh       = q_nb.plus_left(rr_bh);
    Eigen::Matrix3d Jac_q = q_nb.jac_left_plus_wrt_first(r_en);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarN());
    rodrigues SSq_nh      = q_eb.plus_left(Deltar_bh_q);

    check("q jac left plus first 1 taylor    ", Sq_nh()(0), SSq_nh()(0), 1e-6);
    check("q jac left plus first 2 taylor    ", Sq_nh()(1), SSq_nh()(1), 1e-6);
    check("q jac left plus first 3 taylor    ", Sq_nh()(2), SSq_nh()(2), 1e-6);
    check("q jac left plus first 4 taylor    ", Sq_nh()(3), SSq_nh()(3), 1e-6);

    dcm R_eb              = R_nb.plus_left(r_en);
    dcm SR_nh             = R_nb.plus_left(rr_bh);
    Eigen::Matrix3d Jac_R = R_nb.jac_left_plus_wrt_first(r_en);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarN());
    dcm SSR_nh            = R_eb.plus_left(Deltar_bh_R);

    check("R jac left plus first 11 taylor   ", SR_nh()(0,0), SSR_nh()(0,0), 1e-6);
    check("R jac left plus first 12 taylor   ", SR_nh()(0,1), SSR_nh()(0,1), 1e-6);
    check("R jac left plus first 13 taylor   ", SR_nh()(0,2), SSR_nh()(0,2), 1e-6);
    check("R jac left plus first 21 taylor   ", SR_nh()(1,0), SSR_nh()(1,0), 1e-6);
    check("R jac left plus first 22 taylor   ", SR_nh()(1,1), SSR_nh()(1,1), 1e-6);
    check("R jac left plus first 23 taylor   ", SR_nh()(1,2), SSR_nh()(1,2), 1e-6);
    check("R jac left plus first 31 taylor   ", SR_nh()(2,0), SSR_nh()(2,0), 1e-6);
    check("R jac left plus first 32 taylor   ", SR_nh()(2,1), SSR_nh()(2,1), 1e-6);
    check("R jac left plus first 33 taylor   ", SR_nh()(2,2), SSR_nh()(2,2), 1e-6);

    rotv r_eb             = r_nb.plus_left(r_en);
    rotv Sr_nh            = r_nb.plus_left(rr_bh);
    Eigen::Matrix3d Jac_r = r_nb.jac_left_plus_wrt_first(r_en);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarN());
    rotv SSr_nh           = r_eb.plus_left(Deltar_bh_r);

    check("r jac left plus first 1 taylor    ", Sr_nh()(0), SSr_nh()(0), 1e-6);
    check("r jac left plus first 2 taylor    ", Sr_nh()(1), SSr_nh()(1), 1e-6);
    check("r jac left plus first 3 taylor    ", Sr_nh()(2), SSr_nh()(2), 1e-6);

} // closes test_jac_left_plus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_plus_wrt_second() {

    // It verifies the behavior of the rotation classes jac_left_plus_wrt_second method, validating the following:
    // rv1 plus (Deltar2N plus R2) = J * Deltar2N plus (rv1 plus R2) */

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    rotv r_en(0.4, -1.2, 0.72);

    rodrigues q_eb        = q_nb.plus_left(r_en);
    rodrigues Sq_nb       = q_nb.plus_left(DeltarN);
    rodrigues Sq_eb       = Sq_nb.plus_left(r_en);
    Eigen::Matrix3d Jac_q = q_nb.jac_left_plus_wrt_second(r_en);
    rotv Deltar_bh_q      = ang::rotv(Jac_q * DeltarN());
    rodrigues SSq_nh      = q_eb.plus_left(Deltar_bh_q);

    check("q jac left plus second 1 taylor    ", Sq_eb()(0), SSq_nh()(0), 1e-12);
    check("q jac left plus second 2 taylor    ", Sq_eb()(1), SSq_nh()(1), 1e-12);
    check("q jac left plus second 3 taylor    ", Sq_eb()(2), SSq_nh()(2), 1e-12);
    check("q jac left plus second 4 taylor    ", Sq_eb()(3), SSq_nh()(3), 1e-12);

    dcm R_eb              = R_nb.plus_left(r_en);
    dcm SR_nb             = R_nb.plus_left(DeltarN);
    dcm SR_eb             = SR_nb.plus_left(r_en);
    Eigen::Matrix3d Jac_R = R_nb.jac_left_plus_wrt_second(r_en);
    rotv Deltar_bh_R      = ang::rotv(Jac_R * DeltarN());
    dcm SSR_nh            = R_eb.plus_left(Deltar_bh_R);

    check("R jac left plus second 11 taylor   ", SR_eb()(0,0), SSR_nh()(0,0), 1e-12);
    check("R jac left plus second 12 taylor   ", SR_eb()(0,1), SSR_nh()(0,1), 1e-12);
    check("R jac left plus second 13 taylor   ", SR_eb()(0,2), SSR_nh()(0,2), 1e-12);
    check("R jac left plus second 21 taylor   ", SR_eb()(1,0), SSR_nh()(1,0), 1e-12);
    check("R jac left plus second 22 taylor   ", SR_eb()(1,1), SSR_nh()(1,1), 1e-12);
    check("R jac left plus second 23 taylor   ", SR_eb()(1,2), SSR_nh()(1,2), 1e-12);
    check("R jac left plus second 31 taylor   ", SR_eb()(2,0), SSR_nh()(2,0), 1e-12);
    check("R jac left plus second 32 taylor   ", SR_eb()(2,1), SSR_nh()(2,1), 1e-12);
    check("R jac left plus second 33 taylor   ", SR_eb()(2,2), SSR_nh()(2,2), 1e-12);

    rotv r_eb             = r_nb.plus_left(r_en);
    rotv Sr_nb            = r_nb.plus_left(DeltarN);
    rotv Sr_eb            = Sr_nb.plus_left(r_en);
    Eigen::Matrix3d Jac_r = r_nb.jac_left_plus_wrt_second(r_en);
    rotv Deltar_bh_r      = ang::rotv(Jac_r * DeltarN());
    rotv SSr_nh           = r_eb.plus_left(Deltar_bh_r);

    check("r jac left plus second 1 taylor    ", Sr_eb()(0), SSr_nh()(0), 1e-12);
    check("r jac left plus second 2 taylor    ", Sr_eb()(1), SSr_nh()(1), 1e-12);
    check("r jac left plus second 3 taylor    ", Sr_eb()(2), SSr_nh()(2), 1e-12);

} // closes test_jac_left_plus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_minus_wrt_first() {

    // It verifies the behavior of the rotation classes jac_left_minus_wrt_first method, validating the following:
    // (Deltar2N plus R2) minus R1 = (R2 minus R1) + J * Deltar2N

    rotv r1_nb(-0.7, 0.45, 1.0);
    rodrigues q1_nb(r1_nb);
    dcm R1_nb(r1_nb);

    rotv DeltarNA(1e-3, -5e-4, 2e-3);

    rotv r2_nb      = r1_nb.plus_left(DeltarNA);
    rodrigues q2_nb = q1_nb.plus_left(DeltarNA);
    dcm R2_nb       = R1_nb.plus_left(DeltarNA);

    rotv DeltarN(-3.2e-4, -1.7e-4, 2.2e-3);

    rotv DeltarN_q        = q2_nb.minus_left(q1_nb);
    rodrigues Sq_nb       = q2_nb.plus_left(DeltarN);
    rotv Sr_nb_q          = Sq_nb.minus_left(q1_nb);
    Eigen::Matrix3d Jac_q = q2_nb.jac_left_minus_wrt_first(q1_nb);
    rotv Deltar_nb_q      = ang::rotv(Jac_q * DeltarN());
    rotv SSr_nb_q(DeltarN_q() + Deltar_nb_q());

    check("q jac left minus first 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-7);
    check("q jac left minus first 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-7);
    check("q jac left minus first 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-7);

    rotv DeltarN_R        = R2_nb.minus_left(R1_nb);
    dcm SR_nb             = R2_nb.plus_left(DeltarN);
    rotv Sr_nb_R          = SR_nb.minus_left(R1_nb);
    Eigen::Matrix3d Jac_R = R2_nb.jac_left_minus_wrt_first(R1_nb);
    rotv Deltar_nb_R      = ang::rotv(Jac_R * DeltarN());
    rotv SSr_nb_R(DeltarN_R() + Deltar_nb_R());

    check("R jac left minus first 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-7);
    check("R jac left minus first 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-7);
    check("R jac left minus first 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-7);

    rotv DeltarN_r        = r2_nb.minus_left(r1_nb);
    rotv Sr_nb            = r2_nb.plus_left(DeltarN);
    rotv Sr_nb_r          = Sr_nb.minus_left(r1_nb);
    Eigen::Matrix3d Jac_r = r2_nb.jac_left_minus_wrt_first(r1_nb);
    rotv Deltar_nb_r      = ang::rotv(Jac_r * DeltarN());
    rotv SSr_nb_r(DeltarN_r() + Deltar_nb_r());

    check("r jac left minus first 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-7);
    check("r jac left minus first 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-7);
    check("r jac left minus first 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-7);

} // closes test_jac_left_minus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_minus_wrt_second() {

    // It verifies the behavior of the rotation classes jac_right_minus_wrt_first method, validating the following:
    // R2 minus (Deltar1N plus R1) = (R2 minus R1) + J * Deltar1N */

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarNA(1e-3, -5e-4, 2e-3);

    rotv Pr_nb      = r_nb.plus_left(DeltarNA);
    rodrigues Pq_nb = q_nb.plus_left(DeltarNA);
    dcm PR_nb       = R_nb.plus_left(DeltarNA);

    rotv DeltarN(-3.2e-4, -1.7e-4, 2.2e-3);

    rotv DeltarN_q        = Pq_nb.minus_left(q_nb);
    rodrigues Sq_nb       = q_nb.plus_left(DeltarN);
    rotv Sr_nb_q          = Pq_nb.minus_left(Sq_nb);
    Eigen::Matrix3d Jac_q = Pq_nb.jac_left_minus_wrt_second(q_nb);
    rotv Deltar_nb_q      = ang::rotv(Jac_q * DeltarN());
    rotv SSr_nb_q(DeltarN_q() + Deltar_nb_q());

    check("q jac left minus second 1 taylor    ", Sr_nb_q()(0), SSr_nb_q()(0), 1e-9);
    check("q jac left minus second 2 taylor    ", Sr_nb_q()(1), SSr_nb_q()(1), 1e-9);
    check("q jac left minus second 3 taylor    ", Sr_nb_q()(2), SSr_nb_q()(2), 1e-9);

    rotv DeltarN_R        = PR_nb.minus_left(R_nb);
    dcm SR_nb             = R_nb.plus_left(DeltarN);
    rotv Sr_nb_R          = PR_nb.minus_left(SR_nb);
    Eigen::Matrix3d Jac_R = PR_nb.jac_left_minus_wrt_second(R_nb);
    rotv Deltar_nb_R      = ang::rotv(Jac_R * DeltarN());
    rotv SSr_nb_R(DeltarN_R() + Deltar_nb_R());

    check("R jac left minus second 1 taylor    ", Sr_nb_R()(0), SSr_nb_R()(0), 1e-9);
    check("R jac left minus second 2 taylor    ", Sr_nb_R()(1), SSr_nb_R()(1), 1e-9);
    check("R jac left minus second 3 taylor    ", Sr_nb_R()(2), SSr_nb_R()(2), 1e-9);

    rotv DeltarN_r        = Pr_nb.minus_left(r_nb);
    rotv Sr_nb            = r_nb.plus_left(DeltarN);
    rotv Sr_nb_r          = Pr_nb.minus_left(Sr_nb);
    Eigen::Matrix3d Jac_r = Pr_nb.jac_left_minus_wrt_second(r_nb);
    rotv Deltar_nb_r      = ang::rotv(Jac_r * DeltarN());
    rotv SSr_nb_r(DeltarN_r() + Deltar_nb_r());

    check("r jac left minus second 1 taylor    ", Sr_nb_r()(0), SSr_nb_r()(0), 1e-9);
    check("r jac left minus second 2 taylor    ", Sr_nb_r()(1), SSr_nb_r()(1), 1e-9);
    check("r jac left minus second 3 taylor    ", Sr_nb_r()(2), SSr_nb_r()(2), 1e-9);

} // closes test_jac_left_minus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_forward_adjoint_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_left_forward_adjoint_wrt_rotation method, validating the following:
    // Ad(DeltarN plus R) | w = AdR | w + J * DeltarN

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    so3_tangent wB(0.33, -0.65, 0.17);

    so3_tangent wN_q          = q_nb | wB;
    rodrigues Sq_nb           = q_nb.plus_left(DeltarN);
    so3_tangent SwN_q         = Sq_nb | wB;
    Eigen::Matrix3d Jac_q     = q_nb.jac_left_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_q = Jac_q * DeltarN();
    so3_tangent SSwN_q        = wN_q + DeltawN_q;

    check("q jac left forward adjoint rotation 1 taylor    ", SwN_q()(0), SSwN_q()(0), 1e-5);
    check("q jac left forward adjoint rotation 2 taylor    ", SwN_q()(1), SSwN_q()(1), 1e-5);
    check("q jac left forward adjoint rotation 3 taylor    ", SwN_q()(2), SSwN_q()(2), 1e-5);

    so3_tangent wN_R          = R_nb | wB;
    dcm SR_nb                 = R_nb.plus_left(DeltarN);
    so3_tangent SwN_R         = SR_nb | wB;
    Eigen::Matrix3d Jac_R     = R_nb.jac_left_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_R = Jac_R * DeltarN();
    so3_tangent SSwN_R        = wN_R + DeltawN_R;

    check("R jac left forward adjoint rotation 1 taylor    ", SwN_R()(0), SSwN_R()(0), 1e-5);
    check("R jac left forward adjoint rotation 2 taylor    ", SwN_R()(1), SSwN_R()(1), 1e-5);
    check("R jac left forward adjoint rotation 3 taylor    ", SwN_R()(2), SSwN_R()(2), 1e-5);

    so3_tangent wN_r          = r_nb | wB;
    rotv Sr_nb                = r_nb.plus_left(DeltarN);
    so3_tangent SwN_r         = Sr_nb | wB;
    Eigen::Matrix3d Jac_r     = r_nb.jac_left_forward_adjoint_wrt_rotation(wB);
    Eigen::Vector3d DeltawN_r = Jac_r * DeltarN();
    so3_tangent SSwN_r        = wN_r + DeltawN_r;

    check("r jac left forward adjoint rotation 1 taylor    ", SwN_r()(0), SSwN_r()(0), 1e-5);
    check("r jac left forward adjoint rotation 2 taylor    ", SwN_r()(1), SSwN_r()(1), 1e-5);
    check("r jac left forward adjoint rotation 3 taylor    ", SwN_r()(2), SSwN_r()(2), 1e-5);
} // closes test_jac_left_forward_adjoint_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_left_backward_adjoint_wrt_rotation() {

    // It verifies the behavior of the rotation classes jac_left_backward_adjoint_wrt_rotation method, validating the following:
    // Ad(DeltarN plus R) % w = AdR % w + J * DeltarN

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    rotv DeltarN(1e-3, -5e-4, 2e-3);

    so3_tangent wN(0.33, -0.65, 0.17);

    so3_tangent wB_q          = q_nb % wN;
    rodrigues Sq_nb           = q_nb.plus_left(DeltarN);
    so3_tangent SwB_q         = Sq_nb % wN;
    Eigen::Matrix3d Jac_q     = q_nb.jac_left_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_q = Jac_q * DeltarN();
    so3_tangent SSwB_q        = wB_q + DeltawB_q;

    check("q jac left backward adjoint rotation 1 taylor    ", SwB_q()(0), SSwB_q()(0), 1e-5);
    check("q jac left backward adjoint rotation 2 taylor    ", SwB_q()(1), SSwB_q()(1), 1e-5);
    check("q jac left backward adjoint rotation 3 taylor    ", SwB_q()(2), SSwB_q()(2), 1e-5);

    so3_tangent wB_R          = R_nb % wN;
    dcm SR_nb                 = R_nb.plus_left(DeltarN);
    so3_tangent SwB_R         = SR_nb % wN;
    Eigen::Matrix3d Jac_R     = R_nb.jac_left_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_R = Jac_R * DeltarN();
    so3_tangent SSwB_R        = wB_R + DeltawB_R;

    check("R jac left backward adjoint rotation 1 taylor    ", SwB_R()(0), SSwB_R()(0), 1e-5);
    check("R jac left backward adjoint rotation 2 taylor    ", SwB_R()(1), SSwB_R()(1), 1e-5);
    check("R jac left backward adjoint rotation 3 taylor    ", SwB_R()(2), SSwB_R()(2), 1e-5);

    so3_tangent wB_r          = r_nb % wN;
    rotv Sr_nb                = r_nb.plus_left(DeltarN);
    so3_tangent SwB_r         = Sr_nb % wN;
    Eigen::Matrix3d Jac_r     = r_nb.jac_left_backward_adjoint_wrt_rotation(wN);
    Eigen::Vector3d DeltawB_r = Jac_r * DeltarN();
    so3_tangent SSwB_r        = wB_r + DeltawB_r;

    check("r jac left backward adjoint rotation 1 taylor    ", SwB_r()(0), SSwB_r()(0), 1e-5);
    check("r jac left backward adjoint rotation 2 taylor    ", SwB_r()(1), SSwB_r()(1), 1e-5);
    check("r jac left backward adjoint rotation 3 taylor    ", SwB_r()(2), SSwB_r()(2), 1e-5);
} // closes test_jac_left_backward_adjoint_wrt_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_forward_rotation_wrt_vector() {

    // It verifies the behavior of the rotation classes jac_euclidean_forward_rotation_wrt_vector method, validating the following:
    // R * v = J * v
    // R * (v + Deltav) = R * v + J * Deltav

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    Eigen::Vector3d vB(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltavB(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vN_q      = q_nb * vB;
    Eigen::Vector3d SvB_q     = vB + DeltavB;
    Eigen::Vector3d SvN_q     = q_nb * SvB_q;
    Eigen::Matrix3d Jac_q     = q_nb.jac_euclidean_forward_rotation_wrt_vector();
    Eigen::Vector3d DeltavN_q = Jac_q * DeltavB;
    Eigen::Vector3d SSvN_q    = vN_q + DeltavN_q;

    check("q jac euclidean forward rotation motion 1 taylor    ", SvN_q(0), SSvN_q(0), 1e-12);
    check("q jac euclidean forward rotation motion 2 taylor    ", SvN_q(1), SSvN_q(1), 1e-12);
    check("q jac euclidean forward rotation motion 3 taylor    ", SvN_q(2), SSvN_q(2), 1e-12);

    Eigen::Vector3d vN_R      = R_nb * vB;
    Eigen::Vector3d SvB_R     = vB + DeltavB;
    Eigen::Vector3d SvN_R     = R_nb * SvB_R;
    Eigen::Matrix3d Jac_R     = R_nb.jac_euclidean_forward_rotation_wrt_vector();
    Eigen::Vector3d DeltavN_R = Jac_R * DeltavB;
    Eigen::Vector3d SSvN_R    = vN_R + DeltavN_R;

    check("R jac euclidean forward rotation motion 1 taylor    ", SvN_R(0), SSvN_R(0), 1e-12);
    check("R jac euclidean forward rotation motion 2 taylor    ", SvN_R(1), SSvN_R(1), 1e-12);
    check("R jac euclidean forward rotation motion 3 taylor    ", SvN_R(2), SSvN_R(2), 1e-12);

    Eigen::Vector3d vN_r      = r_nb * vB;
    Eigen::Vector3d SvB_r     = vB + DeltavB;
    Eigen::Vector3d SvN_r     = r_nb * SvB_r;
    Eigen::Matrix3d Jac_r     = r_nb.jac_euclidean_forward_rotation_wrt_vector();
    Eigen::Vector3d DeltavN_r = Jac_r * DeltavB;
    Eigen::Vector3d SSvN_r    = vN_r + DeltavN_r;

    check("r jac euclidean forward rotation motion 1 taylor    ", SvN_r(0), SSvN_r(0), 1e-12);
    check("r jac euclidean forward rotation motion 2 taylor    ", SvN_r(1), SSvN_r(1), 1e-12);
    check("r jac euclidean forward rotation motion 3 taylor    ", SvN_r(2), SSvN_r(2), 1e-12);

    Eigen::Vector3d vN_linear_q = q_nb.jac_euclidean_forward_rotation_wrt_vector() * vB;
    Eigen::Vector3d vN_linear_R = R_nb.jac_euclidean_forward_rotation_wrt_vector() * vB;
    Eigen::Vector3d vN_linear_r = r_nb.jac_euclidean_forward_rotation_wrt_vector() * vB;

    check("q jac euclidean forward rotation motion 1 linear    ", vN_q(0), vN_linear_q(0), 1e-12);
    check("q jac euclidean forward rotation motion 2 linear    ", vN_q(1), vN_linear_q(1), 1e-12);
    check("q jac euclidean forward rotation motion 3 linear    ", vN_q(2), vN_linear_q(2), 1e-12);
    check("R jac euclidean forward rotation motion 1 linear    ", vN_R(0), vN_linear_R(0), 1e-12);
    check("R jac euclidean forward rotation motion 2 linear    ", vN_R(1), vN_linear_R(1), 1e-12);
    check("R jac euclidean forward rotation motion 3 linear    ", vN_R(2), vN_linear_R(2), 1e-12);
    check("r jac euclidean forward rotation motion 1 linear    ", vN_r(0), vN_linear_r(0), 1e-12);
    check("r jac euclidean forward rotation motion 2 linear    ", vN_r(1), vN_linear_r(1), 1e-12);
    check("r jac euclidean forward rotation motion 3 linear    ", vN_r(2), vN_linear_r(2), 1e-12);

} // closes test_jac_euclidean_forward_rotation_wrt_vector

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_backward_rotation_wrt_vector() {

    // It verifies the behavior of the rotation classes jac_euclidean_backward_rotation_wrt_vector method, validating the following:
    // R / v = J * v
    // R / (v + Deltav) = R / v + J * Deltav

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    Eigen::Vector3d vN(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltavN(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d vB_q      = q_nb / vN;
    Eigen::Vector3d SvN_q     = vN + DeltavN;
    Eigen::Vector3d SvB_q     = q_nb / SvN_q;
    Eigen::Matrix3d Jac_q     = q_nb.jac_euclidean_backward_rotation_wrt_vector();
    Eigen::Vector3d DeltavB_q = Jac_q * DeltavN;
    Eigen::Vector3d SSvB_q    = vB_q + DeltavB_q;

    check("q jac euclidean backward rotation motion 1 taylor    ", SvB_q(0), SSvB_q(0), 1e-12);
    check("q jac euclidean backward rotation motion 2 taylor    ", SvB_q(1), SSvB_q(1), 1e-12);
    check("q jac euclidean backward rotation motion 3 taylor    ", SvB_q(2), SSvB_q(2), 1e-12);

    Eigen::Vector3d vB_R      = R_nb / vN;
    Eigen::Vector3d SvN_R     = vN + DeltavN;
    Eigen::Vector3d SvB_R     = R_nb / SvN_R;
    Eigen::Matrix3d Jac_R     = R_nb.jac_euclidean_backward_rotation_wrt_vector();
    Eigen::Vector3d DeltavB_R = Jac_R * DeltavN;
    Eigen::Vector3d SSvB_R    = vB_R + DeltavB_R;

    check("R jac euclidean backward rotation motion 1 taylor    ", SvB_R(0), SSvB_R(0), 1e-12);
    check("R jac euclidean backward rotation motion 2 taylor    ", SvB_R(1), SSvB_R(1), 1e-12);
    check("R jac euclidean backward rotation motion 3 taylor    ", SvB_R(2), SSvB_R(2), 1e-12);

    Eigen::Vector3d vB_r      = r_nb / vN;
    Eigen::Vector3d SvN_r     = vN + DeltavN;
    Eigen::Vector3d SvB_r     = r_nb / SvN_r;
    Eigen::Matrix3d Jac_r     = r_nb.jac_euclidean_backward_rotation_wrt_vector();
    Eigen::Vector3d DeltavB_r = Jac_r * DeltavN;
    Eigen::Vector3d SSvB_r    = vB_r + DeltavB_r;

    check("r jac euclidean backward rotation motion 1 taylor    ", SvB_r(0), SSvB_r(0), 1e-12);
    check("r jac euclidean backward rotation motion 2 taylor    ", SvB_r(1), SSvB_r(1), 1e-12);
    check("r jac euclidean backward rotation motion 3 taylor    ", SvB_r(2), SSvB_r(2), 1e-12);

    Eigen::Vector3d vB_linear_q = q_nb.jac_euclidean_backward_rotation_wrt_vector() * vN;
    Eigen::Vector3d vB_linear_R = R_nb.jac_euclidean_backward_rotation_wrt_vector() * vN;
    Eigen::Vector3d vB_linear_r = r_nb.jac_euclidean_backward_rotation_wrt_vector() * vN;

    check("q jac euclidean backward rotation motion 1 linear    ", vB_q(0), vB_linear_q(0), 1e-12);
    check("q jac euclidean backward rotation motion 2 linear    ", vB_q(1), vB_linear_q(1), 1e-12);
    check("q jac euclidean backward rotation motion 3 linear    ", vB_q(2), vB_linear_q(2), 1e-12);
    check("R jac euclidean backward rotation motion 1 linear    ", vB_R(0), vB_linear_R(0), 1e-12);
    check("R jac euclidean backward rotation motion 2 linear    ", vB_R(1), vB_linear_R(1), 1e-12);
    check("R jac euclidean backward rotation motion 3 linear    ", vB_R(2), vB_linear_R(2), 1e-12);
    check("r jac euclidean backward rotation motion 1 linear    ", vB_r(0), vB_linear_r(0), 1e-12);
    check("r jac euclidean backward rotation motion 2 linear    ", vB_r(1), vB_linear_r(1), 1e-12);
    check("r jac euclidean backward rotation motion 3 linear    ", vB_r(2), vB_linear_r(2), 1e-12);

} // closes test_jac_euclidean_backward_rotation_wrt_vector

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_forward_adjoint_wrt_tangent() {

    // It verifies the behavior of the rotation classes jac_euclidean_forward_adjoint_wrt_tangent method, validating the following:
    // AdR | w = J * w
    // AdR | (w + Deltaw) = AdR | w + J * Deltaw

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    so3_tangent wB(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltawB(1e-3, -5e-4, 2e-3);

    so3_tangent wN_q          = q_nb | wB;
    so3_tangent SwB_q         = wB + DeltawB;
    so3_tangent SwN_q         = q_nb | SwB_q;
    Eigen::Matrix3d Jac_q     = q_nb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawN_q = Jac_q * DeltawB;
    so3_tangent SSwN_q        = wN_q + DeltawN_q;

    check("q jac euclidean forward adjoint tangent 1 taylor    ", SwN_q()(0), SSwN_q()(0), 1e-12);
    check("q jac euclidean forward adjoint tangent 2 taylor    ", SwN_q()(1), SSwN_q()(1), 1e-12);
    check("q jac euclidean forward adjoint tangent 3 taylor    ", SwN_q()(2), SSwN_q()(2), 1e-12);

    so3_tangent wN_R          = R_nb | wB;
    so3_tangent SwB_R         = wB + DeltawB;
    so3_tangent SwN_R         = R_nb | SwB_R;
    Eigen::Matrix3d Jac_R     = R_nb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawN_R = Jac_R * DeltawB;
    so3_tangent SSwN_R        = wN_R + DeltawN_R;

    check("R jac euclidean forward adjoint tangent 1 taylor    ", SwN_R()(0), SSwN_R()(0), 1e-12);
    check("R jac euclidean forward adjoint tangent 2 taylor    ", SwN_R()(1), SSwN_R()(1), 1e-12);
    check("R jac euclidean forward adjoint tangent 3 taylor    ", SwN_R()(2), SSwN_R()(2), 1e-12);

    so3_tangent wN_r          = r_nb | wB;
    so3_tangent SwB_r         = wB + DeltawB;
    so3_tangent SwN_r         = r_nb | SwB_r;
    Eigen::Matrix3d Jac_r     = r_nb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawN_r = Jac_r * DeltawB;
    so3_tangent SSwN_r        = wN_r + DeltawN_r;

    check("r jac euclidean forward adjoint tangent 1 taylor    ", SwN_r()(0), SSwN_r()(0), 1e-12);
    check("r jac euclidean forward adjoint tangent 2 taylor    ", SwN_r()(1), SSwN_r()(1), 1e-12);
    check("r jac euclidean forward adjoint tangent 3 taylor    ", SwN_r()(2), SSwN_r()(2), 1e-12);

    Eigen::Vector3d wN_linear_q = q_nb.jac_euclidean_forward_adjoint_wrt_tangent() * wB();
    Eigen::Vector3d wN_linear_R = R_nb.jac_euclidean_forward_adjoint_wrt_tangent() * wB();
    Eigen::Vector3d wN_linear_r = r_nb.jac_euclidean_forward_adjoint_wrt_tangent() * wB();

    check("q jac euclidean forward adjoint tangent 1 linear    ", wN_q()(0), wN_linear_q(0), 1e-12);
    check("q jac euclidean forward adjoint tangent 2 linear    ", wN_q()(1), wN_linear_q(1), 1e-12);
    check("q jac euclidean forward adjoint tangent 3 linear    ", wN_q()(2), wN_linear_q(2), 1e-12);
    check("R jac euclidean forward adjoint tangent 1 linear    ", wN_R()(0), wN_linear_R(0), 1e-12);
    check("R jac euclidean forward adjoint tangent 2 linear    ", wN_R()(1), wN_linear_R(1), 1e-12);
    check("R jac euclidean forward adjoint tangent 3 linear    ", wN_R()(2), wN_linear_R(2), 1e-12);
    check("r jac euclidean forward adjoint tangent 1 linear    ", wN_r()(0), wN_linear_r(0), 1e-12);
    check("r jac euclidean forward adjoint tangent 2 linear    ", wN_r()(1), wN_linear_r(1), 1e-12);
    check("r jac euclidean forward adjoint tangent 3 linear    ", wN_r()(2), wN_linear_r(2), 1e-12);

} // closes test_jac_euclidean_forward_adjoint_wrt_tangent

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_backward_adjoint_wrt_tangent() {

    // It verifies the behavior of the rotation classes jac_euclidean_backward_adjoint_wrt_tangent method, validating the following:
    // AdR % w = J * w
    // AdR % (w + Deltaw) = AdR % w + J * Deltaw

    rotv r_nb(-0.7, 0.45, 1.0);
    rodrigues q_nb(r_nb);
    dcm R_nb(r_nb);

    so3_tangent wN(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltawN(1e-3, -5e-4, 2e-3);

    so3_tangent wB_q          = q_nb % wN;
    so3_tangent SwN_q         = wN + DeltawN;
    so3_tangent SwB_q         = q_nb % SwN_q;
    Eigen::Matrix3d Jac_q     = q_nb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawB_q = Jac_q * DeltawN;
    so3_tangent SSwB_q        = wB_q + DeltawB_q;

    check("q jac euclidean backward adjoint tangent 1 taylor    ", SwB_q()(0), SSwB_q()(0), 1e-12);
    check("q jac euclidean backward adjoint tangent 2 taylor    ", SwB_q()(1), SSwB_q()(1), 1e-12);
    check("q jac euclidean backward adjoint tangent 3 taylor    ", SwB_q()(2), SSwB_q()(2), 1e-12);

    so3_tangent wB_R          = R_nb % wN;
    so3_tangent SwN_R         = wN + DeltawN;
    so3_tangent SwB_R         = R_nb % SwN_R;
    Eigen::Matrix3d Jac_R     = R_nb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawB_R = Jac_R * DeltawN;
    so3_tangent SSwB_R        = wB_R + DeltawB_R;

    check("R jac euclidean backward adjoint tangent 1 taylor    ", SwB_R()(0), SSwB_R()(0), 1e-12);
    check("R jac euclidean backward adjoint tangent 2 taylor    ", SwB_R()(1), SSwB_R()(1), 1e-12);
    check("R jac euclidean backward adjoint tangent 3 taylor    ", SwB_R()(2), SSwB_R()(2), 1e-12);

    so3_tangent wB_r          = r_nb % wN;
    so3_tangent SwN_r         = wN + DeltawN;
    so3_tangent SwB_r         = r_nb % SwN_r;
    Eigen::Matrix3d Jac_r     = r_nb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector3d DeltawB_r = Jac_r * DeltawN;
    so3_tangent SSwB_r        = wB_r + DeltawB_r;

    check("r jac euclidean backward adjoint tangent 1 taylor    ", SwB_r()(0), SSwB_r()(0), 1e-12);
    check("r jac euclidean backward adjoint tangent 2 taylor    ", SwB_r()(1), SSwB_r()(1), 1e-12);
    check("r jac euclidean backward adjoint tangent 3 taylor    ", SwB_r()(2), SSwB_r()(2), 1e-12);

    Eigen::Vector3d wB_linear_q = q_nb.jac_euclidean_backward_adjoint_wrt_tangent() * wN();
    Eigen::Vector3d wB_linear_R = R_nb.jac_euclidean_backward_adjoint_wrt_tangent() * wN();
    Eigen::Vector3d wB_linear_r = r_nb.jac_euclidean_backward_adjoint_wrt_tangent() * wN();

    check("q jac euclidean backward adjoint tangent 1 linear    ", wB_q()(0), wB_linear_q(0), 1e-12);
    check("q jac euclidean backward adjoint tangent 2 linear    ", wB_q()(1), wB_linear_q(1), 1e-12);
    check("q jac euclidean backward adjoint tangent 3 linear    ", wB_q()(2), wB_linear_q(2), 1e-12);
    check("R jac euclidean backward adjoint tangent 1 linear    ", wB_R()(0), wB_linear_R(0), 1e-12);
    check("R jac euclidean backward adjoint tangent 2 linear    ", wB_R()(1), wB_linear_R(1), 1e-12);
    check("R jac euclidean backward adjoint tangent 3 linear    ", wB_R()(2), wB_linear_R(2), 1e-12);
    check("r jac euclidean backward adjoint tangent 1 linear    ", wB_r()(0), wB_linear_r(0), 1e-12);
    check("r jac euclidean backward adjoint tangent 2 linear    ", wB_r()(1), wB_linear_r(1), 1e-12);
    check("r jac euclidean backward adjoint tangent 3 linear    ", wB_r()(2), wB_linear_r(2), 1e-12);

} // closes test_jac_euclidean_backward_adjoint_wrt_tangent

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_wrt_dcm() {

    // It verifies the behavior of the dcm class jac_forward_rotation_wrt_dcm and
    // jac_backward_rotation_wrt_dcm, validating the following equations:
    // [R * exp(Delta r)] * v ~= R * v + d(R * v)/dq |R*v * {[R * exp(Delta r)] - R}
    // [exp(Delta r) * R] * v ~= R * v + d(R * v)/dq |R*v * {[exp(Delta r) * R] - R}
    // [R * exp(Delta r)] / v ~= R / v + d(R / v)/dq |R/v * {[q * exp(Delta r)] - R}
    // [exp(Delta r) * R] / v ~= R / v + d(R / v)/dq |R/v * {[exp(Delta r) * R] - R}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta r

    double d2r = math::constant::D2R();

    euler euler_nb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    dcm R_nb(euler_nb);
    Eigen::Vector9d R_nb_wedge = dcm::wedge(R_nb);

    rotv Deltar_nb(0.001, -0.002, 0.0015); // WATCH OUT - needs to be small
    dcm DeltaR_nb(Deltar_nb);
    Eigen::Vector9d DeltaR_nb_wedge = dcm::wedge(DeltaR_nb);

    Eigen::Vector3d v_b(3.8, -4.1, -2.3);
    Eigen::Vector3d v_n = R_nb * v_b;

    dcm              SR_nb       = R_nb.plus_right(Deltar_nb); // WATCH OUT - same as R_nb * DeltaR_nb
    Eigen::Vector3d  Sv_n        = SR_nb * v_b;
    Eigen::Matrix39d Jac_f_r     = R_nb.jac_euclidean_forward_rotation_wrt_dcm(v_b);
    Eigen::Vector3d  SDelta_v_n  = Jac_f_r * dcm::wedge(SR_nb() - R_nb()); // WATCH OUT - very different from Jac_f_r * DeltaR_nb().get()
    Eigen::Vector3d  SSv_n       = v_n + SDelta_v_n;

    check("q jac forward dcm rotate 0 taylor     ", Sv_n(0), SSv_n(0), 1e-5);
    check("q jac forward dcm rotate 1 taylor     ", Sv_n(1), SSv_n(1), 1e-6);
    check("q jac forward dcm rotate 2 taylor     ", Sv_n(2), SSv_n(2), 1e-4);

    dcm             TR_nb        = R_nb.plus_left(Deltar_nb);
    Eigen::Vector3d Tv_n         = TR_nb * v_b;
    Eigen::Vector3d TDeltaB_v_n  = Jac_f_r * dcm::wedge(TR_nb() - R_nb());
    Eigen::Vector3d TTv_n        = v_n + TDeltaB_v_n;

    check("q jac forward dcm rotate 0 taylor     ", Tv_n(0), TTv_n(0), 1e-5);
    check("q jac forward dcm rotate 1 taylor     ", Tv_n(1), TTv_n(1), 1e-6);
    check("q jac forward dcm rotate 2 taylor     ", Tv_n(2), TTv_n(2), 1e-4);

    Eigen::Vector3d  Sv_b         = SR_nb / v_n;
    Eigen::Matrix39d Jac_b_r      = R_nb.jac_euclidean_backward_rotation_wrt_dcm(v_n);
    Eigen::Vector3d  SDelta_v_b   = Jac_b_r * dcm::wedge(SR_nb() - R_nb());
    Eigen::Vector3d  SSv_b        = v_b + SDelta_v_b;

    check("q jac backward dcm rotate 0 taylor    ", Sv_b(0), SSv_b(0), 1e-5);
    check("q jac backward dcm rotate 1 taylor    ", Sv_b(1), SSv_b(1), 1e-5);
    check("q jac backward dcm rotate 2 taylor    ", Sv_b(2), SSv_b(2), 1e-4);

    Eigen::Vector3d Tv_b            = TR_nb / v_n;
    Eigen::Vector3d TDelta_v_b      = Jac_b_r * dcm::wedge(TR_nb() - R_nb());
    Eigen::Vector3d TTv_b           = v_b + TDelta_v_b;

    check("q jac backward dcm rotate 0 taylor    ", Tv_b(0), TTv_b(0), 1e-5);
    check("q jac backward dcm rotate 1 taylor    ", Tv_b(1), TTv_b(1), 1e-5);
    check("q jac backward dcm rotate 2 taylor    ", Tv_b(2), TTv_b(2), 1e-4);

    // It verifies the behavior of the dcm class jac_forward_adjoint_wrt_dcm and
    // jac_backward_adjoint_wrt_dcm, validating the following equations:
    // [R * exp(Delta r)] | w ~= R | w + d(R | w)/dR |R|w * {[R * exp(Delta r)] - R}
    // [exp(Delta r) * R] | w ~= R | w + d(R | w)/dR |R|w * {[exp(Delta r) * R] - R}
    // [R * exp(Delta r)] % w ~= R % w + d(R % w)/dR |R%w * {[R * exp(Delta r)] - R}
    // [exp(Delta r) * R] % w ~= R % w + d(R % w)/dR |R%w * {[exp(Delta r) * R] - R}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta r

    ang::so3_tangent w_nbb(3.8, -4.1, -2.3);
    ang::so3_tangent w_nbn = R_nb | w_nbb;

    ang::so3_tangent Sw_nbn         = SR_nb | w_nbb;
    Eigen::Matrix39d Jac_f_a        = R_nb.jac_euclidean_forward_adjoint_wrt_dcm(w_nbb);
    ang::so3_tangent SDeltaA_w_nbn(Jac_f_a * dcm::wedge(SR_nb() - R_nb()));
    ang::so3_tangent SSw_nbn        = w_nbn + SDeltaA_w_nbn;

    check("q jac forward dcm adjoint 0 taylor    ", Sw_nbn()(0), SSw_nbn()(0), 1e-5);
    check("q jac forward dcm adjoint 1 taylor    ", Sw_nbn()(1), SSw_nbn()(1), 1e-6);
    check("q jac forward dcm adjoint 2 taylor    ", Sw_nbn()(2), SSw_nbn()(2), 1e-4);

    ang::so3_tangent Tw_nbn      = TR_nb | w_nbb;
    ang::so3_tangent TDelta_w_nbn(Jac_f_a * dcm::wedge(TR_nb() - R_nb()));
    ang::so3_tangent TTw_nbn     = w_nbn + TDelta_w_nbn;

    check("q jac forward dcm adjoint 0 taylor    ", Tw_nbn()(0), TTw_nbn()(0), 1e-5);
    check("q jac forward dcm adjoint 1 taylor    ", Tw_nbn()(1), TTw_nbn()(1), 1e-6);
    check("q jac forward dcm adjoint 2 taylor    ", Tw_nbn()(2), TTw_nbn()(2), 1e-4);

    ang::so3_tangent Sw_nbb         = SR_nb % w_nbn;
    Eigen::Matrix39d Jac_b_a        = R_nb.jac_euclidean_backward_adjoint_wrt_dcm(w_nbn);
    ang::so3_tangent SDelta_w_nbb(Jac_b_a * dcm::wedge(SR_nb() - R_nb()));
    ang::so3_tangent SSw_nbb        = w_nbb + SDelta_w_nbb;

    check("q jac backward dcm adjoint 0 taylor   ", Sw_nbb()(0), SSw_nbb()(0), 1e-5);
    check("q jac backward dcm adjoint 1 taylor   ", Sw_nbb()(1), SSw_nbb()(1), 1e-5);
    check("q jac backward dcm adjoint 2 taylor   ", Sw_nbb()(2), SSw_nbb()(2), 1e-4);

    ang::so3_tangent Tw_nbb         = TR_nb % w_nbn;
    ang::so3_tangent TDelta_w_nbb(Jac_b_a * dcm::wedge(TR_nb() - R_nb()));
    ang::so3_tangent TTw_nbb        = w_nbb + TDelta_w_nbb;

    check("q jac backward dcm adjoint 0 taylor   ", Tw_nbb()(0), TTw_nbb()(0), 1e-5);
    check("q jac backward dcm adjoint 1 taylor   ", Tw_nbb()(1), TTw_nbb()(1), 1e-5);
    check("q jac backward dcm adjoint 2 taylor   ", Tw_nbb()(2), TTw_nbb()(2), 1e-4);

} // closes test_jac_euclidean_wrt_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_wrt_dcm_linear() {

    // It verifies the behavior of the rotation matrix class jac_euclidean_forward_rotation_wrt_dcm and
    // jac_euclidean_backward_rotation_wrt_dcm, validating the following equations:

    // R * v = d(R * v)/dR |R*v * Rwedge
    // (R + DeltaR) * v = R * v + d(R * v)/dR |R*v * DeltaRwedge
    // R / v = d(R / v)/dR |R/v * Rwedge
    // (R + DeltaR) / v = R / v + d(R / v)/dR |R/v * DeltaRwedge

    double d2r = math::constant::D2R();

    euler euler_nb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    dcm R_nb(euler_nb);
    dcm R_bn = R_nb.inverse();
    Eigen::Vector9d R_nb_wedge = dcm::wedge(R_nb);

    Eigen::Vector3d v_b(3.8, -4.1, -2.3);
    Eigen::Vector3d v_n = R_nb * v_b;

    euler Deltaeuler_nb(3.2 * d2r, 8.7 * d2r, -4.5 * d2r); // does NOT need to be small
    dcm DeltaR_nb(Deltaeuler_nb);
    Eigen::Vector9d DeltaR_nb_wedge = dcm::wedge(DeltaR_nb);

    Eigen::Matrix39d J_forward_rotation = R_nb.jac_euclidean_forward_rotation_wrt_dcm(v_b);
    Eigen::Vector3d v_n2                = J_forward_rotation * R_nb_wedge;

    check("R jac forward dcm rotate 0 linear     ", v_n(0), v_n2(0), 1e-14);
    check("R jac forward dcm rotate 1 linear     ", v_n(1), v_n2(1), 1e-14);
    check("R jac forward dcm rotate 2 linear     ", v_n(2), v_n2(2), 1e-14);

    Eigen::Matrix3d SR_nb = R_nb() + DeltaR_nb(); // WATCH OUT - This is a matrix, not a dcm, as otherwise it normalizes
    Eigen::Vector3d Sv_n  = SR_nb * v_b;
    Eigen::Vector3d SSv_n = v_n + J_forward_rotation * DeltaR_nb_wedge;

    check("R jac forward dcm rotate 0 linear     ", Sv_n(0), SSv_n(0), 1e-14);
    check("R jac forward dcm rotate 1 linear     ", Sv_n(1), SSv_n(1), 1e-14);
    check("R jac forward dcm rotate 2 linear     ", Sv_n(2), SSv_n(2), 1e-14);

    Eigen::Matrix39d J_backward_rotation = R_nb.jac_euclidean_backward_rotation_wrt_dcm(v_n);
    Eigen::Vector3d v_b2                 = J_backward_rotation * R_nb_wedge;

    check("R jac backward dcm rotate 0 linear    ", v_b(0), v_b2(0), 1e-14);
    check("R jac backward dcm rotate 1 linear    ", v_b(1), v_b2(1), 1e-14);
    check("R jac backward dcm rotate 2 linear    ", v_b(2), v_b2(2), 1e-14);

    Eigen::Vector3d Sv_b  = SR_nb.transpose() *  v_n; // WATCH OUT - This can not be the inverse rotation operator /
    Eigen::Vector3d SSv_b = v_b + J_backward_rotation * DeltaR_nb_wedge;

    check("R jac backward dcm rotate 0 linear    ", Sv_b(0), SSv_b(0), 1e-14);
    check("R jac backward dcm rotate 1 linear    ", Sv_b(1), SSv_b(1), 1e-14);
    check("R jac backward dcm rotate 2 linear    ", Sv_b(2), SSv_b(2), 1e-14);

    // It verifies the behavior of the rotation matrix class jac_euclidean_forward_adjoint_wrt_dcm and
    // jac_euclidean_backward_adjoint_wrt_dcm, validating the following equations:

    // R | w = d(R | w)/dR |R|w * Rwedge
    // (R + DeltaR) | w = R | w + d(R | w)/dR |R|w * DeltaRwedge
    // R % w = d(R % w)/dR |R%w * Rwedge
    // (R + DeltaR) % w = R % w + d(R % w)/dR |R%w * DeltaRwedge

    so3_tangent w_b(1.8, -2.1, -3.3);
    so3_tangent w_n = R_nb | w_b;

    Eigen::Matrix39d J_forward_adjoint = R_nb.jac_euclidean_forward_adjoint_wrt_dcm(w_b);
    Eigen::Vector3d w_n2               = J_forward_adjoint * R_nb_wedge;

    check("R jac forward dcm adjoint 0 linear    ", w_n()(0), w_n2(0), 1e-14);
    check("R jac forward dcm adjoint 1 linear    ", w_n()(1), w_n2(1), 1e-14);
    check("R jac forward dcm adjoint 2 linear    ", w_n()(2), w_n2(2), 1e-14);

    Eigen::Vector3d Sw_n  = SR_nb * w_b(); // WATCH OUT - This can not be the forward adjoint operator |
    Eigen::Vector3d SSw_n = w_n() + J_forward_adjoint * DeltaR_nb_wedge;

    check("R jac forward dcm adjoint 0 linear    ", Sw_n(0), SSw_n(0), 1e-14);
    check("R jac forward dcm adjoint 1 linear    ", Sw_n(1), SSw_n(1), 1e-14);
    check("R jac forward dcm adjoint 2 linear    ", Sw_n(2), SSw_n(2), 1e-14);

    Eigen::Matrix39d J_backward_adjoint = R_nb.jac_euclidean_backward_adjoint_wrt_dcm(w_n);
    Eigen::Vector3d w_b2                = J_backward_adjoint * R_nb_wedge;

    check("R jac backward dcm adjoint 0 linear   ", w_b()(0), w_b2(0), 1e-14);
    check("R jac backward dcm adjoint 1 linear   ", w_b()(1), w_b2(1), 1e-14);
    check("R jac backward dcm adjoint 2 linear   ", w_b()(2), w_b2(2), 1e-14);

    Eigen::Vector3d Sw_b  = SR_nb.transpose() * w_n(); // WATCH OUT - This can not be the inverse adjoint operator %
    Eigen::Vector3d SSw_b = w_b() + J_backward_adjoint * DeltaR_nb_wedge;

    check("R jac backward dcm adjoint 0 linear   ", Sw_b(0), SSw_b(0), 1e-14);
    check("R jac backward dcm adjoint 1 linear   ", Sw_b(1), SSw_b(1), 1e-14);
    check("R jac backward dcm adjoint 2 linear   ", Sw_b(2), SSw_b(2), 1e-14);

} // closes test_jac_euclidean_wrt_dcm_linear

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_wrt_rodrigues() {

    // It verifies the behavior of the rodrigues class jac_forward_rotation_wrt_rodrigues and
    // jac_backward_rotation_wrt_rodrigues, validating the following equations:
    // [q * exp(Delta r)] * v ~= q * v + d(q * v)/dq |q*v * {[q * exp(Delta r)] - q}
    // [exp(Delta r) * q] * v ~= q * v + d(q * v)/dq |q*v * {[exp(Delta r) * q] - q}
    // [q * exp(Delta r)] / v ~= q / v + d(q / v)/dq |q/v * {[q * exp(Delta r)] - q}
    // [exp(Delta r) * q] / v ~= q / v + d(q / v)/dq |q/v * {[exp(Delta r) * q] - q}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta r

    double d2r = math::constant::D2R();

    euler euler_nb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    rodrigues q_nb(euler_nb);
    quat q_nb_wedge = q_nb();

    rotv Deltar_nb(0.001, -0.002, 0.0015); // WATCH OUT - needs to be small
    rodrigues Deltaq_nb(Deltar_nb);
    quat Deltaq_nb_wedge = Deltaq_nb();

    Eigen::Vector3d v_b(3.8, -4.1, -2.3);
    Eigen::Vector3d v_n = q_nb * v_b;

    rodrigues        Sq_nb       = q_nb.plus_right(Deltar_nb); // WATCH OUT - same as q_nb * Deltaq_nb
    Eigen::Vector3d  Sv_n        = Sq_nb * v_b;
    Eigen::Matrix34d Jac_f_r     = q_nb.jac_euclidean_forward_rotation_wrt_rodrigues(v_b);
    Eigen::Vector3d  SDelta_v_n  = Jac_f_r * (Sq_nb().get() - q_nb().get()); // WATCH OUT - very different from Jac_f_r * Deltaq_nb().get()
    Eigen::Vector3d  SSv_n       = v_n + SDelta_v_n;

    check("q jac forward quat rotate 0 taylor     ", Sv_n(0), SSv_n(0), 1e-5);
    check("q jac forward quat rotate 1 taylor     ", Sv_n(1), SSv_n(1), 1e-6);
    check("q jac forward quat rotate 2 taylor     ", Sv_n(2), SSv_n(2), 1e-4);

    rodrigues       Tq_nb       = q_nb.plus_left(Deltar_nb);
    Eigen::Vector3d Tv_n        = Tq_nb * v_b;
    Eigen::Vector3d TDeltaB_v_n = Jac_f_r * (Tq_nb().get() - q_nb().get());
    Eigen::Vector3d TTv_n      = v_n + TDeltaB_v_n;

    check("q jac forward quat rotate 0 taylor     ", Tv_n(0), TTv_n(0), 1e-5);
    check("q jac forward quat rotate 1 taylor     ", Tv_n(1), TTv_n(1), 1e-6);
    check("q jac forward quat rotate 2 taylor     ", Tv_n(2), TTv_n(2), 1e-4);

    Eigen::Vector3d  Sv_b         = Sq_nb / v_n;
    Eigen::Matrix34d Jac_b_r      = q_nb.jac_euclidean_backward_rotation_wrt_rodrigues(v_n);
    Eigen::Vector3d  SDelta_v_b   = Jac_b_r * (Sq_nb().get() - q_nb().get());
    Eigen::Vector3d  SSv_b        = v_b + SDelta_v_b;

    check("q jac backward quat rotate 0 taylor    ", Sv_b(0), SSv_b(0), 1e-5);
    check("q jac backward quat rotate 1 taylor    ", Sv_b(1), SSv_b(1), 1e-5);
    check("q jac backward quat rotate 2 taylor    ", Sv_b(2), SSv_b(2), 1e-4);

    Eigen::Vector3d Tv_b           = Tq_nb / v_n;
    Eigen::Vector3d TDelta_v_b     = Jac_b_r * (Tq_nb().get() - q_nb().get());
    Eigen::Vector3d TTv_b          = v_b + TDelta_v_b;

    check("q jac backward quat rotate 0 taylor    ", Tv_b(0), TTv_b(0), 1e-5);
    check("q jac backward quat rotate 1 taylor    ", Tv_b(1), TTv_b(1), 1e-5);
    check("q jac backward quat rotate 2 taylor    ", Tv_b(2), TTv_b(2), 1e-4);

    // It verifies the behavior of the rodrigues class jac_forward_adjoint_wrt_rodrigues and
    // jac_backward_adjoint_wrt_rodrigues, validating the following equations:
    // [q * exp(Delta r)] | w ~= q | w + d(q | w)/dq |q|w * {[q * exp(Delta r)] - q}
    // [exp(Delta r) * q] | w ~= q | w + d(q | w)/dq |q|w * {[exp(Delta r) * q] - q}
    // [q * exp(Delta r)] % w ~= q % w + d(q % w)/dq |q%w * {[q * exp(Delta r)] - q}
    // [exp(Delta r) * q] % w ~= q % w + d(q % w)/dq |q%w * {[exp(Delta r) * q] - q}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta r

    ang::so3_tangent w_nbb(3.8, -4.1, -2.3);
    ang::so3_tangent w_nbn = q_nb | w_nbb;

    ang::so3_tangent Sw_nbn         = Sq_nb | w_nbb;
    Eigen::Matrix34d Jac_f_a        = q_nb.jac_euclidean_forward_adjoint_wrt_rodrigues(w_nbb);
    ang::so3_tangent SDeltaA_w_nbn(Jac_f_a * (Sq_nb().get() - q_nb().get()));
    ang::so3_tangent SSw_nbn        = w_nbn + SDeltaA_w_nbn;

    check("q jac forward quat adjoint 0 taylor    ", Sw_nbn()(0), SSw_nbn()(0), 1e-5);
    check("q jac forward quat adjoint 1 taylor    ", Sw_nbn()(1), SSw_nbn()(1), 1e-6);
    check("q jac forward quat adjoint 2 taylor    ", Sw_nbn()(2), SSw_nbn()(2), 1e-4);

    ang::so3_tangent Tw_nbn      = Tq_nb | w_nbb;
    ang::so3_tangent TDelta_w_nbn(Jac_f_a * (Tq_nb().get() - q_nb().get()));
    ang::so3_tangent TTw_nbn     = w_nbn + TDelta_w_nbn;

    check("q jac forward quat adjoint 0 taylor    ", Tw_nbn()(0), TTw_nbn()(0), 1e-5);
    check("q jac forward quat adjoint 1 taylor    ", Tw_nbn()(1), TTw_nbn()(1), 1e-6);
    check("q jac forward quat adjoint 2 taylor    ", Tw_nbn()(2), TTw_nbn()(2), 1e-4);

    ang::so3_tangent Sw_nbb         = Sq_nb % w_nbn;
    Eigen::Matrix34d Jac_b_a        = q_nb.jac_euclidean_backward_adjoint_wrt_rodrigues(w_nbn);
    ang::so3_tangent SDelta_w_nbb(Jac_b_a * (Sq_nb().get() - q_nb().get()));
    ang::so3_tangent SSw_nbb        = w_nbb + SDelta_w_nbb;

    check("q jac backward quat adjoint 0 taylor   ", Sw_nbb()(0), SSw_nbb()(0), 1e-5);
    check("q jac backward quat adjoint 1 taylor   ", Sw_nbb()(1), SSw_nbb()(1), 1e-5);
    check("q jac backward quat adjoint 2 taylor   ", Sw_nbb()(2), SSw_nbb()(2), 1e-4);

    ang::so3_tangent Tw_nbb         = Tq_nb % w_nbn;
    ang::so3_tangent TDelta_w_nbb(Jac_b_a * (Tq_nb().get() - q_nb().get()));
    ang::so3_tangent TTw_nbb        = w_nbb + TDelta_w_nbb;

    check("q jac backward quat adjoint 0 taylor   ", Tw_nbb()(0), TTw_nbb()(0), 1e-5);
    check("q jac backward quat adjoint 1 taylor   ", Tw_nbb()(1), TTw_nbb()(1), 1e-5);
    check("q jac backward quat adjoint 2 taylor   ", Tw_nbb()(2), TTw_nbb()(2), 1e-4);

} // closes test_jac_euclidean_wrt_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_wrt_rodrigues_bilinear() {

    // It verifies the behavior of the rodrigues class jac_euclidean_forward_rotation_wrt_rodrigues and
    // jac_euclidean_backward_rotation_wrt_rodrigues, validating the following equations:

    // q * v = 0.5 * d(q * v)/dq |q*v * qwedge
    // q / v = 0.5 * d(q / v)/dq |q/v * qwedge

    // NOTICE the 0.5 FACTOR because it is BILINEAR

    double d2r = math::constant::D2R();

    euler euler_nb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    rodrigues q_nb(euler_nb);
    rodrigues q_bn = q_nb.inverse();
    quat q_nb_wedge = q_nb();

    Eigen::Vector3d v_b(3.8, -4.1, -2.3);
    Eigen::Vector3d v_n = q_nb * v_b;

    euler Deltaeuler_nb(3.2 * d2r, 8.7 * d2r, -4.5 * d2r); // does NOT need to be small
    rodrigues Deltaq_nb(Deltaeuler_nb);
    quat Deltaq_nb_wedge = Deltaq_nb();

    Eigen::Matrix34d J_forward_rotation = q_nb.jac_euclidean_forward_rotation_wrt_rodrigues(v_b);
    Eigen::Vector3d v_n2                = 0.5 * J_forward_rotation * q_nb_wedge; // WATCH OUT - Note the 0.5 factor

    check("q jac forward quat rotate 0 linear     ", v_n(0), v_n2(0), 1e-14);
    check("q jac forward quat rotate 1 linear     ", v_n(1), v_n2(1), 1e-14);
    check("q jac forward quat rotate 2 linear     ", v_n(2), v_n2(2), 1e-14);

    Eigen::Matrix<double,3,4> jac_bwd_rotate1 = q_nb.jac_euclidean_backward_rotation_wrt_rodrigues(v_n);
    Eigen::Vector3d v_b2 = 0.5 * jac_bwd_rotate1 * q_nb().get();

    check("q jac backward quat rotate 0 bilinear  ", v_b(0), v_b2(0), 1e-14);
    check("q jac backward quat rotate 1 bilinear  ", v_b(1), v_b2(1), 1e-14);
    check("q jac backward quat rotate 2 bilinear  ", v_b(2), v_b2(2), 1e-14);

    Eigen::Matrix<double,3,4> jac_fwd_rotate2 = q_bn.jac_euclidean_backward_rotation_wrt_rodrigues(v_b);

    check("q jac_quat rotate 0 - 0     ", J_forward_rotation(0,0), + jac_fwd_rotate2(0,0), 1e-12);
    check("q jac_quat rotate 0 - 1     ", J_forward_rotation(0,1), - jac_fwd_rotate2(0,1), 1e-12);
    check("q jac_quat rotate 0 - 2     ", J_forward_rotation(0,2), - jac_fwd_rotate2(0,2), 1e-12);
    check("q jac_quat rotate 0 - 3     ", J_forward_rotation(0,3), - jac_fwd_rotate2(0,3), 1e-12);
    check("q jac_quat rotate 1 - 0     ", J_forward_rotation(1,0), + jac_fwd_rotate2(1,0), 1e-12);
    check("q jac_quat rotate 1 - 1     ", J_forward_rotation(1,1), - jac_fwd_rotate2(1,1), 1e-12);
    check("q jac_quat rotate 1 - 2     ", J_forward_rotation(1,2), - jac_fwd_rotate2(1,2), 1e-12);
    check("q jac_quat rotate 1 - 3     ", J_forward_rotation(1,3), - jac_fwd_rotate2(1,3), 1e-12);
    check("q jac_quat rotate 2 - 0     ", J_forward_rotation(2,0), + jac_fwd_rotate2(2,0), 1e-12);
    check("q jac_quat rotate 2 - 1     ", J_forward_rotation(2,1), - jac_fwd_rotate2(2,1), 1e-12);
    check("q jac_quat rotate 2 - 2     ", J_forward_rotation(2,2), - jac_fwd_rotate2(2,2), 1e-12);
    check("q jac_quat rotate 2 - 3     ", J_forward_rotation(2,3), - jac_fwd_rotate2(2,3), 1e-12);

    Eigen::Matrix<double,3,4> jac_bwd_rotate2 = q_bn.jac_euclidean_forward_rotation_wrt_rodrigues(v_n);

    check("q jac_quat rotate 0 - 0     ", jac_bwd_rotate1(0,0), + jac_bwd_rotate2(0,0), 1e-12);
    check("q jac_quat rotate 0 - 1     ", jac_bwd_rotate1(0,1), - jac_bwd_rotate2(0,1), 1e-12);
    check("q jac_quat rotate 0 - 2     ", jac_bwd_rotate1(0,2), - jac_bwd_rotate2(0,2), 1e-12);
    check("q jac_quat rotate 0 - 3     ", jac_bwd_rotate1(0,3), - jac_bwd_rotate2(0,3), 1e-12);
    check("q jac_quat rotate 1 - 0     ", jac_bwd_rotate1(1,0), + jac_bwd_rotate2(1,0), 1e-12);
    check("q jac_quat rotate 1 - 1     ", jac_bwd_rotate1(1,1), - jac_bwd_rotate2(1,1), 1e-12);
    check("q jac_quat rotate 1 - 2     ", jac_bwd_rotate1(1,2), - jac_bwd_rotate2(1,2), 1e-12);
    check("q jac_quat rotate 1 - 3     ", jac_bwd_rotate1(1,3), - jac_bwd_rotate2(1,3), 1e-12);
    check("q jac_quat rotate 2 - 0     ", jac_bwd_rotate1(2,0), + jac_bwd_rotate2(2,0), 1e-12);
    check("q jac_quat rotate 2 - 1     ", jac_bwd_rotate1(2,1), - jac_bwd_rotate2(2,1), 1e-12);
    check("q jac_quat rotate 2 - 2     ", jac_bwd_rotate1(2,2), - jac_bwd_rotate2(2,2), 1e-12);
    check("q jac_quat rotate 2 - 3     ", jac_bwd_rotate1(2,3), - jac_bwd_rotate2(2,3), 1e-12);

    // It verifies the behavior of the rodrigues class jacobian_rodrigues_forward_adjoint and
    // jacobian_rodrigues_backward_adjoint, validating the following equations:
    // q | w = 0.5 * d(q | w)/dq |q|w * q
    // q % w = 0.5 * d(q % w)/dq |q%w * q

    ang::so3_tangent w_nbb(3.8, -4.1, -2.3);
    ang::so3_tangent w_nbn = q_nb | w_nbb;

    Eigen::Matrix<double,3,4> jac_fwd_adjoint1 = q_nb.jac_euclidean_forward_adjoint_wrt_rodrigues(w_nbb);
    ang::so3_tangent w_nbn2(0.5 * jac_fwd_adjoint1 * q_nb().get());

    check("q jac forward quat adjoint 0 bilinear  ", w_nbn()(0), w_nbn2()(0), 1e-14);
    check("q jac forward quat adjoint 1 bilinear  ", w_nbn()(1), w_nbn2()(1), 1e-14);
    check("q jac forward quat adjoint 2 bilinear  ", w_nbn()(2), w_nbn2()(2), 1e-14);


    Eigen::Matrix<double,3,4> jac_bwd_adjoint1 = q_nb.jac_euclidean_backward_adjoint_wrt_rodrigues(w_nbn);
    ang::so3_tangent w_nbb2(0.5 * jac_bwd_adjoint1 * q_nb().get());

    check("q jac backward quat adjoint 0 bilinear ", w_nbb()(0), w_nbb2()(0), 1e-14);
    check("q jac backward quat adjoint 1 bilinear ", w_nbb()(1), w_nbb2()(1), 1e-14);
    check("q jac backward quat adjoint 2 bilinear ", w_nbb()(2), w_nbb2()(2), 1e-14);

    Eigen::Matrix<double,3,4> jac_fwd_adjoint2 = q_bn.jac_euclidean_backward_adjoint_wrt_rodrigues(w_nbb);

    check("q jac_quat adjoint 0 - 0    ", jac_fwd_adjoint1(0,0), + jac_fwd_adjoint2(0,0), 1e-12);
    check("q jac_quat adjoint 0 - 1    ", jac_fwd_adjoint1(0,1), - jac_fwd_adjoint2(0,1), 1e-12);
    check("q jac_quat adjoint 0 - 2    ", jac_fwd_adjoint1(0,2), - jac_fwd_adjoint2(0,2), 1e-12);
    check("q jac_quat adjoint 0 - 3    ", jac_fwd_adjoint1(0,3), - jac_fwd_adjoint2(0,3), 1e-12);
    check("q jac_quat adjoint 1 - 0    ", jac_fwd_adjoint1(1,0), + jac_fwd_adjoint2(1,0), 1e-12);
    check("q jac_quat adjoint 1 - 1    ", jac_fwd_adjoint1(1,1), - jac_fwd_adjoint2(1,1), 1e-12);
    check("q jac_quat adjoint 1 - 2    ", jac_fwd_adjoint1(1,2), - jac_fwd_adjoint2(1,2), 1e-12);
    check("q jac_quat adjoint 1 - 3    ", jac_fwd_adjoint1(1,3), - jac_fwd_adjoint2(1,3), 1e-12);
    check("q jac_quat adjoint 2 - 0    ", jac_fwd_adjoint1(2,0), + jac_fwd_adjoint2(2,0), 1e-12);
    check("q jac_quat adjoint 2 - 1    ", jac_fwd_adjoint1(2,1), - jac_fwd_adjoint2(2,1), 1e-12);
    check("q jac_quat adjoint 2 - 2    ", jac_fwd_adjoint1(2,2), - jac_fwd_adjoint2(2,2), 1e-12);
    check("q jac_quat adjoint 2 - 3    ", jac_fwd_adjoint1(2,3), - jac_fwd_adjoint2(2,3), 1e-12);

    Eigen::Matrix<double,3,4> jac_bwd_adjoint2 = q_bn.jac_euclidean_forward_adjoint_wrt_rodrigues(w_nbn);

    check("q jac_quat adjoint 0 - 0    ", jac_bwd_adjoint1(0,0), + jac_bwd_adjoint2(0,0), 1e-12);
    check("q jac_quat adjoint 0 - 1    ", jac_bwd_adjoint1(0,1), - jac_bwd_adjoint2(0,1), 1e-12);
    check("q jac_quat adjoint 0 - 2    ", jac_bwd_adjoint1(0,2), - jac_bwd_adjoint2(0,2), 1e-12);
    check("q jac_quat adjoint 0 - 3    ", jac_bwd_adjoint1(0,3), - jac_bwd_adjoint2(0,3), 1e-12);
    check("q jac_quat adjoint 1 - 0    ", jac_bwd_adjoint1(1,0), + jac_bwd_adjoint2(1,0), 1e-12);
    check("q jac_quat adjoint 1 - 1    ", jac_bwd_adjoint1(1,1), - jac_bwd_adjoint2(1,1), 1e-12);
    check("q jac_quat adjoint 1 - 2    ", jac_bwd_adjoint1(1,2), - jac_bwd_adjoint2(1,2), 1e-12);
    check("q jac_quat adjoint 1 - 3    ", jac_bwd_adjoint1(1,3), - jac_bwd_adjoint2(1,3), 1e-12);
    check("q jac_quat adjoint 2 - 0    ", jac_bwd_adjoint1(2,0), + jac_bwd_adjoint2(2,0), 1e-12);
    check("q jac_quat adjoint 2 - 1    ", jac_bwd_adjoint1(2,1), - jac_bwd_adjoint2(2,1), 1e-12);
    check("q jac_quat adjoint 2 - 2    ", jac_bwd_adjoint1(2,2), - jac_bwd_adjoint2(2,2), 1e-12);
    check("q jac_quat adjoint 2 - 3    ", jac_bwd_adjoint1(2,3), - jac_bwd_adjoint2(2,3), 1e-12);

} // closes test_jac_euclidean_wrt_rodrigues_bilinear

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3_jacobian::test_jac_euclidean_wrt_rotv() {

    // It checks the proper behavior of the "rotv" class methods related with the euclidean action jacobians,
    // validating the following equations:
    // (r + Delta r) * v ~= r * v + d(r * v)/dr |r*v * Delta r
    // (r + Delta r) / v ~= r / v + d(r / v)/dr |r/v * Delta r

    // It checks the proper behavior of the "rotv" class methods related with the euclidean adjoint jacobians,
    // validating the following equations:
    // (r + Delta r) | w ~= r | w + d(r | w)/dr |r|w * Delta r
    // (r + Delta r) % w ~= r % w + d(r % w)/dr |r%w * Delta r

    ang::rotv rv(0.4, -0.7, 0.3);
    ang::rotv Delta_rv(0.0002, 0.00025, -0.00015); // needs to be small
    ang::rotv Srv(rv() + Delta_rv());

    Eigen::Vector3d v_b(3.8, -4.1, -2.3);
    Eigen::Vector3d v_n   = rv * v_b;
    Eigen::Vector3d Sv_n  = Srv * v_b;
    Eigen::Vector3d Sv_n2 = v_n + rv.jac_euclidean_forward_rotation_wrt_rotv(v_b) * Delta_rv();
    Eigen::Vector3d Sv_n3 = v_n + rv.jac_euclidean_forward_rotation_wrt_rotv_bis(v_b) * Delta_rv();

    check("jac sum rotv rotate 0     ", Sv_n(0), Sv_n2(0), 1e-6);
    check("jac sum rotv rotate 1     ", Sv_n(1), Sv_n2(1), 1e-6);
    check("jac sum rotv rotate 2     ", Sv_n(2), Sv_n2(2), 1e-8);
    check("jac sum rotv rotate 0     ", Sv_n(0), Sv_n3(0), 1e-6);
    check("jac sum rotv rotate 1     ", Sv_n(1), Sv_n3(1), 1e-6);
    check("jac sum rotv rotate 2     ", Sv_n(2), Sv_n3(2), 1e-8);

    ang::so3_tangent w_nbb(3.8, -4.1, -2.3);
    ang::so3_tangent w_nbn   = rv | w_nbb;
    ang::so3_tangent Sw_nbn  = Srv | w_nbb;
    ang::so3_tangent Sw_nbn2 = w_nbn + rv.jac_euclidean_forward_adjoint_wrt_rotv(w_nbb) * Delta_rv();
    ang::so3_tangent Sw_nbn3 = w_nbn + rv.jac_euclidean_forward_adjoint_wrt_rotv_bis(w_nbb) * Delta_rv();

    check("jac sum rotv adjoint 0    ", Sw_nbn()(0), Sw_nbn2()(0), 1e-6);
    check("jac sum rotv adjoint 1    ", Sw_nbn()(1), Sw_nbn2()(1), 1e-6);
    check("jac sum rotv adjoint 2    ", Sw_nbn()(2), Sw_nbn2()(2), 1e-8);
    check("jac sum rotv adjoint 0    ", Sw_nbn()(0), Sw_nbn3()(0), 1e-6);
    check("jac sum rotv adjoint 1    ", Sw_nbn()(1), Sw_nbn3()(1), 1e-6);
    check("jac sum rotv adjoint 2    ", Sw_nbn()(2), Sw_nbn3()(2), 1e-8);

    Eigen::Vector3d Sv_b  = Srv / v_n;
    Eigen::Vector3d Sv_b2 = v_b + rv.jac_euclidean_backward_rotation_wrt_rotv(v_n) * Delta_rv();
    Eigen::Vector3d Sv_b3 = v_b + rv.jac_euclidean_backward_rotation_wrt_rotv_bis(v_n) * Delta_rv();

    check("jac sum rotv rotate 0     ", Sv_b(0), Sv_b2(0), 1e-6);
    check("jac sum rotv rotate 1     ", Sv_b(1), Sv_b2(1), 1e-6);
    check("jac sum rotv rotate 2     ", Sv_b(2), Sv_b2(2), 1e-7);
    check("jac sum rotv rotate 0     ", Sv_b(0), Sv_b3(0), 1e-6);
    check("jac sum rotv rotate 1     ", Sv_b(1), Sv_b3(1), 1e-6);
    check("jac sum rotv rotate 2     ", Sv_b(2), Sv_b3(2), 1e-7);

    ang::so3_tangent Sw_nbb = Srv % w_nbn;
    ang::so3_tangent Sw_nbb2 = w_nbb + rv.jac_euclidean_backward_adjoint_wrt_rotv(w_nbn) * Delta_rv();
    ang::so3_tangent Sw_nbb3 = w_nbb + rv.jac_euclidean_backward_adjoint_wrt_rotv_bis(w_nbn) * Delta_rv();

    check("jac sum rotv adjoint 0    ", Sw_nbb()(0), Sw_nbb2()(0), 1e-6);
    check("jac sum rotv adjoint 1    ", Sw_nbb()(1), Sw_nbb2()(1), 1e-6);
    check("jac sum rotv adjoint 2    ", Sw_nbb()(2), Sw_nbb2()(2), 1e-7);
    check("jac sum rotv adjoint 0    ", Sw_nbb()(0), Sw_nbb3()(0), 1e-6);
    check("jac sum rotv adjoint 1    ", Sw_nbb()(1), Sw_nbb3()(1), 1e-6);
    check("jac sum rotv adjoint 2    ", Sw_nbb()(2), Sw_nbb3()(2), 1e-7);

} // closes test_jac_euclidean_wrt_rotv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

















