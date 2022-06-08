#include "Tse3_jacobian.h"
#include "../ang/rotate/dcm.h"
#include "../ang/rotate/rodrigues.h"
#include "../ang/rotate/euler.h"
#include "../ang/rotate/rotv.h"
#include "../ang/transform/speu_rodrigues.h"
#include "../ang/transform/speu_dcm.h"
#include "../ang/transform/homogeneous.h"
#include "../ang/transform/trfv.h"
#include "../ang/transform/screw.h"
#include "../ang/transform/dual.h"
#include "../ang/transform/se3_tangent.h"
#include "../ang/tools.h"
#include <iostream>

using namespace std;

ang::test::Tse3_jacobian::Tse3_jacobian(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tse3_jacobian::run() {
	::jail::unit_test::run();
    test_jac_right_left();                              cout << endl << endl;

    test_jac_right_inverse();                           cout << endl << endl;
    test_jac_right_composition_wrt_first();             cout << endl << endl;
    test_jac_right_composition_wrt_second();            cout << endl << endl;
    test_jac_right_forward_motion_wrt_motion();         cout << endl << endl;
    test_jac_right_backward_motion_wrt_motion();        cout << endl << endl;
    test_jac_right_log();                               cout << endl << endl;
    test_jac_right_plus_wrt_first();                    cout << endl << endl;
    test_jac_right_plus_wrt_second();                   cout << endl << endl;
    test_jac_right_minus_wrt_first();                   cout << endl << endl;
    test_jac_right_minus_wrt_second();                  cout << endl << endl;
    test_jac_right_forward_adjoint_wrt_motion();        cout << endl << endl;
    test_jac_right_backward_adjoint_wrt_motion();       cout << endl << endl;

    test_jac_left_inverse();                            cout << endl << endl;
    test_jac_left_composition_wrt_first();              cout << endl << endl;
    test_jac_left_composition_wrt_second();             cout << endl << endl;
    test_jac_left_forward_motion_wrt_motion();          cout << endl << endl;
    test_jac_left_backward_motion_wrt_motion();         cout << endl << endl;
    test_jac_left_log();                                cout << endl << endl;
    test_jac_left_plus_wrt_first();                     cout << endl << endl;
    test_jac_left_plus_wrt_second();                    cout << endl << endl;
    test_jac_left_minus_wrt_first();                    cout << endl << endl;
    test_jac_left_minus_wrt_second();                   cout << endl << endl;
    test_jac_left_forward_adjoint_wrt_motion();         cout << endl << endl;
    test_jac_left_backward_adjoint_wrt_motion();        cout << endl << endl;

    test_jac_euclidean_forward_motion_wrt_point();      cout << endl << endl;
    test_jac_euclidean_backward_motion_wrt_point();     cout << endl << endl;
    test_jac_euclidean_forward_adjoint_wrt_tangent();   cout << endl << endl;
    test_jac_euclidean_backward_adjoint_wrt_tangent();  cout << endl << endl;

    test_jac_euclidean_wrt_speu_dcm();                  cout << endl << endl;
    test_jac_euclidean_wrt_speu_dcm_linear();           cout << endl << endl;
    test_jac_euclidean_wrt_speu_rodrigues();            cout << endl << endl;
    test_jac_euclidean_wrt_homogeneous();               cout << endl << endl;
    test_jac_euclidean_wrt_homogeneous_linear();        cout << endl << endl;
    test_jac_euclidean_wrt_trfv();                      cout << endl << endl;
    test_jac_euclidean_wrt_trfv_bis();                  cout << endl << endl;

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_left() {

    // It checks the proper behavior of the "trfv" right and left jacobian methods,
    // validating the following equations:
    // exp(tau + Deltatau) ~= exp(tau) plus (JR * Deltatau)
    // exp(tau + Deltatau) ~= (JL * Deltatau) plus_left exp(tau) */
    // tau + jacRinv * Deltatau ~= log[exp(tau) plus Deltatau]
    // tau + JLinv * Deltatau ~= log[Deltatau plus_left exp(tau)] */

    Eigen::Vector6d Itau; Itau << 4.0, 11.0, 7.0, 0.3, 0.5, 0.2;
    trfv tau(Itau);
    speu_dcm gR = tau.exp_map_speu_dcm();

    Eigen::Vector6d IDelta_tau; IDelta_tau << 0.0003, -0.0008, 0.0002, -0.0002, -0.00031, 0.00013; // needs to be small
    trfv Delta_tau(IDelta_tau);
    speu_dcm Delta_GR = Delta_tau.exp_map_speu_dcm();

    // I do an algebraic sum of the transform vectors (not a concatenation - composition).
    // Then I compute the perturbation as the product of the jacobian by the transform vector
    trfv Stau(tau() + Delta_tau());
    speu_dcm SgR = Stau.exp_map_speu_dcm();

    trfv Delta_tau_right(tau.jac_right() * Delta_tau());
    speu_dcm SgR_right = gR.plus_right(Delta_tau_right);

    check("jac right 0 - 0  ", SgR.get_dcm()()(0,0), SgR_right.get_dcm()()(0,0), 1e-8);
    check("jac right 0 - 1  ", SgR.get_dcm()()(0,1), SgR_right.get_dcm()()(0,1), 1e-8);
    check("jac right 0 - 2  ", SgR.get_dcm()()(0,2), SgR_right.get_dcm()()(0,2), 1e-8);
    check("jac right 1 - 0  ", SgR.get_dcm()()(1,0), SgR_right.get_dcm()()(1,0), 1e-8);
    check("jac right 1 - 1  ", SgR.get_dcm()()(1,1), SgR_right.get_dcm()()(1,1), 1e-8);
    check("jac right 1 - 2  ", SgR.get_dcm()()(1,2), SgR_right.get_dcm()()(1,2), 1e-8);
    check("jac right 2 - 0  ", SgR.get_dcm()()(2,0), SgR_right.get_dcm()()(2,0), 1e-8);
    check("jac right 2 - 1  ", SgR.get_dcm()()(2,1), SgR_right.get_dcm()()(2,1), 1e-8);
    check("jac right 2 - 2  ", SgR.get_dcm()()(2,2), SgR_right.get_dcm()()(2,2), 1e-8);
    check("jac right 0      ", SgR.get_T()(0),       SgR_right.get_T()(0),       1e-7);
    check("jac right 1      ", SgR.get_T()(1),       SgR_right.get_T()(1),       1e-7);
    check("jac right 2      ", SgR.get_T()(2),       SgR_right.get_T()(2),       1e-6);

    trfv Delta_tau_left(tau.jac_left() * Delta_tau());
    speu_dcm SgR_left = gR.plus_left(Delta_tau_left);

    check("jac left 0 - 0   ", SgR.get_dcm()()(0,0), SgR_left.get_dcm()()(0,0), 1e-8);
    check("jac left 0 - 1   ", SgR.get_dcm()()(0,1), SgR_left.get_dcm()()(0,1), 1e-8);
    check("jac left 0 - 2   ", SgR.get_dcm()()(0,2), SgR_left.get_dcm()()(0,2), 1e-8);
    check("jac left 1 - 0   ", SgR.get_dcm()()(1,0), SgR_left.get_dcm()()(1,0), 1e-8);
    check("jac left 1 - 1   ", SgR.get_dcm()()(1,1), SgR_left.get_dcm()()(1,1), 1e-8);
    check("jac left 1 - 2   ", SgR.get_dcm()()(1,2), SgR_left.get_dcm()()(1,2), 1e-8);
    check("jac left 2 - 0   ", SgR.get_dcm()()(2,0), SgR_left.get_dcm()()(2,0), 1e-8);
    check("jac left 2 - 1   ", SgR.get_dcm()()(2,1), SgR_left.get_dcm()()(2,1), 1e-8);
    check("jac left 2 - 2   ", SgR.get_dcm()()(2,2), SgR_left.get_dcm()()(2,2), 1e-8);
    check("jac left 0       ", SgR.get_T()(0),       SgR_left.get_T()(0),       1e-7);
    check("jac left 1       ", SgR.get_T()(1),       SgR_left.get_T()(1),       1e-7);
    check("jac left 2       ", SgR.get_T()(2),       SgR_left.get_T()(2),       1e-6);

    Eigen::Vector6d tau_right_add  = tau() + tau.jac_right_inv() * Delta_tau();
    Eigen::Vector6d tau_right_add2 = (tau.exp_map_speu_dcm().plus_right( Delta_tau)).log_map_trfv()();

    check("jac right inv 0  ", tau_right_add(0), tau_right_add2(0), 1e-7);
    check("jac right inv 1  ", tau_right_add(1), tau_right_add2(1), 1e-7);
    check("jac right inv 2  ", tau_right_add(2), tau_right_add2(2), 1e-6);
    check("jac right inv 3  ", tau_right_add(3), tau_right_add2(3), 1e-8);
    check("jac right inv 4  ", tau_right_add(4), tau_right_add2(4), 1e-8);
    check("jac right inv 5  ", tau_right_add(5), tau_right_add2(5), 1e-8);

    Eigen::Vector6d tau_left_add = tau() + tau.jac_left_inv() * Delta_tau();
    Eigen::Vector6d tau_left_add2 = (tau.exp_map_speu_dcm().plus_left(Delta_tau)).log_map_trfv()();

    check("jac left inv 0   ", tau_left_add(0), tau_left_add2(0), 1e-7);
    check("jac left inv 1   ", tau_left_add(1), tau_left_add2(1), 1e-7);
    check("jac left inv 2   ", tau_left_add(2), tau_left_add2(2), 1e-6);
    check("jac left inv 3   ", tau_left_add(3), tau_left_add2(3), 1e-8);
    check("jac left inv 4   ", tau_left_add(4), tau_left_add2(4), 1e-8);
    check("jac left inv 5   ", tau_left_add(5), tau_left_add2(5), 1e-8);

    // check that what is called inverse is in fact the inverse
    Eigen::Matrix6d J_right      = tau.jac_right();
    Eigen::Matrix6d J_right_inv  = tau.jac_right_inv();
    Eigen::Matrix6d J_right_inv2 = J_right.inverse();

    check("jac right inverse 0 - 0  ", J_right_inv(0,0), J_right_inv2(0,0), 1e-12);
    check("jac right inverse 0 - 1  ", J_right_inv(0,1), J_right_inv2(0,1), 1e-12);
    check("jac right inverse 0 - 2  ", J_right_inv(0,2), J_right_inv2(0,2), 1e-12);
    check("jac right inverse 0 - 3  ", J_right_inv(0,3), J_right_inv2(0,3), 1e-12);
    check("jac right inverse 0 - 4  ", J_right_inv(0,4), J_right_inv2(0,4), 1e-12);
    check("jac right inverse 0 - 5  ", J_right_inv(0,5), J_right_inv2(0,5), 1e-12);

    check("jac right inverse 1 - 0  ", J_right_inv(1,0), J_right_inv2(1,0), 1e-12);
    check("jac right inverse 1 - 1  ", J_right_inv(1,1), J_right_inv2(1,1), 1e-12);
    check("jac right inverse 1 - 2  ", J_right_inv(1,2), J_right_inv2(1,2), 1e-12);
    check("jac right inverse 1 - 3  ", J_right_inv(1,3), J_right_inv2(1,3), 1e-12);
    check("jac right inverse 1 - 4  ", J_right_inv(1,4), J_right_inv2(1,4), 1e-12);
    check("jac right inverse 1 - 5  ", J_right_inv(1,5), J_right_inv2(1,5), 1e-12);

    check("jac right inverse 2 - 0  ", J_right_inv(2,0), J_right_inv2(2,0), 1e-12);
    check("jac right inverse 2 - 1  ", J_right_inv(2,1), J_right_inv2(2,1), 1e-12);
    check("jac right inverse 2 - 2  ", J_right_inv(2,2), J_right_inv2(2,2), 1e-12);
    check("jac right inverse 2 - 3  ", J_right_inv(2,3), J_right_inv2(2,3), 1e-12);
    check("jac right inverse 2 - 4  ", J_right_inv(2,4), J_right_inv2(2,4), 1e-12);
    check("jac right inverse 2 - 5  ", J_right_inv(2,5), J_right_inv2(2,5), 1e-12);

    check("jac right inverse 3 - 0  ", J_right_inv(3,0), J_right_inv2(3,0), 1e-12);
    check("jac right inverse 3 - 1  ", J_right_inv(3,1), J_right_inv2(3,1), 1e-12);
    check("jac right inverse 3 - 2  ", J_right_inv(3,2), J_right_inv2(3,2), 1e-12);
    check("jac right inverse 3 - 3  ", J_right_inv(3,3), J_right_inv2(3,3), 1e-12);
    check("jac right inverse 3 - 4  ", J_right_inv(3,4), J_right_inv2(3,4), 1e-12);
    check("jac right inverse 3 - 5  ", J_right_inv(3,5), J_right_inv2(3,5), 1e-12);

    check("jac right inverse 4 - 0  ", J_right_inv(4,0), J_right_inv2(4,0), 1e-12);
    check("jac right inverse 4 - 1  ", J_right_inv(4,1), J_right_inv2(4,1), 1e-12);
    check("jac right inverse 4 - 2  ", J_right_inv(4,2), J_right_inv2(4,2), 1e-12);
    check("jac right inverse 4 - 3  ", J_right_inv(4,3), J_right_inv2(4,3), 1e-12);
    check("jac right inverse 4 - 4  ", J_right_inv(4,4), J_right_inv2(4,4), 1e-12);
    check("jac right inverse 4 - 5  ", J_right_inv(4,5), J_right_inv2(4,5), 1e-12);

    check("jac right inverse 5 - 0  ", J_right_inv(5,0), J_right_inv2(5,0), 1e-12);
    check("jac right inverse 5 - 1  ", J_right_inv(5,1), J_right_inv2(5,1), 1e-12);
    check("jac right inverse 5 - 2  ", J_right_inv(5,2), J_right_inv2(5,2), 1e-12);
    check("jac right inverse 5 - 3  ", J_right_inv(5,3), J_right_inv2(5,3), 1e-12);
    check("jac right inverse 5 - 4  ", J_right_inv(5,4), J_right_inv2(5,4), 1e-12);
    check("jac right inverse 5 - 5  ", J_right_inv(5,5), J_right_inv2(5,5), 1e-12);

    Eigen::Matrix6d J_left      = tau.jac_left();
    Eigen::Matrix6d J_left_inv  = tau.jac_left_inv();
    Eigen::Matrix6d J_left_inv2 = J_left.inverse();

    check("jac left inverse 0 - 0  ", J_left_inv(0,0), J_left_inv2(0,0), 1e-12);
    check("jac left inverse 0 - 1  ", J_left_inv(0,1), J_left_inv2(0,1), 1e-12);
    check("jac left inverse 0 - 2  ", J_left_inv(0,2), J_left_inv2(0,2), 1e-12);
    check("jac left inverse 0 - 3  ", J_left_inv(0,3), J_left_inv2(0,3), 1e-12);
    check("jac left inverse 0 - 4  ", J_left_inv(0,4), J_left_inv2(0,4), 1e-12);
    check("jac left inverse 0 - 5  ", J_left_inv(0,5), J_left_inv2(0,5), 1e-12);

    check("jac left inverse 1 - 0  ", J_left_inv(1,0), J_left_inv2(1,0), 1e-12);
    check("jac left inverse 1 - 1  ", J_left_inv(1,1), J_left_inv2(1,1), 1e-12);
    check("jac left inverse 1 - 2  ", J_left_inv(1,2), J_left_inv2(1,2), 1e-12);
    check("jac left inverse 1 - 3  ", J_left_inv(1,3), J_left_inv2(1,3), 1e-12);
    check("jac left inverse 1 - 4  ", J_left_inv(1,4), J_left_inv2(1,4), 1e-12);
    check("jac left inverse 1 - 5  ", J_left_inv(1,5), J_left_inv2(1,5), 1e-12);

    check("jac left inverse 2 - 0  ", J_left_inv(2,0), J_left_inv2(2,0), 1e-12);
    check("jac left inverse 2 - 1  ", J_left_inv(2,1), J_left_inv2(2,1), 1e-12);
    check("jac left inverse 2 - 2  ", J_left_inv(2,2), J_left_inv2(2,2), 1e-12);
    check("jac left inverse 2 - 3  ", J_left_inv(2,3), J_left_inv2(2,3), 1e-12);
    check("jac left inverse 2 - 4  ", J_left_inv(2,4), J_left_inv2(2,4), 1e-12);
    check("jac left inverse 2 - 5  ", J_left_inv(2,5), J_left_inv2(2,5), 1e-12);

    check("jac left inverse 3 - 0  ", J_left_inv(3,0), J_left_inv2(3,0), 1e-12);
    check("jac left inverse 3 - 1  ", J_left_inv(3,1), J_left_inv2(3,1), 1e-12);
    check("jac left inverse 3 - 2  ", J_left_inv(3,2), J_left_inv2(3,2), 1e-12);
    check("jac left inverse 3 - 3  ", J_left_inv(3,3), J_left_inv2(3,3), 1e-12);
    check("jac left inverse 3 - 4  ", J_left_inv(3,4), J_left_inv2(3,4), 1e-12);
    check("jac left inverse 3 - 5  ", J_left_inv(3,5), J_left_inv2(3,5), 1e-12);

    check("jac left inverse 4 - 0  ", J_left_inv(4,0), J_left_inv2(4,0), 1e-12);
    check("jac left inverse 4 - 1  ", J_left_inv(4,1), J_left_inv2(4,1), 1e-12);
    check("jac left inverse 4 - 2  ", J_left_inv(4,2), J_left_inv2(4,2), 1e-12);
    check("jac left inverse 4 - 3  ", J_left_inv(4,3), J_left_inv2(4,3), 1e-12);
    check("jac left inverse 4 - 4  ", J_left_inv(4,4), J_left_inv2(4,4), 1e-12);
    check("jac left inverse 4 - 5  ", J_left_inv(4,5), J_left_inv2(4,5), 1e-12);

    check("jac left inverse 5 - 0  ", J_left_inv(5,0), J_left_inv2(5,0), 1e-12);
    check("jac left inverse 5 - 1  ", J_left_inv(5,1), J_left_inv2(5,1), 1e-12);
    check("jac left inverse 5 - 2  ", J_left_inv(5,2), J_left_inv2(5,2), 1e-12);
    check("jac left inverse 5 - 3  ", J_left_inv(5,3), J_left_inv2(5,3), 1e-12);
    check("jac left inverse 5 - 4  ", J_left_inv(5,4), J_left_inv2(5,4), 1e-12);
    check("jac left inverse 5 - 5  ", J_left_inv(5,5), J_left_inv2(5,5), 1e-12);
} // closes test_jac_right_left

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_inverse() {

    // It verifies the behavior of the rotation classes jac_right_inverse method, validating the following:
    // (M plus DeltatauB)^-1 = M^-1 plus J * DeltatauB

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    speu_rodrigues Sgq_eb  = gq_eb.plus_right(DeltatauB);
    speu_rodrigues gq_be   = gq_eb.inverse();
    speu_rodrigues Sgq_be  = Sgq_eb.inverse();
    Eigen::Matrix6d Jac_gq = gq_eb.jac_right_inverse();
    trfv Deltatau_be_gq    = trfv(Jac_gq * DeltatauB());
    speu_rodrigues SSgq_be = gq_be.plus_right(Deltatau_be_gq);

    check("gq jac right inverse 1 taylor   ", Sgq_be.get_rotv()()(0), SSgq_be.get_rotv()()(0), 1e-12);
    check("gq jac right inverse 2 taylor   ", Sgq_be.get_rotv()()(1), SSgq_be.get_rotv()()(1), 1e-12);
    check("gq jac right inverse 3 taylor   ", Sgq_be.get_rotv()()(2), SSgq_be.get_rotv()()(2), 1e-12);
    check("gq jac right inverse 4 taylor   ", Sgq_be.get_T()(0), SSgq_be.get_T()(0), 1e-12);
    check("gq jac right inverse 5 taylor   ", Sgq_be.get_T()(1), SSgq_be.get_T()(1), 1e-12);
    check("gq jac right inverse 6 taylor   ", Sgq_be.get_T()(2), SSgq_be.get_T()(2), 1e-12);

    speu_dcm SgR_eb        = gR_eb.plus_right(DeltatauB);
    speu_dcm gR_be         = gR_eb.inverse();
    speu_dcm SgR_be        = SgR_eb.inverse();
    Eigen::Matrix6d Jac_gR = gR_eb.jac_right_inverse();
    trfv Deltatau_be_gR    = trfv(Jac_gR * DeltatauB());
    speu_dcm SSgR_be       = gR_be.plus_right(Deltatau_be_gR);

    check("gR jac right inverse 1 taylor   ", SgR_be.get_rotv()()(0), SSgR_be.get_rotv()()(0), 1e-12);
    check("gR jac right inverse 2 taylor   ", SgR_be.get_rotv()()(1), SSgR_be.get_rotv()()(1), 1e-12);
    check("gR jac right inverse 3 taylor   ", SgR_be.get_rotv()()(2), SSgR_be.get_rotv()()(2), 1e-12);
    check("gR jac right inverse 4 taylor   ", SgR_be.get_T()(0), SSgR_be.get_T()(0), 1e-12);
    check("gR jac right inverse 5 taylor   ", SgR_be.get_T()(1), SSgR_be.get_T()(1), 1e-12);
    check("gR jac right inverse 6 taylor   ", SgR_be.get_T()(2), SSgR_be.get_T()(2), 1e-12);

    homogeneous SM_eb     = M_eb.plus_right(DeltatauB);
    homogeneous M_be      = M_eb.inverse();
    homogeneous SM_be     = SM_eb.inverse();
    Eigen::Matrix6d Jac_M = M_eb.jac_right_inverse();
    trfv Deltatau_be_M    = trfv(Jac_M * DeltatauB());
    homogeneous SSM_be    = M_be.plus_right(Deltatau_be_M);

    check("M jac right inverse 1 taylor    ", SM_be.get_rotv()()(0), SSM_be.get_rotv()()(0), 1e-12);
    check("M jac right inverse 2 taylor    ", SM_be.get_rotv()()(1), SSM_be.get_rotv()()(1), 1e-12);
    check("M jac right inverse 3 taylor    ", SM_be.get_rotv()()(2), SSM_be.get_rotv()()(2), 1e-12);
    check("M jac right inverse 4 taylor    ", SM_be.get_T()(0), SSM_be.get_T()(0), 1e-12);
    check("M jac right inverse 5 taylor    ", SM_be.get_T()(1), SSM_be.get_T()(1), 1e-12);
    check("M jac right inverse 6 taylor    ", SM_be.get_T()(2), SSM_be.get_T()(2), 1e-12);

    dual SZ_eb            = Z_eb.plus_right(DeltatauB);
    dual Z_be             = Z_eb.inverse();
    dual SZ_be            = SZ_eb.inverse();
    Eigen::Matrix6d Jac_Z = Z_eb.jac_right_inverse();
    trfv Deltatau_be_Z    = trfv(Jac_Z * DeltatauB());
    dual SSZ_be           = Z_be.plus_right(Deltatau_be_Z);

    check("Z jac right inverse 1 taylor    ", SZ_be.get_rotv()()(0), SSZ_be.get_rotv()()(0), 1e-12);
    check("Z jac right inverse 2 taylor    ", SZ_be.get_rotv()()(1), SSZ_be.get_rotv()()(1), 1e-12);
    check("Z jac right inverse 3 taylor    ", SZ_be.get_rotv()()(2), SSZ_be.get_rotv()()(2), 1e-12);
    check("Z jac right inverse 4 taylor    ", SZ_be.get_T()(0), SSZ_be.get_T()(0), 1e-12);
    check("Z jac right inverse 5 taylor    ", SZ_be.get_T()(1), SSZ_be.get_T()(1), 1e-12);
    check("Z jac right inverse 6 taylor    ", SZ_be.get_T()(2), SSZ_be.get_T()(2), 1e-12);

    trfv tau_be             = tau_eb.inverse();
    trfv Stau_eb            = tau_eb.plus_right(DeltatauB);
    trfv Stau_be            = Stau_eb.inverse();
    Eigen::Matrix6d Jac_tau = tau_eb.jac_right_inverse();
    trfv Deltatau_be_tau    = trfv(Jac_tau * DeltatauB());
    trfv SStau_be           = tau_be.plus_right(Deltatau_be_tau);

    check("tau jac right inverse 1 taylor  ", Stau_be()(0), SStau_be()(0), 1e-12);
    check("tau jac right inverse 2 taylor  ", Stau_be()(1), SStau_be()(1), 1e-12);
    check("tau jac right inverse 3 taylor  ", Stau_be()(2), SStau_be()(2), 1e-12);
    check("tau jac right inverse 4 taylor  ", Stau_be()(3), SStau_be()(3), 1e-12);
    check("tau jac right inverse 5 taylor  ", Stau_be()(4), SStau_be()(4), 1e-12);
    check("tau jac right inverse 6 taylor  ", Stau_be()(5), SStau_be()(5), 1e-12);

} // closes test_jac_right_inverse

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_composition_wrt_first() {

    // It verifies the behavior of the motion classes jac_right_composition_wrt_first method, validating the following:
    // (M1 plus DeltatauB1) * M2 = M1 * M2 plus J * DeltatauB1

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);

    trfv tau_bh(r_bh, T_bhb);
    speu_dcm gR_bh(tau_bh);
    speu_rodrigues gq_bh(tau_bh);
    homogeneous M_bh(tau_bh);
    dual Z_bh(tau_bh);

    speu_rodrigues Sgq_eb    = gq_eb.plus_right(DeltatauB);
    speu_rodrigues gq_eh     = gq_eb * gq_bh;
    speu_rodrigues Sgq_eh    = Sgq_eb * gq_bh;
    Eigen::Matrix6d Jac_gq   = gq_eb.jac_right_composition_wrt_first(gq_bh);
    trfv Deltatau_bh_gq      = ang::trfv(Jac_gq * DeltatauB());
    speu_rodrigues SSgq_eh   = gq_eh.plus_right(Deltatau_bh_gq);

    check("gq jac right composition first 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-12);
    check("gq jac right composition first 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-12);
    check("gq jac right composition first 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-12);
    check("gq jac right composition first 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-12);
    check("gq jac right composition first 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-12);
    check("gq jac right composition first 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-12);

    speu_dcm SgR_eb        = gR_eb.plus_right(DeltatauB);
    speu_dcm gR_eh         = gR_eb * gR_bh;
    speu_dcm SgR_eh        = SgR_eb * gR_bh;
    Eigen::Matrix6d Jac_gR = gR_eb.jac_right_composition_wrt_first(gR_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauB());
    speu_dcm SSgR_eh       = gR_eh.plus_right(Deltatau_bh_gR);

    check("gR jac right composition first 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-12);
    check("gR jac right composition first 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-12);
    check("gR jac right composition first 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-12);
    check("gR jac right composition first 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-12);
    check("gR jac right composition first 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-12);
    check("gR jac right composition first 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-12);

    homogeneous SM_eb     = M_eb.plus_right(DeltatauB);
    homogeneous M_eh      = M_eb * M_bh;
    homogeneous SM_eh     = SM_eb * M_bh;
    Eigen::Matrix6d Jac_M = M_eb.jac_right_composition_wrt_first(M_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauB());
    homogeneous SSM_eh    = M_eh.plus_right(Deltatau_bh_M);

    check("M jac right composition first 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-12);
    check("M jac right composition first 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-12);
    check("M jac right composition first 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-12);
    check("M jac right composition first 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-12);
    check("M jac right composition first 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-12);
    check("M jac right composition first 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-12);

    dual SZ_eb            = Z_eb.plus_right(DeltatauB);
    dual Z_eh             = Z_eb * Z_bh;
    dual SZ_eh            = SZ_eb * Z_bh;
    Eigen::Matrix6d Jac_Z = Z_eb.jac_right_composition_wrt_first(Z_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauB());
    dual SSZ_eh           = Z_eh.plus_right(Deltatau_bh_Z);

    check("Z jac right composition first 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-12);
    check("Z jac right composition first 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-12);
    check("Z jac right composition first 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-12);
    check("Z jac right composition first 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-12);
    check("Z jac right composition first 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-12);
    check("Z jac right composition first 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-12);

    trfv Stau_eb            = tau_eb.plus_right(DeltatauB);
    trfv tau_eh             = tau_eb * tau_bh;
    trfv Stau_eh            = Stau_eb * tau_bh;
    Eigen::Matrix6d Jac_tau = tau_eb.jac_right_composition_wrt_first(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauB());
    trfv SStau_eh           = tau_eh.plus_right(Deltatau_bh_tau);

    check("tau jac right composition first 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-12);
    check("tau jac right composition first 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-12);
    check("tau jac right composition first 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-12);
    check("tau jac right composition first 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-12);
    check("tau jac right composition first 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-12);
    check("tau jac right composition first 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-12);

} // closes test_jac_right_composition_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_composition_wrt_second() {

    // It verifies the behavior of the motion classes jac_right_composition_wrt_second method, validating the following:
    // M1 * (M2 plus DeltatauB2) = M1 * M2 plus J * DeltatauB2

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);

    trfv tau_bh(r_bh, T_bhb);
    speu_dcm gR_bh(tau_bh);
    speu_rodrigues gq_bh(tau_bh);
    homogeneous M_bh(tau_bh);
    dual Z_bh(tau_bh);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauH(v);

    speu_rodrigues Sgq_bh    = gq_bh.plus_right(DeltatauH);
    speu_rodrigues gq_eh     = gq_eb * gq_bh;
    speu_rodrigues Sgq_eh    = gq_eb * Sgq_bh;
    Eigen::Matrix6d Jac_gq   = gq_eb.jac_right_composition_wrt_second(gq_bh);
    trfv Deltatau_bh_gq      = ang::trfv(Jac_gq * DeltatauH());
    speu_rodrigues SSgq_eh   = gq_eh.plus_right(Deltatau_bh_gq);

    check("gq jac right composition second 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-12);
    check("gq jac right composition second 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-12);
    check("gq jac right composition second 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-12);
    check("gq jac right composition second 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-12);
    check("gq jac right composition second 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-12);
    check("gq jac right composition second 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-12);

    speu_dcm SgR_bh        = gR_bh.plus_right(DeltatauH);
    speu_dcm gR_eh         = gR_eb * gR_bh;
    speu_dcm SgR_eh        = gR_eb * SgR_bh;
    Eigen::Matrix6d Jac_gR = gR_eb.jac_right_composition_wrt_second(gR_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauH());
    speu_dcm SSgR_eh       = gR_eh.plus_right(Deltatau_bh_gR);

    check("gR jac right composition second 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-12);
    check("gR jac right composition second 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-12);
    check("gR jac right composition second 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-12);
    check("gR jac right composition second 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-12);
    check("gR jac right composition second 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-12);
    check("gR jac right composition second 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-12);

    homogeneous SM_bh     = M_bh.plus_right(DeltatauH);
    homogeneous M_eh      = M_eb * M_bh;
    homogeneous SM_eh     = M_eb * SM_bh;
    Eigen::Matrix6d Jac_M = M_eb.jac_right_composition_wrt_second(M_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauH());
    homogeneous SSM_eh    = M_eh.plus_right(Deltatau_bh_M);

    check("M jac right composition second 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-12);
    check("M jac right composition second 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-12);
    check("M jac right composition second 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-12);
    check("M jac right composition second 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-12);
    check("M jac right composition second 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-12);
    check("M jac right composition second 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-12);

    dual SZ_bh            = Z_bh.plus_right(DeltatauH);
    dual Z_eh             = Z_eb * Z_bh;
    dual SZ_eh            = Z_eb * SZ_bh;
    Eigen::Matrix6d Jac_Z = Z_eb.jac_right_composition_wrt_second(Z_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauH());
    dual SSZ_eh           = Z_eh.plus_right(Deltatau_bh_Z);

    check("Z jac right composition second 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-12);
    check("Z jac right composition second 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-12);
    check("Z jac right composition second 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-12);
    check("Z jac right composition second 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-12);
    check("Z jac right composition second 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-12);
    check("Z jac right composition second 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-12);

    trfv Stau_bh            = tau_bh.plus_right(DeltatauH);
    trfv tau_eh             = tau_eb * tau_bh;
    trfv Stau_eh            = tau_eb * Stau_bh;
    Eigen::Matrix6d Jac_tau = tau_eb.jac_right_composition_wrt_second(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauH());
    trfv SStau_eh           = tau_eh.plus_right(Deltatau_bh_tau);

    check("tau jac right composition second 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-12);
    check("tau jac right composition second 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-12);
    check("tau jac right composition second 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-12);
    check("tau jac right composition second 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-12);
    check("tau jac right composition second 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-12);
    check("tau jac right composition second 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-12);

} // closes test_jac_right_composition_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_forward_motion_wrt_motion() {

    // It verifies the behavior of the motion classes jac_right_forward_motion_wrt_motion method, validating the following:
    // (M plus DeltatauB) * p = M * p + J * DeltatauB

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    Eigen::Vector3d pB(0.33, -0.65, 0.17);

    Eigen::Vector3d pE_gq      = gq_eb * pB;
    speu_rodrigues Sgq_eb      = gq_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpE_gq     = Sgq_eb * pB;
    Eigen::Matrix36d Jac_gq    = gq_eb.jac_right_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_gq = Jac_gq * DeltatauB();
    Eigen::Vector3d SSpE_gq    = pE_gq + DeltapE_gq;

    check("gq jac right forward motion motion 1 taylor    ", SpE_gq(0), SSpE_gq(0), 1e-5);
    check("gq jac right forward motion motion 2 taylor    ", SpE_gq(1), SSpE_gq(1), 1e-4);
    check("gq jac right forward motion motion 3 taylor    ", SpE_gq(2), SSpE_gq(2), 1e-4);

    Eigen::Vector3d pE_gR      = gR_eb * pB;
    speu_dcm SgR_eb            = gR_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpE_gR     = SgR_eb * pB;
    Eigen::Matrix36d Jac_gR    = gR_eb.jac_right_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_gR = Jac_gR * DeltatauB();
    Eigen::Vector3d SSpE_gR    = pE_gR + DeltapE_gR;

    check("gR jac right forward motion motion 1 taylor    ", SpE_gR(0), SSpE_gR(0), 1e-5);
    check("gR jac right forward motion motion 2 taylor    ", SpE_gR(1), SSpE_gR(1), 1e-4);
    check("gR jac right forward motion motion 3 taylor    ", SpE_gR(2), SSpE_gR(2), 1e-4);

    Eigen::Vector3d pE_M      = M_eb * pB;
    homogeneous SM_eb         = M_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpE_M     = SM_eb * pB;
    Eigen::Matrix36d Jac_M    = M_eb.jac_right_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_M = Jac_M * DeltatauB();
    Eigen::Vector3d SSpE_M    = pE_M + DeltapE_M;

    check("M jac right forward motion motion 1 taylor     ", SpE_M(0), SSpE_M(0), 1e-5);
    check("M jac right forward motion motion 2 taylor     ", SpE_M(1), SSpE_M(1), 1e-4);
    check("M jac right forward motion motion 3 taylor     ", SpE_M(2), SSpE_M(2), 1e-4);

    Eigen::Vector3d pE_Z      = Z_eb * pB;
    dual SZ_eb                = Z_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpE_Z     = SZ_eb * pB;
    Eigen::Matrix36d Jac_Z    = Z_eb.jac_right_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_Z = Jac_Z * DeltatauB();
    Eigen::Vector3d SSpE_Z    = pE_Z + DeltapE_Z;

    check("Z jac right forward motion motion 1 taylor     ", SpE_Z(0), SSpE_Z(0), 1e-5);
    check("Z jac right forward motion motion 2 taylor     ", SpE_Z(1), SSpE_Z(1), 1e-4);
    check("Z jac right forward motion motion 3 taylor     ", SpE_Z(2), SSpE_Z(2), 1e-4);

    Eigen::Vector3d pE_tau      = tau_eb * pB;
    trfv Stau_eb                = tau_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpE_tau     = Stau_eb * pB;
    Eigen::Matrix36d Jac_tau    = tau_eb.jac_right_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_tau = Jac_tau * DeltatauB();
    Eigen::Vector3d SSpE_tau    = pE_tau + DeltapE_tau;

    check("tau jac right forward motion motion 1 taylor   ", SpE_tau(0), SSpE_tau(0), 1e-5);
    check("tau jac right forward motion motion 2 taylor   ", SpE_tau(1), SSpE_tau(1), 1e-4);
    check("tau jac right forward motion motion 3 taylor   ", SpE_tau(2), SSpE_tau(2), 1e-4);

} // closes test_jac_right_forward_motion_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_backward_motion_wrt_motion() {

    // It verifies the behavior of the motion classes jac_right_backward_motion_wrt_motion method, validating the following:
    // (M plus DeltatauB) / p = M / p + J * DeltatauB

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    Eigen::Vector3d pE(0.33, -0.65, 0.17);

    Eigen::Vector3d pB_gq      = gq_eb / pE;
    speu_rodrigues Sgq_eb      = gq_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpB_gq     = Sgq_eb / pE;
    Eigen::Matrix36d Jac_gq    = gq_eb.jac_right_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_gq = Jac_gq * DeltatauB();
    Eigen::Vector3d SSpB_gq    = pB_gq + DeltapB_gq;

    check("gq jac right backward motion motion 1 taylor    ", SpB_gq(0), SSpB_gq(0), 1e-4);
    check("gq jac right backward motion motion 2 taylor    ", SpB_gq(1), SSpB_gq(1), 1e-4);
    check("gq jac right backward motion motion 3 taylor    ", SpB_gq(2), SSpB_gq(2), 1e-4);

    Eigen::Vector3d pB_gR      = gR_eb / pE;
    speu_dcm SgR_eb            = gR_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpB_gR     = SgR_eb / pE;
    Eigen::Matrix36d Jac_gR    = gR_eb.jac_right_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_gR = Jac_gR * DeltatauB();
    Eigen::Vector3d SSpB_gR    = pB_gR + DeltapB_gR;

    check("gR jac right backward motion motion 1 taylor    ", SpB_gR(0), SSpB_gR(0), 1e-4);
    check("gR jac right backward motion motion 2 taylor    ", SpB_gR(1), SSpB_gR(1), 1e-4);
    check("gR jac right backward motion motion 3 taylor    ", SpB_gR(2), SSpB_gR(2), 1e-4);

    Eigen::Vector3d pB_M      = M_eb / pE;
    homogeneous SM_eb         = M_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpB_M     = SM_eb / pE;
    Eigen::Matrix36d Jac_M    = M_eb.jac_right_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_M = Jac_M * DeltatauB();
    Eigen::Vector3d SSpB_M    = pB_M + DeltapB_M;

    check("M jac right backward motion motion 1 taylor     ", SpB_M(0), SSpB_M(0), 1e-4);
    check("M jac right backward motion motion 2 taylor     ", SpB_M(1), SSpB_M(1), 1e-4);
    check("M jac right backward motion motion 3 taylor     ", SpB_M(2), SSpB_M(2), 1e-4);

    Eigen::Vector3d pB_Z      = Z_eb / pE;
    dual SZ_eb                = Z_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpB_Z     = SZ_eb / pE;
    Eigen::Matrix36d Jac_Z    = Z_eb.jac_right_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_Z = Jac_Z * DeltatauB();
    Eigen::Vector3d SSpB_Z    = pB_Z + DeltapB_Z;

    check("Z jac right backward motion motion 1 taylor     ", SpB_Z(0), SSpB_Z(0), 1e-4);
    check("Z jac right backward motion motion 2 taylor     ", SpB_Z(1), SSpB_Z(1), 1e-4);
    check("Z jac right backward motion motion 3 taylor     ", SpB_Z(2), SSpB_Z(2), 1e-4);

    Eigen::Vector3d pB_tau      = tau_eb / pE;
    trfv Stau_eb                = tau_eb.plus_right(DeltatauB);
    Eigen::Vector3d SpB_tau     = Stau_eb / pE;
    Eigen::Matrix36d Jac_tau    = tau_eb.jac_right_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_tau = Jac_tau * DeltatauB();
    Eigen::Vector3d SSpB_tau    = pB_tau + DeltapB_tau;

    check("tau jac right backward motion motion 1 taylor   ", SpB_tau(0), SSpB_tau(0), 1e-4);
    check("tau jac right backward motion motion 2 taylor   ", SpB_tau(1), SSpB_tau(1), 1e-4);
    check("tau jac right backward motion motion 3 taylor   ", SpB_tau(2), SSpB_tau(2), 1e-4);

} // closes test_jac_right_backward_motion_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_log() {

    // It checks the proper behavior of the "rotv" right logarithmic jacobian method,
    // validating the following equations:
    // Log(gR plus Deltatau) ~= Log(gR) + (J * Deltatau)
    // Log(gq plus Deltatau) ~= Log(gq) + (J * Deltatau)
    // Log(M plus Deltatau) ~= Log(M) + (J * Deltatau)
    // Log(zeta plus Deltatau) ~= Log(zeta) + (J * Deltatau)
    // Log(tau plus Deltatau) ~= Log(tau) + (J * Deltatau)

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    speu_rodrigues Sgq_eb       = gq_eb.plus_right(DeltatauB);
    trfv Stau_eb_gq             = Sgq_eb.log_map_trfv();
    trfv tau_eb_gq              = gq_eb.log_map_trfv();
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_right_log();
    trfv Deltatau_eb_gq         = trfv(Jac_gq * DeltatauB());
    trfv SStau_eb_gq(tau_eb_gq() + Deltatau_eb_gq());

    check("gq jac right log 1 taylor    ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-5);
    check("gq jac right log 2 taylor    ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-5);
    check("gq jac right log 3 taylor    ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-5);
    check("gq jac right log 4 taylor    ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-6);
    check("gq jac right log 5 taylor    ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-6);
    check("gq jac right log 6 taylor    ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-6);

    speu_dcm SgR_eb             = gR_eb.plus_right(DeltatauB);
    trfv Stau_eb_gR             = SgR_eb.log_map_trfv();
    trfv tau_eb_gR              = gR_eb.log_map_trfv();
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_right_log();
    trfv Deltatau_eb_gR         = trfv(Jac_gR * DeltatauB());
    trfv SStau_eb_gR(tau_eb_gR() + Deltatau_eb_gR());

    check("gR jac right log 1 taylor    ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-5);
    check("gR jac right log 2 taylor    ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-5);
    check("gR jac right log 3 taylor    ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-5);
    check("gR jac right log 4 taylor    ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-6);
    check("gR jac right log 5 taylor    ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-6);
    check("gR jac right log 6 taylor    ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-6);

    homogeneous SM_eb          = M_eb.plus_right(DeltatauB);
    trfv Stau_eb_M             = SM_eb.log_map_trfv();
    trfv tau_eb_M              = M_eb.log_map_trfv();
    Eigen::Matrix6d Jac_M      = M_eb.jac_right_log();
    trfv Deltatau_eb_M         = trfv(Jac_M * DeltatauB());
    trfv SStau_eb_M(tau_eb_M() + Deltatau_eb_M());

    check("M jac right log 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-5);
    check("M jac right log 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-5);
    check("M jac right log 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-5);
    check("M jac right log 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-6);
    check("M jac right log 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-6);
    check("M jac right log 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-6);

    dual SZ_eb                 = Z_eb.plus_right(DeltatauB);
    trfv Stau_eb_Z             = SZ_eb.log_map_trfv();
    trfv tau_eb_Z              = Z_eb.log_map_trfv();
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_right_log();
    trfv Deltatau_eb_Z         = trfv(Jac_Z * DeltatauB());
    trfv SStau_eb_Z(tau_eb_Z() + Deltatau_eb_Z());

    check("Z jac right log 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-5);
    check("Z jac right log 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-5);
    check("Z jac right log 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-5);
    check("Z jac right log 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-6);
    check("Z jac right log 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-6);
    check("Z jac right log 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-6);

    trfv Stau_eb                 = tau_eb.plus_right(DeltatauB);
    trfv Stau_eb_tau             = Stau_eb;
    trfv tau_eb_tau              = tau_eb;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_right_log();
    trfv Deltatau_eb_tau         = trfv(Jac_tau * DeltatauB());
    trfv SStau_eb_tau(tau_eb_tau() + Deltatau_eb_tau());

    check("tau jac right log 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-5);
    check("tau jac right log 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-5);
    check("tau jac right log 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-5);
    check("tau jac right log 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-6);
    check("tau jac right log 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-6);
    check("tau jac right log 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-6);

} // closes test_jacobian_right_log

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_plus_wrt_first() {

    // It verifies the behavior of the rotation classes jac_right_plus_wrt_first method, validating the following:
    // (M1 plus Deltatau1B) plus tau2 = (M1 plus tau2) plus J * Deltatau1B

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);
    trfv tau_bh(r_bh, T_bhb);

    speu_rodrigues gq_eh   = gq_eb.plus_right(tau_bh);
    speu_rodrigues Sgq_eb  = gq_eb.plus_right(DeltatauB);
    speu_rodrigues Sgq_eh  = Sgq_eb.plus_right(tau_bh);
    Eigen::Matrix6d Jac_gq = gq_eb.jac_right_plus_wrt_first(tau_bh);
    trfv Deltatau_bh_gq    = ang::trfv(Jac_gq * DeltatauB());
    speu_rodrigues SSgq_eh = gq_eh.plus_right(Deltatau_bh_gq);

    check("gq jac right plus first 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-12);
    check("gq jac right plus first 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-12);
    check("gq jac right plus first 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-12);
    check("gq jac right plus first 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-12);
    check("gq jac right plus first 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-12);
    check("gq jac right plus first 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-12);

    speu_dcm gR_eh         = gR_eb.plus_right(tau_bh);
    speu_dcm SgR_eb        = gR_eb.plus_right(DeltatauB);
    speu_dcm SgR_eh        = SgR_eb.plus_right(tau_bh);
    Eigen::Matrix6d Jac_gR = gR_eb.jac_right_plus_wrt_first(tau_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauB());
    speu_dcm SSgR_eh       = gR_eh.plus_right(Deltatau_bh_gR);

    check("gR jac right plus first 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-12);
    check("gR jac right plus first 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-12);
    check("gR jac right plus first 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-12);
    check("gR jac right plus first 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-12);
    check("gR jac right plus first 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-12);
    check("gR jac right plus first 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-12);

    homogeneous M_eh      = M_eb.plus_right(tau_bh);
    homogeneous SM_eb     = M_eb.plus_right(DeltatauB);
    homogeneous SM_eh     = SM_eb.plus_right(tau_bh);
    Eigen::Matrix6d Jac_M = M_eb.jac_right_plus_wrt_first(tau_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauB());
    homogeneous SSM_eh    = M_eh.plus_right(Deltatau_bh_M);

    check("M jac right plus first 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-12);
    check("M jac right plus first 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-12);
    check("M jac right plus first 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-12);
    check("M jac right plus first 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-12);
    check("M jac right plus first 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-12);
    check("M jac right plus first 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-12);

    dual Z_eh             = Z_eb.plus_right(tau_bh);
    dual SZ_eb            = Z_eb.plus_right(DeltatauB);
    dual SZ_eh            = SZ_eb.plus_right(tau_bh);
    Eigen::Matrix6d Jac_Z = Z_eb.jac_right_plus_wrt_first(tau_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauB());
    dual SSZ_eh           = Z_eh.plus_right(Deltatau_bh_Z);

    check("Z jac right plus first 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-12);
    check("Z jac right plus first 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-12);
    check("Z jac right plus first 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-12);
    check("Z jac right plus first 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-12);
    check("Z jac right plus first 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-12);
    check("Z jac right plus first 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-12);

    trfv tau_eh             = tau_eb.plus_right(tau_bh);
    trfv Stau_eb            = tau_eb.plus_right(DeltatauB);
    trfv Stau_eh            = Stau_eb.plus_right(tau_bh);
    Eigen::Matrix6d Jac_tau = tau_eb.jac_right_plus_wrt_first(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauB());
    trfv SStau_eh           = tau_eh.plus_right(Deltatau_bh_tau);

    check("tau jac right plus first 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-12);
    check("tau jac right plus first 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-12);
    check("tau jac right plus first 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-12);
    check("tau jac right plus first 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-12);
    check("tau jac right plus first 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-12);
    check("tau jac right plus first 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-12);

} // closes test_jac_right_plus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_plus_wrt_second() {

    // It verifies the behavior of the motion classes jac_right_plus_wrt_second method, validating the following:
    // M1 plus (tau2 + Deltatau2B) = (M1 plus tau2) plus J * Deltatau2B

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauH(v);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);
    trfv tau_bh(r_bh, T_bhb);

    trfv tautau_bh(tau_bh() + DeltatauH());

    speu_rodrigues gq_eh   = gq_eb.plus_right(tau_bh);
    speu_rodrigues Sgq_eh  = gq_eb.plus_right(tautau_bh);
    Eigen::Matrix6d Jac_gq = gq_eb.jac_right_plus_wrt_second(tau_bh);
    trfv Deltatau_bh_gq    = ang::trfv(Jac_gq * DeltatauH());
    speu_rodrigues SSgq_eh = gq_eh.plus_right(Deltatau_bh_gq);

    check("gq jac right plus second 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-6);
    check("gq jac right plus second 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-6);
    check("gq jac right plus second 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-6);
    check("gq jac right plus second 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-5);
    check("gq jac right plus second 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-5);
    check("gq jac right plus second 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-5);

    speu_dcm gR_eh         = gR_eb.plus_right(tau_bh);
    speu_dcm SgR_eh        = gR_eb.plus_right(tautau_bh);
    Eigen::Matrix6d Jac_gR = gR_eb.jac_right_plus_wrt_second(tau_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauH());
    speu_dcm SSgR_eh       = gR_eh.plus_right(Deltatau_bh_gR);

    check("gR jac right plus second 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-6);
    check("gR jac right plus second 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-6);
    check("gR jac right plus second 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-6);
    check("gR jac right plus second 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-5);
    check("gR jac right plus second 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-5);
    check("gR jac right plus second 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-5);

    homogeneous M_eh      = M_eb.plus_right(tau_bh);
    homogeneous SM_eh     = M_eb.plus_right(tautau_bh);
    Eigen::Matrix6d Jac_M = M_eb.jac_right_plus_wrt_second(tau_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauH());
    homogeneous SSM_eh    = M_eh.plus_right(Deltatau_bh_M);

    check("M jac right plus second 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-6);
    check("M jac right plus second 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-6);
    check("M jac right plus second 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-6);
    check("M jac right plus second 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-5);
    check("M jac right plus second 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-5);
    check("M jac right plus second 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-5);

    dual Z_eh             = Z_eb.plus_right(tau_bh);
    dual SZ_eh            = Z_eb.plus_right(tautau_bh);
    Eigen::Matrix6d Jac_Z = Z_eb.jac_right_plus_wrt_second(tau_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauH());
    dual SSZ_eh           = Z_eh.plus_right(Deltatau_bh_Z);

    check("Z jac right plus second 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-6);
    check("Z jac right plus second 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-6);
    check("Z jac right plus second 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-6);
    check("Z jac right plus second 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-5);
    check("Z jac right plus second 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-5);
    check("Z jac right plus second 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-5);

    trfv tau_eh             = tau_eb.plus_right(tau_bh);
    trfv Stau_eh            = tau_eb.plus_right(tautau_bh);
    Eigen::Matrix6d Jac_tau = tau_eb.jac_right_plus_wrt_second(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauH());
    trfv SStau_eh           = tau_eh.plus_right(Deltatau_bh_tau);

    check("tau jac right plus second 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-6);
    check("tau jac right plus second 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-6);
    check("tau jac right plus second 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-6);
    check("tau jac right plus second 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-5);
    check("tau jac right plus second 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-5);
    check("tau jac right plus second 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-5);

} // closes test_jac_right_plus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_minus_wrt_first() {

    // It verifies the behavior of the motion classes jac_right_minus_wrt_first method, validating the following:
    // (M2 plus Deltatau2B) minus M1 = (M2 minus M1) + J * Deltatau2B

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    trfv Ptau_eb          = tau_eb.plus_right(DeltatauB);
    speu_dcm PgR_eb       = gR_eb.plus_right(DeltatauB);
    speu_rodrigues Pgq_eb = gq_eb.plus_right(DeltatauB);
    homogeneous PM_eb     = M_eb.plus_right(DeltatauB);
    dual PZ_eb            = Z_eb.plus_right(DeltatauB);

    trfv DeltatauB_gq      = Pgq_eb.minus_right_trfv(gq_eb);
    speu_rodrigues Sgq_eb  = Pgq_eb.plus_right(DeltatauB);
    trfv Stau_eb_gq        = Sgq_eb.minus_right_trfv(gq_eb);
    Eigen::Matrix6d Jac_gq = Pgq_eb.jac_right_minus_wrt_first(gq_eb);
    trfv Deltatau_eb_gq    = ang::trfv(Jac_gq * DeltatauB());
    trfv SStau_eb_gq(DeltatauB_gq() + Deltatau_eb_gq());

    check("gq jac right minus first 1 taylor   ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-7);
    check("gq jac right minus first 2 taylor   ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-7);
    check("gq jac right minus first 3 taylor   ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-7);
    check("gq jac right minus first 4 taylor   ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-7);
    check("gq jac right minus first 5 taylor   ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-7);
    check("gq jac right minus first 6 taylor   ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-7);

    trfv DeltatauB_gR      = PgR_eb.minus_right_trfv(gR_eb);
    speu_dcm SgR_eb        = PgR_eb.plus_right(DeltatauB);
    trfv Stau_eb_gR        = SgR_eb.minus_right_trfv(gR_eb);
    Eigen::Matrix6d Jac_gR = PgR_eb.jac_right_minus_wrt_first(gR_eb);
    trfv Deltatau_eb_gR    = ang::trfv(Jac_gR * DeltatauB());
    trfv SStau_eb_gR(DeltatauB_gR() + Deltatau_eb_gR());

    check("gR jac right minus first 1 taylor   ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-7);
    check("gR jac right minus first 2 taylor   ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-7);
    check("gR jac right minus first 3 taylor   ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-7);
    check("gR jac right minus first 4 taylor   ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-7);
    check("gR jac right minus first 5 taylor   ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-7);
    check("gR jac right minus first 6 taylor   ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-7);

    trfv DeltatauB_M      = PM_eb.minus_right_trfv(M_eb);
    homogeneous SM_eb     = PM_eb.plus_right(DeltatauB);
    trfv Stau_eb_M        = SM_eb.minus_right_trfv(M_eb);
    Eigen::Matrix6d Jac_M = PM_eb.jac_right_minus_wrt_first(M_eb);
    trfv Deltatau_eb_M    = ang::trfv(Jac_M * DeltatauB());
    trfv SStau_eb_M(DeltatauB_M() + Deltatau_eb_M());

    check("M jac right minus first 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-7);
    check("M jac right minus first 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-7);
    check("M jac right minus first 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-7);
    check("M jac right minus first 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-7);
    check("M jac right minus first 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-7);
    check("M jac right minus first 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-7);

    trfv DeltatauB_Z      = PZ_eb.minus_right_trfv(Z_eb);
    dual SZ_eb            = PZ_eb.plus_right(DeltatauB);
    trfv Stau_eb_Z        = SZ_eb.minus_right_trfv(Z_eb);
    Eigen::Matrix6d Jac_Z = PZ_eb.jac_right_minus_wrt_first(Z_eb);
    trfv Deltatau_eb_Z    = ang::trfv(Jac_Z * DeltatauB());
    trfv SStau_eb_Z(DeltatauB_Z() + Deltatau_eb_Z());

    check("Z jac right minus first 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-7);
    check("Z jac right minus first 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-7);
    check("Z jac right minus first 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-7);
    check("Z jac right minus first 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-7);
    check("Z jac right minus first 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-7);
    check("Z jac right minus first 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-7);

    trfv DeltatauB_tau      = Ptau_eb.minus_right_trfv(tau_eb);
    trfv Stau_eb            = Ptau_eb.plus_right(DeltatauB);
    trfv Stau_eb_tau        = Stau_eb.minus_right_trfv(tau_eb);
    Eigen::Matrix6d Jac_tau = Ptau_eb.jac_right_minus_wrt_first(tau_eb);
    trfv Deltatau_eb_tau    = ang::trfv(Jac_tau * DeltatauB());
    trfv SStau_eb_tau(DeltatauB_tau() + Deltatau_eb_tau());

    check("tau jac right minus first 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-7);
    check("tau jac right minus first 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-7);
    check("tau jac right minus first 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-7);
    check("tau jac right minus first 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-7);
    check("tau jac right minus first 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-7);
    check("tau jac right minus first 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-7);

} // closes test_jac_right_minus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_minus_wrt_second() {

    // It verifies the behavior of the motion classes jac_right_minus_wrt_second method, validating the following:
    // M2 minus(M1 plus Deltatau1B) = (M2 minus M1) + J * Deltatau1B

    rotv r1_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T1_ebe(7.5, -12.8, -3.2);

    trfv tau1_eb(r1_eb, T1_ebe);
    speu_dcm gR1_eb(tau1_eb);
    speu_rodrigues gq1_eb(tau1_eb);
    homogeneous M1_eb(tau1_eb);
    dual Z1_eb(tau1_eb);

    Eigen::Vector6d vA; vA << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauBA(vA);

    trfv tau2_eb          = tau1_eb.plus_right(DeltatauBA);
    speu_dcm gR2_eb       = gR1_eb.plus_right(DeltatauBA);
    speu_rodrigues gq2_eb = gq1_eb.plus_right(DeltatauBA);
    homogeneous M2_eb     = M1_eb.plus_right(DeltatauBA);
    dual Z2_eb            = Z1_eb.plus_right(DeltatauBA);

    Eigen::Vector6d v; v << -0.002, -0.003, -0.002, -1.5e-3, -9e-4, 2.8e-3;
    trfv DeltatauB(v);

    trfv DeltatauB_gq        = gq2_eb.minus_right_trfv(gq1_eb);
    speu_rodrigues Sgq_eb    = gq1_eb.plus_right(DeltatauB);
    trfv Stau_eb_gq          = gq2_eb.minus_right_trfv(Sgq_eb);
    Eigen::Matrix6d Jac_gq   = gq2_eb.jac_right_minus_wrt_second(gq1_eb);
    trfv Deltatau_eb_gq      = ang::trfv(Jac_gq * DeltatauB());
    trfv SStau_eb_gq(DeltatauB_gq() + Deltatau_eb_gq());

    check("gq jac right minus second 1 taylor   ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-7);
    check("gq jac right minus second 2 taylor   ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-7);
    check("gq jac right minus second 3 taylor   ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-7);
    check("gq jac right minus second 4 taylor   ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-7);
    check("gq jac right minus second 5 taylor   ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-7);
    check("gq jac right minus second 6 taylor   ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-7);

    trfv DeltatauB_gR        = gR2_eb.minus_right_trfv(gR1_eb);
    speu_dcm SgR_eb             = gR1_eb.plus_right(DeltatauB);
    trfv Stau_eb_gR          = gR2_eb.minus_right_trfv(SgR_eb);
    Eigen::Matrix6d Jac_gR = gR2_eb.jac_right_minus_wrt_second(gR1_eb);
    trfv Deltatau_eb_gR      = ang::trfv(Jac_gR * DeltatauB());
    trfv SStau_eb_gR(DeltatauB_gR() + Deltatau_eb_gR());

    check("gR jac right minus second 1 taylor   ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-7);
    check("gR jac right minus second 2 taylor   ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-7);
    check("gR jac right minus second 3 taylor   ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-7);
    check("gR jac right minus second 4 taylor   ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-7);
    check("gR jac right minus second 5 taylor   ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-7);
    check("gR jac right minus second 6 taylor   ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-7);

    trfv DeltatauB_M      = M2_eb.minus_right_trfv(M1_eb);
    homogeneous SM_eb     = M1_eb.plus_right(DeltatauB);
    trfv Stau_eb_M        = M2_eb.minus_right_trfv(SM_eb);
    Eigen::Matrix6d Jac_M = M2_eb.jac_right_minus_wrt_second(M1_eb);
    trfv Deltatau_eb_M    = ang::trfv(Jac_M * DeltatauB());
    trfv SStau_eb_M(DeltatauB_M() + Deltatau_eb_M());

    check("M jac right minus second 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-7);
    check("M jac right minus second 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-7);
    check("M jac right minus second 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-7);
    check("M jac right minus second 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-7);
    check("M jac right minus second 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-7);
    check("M jac right minus second 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-7);

    trfv DeltatauB_Z      = Z2_eb.minus_right_trfv(Z1_eb);
    dual SZ_eb            = Z1_eb.plus_right(DeltatauB);
    trfv Stau_eb_Z        = Z2_eb.minus_right_trfv(SZ_eb);
    Eigen::Matrix6d Jac_Z = Z2_eb.jac_right_minus_wrt_second(Z1_eb);
    trfv Deltatau_eb_Z    = ang::trfv(Jac_Z * DeltatauB());
    trfv SStau_eb_Z(DeltatauB_Z() + Deltatau_eb_Z());

    check("Z jac right minus second 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-7);
    check("Z jac right minus second 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-7);
    check("Z jac right minus second 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-7);
    check("Z jac right minus second 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-7);
    check("Z jac right minus second 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-7);
    check("Z jac right minus second 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-7);

    trfv DeltatauB_tau      = tau2_eb.minus_right_trfv(tau1_eb);
    trfv Stau_eb            = tau1_eb.plus_right(DeltatauB);
    trfv Stau_eb_tau        = tau2_eb.minus_right_trfv(Stau_eb);
    Eigen::Matrix6d Jac_tau = tau2_eb.jac_right_minus_wrt_second(tau1_eb);
    trfv Deltatau_eb_tau    = ang::trfv(Jac_tau * DeltatauB());
    trfv SStau_eb_tau(DeltatauB_tau() + Deltatau_eb_tau());

    check("tau jac right minus second 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-7);
    check("tau jac right minus second 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-7);
    check("tau jac right minus second 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-7);
    check("tau jac right minus second 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-7);
    check("tau jac right minus second 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-7);
    check("tau jac right minus second 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-7);

} // closes test_jac_right_minus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_forward_adjoint_wrt_motion() {

    // It verifies the behavior of the motion classes jac_right_forward_adjoint_wrt_motion method, validating the following:
    // Ad(M plus DeltatauB) | xi = AdM | xi + J * DeltatauB *

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d tau; tau << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(tau);

    se3_tangent xiB(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    se3_tangent xiE_gq          = gq_eb | xiB;
    speu_rodrigues Sgq_eb       = gq_eb.plus_right(DeltatauB);
    se3_tangent SxiE_gq         = Sgq_eb | xiB;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_right_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_gq = Jac_gq * DeltatauB();
    se3_tangent SSxiE_gq        = xiE_gq + DeltaxiE_gq;

    check("gq jac right forward adjoint motion 1 taylor    ", SxiE_gq()(0), SSxiE_gq()(0), 1e-4);
    check("gq jac right forward adjoint motion 2 taylor    ", SxiE_gq()(1), SSxiE_gq()(1), 1e-4);
    check("gq jac right forward adjoint motion 3 taylor    ", SxiE_gq()(2), SSxiE_gq()(2), 1e-4);
    check("gq jac right forward adjoint motion 4 taylor    ", SxiE_gq()(3), SSxiE_gq()(3), 1e-4);
    check("gq jac right forward adjoint motion 5 taylor    ", SxiE_gq()(4), SSxiE_gq()(4), 1e-4);
    check("gq jac right forward adjoint motion 6 taylor    ", SxiE_gq()(5), SSxiE_gq()(5), 1e-4);

    se3_tangent xiE_gR          = gR_eb | xiB;
    speu_dcm SgR_eb             = gR_eb.plus_right(DeltatauB);
    se3_tangent SxiE_gR         = SgR_eb | xiB;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_right_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_gR = Jac_gR * DeltatauB();
    se3_tangent SSxiE_gR        = xiE_gR + DeltaxiE_gR;

    check("gR jac right forward adjoint motion 1 taylor    ", SxiE_gR()(0), SSxiE_gR()(0), 1e-4);
    check("gR jac right forward adjoint motion 2 taylor    ", SxiE_gR()(1), SSxiE_gR()(1), 1e-4);
    check("gR jac right forward adjoint motion 3 taylor    ", SxiE_gR()(2), SSxiE_gR()(2), 1e-4);
    check("gR jac right forward adjoint motion 4 taylor    ", SxiE_gR()(3), SSxiE_gR()(3), 1e-4);
    check("gR jac right forward adjoint motion 5 taylor    ", SxiE_gR()(4), SSxiE_gR()(4), 1e-4);
    check("gR jac right forward adjoint motion 6 taylor    ", SxiE_gR()(5), SSxiE_gR()(5), 1e-4);

    se3_tangent xiE_M          = M_eb | xiB;
    homogeneous SM_eb          = M_eb.plus_right(DeltatauB);
    se3_tangent SxiE_M         = SM_eb | xiB;
    Eigen::Matrix6d Jac_M      = M_eb.jac_right_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_M = Jac_M * DeltatauB();
    se3_tangent SSxiE_M        = xiE_M + DeltaxiE_M;

    check("M jac right forward adjoint motion 1 taylor     ", SxiE_M()(0), SSxiE_M()(0), 1e-4);
    check("M jac right forward adjoint motion 2 taylor     ", SxiE_M()(1), SSxiE_M()(1), 1e-4);
    check("M jac right forward adjoint motion 3 taylor     ", SxiE_M()(2), SSxiE_M()(2), 1e-4);
    check("M jac right forward adjoint motion 4 taylor     ", SxiE_M()(3), SSxiE_M()(3), 1e-4);
    check("M jac right forward adjoint motion 5 taylor     ", SxiE_M()(4), SSxiE_M()(4), 1e-4);
    check("M jac right forward adjoint motion 6 taylor     ", SxiE_M()(5), SSxiE_M()(5), 1e-4);
    se3_tangent xiE_Z          = Z_eb | xiB;
    dual SZ_eb                 = Z_eb.plus_right(DeltatauB);
    se3_tangent SxiE_Z         = SZ_eb | xiB;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_right_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_Z = Jac_Z * DeltatauB();
    se3_tangent SSxiE_Z        = xiE_Z + DeltaxiE_Z;

    check("Z jac right forward adjoint motion 1 taylor     ", SxiE_Z()(0), SSxiE_Z()(0), 1e-4);
    check("Z jac right forward adjoint motion 2 taylor     ", SxiE_Z()(1), SSxiE_Z()(1), 1e-4);
    check("Z jac right forward adjoint motion 3 taylor     ", SxiE_Z()(2), SSxiE_Z()(2), 1e-4);
    check("Z jac right forward adjoint motion 4 taylor     ", SxiE_Z()(3), SSxiE_Z()(3), 1e-4);
    check("Z jac right forward adjoint motion 5 taylor     ", SxiE_Z()(4), SSxiE_Z()(4), 1e-4);
    check("Z jac right forward adjoint motion 6 taylor     ", SxiE_Z()(5), SSxiE_Z()(5), 1e-4);

    se3_tangent xiE_tau          = tau_eb | xiB;
    trfv Stau_eb                 = tau_eb.plus_right(DeltatauB);
    se3_tangent SxiE_tau         = Stau_eb | xiB;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_right_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_tau = Jac_tau * DeltatauB();
    se3_tangent SSxiE_tau        = xiE_tau + DeltaxiE_tau;

    check("tau jac right forward adjoint motion 1 taylor   ", SxiE_tau()(0), SSxiE_tau()(0), 1e-4);
    check("tau jac right forward adjoint motion 2 taylor   ", SxiE_tau()(1), SSxiE_tau()(1), 1e-4);
    check("tau jac right forward adjoint motion 3 taylor   ", SxiE_tau()(2), SSxiE_tau()(2), 1e-4);
    check("tau jac right forward adjoint motion 4 taylor   ", SxiE_tau()(3), SSxiE_tau()(3), 1e-4);
    check("tau jac right forward adjoint motion 5 taylor   ", SxiE_tau()(4), SSxiE_tau()(4), 1e-4);
    check("tau jac right forward adjoint motion 6 taylor   ", SxiE_tau()(5), SSxiE_tau()(5), 1e-4);

} // closes test_jac_right_forward_adjoint_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_right_backward_adjoint_wrt_motion() {

    // It verifies the behavior of the motion classes jac_right_backward_adjoint_wrt_motion method, validating the following:
    // Ad(M plus DeltatauB) % xi = AdM % xi + J * DeltatauB *

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d tau; tau << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(tau);

    se3_tangent xiE(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    se3_tangent xiB_gq          = gq_eb % xiE;
    speu_rodrigues Sgq_eb       = gq_eb.plus_right(DeltatauB);
    se3_tangent SxiB_gq         = Sgq_eb % xiE;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_right_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_gq = Jac_gq * DeltatauB();
    se3_tangent SSxiB_gq        = xiB_gq + DeltaxiB_gq;

    check("gq jac right backward adjoint motion 1 taylor    ", SxiB_gq()(0), SSxiB_gq()(0), 1e-4);
    check("gq jac right backward adjoint motion 2 taylor    ", SxiB_gq()(1), SSxiB_gq()(1), 1e-4);
    check("gq jac right backward adjoint motion 3 taylor    ", SxiB_gq()(2), SSxiB_gq()(2), 1e-4);
    check("gq jac right backward adjoint motion 4 taylor    ", SxiB_gq()(3), SSxiB_gq()(3), 1e-4);
    check("gq jac right backward adjoint motion 5 taylor    ", SxiB_gq()(4), SSxiB_gq()(4), 1e-4);
    check("gq jac right backward adjoint motion 6 taylor    ", SxiB_gq()(5), SSxiB_gq()(5), 1e-4);

    se3_tangent xiB_gR          = gR_eb % xiE;
    speu_dcm SgR_eb             = gR_eb.plus_right(DeltatauB);
    se3_tangent SxiB_gR         = SgR_eb % xiE;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_right_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_gR = Jac_gR * DeltatauB();
    se3_tangent SSxiB_gR        = xiB_gR + DeltaxiB_gR;

    check("gR jac right backward adjoint motion 1 taylor    ", SxiB_gR()(0), SSxiB_gR()(0), 1e-4);
    check("gR jac right backward adjoint motion 2 taylor    ", SxiB_gR()(1), SSxiB_gR()(1), 1e-4);
    check("gR jac right backward adjoint motion 3 taylor    ", SxiB_gR()(2), SSxiB_gR()(2), 1e-4);
    check("gR jac right backward adjoint motion 4 taylor    ", SxiB_gR()(3), SSxiB_gR()(3), 1e-4);
    check("gR jac right backward adjoint motion 5 taylor    ", SxiB_gR()(4), SSxiB_gR()(4), 1e-4);
    check("gR jac right backward adjoint motion 6 taylor    ", SxiB_gR()(5), SSxiB_gR()(5), 1e-4);

    se3_tangent xiB_M          = M_eb % xiE;
    homogeneous SM_eb          = M_eb.plus_right(DeltatauB);
    se3_tangent SxiB_M         = SM_eb % xiE;
    Eigen::Matrix6d Jac_M      = M_eb.jac_right_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_M = Jac_M * DeltatauB();
    se3_tangent SSxiB_M        = xiB_M + DeltaxiB_M;

    check("M jac right backward adjoint motion 1 taylor     ", SxiB_M()(0), SSxiB_M()(0), 1e-4);
    check("M jac right backward adjoint motion 2 taylor     ", SxiB_M()(1), SSxiB_M()(1), 1e-4);
    check("M jac right backward adjoint motion 3 taylor     ", SxiB_M()(2), SSxiB_M()(2), 1e-4);
    check("M jac right backward adjoint motion 4 taylor     ", SxiB_M()(3), SSxiB_M()(3), 1e-4);
    check("M jac right backward adjoint motion 5 taylor     ", SxiB_M()(4), SSxiB_M()(4), 1e-4);
    check("M jac right backward adjoint motion 6 taylor     ", SxiB_M()(5), SSxiB_M()(5), 1e-4);
    se3_tangent xiB_Z          = Z_eb % xiE;
    dual SZ_eb                 = Z_eb.plus_right(DeltatauB);
    se3_tangent SxiB_Z         = SZ_eb % xiE;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_right_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_Z = Jac_Z * DeltatauB();
    se3_tangent SSxiB_Z        = xiB_Z + DeltaxiB_Z;

    check("Z jac right backward adjoint motion 1 taylor     ", SxiB_Z()(0), SSxiB_Z()(0), 1e-4);
    check("Z jac right backward adjoint motion 2 taylor     ", SxiB_Z()(1), SSxiB_Z()(1), 1e-4);
    check("Z jac right backward adjoint motion 3 taylor     ", SxiB_Z()(2), SSxiB_Z()(2), 1e-4);
    check("Z jac right backward adjoint motion 4 taylor     ", SxiB_Z()(3), SSxiB_Z()(3), 1e-4);
    check("Z jac right backward adjoint motion 5 taylor     ", SxiB_Z()(4), SSxiB_Z()(4), 1e-4);
    check("Z jac right backward adjoint motion 6 taylor     ", SxiB_Z()(5), SSxiB_Z()(5), 1e-4);

    se3_tangent xiB_tau          = tau_eb % xiE;
    trfv Stau_eb                 = tau_eb.plus_right(DeltatauB);
    se3_tangent SxiB_tau         = Stau_eb % xiE;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_right_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_tau = Jac_tau * DeltatauB();
    se3_tangent SSxiB_tau        = xiB_tau + DeltaxiB_tau;

    check("tau jac right backward adjoint motion 1 taylor   ", SxiB_tau()(0), SSxiB_tau()(0), 1e-4);
    check("tau jac right backward adjoint motion 2 taylor   ", SxiB_tau()(1), SSxiB_tau()(1), 1e-4);
    check("tau jac right backward adjoint motion 3 taylor   ", SxiB_tau()(2), SSxiB_tau()(2), 1e-4);
    check("tau jac right backward adjoint motion 4 taylor   ", SxiB_tau()(3), SSxiB_tau()(3), 1e-4);
    check("tau jac right backward adjoint motion 5 taylor   ", SxiB_tau()(4), SSxiB_tau()(4), 1e-4);
    check("tau jac right backward adjoint motion 6 taylor   ", SxiB_tau()(5), SSxiB_tau()(5), 1e-4);

} // closes test_jac_right_backward_adjoint_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_inverse() {

    // It verifies the behavior of the motion classes jac_right_inverse method, validating the following:
    // (DeltatauE plus M)^-1 = J * DeltatauB plus M^-1

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(v);

    speu_rodrigues Sgq_eb  = gq_eb.plus_left(DeltatauE);
    speu_rodrigues gq_be   = gq_eb.inverse();
    speu_rodrigues Sgq_be  = Sgq_eb.inverse();
    Eigen::Matrix6d Jac_gq = gq_eb.jac_left_inverse();
    trfv Deltatau_be_gq    = trfv(Jac_gq * DeltatauE());
    speu_rodrigues SSgq_be = gq_be.plus_left(Deltatau_be_gq);

    check("gq jac left inverse 1 taylor   ", Sgq_be.get_rotv()()(0), SSgq_be.get_rotv()()(0), 1e-12);
    check("gq jac left inverse 2 taylor   ", Sgq_be.get_rotv()()(1), SSgq_be.get_rotv()()(1), 1e-12);
    check("gq jac left inverse 3 taylor   ", Sgq_be.get_rotv()()(2), SSgq_be.get_rotv()()(2), 1e-12);
    check("gq jac left inverse 4 taylor   ", Sgq_be.get_T()(0), SSgq_be.get_T()(0), 1e-12);
    check("gq jac left inverse 5 taylor   ", Sgq_be.get_T()(1), SSgq_be.get_T()(1), 1e-12);
    check("gq jac left inverse 6 taylor   ", Sgq_be.get_T()(2), SSgq_be.get_T()(2), 1e-12);

    speu_dcm SgR_eb        = gR_eb.plus_left(DeltatauE);
    speu_dcm gR_be         = gR_eb.inverse();
    speu_dcm SgR_be        = SgR_eb.inverse();
    Eigen::Matrix6d Jac_gR = gR_eb.jac_left_inverse();
    trfv Deltatau_be_gR    = trfv(Jac_gR * DeltatauE());
    speu_dcm SSgR_be       = gR_be.plus_left(Deltatau_be_gR);

    check("gR jac left inverse 1 taylor   ", SgR_be.get_rotv()()(0), SSgR_be.get_rotv()()(0), 1e-12);
    check("gR jac left inverse 2 taylor   ", SgR_be.get_rotv()()(1), SSgR_be.get_rotv()()(1), 1e-12);
    check("gR jac left inverse 3 taylor   ", SgR_be.get_rotv()()(2), SSgR_be.get_rotv()()(2), 1e-12);
    check("gR jac left inverse 4 taylor   ", SgR_be.get_T()(0), SSgR_be.get_T()(0), 1e-12);
    check("gR jac left inverse 5 taylor   ", SgR_be.get_T()(1), SSgR_be.get_T()(1), 1e-12);
    check("gR jac left inverse 6 taylor   ", SgR_be.get_T()(2), SSgR_be.get_T()(2), 1e-12);

    homogeneous SM_eb     = M_eb.plus_left(DeltatauE);
    homogeneous M_be      = M_eb.inverse();
    homogeneous SM_be     = SM_eb.inverse();
    Eigen::Matrix6d Jac_M = M_eb.jac_left_inverse();
    trfv Deltatau_be_M    = trfv(Jac_M * DeltatauE());
    homogeneous SSM_be    = M_be.plus_left(Deltatau_be_M);

    check("M jac left inverse 1 taylor    ", SM_be.get_rotv()()(0), SSM_be.get_rotv()()(0), 1e-12);
    check("M jac left inverse 2 taylor    ", SM_be.get_rotv()()(1), SSM_be.get_rotv()()(1), 1e-12);
    check("M jac left inverse 3 taylor    ", SM_be.get_rotv()()(2), SSM_be.get_rotv()()(2), 1e-12);
    check("M jac left inverse 4 taylor    ", SM_be.get_T()(0), SSM_be.get_T()(0), 1e-12);
    check("M jac left inverse 5 taylor    ", SM_be.get_T()(1), SSM_be.get_T()(1), 1e-12);
    check("M jac left inverse 6 taylor    ", SM_be.get_T()(2), SSM_be.get_T()(2), 1e-12);

    dual SZ_eb            = Z_eb.plus_left(DeltatauE);
    dual Z_be             = Z_eb.inverse();
    dual SZ_be            = SZ_eb.inverse();
    Eigen::Matrix6d Jac_Z = Z_eb.jac_left_inverse();
    trfv Deltatau_be_Z    = trfv(Jac_Z * DeltatauE());
    dual SSZ_be           = Z_be.plus_left(Deltatau_be_Z);

    check("Z jac left inverse 1 taylor    ", SZ_be.get_rotv()()(0), SSZ_be.get_rotv()()(0), 1e-12);
    check("Z jac left inverse 2 taylor    ", SZ_be.get_rotv()()(1), SSZ_be.get_rotv()()(1), 1e-12);
    check("Z jac left inverse 3 taylor    ", SZ_be.get_rotv()()(2), SSZ_be.get_rotv()()(2), 1e-12);
    check("Z jac left inverse 4 taylor    ", SZ_be.get_T()(0), SSZ_be.get_T()(0), 1e-12);
    check("Z jac left inverse 5 taylor    ", SZ_be.get_T()(1), SSZ_be.get_T()(1), 1e-12);
    check("Z jac left inverse 6 taylor    ", SZ_be.get_T()(2), SSZ_be.get_T()(2), 1e-12);

    trfv tau_be             = tau_eb.inverse();
    trfv Stau_eb            = tau_eb.plus_left(DeltatauE);
    trfv Stau_be            = Stau_eb.inverse();
    Eigen::Matrix6d Jac_tau = tau_eb.jac_left_inverse();
    trfv Deltatau_be_tau    = trfv(Jac_tau * DeltatauE());
    trfv SStau_be           = tau_be.plus_left(Deltatau_be_tau);

    check("tau jac left inverse 1 taylor  ", Stau_be()(0), SStau_be()(0), 1e-12);
    check("tau jac left inverse 2 taylor  ", Stau_be()(1), SStau_be()(1), 1e-12);
    check("tau jac left inverse 3 taylor  ", Stau_be()(2), SStau_be()(2), 1e-12);
    check("tau jac left inverse 4 taylor  ", Stau_be()(3), SStau_be()(3), 1e-12);
    check("tau jac left inverse 5 taylor  ", Stau_be()(4), SStau_be()(4), 1e-12);
    check("tau jac left inverse 6 taylor  ", Stau_be()(5), SStau_be()(5), 1e-12);

} // closes test_jac_left_inverse

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_composition_wrt_first() {

    // It verifies the behavior of the motion classes jac_left_composition_wrt_first method, validating the following:
    // (DeltatauE1 plus M1) * M2 = J * DeltatauE1 plus M1 * M2

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(v);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);

    trfv tau_bh(r_bh, T_bhb);
    speu_dcm gR_bh(tau_bh);
    speu_rodrigues gq_bh(tau_bh);
    homogeneous M_bh(tau_bh);
    dual Z_bh(tau_bh);

    speu_rodrigues Sgq_eb    = gq_eb.plus_left(DeltatauE);
    speu_rodrigues gq_eh     = gq_eb * gq_bh;
    speu_rodrigues Sgq_eh    = Sgq_eb * gq_bh;
    Eigen::Matrix6d Jac_gq   = gq_eb.jac_left_composition_wrt_first(gq_bh);
    trfv Deltatau_bh_gq      = ang::trfv(Jac_gq * DeltatauE());
    speu_rodrigues SSgq_eh   = gq_eh.plus_left(Deltatau_bh_gq);

    check("gq jac left composition first 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-12);
    check("gq jac left composition first 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-12);
    check("gq jac left composition first 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-12);
    check("gq jac left composition first 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-12);
    check("gq jac left composition first 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-12);
    check("gq jac left composition first 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-12);

    speu_dcm SgR_eb        = gR_eb.plus_left(DeltatauE);
    speu_dcm gR_eh         = gR_eb * gR_bh;
    speu_dcm SgR_eh        = SgR_eb * gR_bh;
    Eigen::Matrix6d Jac_gR = gR_eb.jac_left_composition_wrt_first(gR_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauE());
    speu_dcm SSgR_eh       = gR_eh.plus_left(Deltatau_bh_gR);

    check("gR jac left composition first 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-12);
    check("gR jac left composition first 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-12);
    check("gR jac left composition first 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-12);
    check("gR jac left composition first 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-12);
    check("gR jac left composition first 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-12);
    check("gR jac left composition first 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-12);

    homogeneous SM_eb     = M_eb.plus_left(DeltatauE);
    homogeneous M_eh      = M_eb * M_bh;
    homogeneous SM_eh     = SM_eb * M_bh;
    Eigen::Matrix6d Jac_M = M_eb.jac_left_composition_wrt_first(M_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauE());
    homogeneous SSM_eh    = M_eh.plus_left(Deltatau_bh_M);

    check("M jac left composition first 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-12);
    check("M jac left composition first 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-12);
    check("M jac left composition first 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-12);
    check("M jac left composition first 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-12);
    check("M jac left composition first 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-12);
    check("M jac left composition first 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-12);

    dual SZ_eb            = Z_eb.plus_left(DeltatauE);
    dual Z_eh             = Z_eb * Z_bh;
    dual SZ_eh            = SZ_eb * Z_bh;
    Eigen::Matrix6d Jac_Z = Z_eb.jac_left_composition_wrt_first(Z_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauE());
    dual SSZ_eh           = Z_eh.plus_left(Deltatau_bh_Z);

    check("Z jac left composition first 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-12);
    check("Z jac left composition first 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-12);
    check("Z jac left composition first 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-12);
    check("Z jac left composition first 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-12);
    check("Z jac left composition first 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-12);
    check("Z jac left composition first 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-12);

    trfv Stau_eb            = tau_eb.plus_left(DeltatauE);
    trfv tau_eh             = tau_eb * tau_bh;
    trfv Stau_eh            = Stau_eb * tau_bh;
    Eigen::Matrix6d Jac_tau = tau_eb.jac_left_composition_wrt_first(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauE());
    trfv SStau_eh           = tau_eh.plus_left(Deltatau_bh_tau);

    check("tau jac left composition first 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-12);
    check("tau jac left composition first 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-12);
    check("tau jac left composition first 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-12);
    check("tau jac left composition first 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-12);
    check("tau jac left composition first 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-12);
    check("tau jac left composition first 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-12);

} // closes test_jac_left_composition_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_composition_wrt_second() {

    // It verifies the behavior of the motion classes jac_left_composition_wrt_second method, validating the following:
    // M1 * (DeltatauE2 plus M2) = J * DeltatauE2 plus M1 * M2

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    rotv r_bh(0.4, -1.2, 0.72);
    Eigen::Vector3d T_bhb(-2.8, -3.9, 4.2);

    trfv tau_bh(r_bh, T_bhb);
    speu_dcm gR_bh(tau_bh);
    speu_rodrigues gq_bh(tau_bh);
    homogeneous M_bh(tau_bh);
    dual Z_bh(tau_bh);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauB(v);

    speu_rodrigues Sgq_bh    = gq_bh.plus_left(DeltatauB);
    speu_rodrigues gq_eh     = gq_eb * gq_bh;
    speu_rodrigues Sgq_eh    = gq_eb * Sgq_bh;
    Eigen::Matrix6d Jac_gq   = gq_eb.jac_left_composition_wrt_second(gq_bh);
    trfv Deltatau_bh_gq      = ang::trfv(Jac_gq * DeltatauB());
    speu_rodrigues SSgq_eh   = gq_eh.plus_left(Deltatau_bh_gq);

    check("gq jac left composition second 1 taylor   ", Sgq_eh.get_rotv()()(0), SSgq_eh.get_rotv()()(0), 1e-12);
    check("gq jac left composition second 2 taylor   ", Sgq_eh.get_rotv()()(1), SSgq_eh.get_rotv()()(1), 1e-12);
    check("gq jac left composition second 3 taylor   ", Sgq_eh.get_rotv()()(2), SSgq_eh.get_rotv()()(2), 1e-12);
    check("gq jac left composition second 4 taylor   ", Sgq_eh.get_T()(0), SSgq_eh.get_T()(0), 1e-12);
    check("gq jac left composition second 5 taylor   ", Sgq_eh.get_T()(1), SSgq_eh.get_T()(1), 1e-12);
    check("gq jac left composition second 6 taylor   ", Sgq_eh.get_T()(2), SSgq_eh.get_T()(2), 1e-12);

    speu_dcm SgR_bh        = gR_bh.plus_left(DeltatauB);
    speu_dcm gR_eh         = gR_eb * gR_bh;
    speu_dcm SgR_eh        = gR_eb * SgR_bh;
    Eigen::Matrix6d Jac_gR = gR_eb.jac_left_composition_wrt_second(gR_bh);
    trfv Deltatau_bh_gR    = ang::trfv(Jac_gR * DeltatauB());
    speu_dcm SSgR_eh       = gR_eh.plus_left(Deltatau_bh_gR);

    check("gR jac left composition second 1 taylor   ", SgR_eh.get_rotv()()(0), SSgR_eh.get_rotv()()(0), 1e-12);
    check("gR jac left composition second 2 taylor   ", SgR_eh.get_rotv()()(1), SSgR_eh.get_rotv()()(1), 1e-12);
    check("gR jac left composition second 3 taylor   ", SgR_eh.get_rotv()()(2), SSgR_eh.get_rotv()()(2), 1e-12);
    check("gR jac left composition second 4 taylor   ", SgR_eh.get_T()(0), SSgR_eh.get_T()(0), 1e-12);
    check("gR jac left composition second 5 taylor   ", SgR_eh.get_T()(1), SSgR_eh.get_T()(1), 1e-12);
    check("gR jac left composition second 6 taylor   ", SgR_eh.get_T()(2), SSgR_eh.get_T()(2), 1e-12);

    homogeneous SM_bh     = M_bh.plus_left(DeltatauB);
    homogeneous M_eh      = M_eb * M_bh;
    homogeneous SM_eh     = M_eb * SM_bh;
    Eigen::Matrix6d Jac_M = M_eb.jac_left_composition_wrt_second(M_bh);
    trfv Deltatau_bh_M    = ang::trfv(Jac_M * DeltatauB());
    homogeneous SSM_eh    = M_eh.plus_left(Deltatau_bh_M);

    check("M jac left composition second 1 taylor    ", SM_eh.get_rotv()()(0), SSM_eh.get_rotv()()(0), 1e-12);
    check("M jac left composition second 2 taylor    ", SM_eh.get_rotv()()(1), SSM_eh.get_rotv()()(1), 1e-12);
    check("M jac left composition second 3 taylor    ", SM_eh.get_rotv()()(2), SSM_eh.get_rotv()()(2), 1e-12);
    check("M jac left composition second 4 taylor    ", SM_eh.get_T()(0), SSM_eh.get_T()(0), 1e-12);
    check("M jac left composition second 5 taylor    ", SM_eh.get_T()(1), SSM_eh.get_T()(1), 1e-12);
    check("M jac left composition second 6 taylor    ", SM_eh.get_T()(2), SSM_eh.get_T()(2), 1e-12);

    dual SZ_bh            = Z_bh.plus_left(DeltatauB);
    dual Z_eh             = Z_eb * Z_bh;
    dual SZ_eh            = Z_eb * SZ_bh;
    Eigen::Matrix6d Jac_Z = Z_eb.jac_left_composition_wrt_second(Z_bh);
    trfv Deltatau_bh_Z    = ang::trfv(Jac_Z * DeltatauB());
    dual SSZ_eh           = Z_eh.plus_left(Deltatau_bh_Z);

    check("Z jac left composition second 1 taylor    ", SZ_eh.get_rotv()()(0), SSZ_eh.get_rotv()()(0), 1e-12);
    check("Z jac left composition second 2 taylor    ", SZ_eh.get_rotv()()(1), SSZ_eh.get_rotv()()(1), 1e-12);
    check("Z jac left composition second 3 taylor    ", SZ_eh.get_rotv()()(2), SSZ_eh.get_rotv()()(2), 1e-12);
    check("Z jac left composition second 4 taylor    ", SZ_eh.get_T()(0), SSZ_eh.get_T()(0), 1e-12);
    check("Z jac left composition second 5 taylor    ", SZ_eh.get_T()(1), SSZ_eh.get_T()(1), 1e-12);
    check("Z jac left composition second 6 taylor    ", SZ_eh.get_T()(2), SSZ_eh.get_T()(2), 1e-12);

    trfv Stau_bh            = tau_bh.plus_left(DeltatauB);
    trfv tau_eh             = tau_eb * tau_bh;
    trfv Stau_eh            = tau_eb * Stau_bh;
    Eigen::Matrix6d Jac_tau = tau_eb.jac_left_composition_wrt_second(tau_bh);
    trfv Deltatau_bh_tau    = ang::trfv(Jac_tau * DeltatauB());
    trfv SStau_eh           = tau_eh.plus_left(Deltatau_bh_tau);

    check("tau jac left composition second 1 taylor  ", Stau_eh.get_rotv()()(0), SStau_eh.get_rotv()()(0), 1e-12);
    check("tau jac left composition second 2 taylor  ", Stau_eh.get_rotv()()(1), SStau_eh.get_rotv()()(1), 1e-12);
    check("tau jac left composition second 3 taylor  ", Stau_eh.get_rotv()()(2), SStau_eh.get_rotv()()(2), 1e-12);
    check("tau jac left composition second 4 taylor  ", Stau_eh.get_T()(0), SStau_eh.get_T()(0), 1e-12);
    check("tau jac left composition second 5 taylor  ", Stau_eh.get_T()(1), SStau_eh.get_T()(1), 1e-12);
    check("tau jac left composition second 6 taylor  ", Stau_eh.get_T()(2), SStau_eh.get_T()(2), 1e-12);

} // closes test_jac_left_composition_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_forward_motion_wrt_motion() {

    // It verifies the behavior of the motion classes jac_left_forward_motion_wrt_motion method, validating the following:
    // (DeltatauE plus M) * p = M * p + J * DeltatauE

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(v);

    Eigen::Vector3d pB(0.33, -0.65, 0.17);

    Eigen::Vector3d pE_gq      = gq_eb * pB;
    speu_rodrigues Sgq_eb      = gq_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpE_gq     = Sgq_eb * pB;
    Eigen::Matrix36d Jac_gq    = gq_eb.jac_left_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_gq = Jac_gq * DeltatauE();
    Eigen::Vector3d SSpE_gq    = pE_gq + DeltapE_gq;

    check("gq jac left forward motion motion 1 taylor    ", SpE_gq(0), SSpE_gq(0), 1e-4);
    check("gq jac left forward motion motion 2 taylor    ", SpE_gq(1), SSpE_gq(1), 1e-4);
    check("gq jac left forward motion motion 3 taylor    ", SpE_gq(2), SSpE_gq(2), 1e-4);

    Eigen::Vector3d pE_gR      = gR_eb * pB;
    speu_dcm SgR_eb            = gR_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpE_gR     = SgR_eb * pB;
    Eigen::Matrix36d Jac_gR    = gR_eb.jac_left_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_gR = Jac_gR * DeltatauE();
    Eigen::Vector3d SSpE_gR    = pE_gR + DeltapE_gR;

    check("gR jac left forward motion motion 1 taylor    ", SpE_gR(0), SSpE_gR(0), 1e-4);
    check("gR jac left forward motion motion 2 taylor    ", SpE_gR(1), SSpE_gR(1), 1e-4);
    check("gR jac left forward motion motion 3 taylor    ", SpE_gR(2), SSpE_gR(2), 1e-4);

    Eigen::Vector3d pE_M      = M_eb * pB;
    homogeneous SM_eb         = M_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpE_M     = SM_eb * pB;
    Eigen::Matrix36d Jac_M    = M_eb.jac_left_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_M = Jac_M * DeltatauE();
    Eigen::Vector3d SSpE_M    = pE_M + DeltapE_M;

    check("M jac left forward motion motion 1 taylor     ", SpE_M(0), SSpE_M(0), 1e-4);
    check("M jac left forward motion motion 2 taylor     ", SpE_M(1), SSpE_M(1), 1e-4);
    check("M jac left forward motion motion 3 taylor     ", SpE_M(2), SSpE_M(2), 1e-4);

    Eigen::Vector3d pE_Z      = Z_eb * pB;
    dual SZ_eb                = Z_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpE_Z     = SZ_eb * pB;
    Eigen::Matrix36d Jac_Z    = Z_eb.jac_left_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_Z = Jac_Z * DeltatauE();
    Eigen::Vector3d SSpE_Z    = pE_Z + DeltapE_Z;

    check("Z jac left forward motion motion 1 taylor     ", SpE_Z(0), SSpE_Z(0), 1e-4);
    check("Z jac left forward motion motion 2 taylor     ", SpE_Z(1), SSpE_Z(1), 1e-4);
    check("Z jac left forward motion motion 3 taylor     ", SpE_Z(2), SSpE_Z(2), 1e-4);

    Eigen::Vector3d pE_tau      = tau_eb * pB;
    trfv Stau_eb                = tau_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpE_tau     = Stau_eb * pB;
    Eigen::Matrix36d Jac_tau    = tau_eb.jac_left_forward_motion_wrt_motion(pB);
    Eigen::Vector3d DeltapE_tau = Jac_tau * DeltatauE();
    Eigen::Vector3d SSpE_tau    = pE_tau + DeltapE_tau;

    check("tau jac left forward motion motion 1 taylor   ", SpE_tau(0), SSpE_tau(0), 1e-4);
    check("tau jac left forward motion motion 2 taylor   ", SpE_tau(1), SSpE_tau(1), 1e-4);
    check("tau jac left forward motion motion 3 taylor   ", SpE_tau(2), SSpE_tau(2), 1e-4);

} // closes test_jac_left_forward_motion_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_backward_motion_wrt_motion() {

    // It verifies the behavior of the motion classes jac_left_backward_motion_wrt_motion method, validating the following:
    // (DeltatauE plus M) / p = M / p + J * DeltatauE

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(v);

    Eigen::Vector3d pE(0.33, -0.65, 0.17);

    Eigen::Vector3d pB_gq      = gq_eb / pE;
    speu_rodrigues Sgq_eb      = gq_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpB_gq     = Sgq_eb / pE;
    Eigen::Matrix36d Jac_gq    = gq_eb.jac_left_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_gq = Jac_gq * DeltatauE();
    Eigen::Vector3d SSpB_gq    = pB_gq + DeltapB_gq;

    check("gq jac left backward motion motion 1 taylor    ", SpB_gq(0), SSpB_gq(0), 1e-4);
    check("gq jac left backward motion motion 2 taylor    ", SpB_gq(1), SSpB_gq(1), 1e-4);
    check("gq jac left backward motion motion 3 taylor    ", SpB_gq(2), SSpB_gq(2), 1e-4);

    Eigen::Vector3d pB_gR      = gR_eb / pE;
    speu_dcm SgR_eb            = gR_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpB_gR     = SgR_eb / pE;
    Eigen::Matrix36d Jac_gR    = gR_eb.jac_left_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_gR = Jac_gR * DeltatauE();
    Eigen::Vector3d SSpB_gR    = pB_gR + DeltapB_gR;

    check("gR jac left backward motion motion 1 taylor    ", SpB_gR(0), SSpB_gR(0), 1e-4);
    check("gR jac left backward motion motion 2 taylor    ", SpB_gR(1), SSpB_gR(1), 1e-4);
    check("gR jac left backward motion motion 3 taylor    ", SpB_gR(2), SSpB_gR(2), 1e-4);

    Eigen::Vector3d pB_M      = M_eb / pE;
    homogeneous SM_eb         = M_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpB_M     = SM_eb / pE;
    Eigen::Matrix36d Jac_M    = M_eb.jac_left_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_M = Jac_M * DeltatauE();
    Eigen::Vector3d SSpB_M    = pB_M + DeltapB_M;

    check("M jac left backward motion motion 1 taylor     ", SpB_M(0), SSpB_M(0), 1e-4);
    check("M jac left backward motion motion 2 taylor     ", SpB_M(1), SSpB_M(1), 1e-4);
    check("M jac left backward motion motion 3 taylor     ", SpB_M(2), SSpB_M(2), 1e-4);

    Eigen::Vector3d pB_Z      = Z_eb / pE;
    dual SZ_eb                = Z_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpB_Z     = SZ_eb / pE;
    Eigen::Matrix36d Jac_Z    = Z_eb.jac_left_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_Z = Jac_Z * DeltatauE();
    Eigen::Vector3d SSpB_Z    = pB_Z + DeltapB_Z;

    check("Z jac left backward motion motion 1 taylor     ", SpB_Z(0), SSpB_Z(0), 1e-4);
    check("Z jac left backward motion motion 2 taylor     ", SpB_Z(1), SSpB_Z(1), 1e-4);
    check("Z jac left backward motion motion 3 taylor     ", SpB_Z(2), SSpB_Z(2), 1e-4);

    Eigen::Vector3d pB_tau      = tau_eb / pE;
    trfv Stau_eb                = tau_eb.plus_left(DeltatauE);
    Eigen::Vector3d SpB_tau     = Stau_eb / pE;
    Eigen::Matrix36d Jac_tau    = tau_eb.jac_left_backward_motion_wrt_motion(pE);
    Eigen::Vector3d DeltapB_tau = Jac_tau * DeltatauE();
    Eigen::Vector3d SSpB_tau    = pB_tau + DeltapB_tau;

    check("tau jac left backward motion motion 1 taylor   ", SpB_tau(0), SSpB_tau(0), 1e-4);
    check("tau jac left backward motion motion 2 taylor   ", SpB_tau(1), SSpB_tau(1), 1e-4);
    check("tau jac left backward motion motion 3 taylor   ", SpB_tau(2), SSpB_tau(2), 1e-4);

} // closes test_jac_left_backward_motion_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_log() {

    // It checks the proper behavior of the motion left logarithmic jacobian method,
    // validating the following equations:
    // Log(Deltatau plus M) ~= Log(M) + (J * Deltatau)

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d v; v << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(v);

    speu_rodrigues Sgq_eb       = gq_eb.plus_left(DeltatauE);
    trfv Stau_eb_gq             = Sgq_eb.log_map_trfv();
    trfv tau_eb_gq              = gq_eb.log_map_trfv();
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_left_log();
    trfv Deltatau_eb_gq         = trfv(Jac_gq * DeltatauE());
    trfv SStau_eb_gq(tau_eb_gq() + Deltatau_eb_gq());

    check("gq jac left log 1 taylor    ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-5);
    check("gq jac left log 2 taylor    ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-5);
    check("gq jac left log 3 taylor    ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-5);
    check("gq jac left log 4 taylor    ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-6);
    check("gq jac left log 5 taylor    ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-6);
    check("gq jac left log 6 taylor    ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-6);

    speu_dcm SgR_eb             = gR_eb.plus_left(DeltatauE);
    trfv Stau_eb_gR             = SgR_eb.log_map_trfv();
    trfv tau_eb_gR              = gR_eb.log_map_trfv();
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_left_log();
    trfv Deltatau_eb_gR         = trfv(Jac_gR * DeltatauE());
    trfv SStau_eb_gR(tau_eb_gR() + Deltatau_eb_gR());

    check("gR jac left log 1 taylor    ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-5);
    check("gR jac left log 2 taylor    ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-5);
    check("gR jac left log 3 taylor    ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-5);
    check("gR jac left log 4 taylor    ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-6);
    check("gR jac left log 5 taylor    ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-6);
    check("gR jac left log 6 taylor    ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-6);

    homogeneous SM_eb          = M_eb.plus_left(DeltatauE);
    trfv Stau_eb_M             = SM_eb.log_map_trfv();
    trfv tau_eb_M              = M_eb.log_map_trfv();
    Eigen::Matrix6d Jac_M      = M_eb.jac_left_log();
    trfv Deltatau_eb_M         = trfv(Jac_M * DeltatauE());
    trfv SStau_eb_M(tau_eb_M() + Deltatau_eb_M());

    check("M jac left log 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-5);
    check("M jac left log 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-5);
    check("M jac left log 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-5);
    check("M jac left log 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-6);
    check("M jac left log 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-6);
    check("M jac left log 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-6);

    dual SZ_eb                 = Z_eb.plus_left(DeltatauE);
    trfv Stau_eb_Z             = SZ_eb.log_map_trfv();
    trfv tau_eb_Z              = Z_eb.log_map_trfv();
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_left_log();
    trfv Deltatau_eb_Z         = trfv(Jac_Z * DeltatauE());
    trfv SStau_eb_Z(tau_eb_Z() + Deltatau_eb_Z());

    check("Z jac left log 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-5);
    check("Z jac left log 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-5);
    check("Z jac left log 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-5);
    check("Z jac left log 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-6);
    check("Z jac left log 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-6);
    check("Z jac left log 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-6);

    trfv Stau_eb                 = tau_eb.plus_left(DeltatauE);
    trfv Stau_eb_tau             = Stau_eb;
    trfv tau_eb_tau              = tau_eb;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_left_log();
    trfv Deltatau_eb_tau         = trfv(Jac_tau * DeltatauE());
    trfv SStau_eb_tau(tau_eb_tau() + Deltatau_eb_tau());

    check("tau jac left log 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-5);
    check("tau jac left log 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-5);
    check("tau jac left log 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-5);
    check("tau jac left log 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-6);
    check("tau jac left log 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-6);
    check("tau jac left log 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-6);

} // closes test_jacobian_left_log

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_plus_wrt_first() {

    // It verifies the behavior of the motion classes jac_left_plus_wrt_first method, validating the following:
    // (tau1 + Deltatau1E) plus M2 = J * Deltatau1E plus (tau1 plus M2)

    rotv r_nb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_nbn(7.5, -12.8, -3.2);

    trfv tau_nb(r_nb, T_nbn);
    speu_dcm gR_nb(tau_nb);
    speu_rodrigues gq_nb(tau_nb);
    homogeneous M_nb(tau_nb);
    dual Z_nb(tau_nb);

    Eigen::Vector6d v; v << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauN(v);

    rotv r_en(0.4, -1.2, 0.72);
    Eigen::Vector3d T_ene(-2.8, -3.9, 4.2);
    trfv tau_en(r_en, T_ene);

    trfv tautau_en(tau_en() + DeltatauN());

    speu_rodrigues gq_eb   = gq_nb.plus_left(tau_en);
    speu_rodrigues Sgq_eb  = gq_nb.plus_left(tautau_en);
    Eigen::Matrix6d Jac_gq = gq_nb.jac_left_plus_wrt_first(tau_en);
    trfv Deltatau_en_gq    = ang::trfv(Jac_gq * DeltatauN());
    speu_rodrigues SSgq_eb = gq_eb.plus_left(Deltatau_en_gq);

    check("gq jac left plus first 1 taylor   ", Sgq_eb.get_rotv()()(0), SSgq_eb.get_rotv()()(0), 1e-6);
    check("gq jac left plus first 2 taylor   ", Sgq_eb.get_rotv()()(1), SSgq_eb.get_rotv()()(1), 1e-6);
    check("gq jac left plus first 3 taylor   ", Sgq_eb.get_rotv()()(2), SSgq_eb.get_rotv()()(2), 1e-6);
    check("gq jac left plus first 4 taylor   ", Sgq_eb.get_T()(0), SSgq_eb.get_T()(0), 1e-5);
    check("gq jac left plus first 5 taylor   ", Sgq_eb.get_T()(1), SSgq_eb.get_T()(1), 1e-5);
    check("gq jac left plus first 6 taylor   ", Sgq_eb.get_T()(2), SSgq_eb.get_T()(2), 1e-5);

    speu_dcm gR_eb         = gR_nb.plus_left(tau_en);
    speu_dcm SgR_eb        = gR_nb.plus_left(tautau_en);
    Eigen::Matrix6d Jac_gR = gR_nb.jac_left_plus_wrt_first(tau_en);
    trfv Deltatau_en_gR    = ang::trfv(Jac_gR * DeltatauN());
    speu_dcm SSgR_eb       = gR_eb.plus_left(Deltatau_en_gR);

    check("gR jac left plus first 1 taylor   ", SgR_eb.get_rotv()()(0), SSgR_eb.get_rotv()()(0), 1e-6);
    check("gR jac left plus first 2 taylor   ", SgR_eb.get_rotv()()(1), SSgR_eb.get_rotv()()(1), 1e-6);
    check("gR jac left plus first 3 taylor   ", SgR_eb.get_rotv()()(2), SSgR_eb.get_rotv()()(2), 1e-6);
    check("gR jac left plus first 4 taylor   ", SgR_eb.get_T()(0), SSgR_eb.get_T()(0), 1e-5);
    check("gR jac left plus first 5 taylor   ", SgR_eb.get_T()(1), SSgR_eb.get_T()(1), 1e-5);
    check("gR jac left plus first 6 taylor   ", SgR_eb.get_T()(2), SSgR_eb.get_T()(2), 1e-5);

    homogeneous M_eb      = M_nb.plus_left(tau_en);
    homogeneous SM_eb     = M_nb.plus_left(tautau_en);
    Eigen::Matrix6d Jac_M = M_nb.jac_left_plus_wrt_first(tau_en);
    trfv Deltatau_en_M    = ang::trfv(Jac_M * DeltatauN());
    homogeneous SSM_eb    = M_eb.plus_left(Deltatau_en_M);

    check("M jac left plus first 1 taylor    ", SM_eb.get_rotv()()(0), SSM_eb.get_rotv()()(0), 1e-6);
    check("M jac left plus first 2 taylor    ", SM_eb.get_rotv()()(1), SSM_eb.get_rotv()()(1), 1e-6);
    check("M jac left plus first 3 taylor    ", SM_eb.get_rotv()()(2), SSM_eb.get_rotv()()(2), 1e-6);
    check("M jac left plus first 4 taylor    ", SM_eb.get_T()(0), SSM_eb.get_T()(0), 1e-5);
    check("M jac left plus first 5 taylor    ", SM_eb.get_T()(1), SSM_eb.get_T()(1), 1e-5);
    check("M jac left plus first 6 taylor    ", SM_eb.get_T()(2), SSM_eb.get_T()(2), 1e-5);

    dual Z_eb             = Z_nb.plus_left(tau_en);
    dual SZ_eb            = Z_nb.plus_left(tautau_en);
    Eigen::Matrix6d Jac_Z = Z_nb.jac_left_plus_wrt_first(tau_en);
    trfv Deltatau_en_Z    = ang::trfv(Jac_Z * DeltatauN());
    dual SSZ_eb           = Z_eb.plus_left(Deltatau_en_Z);

    check("Z jac left plus first 1 taylor    ", SZ_eb.get_rotv()()(0), SSZ_eb.get_rotv()()(0), 1e-6);
    check("Z jac left plus first 2 taylor    ", SZ_eb.get_rotv()()(1), SSZ_eb.get_rotv()()(1), 1e-6);
    check("Z jac left plus first 3 taylor    ", SZ_eb.get_rotv()()(2), SSZ_eb.get_rotv()()(2), 1e-6);
    check("Z jac left plus first 4 taylor    ", SZ_eb.get_T()(0), SSZ_eb.get_T()(0), 1e-5);
    check("Z jac left plus first 5 taylor    ", SZ_eb.get_T()(1), SSZ_eb.get_T()(1), 1e-5);
    check("Z jac left plus first 6 taylor    ", SZ_eb.get_T()(2), SSZ_eb.get_T()(2), 1e-5);

    trfv tau_eb             = tau_nb.plus_left(tau_en);
    trfv Stau_eb            = tau_nb.plus_left(tautau_en);
    Eigen::Matrix6d Jac_tau = tau_nb.jac_left_plus_wrt_first(tau_en);
    trfv Deltatau_en_tau    = ang::trfv(Jac_tau * DeltatauN());
    trfv SStau_eb           = tau_eb.plus_left(Deltatau_en_tau);

    check("tau jac left plus first 1 taylor  ", Stau_eb.get_rotv()()(0), SStau_eb.get_rotv()()(0), 1e-6);
    check("tau jac left plus first 2 taylor  ", Stau_eb.get_rotv()()(1), SStau_eb.get_rotv()()(1), 1e-6);
    check("tau jac left plus first 3 taylor  ", Stau_eb.get_rotv()()(2), SStau_eb.get_rotv()()(2), 1e-6);
    check("tau jac left plus first 4 taylor  ", Stau_eb.get_T()(0), SStau_eb.get_T()(0), 1e-5);
    check("tau jac left plus first 5 taylor  ", Stau_eb.get_T()(1), SStau_eb.get_T()(1), 1e-5);
    check("tau jac left plus first 6 taylor  ", Stau_eb.get_T()(2), SStau_eb.get_T()(2), 1e-5);

} // closes test_jac_left_plus_wrt_first

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_plus_wrt_second() {

    // It verifies the behavior of the rotation classes jac_left_plus_wrt_second method, validating the following:
    // tau1 plus (Deltatau2E plus M2) = J * Deltatau2E plus (tau1 plus M2)

    rotv r_nb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_nbn(7.5, -12.8, -3.2);

    trfv tau_nb(r_nb, T_nbn);
    speu_dcm gR_nb(tau_nb);
    speu_rodrigues gq_nb(tau_nb);
    homogeneous M_nb(tau_nb);
    dual Z_nb(tau_nb);

    Eigen::Vector6d v; v << 0.04, -0.05, 0.07, 1e-3, -5e-4, 2e-3;
    trfv DeltatauN(v);

    rotv r_en(0.4, -1.2, 0.72);
    Eigen::Vector3d T_ene(-2.8, -3.9, 4.2);
    trfv tau_en(r_en, T_ene);

    speu_rodrigues gq_eb   = gq_nb.plus_left(tau_en);
    speu_rodrigues Sgq_nb  = gq_nb.plus_left(DeltatauN);
    speu_rodrigues Sgq_eb  = Sgq_nb.plus_left(tau_en);
    Eigen::Matrix6d Jac_gq = gq_nb.jac_left_plus_wrt_second(tau_en);
    trfv Deltatau_eb_gq    = ang::trfv(Jac_gq * DeltatauN());
    speu_rodrigues SSgq_eb = gq_eb.plus_left(Deltatau_eb_gq);

    check("gq jac left plus second 1 taylor   ", Sgq_eb.get_rotv()()(0), SSgq_eb.get_rotv()()(0), 1e-12);
    check("gq jac left plus second 2 taylor   ", Sgq_eb.get_rotv()()(1), SSgq_eb.get_rotv()()(1), 1e-12);
    check("gq jac left plus second 3 taylor   ", Sgq_eb.get_rotv()()(2), SSgq_eb.get_rotv()()(2), 1e-12);
    check("gq jac left plus second 4 taylor   ", Sgq_eb.get_T()(0), SSgq_eb.get_T()(0), 1e-12);
    check("gq jac left plus second 5 taylor   ", Sgq_eb.get_T()(1), SSgq_eb.get_T()(1), 1e-12);
    check("gq jac left plus second 6 taylor   ", Sgq_eb.get_T()(2), SSgq_eb.get_T()(2), 1e-12);

    speu_dcm gR_eb         = gR_nb.plus_left(tau_en);
    speu_dcm SgR_nb        = gR_nb.plus_left(DeltatauN);
    speu_dcm SgR_eb        = SgR_nb.plus_left(tau_en);
    Eigen::Matrix6d Jac_gR = gR_nb.jac_left_plus_wrt_second(tau_en);
    trfv Deltatau_eb_gR    = ang::trfv(Jac_gR * DeltatauN());
    speu_dcm SSgR_eb       = gR_eb.plus_left(Deltatau_eb_gR);

    check("gR jac left plus second 1 taylor   ", SgR_eb.get_rotv()()(0), SSgR_eb.get_rotv()()(0), 1e-12);
    check("gR jac left plus second 2 taylor   ", SgR_eb.get_rotv()()(1), SSgR_eb.get_rotv()()(1), 1e-12);
    check("gR jac left plus second 3 taylor   ", SgR_eb.get_rotv()()(2), SSgR_eb.get_rotv()()(2), 1e-12);
    check("gR jac left plus second 4 taylor   ", SgR_eb.get_T()(0), SSgR_eb.get_T()(0), 1e-12);
    check("gR jac left plus second 5 taylor   ", SgR_eb.get_T()(1), SSgR_eb.get_T()(1), 1e-12);
    check("gR jac left plus second 6 taylor   ", SgR_eb.get_T()(2), SSgR_eb.get_T()(2), 1e-12);

    homogeneous M_eb      = M_nb.plus_left(tau_en);
    homogeneous SM_nb     = M_nb.plus_left(DeltatauN);
    homogeneous SM_eb     = SM_nb.plus_left(tau_en);
    Eigen::Matrix6d Jac_M = M_nb.jac_left_plus_wrt_second(tau_en);
    trfv Deltatau_eb_M    = ang::trfv(Jac_M * DeltatauN());
    homogeneous SSM_eb    = M_eb.plus_left(Deltatau_eb_M);

    check("M jac left plus second 1 taylor    ", SM_eb.get_rotv()()(0), SSM_eb.get_rotv()()(0), 1e-12);
    check("M jac left plus second 2 taylor    ", SM_eb.get_rotv()()(1), SSM_eb.get_rotv()()(1), 1e-12);
    check("M jac left plus second 3 taylor    ", SM_eb.get_rotv()()(2), SSM_eb.get_rotv()()(2), 1e-12);
    check("M jac left plus second 4 taylor    ", SM_eb.get_T()(0), SSM_eb.get_T()(0), 1e-12);
    check("M jac left plus second 5 taylor    ", SM_eb.get_T()(1), SSM_eb.get_T()(1), 1e-12);
    check("M jac left plus second 6 taylor    ", SM_eb.get_T()(2), SSM_eb.get_T()(2), 1e-12);

    dual Z_eb             = Z_nb.plus_left(tau_en);
    dual SZ_nb            = Z_nb.plus_left(DeltatauN);
    dual SZ_eb            = SZ_nb.plus_left(tau_en);
    Eigen::Matrix6d Jac_Z = Z_nb.jac_left_plus_wrt_second(tau_en);
    trfv Deltatau_eb_Z    = ang::trfv(Jac_Z * DeltatauN());
    dual SSZ_eb           = Z_eb.plus_left(Deltatau_eb_Z);

    check("Z jac left plus second 1 taylor    ", SZ_eb.get_rotv()()(0), SSZ_eb.get_rotv()()(0), 1e-12);
    check("Z jac left plus second 2 taylor    ", SZ_eb.get_rotv()()(1), SSZ_eb.get_rotv()()(1), 1e-12);
    check("Z jac left plus second 3 taylor    ", SZ_eb.get_rotv()()(2), SSZ_eb.get_rotv()()(2), 1e-12);
    check("Z jac left plus second 4 taylor    ", SZ_eb.get_T()(0), SSZ_eb.get_T()(0), 1e-12);
    check("Z jac left plus second 5 taylor    ", SZ_eb.get_T()(1), SSZ_eb.get_T()(1), 1e-12);
    check("Z jac left plus second 6 taylor    ", SZ_eb.get_T()(2), SSZ_eb.get_T()(2), 1e-12);

    trfv tau_eb             = tau_nb.plus_left(tau_en);
    trfv Stau_nb            = tau_nb.plus_left(DeltatauN);
    trfv Stau_eb            = Stau_nb.plus_left(tau_en);
    Eigen::Matrix6d Jac_tau = tau_nb.jac_left_plus_wrt_second(tau_en);
    trfv Deltatau_eb_tau    = ang::trfv(Jac_tau * DeltatauN());
    trfv SStau_eb           = tau_eb.plus_left(Deltatau_eb_tau);

    check("tau jac left plus second 1 taylor  ", Stau_eb.get_rotv()()(0), SStau_eb.get_rotv()()(0), 1e-12);
    check("tau jac left plus second 2 taylor  ", Stau_eb.get_rotv()()(1), SStau_eb.get_rotv()()(1), 1e-12);
    check("tau jac left plus second 3 taylor  ", Stau_eb.get_rotv()()(2), SStau_eb.get_rotv()()(2), 1e-12);
    check("tau jac left plus second 4 taylor  ", Stau_eb.get_T()(0), SStau_eb.get_T()(0), 1e-12);
    check("tau jac left plus second 5 taylor  ", Stau_eb.get_T()(1), SStau_eb.get_T()(1), 1e-12);
    check("tau jac left plus second 6 taylor  ", Stau_eb.get_T()(2), SStau_eb.get_T()(2), 1e-12);

} // closes test_jac_left_plus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_minus_wrt_first() {

    // It verifies the behavior of the motion classes jac_left_minus_wrt_first method, validating the following:
    // (Deltatau2E plus M2) minus M1 = (M2 minus M1) + J * Deltatau2E */

    rotv r1_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T1_ebe(7.5, -12.8, -3.2);

    trfv tau1_eb(r1_eb, T1_ebe);
    speu_dcm gR1_eb(tau1_eb);
    speu_rodrigues gq1_eb(tau1_eb);
    homogeneous M1_eb(tau1_eb);
    dual Z1_eb(tau1_eb);

    Eigen::Vector6d vA; vA << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauEA(vA);

    trfv tau2_eb          = tau1_eb.plus_left(DeltatauEA);
    speu_dcm gR2_eb       = gR1_eb.plus_left(DeltatauEA);
    speu_rodrigues gq2_eb = gq1_eb.plus_left(DeltatauEA);
    homogeneous M2_eb     = M1_eb.plus_left(DeltatauEA);
    dual Z2_eb            = Z1_eb.plus_left(DeltatauEA);

    Eigen::Vector6d v; v << -0.002, -0.003, -0.002, -1.5e-3, -9e-4, 2.8e-3;
    trfv DeltatauE(v);

    trfv DeltatauE_gq        = gq2_eb.minus_left_trfv(gq1_eb);
    speu_rodrigues Sgq_eb    = gq2_eb.plus_left(DeltatauE);
    trfv Stau_eb_gq          = Sgq_eb.minus_left_trfv(gq1_eb);
    Eigen::Matrix6d Jac_gq   = gq2_eb.jac_left_minus_wrt_first(gq1_eb);
    trfv Deltatau_eb_gq      = ang::trfv(Jac_gq * DeltatauE());
    trfv SStau_eb_gq(DeltatauE_gq() + Deltatau_eb_gq());

    check("gq jac left minus first 1 taylor   ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-7);
    check("gq jac left minus first 2 taylor   ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-7);
    check("gq jac left minus first 3 taylor   ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-7);
    check("gq jac left minus first 4 taylor   ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-7);
    check("gq jac left minus first 5 taylor   ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-7);
    check("gq jac left minus first 6 taylor   ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-7);

    trfv DeltatauE_gR        = gR2_eb.minus_left_trfv(gR1_eb);
    speu_dcm SgR_eb          = gR2_eb.plus_left(DeltatauE);
    trfv Stau_eb_gR          = SgR_eb.minus_left_trfv(gR1_eb);
    Eigen::Matrix6d Jac_gR = gR2_eb.jac_left_minus_wrt_first(gR1_eb);
    trfv Deltatau_eb_gR      = ang::trfv(Jac_gR * DeltatauE());
    trfv SStau_eb_gR(DeltatauE_gR() + Deltatau_eb_gR());

    check("gR jac left minus first 1 taylor   ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-7);
    check("gR jac left minus first 2 taylor   ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-7);
    check("gR jac left minus first 3 taylor   ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-7);
    check("gR jac left minus first 4 taylor   ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-7);
    check("gR jac left minus first 5 taylor   ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-7);
    check("gR jac left minus first 6 taylor   ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-7);

    trfv DeltatauE_M      = M2_eb.minus_left_trfv(M1_eb);
    homogeneous SM_eb     = M2_eb.plus_left(DeltatauE);
    trfv Stau_eb_M        = SM_eb.minus_left_trfv(M1_eb);
    Eigen::Matrix6d Jac_M = M2_eb.jac_left_minus_wrt_first(M1_eb);
    trfv Deltatau_eb_M    = ang::trfv(Jac_M * DeltatauE());
    trfv SStau_eb_M(DeltatauE_M() + Deltatau_eb_M());

    check("M jac left minus first 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-7);
    check("M jac left minus first 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-7);
    check("M jac left minus first 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-7);
    check("M jac left minus first 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-7);
    check("M jac left minus first 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-7);
    check("M jac left minus first 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-7);

    trfv DeltatauE_Z      = Z2_eb.minus_left_trfv(Z1_eb);
    dual SZ_eb            = Z2_eb.plus_left(DeltatauE);
    trfv Stau_eb_Z        = SZ_eb.minus_left_trfv(Z1_eb);
    Eigen::Matrix6d Jac_Z = Z2_eb.jac_left_minus_wrt_first(Z1_eb);
    trfv Deltatau_eb_Z    = ang::trfv(Jac_Z * DeltatauE());
    trfv SStau_eb_Z(DeltatauE_Z() + Deltatau_eb_Z());

    check("Z jac left minus first 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-7);
    check("Z jac left minus first 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-7);
    check("Z jac left minus first 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-7);
    check("Z jac left minus first 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-7);
    check("Z jac left minus first 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-7);
    check("Z jac left minus first 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-7);

    trfv DeltatauE_tau      = tau2_eb.minus_left_trfv(tau1_eb);
    trfv Stau_eb            = tau2_eb.plus_left(DeltatauE);
    trfv Stau_eb_tau        = Stau_eb.minus_left_trfv(tau1_eb);
    Eigen::Matrix6d Jac_tau = tau2_eb.jac_left_minus_wrt_first(tau1_eb);
    trfv Deltatau_eb_tau    = ang::trfv(Jac_tau * DeltatauE());
    trfv SStau_eb_tau(DeltatauE_tau() + Deltatau_eb_tau());

    check("tau jac left minus first 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-7);
    check("tau jac left minus first 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-7);
    check("tau jac left minus first 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-7);
    check("tau jac left minus first 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-7);
    check("tau jac left minus first 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-7);
    check("tau jac left minus first 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-7);

} // closes test_jac_left_minus_wrt_first.

// ///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_minus_wrt_second() {

    // It verifies the behavior of the motion classes jac_left_minus_wrt_second method, validating the following:
    // M2 minus (Deltatau1E plus M1) = (M2 minus M1) + J * Deltatau1E *

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d vA; vA << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauEA(vA);

    trfv Ptau_eb          = tau_eb.plus_left(DeltatauEA);
    speu_dcm PgR_eb       = gR_eb.plus_left(DeltatauEA);
    speu_rodrigues Pgq_eb = gq_eb.plus_left(DeltatauEA);
    homogeneous PM_eb     = M_eb.plus_left(DeltatauEA);
    dual PZ_eb            = Z_eb.plus_left(DeltatauEA);

    Eigen::Vector6d v; v << 0.002, -0.003, -0.002, 3e-4, 7e-4, 8e-4;
    trfv DeltatauE(v);

    trfv DeltatauE_gq      = Pgq_eb.minus_left_trfv(gq_eb);
    speu_rodrigues Sgq_eb  = gq_eb.plus_left(DeltatauE);
    trfv Stau_eb_gq        = Pgq_eb.minus_left_trfv(Sgq_eb);
    Eigen::Matrix6d Jac_gq = Pgq_eb.jac_left_minus_wrt_second(gq_eb);
    trfv Deltatau_eb_gq    = ang::trfv(Jac_gq * DeltatauE());
    trfv SStau_eb_gq(DeltatauE_gq() + Deltatau_eb_gq());

    check("gq jac left minus second 1 taylor   ", Stau_eb_gq()(0), SStau_eb_gq()(0), 1e-7);
    check("gq jac left minus second 2 taylor   ", Stau_eb_gq()(1), SStau_eb_gq()(1), 1e-7);
    check("gq jac left minus second 3 taylor   ", Stau_eb_gq()(2), SStau_eb_gq()(2), 1e-7);
    check("gq jac left minus second 4 taylor   ", Stau_eb_gq()(3), SStau_eb_gq()(3), 1e-7);
    check("gq jac left minus second 5 taylor   ", Stau_eb_gq()(4), SStau_eb_gq()(4), 1e-7);
    check("gq jac left minus second 6 taylor   ", Stau_eb_gq()(5), SStau_eb_gq()(5), 1e-7);

    trfv DeltatauE_gR      = PgR_eb.minus_left_trfv(gR_eb);
    speu_dcm SgR_eb        = gR_eb.plus_left(DeltatauE);
    trfv Stau_eb_gR        = PgR_eb.minus_left_trfv(SgR_eb);
    Eigen::Matrix6d Jac_gR = PgR_eb.jac_left_minus_wrt_second(gR_eb);
    trfv Deltatau_eb_gR    = ang::trfv(Jac_gR * DeltatauE());
    trfv SStau_eb_gR(DeltatauE_gR() + Deltatau_eb_gR());

    check("gR jac left minus second 1 taylor   ", Stau_eb_gR()(0), SStau_eb_gR()(0), 1e-7);
    check("gR jac left minus second 2 taylor   ", Stau_eb_gR()(1), SStau_eb_gR()(1), 1e-7);
    check("gR jac left minus second 3 taylor   ", Stau_eb_gR()(2), SStau_eb_gR()(2), 1e-7);
    check("gR jac left minus second 4 taylor   ", Stau_eb_gR()(3), SStau_eb_gR()(3), 1e-7);
    check("gR jac left minus second 5 taylor   ", Stau_eb_gR()(4), SStau_eb_gR()(4), 1e-7);
    check("gR jac left minus second 6 taylor   ", Stau_eb_gR()(5), SStau_eb_gR()(5), 1e-7);

    trfv DeltatauE_M      = PM_eb.minus_left_trfv(M_eb);
    homogeneous SM_eb     = M_eb.plus_left(DeltatauE);
    trfv Stau_eb_M        = PM_eb.minus_left_trfv(SM_eb);
    Eigen::Matrix6d Jac_M = PM_eb.jac_left_minus_wrt_second(M_eb);
    trfv Deltatau_eb_M    = ang::trfv(Jac_M * DeltatauE());
    trfv SStau_eb_M(DeltatauE_M() + Deltatau_eb_M());

    check("M jac left minus second 1 taylor    ", Stau_eb_M()(0), SStau_eb_M()(0), 1e-7);
    check("M jac left minus second 2 taylor    ", Stau_eb_M()(1), SStau_eb_M()(1), 1e-7);
    check("M jac left minus second 3 taylor    ", Stau_eb_M()(2), SStau_eb_M()(2), 1e-7);
    check("M jac left minus second 4 taylor    ", Stau_eb_M()(3), SStau_eb_M()(3), 1e-7);
    check("M jac left minus second 5 taylor    ", Stau_eb_M()(4), SStau_eb_M()(4), 1e-7);
    check("M jac left minus second 6 taylor    ", Stau_eb_M()(5), SStau_eb_M()(5), 1e-7);

    trfv DeltatauE_Z      = PZ_eb.minus_left_trfv(Z_eb);
    dual SZ_eb            = Z_eb.plus_left(DeltatauE);
    trfv Stau_eb_Z        = PZ_eb.minus_left_trfv(SZ_eb);
    Eigen::Matrix6d Jac_Z = PZ_eb.jac_left_minus_wrt_second(Z_eb);
    trfv Deltatau_eb_Z    = ang::trfv(Jac_Z * DeltatauE());
    trfv SStau_eb_Z(DeltatauE_Z() + Deltatau_eb_Z());

    check("Z jac left minus second 1 taylor    ", Stau_eb_Z()(0), SStau_eb_Z()(0), 1e-7);
    check("Z jac left minus second 2 taylor    ", Stau_eb_Z()(1), SStau_eb_Z()(1), 1e-7);
    check("Z jac left minus second 3 taylor    ", Stau_eb_Z()(2), SStau_eb_Z()(2), 1e-7);
    check("Z jac left minus second 4 taylor    ", Stau_eb_Z()(3), SStau_eb_Z()(3), 1e-7);
    check("Z jac left minus second 5 taylor    ", Stau_eb_Z()(4), SStau_eb_Z()(4), 1e-7);
    check("Z jac left minus second 6 taylor    ", Stau_eb_Z()(5), SStau_eb_Z()(5), 1e-7);

    trfv DeltatauE_tau      = Ptau_eb.minus_left_trfv(tau_eb);
    trfv Stau_eb            = tau_eb.plus_left(DeltatauE);
    trfv Stau_eb_tau        = Ptau_eb.minus_left_trfv(Stau_eb);
    Eigen::Matrix6d Jac_tau = Ptau_eb.jac_left_minus_wrt_second(tau_eb);
    trfv Deltatau_eb_tau    = ang::trfv(Jac_tau * DeltatauE());
    trfv SStau_eb_tau(DeltatauE_tau() + Deltatau_eb_tau());

    check("tau jac left minus second 1 taylor  ", Stau_eb_tau()(0), SStau_eb_tau()(0), 1e-7);
    check("tau jac left minus second 2 taylor  ", Stau_eb_tau()(1), SStau_eb_tau()(1), 1e-7);
    check("tau jac left minus second 3 taylor  ", Stau_eb_tau()(2), SStau_eb_tau()(2), 1e-7);
    check("tau jac left minus second 4 taylor  ", Stau_eb_tau()(3), SStau_eb_tau()(3), 1e-7);
    check("tau jac left minus second 5 taylor  ", Stau_eb_tau()(4), SStau_eb_tau()(4), 1e-7);
    check("tau jac left minus second 6 taylor  ", Stau_eb_tau()(5), SStau_eb_tau()(5), 1e-7);

} // closes test_jac_left_minus_wrt_second

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_forward_adjoint_wrt_motion() {

    // It verifies the behavior of the motion classes jac_left_forward_adjoint_wrt_motion method, validating the following:
    // Ad(DeltatauE plus M) | xi = AdM | xi + J * DeltatauE

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d tau; tau << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(tau);

    se3_tangent xiB(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    se3_tangent xiE_gq          = gq_eb | xiB;
    speu_rodrigues Sgq_eb       = gq_eb.plus_left(DeltatauE);
    se3_tangent SxiE_gq         = Sgq_eb | xiB;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_left_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_gq = Jac_gq * DeltatauE();
    se3_tangent SSxiE_gq        = xiE_gq + DeltaxiE_gq;

    check("gq jac left forward adjoint motion 1 taylor    ", SxiE_gq()(0), SSxiE_gq()(0), 1e-4);
    check("gq jac left forward adjoint motion 2 taylor    ", SxiE_gq()(1), SSxiE_gq()(1), 1e-4);
    check("gq jac left forward adjoint motion 3 taylor    ", SxiE_gq()(2), SSxiE_gq()(2), 1e-4);
    check("gq jac left forward adjoint motion 4 taylor    ", SxiE_gq()(3), SSxiE_gq()(3), 1e-4);
    check("gq jac left forward adjoint motion 5 taylor    ", SxiE_gq()(4), SSxiE_gq()(4), 1e-4);
    check("gq jac left forward adjoint motion 6 taylor    ", SxiE_gq()(5), SSxiE_gq()(5), 1e-4);

    se3_tangent xiE_gR          = gR_eb | xiB;
    speu_dcm SgR_eb             = gR_eb.plus_left(DeltatauE);
    se3_tangent SxiE_gR         = SgR_eb | xiB;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_left_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_gR = Jac_gR * DeltatauE();
    se3_tangent SSxiE_gR        = xiE_gR + DeltaxiE_gR;

    check("gR jac left forward adjoint motion 1 taylor    ", SxiE_gR()(0), SSxiE_gR()(0), 1e-4);
    check("gR jac left forward adjoint motion 2 taylor    ", SxiE_gR()(1), SSxiE_gR()(1), 1e-4);
    check("gR jac left forward adjoint motion 3 taylor    ", SxiE_gR()(2), SSxiE_gR()(2), 1e-4);
    check("gR jac left forward adjoint motion 4 taylor    ", SxiE_gR()(3), SSxiE_gR()(3), 1e-4);
    check("gR jac left forward adjoint motion 5 taylor    ", SxiE_gR()(4), SSxiE_gR()(4), 1e-4);
    check("gR jac left forward adjoint motion 6 taylor    ", SxiE_gR()(5), SSxiE_gR()(5), 1e-4);

    se3_tangent xiE_M          = M_eb | xiB;
    homogeneous SM_eb          = M_eb.plus_left(DeltatauE);
    se3_tangent SxiE_M         = SM_eb | xiB;
    Eigen::Matrix6d Jac_M      = M_eb.jac_left_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_M = Jac_M * DeltatauE();
    se3_tangent SSxiE_M        = xiE_M + DeltaxiE_M;

    check("M jac left forward adjoint motion 1 taylor     ", SxiE_M()(0), SSxiE_M()(0), 1e-4);
    check("M jac left forward adjoint motion 2 taylor     ", SxiE_M()(1), SSxiE_M()(1), 1e-4);
    check("M jac left forward adjoint motion 3 taylor     ", SxiE_M()(2), SSxiE_M()(2), 1e-4);
    check("M jac left forward adjoint motion 4 taylor     ", SxiE_M()(3), SSxiE_M()(3), 1e-4);
    check("M jac left forward adjoint motion 5 taylor     ", SxiE_M()(4), SSxiE_M()(4), 1e-4);
    check("M jac left forward adjoint motion 6 taylor     ", SxiE_M()(5), SSxiE_M()(5), 1e-4);
    se3_tangent xiE_Z          = Z_eb | xiB;
    dual SZ_eb                 = Z_eb.plus_left(DeltatauE);
    se3_tangent SxiE_Z         = SZ_eb | xiB;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_left_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_Z = Jac_Z * DeltatauE();
    se3_tangent SSxiE_Z        = xiE_Z + DeltaxiE_Z;

    check("Z jac left forward adjoint motion 1 taylor     ", SxiE_Z()(0), SSxiE_Z()(0), 1e-4);
    check("Z jac left forward adjoint motion 2 taylor     ", SxiE_Z()(1), SSxiE_Z()(1), 1e-4);
    check("Z jac left forward adjoint motion 3 taylor     ", SxiE_Z()(2), SSxiE_Z()(2), 1e-4);
    check("Z jac left forward adjoint motion 4 taylor     ", SxiE_Z()(3), SSxiE_Z()(3), 1e-4);
    check("Z jac left forward adjoint motion 5 taylor     ", SxiE_Z()(4), SSxiE_Z()(4), 1e-4);
    check("Z jac left forward adjoint motion 6 taylor     ", SxiE_Z()(5), SSxiE_Z()(5), 1e-4);

    se3_tangent xiE_tau          = tau_eb | xiB;
    trfv Stau_eb                 = tau_eb.plus_left(DeltatauE);
    se3_tangent SxiE_tau         = Stau_eb | xiB;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_left_forward_adjoint_wrt_motion(xiB);
    Eigen::Vector6d DeltaxiE_tau = Jac_tau * DeltatauE();
    se3_tangent SSxiE_tau        = xiE_tau + DeltaxiE_tau;

    check("tau jac left forward adjoint motion 1 taylor   ", SxiE_tau()(0), SSxiE_tau()(0), 1e-4);
    check("tau jac left forward adjoint motion 2 taylor   ", SxiE_tau()(1), SSxiE_tau()(1), 1e-4);
    check("tau jac left forward adjoint motion 3 taylor   ", SxiE_tau()(2), SSxiE_tau()(2), 1e-4);
    check("tau jac left forward adjoint motion 4 taylor   ", SxiE_tau()(3), SSxiE_tau()(3), 1e-4);
    check("tau jac left forward adjoint motion 5 taylor   ", SxiE_tau()(4), SSxiE_tau()(4), 1e-4);
    check("tau jac left forward adjoint motion 6 taylor   ", SxiE_tau()(5), SSxiE_tau()(5), 1e-4);

} // closes test_jac_left_forward_adjoint_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_left_backward_adjoint_wrt_motion() {

    // It verifies the behavior of the motion classes jac_left_backward_adjoint_wrt_motion method, validating the following:
    // Ad(DeltatauE plus M) % xi = AdM % xi + J * DeltatauE

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector6d tau; tau << 0.004, -0.005, 0.007, 1e-3, -5e-4, 2e-3;
    trfv DeltatauE(tau);

    se3_tangent xiE(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    se3_tangent xiB_gq          = gq_eb % xiE;
    speu_rodrigues Sgq_eb       = gq_eb.plus_left(DeltatauE);
    se3_tangent SxiB_gq         = Sgq_eb % xiE;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_left_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_gq = Jac_gq * DeltatauE();
    se3_tangent SSxiB_gq        = xiB_gq + DeltaxiB_gq;

    check("gq jac left backward adjoint motion 1 taylor    ", SxiB_gq()(0), SSxiB_gq()(0), 1e-4);
    check("gq jac left backward adjoint motion 2 taylor    ", SxiB_gq()(1), SSxiB_gq()(1), 1e-4);
    check("gq jac left backward adjoint motion 3 taylor    ", SxiB_gq()(2), SSxiB_gq()(2), 1e-4);
    check("gq jac left backward adjoint motion 4 taylor    ", SxiB_gq()(3), SSxiB_gq()(3), 1e-4);
    check("gq jac left backward adjoint motion 5 taylor    ", SxiB_gq()(4), SSxiB_gq()(4), 1e-4);
    check("gq jac left backward adjoint motion 6 taylor    ", SxiB_gq()(5), SSxiB_gq()(5), 1e-4);

    se3_tangent xiB_gR          = gR_eb % xiE;
    speu_dcm SgR_eb             = gR_eb.plus_left(DeltatauE);
    se3_tangent SxiB_gR         = SgR_eb % xiE;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_left_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_gR = Jac_gR * DeltatauE();
    se3_tangent SSxiB_gR        = xiB_gR + DeltaxiB_gR;

    check("gR jac left backward adjoint motion 1 taylor    ", SxiB_gR()(0), SSxiB_gR()(0), 1e-4);
    check("gR jac left backward adjoint motion 2 taylor    ", SxiB_gR()(1), SSxiB_gR()(1), 1e-4);
    check("gR jac left backward adjoint motion 3 taylor    ", SxiB_gR()(2), SSxiB_gR()(2), 1e-4);
    check("gR jac left backward adjoint motion 4 taylor    ", SxiB_gR()(3), SSxiB_gR()(3), 1e-4);
    check("gR jac left backward adjoint motion 5 taylor    ", SxiB_gR()(4), SSxiB_gR()(4), 1e-4);
    check("gR jac left backward adjoint motion 6 taylor    ", SxiB_gR()(5), SSxiB_gR()(5), 1e-4);

    se3_tangent xiB_M          = M_eb % xiE;
    homogeneous SM_eb          = M_eb.plus_left(DeltatauE);
    se3_tangent SxiB_M         = SM_eb % xiE;
    Eigen::Matrix6d Jac_M      = M_eb.jac_left_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_M = Jac_M * DeltatauE();
    se3_tangent SSxiB_M        = xiB_M + DeltaxiB_M;

    check("M jac left backward adjoint motion 1 taylor     ", SxiB_M()(0), SSxiB_M()(0), 1e-4);
    check("M jac left backward adjoint motion 2 taylor     ", SxiB_M()(1), SSxiB_M()(1), 1e-4);
    check("M jac left backward adjoint motion 3 taylor     ", SxiB_M()(2), SSxiB_M()(2), 1e-4);
    check("M jac left backward adjoint motion 4 taylor     ", SxiB_M()(3), SSxiB_M()(3), 1e-4);
    check("M jac left backward adjoint motion 5 taylor     ", SxiB_M()(4), SSxiB_M()(4), 1e-4);
    check("M jac left backward adjoint motion 6 taylor     ", SxiB_M()(5), SSxiB_M()(5), 1e-4);
    se3_tangent xiB_Z          = Z_eb % xiE;
    dual SZ_eb                 = Z_eb.plus_left(DeltatauE);
    se3_tangent SxiB_Z         = SZ_eb % xiE;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_left_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_Z = Jac_Z * DeltatauE();
    se3_tangent SSxiB_Z        = xiB_Z + DeltaxiB_Z;

    check("Z jac left backward adjoint motion 1 taylor     ", SxiB_Z()(0), SSxiB_Z()(0), 1e-4);
    check("Z jac left backward adjoint motion 2 taylor     ", SxiB_Z()(1), SSxiB_Z()(1), 1e-4);
    check("Z jac left backward adjoint motion 3 taylor     ", SxiB_Z()(2), SSxiB_Z()(2), 1e-4);
    check("Z jac left backward adjoint motion 4 taylor     ", SxiB_Z()(3), SSxiB_Z()(3), 1e-4);
    check("Z jac left backward adjoint motion 5 taylor     ", SxiB_Z()(4), SSxiB_Z()(4), 1e-4);
    check("Z jac left backward adjoint motion 6 taylor     ", SxiB_Z()(5), SSxiB_Z()(5), 1e-4);

    se3_tangent xiB_tau          = tau_eb % xiE;
    trfv Stau_eb                 = tau_eb.plus_left(DeltatauE);
    se3_tangent SxiB_tau         = Stau_eb % xiE;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_left_backward_adjoint_wrt_motion(xiE);
    Eigen::Vector6d DeltaxiB_tau = Jac_tau * DeltatauE();
    se3_tangent SSxiB_tau        = xiB_tau + DeltaxiB_tau;

    check("tau jac left backward adjoint motion 1 taylor   ", SxiB_tau()(0), SSxiB_tau()(0), 1e-4);
    check("tau jac left backward adjoint motion 2 taylor   ", SxiB_tau()(1), SSxiB_tau()(1), 1e-4);
    check("tau jac left backward adjoint motion 3 taylor   ", SxiB_tau()(2), SSxiB_tau()(2), 1e-4);
    check("tau jac left backward adjoint motion 4 taylor   ", SxiB_tau()(3), SSxiB_tau()(3), 1e-4);
    check("tau jac left backward adjoint motion 5 taylor   ", SxiB_tau()(4), SSxiB_tau()(4), 1e-4);
    check("tau jac left backward adjoint motion 6 taylor   ", SxiB_tau()(5), SSxiB_tau()(5), 1e-4);

} // closes test_jac_left_backward_adjoint_wrt_motion

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_forward_motion_wrt_point() {

    // It verifies the behavior of the motion classes jac_euclidean_forward_motion_wrt_point method, validating the following:
    // M * (p + Deltap) = M * p + J * Deltap

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector3d pB(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltapB(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d pE_gq      = gq_eb * pB;
    Eigen::Vector3d SpB_gq     = pB + DeltapB;
    Eigen::Vector3d SpE_gq     = gq_eb * SpB_gq;
    Eigen::Matrix3d Jac_gq     = gq_eb.jac_euclidean_forward_motion_wrt_point();
    Eigen::Vector3d DeltapE_gq = Jac_gq * DeltapB;
    Eigen::Vector3d SSpE_gq    = pE_gq + DeltapE_gq;

    check("gq jac euclidean forward motion point 1 taylor    ", SpE_gq(0), SSpE_gq(0), 1e-12);
    check("gq jac euclidean forward motion point 2 taylor    ", SpE_gq(1), SSpE_gq(1), 1e-12);
    check("gq jac euclidean forward motion point 3 taylor    ", SpE_gq(2), SSpE_gq(2), 1e-12);

    Eigen::Vector3d pE_gR      = gR_eb * pB;
    Eigen::Vector3d SpB_gR     = pB + DeltapB;
    Eigen::Vector3d SpE_gR     = gR_eb * SpB_gR;
    Eigen::Matrix3d Jac_gR     = gR_eb.jac_euclidean_forward_motion_wrt_point();
    Eigen::Vector3d DeltapE_gR = Jac_gR * DeltapB;
    Eigen::Vector3d SSpE_gR    = pE_gR + DeltapE_gR;

    check("gR jac euclidean forward motion point 1 taylor    ", SpE_gR(0), SSpE_gR(0), 1e-12);
    check("gR jac euclidean forward motion point 2 taylor    ", SpE_gR(1), SSpE_gR(1), 1e-12);
    check("gR jac euclidean forward motion point 3 taylor    ", SpE_gR(2), SSpE_gR(2), 1e-12);

    Eigen::Vector3d pE_M      = M_eb * pB;
    Eigen::Vector3d SpB_M     = pB + DeltapB;
    Eigen::Vector3d SpE_M     = M_eb * SpB_M;
    Eigen::Matrix3d Jac_M     = M_eb.jac_euclidean_forward_motion_wrt_point();
    Eigen::Vector3d DeltapE_M = Jac_M * DeltapB;
    Eigen::Vector3d SSpE_M    = pE_M + DeltapE_M;

    check("M jac euclidean forward motion point 1 taylor     ", SpE_M(0), SSpE_M(0), 1e-12);
    check("M jac euclidean forward motion point 2 taylor     ", SpE_M(1), SSpE_M(1), 1e-12);
    check("M jac euclidean forward motion point 3 taylor     ", SpE_M(2), SSpE_M(2), 1e-12);


    Eigen::Vector3d pE_Z      = Z_eb * pB;
    Eigen::Vector3d SpB_Z     = pB + DeltapB;
    Eigen::Vector3d SpE_Z     = Z_eb * SpB_Z;
    Eigen::Matrix3d Jac_Z     = Z_eb.jac_euclidean_forward_motion_wrt_point();
    Eigen::Vector3d DeltapE_Z = Jac_Z * DeltapB;
    Eigen::Vector3d SSpE_Z    = pE_Z + DeltapE_Z;

    check("Z jac euclidean forward motion point 1 taylor     ", SpE_Z(0), SSpE_Z(0), 1e-12);
    check("Z jac euclidean forward motion point 2 taylor     ", SpE_Z(1), SSpE_Z(1), 1e-12);
    check("Z jac euclidean forward motion point 3 taylor     ", SpE_Z(2), SSpE_Z(2), 1e-12);

    Eigen::Vector3d pE_tau      = tau_eb * pB;
    Eigen::Vector3d SpB_tau     = pB + DeltapB;
    Eigen::Vector3d SpE_tau     = tau_eb * SpB_tau;
    Eigen::Matrix3d Jac_tau     = tau_eb.jac_euclidean_forward_motion_wrt_point();
    Eigen::Vector3d DeltapE_tau = Jac_tau * DeltapB;
    Eigen::Vector3d SSpE_tau    = pE_tau + DeltapE_tau;

    check("tau jac euclidean forward motion point 1 taylor   ", SpE_tau(0), SSpE_tau(0), 1e-12);
    check("tau jac euclidean forward motion point 2 taylor   ", SpE_tau(1), SSpE_tau(1), 1e-12);
    check("tau jac euclidean forward motion point 3 taylor   ", SpE_tau(2), SSpE_tau(2), 1e-12);
} // closes test_jac_euclidean_forward_motion_wrt_point

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_backward_motion_wrt_point() {

    // It verifies the behavior of the motion classes jac_euclidean_backward_motion_wrt_point method, validating the following:
    // M / (p + Deltap) = M / p + J * Deltap

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    Eigen::Vector3d pE(0.33, -0.65, 0.17);

    Eigen::Vector3d DeltapE(1e-3, -5e-4, 2e-3);

    Eigen::Vector3d pB_gq      = gq_eb / pE;
    Eigen::Vector3d SpE_gq     = pE + DeltapE;
    Eigen::Vector3d SpB_gq     = gq_eb / SpE_gq;
    Eigen::Matrix3d Jac_gq     = gq_eb.jac_euclidean_backward_motion_wrt_point();
    Eigen::Vector3d DeltapB_gq = Jac_gq * DeltapE;
    Eigen::Vector3d SSpB_gq    = pB_gq + DeltapB_gq;

    check("gq jac euclidean backward motion point 1 taylor    ", SpB_gq(0), SSpB_gq(0), 1e-12);
    check("gq jac euclidean backward motion point 2 taylor    ", SpB_gq(1), SSpB_gq(1), 1e-12);
    check("gq jac euclidean backward motion point 3 taylor    ", SpB_gq(2), SSpB_gq(2), 1e-12);

    Eigen::Vector3d pB_gR      = gR_eb / pE;
    Eigen::Vector3d SpE_gR     = pE + DeltapE;
    Eigen::Vector3d SpB_gR     = gR_eb / SpE_gR;
    Eigen::Matrix3d Jac_gR     = gR_eb.jac_euclidean_backward_motion_wrt_point();
    Eigen::Vector3d DeltapB_gR = Jac_gR * DeltapE;
    Eigen::Vector3d SSpB_gR    = pB_gR + DeltapB_gR;

    check("gR jac euclidean backward motion point 1 taylor    ", SpB_gR(0), SSpB_gR(0), 1e-12);
    check("gR jac euclidean backward motion point 2 taylor    ", SpB_gR(1), SSpB_gR(1), 1e-12);
    check("gR jac euclidean backward motion point 3 taylor    ", SpB_gR(2), SSpB_gR(2), 1e-12);

    Eigen::Vector3d pB_M      = M_eb / pE;
    Eigen::Vector3d SpE_M     = pE + DeltapE;
    Eigen::Vector3d SpB_M     = M_eb / SpE_M;
    Eigen::Matrix3d Jac_M     = M_eb.jac_euclidean_backward_motion_wrt_point();
    Eigen::Vector3d DeltapB_M = Jac_M * DeltapE;
    Eigen::Vector3d SSpB_M    = pB_M + DeltapB_M;

    check("M jac euclidean backward motion point 1 taylor     ", SpB_M(0), SSpB_M(0), 1e-12);
    check("M jac euclidean backward motion point 2 taylor     ", SpB_M(1), SSpB_M(1), 1e-12);
    check("M jac euclidean backward motion point 3 taylor     ", SpB_M(2), SSpB_M(2), 1e-12);

    Eigen::Vector3d pB_Z      = Z_eb / pE;
    Eigen::Vector3d SpE_Z     = pE + DeltapE;
    Eigen::Vector3d SpB_Z     = Z_eb / SpE_Z;
    Eigen::Matrix3d Jac_Z     = Z_eb.jac_euclidean_backward_motion_wrt_point();
    Eigen::Vector3d DeltapB_Z = Jac_Z * DeltapE;
    Eigen::Vector3d SSpB_Z    = pB_Z + DeltapB_Z;

    check("Z jac euclidean backward motion point 1 taylor     ", SpB_Z(0), SSpB_Z(0), 1e-12);
    check("Z jac euclidean backward motion point 2 taylor     ", SpB_Z(1), SSpB_Z(1), 1e-12);
    check("Z jac euclidean backward motion point 3 taylor     ", SpB_Z(2), SSpB_Z(2), 1e-12);

    Eigen::Vector3d pB_tau      = tau_eb / pE;
    Eigen::Vector3d SpE_tau     = pE + DeltapE;
    Eigen::Vector3d SpB_tau     = tau_eb / SpE_tau;
    Eigen::Matrix3d Jac_tau     = tau_eb.jac_euclidean_backward_motion_wrt_point();
    Eigen::Vector3d DeltapB_tau = Jac_tau * DeltapE;
    Eigen::Vector3d SSpB_tau    = pB_tau + DeltapB_tau;

    check("tau jac euclidean backward motion point 1 taylor   ", SpB_tau(0), SSpB_tau(0), 1e-12);
    check("tau jac euclidean backward motion point 2 taylor   ", SpB_tau(1), SSpB_tau(1), 1e-12);
    check("tau jac euclidean backward motion point 3 taylor   ", SpB_tau(2), SSpB_tau(2), 1e-12);
} // closes test_jac_euclidean_backward_motion_wrt_point

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_forward_adjoint_wrt_tangent() {

    // It verifies the behavior of the motion classes jac_euclidean_forward_adjoint_wrt_tangent method, validating the following:
    // AdM | (xi + Deltaxi) = AdM | xi + J * Deltaxi

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    se3_tangent xiB(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    Eigen::Vector6d xi; xi << 1e-3, -5e-4, 2e-3, -6e-4, 3e-3, 8e-3;
    Eigen::Vector6d DeltaxiB(xi);

    se3_tangent xiE_gq          = gq_eb | xiB;
    se3_tangent SxiB_gq         = xiB + DeltaxiB;
    se3_tangent SxiE_gq         = gq_eb | SxiB_gq;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiE_gq = Jac_gq * DeltaxiB;
    se3_tangent  SSxiE_gq       = xiE_gq + DeltaxiE_gq;

    check("gq jac euclidean forward adjoint tangent 1 taylor    ", SxiE_gq()(0), SSxiE_gq()(0), 1e-12);
    check("gq jac euclidean forward adjoint tangent 2 taylor    ", SxiE_gq()(1), SSxiE_gq()(1), 1e-12);
    check("gq jac euclidean forward adjoint tangent 3 taylor    ", SxiE_gq()(2), SSxiE_gq()(2), 1e-12);
    check("gq jac euclidean forward adjoint tangent 4 taylor    ", SxiE_gq()(3), SSxiE_gq()(3), 1e-12);
    check("gq jac euclidean forward adjoint tangent 5 taylor    ", SxiE_gq()(4), SSxiE_gq()(4), 1e-12);
    check("gq jac euclidean forward adjoint tangent 6 taylor    ", SxiE_gq()(5), SSxiE_gq()(5), 1e-12);

    se3_tangent xiE_gR          = gR_eb | xiB;
    se3_tangent SxiB_gR         = xiB + DeltaxiB;
    se3_tangent SxiE_gR         = gR_eb | SxiB_gR;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiE_gR = Jac_gR * DeltaxiB;
    se3_tangent SSxiE_gR        = xiE_gR + DeltaxiE_gR;

    check("gR jac euclidean forward adjoint tangent 1 taylor    ", SxiE_gR()(0), SSxiE_gR()(0), 1e-12);
    check("gR jac euclidean forward adjoint tangent 2 taylor    ", SxiE_gR()(1), SSxiE_gR()(1), 1e-12);
    check("gR jac euclidean forward adjoint tangent 3 taylor    ", SxiE_gR()(2), SSxiE_gR()(2), 1e-12);
    check("gR jac euclidean forward adjoint tangent 4 taylor    ", SxiE_gR()(3), SSxiE_gR()(3), 1e-12);
    check("gR jac euclidean forward adjoint tangent 5 taylor    ", SxiE_gR()(4), SSxiE_gR()(4), 1e-12);
    check("gR jac euclidean forward adjoint tangent 6 taylor    ", SxiE_gR()(5), SSxiE_gR()(5), 1e-12);

    se3_tangent xiE_M          = M_eb | xiB;
    se3_tangent SxiB_M         = xiB + DeltaxiB;
    se3_tangent SxiE_M         = M_eb | SxiB_M;
    Eigen::Matrix6d Jac_M      = M_eb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiE_M = Jac_M * DeltaxiB;
    se3_tangent SSxiE_M        = xiE_M + DeltaxiE_M;

    check("M jac euclidean forward adjoint tangent 1 taylor     ", SxiE_M()(0), SSxiE_M()(0), 1e-12);
    check("M jac euclidean forward adjoint tangent 2 taylor     ", SxiE_M()(1), SSxiE_M()(1), 1e-12);
    check("M jac euclidean forward adjoint tangent 3 taylor     ", SxiE_M()(2), SSxiE_M()(2), 1e-12);
    check("M jac euclidean forward adjoint tangent 4 taylor     ", SxiE_M()(3), SSxiE_M()(3), 1e-12);
    check("M jac euclidean forward adjoint tangent 5 taylor     ", SxiE_M()(4), SSxiE_M()(4), 1e-12);
    check("M jac euclidean forward adjoint tangent 6 taylor     ", SxiE_M()(5), SSxiE_M()(5), 1e-12);

    se3_tangent xiE_Z          = Z_eb | xiB;
    se3_tangent SxiB_Z         = xiB + DeltaxiB;
    se3_tangent SxiE_Z         = Z_eb | SxiB_Z;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiE_Z = Jac_Z * DeltaxiB;
    se3_tangent SSxiE_Z        = xiE_Z + DeltaxiE_Z;

    check("Z jac euclidean forward adjoint tangent 1 taylor     ", SxiE_Z()(0), SSxiE_Z()(0), 1e-12);
    check("Z jac euclidean forward adjoint tangent 2 taylor     ", SxiE_Z()(1), SSxiE_Z()(1), 1e-12);
    check("Z jac euclidean forward adjoint tangent 3 taylor     ", SxiE_Z()(2), SSxiE_Z()(2), 1e-12);
    check("Z jac euclidean forward adjoint tangent 4 taylor     ", SxiE_Z()(3), SSxiE_Z()(3), 1e-12);
    check("Z jac euclidean forward adjoint tangent 5 taylor     ", SxiE_Z()(4), SSxiE_Z()(4), 1e-12);
    check("Z jac euclidean forward adjoint tangent 6 taylor     ", SxiE_Z()(5), SSxiE_Z()(5), 1e-12);

    se3_tangent xiE_tau          = tau_eb | xiB;
    se3_tangent SxiB_tau         = xiB + DeltaxiB;
    se3_tangent SxiE_tau         = tau_eb | SxiB_tau;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_euclidean_forward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiE_tau = Jac_tau * DeltaxiB;
    se3_tangent SSxiE_tau        = xiE_tau + DeltaxiE_tau;

    check("tau jac euclidean forward adjoint tangent 1 taylor   ", SxiE_tau()(0), SSxiE_tau()(0), 1e-12);
    check("tau jac euclidean forward adjoint tangent 2 taylor   ", SxiE_tau()(1), SSxiE_tau()(1), 1e-12);
    check("tau jac euclidean forward adjoint tangent 3 taylor   ", SxiE_tau()(2), SSxiE_tau()(2), 1e-12);
    check("tau jac euclidean forward adjoint tangent 4 taylor   ", SxiE_tau()(3), SSxiE_tau()(3), 1e-12);
    check("tau jac euclidean forward adjoint tangent 5 taylor   ", SxiE_tau()(4), SSxiE_tau()(4), 1e-12);
    check("tau jac euclidean forward adjoint tangent 6 taylor   ", SxiE_tau()(5), SSxiE_tau()(5), 1e-12);

} // closes test_jac_euclidean_forward_adjoint_wrt_tangent

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_backward_adjoint_wrt_tangent() {

    // It verifies the behavior of the motion classes jac_euclidean_backward_adjoint_wrt_tangent method, validating the following:
    // AdM % (xi + Deltaxi) = AdM % xi + J * Deltaxi

    rotv r_eb(-0.7, 0.45, 1.0);
    Eigen::Vector3d T_ebe(7.5, -12.8, -3.2);

    trfv tau_eb(r_eb, T_ebe);
    speu_dcm gR_eb(tau_eb);
    speu_rodrigues gq_eb(tau_eb);
    homogeneous M_eb(tau_eb);
    dual Z_eb(tau_eb);

    se3_tangent xiE(0.33, -0.65, 0.17, 1.11, 0.17, -0.58);

    Eigen::Vector6d xi; xi << 1e-3, -5e-4, 2e-3, -6e-4, 3e-3, 8e-3;
    Eigen::Vector6d DeltaxiE(xi);

    se3_tangent xiB_gq          = gq_eb % xiE;
    se3_tangent SxiE_gq         = xiE + DeltaxiE;
    se3_tangent SxiB_gq         = gq_eb % SxiE_gq;
    Eigen::Matrix6d Jac_gq      = gq_eb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiB_gq = Jac_gq * DeltaxiE;
    se3_tangent  SSxiB_gq       = xiB_gq + DeltaxiB_gq;

    check("gq jac euclidean backward adjoint tangent 1 taylor    ", SxiB_gq()(0), SSxiB_gq()(0), 1e-12);
    check("gq jac euclidean backward adjoint tangent 2 taylor    ", SxiB_gq()(1), SSxiB_gq()(1), 1e-12);
    check("gq jac euclidean backward adjoint tangent 3 taylor    ", SxiB_gq()(2), SSxiB_gq()(2), 1e-12);
    check("gq jac euclidean backward adjoint tangent 4 taylor    ", SxiB_gq()(3), SSxiB_gq()(3), 1e-12);
    check("gq jac euclidean backward adjoint tangent 5 taylor    ", SxiB_gq()(4), SSxiB_gq()(4), 1e-12);
    check("gq jac euclidean backward adjoint tangent 6 taylor    ", SxiB_gq()(5), SSxiB_gq()(5), 1e-12);

    se3_tangent xiB_gR          = gR_eb % xiE;
    se3_tangent SxiE_gR         = xiE + DeltaxiE;
    se3_tangent SxiB_gR         = gR_eb % SxiE_gR;
    Eigen::Matrix6d Jac_gR      = gR_eb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiB_gR = Jac_gR * DeltaxiE;
    se3_tangent SSxiB_gR        = xiB_gR + DeltaxiB_gR;

    check("gR jac euclidean backward adjoint tangent 1 taylor    ", SxiB_gR()(0), SSxiB_gR()(0), 1e-12);
    check("gR jac euclidean backward adjoint tangent 2 taylor    ", SxiB_gR()(1), SSxiB_gR()(1), 1e-12);
    check("gR jac euclidean backward adjoint tangent 3 taylor    ", SxiB_gR()(2), SSxiB_gR()(2), 1e-12);
    check("gR jac euclidean backward adjoint tangent 4 taylor    ", SxiB_gR()(3), SSxiB_gR()(3), 1e-12);
    check("gR jac euclidean backward adjoint tangent 5 taylor    ", SxiB_gR()(4), SSxiB_gR()(4), 1e-12);
    check("gR jac euclidean backward adjoint tangent 6 taylor    ", SxiB_gR()(5), SSxiB_gR()(5), 1e-12);

    se3_tangent xiB_M          = M_eb % xiE;
    se3_tangent SxiE_M         = xiE + DeltaxiE;
    se3_tangent SxiB_M         = M_eb % SxiE_M;
    Eigen::Matrix6d Jac_M      = M_eb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiB_M = Jac_M * DeltaxiE;
    se3_tangent SSxiB_M        = xiB_M + DeltaxiB_M;

    check("M jac euclidean backward adjoint tangent 1 taylor     ", SxiB_M()(0), SSxiB_M()(0), 1e-12);
    check("M jac euclidean backward adjoint tangent 2 taylor     ", SxiB_M()(1), SSxiB_M()(1), 1e-12);
    check("M jac euclidean backward adjoint tangent 3 taylor     ", SxiB_M()(2), SSxiB_M()(2), 1e-12);
    check("M jac euclidean backward adjoint tangent 4 taylor     ", SxiB_M()(3), SSxiB_M()(3), 1e-12);
    check("M jac euclidean backward adjoint tangent 5 taylor     ", SxiB_M()(4), SSxiB_M()(4), 1e-12);
    check("M jac euclidean backward adjoint tangent 6 taylor     ", SxiB_M()(5), SSxiB_M()(5), 1e-12);

    se3_tangent xiB_Z          = Z_eb % xiE;
    se3_tangent SxiE_Z         = xiE + DeltaxiE;
    se3_tangent SxiB_Z         = Z_eb % SxiE_Z;
    Eigen::Matrix6d Jac_Z      = Z_eb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiB_Z = Jac_Z * DeltaxiE;
    se3_tangent SSxiB_Z        = xiB_Z + DeltaxiB_Z;

    check("Z jac euclidean backward adjoint tangent 1 taylor     ", SxiB_Z()(0), SSxiB_Z()(0), 1e-12);
    check("Z jac euclidean backward adjoint tangent 2 taylor     ", SxiB_Z()(1), SSxiB_Z()(1), 1e-12);
    check("Z jac euclidean backward adjoint tangent 3 taylor     ", SxiB_Z()(2), SSxiB_Z()(2), 1e-12);
    check("Z jac euclidean backward adjoint tangent 4 taylor     ", SxiB_Z()(3), SSxiB_Z()(3), 1e-12);
    check("Z jac euclidean backward adjoint tangent 5 taylor     ", SxiB_Z()(4), SSxiB_Z()(4), 1e-12);
    check("Z jac euclidean backward adjoint tangent 6 taylor     ", SxiB_Z()(5), SSxiB_Z()(5), 1e-12);

    se3_tangent xiB_tau          = tau_eb % xiE;
    se3_tangent SxiE_tau         = xiE + DeltaxiE;
    se3_tangent SxiB_tau         = tau_eb % SxiE_tau;
    Eigen::Matrix6d Jac_tau      = tau_eb.jac_euclidean_backward_adjoint_wrt_tangent();
    Eigen::Vector6d DeltaxiB_tau = Jac_tau * DeltaxiE;
    se3_tangent SSxiB_tau        = xiB_tau + DeltaxiB_tau;

    check("tau jac euclidean backward adjoint tangent 1 taylor   ", SxiB_tau()(0), SSxiB_tau()(0), 1e-12);
    check("tau jac euclidean backward adjoint tangent 2 taylor   ", SxiB_tau()(1), SSxiB_tau()(1), 1e-12);
    check("tau jac euclidean backward adjoint tangent 3 taylor   ", SxiB_tau()(2), SSxiB_tau()(2), 1e-12);
    check("tau jac euclidean backward adjoint tangent 4 taylor   ", SxiB_tau()(3), SSxiB_tau()(3), 1e-12);
    check("tau jac euclidean backward adjoint tangent 5 taylor   ", SxiB_tau()(4), SSxiB_tau()(4), 1e-12);
    check("tau jac euclidean backward adjoint tangent 6 taylor   ", SxiB_tau()(5), SSxiB_tau()(5), 1e-12);

} // closes test_jac_euclidean_backward_adjoint_wrt_tangent

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_speu_dcm() {

    // It verifies the behavior of the speu_rodrigues class jac_euclidean_forward_motion_wrt_speu_rodrigues and
    // jac_backward_motion_wrt_speu_rodrigues, validating the following equations:
    // [Gr * exp(Delta tau)] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[Gr * exp(Delta tau)] - Gr}
    // [exp(Delta tau) * Gr] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[exp(Delta tau) * Gr] - Gr}
    // [Gr * exp(Delta tau)] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[Gr * exp(Delta tau)] - Gr}
    // [exp(Delta tau) * Gr] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[exp(Delta tau) * Gr] - Gr}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta tau

    double d2r = math::constant::D2R();

    euler euler_eb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    dcm R_eb(euler_eb);
    Eigen::Vector3d T_ebe(7.2, -8.1, -3.4);
    speu_dcm GR_eb(R_eb, T_ebe);
    Eigen::Vector12d GR_eb_wedge = ang::speu_dcm::wedge(GR_eb);

    rotv Deltar_eb( 0.001, -0.002, 0.0015); // WATCH OUT - needs to be small
    Eigen::Vector3d DeltaT_ebe(-0.003, -0.0012, 0.0018);
    ang::trfv Deltatau_eb(Deltar_eb, DeltaT_ebe);
    speu_dcm DeltaGR_eb(Deltatau_eb);
    Eigen::Vector12d DeltaGR_eb_wedge = ang::speu_dcm::wedge(DeltaGR_eb);

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_e = GR_eb * p_b;

    speu_dcm          SGR_eb       = GR_eb.plus_right(Deltatau_eb); // WATCH OUT - same as GR_eb * DeltaGR_eb
    Eigen::Vector12d  SGR_eb_wedge = ang::speu_dcm::wedge(SGR_eb);
    Eigen::Vector3d   Sp_e         = SGR_eb * p_b;
    Eigen::Matrix312d Jac_f        = GR_eb.jac_euclidean_forward_motion_wrt_speu_dcm(p_b);
    Eigen::Vector3d   SDelta_p_e   = Jac_f * (SGR_eb_wedge - GR_eb_wedge); // WATCH OUT - very different from Jac_f * speu_dcm::wedge(DeltaGR_eb)
    Eigen::Vector3d   SSp_e        = p_e + SDelta_p_e;

    check("Gq jac forward speu dcm 0 taylor  ", Sp_e(0), SSp_e(0), 1e-5);
    check("Gq jac forward speu dcm 1 taylor  ", Sp_e(1), SSp_e(1), 1e-6);
    check("Gq jac forward speu dcm 2 taylor  ", Sp_e(2), SSp_e(2), 1e-4);

    speu_dcm         TGR_eb       = GR_eb.plus_left(Deltatau_eb);
    Eigen::Vector12d TGR_eb_wedge = ang::speu_dcm::wedge(TGR_eb);
    Eigen::Vector3d  Tp_e         = TGR_eb * p_b;
    Eigen::Vector3d  TDelta_p_e   = Jac_f * (TGR_eb_wedge - GR_eb_wedge);
    Eigen::Vector3d  TTp_e        = p_e + TDelta_p_e;

    check("Gq jac forward speu dcm 0 taylor  ", Tp_e(0), TTp_e(0), 1e-5);
    check("Gq jac forward speu dcm 1 taylor  ", Tp_e(1), TTp_e(1), 1e-6);
    check("Gq jac forward speu dcm 2 taylor  ", Tp_e(2), TTp_e(2), 1e-4);

    Eigen::Vector3d   Sp_b       = SGR_eb / p_e;
    Eigen::Matrix312d Jac_b      = GR_eb.jac_euclidean_backward_motion_wrt_speu_dcm(p_e);
    Eigen::Vector3d   SDelta_p_b = Jac_b * (SGR_eb_wedge - GR_eb_wedge);
    Eigen::Vector3d   SSp_b      = p_b + SDelta_p_b;

    check("Gq jac backward speu dcm 0 taylor ", Sp_b(0), SSp_b(0), 1e-5);
    check("Gq jac backward speu dcm 1 taylor ", Sp_b(1), SSp_b(1), 1e-5);
    check("Gq jac backward speu dcm 2 taylor ", Sp_b(2), SSp_b(2), 1e-5);

    Eigen::Vector3d  Tp_b       = TGR_eb / p_e;
    Eigen::Vector3d  TDelta_p_b = Jac_b * (TGR_eb_wedge - GR_eb_wedge);
    Eigen::Vector3d  TTp_b      = p_b + TDelta_p_b;

    check("Gq jac backward speu dcm 0 taylor ", Tp_b(0), TTp_b(0), 1e-4);
    check("Gq jac backward speu dcm 1 taylor ", Tp_b(1), TTp_b(1), 1e-4);
    check("Gq jac backward speu dcm 2 taylor ", Tp_b(2), TTp_b(2), 1e-4);
} // closes test_jac_euclidean_wrt_speu_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_speu_dcm_linear() {

    // It verifies the behavior of the speu_dcm class jac_euclidean_forward_motion_wrt_speu_dcm method,
    // validating the following equation:
    // Gr * p = d(Gr * p)/dGr |Gr*p * p
    // Note that unlike the forward transform, the backward transform is not linear, so
    // no similar expression exists

    double d2r = math::constant::D2R();

    ang::euler euler_eb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    ang::dcm R_eb(euler_eb);
    Eigen::Vector3d T_ebe(7.2, -8.1, -3.4);
    ang::speu_dcm GR_eb(R_eb, T_ebe);
    Eigen::Vector12d GR_eb_wedge = ang::speu_dcm::wedge(GR_eb);

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_e = GR_eb * p_b;

    Eigen::Matrix312d Jac_f   = GR_eb.jac_euclidean_forward_motion_wrt_speu_dcm(p_b);
    Eigen::Vector3d   p_e2    = Jac_f * GR_eb_wedge;

    check("Gr jac forward speu dcm 0 linear    ", p_e(0), p_e2(0), 1e-14);
    check("Gr jac forward speu dcm 1 linear    ", p_e(1), p_e2(1), 1e-14);
    check("Gr jac forward speu dcm 2 linear    ", p_e(2), p_e2(2), 1e-14);

} // closes test_jac_euclidean_wrt_speu_dcm_Linear

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_speu_rodrigues() {

    // It verifies the behavior of the speu_rodrigues class jac_euclidean_forward_motion_wrt_speu_rodrigues and
    // jac_backward_motion_wrt_speu_rodrigues, validating the following equations:
    // [Gq * exp(Delta tau)] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[Gq * exp(Delta tau)] - Gq}
    // [exp(Delta tau) * Gq] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[exp(Delta tau) * Gq] - Gq}
    // [Gq * exp(Delta tau)] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[Gq * exp(Delta tau)] - Gq}
    // [exp(Delta tau) * Gq] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[exp(Delta tau) * Gq] - Gq}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta tau

    double d2r = math::constant::D2R();

    euler euler_eb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    rodrigues q_eb(euler_eb);
    Eigen::Vector3d T_ebe(7.2, -8.1, -3.4);
    speu_rodrigues Gq_eb(q_eb, T_ebe);
    Eigen::Vector7d Gq_eb_wedge = ang::speu_rodrigues::wedge(Gq_eb);

    rotv Deltar_eb( 0.001, -0.002, 0.0015); // WATCH OUT - needs to be small
    Eigen::Vector3d DeltaT_ebe(-0.003, -0.0012, 0.0018);
    ang::trfv Deltatau_eb(Deltar_eb, DeltaT_ebe);
    speu_rodrigues DeltaGq_eb(Deltatau_eb);
    Eigen::Vector7d DeltaGq_eb_wedge = ang::speu_rodrigues::wedge(DeltaGq_eb);

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_e = Gq_eb * p_b;

    speu_rodrigues   SGq_eb       = Gq_eb.plus_right(Deltatau_eb); // WATCH OUT - same as Gq_eb * DeltaGq_eb
    Eigen::Vector7d  SGq_eb_wedge = ang::speu_rodrigues::wedge(SGq_eb);
    Eigen::Vector3d  Sp_e         = SGq_eb * p_b;
    Eigen::Matrix37d Jac_f        = Gq_eb.jac_euclidean_forward_motion_wrt_speu_rodrigues(p_b);
    Eigen::Vector3d  SDelta_p_e   = Jac_f * (SGq_eb_wedge - Gq_eb_wedge); // WATCH OUT - very different from Jac_f * speu_rodrigues::wedge(DeltaGq_eb)
    Eigen::Vector3d  SSp_e        = p_e + SDelta_p_e;

    check("Gq jac forward speu rodrigues 0 taylor  ", Sp_e(0), SSp_e(0), 1e-5);
    check("Gq jac forward speu rodrigues 1 taylor  ", Sp_e(1), SSp_e(1), 1e-6);
    check("Gq jac forward speu rodrigues 2 taylor  ", Sp_e(2), SSp_e(2), 1e-4);

    speu_rodrigues   TGq_eb       = Gq_eb.plus_left(Deltatau_eb);
    Eigen::Vector7d  TGq_eb_wedge = ang::speu_rodrigues::wedge(TGq_eb);
    Eigen::Vector3d  Tp_e         = TGq_eb * p_b;
    Eigen::Vector3d  TDelta_p_e   = Jac_f * (TGq_eb_wedge - Gq_eb_wedge);
    Eigen::Vector3d  TTp_e        = p_e + TDelta_p_e;

    check("Gq jac forward speu rodrigues 0 taylor  ", Tp_e(0), TTp_e(0), 1e-5);
    check("Gq jac forward speu rodrigues 1 taylor  ", Tp_e(1), TTp_e(1), 1e-6);
    check("Gq jac forward speu rodrigues 2 taylor  ", Tp_e(2), TTp_e(2), 1e-4);

    Eigen::Vector3d  Sp_b       = SGq_eb / p_e;
    Eigen::Matrix37d Jac_b      = Gq_eb.jac_euclidean_backward_motion_wrt_speu_rodrigues(p_e);
    Eigen::Vector3d  SDelta_p_b = Jac_b * (SGq_eb_wedge - Gq_eb_wedge);
    Eigen::Vector3d  SSp_b      = p_b + SDelta_p_b;

    check("Gq jac backward speu rodrigues 0 taylor ", Sp_b(0), SSp_b(0), 1e-5);
    check("Gq jac backward speu rodrigues 1 taylor ", Sp_b(1), SSp_b(1), 1e-5);
    check("Gq jac backward speu rodrigues 2 taylor ", Sp_b(2), SSp_b(2), 1e-5);

    Eigen::Vector3d  Tp_b       = TGq_eb / p_e;
    Eigen::Vector3d  TDelta_p_b = Jac_b * (TGq_eb_wedge - Gq_eb_wedge);
    Eigen::Vector3d  TTp_b      = p_b + TDelta_p_b;

    check("Gq jac backward speu rodrigues 0 taylor ", Tp_b(0), TTp_b(0), 1e-4);
    check("Gq jac backward speu rodrigues 1 taylor ", Tp_b(1), TTp_b(1), 1e-5);
    check("Gq jac backward speu rodrigues 2 taylor ", Tp_b(2), TTp_b(2), 1e-4);

} // closes test_jac_euclidean_wrt_speu_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_homogeneous() {

    // It verifies the behavior of the homogeneous class jac_euclidean_forward_motion_wrt_homogeneous and
    // jac_backward_motion_wrt_homogeneous, validating the following equations:
    // [M * exp(Delta tau)] * p ~= M * p + d(M * p)/dM |M*p * {[M * exp(Delta tau)] - M}
    // [exp(Delta tau) * M] * p ~= M * p + d(M * p)/dM |M*p * {[exp(Delta tau) * M] - M}
    // [M * exp(Delta tau)] / p ~= M / p + d(M / p)/dM |M/p * {[M * exp(Delta tau)] - M}
    // [exp(Delta tau) * M] / p ~= M / p + d(M / p)/dM |M/p * {[exp(Delta tau) * M] - M}
    // Note that ~= means that in the first order Taylor expansion these are only valid for small Delta tau

    double d2r = math::constant::D2R();

    euler euler_eb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    dcm R_eb(euler_eb);
    Eigen::Vector3d T_ebe(7.2, -8.1, -3.4);
    homogeneous M_eb(R_eb, T_ebe);
    Eigen::Vector12d M_eb_wedge = ang::homogeneous::wedge(M_eb);

    rotv Deltar_eb( 0.001, -0.002, 0.0015); // WATCH OUT - needs to be small
    Eigen::Vector3d DeltaT_ebe(-0.003, -0.0012, 0.0018); // WATCH OUT - needs to be small
    ang::trfv Deltatau_eb(Deltar_eb, DeltaT_ebe);
    homogeneous DeltaM_eb(Deltatau_eb);
    Eigen::Vector12d DeltaM_eb_wedge = ang::homogeneous::wedge(DeltaM_eb);

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_e = M_eb * p_b;

    homogeneous       SM_eb        = M_eb.plus_right(Deltatau_eb); // WATCH OUT - same as M_eb * DeltaM_eb
    Eigen::Vector12d  SM_eb_wedge  = ang::homogeneous::wedge(SM_eb);
    Eigen::Vector3d   Sp_e         = SM_eb * p_b;
    Eigen::Matrix312d Jac_f        = M_eb.jac_euclidean_forward_motion_wrt_homogeneous(p_b);
    Eigen::Vector3d   SDelta_p_e   = Jac_f * (SM_eb_wedge - M_eb_wedge); // WATCH OUT - very different from Jac_f * homogeneous::wedge(DeltaM_eb)
    Eigen::Vector3d   SSp_e        = p_e + SDelta_p_e;

    check("M jac forward homogeneous 0 taylor  ", Sp_e(0), SSp_e(0), 1e-5);
    check("M jac forward homogeneous 1 taylor  ", Sp_e(1), SSp_e(1), 1e-6);
    check("M jac forward homogeneous 2 taylor  ", Sp_e(2), SSp_e(2), 1e-4);

    homogeneous      TM_eb        = M_eb.plus_left(Deltatau_eb);
    Eigen::Vector12d TM_eb_wedge  = ang::homogeneous::wedge(TM_eb);
    Eigen::Vector3d  Tp_e         = TM_eb * p_b;
    Eigen::Vector3d  TDelta_p_e   = Jac_f * (TM_eb_wedge - M_eb_wedge);
    Eigen::Vector3d  TTp_e        = p_e + TDelta_p_e;

    check("M jac forward homogeneous 0 taylor  ", Tp_e(0), TTp_e(0), 1e-5);
    check("M jac forward homogeneous 1 taylor  ", Tp_e(1), TTp_e(1), 1e-6);
    check("M jac forward homogeneous 2 taylor  ", Tp_e(2), TTp_e(2), 1e-4);

    Eigen::Vector3d   Sp_b       = SM_eb / p_e;
    Eigen::Matrix312d Jac_b      = M_eb.jac_euclidean_backward_motion_wrt_homogeneous(p_e);
    Eigen::Vector3d   SDelta_p_b = Jac_b * (SM_eb_wedge - M_eb_wedge);
    Eigen::Vector3d   SSp_b      = p_b + SDelta_p_b;

    check("M jac backward homogeneous 0 taylor ", Sp_b(0), SSp_b(0), 1e-5);
    check("M jac backward homogeneous 1 taylor ", Sp_b(1), SSp_b(1), 1e-5);
    check("M jac backward homogeneous 2 taylor ", Sp_b(2), SSp_b(2), 1e-5);

    Eigen::Vector3d  Tp_b       = TM_eb / p_e;
    Eigen::Vector3d  TDelta_p_b = Jac_b * (TM_eb_wedge - M_eb_wedge);
    Eigen::Vector3d  TTp_b      = p_b + TDelta_p_b;

    check("Mjac backward homogeneous 0 taylor ", Tp_b(0), TTp_b(0), 1e-4);
    check("Mjac backward homogeneous 1 taylor ", Tp_b(1), TTp_b(1), 1e-4);
    check("Mjac backward homogeneous 2 taylor ", Tp_b(2), TTp_b(2), 1e-4);

} // closes test_jac_euclidean_wrt_homogeneous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_homogeneous_linear() {

    // It verifies the behavior of the homogeneous class jac_euclidean_forward_motion_wrt_homogeneous method,
    // validating the following equation:
    // M * p = d(M * p)/dM |M*p * M
    // Note that unlike the forward transform, the backward transform is not linear, so
    // no similar expression exists

    double d2r = math::constant::D2R();

    euler euler_eb(-10.7 * d2r, -5.5 * d2r, 2.1 * d2r);
    dcm R_eb(euler_eb);
    Eigen::Vector3d T_ebe(7.2, -8.1, -3.4);
    homogeneous M_eb(R_eb, T_ebe);
    Eigen::Vector12d M_eb_wedge = ang::homogeneous::wedge(M_eb);

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_e = M_eb * p_b;

    Eigen::Matrix312d Jac_f   = M_eb.jac_euclidean_forward_motion_wrt_homogeneous(p_b);
    Eigen::Vector3d   p_e2    = Jac_f * M_eb_wedge;

    check("M jac forward homogeneous 0 linear    ", p_e(0), p_e2(0), 1e-14);
    check("M jac forward homogeneous 1 linear    ", p_e(1), p_e2(1), 1e-14);
    check("M jac forward homogeneous 2 linear    ", p_e(2), p_e2(2), 1e-14);

} // closes test_jac_euclidean_wrt_homogeneous_linear

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_trfv() {

    // It checks the proper behavior of the "trfv" class methods related with the euclidean action jacobians,
    // validating the following equations:
    // (tau + Delta tau) * p ~= tau * p + d(tau * p)/dtau |tau*p * Delta tau
    // (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau

    Eigen::Vector6d Itau; Itau << 4.0, 11.0, 7.0, 0.3, 0.5, 0.2;
    ang::trfv tau(Itau);

    Eigen::Vector6d IDelta_tau; IDelta_tau << 0.0003, -0.0008, 0.0002, -0.0002, -0.00031, 0.00013; // needs to be small
    ang::trfv Delta_tau(IDelta_tau);

    ang::trfv Stau(tau() + Delta_tau());

    Eigen::Vector3d p_b(3.8, -4.1, -2.3);
    Eigen::Vector3d p_x = tau * p_b;
    Eigen::Vector3d Sp_x  = Stau * p_b;
    Eigen::Vector3d Sp_x2 = p_x + tau.jac_euclidean_forward_motion_wrt_trfv(p_b) * Delta_tau();
    Eigen::Vector3d Sp_x3 = p_x + tau.jac_euclidean_forward_motion_wrt_trfv_bis(p_b) * Delta_tau();
    Eigen::Vector3d Sp_x4 = p_x + tau.jac_euclidean_forward_motion_wrt_trfv_tri(p_b) * Delta_tau();

    check("jac sum trfv 0     ", Sp_x(0), Sp_x2(0), 1e-6);
    check("jac sum trfv 1     ", Sp_x(1), Sp_x2(1), 1e-6);
    check("jac sum trfv 2     ", Sp_x(2), Sp_x2(2), 1e-6);
    check("jac sum trfv 0     ", Sp_x(0), Sp_x3(0), 1e-6);
    check("jac sum trfv 1     ", Sp_x(1), Sp_x3(1), 1e-6);
    check("jac sum trfv 2     ", Sp_x(2), Sp_x3(2), 1e-6);
    check("jac sum trfv 0     ", Sp_x(0), Sp_x4(0), 1e-6);
    check("jac sum trfv 1     ", Sp_x(1), Sp_x4(1), 1e-6);
    check("jac sum trfv 2     ", Sp_x(2), Sp_x4(2), 1e-6);

    Eigen::Vector3d Sp_b  = Stau / p_x;
    Eigen::Vector3d Sp_b1 = p_b + tau.jac_euclidean_backward_motion_wrt_trfv_zero(p_x) * Delta_tau();
    Eigen::Vector3d Sp_b2 = p_b + tau.jac_euclidean_backward_motion_wrt_trfv(p_x) * Delta_tau();
    Eigen::Vector3d Sp_b3 = p_b + tau.jac_euclidean_backward_motion_wrt_trfv_bis(p_x) * Delta_tau();
    Eigen::Vector3d Sp_b4 = p_b + tau.jac_euclidean_backward_motion_wrt_trfv_tri(p_x) * Delta_tau();

    check("jac sum trfv 0     ", Sp_b(0), Sp_b1(0), 1e-6);
    check("jac sum trfv 1     ", Sp_b(1), Sp_b1(1), 1e-6);
    check("jac sum trfv 2     ", Sp_b(2), Sp_b1(2), 1e-6);
    check("jac sum trfv 0     ", Sp_b(0), Sp_b2(0), 1e-6);
    check("jac sum trfv 1     ", Sp_b(1), Sp_b2(1), 1e-6);
    check("jac sum trfv 2     ", Sp_b(2), Sp_b2(2), 1e-6);
    check("jac sum trfv 0     ", Sp_b(0), Sp_b3(0), 1e-6);
    check("jac sum trfv 1     ", Sp_b(1), Sp_b3(1), 1e-7);
    check("jac sum trfv 2     ", Sp_b(2), Sp_b3(2), 1e-6);
    check("jac sum trfv 0     ", Sp_b(0), Sp_b4(0), 1e-6);
    check("jac sum trfv 1     ", Sp_b(1), Sp_b4(1), 1e-7);
    check("jac sum trfv 2     ", Sp_b(2), Sp_b4(2), 1e-6);

} // closes test_jac_euclidean_wrt_trfv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3_jacobian::test_jac_euclidean_wrt_trfv_bis() {

    // It verifies the behavior of the transform vector class jacobian methods, validating the following expressions:
    // exp(tau + Deltatau).get_T() ~= exp(tau).get_T + d(exp(tau).get_T())/drv |[exp(tau)] * Deltarv
    // exp(tau.inv() + Deltatau).get_T() ~= exp(tau.inv()).get_T + d(exp(tau.inv()).get_T())/drv |[exp(tau.inv())] * Deltarv
    // exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/drv |[exp(tau) * p] * Deltarv
    // exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) / p)/drv |[exp(tau)/p] * Deltarv
    // exp(tau + Deltatau) * p ==  exp(tau) * p + d(exp(tau) * p)/ds |[exp(tau)*p] * Deltas
    // exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) * p)/ds |[exp(tau)/p] * Deltas
    // exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/dtau |[exp(tau)*p] * Deltatau
    // exp(tau + Deltatau) / p ~=  exp(tau) / p + d(exp(tau) / p)/dtau |[exp(tau)/p] * Deltatau

    Eigen::Vector6d Itau_xb; Itau_xb << 4.0, 11.0, 7.0, 1.0, 0.5, 0.2;
    ang::trfv tau_xb(Itau_xb);

    Eigen::Vector3d Delta_rv_xb(-0.8e-5, 1.1e-5, -1.2e-5);
    Eigen::Vector6d Itau1_xb; Itau1_xb << Itau_xb(0), Itau_xb(1), Itau_xb(2), Itau_xb(3) + Delta_rv_xb(0), Itau_xb(4) + Delta_rv_xb(1), Itau_xb(5) + Delta_rv_xb(2);
    ang::trfv tau1_xb(Itau1_xb);

    Eigen::Vector3d Delta_s_xb(-0.3 * 1e-3, 0.4 * 1e-3, 0.7 * 1e-3);
    Eigen::Vector6d Itau2_xb; Itau2_xb << Itau_xb(0) + Delta_s_xb(0), Itau_xb(1) + Delta_s_xb(1), Itau_xb(2) + Delta_s_xb(2), Itau_xb(3), Itau_xb(4), Itau_xb(5);
    ang::trfv tau2_xb(Itau2_xb);

    Eigen::Vector6d Itau3_xb; Itau3_xb << Itau_xb(0) + Delta_s_xb(0), Itau_xb(1) + Delta_s_xb(1), Itau_xb(2) + Delta_s_xb(2), Itau_xb(3) + Delta_rv_xb(0), Itau_xb(4) + Delta_rv_xb(1), Itau_xb(5) + Delta_rv_xb(2);
    ang::trfv tau3_xb(Itau3_xb);
    Eigen::Vector6d Delta_tau_xb; Delta_tau_xb << Delta_s_xb, Delta_rv_xb;

    Eigen::Vector3d p_b(7.1, -10.2, 4.3);
    Eigen::Vector3d p_x = tau_xb * p_b;
    Eigen::Vector3d p_b2 = tau_xb / p_x;

    Eigen::Vector3d T11   = tau_xb.get_T();
    Eigen::Matrix3d J11   = tau_xb.jac_euclidean_forward_motion_T_wrt_rotv();
    Eigen::Vector3d TT11  = tau1_xb.get_T();
    Eigen::Vector3d TTT11 = T11 + J11 * Delta_rv_xb;

    check("jacobian T 1       ", TT11(0), TTT11(0), 1e-9);
    check("jacobian T 2       ", TT11(1), TTT11(1), 1e-9);
    check("jacobian T 3       ", TT11(2), TTT11(2), 1e-9);

    Eigen::Matrix3d JX11  = tau_xb.jac_euclidean_forward_motion_wrt_rotv(p_b);
    Eigen::Vector3d OpXX  = tau1_xb * p_b;
    Eigen::Vector3d OpXXX = p_x + JX11 * Delta_rv_xb;

    check("jacobian tau 1     ", OpXX(0), OpXXX(0), 1e-9);
    check("jacobian tau 2     ", OpXX(1), OpXXX(1), 1e-9);
    check("jacobian tau 3     ", OpXX(2), OpXXX(2), 1e-9);

    Eigen::Vector3d Tinv11   = tau_xb.get_inverse_T();
    Eigen::Matrix3d Jinv11   = tau_xb.jac_euclidean_backward_motion_T_wrt_rotv();
    Eigen::Vector3d TTinv11  = tau1_xb.get_inverse_T();
    Eigen::Vector3d TTTinv11 = Tinv11 + Jinv11 * Delta_rv_xb;

    check("jacobian T inv 1   ", TTinv11(0), TTTinv11(0), 1e-9);
    check("jacobian T inv 2   ", TTinv11(1), TTTinv11(1), 1e-9);
    check("jacobian T inv 3   ", TTinv11(2), TTTinv11(2), 1e-9);

    Eigen::Matrix3d _JX11  = tau_xb.jac_euclidean_backward_motion_wrt_rotv(p_x);
    Eigen::Vector3d _OpBB  = tau1_xb / p_x;
    Eigen::Vector3d _OpBBB = p_b2 + _JX11 * Delta_rv_xb;

    check("jacobian tau inv 1 ", _OpBB(0), _OpBBB(0), 1e-8);
    check("jacobian tau inv 2 ", _OpBB(1), _OpBBB(1), 1e-9);
    check("jacobian tau inv 3 ", _OpBB(2), _OpBBB(2), 1e-9);

    Eigen::Matrix3d AJX11  = tau_xb.jac_euclidean_forward_motion_wrt_s(p_b);
    Eigen::Vector3d ApXX   = tau2_xb * p_b;
    Eigen::Vector3d ApXXX  = p_x + AJX11 * Delta_s_xb;

    check("jacobian tau 1     ", ApXX(0), ApXXX(0), 1e-14);
    check("jacobian tau 2     ", ApXX(1), ApXXX(1), 1e-14);
    check("jacobian tau 3     ", ApXX(2), ApXXX(2), 1e-14);

    Eigen::Matrix3d _AJX11 = tau_xb.jac_euclidean_backward_motion_wrt_s(p_x);
    Eigen::Vector3d _ApBB  = tau2_xb / p_x;
    Eigen::Vector3d _ApBBB = p_b2 + _AJX11 * Delta_s_xb;

    check("jacobian tau inv 1 ", _ApBB(0), _ApBBB(0), 1e-14);
    check("jacobian tau inv 2 ", _ApBB(1), _ApBBB(1), 1e-14);
    check("jacobian tau inv 3 ", _ApBB(2), _ApBBB(2), 1e-14);

    Eigen::Vector3d BpXX            = tau3_xb * p_b;
    Eigen::Vector3d p_x_sum4 = p_x + tau_xb.jac_euclidean_forward_motion_wrt_trfv_tri(p_b) * Delta_tau_xb;

    check("jacobian tau 1     ", BpXX(0), p_x_sum4(0), 1e-8);
    check("jacobian tau 2     ", BpXX(1), p_x_sum4(1), 1e-8);
    check("jacobian tau 3     ", BpXX(2), p_x_sum4(2), 1e-8);

    Eigen::Vector3d _BpBB            = tau3_xb / p_x;
    Eigen::Vector3d p_b_sum4 = p_b2 + tau_xb.jac_euclidean_backward_motion_wrt_trfv_tri(p_x) * Delta_tau_xb;

    check("jacobian tau inv 1 ", _BpBB(0), p_b_sum4(0), 1e-8);
    check("jacobian tau inv 2 ", _BpBB(1), p_b_sum4(1), 1e-8);
    check("jacobian tau inv 3 ", _BpBB(2), p_b_sum4(2), 1e-8);

} // closes test_jac_euclidean_wrt_trfv_bis

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
























