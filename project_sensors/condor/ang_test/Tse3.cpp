#include "Tse3.h"
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

ang::test::Tse3::Tse3(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tse3::run() {
	::jail::unit_test::run();

    test_se3();		                       cout << endl << endl;
    test_exp_log_maps();                   cout << endl << endl;
    test_power();                          cout << endl << endl;
    test_sclerp();                         cout << endl << endl;
    test_plus_minus();                     cout << endl << endl;
    test_adjoint();                        cout << endl << endl;
    test_velocity();                       cout << endl << endl;
    test_dual_set();                       cout << endl << endl;

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_se3() {
    double d2r = math::constant::D2R();

    // different rotation directions (not unitary)
    vector<Eigen::Vector3d> Vdir_full(26);
    Vdir_full[0]  << +1.,  0.,  0.;
    Vdir_full[1]  << -1.,  0.,  0.;
    Vdir_full[2]  <<  0., +1.,  0.;
    Vdir_full[3]  <<  0., -1.,  0.;
    Vdir_full[4]  <<  0.,  0., +1.;
    Vdir_full[5]  <<  0.,  0., -1.;
    Vdir_full[6]  << +1., +2.,  0.;
    Vdir_full[7]  << +1., -2.,  0.;
    Vdir_full[8]  << -1., +2.,  0.;
    Vdir_full[9]  << -1., -2.,  0.;
    Vdir_full[10] << +1.,  0., +2.;
    Vdir_full[11] << +1.,  0., -2.;
    Vdir_full[12] << -1.,  0., +2.;
    Vdir_full[13] << -1.,  0., -2.;
    Vdir_full[14] <<  0., +1., +2.;
    Vdir_full[15] <<  0., +1., -2.;
    Vdir_full[16] <<  0., -1., +2.;
    Vdir_full[17] <<  0., -1., -2.;
    Vdir_full[18] << +1., +2., +3.;
    Vdir_full[19] << +1., +2., -3.;
    Vdir_full[20] << +1., -2., +3.;
    Vdir_full[21] << +1., -2., -3.;
    Vdir_full[22] << -1., +2., +3.;
    Vdir_full[23] << -1., +2., -3.;
    Vdir_full[24] << -1., -2., +3.;
    Vdir_full[25] << -1., -2., -3.;

    // different rotation angles (less than 180 [deg] as otherwise shortest path is chosen)
    vector<double> Vangle_rad(11);
    Vangle_rad[0]  = 0.7 * math::constant::PI(); // more than  90 deg
    Vangle_rad[1]  = 0.2 * math::constant::PI(); // less than 90 deg
    Vangle_rad[2]  = 1e-3;
    Vangle_rad[3]  = 1e-5;
    Vangle_rad[4]  = 1e-7;
    Vangle_rad[5]  = 1e-9;
    Vangle_rad[6]  = 1e-11;
    Vangle_rad[7]  = 1e-13;
    Vangle_rad[8]  = 1e-15;
    Vangle_rad[9]  = 1e-16;
    Vangle_rad[10] = 0.;

    Eigen::Vector3d T_xbx(50.0, 75.0, -33.0); // translation
    Eigen::Vector3d p_x(-4.0, 3.0, 7.0); // input point
    Eigen::Vector3d v_x(3.0, -8.0, 2.0); // input vector
    ang::se3_tangent xi_xbb_mrps(0.5, -0.3, 0.8, 0.2*d2r, -0.1*d2r, 0.3*d2r);

    Eigen::Vector4d p4_x; p4_x << p_x, 1.;
    Eigen::Vector4d p4_b_M, p4P_b_M;
    Eigen::Vector4d p4_x_M, p4P_x_M;

    Eigen::Vector3d  dir_unit;
    ang::rodrigues gq1X_xb, gq2X_xb, gq3X_xb, gq4X_xb, gq5X_xb, gq6X_xb;
    Eigen::Vector3d p_b_q, p_b_R, p_b_M, p_b_tau, p_b_S, p_b_Z;
    Eigen::Vector3d pP_b_q, pP_b_R, pP_b_M, pP_b_tau, pP_b_S, pP_b_Z;
    Eigen::Vector3d p_x_q, p_x_R, p_x_M, p_x_tau, p_x_S, p_x_Z;
    Eigen::Vector3d pP_x_q, pP_x_R, pP_x_M, pP_x_tau, pP_x_S, pP_x_Z;

    Eigen::Vector3d vR_b_q, vR_b_R, vR_b_M, vR_b_tau, vR_b_S, vR_b_Z;
    Eigen::Vector3d vPR_b_q, vPR_b_R, vPR_b_M, vPR_b_tau, vPR_b_S, vPR_b_Z;
    Eigen::Vector3d vR_x_q, vR_x_R, vR_x_M, vR_x_tau, vR_x_S, vR_x_Z;
    Eigen::Vector3d vPR_x_q, vPR_x_R, vPR_x_M, vPR_x_tau, vPR_x_S, vPR_x_Z;

    Eigen::Matrix<double,7,1> gqdot1_xb, gqdot2_xb;
    Eigen::Matrix<double,3,4> gRdot1_xb, gRdot2_xb;
    Eigen::Matrix4d gMdot1_xb, gMdot2_xb;
    ang::dual_quat gZdot1_xb, gZdot2_xb;
    ang::se3_tangent xi_xbx_mrps;

    ang::se3_tangent xi_xbb_mrps_homogeneous, xi_xbx_mrps_homogeneous, xi_xbb_mrps_speu_dcm, xi_xbx_mrps_speu_dcm, xi_xbb_mrps_speu_rodrigues, xi_xbx_mrps_speu_rodrigues, xi_xbb_mrps_dual, xi_xbx_mrps_dual;

    ang::euler      euler_bc(18 * d2r, -27 * d2r, +44 * d2r);
    Eigen::Vector3d T_bcb(-17.0, 24.0, +34.0);
    ang::rodrigues  q_bc(euler_bc);
    ang::dcm        R_bc(euler_bc);
    ang::rotv       rv_bc(euler_bc);

    ang::speu_rodrigues gq_bc(q_bc, T_bcb);
    ang::speu_dcm       gR_bc(R_bc, T_bcb);
    ang::homogeneous    gM_bc(R_bc, T_bcb);
    ang::trfv           gtau_bc(rv_bc, T_bcb);
    ang::screw          gS_bc(rv_bc, T_bcb);
    ang::dual           gZ_bc(rv_bc, T_bcb);

    ang::speu_rodrigues gq_xc, gqX_xc, gq_bx;
    ang::speu_dcm       gR_xc, gRX_xc, gR_bx;
    ang::homogeneous    gM_xc, gMX_xc, gM_bx;
    ang::trfv           gtau_xc, gtauX_xc, gtau_bx;
    ang::screw          gS_xc, gSX_xc, gS_bx;
    ang::dual           gZ_xc, gZX_xc, gZ_bx;
    ang::rodrigues gqq_xc, gqR_xc, gqM_xc, gqtau_xc, gqS_xc, gqZ_xc, gqs_xc, gqqX_xc, gqRX_xc, gqMX_xc, gqtauX_xc, gqSX_xc, gqZX_xc, gqsX_xc;

    ang::rodrigues gZ1, gZ2, gZ3, gZ4, gZ5, gZ6;
    ang::rotv rvZ7, rvZ8, rvZ9, rvZ10, rvZ11, rvZ12, rvZ13, rvZ14;
    Eigen::Vector3d TZ1, TZ2, TZ3, TZ4, TZ5, TZ6, TZ7, TZ8, TZ9, TZ10, TZ11, TZ12, TZ13, TZ14;

    for (unsigned int i = 0; i != Vdir_full.size(); ++i) {
        dir_unit = Vdir_full[i] / Vdir_full[i].norm();
        for (unsigned int j = 0; j != Vangle_rad.size(); ++j) {

            cout << endl << "Iteration i = " << i << "  j = " << j << endl;

            ang::rotv rv_xb(dir_unit * Vangle_rad[j]);
            ang::dcm R_xb(rv_xb);
            ang::rodrigues q_xb(rv_xb);

            // create the different SE3 representations among each other
            ang::speu_rodrigues gq1_xb(q_xb, T_xbx);
            ang::speu_dcm gR1_xb(R_xb, T_xbx);
            ang::homogeneous gM1_xb(R_xb, T_xbx);
            ang::trfv gtau1_xb(rv_xb, T_xbx);
            ang::dual gZ1_xb(rv_xb, T_xbx);
            ang::screw gS1_xb(rv_xb, T_xbx);

            ang::speu_rodrigues gq2_xb(gR1_xb);
            ang::speu_rodrigues gq3_xb(gM1_xb);
            ang::speu_rodrigues gq4_xb(gtau1_xb);
            ang::speu_rodrigues gq5_xb(gZ1_xb);
            ang::speu_rodrigues gq6_xb(gS1_xb);

            ang::speu_dcm gR2_xb(gq1_xb);
            ang::speu_dcm gR3_xb(gM1_xb);
            ang::speu_dcm gR4_xb(gtau1_xb);
            ang::speu_dcm gR5_xb(gZ1_xb);
            ang::speu_dcm gR6_xb(gS1_xb);

            ang::homogeneous gM2_xb(gq1_xb);
            ang::homogeneous gM3_xb(gR1_xb);
            ang::homogeneous gM4_xb(gtau1_xb);
            ang::homogeneous gM5_xb(gZ1_xb);
            ang::homogeneous gM6_xb(gS1_xb);

            ang::trfv gtau2_xb(gq1_xb);
            ang::trfv gtau3_xb(gR1_xb);
            ang::trfv gtau4_xb(gM1_xb);
            ang::trfv gtau5_xb(gS1_xb);
            ang::trfv gtau6_xb(gZ1_xb);

            ang::dual gZ2_xb(gq1_xb);
            ang::dual gZ3_xb(gR1_xb);
            ang::dual gZ4_xb(gM1_xb);
            ang::dual gZ5_xb(gtau1_xb);
            ang::dual gZ6_xb(gS1_xb);

            ang::screw gS2_xb(gq1_xb);
            ang::screw gS3_xb(gR1_xb);
            ang::screw gS4_xb(gM1_xb);
            ang::screw gS5_xb(gtau1_xb);
            ang::screw gS6_xb(gZ1_xb);

            check("rot vector gR 1a   :", rv_xb()(0), gR1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2a   :", rv_xb()(1), gR1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3a   :", rv_xb()(2), gR1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gR 1b   :", rv_xb()(0), gR2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2b   :", rv_xb()(1), gR2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3b   :", rv_xb()(2), gR2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gR 1c   :", rv_xb()(0), gR3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2c   :", rv_xb()(1), gR3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3c   :", rv_xb()(2), gR3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gR 1d   :", rv_xb()(0), gR4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2d   :", rv_xb()(1), gR4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3d   :", rv_xb()(2), gR4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gR 1e   :", rv_xb()(0), gR5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2e   :", rv_xb()(1), gR5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3e   :", rv_xb()(2), gR5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gR 1f   :", rv_xb()(0), gR6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gR 2f   :", rv_xb()(1), gR6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gR 3f   :", rv_xb()(2), gR6_xb.get_rotv()()(2), 1e-12);

            check("rot vector gq 1a   :", rv_xb()(0), gq1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2a   :", rv_xb()(1), gq1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3a   :", rv_xb()(2), gq1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gq 1b   :", rv_xb()(0), gq2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2b   :", rv_xb()(1), gq2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3b   :", rv_xb()(2), gq2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gq 1c   :", rv_xb()(0), gq3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2c   :", rv_xb()(1), gq3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3c   :", rv_xb()(2), gq3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gq 1d   :", rv_xb()(0), gq4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2d   :", rv_xb()(1), gq4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3d   :", rv_xb()(2), gq4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gq 1e   :", rv_xb()(0), gq5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2e   :", rv_xb()(1), gq5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3e   :", rv_xb()(2), gq5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gq 1f   :", rv_xb()(0), gq6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gq 2f   :", rv_xb()(1), gq6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gq 3f   :", rv_xb()(2), gq6_xb.get_rotv()()(2), 1e-12);

            check("rot vector gM 1a   :", rv_xb()(0), gM1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2a   :", rv_xb()(1), gM1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3a   :", rv_xb()(2), gM1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gM 1b   :", rv_xb()(0), gM2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2b   :", rv_xb()(1), gM2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3b   :", rv_xb()(2), gM2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gM 1c   :", rv_xb()(0), gM3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2c   :", rv_xb()(1), gM3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3c   :", rv_xb()(2), gM3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gM 1d   :", rv_xb()(0), gM4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2d   :", rv_xb()(1), gM4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3d   :", rv_xb()(2), gM4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gM 1e   :", rv_xb()(0), gM5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2e   :", rv_xb()(1), gM5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3e   :", rv_xb()(2), gM5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gM 1f   :", rv_xb()(0), gM6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gM 2f   :", rv_xb()(1), gM6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gM 3f   :", rv_xb()(2), gM6_xb.get_rotv()()(2), 1e-12);

            check("rot vector gtau 1a   :", rv_xb()(0), gtau1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2a   :", rv_xb()(1), gtau1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3a   :", rv_xb()(2), gtau1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gtau 1b   :", rv_xb()(0), gtau2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2b   :", rv_xb()(1), gtau2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3b   :", rv_xb()(2), gtau2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gtau 1c   :", rv_xb()(0), gtau3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2c   :", rv_xb()(1), gtau3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3c   :", rv_xb()(2), gtau3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gtau 1d   :", rv_xb()(0), gtau4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2d   :", rv_xb()(1), gtau4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3d   :", rv_xb()(2), gtau4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gtau 1e   :", rv_xb()(0), gtau5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2e   :", rv_xb()(1), gtau5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3e   :", rv_xb()(2), gtau5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gtau 1f   :", rv_xb()(0), gtau6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gtau 2f   :", rv_xb()(1), gtau6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gtau 3f   :", rv_xb()(2), gtau6_xb.get_rotv()()(2), 1e-12);

            check("rot vector gZ 1a   :", rv_xb()(0), gZ1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2a   :", rv_xb()(1), gZ1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3a   :", rv_xb()(2), gZ1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gZ 1b   :", rv_xb()(0), gZ2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2b   :", rv_xb()(1), gZ2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3b   :", rv_xb()(2), gZ2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gZ 1c   :", rv_xb()(0), gZ3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2c   :", rv_xb()(1), gZ3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3c   :", rv_xb()(2), gZ3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gZ 1d   :", rv_xb()(0), gZ4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2d   :", rv_xb()(1), gZ4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3d   :", rv_xb()(2), gZ4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gZ 1e   :", rv_xb()(0), gZ5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2e   :", rv_xb()(1), gZ5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3e   :", rv_xb()(2), gZ5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gZ 1f   :", rv_xb()(0), gZ6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gZ 2f   :", rv_xb()(1), gZ6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gZ 3f   :", rv_xb()(2), gZ6_xb.get_rotv()()(2), 1e-12);

            check("rot vector gS 1a   :", rv_xb()(0), gS1_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2a   :", rv_xb()(1), gS1_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3a   :", rv_xb()(2), gS1_xb.get_rotv()()(2), 1e-12);
            check("rot vector gS 1b   :", rv_xb()(0), gS2_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2b   :", rv_xb()(1), gS2_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3b   :", rv_xb()(2), gS2_xb.get_rotv()()(2), 1e-12);
            check("rot vector gS 1c   :", rv_xb()(0), gS3_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2c   :", rv_xb()(1), gS3_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3c   :", rv_xb()(2), gS3_xb.get_rotv()()(2), 1e-12);
            check("rot vector gS 1d   :", rv_xb()(0), gS4_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2d   :", rv_xb()(1), gS4_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3d   :", rv_xb()(2), gS4_xb.get_rotv()()(2), 1e-12);
            check("rot vector gS 1e   :", rv_xb()(0), gS5_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2e   :", rv_xb()(1), gS5_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3e   :", rv_xb()(2), gS5_xb.get_rotv()()(2), 1e-12);
            check("rot vector gS 1f   :", rv_xb()(0), gS6_xb.get_rotv()()(0), 1e-12);
            check("rot vector gS 2f   :", rv_xb()(1), gS6_xb.get_rotv()()(1), 1e-12);
            check("rot vector gS 3f   :", rv_xb()(2), gS6_xb.get_rotv()()(2), 1e-12);

            // check translations
            check("translation gR 1a  :", T_xbx(0), gR1_xb.get_T()(0), 1e-12);
            check("translation gR 2a  :", T_xbx(1), gR1_xb.get_T()(1), 1e-12);
            check("translation gR 3a  :", T_xbx(2), gR1_xb.get_T()(2), 1e-12);
            check("translation gR 1b  :", T_xbx(0), gR2_xb.get_T()(0), 1e-12);
            check("translation gR 2b  :", T_xbx(1), gR2_xb.get_T()(1), 1e-12);
            check("translation gR 3b  :", T_xbx(2), gR2_xb.get_T()(2), 1e-12);
            check("translation gR 1c  :", T_xbx(0), gR3_xb.get_T()(0), 1e-12);
            check("translation gR 2c  :", T_xbx(1), gR3_xb.get_T()(1), 1e-12);
            check("translation gR 3c  :", T_xbx(2), gR3_xb.get_T()(2), 1e-12);
            check("translation gR 1d  :", T_xbx(0), gR4_xb.get_T()(0), 1e-12);
            check("translation gR 2d  :", T_xbx(1), gR4_xb.get_T()(1), 1e-12);
            check("translation gR 3d  :", T_xbx(2), gR4_xb.get_T()(2), 1e-12);
            check("translation gR 1e  :", T_xbx(0), gR5_xb.get_T()(0), 1e-12);
            check("translation gR 2e  :", T_xbx(1), gR5_xb.get_T()(1), 1e-12);
            check("translation gR 3e  :", T_xbx(2), gR5_xb.get_T()(2), 1e-12);
            check("translation gR 1f  :", T_xbx(0), gR6_xb.get_T()(0), 1e-12);
            check("translation gR 2f  :", T_xbx(1), gR6_xb.get_T()(1), 1e-12);
            check("translation gR 3f  :", T_xbx(2), gR6_xb.get_T()(2), 1e-12);

            check("translation gq 1a  :", T_xbx(0), gq1_xb.get_T()(0), 1e-12);
            check("translation gq 2a  :", T_xbx(1), gq1_xb.get_T()(1), 1e-12);
            check("translation gq 3a  :", T_xbx(2), gq1_xb.get_T()(2), 1e-12);
            check("translation gq 1b  :", T_xbx(0), gq2_xb.get_T()(0), 1e-12);
            check("translation gq 2b  :", T_xbx(1), gq2_xb.get_T()(1), 1e-12);
            check("translation gq 3b  :", T_xbx(2), gq2_xb.get_T()(2), 1e-12);
            check("translation gq 1c  :", T_xbx(0), gq3_xb.get_T()(0), 1e-12);
            check("translation gq 2c  :", T_xbx(1), gq3_xb.get_T()(1), 1e-12);
            check("translation gq 3c  :", T_xbx(2), gq3_xb.get_T()(2), 1e-12);
            check("translation gq 1d  :", T_xbx(0), gq4_xb.get_T()(0), 1e-12);
            check("translation gq 2d  :", T_xbx(1), gq4_xb.get_T()(1), 1e-12);
            check("translation gq 3d  :", T_xbx(2), gq4_xb.get_T()(2), 1e-12);
            check("translation gq 1e  :", T_xbx(0), gq5_xb.get_T()(0), 1e-12);
            check("translation gq 2e  :", T_xbx(1), gq5_xb.get_T()(1), 1e-12);
            check("translation gq 3e  :", T_xbx(2), gq5_xb.get_T()(2), 1e-12);
            check("translation gq 1f  :", T_xbx(0), gq6_xb.get_T()(0), 1e-12);
            check("translation gq 2f  :", T_xbx(1), gq6_xb.get_T()(1), 1e-12);
            check("translation gq 3f  :", T_xbx(2), gq6_xb.get_T()(2), 1e-12);

            check("translation gM 1a  :", T_xbx(0), gM1_xb.get_T()(0), 1e-12);
            check("translation gM 2a  :", T_xbx(1), gM1_xb.get_T()(1), 1e-12);
            check("translation gM 3a  :", T_xbx(2), gM1_xb.get_T()(2), 1e-12);
            check("translation gM 1b  :", T_xbx(0), gM2_xb.get_T()(0), 1e-12);
            check("translation gM 2b  :", T_xbx(1), gM2_xb.get_T()(1), 1e-12);
            check("translation gM 3b  :", T_xbx(2), gM2_xb.get_T()(2), 1e-12);
            check("translation gM 1c  :", T_xbx(0), gM3_xb.get_T()(0), 1e-12);
            check("translation gM 2c  :", T_xbx(1), gM3_xb.get_T()(1), 1e-12);
            check("translation gM 3c  :", T_xbx(2), gM3_xb.get_T()(2), 1e-12);
            check("translation gM 1d  :", T_xbx(0), gM4_xb.get_T()(0), 1e-12);
            check("translation gM 2d  :", T_xbx(1), gM4_xb.get_T()(1), 1e-12);
            check("translation gM 3d  :", T_xbx(2), gM4_xb.get_T()(2), 1e-12);
            check("translation gM 1e  :", T_xbx(0), gM5_xb.get_T()(0), 1e-12);
            check("translation gM 2e  :", T_xbx(1), gM5_xb.get_T()(1), 1e-12);
            check("translation gM 3e  :", T_xbx(2), gM5_xb.get_T()(2), 1e-12);
            check("translation gM 1f  :", T_xbx(0), gM6_xb.get_T()(0), 1e-12);
            check("translation gM 2f  :", T_xbx(1), gM6_xb.get_T()(1), 1e-12);
            check("translation gM 3f  :", T_xbx(2), gM6_xb.get_T()(2), 1e-12);

            check("translation gtau 1a  :", T_xbx(0), gtau1_xb.get_T()(0), 1e-12);
            check("translation gtau 2a  :", T_xbx(1), gtau1_xb.get_T()(1), 1e-12);
            check("translation gtau 3a  :", T_xbx(2), gtau1_xb.get_T()(2), 1e-12);
            check("translation gtau 1b  :", T_xbx(0), gtau2_xb.get_T()(0), 1e-12);
            check("translation gtau 2b  :", T_xbx(1), gtau2_xb.get_T()(1), 1e-12);
            check("translation gtau 3b  :", T_xbx(2), gtau2_xb.get_T()(2), 1e-12);
            check("translation gtau 1c  :", T_xbx(0), gtau3_xb.get_T()(0), 1e-12);
            check("translation gtau 2c  :", T_xbx(1), gtau3_xb.get_T()(1), 1e-12);
            check("translation gtau 3c  :", T_xbx(2), gtau3_xb.get_T()(2), 1e-12);
            check("translation gtau 1d  :", T_xbx(0), gtau4_xb.get_T()(0), 1e-12);
            check("translation gtau 2d  :", T_xbx(1), gtau4_xb.get_T()(1), 1e-12);
            check("translation gtau 3d  :", T_xbx(2), gtau4_xb.get_T()(2), 1e-12);
            check("translation gtau 1e  :", T_xbx(0), gtau5_xb.get_T()(0), 1e-12);
            check("translation gtau 2e  :", T_xbx(1), gtau5_xb.get_T()(1), 1e-12);
            check("translation gtau 3e  :", T_xbx(2), gtau5_xb.get_T()(2), 1e-12);
            check("translation gtau 1f  :", T_xbx(0), gtau6_xb.get_T()(0), 1e-12);
            check("translation gtau 2f  :", T_xbx(1), gtau6_xb.get_T()(1), 1e-12);
            check("translation gtau 3f  :", T_xbx(2), gtau6_xb.get_T()(2), 1e-12);

            check("translation gZ 1a  :", T_xbx(0), gZ1_xb.get_T()(0), 1e-12);
            check("translation gZ 2a  :", T_xbx(1), gZ1_xb.get_T()(1), 1e-12);
            check("translation gZ 3a  :", T_xbx(2), gZ1_xb.get_T()(2), 1e-12);
            check("translation gZ 1b  :", T_xbx(0), gZ2_xb.get_T()(0), 1e-12);
            check("translation gZ 2b  :", T_xbx(1), gZ2_xb.get_T()(1), 1e-12);
            check("translation gZ 3b  :", T_xbx(2), gZ2_xb.get_T()(2), 1e-12);
            check("translation gZ 1c  :", T_xbx(0), gZ3_xb.get_T()(0), 1e-12);
            check("translation gZ 2c  :", T_xbx(1), gZ3_xb.get_T()(1), 1e-12);
            check("translation gZ 3c  :", T_xbx(2), gZ3_xb.get_T()(2), 1e-12);
            check("translation gZ 1d  :", T_xbx(0), gZ4_xb.get_T()(0), 1e-12);
            check("translation gZ 2d  :", T_xbx(1), gZ4_xb.get_T()(1), 1e-12);
            check("translation gZ 3d  :", T_xbx(2), gZ4_xb.get_T()(2), 1e-12);
            check("translation gZ 1e  :", T_xbx(0), gZ5_xb.get_T()(0), 1e-12);
            check("translation gZ 2e  :", T_xbx(1), gZ5_xb.get_T()(1), 1e-12);
            check("translation gZ 3e  :", T_xbx(2), gZ5_xb.get_T()(2), 1e-12);
            check("translation gZ 1f  :", T_xbx(0), gZ6_xb.get_T()(0), 1e-12);
            check("translation gZ 2f  :", T_xbx(1), gZ6_xb.get_T()(1), 1e-12);
            check("translation gZ 3f  :", T_xbx(2), gZ6_xb.get_T()(2), 1e-12);

            check("translation gS 1a  :", T_xbx(0), gS1_xb.get_T()(0), 1e-12);
            check("translation gS 2a  :", T_xbx(1), gS1_xb.get_T()(1), 1e-12);
            check("translation gS 3a  :", T_xbx(2), gS1_xb.get_T()(2), 1e-12);
            check("translation gS 1b  :", T_xbx(0), gS2_xb.get_T()(0), 1e-12);
            check("translation gS 2b  :", T_xbx(1), gS2_xb.get_T()(1), 1e-12);
            check("translation gS 3b  :", T_xbx(2), gS2_xb.get_T()(2), 1e-12);
            check("translation gS 1c  :", T_xbx(0), gS3_xb.get_T()(0), 1e-12);
            check("translation gS 2c  :", T_xbx(1), gS3_xb.get_T()(1), 1e-12);
            check("translation gS 3c  :", T_xbx(2), gS3_xb.get_T()(2), 1e-12);
            check("translation gS 1d  :", T_xbx(0), gS4_xb.get_T()(0), 1e-12);
            check("translation gS 2d  :", T_xbx(1), gS4_xb.get_T()(1), 1e-12);
            check("translation gS 3d  :", T_xbx(2), gS4_xb.get_T()(2), 1e-12);
            check("translation gS 1e  :", T_xbx(0), gS5_xb.get_T()(0), 1e-12);
            check("translation gS 2e  :", T_xbx(1), gS5_xb.get_T()(1), 1e-12);
            check("translation gS 3e  :", T_xbx(2), gS5_xb.get_T()(2), 1e-12);
            check("translation gS 1f  :", T_xbx(0), gS6_xb.get_T()(0), 1e-12);
            check("translation gS 2f  :", T_xbx(1), gS6_xb.get_T()(1), 1e-12);
            check("translation gS 3f  :", T_xbx(2), gS6_xb.get_T()(2), 1e-12);

            // inverse transformations
            p_b_q     = gq1_xb / p_x;
            p_b_R     = gR1_xb / p_x;
            p_b_M     = gM1_xb / p_x;
            p_b_tau   = gtau1_xb / p_x;
            p_b_Z     = gZ1_xb / p_x;
            p_b_S     = gS1_xb / p_x;
            p4_b_M    = gM1_xb / p4_x;
            pP_b_q    = gq1_xb.inverse() * p_x;
            pP_b_R    = gR1_xb.inverse() * p_x;
            pP_b_M    = gM1_xb.inverse() * p_x;
            pP_b_tau  = gtau1_xb.inverse() * p_x;
            pP_b_Z    = gZ1_xb.inverse() * p_x;
            pP_b_S    = gS1_xb.inverse() * p_x;
            p4P_b_M   = gM1_xb.inverse() * p4_x;

            check("method / - 1a             ", p_b_q(0), p_b_R(0),   1e-12);
            check("method / - 2a             ", p_b_q(1), p_b_R(1),   1e-12);
            check("method / - 3a             ", p_b_q(2), p_b_R(2),   1e-12);
            check("method / - 1b             ", p_b_q(0), p_b_M(0),   1e-12);
            check("method / - 2b             ", p_b_q(1), p_b_M(1),   1e-12);
            check("method / - 3b             ", p_b_q(2), p_b_M(2),   1e-12);
            check("method / - 1c             ", p_b_q(0), p_b_tau(0), 1e-12);
            check("method / - 2c             ", p_b_q(1), p_b_tau(1), 1e-12);
            check("method / - 3c             ", p_b_q(2), p_b_tau(2), 1e-12);
            check("method / - 1d             ", p_b_q(0), p_b_Z(0),   1e-12);
            check("method / - 2d             ", p_b_q(1), p_b_Z(1),   1e-12);
            check("method / - 3d             ", p_b_q(2), p_b_Z(2),   1e-12);
            check("method / - 1e             ", p_b_q(0), p_b_S(0),   1e-12);
            check("method / - 2e             ", p_b_q(1), p_b_S(1),   1e-12);
            check("method / - 3e             ", p_b_q(2), p_b_S(2),   1e-12);
            check("method / - 1g             ", p_b_q(0), p4_b_M(0),  1e-12);
            check("method / - 2g             ", p_b_q(1), p4_b_M(1),  1e-12);
            check("method / - 3g             ", p_b_q(2), p4_b_M(2),  1e-12);
            check("method / - 4g             ", 1.,         p4_b_M(3),  1e-12);
            check("method inverse and * - 1a ", p_b_q(0), pP_b_q(0),  1e-12);
            check("method inverse and * - 2a ", p_b_q(1), pP_b_q(1),  1e-12);
            check("method inverse and * - 3a ", p_b_q(2), pP_b_q(2),  1e-12);
            check("method inverse and * - 1b ", p_b_q(0), pP_b_R(0),  1e-12);
            check("method inverse and * - 2b ", p_b_q(1), pP_b_R(1),  1e-12);
            check("method inverse and * - 3b ", p_b_q(2), pP_b_R(2),  1e-12);
            check("method inverse and * - 1c ", p_b_q(0), pP_b_M(0),  1e-12);
            check("method inverse and * - 2c ", p_b_q(1), pP_b_M(1),  1e-12);
            check("method inverse and * - 3c ", p_b_q(2), pP_b_M(2),  1e-12);
            check("method inverse and * - 1d ", p_b_q(0), pP_b_tau(0),1e-12);
            check("method inverse and * - 2d ", p_b_q(1), pP_b_tau(1),1e-12);
            check("method inverse and * - 3d ", p_b_q(2), pP_b_tau(2),1e-12);
            check("method inverse and * - 1e ", p_b_q(0), pP_b_Z(0),  1e-12);
            check("method inverse and * - 2e ", p_b_q(1), pP_b_Z(1),  1e-12);
            check("method inverse and * - 3e ", p_b_q(2), pP_b_Z(2),  1e-12);
            check("method inverse and * - 1f ", p_b_q(0), pP_b_S(0),  1e-12);
            check("method inverse and * - 2f ", p_b_q(1), pP_b_S(1),  1e-12);
            check("method inverse and * - 3f ", p_b_q(2), pP_b_S(2),  1e-12);
            check("method inverse and * - 1h ", p_b_q(0), p4P_b_M(0), 1e-12);
            check("method inverse and * - 2h ", p_b_q(1), p4P_b_M(1), 1e-12);
            check("method inverse and * - 3h ", p_b_q(2), p4P_b_M(2), 1e-12);
            check("method inverse and * - 4h ", 1.,         p4P_b_M(3), 1e-12);

            // direct transformations
            p_x_q     = gq1_xb * p_b_q;
            p_x_R     = gR1_xb * p_b_R;
            p_x_M     = gM1_xb * p_b_M;
            p_x_tau   = gtau1_xb * p_b_tau;
            p_x_Z     = gZ1_xb * p_b_Z;
            p_x_S     = gS1_xb * p_b_S;
            p4_x_M    = gM1_xb * p4_b_M;
            pP_x_q    = gq1_xb.inverse() / pP_b_q;
            pP_x_R    = gR1_xb.inverse() / pP_b_R;
            pP_x_M    = gM1_xb.inverse() / pP_b_M;
            pP_x_tau  = gtau1_xb.inverse() / pP_b_tau;
            pP_x_Z    = gZ1_xb.inverse() / pP_b_Z;
            pP_x_S    = gS1_xb.inverse() / pP_b_S;
            p4P_x_M   = gM1_xb.inverse() / p4P_b_M;

            check("method * - 1a             ", p_x(0), p_x_q(0),    1e-12);
            check("method * - 2a             ", p_x(1), p_x_q(1),    1e-12);
            check("method * - 3a             ", p_x(2), p_x_q(2),    1e-12);
            check("method * - 1b             ", p_x(0), p_x_R(0),    1e-12);
            check("method * - 2b             ", p_x(1), p_x_R(1),    1e-12);
            check("method * - 3b             ", p_x(2), p_x_R(2),    1e-12);
            check("method * - 1c             ", p_x(0), p_x_M(0),    1e-12);
            check("method * - 2c             ", p_x(1), p_x_M(1),    1e-12);
            check("method * - 3c             ", p_x(2), p_x_M(2),    1e-12);
            check("method * - 1d             ", p_x(0), p_x_tau(0),  1e-12);
            check("method * - 2d             ", p_x(1), p_x_tau(1),  1e-12);
            check("method * - 3d             ", p_x(2), p_x_tau(2),  1e-12);
            check("method * - 1e             ", p_x(0), p_x_Z(0),    1e-12);
            check("method * - 2e             ", p_x(1), p_x_Z(1),    1e-12);
            check("method * - 3e             ", p_x(2), p_x_Z(2),    1e-12);
            check("method * - 1f             ", p_x(0), p_x_S(0),    1e-12);
            check("method * - 2f             ", p_x(1), p_x_S(1),    1e-12);
            check("method * - 3f             ", p_x(2), p_x_S(2),    1e-12);
            check("method * - 1h             ", p4_x(0), p4_x_M(0),  1e-12);
            check("method * - 2h             ", p4_x(1), p4_x_M(1),  1e-12);
            check("method * - 3h             ", p4_x(2), p4_x_M(2),  1e-12);
            check("method * - 4h             ", p4_x(3), p4_x_M(3),  1e-12);
            check("method inverse and / - 1a ", p_x(0), pP_x_q(0),   1e-12);
            check("method inverse and / - 2a ", p_x(1), pP_x_q(1),   1e-12);
            check("method inverse and / - 3a ", p_x(2), pP_x_q(2),   1e-12);
            check("method inverse and / - 1b ", p_x(0), pP_x_R(0),   1e-12);
            check("method inverse and / - 2b ", p_x(1), pP_x_R(1),   1e-12);
            check("method inverse and / - 3b ", p_x(2), pP_x_R(2),   1e-12);
            check("method inverse and / - 1c ", p_x(0), pP_x_M(0),   1e-12);
            check("method inverse and / - 2c ", p_x(1), pP_x_M(1),   1e-12);
            check("method inverse and / - 3c ", p_x(2), pP_x_M(2),   1e-12);
            check("method inverse and / - 1d ", p_x(0), pP_x_tau(0), 1e-12);
            check("method inverse and / - 2d ", p_x(1), pP_x_tau(1), 1e-12);
            check("method inverse and / - 3d ", p_x(2), pP_x_tau(2), 1e-12);
            check("method inverse and / - 1e ", p_x(0), pP_x_Z(0),   1e-12);
            check("method inverse and / - 2e ", p_x(1), pP_x_Z(1),   1e-12);
            check("method inverse and / - 3e ", p_x(2), pP_x_Z(2),   1e-12);
            check("method inverse and / - 1f ", p_x(0), pP_x_S(0),   1e-12);
            check("method inverse and / - 2f ", p_x(1), pP_x_S(1),   1e-12);
            check("method inverse and / - 3f ", p_x(2), pP_x_S(2),   1e-12);
            check("method inverse and / - 1h ", p4_x(0), p4P_x_M(0), 1e-12);
            check("method inverse and / - 2h ", p4_x(1), p4P_x_M(1), 1e-12);
            check("method inverse and / - 3h ", p4_x(2), p4P_x_M(2), 1e-12);
            check("method inverse and / - 4h ", p4_x(3), p4P_x_M(3), 1e-12);

            // inverse rotations only (no translations)
            vR_b_q    = gq1_xb & v_x;
            vR_b_R    = gR1_xb & v_x;
            vR_b_M    = gM1_xb & v_x;
            vR_b_tau  = gtau1_xb & v_x;
            vR_b_Z    = gZ1_xb & v_x;
            vR_b_S    = gS1_xb & v_x;
            vPR_b_q   = gq1_xb.inverse() ^ v_x;
            vPR_b_R   = gR1_xb.inverse() ^ v_x;
            vPR_b_M   = gM1_xb.inverse() ^ v_x;
            vPR_b_tau = gtau1_xb.inverse() ^ v_x;
            vPR_b_Z   = gZ1_xb.inverse() ^ v_x;
            vPR_b_S   = gS1_xb.inverse() ^ v_x;

            check("method & - 1a             ", vR_b_q(0), vR_b_R(0),   1e-12);
            check("method & - 2a             ", vR_b_q(1), vR_b_R(1),   1e-12);
            check("method & - 3a             ", vR_b_q(2), vR_b_R(2),   1e-12);
            check("method & - 1b             ", vR_b_q(0), vR_b_M(0),   1e-12);
            check("method & - 2b             ", vR_b_q(1), vR_b_M(1),   1e-12);
            check("method & - 3b             ", vR_b_q(2), vR_b_M(2),   1e-12);
            check("method & - 1c             ", vR_b_q(0), vR_b_tau(0), 1e-12);
            check("method & - 2c             ", vR_b_q(1), vR_b_tau(1), 1e-12);
            check("method & - 3c             ", vR_b_q(2), vR_b_tau(2), 1e-12);
            check("method & - 1d             ", vR_b_q(0), vR_b_Z(0),   1e-12);
            check("method & - 2d             ", vR_b_q(1), vR_b_Z(1),   1e-12);
            check("method & - 3d             ", vR_b_q(2), vR_b_Z(2),   1e-12);
            check("method & - 1e             ", vR_b_q(0), vR_b_S(0),   1e-12);
            check("method & - 2e             ", vR_b_q(1), vR_b_S(1),   1e-12);
            check("method & - 3e             ", vR_b_q(2), vR_b_S(2),   1e-12);
            check("method inverse and ^ - 1a ", vR_b_q(0), vPR_b_q(0),  1e-12);
            check("method inverse and ^ - 2a ", vR_b_q(1), vPR_b_q(1),  1e-12);
            check("method inverse and ^ - 3a ", vR_b_q(2), vPR_b_q(2),  1e-12);
            check("method inverse and ^ - 1b ", vR_b_q(0), vPR_b_R(0),  1e-12);
            check("method inverse and ^ - 2b ", vR_b_q(1), vPR_b_R(1),  1e-12);
            check("method inverse and ^ - 3b ", vR_b_q(2), vPR_b_R(2),  1e-12);
            check("method inverse and ^ - 1c ", vR_b_q(0), vPR_b_M(0),  1e-12);
            check("method inverse and ^ - 2c ", vR_b_q(1), vPR_b_M(1),  1e-12);
            check("method inverse and ^ - 3c ", vR_b_q(2), vPR_b_M(2),  1e-12);
            check("method inverse and ^ - 1d ", vR_b_q(0), vPR_b_tau(0),1e-12);
            check("method inverse and ^ - 2d ", vR_b_q(1), vPR_b_tau(1),1e-12);
            check("method inverse and ^ - 3d ", vR_b_q(2), vPR_b_tau(2),1e-12);
            check("method inverse and ^ - 1e ", vR_b_q(0), vPR_b_Z(0),  1e-12);
            check("method inverse and ^ - 2e ", vR_b_q(1), vPR_b_Z(1),  1e-12);
            check("method inverse and ^ - 3e ", vR_b_q(2), vPR_b_Z(2),  1e-12);
            check("method inverse and ^ - 1f ", vR_b_q(0), vPR_b_S(0),  1e-12);
            check("method inverse and ^ - 2f ", vR_b_q(1), vPR_b_S(1),  1e-12);
            check("method inverse and ^ - 3f ", vR_b_q(2), vPR_b_S(2),  1e-12);

            // direct rotations (no translations)
            vR_x_q    = gq1_xb ^ vR_b_q;
            vR_x_R    = gR1_xb ^ vR_b_R;
            vR_x_M    = gM1_xb ^ vR_b_M;
            vR_x_tau  = gtau1_xb ^ vR_b_tau;
            vR_x_Z    = gZ1_xb ^ vR_b_Z;
            vR_x_S    = gS1_xb ^ vR_b_S;
            vPR_x_q   = gq1_xb.inverse() & vPR_b_q;
            vPR_x_R   = gR1_xb.inverse() & vPR_b_R;
            vPR_x_M   = gM1_xb.inverse() & vPR_b_M;
            vPR_x_tau = gtau1_xb.inverse() & vPR_b_tau;
            vPR_x_Z   = gZ1_xb.inverse() & vPR_b_Z;
            vPR_x_S   = gS1_xb.inverse() & vPR_b_S;

            check("method ^ - 1a             ", v_x(0), vR_x_q(0),   1e-12);
            check("method ^ - 2a             ", v_x(1), vR_x_q(1),   1e-12);
            check("method ^ - 3a             ", v_x(2), vR_x_q(2),   1e-12);
            check("method ^ - 1b             ", v_x(0), vR_x_R(0),   1e-12);
            check("method ^ - 2b             ", v_x(1), vR_x_R(1),   1e-12);
            check("method ^ - 3b             ", v_x(2), vR_x_R(2),   1e-12);
            check("method ^ - 1c             ", v_x(0), vR_x_M(0),   1e-12);
            check("method ^ - 2c             ", v_x(1), vR_x_M(1),   1e-12);
            check("method ^ - 3c             ", v_x(2), vR_x_M(2),   1e-12);
            check("method ^ - 1d             ", v_x(0), vR_x_tau(0),  1e-12);
            check("method ^ - 2d             ", v_x(1), vR_x_tau(1),  1e-12);
            check("method ^ - 3d             ", v_x(2), vR_x_tau(2),  1e-12);
            check("method ^ - 1e             ", v_x(0), vR_x_Z(0),   1e-12);
            check("method ^ - 2e             ", v_x(1), vR_x_Z(1),   1e-12);
            check("method ^ - 3e             ", v_x(2), vR_x_Z(2),   1e-12);
            check("method ^ - 1f             ", v_x(0), vR_x_S(0),   1e-12);
            check("method ^ - 2f             ", v_x(1), vR_x_S(1),   1e-12);
            check("method ^ - 3f             ", v_x(2), vR_x_S(2),   1e-12);
            check("method inverse and & - 1a ", v_x(0), vPR_x_q(0),  1e-12);
            check("method inverse and & - 2a ", v_x(1), vPR_x_q(1),  1e-12);
            check("method inverse and & - 3a ", v_x(2), vPR_x_q(2),  1e-12);
            check("method inverse and & - 1b ", v_x(0), vPR_x_R(0),  1e-12);
            check("method inverse and & - 2b ", v_x(1), vPR_x_R(1),  1e-12);
            check("method inverse and & - 3b ", v_x(2), vPR_x_R(2),  1e-12);
            check("method inverse and & - 1c ", v_x(0), vPR_x_M(0),  1e-12);
            check("method inverse and & - 2c ", v_x(1), vPR_x_M(1),  1e-12);
            check("method inverse and & - 3c ", v_x(2), vPR_x_M(2),  1e-12);
            check("method inverse and & - 1d ", v_x(0), vPR_x_tau(0), 1e-12);
            check("method inverse and & - 2d ", v_x(1), vPR_x_tau(1), 1e-12);
            check("method inverse and & - 3d ", v_x(2), vPR_x_tau(2), 1e-12);
            check("method inverse and & - 1e ", v_x(0), vPR_x_Z(0),  1e-12);
            check("method inverse and & - 2e ", v_x(1), vPR_x_Z(1),  1e-12);
            check("method inverse and & - 3e ", v_x(2), vPR_x_Z(2),  1e-12);
            check("method inverse and & - 1f ", v_x(0), vPR_x_S(0),  1e-12);
            check("method inverse and & - 2f ", v_x(1), vPR_x_S(1),  1e-12);
            check("method inverse and & - 3f ", v_x(2), vPR_x_S(2),  1e-12);

            // body angular and linear velocities
            gqdot1_xb = gq1_xb.xibody2dot(xi_xbb_mrps);
            gRdot1_xb = gR1_xb.xibody2dot(xi_xbb_mrps);
            gMdot1_xb = gM1_xb.xibody2dot(xi_xbb_mrps);
            gZdot1_xb = gZ1_xb.xibody2dot(xi_xbb_mrps);

            xi_xbb_mrps_speu_rodrigues = gq1_xb.dot2xibody(gqdot1_xb);
            xi_xbb_mrps_speu_dcm       = gR1_xb.dot2xibody(gRdot1_xb);
            xi_xbb_mrps_homogeneous    = gM1_xb.dot2xibody(gMdot1_xb);
            xi_xbb_mrps_dual           = gZ1_xb.dot2xibody(gZdot1_xb);

//            cout << gRdot1_xb << endl << endl;
//            cout << gMdot1_xb << endl << endl;
//            cout << xi_xbb_mrps_speu_dcm() << endl << endl;
//            cout << xi_xbb_mrps_homogeneous() << endl << endl;

            check("vi body - 1a             ", xi_xbb_mrps()(0), xi_xbb_mrps_speu_rodrigues()(0), 1e-12);
            check("vi body - 2a             ", xi_xbb_mrps()(1), xi_xbb_mrps_speu_rodrigues()(1), 1e-12);
            check("vi body - 3a             ", xi_xbb_mrps()(2), xi_xbb_mrps_speu_rodrigues()(2), 1e-12);
            check("omega body - 1a          ", xi_xbb_mrps()(3), xi_xbb_mrps_speu_rodrigues()(3), 1e-12);
            check("omega body - 2a          ", xi_xbb_mrps()(4), xi_xbb_mrps_speu_rodrigues()(4), 1e-12);
            check("omega body - 3a          ", xi_xbb_mrps()(5), xi_xbb_mrps_speu_rodrigues()(5), 1e-12);
            check("vi body - 1b             ", xi_xbb_mrps()(0), xi_xbb_mrps_speu_dcm()(0), 1e-12);
            check("vi body - 2b             ", xi_xbb_mrps()(1), xi_xbb_mrps_speu_dcm()(1), 1e-12);
            check("vi body - 3b             ", xi_xbb_mrps()(2), xi_xbb_mrps_speu_dcm()(2), 1e-12);
            check("omega body - 1b          ", xi_xbb_mrps()(3), xi_xbb_mrps_speu_dcm()(3), 1e-12);
            check("omega body - 2b          ", xi_xbb_mrps()(4), xi_xbb_mrps_speu_dcm()(4), 1e-12);
            check("omega body - 3b          ", xi_xbb_mrps()(5), xi_xbb_mrps_speu_dcm()(5), 1e-12);
            check("vi body - 1c             ", xi_xbb_mrps()(0), xi_xbb_mrps_homogeneous()(0), 1e-12);
            check("vi body - 2c             ", xi_xbb_mrps()(1), xi_xbb_mrps_homogeneous()(1), 1e-12);
            check("vi body - 3c             ", xi_xbb_mrps()(2), xi_xbb_mrps_homogeneous()(2), 1e-12);
            check("omega body - 1c          ", xi_xbb_mrps()(3), xi_xbb_mrps_homogeneous()(3), 1e-12);
            check("omega body - 2c          ", xi_xbb_mrps()(4), xi_xbb_mrps_homogeneous()(4), 1e-12);
            check("omega body - 3c          ", xi_xbb_mrps()(5), xi_xbb_mrps_homogeneous()(5), 1e-12);
            check("vi body - 1d             ", xi_xbb_mrps()(0), xi_xbb_mrps_dual()(0), 1e-12);
            check("vi body - 2d             ", xi_xbb_mrps()(1), xi_xbb_mrps_dual()(1), 1e-12);
            check("vi body - 3d             ", xi_xbb_mrps()(2), xi_xbb_mrps_dual()(2), 1e-12);
            check("omega body - 1d          ", xi_xbb_mrps()(3), xi_xbb_mrps_dual()(3), 1e-12);
            check("omega body - 2d          ", xi_xbb_mrps()(4), xi_xbb_mrps_dual()(4), 1e-12);
            check("omega body - 3d          ", xi_xbb_mrps()(5), xi_xbb_mrps_dual()(5), 1e-12);

            // space angular and linear velocities
            xi_xbx_mrps = gM1_xb | xi_xbb_mrps;
            gqdot2_xb   = gq1_xb.xispace2dot(xi_xbx_mrps);
            gRdot2_xb   = gR1_xb.xispace2dot(xi_xbx_mrps);
            gMdot2_xb   = gM1_xb.xispace2dot(xi_xbx_mrps);
            gZdot2_xb   = gZ1_xb.xispace2dot(xi_xbx_mrps);

            xi_xbx_mrps_speu_rodrigues = gq1_xb.dot2xispace(gqdot2_xb);
            xi_xbx_mrps_speu_dcm       = gR1_xb.dot2xispace(gRdot2_xb);
            xi_xbx_mrps_homogeneous    = gM1_xb.dot2xispace(gMdot2_xb);
            xi_xbx_mrps_dual           = gZ1_xb.dot2xispace(gZdot2_xb);

//            cout << gRdot2_xb << endl << endl;
//            cout << gMdot2_xb << endl << endl;

            check("omega body-space speu_rodrigues0", gqdot1_xb(0), gqdot2_xb(0), 1e-12);
            check("omega body-space speu_rodrigues1", gqdot1_xb(1), gqdot2_xb(1), 1e-12);
            check("omega body-space speu_rodrigues2", gqdot1_xb(2), gqdot2_xb(2), 1e-12);
            check("omega body-space speu_rodrigues3", gqdot1_xb(3), gqdot2_xb(3), 1e-12);
            check("omega body-space speu_rodrigues4", gqdot1_xb(4), gqdot2_xb(4), 1e-12);
            check("omega body-space speu_rodrigues5", gqdot1_xb(5), gqdot2_xb(5), 1e-12);
            check("omega body-space speu_rodrigues6", gqdot1_xb(6), gqdot2_xb(6), 1e-12);
            check("omega body-space speu_dcm00     ", gRdot1_xb(0,0), gRdot2_xb(0,0), 1e-12);
            check("omega body-space speu_dcm01     ", gRdot1_xb(0,1), gRdot2_xb(0,1), 1e-12);
            check("omega body-space speu_dcm02     ", gRdot1_xb(0,2), gRdot2_xb(0,2), 1e-12);
            check("omega body-space speu_dcm03     ", gRdot1_xb(0,3), gRdot2_xb(0,3), 1e-12);
            check("omega body-space speu_dcm10     ", gRdot1_xb(1,0), gRdot2_xb(1,0), 1e-12);
            check("omega body-space speu_dcm11     ", gRdot1_xb(1,1), gRdot2_xb(1,1), 1e-12);
            check("omega body-space speu_dcm12     ", gRdot1_xb(1,2), gRdot2_xb(1,2), 1e-12);
            check("omega body-space speu_dcm13     ", gRdot1_xb(1,3), gRdot2_xb(1,3), 1e-12);
            check("omega body-space speu_dcm20     ", gRdot1_xb(2,0), gRdot2_xb(2,0), 1e-12);
            check("omega body-space speu_dcm21     ", gRdot1_xb(2,1), gRdot2_xb(2,1), 1e-12);
            check("omega body-space speu_dcm22     ", gRdot1_xb(2,2), gRdot2_xb(2,2), 1e-12);
            check("omega body-space speu_dcm23     ", gRdot1_xb(2,3), gRdot2_xb(2,3), 1e-12);
            check("omega body-space homogeneous00  ", gMdot1_xb(0,0), gMdot2_xb(0,0), 1e-12);
            check("omega body-space homogeneous01  ", gMdot1_xb(0,1), gMdot2_xb(0,1), 1e-12);
            check("omega body-space homogeneous02  ", gMdot1_xb(0,2), gMdot2_xb(0,2), 1e-12);
            check("omega body-space homogeneous03  ", gMdot1_xb(0,3), gMdot2_xb(0,3), 1e-12);
            check("omega body-space homogeneous10  ", gMdot1_xb(1,0), gMdot2_xb(1,0), 1e-12);
            check("omega body-space homogeneous11  ", gMdot1_xb(1,1), gMdot2_xb(1,1), 1e-12);
            check("omega body-space homogeneous12  ", gMdot1_xb(1,2), gMdot2_xb(1,2), 1e-12);
            check("omega body-space homogeneous13  ", gMdot1_xb(1,3), gMdot2_xb(1,3), 1e-12);
            check("omega body-space homogeneous20  ", gMdot1_xb(2,0), gMdot2_xb(2,0), 1e-12);
            check("omega body-space homogeneous21  ", gMdot1_xb(2,1), gMdot2_xb(2,1), 1e-12);
            check("omega body-space homogeneous22  ", gMdot1_xb(2,2), gMdot2_xb(2,2), 1e-12);
            check("omega body-space homogeneous23  ", gMdot1_xb(2,3), gMdot2_xb(2,3), 1e-12);
            check("omega body-space homogeneous30  ", gMdot1_xb(3,0), gMdot2_xb(3,0), 1e-12);
            check("omega body-space homogeneous31  ", gMdot1_xb(3,1), gMdot2_xb(3,1), 1e-12);
            check("omega body-space homogeneous32  ", gMdot1_xb(3,2), gMdot2_xb(3,2), 1e-12);
            check("omega body-space homogeneous33  ", gMdot1_xb(3,3), gMdot2_xb(3,3), 1e-12);
            check("omega body-space dual r0        ", gZdot1_xb.get_qr()(0), gZdot2_xb.get_qr()(0), 1e-12);
            check("omega body-space dual r1        ", gZdot1_xb.get_qr()(1), gZdot2_xb.get_qr()(1), 1e-12);
            check("omega body-space dual r2        ", gZdot1_xb.get_qr()(2), gZdot2_xb.get_qr()(2), 1e-12);
            check("omega body-space dual r3        ", gZdot1_xb.get_qr()(3), gZdot2_xb.get_qr()(3), 1e-12);
            check("omega body-space dual d0        ", gZdot1_xb.get_qd()(0), gZdot2_xb.get_qd()(0), 1e-12);
            check("omega body-space dual d1        ", gZdot1_xb.get_qd()(1), gZdot2_xb.get_qd()(1), 1e-12);
            check("omega body-space dual d2        ", gZdot1_xb.get_qd()(2), gZdot2_xb.get_qd()(2), 1e-12);
            check("omega body-space dual d3        ", gZdot1_xb.get_qd()(3), gZdot2_xb.get_qd()(3), 1e-12);

            check("vi space - 1a            ", xi_xbx_mrps()(0), xi_xbx_mrps_speu_rodrigues()(0), 1e-12);
            check("vi space - 2a            ", xi_xbx_mrps()(1), xi_xbx_mrps_speu_rodrigues()(1), 1e-12);
            check("vi space - 3a            ", xi_xbx_mrps()(2), xi_xbx_mrps_speu_rodrigues()(2), 1e-12);
            check("omega space - 1a         ", xi_xbx_mrps()(3), xi_xbx_mrps_speu_rodrigues()(3), 1e-12);
            check("omega space - 2a         ", xi_xbx_mrps()(4), xi_xbx_mrps_speu_rodrigues()(4), 1e-12);
            check("omega space - 3a         ", xi_xbx_mrps()(5), xi_xbx_mrps_speu_rodrigues()(5), 1e-12);
            check("vi space - 1b            ", xi_xbx_mrps()(0), xi_xbx_mrps_speu_dcm()(0), 1e-12);
            check("vi space - 2b            ", xi_xbx_mrps()(1), xi_xbx_mrps_speu_dcm()(1), 1e-12);
            check("vi space - 3b            ", xi_xbx_mrps()(2), xi_xbx_mrps_speu_dcm()(2), 1e-12);
            check("omega space - 1b         ", xi_xbx_mrps()(3), xi_xbx_mrps_speu_dcm()(3), 1e-12);
            check("omega space - 2b         ", xi_xbx_mrps()(4), xi_xbx_mrps_speu_dcm()(4), 1e-12);
            check("omega space - 3b         ", xi_xbx_mrps()(5), xi_xbx_mrps_speu_dcm()(5), 1e-12);
            check("vi space - 1c            ", xi_xbx_mrps()(0), xi_xbx_mrps_homogeneous()(0), 1e-12);
            check("vi space - 2c            ", xi_xbx_mrps()(1), xi_xbx_mrps_homogeneous()(1), 1e-12);
            check("vi space - 3c            ", xi_xbx_mrps()(2), xi_xbx_mrps_homogeneous()(2), 1e-12);
            check("omega space - 1c         ", xi_xbx_mrps()(3), xi_xbx_mrps_homogeneous()(3), 1e-12);
            check("omega space - 2c         ", xi_xbx_mrps()(4), xi_xbx_mrps_homogeneous()(4), 1e-12);
            check("omega space - 3c         ", xi_xbx_mrps()(5), xi_xbx_mrps_homogeneous()(5), 1e-12);
            check("vi space - 1d            ", xi_xbx_mrps()(0), xi_xbx_mrps_dual()(0), 1e-12);
            check("vi space - 2d            ", xi_xbx_mrps()(1), xi_xbx_mrps_dual()(1), 1e-12);
            check("vi space - 3d            ", xi_xbx_mrps()(2), xi_xbx_mrps_dual()(2), 1e-12);
            check("omega space - 1d         ", xi_xbx_mrps()(3), xi_xbx_mrps_dual()(3), 1e-12);
            check("omega space - 2d         ", xi_xbx_mrps()(4), xi_xbx_mrps_dual()(4), 1e-12);
            check("omega space - 3d         ", xi_xbx_mrps()(5), xi_xbx_mrps_dual()(5), 1e-12);

            // direct transformation concatenation
            gq_xc = gq1_xb * gq_bc;
            gR_xc = gR1_xb * gR_bc;
            gM_xc = gM1_xb * gM_bc;
            gtau_xc  = gtau1_xb * gtau_bc;
            gZ_xc  = gZ1_xb * gZ_bc;
            gS_xc  = gS1_xb * gS_bc;

            check("concatenation rot vector 1a   :", gq_xc.get_rotv()()(0), gR_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2a   :", gq_xc.get_rotv()()(1), gR_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3a   :", gq_xc.get_rotv()()(2), gR_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1b   :", gq_xc.get_rotv()()(0), gM_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2b   :", gq_xc.get_rotv()()(1), gM_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3b   :", gq_xc.get_rotv()()(2), gM_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1c   :", gq_xc.get_rotv()()(0), gtau_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2c   :", gq_xc.get_rotv()()(1), gtau_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3c   :", gq_xc.get_rotv()()(2), gtau_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1d   :", gq_xc.get_rotv()()(0), gZ_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2d   :", gq_xc.get_rotv()()(1), gZ_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3d   :", gq_xc.get_rotv()()(2), gZ_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1e   :", gq_xc.get_rotv()()(0), gS_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2e   :", gq_xc.get_rotv()()(1), gS_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3e   :", gq_xc.get_rotv()()(2), gS_xc.get_rotv()()(2), 1e-12);

            check("concatenation translation 1a   :", gq_xc.get_T()(0), gR_xc.get_T()(0), 1e-12);
            check("concatenation translation 2a   :", gq_xc.get_T()(1), gR_xc.get_T()(1), 1e-12);
            check("concatenation translation 3a   :", gq_xc.get_T()(2), gR_xc.get_T()(2), 1e-12);
            check("concatenation translation 1b   :", gq_xc.get_T()(0), gM_xc.get_T()(0), 1e-12);
            check("concatenation translation 2b   :", gq_xc.get_T()(1), gM_xc.get_T()(1), 1e-12);
            check("concatenation translation 3b   :", gq_xc.get_T()(2), gM_xc.get_T()(2), 1e-12);
            check("concatenation translation 1c   :", gq_xc.get_T()(0), gtau_xc.get_T()(0), 1e-12);
            check("concatenation translation 2c   :", gq_xc.get_T()(1), gtau_xc.get_T()(1), 1e-12);
            check("concatenation translation 3c   :", gq_xc.get_T()(2), gtau_xc.get_T()(2), 1e-12);
            check("concatenation translation 1d   :", gq_xc.get_T()(0), gZ_xc.get_T()(0), 1e-12);
            check("concatenation translation 2d   :", gq_xc.get_T()(1), gZ_xc.get_T()(1), 1e-12);
            check("concatenation translation 3d   :", gq_xc.get_T()(2), gZ_xc.get_T()(2), 1e-12);
            check("concatenation translation 1e   :", gq_xc.get_T()(0), gS_xc.get_T()(0), 1e-12);
            check("concatenation translation 2e   :", gq_xc.get_T()(1), gS_xc.get_T()(1), 1e-12);
            check("concatenation translation 3e   :", gq_xc.get_T()(2), gS_xc.get_T()(2), 1e-12);

            // inverse transformation concatenation
            gqX_xc = gq1_xb / gq_bc;
            gRX_xc = gR1_xb / gR_bc;
            gMX_xc = gM1_xb / gM_bc;
            gtauX_xc  = gtau1_xb / gtau_bc;
            gZX_xc  = gZ1_xb / gZ_bc;
            gSX_xc  = gS1_xb / gS_bc;

            check("concatenation rot vector 1a   :", gqX_xc.get_rotv()()(0), gRX_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2a   :", gqX_xc.get_rotv()()(1), gRX_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3a   :", gqX_xc.get_rotv()()(2), gRX_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1b   :", gqX_xc.get_rotv()()(0), gMX_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2b   :", gqX_xc.get_rotv()()(1), gMX_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3b   :", gqX_xc.get_rotv()()(2), gMX_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1c   :", gqX_xc.get_rotv()()(0), gtauX_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2c   :", gqX_xc.get_rotv()()(1), gtauX_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3c   :", gqX_xc.get_rotv()()(2), gtauX_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1d   :", gqX_xc.get_rotv()()(0), gZX_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2d   :", gqX_xc.get_rotv()()(1), gZX_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3d   :", gqX_xc.get_rotv()()(2), gZX_xc.get_rotv()()(2), 1e-12);
            check("concatenation rot vector 1e   :", gqX_xc.get_rotv()()(0), gSX_xc.get_rotv()()(0), 1e-12);
            check("concatenation rot vector 2e   :", gqX_xc.get_rotv()()(1), gSX_xc.get_rotv()()(1), 1e-12);
            check("concatenation rot vector 3e   :", gqX_xc.get_rotv()()(2), gSX_xc.get_rotv()()(2), 1e-12);

            check("concatenation translation 1a   :", gqX_xc.get_T()(0), gRX_xc.get_T()(0), 1e-12);
            check("concatenation translation 2a   :", gqX_xc.get_T()(1), gRX_xc.get_T()(1), 1e-12);
            check("concatenation translation 3a   :", gqX_xc.get_T()(2), gRX_xc.get_T()(2), 1e-12);
            check("concatenation translation 1b   :", gqX_xc.get_T()(0), gMX_xc.get_T()(0), 1e-12);
            check("concatenation translation 2b   :", gqX_xc.get_T()(1), gMX_xc.get_T()(1), 1e-12);
            check("concatenation translation 3b   :", gqX_xc.get_T()(2), gMX_xc.get_T()(2), 1e-12);
            check("concatenation translation 1c   :", gqX_xc.get_T()(0), gtauX_xc.get_T()(0), 1e-12);
            check("concatenation translation 2c   :", gqX_xc.get_T()(1), gtauX_xc.get_T()(1), 1e-12);
            check("concatenation translation 3c   :", gqX_xc.get_T()(2), gtauX_xc.get_T()(2), 1e-12);
            check("concatenation translation 1d   :", gqX_xc.get_T()(0), gZX_xc.get_T()(0), 1e-12);
            check("concatenation translation 2d   :", gqX_xc.get_T()(1), gZX_xc.get_T()(1), 1e-12);
            check("concatenation translation 3d   :", gqX_xc.get_T()(2), gZX_xc.get_T()(2), 1e-12);
            check("concatenation translation 1e   :", gqX_xc.get_T()(0), gSX_xc.get_T()(0), 1e-12);
            check("concatenation translation 2e   :", gqX_xc.get_T()(1), gSX_xc.get_T()(1), 1e-12);
            check("concatenation translation 3e   :", gqX_xc.get_T()(2), gSX_xc.get_T()(2), 1e-12);

            // inverse
            gq_bx    = gq1_xb.inverse();
            gR_bx    = gR1_xb.inverse();
            gM_bx    = gM1_xb.inverse();
            gtau_bx  = gtau1_xb.inverse();
            gZ_bx    = gZ1_xb.inverse();
            gS_bx    = gS1_xb.inverse();

            check("inverse rot vector 1a   :", gq_bx.get_rotv()()(0), gR_bx.get_rotv()()(0), 1e-12);
            check("inverse rot vector 2a   :", gq_bx.get_rotv()()(1), gR_bx.get_rotv()()(1), 1e-12);
            check("inverse rot vector 3a   :", gq_bx.get_rotv()()(2), gR_bx.get_rotv()()(2), 1e-12);
            check("inverse rot vector 1b   :", gq_bx.get_rotv()()(0), gM_bx.get_rotv()()(0), 1e-12);
            check("inverse rot vector 2b   :", gq_bx.get_rotv()()(1), gM_bx.get_rotv()()(1), 1e-12);
            check("inverse rot vector 3b   :", gq_bx.get_rotv()()(2), gM_bx.get_rotv()()(2), 1e-12);
            check("inverse rot vector 1c   :", gq_bx.get_rotv()()(0), gtau_bx.get_rotv()()(0), 1e-12);
            check("inverse rot vector 2c   :", gq_bx.get_rotv()()(1), gtau_bx.get_rotv()()(1), 1e-12);
            check("inverse rot vector 3c   :", gq_bx.get_rotv()()(2), gtau_bx.get_rotv()()(2), 1e-12);
            check("inverse rot vector 1d   :", gq_bx.get_rotv()()(0), gZ_bx.get_rotv()()(0), 1e-12);
            check("inverse rot vector 2d   :", gq_bx.get_rotv()()(1), gZ_bx.get_rotv()()(1), 1e-12);
            check("inverse rot vector 3d   :", gq_bx.get_rotv()()(2), gZ_bx.get_rotv()()(2), 1e-12);
            check("inverse rot vector 1e   :", gq_bx.get_rotv()()(0), gS_bx.get_rotv()()(0), 1e-12);
            check("inverse rot vector 2e   :", gq_bx.get_rotv()()(1), gS_bx.get_rotv()()(1), 1e-12);
            check("inverse rot vector 3e   :", gq_bx.get_rotv()()(2), gS_bx.get_rotv()()(2), 1e-12);

            check("inverse translation 1a   :", gq_bx.get_T()(0), gR_bx.get_T()(0), 1e-12);
            check("inverse translation 2a   :", gq_bx.get_T()(1), gR_bx.get_T()(1), 1e-12);
            check("inverse translation 3a   :", gq_bx.get_T()(2), gR_bx.get_T()(2), 1e-12);
            check("inverse translation 1b   :", gq_bx.get_T()(0), gM_bx.get_T()(0), 1e-12);
            check("inverse translation 2b   :", gq_bx.get_T()(1), gM_bx.get_T()(1), 1e-12);
            check("inverse translation 3b   :", gq_bx.get_T()(2), gM_bx.get_T()(2), 1e-12);
            check("inverse translation 1c   :", gq_bx.get_T()(0), gtau_bx.get_T()(0), 1e-12);
            check("inverse translation 2c   :", gq_bx.get_T()(1), gtau_bx.get_T()(1), 1e-12);
            check("inverse translation 3c   :", gq_bx.get_T()(2), gtau_bx.get_T()(2), 1e-12);
            check("inverse translation 1d   :", gq_bx.get_T()(0), gZ_bx.get_T()(0), 1e-12);
            check("inverse translation 2d   :", gq_bx.get_T()(1), gZ_bx.get_T()(1), 1e-12);
            check("inverse translation 3d   :", gq_bx.get_T()(2), gZ_bx.get_T()(2), 1e-12);
            check("inverse translation 1e   :", gq_bx.get_T()(0), gS_bx.get_T()(0), 1e-12);
            check("inverse translation 2e   :", gq_bx.get_T()(1), gS_bx.get_T()(1), 1e-12);
            check("inverse translation 3e   :", gq_bx.get_T()(2), gS_bx.get_T()(2), 1e-12);

            // log and exp
            ang::trfv Ztau2 = gtau1_xb.exp_map_speu_rodrigues().log_map_trfv();
            ang::trfv Ztau3 = gtau1_xb.exp_map_speu_dcm().log_map_trfv();
            ang::trfv Ztau4 = gtau1_xb.exp_map_homogeneous().log_map_trfv();

            check("exp is the opposite of log 1 bb ", gtau1_xb.get_rotv()()(0), Ztau2.get_rotv()()(0), 1e-12);
            check("exp is the opposite of log 2 bb ", gtau1_xb.get_rotv()()(1), Ztau2.get_rotv()()(1), 1e-12);
            check("exp is the opposite of log 3 bb ", gtau1_xb.get_rotv()()(2), Ztau2.get_rotv()()(2), 1e-12);
            check("exp is the opposite of log 1 cc ", gtau1_xb.get_rotv()()(0), Ztau3.get_rotv()()(0), 1e-12);
            check("exp is the opposite of log 2 cc ", gtau1_xb.get_rotv()()(1), Ztau3.get_rotv()()(1), 1e-12);
            check("exp is the opposite of log 3 cc ", gtau1_xb.get_rotv()()(2), Ztau3.get_rotv()()(2), 1e-12);
            check("exp is the opposite of log 1 dd ", gtau1_xb.get_rotv()()(0), Ztau4.get_rotv()()(0), 1e-12);
            check("exp is the opposite of log 2 dd ", gtau1_xb.get_rotv()()(1), Ztau4.get_rotv()()(1), 1e-12);
            check("exp is the opposite of log 3 dd ", gtau1_xb.get_rotv()()(2), Ztau4.get_rotv()()(2), 1e-12);

            check("exp is the opposite of log 1 bb ", gtau1_xb.get_T()(0), Ztau2.get_T()(0), 1e-12);
            check("exp is the opposite of log 2 bb ", gtau1_xb.get_T()(1), Ztau2.get_T()(1), 1e-12);
            check("exp is the opposite of log 3 bb ", gtau1_xb.get_T()(2), Ztau2.get_T()(2), 1e-12);
            check("exp is the opposite of log 1 cc ", gtau1_xb.get_T()(0), Ztau3.get_T()(0), 1e-12);
            check("exp is the opposite of log 2 cc ", gtau1_xb.get_T()(1), Ztau3.get_T()(1), 1e-12);
            check("exp is the opposite of log 3 cc ", gtau1_xb.get_T()(2), Ztau3.get_T()(2), 1e-12);
            check("exp is the opposite of log 1 dd ", gtau1_xb.get_T()(0), Ztau4.get_T()(0), 1e-12);
            check("exp is the opposite of log 2 dd ", gtau1_xb.get_T()(1), Ztau4.get_T()(1), 1e-12);
            check("exp is the opposite of log 3 dd ", gtau1_xb.get_T()(2), Ztau4.get_T()(2), 1e-12);

        }
    }
} // closes test_se3

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_exp_log_maps() {

    // the purpose of this test is to verify that the exponential and logarithmic functions relating
    // the transform vector with other representations with the common
    // arithmetic rules for exponential and logarithmic functions.

    // different rotation directions (not unitary)
    vector<Eigen::Vector3d> Vdir_full(26);
    Vdir_full[0]  << +1.,  0.,  0.;
    Vdir_full[1]  << -1.,  0.,  0.;
    Vdir_full[2]  <<  0., +1.,  0.;
    Vdir_full[3]  <<  0., -1.,  0.;
    Vdir_full[4]  <<  0.,  0., +1.;
    Vdir_full[5]  <<  0.,  0., -1.;
    Vdir_full[6]  << +1., +2.,  0.;
    Vdir_full[7]  << +1., -2.,  0.;
    Vdir_full[8]  << -1., +2.,  0.;
    Vdir_full[9]  << -1., -2.,  0.;
    Vdir_full[10] << +1.,  0., +2.;
    Vdir_full[11] << +1.,  0., -2.;
    Vdir_full[12] << -1.,  0., +2.;
    Vdir_full[13] << -1.,  0., -2.;
    Vdir_full[14] <<  0., +1., +2.;
    Vdir_full[15] <<  0., +1., -2.;
    Vdir_full[16] <<  0., -1., +2.;
    Vdir_full[17] <<  0., -1., -2.;
    Vdir_full[18] << +1., +2., +3.;
    Vdir_full[19] << +1., +2., -3.;
    Vdir_full[20] << +1., -2., +3.;
    Vdir_full[21] << +1., -2., -3.;
    Vdir_full[22] << -1., +2., +3.;
    Vdir_full[23] << -1., +2., -3.;
    Vdir_full[24] << -1., -2., +3.;
    Vdir_full[25] << -1., -2., -3.;

    // different rotation angles (less than 180 [deg] as otherwise shortest path is chosen)
    vector<double> Vangle_rad(11);
    Vangle_rad[0]  = 0.7 * math::constant::PI(); // more than  90 deg
    Vangle_rad[1]  = 0.2 * math::constant::PI(); // less than 90 deg
    Vangle_rad[2]  = 1e-3;
    Vangle_rad[3]  = 1e-5;
    Vangle_rad[4]  = 1e-7;
    Vangle_rad[5]  = 1e-9;
    Vangle_rad[6]  = 1e-11;
    Vangle_rad[7]  = 1e-13;
    Vangle_rad[8]  = 1e-15;
    Vangle_rad[9]  = 1e-16;
    Vangle_rad[10] = 0.;

    Eigen::Vector3d dir_unit;
    for (unsigned int i = 0; i != Vdir_full.size(); ++i) {
        dir_unit = Vdir_full[i] / Vdir_full[i].norm();
        for (unsigned int j = 0; j != Vangle_rad.size(); ++j) {

            cout << endl << "Iteration i = " << i << "  j = " << j << endl;

            ang::rotv Orotv_ba(dir_unit * Vangle_rad[j]);

            Eigen::Vector3d T_bab(0.41, -0.37, 0.11);
            ang::trfv Otrfv_ba(Orotv_ba, T_bab);
            ang::rotv Orotv_cb(-0.02, 0.12, 0.07);
            Eigen::Vector3d T_cbc(0.03, 0.05, -0.08);
            ang::trfv Otrfv_cb(Orotv_cb, T_cbc);
            ang::trfv Otrfv_ca = Otrfv_cb * Otrfv_ba;
            ang::trfv Otrfv_3ba = Otrfv_ba * Otrfv_ba * Otrfv_ba;

            // check that exponential equals the inverse of logarithm --> à = log(exp(a))
            ang::speu_rodrigues Osq_ba = Otrfv_ba.exp_map_speu_rodrigues();
            ang::speu_dcm Osdcm_ba     = Otrfv_ba.exp_map_speu_dcm();
            ang::homogeneous OM_ba     = Otrfv_ba.exp_map_homogeneous();
            ang::screw OS_ba           = Otrfv_ba.explog_map_screw();
            ang::dual OZ_ba            = Otrfv_ba.exp_map_dual();
            ang::trfv Otrfv1_ba        = Osq_ba.log_map_trfv();
            ang::trfv Otrfv2_ba        = Osdcm_ba.log_map_trfv();
            ang::trfv Otrfv3_ba        = OM_ba.log_map_trfv();
            ang::trfv Otrfv4_ba        = OS_ba.explog_map_trfv();
            ang::trfv Otrfv5_ba        = OZ_ba.log_map_trfv();
            check("exp-log quat 0       ", Otrfv_ba.get_rotv()()(0), Otrfv1_ba.get_rotv()()(0), 1e-12);
            check("exp-log quat 1       ", Otrfv_ba.get_rotv()()(1), Otrfv1_ba.get_rotv()()(1), 1e-12);
            check("exp-log quat 2       ", Otrfv_ba.get_rotv()()(2), Otrfv1_ba.get_rotv()()(2), 1e-12);
            check("exp-log quat 0       ", Otrfv_ba.get_T()(0), Otrfv1_ba.get_T()(0), 1e-12);
            check("exp-log quat 1       ", Otrfv_ba.get_T()(1), Otrfv1_ba.get_T()(1), 1e-12);
            check("exp-log quat 2       ", Otrfv_ba.get_T()(2), Otrfv1_ba.get_T()(2), 1e-12);
            check("exp-log dcm 0        ", Otrfv_ba.get_rotv()()(0), Otrfv2_ba.get_rotv()()(0), 1e-12);
            check("exp-log dcm 1        ", Otrfv_ba.get_rotv()()(1), Otrfv2_ba.get_rotv()()(1), 1e-12);
            check("exp-log dcm 2        ", Otrfv_ba.get_rotv()()(2), Otrfv2_ba.get_rotv()()(2), 1e-12);
            check("exp-log dcm 0        ", Otrfv_ba.get_T()(0), Otrfv2_ba.get_T()(0), 1e-12);
            check("exp-log dcm 1        ", Otrfv_ba.get_T()(1), Otrfv2_ba.get_T()(1), 1e-12);
            check("exp-log dcm 2        ", Otrfv_ba.get_T()(2), Otrfv2_ba.get_T()(2), 1e-12);
            check("exp-log homo 0       ", Otrfv_ba.get_rotv()()(0), Otrfv3_ba.get_rotv()()(0), 1e-12);
            check("exp-log homo 1       ", Otrfv_ba.get_rotv()()(1), Otrfv3_ba.get_rotv()()(1), 1e-12);
            check("exp-log homo 2       ", Otrfv_ba.get_rotv()()(2), Otrfv3_ba.get_rotv()()(2), 1e-12);
            check("exp-log homo 0       ", Otrfv_ba.get_T()(0), Otrfv3_ba.get_T()(0), 1e-12);
            check("exp-log homo 1       ", Otrfv_ba.get_T()(1), Otrfv3_ba.get_T()(1), 1e-12);
            check("exp-log homo 2       ", Otrfv_ba.get_T()(2), Otrfv3_ba.get_T()(2), 1e-12);
            check("exp-log screw 0      ", Otrfv_ba.get_rotv()()(0), Otrfv4_ba.get_rotv()()(0), 1e-12);
            check("exp-log screw 1      ", Otrfv_ba.get_rotv()()(1), Otrfv4_ba.get_rotv()()(1), 1e-12);
            check("exp-log screw 2      ", Otrfv_ba.get_rotv()()(2), Otrfv4_ba.get_rotv()()(2), 1e-12);
            check("exp-log screw 0      ", Otrfv_ba.get_T()(0), Otrfv4_ba.get_T()(0), 1e-12);
            check("exp-log screw 1      ", Otrfv_ba.get_T()(1), Otrfv4_ba.get_T()(1), 1e-12);
            check("exp-log screw 2      ", Otrfv_ba.get_T()(2), Otrfv4_ba.get_T()(2), 1e-12);
            check("exp-log dual 0       ", Otrfv_ba.get_rotv()()(0), Otrfv5_ba.get_rotv()()(0), 1e-12);
            check("exp-log dual 1       ", Otrfv_ba.get_rotv()()(1), Otrfv5_ba.get_rotv()()(1), 1e-12);
            check("exp-log dual 2       ", Otrfv_ba.get_rotv()()(2), Otrfv5_ba.get_rotv()()(2), 1e-12);
            check("exp-log dual 0       ", Otrfv_ba.get_T()(0), Otrfv5_ba.get_T()(0), 1e-12);
            check("exp-log dual 1       ", Otrfv_ba.get_T()(1), Otrfv5_ba.get_T()(1), 1e-12);
            check("exp-log dual 2       ", Otrfv_ba.get_T()(2), Otrfv5_ba.get_T()(2), 1e-12);
            check("exp-log SE3 0        ", Otrfv_ba.get_rotv()()(0), Otrfv4_ba()(3), 1e-12);
            check("exp-log SE3 1        ", Otrfv_ba.get_rotv()()(1), Otrfv4_ba()(4), 1e-12);
            check("exp-log SE3 2        ", Otrfv_ba.get_rotv()()(2), Otrfv4_ba()(5), 1e-12);
            check("exp-log SE3 0        ", Otrfv_ba.get_s()(0), Otrfv4_ba()(0), 1e-12);
            check("exp-log SE3 1        ", Otrfv_ba.get_s()(1), Otrfv4_ba()(1), 1e-12);
            check("exp-log SE3 2        ", Otrfv_ba.get_s()(2), Otrfv4_ba()(2), 1e-12);
        }
    }
} // closes test_exp_log_maps

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_power() {

    // the purpose of this test is to verify that the power function of a transformation

    // different rotation directions (not unitary)
    vector<Eigen::Vector3d> Vdir_full(26);
    Vdir_full[0]  << +1.,  0.,  0.;
    Vdir_full[1]  << -1.,  0.,  0.;
    Vdir_full[2]  <<  0., +1.,  0.;
    Vdir_full[3]  <<  0., -1.,  0.;
    Vdir_full[4]  <<  0.,  0., +1.;
    Vdir_full[5]  <<  0.,  0., -1.;
    Vdir_full[6]  << +1., +2.,  0.;
    Vdir_full[7]  << +1., -2.,  0.;
    Vdir_full[8]  << -1., +2.,  0.;
    Vdir_full[9]  << -1., -2.,  0.;
    Vdir_full[10] << +1.,  0., +2.;
    Vdir_full[11] << +1.,  0., -2.;
    Vdir_full[12] << -1.,  0., +2.;
    Vdir_full[13] << -1.,  0., -2.;
    Vdir_full[14] <<  0., +1., +2.;
    Vdir_full[15] <<  0., +1., -2.;
    Vdir_full[16] <<  0., -1., +2.;
    Vdir_full[17] <<  0., -1., -2.;
    Vdir_full[18] << +1., +2., +3.;
    Vdir_full[19] << +1., +2., -3.;
    Vdir_full[20] << +1., -2., +3.;
    Vdir_full[21] << +1., -2., -3.;
    Vdir_full[22] << -1., +2., +3.;
    Vdir_full[23] << -1., +2., -3.;
    Vdir_full[24] << -1., -2., +3.;
    Vdir_full[25] << -1., -2., -3.;

    // different rotation angles (less than 180 [deg] as otherwise shortest path is chosen)
    vector<double> Vangle_rad(11);
    Vangle_rad[0]  = 0.7 * math::constant::PI(); // more than  90 deg
    Vangle_rad[1]  = 0.2 * math::constant::PI(); // less than 90 deg
    Vangle_rad[2]  = 1e-3;
    Vangle_rad[3]  = 1e-5;
    Vangle_rad[4]  = 1e-7;
    Vangle_rad[5]  = 1e-9;
    Vangle_rad[6]  = 1e-11;
    Vangle_rad[7]  = 1e-13;
    Vangle_rad[8]  = 1e-15;
    Vangle_rad[9]  = 1e-16;
    Vangle_rad[10] = 0.;

    Eigen::Vector3d dir_unit;
    for (unsigned int i = 0; i != Vdir_full.size(); ++i) {
        dir_unit = Vdir_full[i] / Vdir_full[i].norm();
        for (unsigned int j = 0; j != Vangle_rad.size(); ++j) {

            cout << endl << "Iteration i = " << i << "  j = " << j << endl;

            ang::rotv rv(dir_unit * Vangle_rad[j]);

            Eigen::Vector3d T(3., -4., -1.);

            ang::rodrigues q(rv);
            ang::dcm R(rv);
            ang::speu_rodrigues gq(q, T);
            ang::speu_dcm gR(R, T);
            ang::homogeneous M(R, T);
            ang::trfv tau(rv, T);
            ang::dual Z(rv, T);
            ang::screw S(rv, T);

            ang::speu_rodrigues gq2 = gq.pow(0.3);
            ang::speu_dcm gR2 = gR.pow(0.3);
            ang::homogeneous M2 = M.pow(0.3);
            ang::trfv tau2 = tau.pow(0.3);
            ang::dual Z2 = Z.pow(0.3);
            ang::screw S2 = S.pow(0.3);

            std::cout << S.get_rotv()() << std::endl;
            std::cout << tau2.get_rotv()() << std::endl;
            std::cout << S2.get_rotv()() << std::endl;
            std::cout << S.get_T() << std::endl;
            std::cout << tau2.get_T() << std::endl;
            std::cout << S2.get_T() << std::endl;

            check("exp1 rotv // gq - 1     ", S2.get_rotv()()(0), gq2.get_rotv()()(0), 1e-12);
            check("exp1 rotv // gq - 2     ", S2.get_rotv()()(1), gq2.get_rotv()()(1), 1e-12);
            check("exp1 rotv // gq - 3     ", S2.get_rotv()()(2), gq2.get_rotv()()(2), 1e-12);
            check("exp1 rotv // gR - 1     ", S2.get_rotv()()(0), gR2.get_rotv()()(0), 1e-12);
            check("exp1 rotv // gR - 2     ", S2.get_rotv()()(1), gR2.get_rotv()()(1), 1e-12);
            check("exp1 rotv // gR - 3     ", S2.get_rotv()()(2), gR2.get_rotv()()(2), 1e-12);
            check("exp1 rotv // homo - 1   ", S2.get_rotv()()(0), M2.get_rotv()()(0), 1e-12);
            check("exp1 rotv // homo - 2   ", S2.get_rotv()()(1), M2.get_rotv()()(1), 1e-12);
            check("exp1 rotv // homo - 3   ", S2.get_rotv()()(2), M2.get_rotv()()(2), 1e-12);
            check("exp1 rotv // trfv - 1   ", S2.get_rotv()()(0), tau2.get_rotv()()(0), 1e-12);
            check("exp1 rotv // trfv - 2   ", S2.get_rotv()()(1), tau2.get_rotv()()(1), 1e-12);
            check("exp1 rotv // trfv - 3   ", S2.get_rotv()()(2), tau2.get_rotv()()(2), 1e-12);
            check("exp1 rotv // dual - 1   ", S2.get_rotv()()(0), Z2.get_rotv()()(0), 1e-12);
            check("exp1 rotv // dual - 2   ", S2.get_rotv()()(1), Z2.get_rotv()()(1), 1e-12);
            check("exp1 rotv // dual - 3   ", S2.get_rotv()()(2), Z2.get_rotv()()(2), 1e-12);

            check("exp1 trans // gq - 1     ", S2.get_T()(0), gq2.get_T()(0), 1e-12);
            check("exp1 trans // gq - 2     ", S2.get_T()(1), gq2.get_T()(1), 1e-12);
            check("exp1 trans // gq - 3     ", S2.get_T()(2), gq2.get_T()(2), 1e-12);
            check("exp1 trans // gR - 1     ", S2.get_T()(0), gR2.get_T()(0), 1e-12);
            check("exp1 trans // gR - 2     ", S2.get_T()(1), gR2.get_T()(1), 1e-12);
            check("exp1 trans // gR - 3     ", S2.get_T()(2), gR2.get_T()(2), 1e-12);
            check("exp1 trans // homo - 1   ", S2.get_T()(0), M2.get_T()(0), 1e-12);
            check("exp1 trans // homo - 2   ", S2.get_T()(1), M2.get_T()(1), 1e-12);
            check("exp1 trans // homo - 3   ", S2.get_T()(2), M2.get_T()(2), 1e-12);
            check("exp1 trans // trfv - 1   ", S2.get_T()(0), tau2.get_T()(0), 1e-12);
            check("exp1 trans // trfv - 2   ", S2.get_T()(1), tau2.get_T()(1), 1e-12);
            check("exp1 trans // trfv - 3   ", S2.get_T()(2), tau2.get_T()(2), 1e-12);
            check("exp1 trans // dual - 1   ", S2.get_T()(0), Z2.get_T()(0), 1e-12);
            check("exp1 trans // dual - 2   ", S2.get_T()(1), Z2.get_T()(1), 1e-12);
            check("exp1 trans // dual - 3   ", S2.get_T()(2), Z2.get_T()(2), 1e-12);

            ang::speu_rodrigues gq3 = gq.pow(3.3);
            ang::speu_dcm gR3 = gR.pow(3.3);
            ang::homogeneous M3 = M.pow(3.3);
            ang::trfv tau3 = tau.pow(3.3);
            ang::dual Z3 = Z.pow(3.3);
            ang::screw S3 = S.pow(3.3);

            check("exp2 rotv // gq - 1     ", S3.get_rotv()()(0), gq3.get_rotv()()(0), 1e-12);
            check("exp2 rotv // gq - 2     ", S3.get_rotv()()(1), gq3.get_rotv()()(1), 1e-12);
            check("exp2 rotv // gq - 3     ", S3.get_rotv()()(2), gq3.get_rotv()()(2), 1e-12);
            check("exp2 rotv // gR - 1     ", S3.get_rotv()()(0), gR3.get_rotv()()(0), 1e-12);
            check("exp2 rotv // gR - 2     ", S3.get_rotv()()(1), gR3.get_rotv()()(1), 1e-12);
            check("exp2 rotv // gR - 3     ", S3.get_rotv()()(2), gR3.get_rotv()()(2), 1e-12);
            check("exp2 rotv // homo - 1   ", S3.get_rotv()()(0), M3.get_rotv()()(0), 1e-12);
            check("exp2 rotv // homo - 2   ", S3.get_rotv()()(1), M3.get_rotv()()(1), 1e-12);
            check("exp2 rotv // homo - 3   ", S3.get_rotv()()(2), M3.get_rotv()()(2), 1e-12);
            check("exp2 rotv // trfv - 1   ", S3.get_rotv()()(0), tau3.get_rotv()()(0), 1e-12);
            check("exp2 rotv // trfv - 2   ", S3.get_rotv()()(1), tau3.get_rotv()()(1), 1e-12);
            check("exp2 rotv // trfv - 3   ", S3.get_rotv()()(2), tau3.get_rotv()()(2), 1e-12);
            check("exp2 rotv // dual - 1   ", S3.get_rotv()()(0), Z3.get_rotv()()(0), 1e-12);
            check("exp2 rotv // dual - 2   ", S3.get_rotv()()(1), Z3.get_rotv()()(1), 1e-12);
            check("exp2 rotv // dual - 3   ", S3.get_rotv()()(2), Z3.get_rotv()()(2), 1e-12);

            check("exp2 trans // gq - 1     ", S3.get_T()(0), gq3.get_T()(0), 1e-12);
            check("exp2 trans // gq - 2     ", S3.get_T()(1), gq3.get_T()(1), 1e-12);
            check("exp2 trans // gq - 3     ", S3.get_T()(2), gq3.get_T()(2), 1e-12);
            check("exp2 trans // gR - 1     ", S3.get_T()(0), gR3.get_T()(0), 1e-12);
            check("exp2 trans // gR - 2     ", S3.get_T()(1), gR3.get_T()(1), 1e-12);
            check("exp2 trans // gR - 3     ", S3.get_T()(2), gR3.get_T()(2), 1e-12);
            check("exp2 trans // homo - 1   ", S3.get_T()(0), M3.get_T()(0), 1e-12);
            check("exp2 trans // homo - 2   ", S3.get_T()(1), M3.get_T()(1), 1e-12);
            check("exp2 trans // homo - 3   ", S3.get_T()(2), M3.get_T()(2), 1e-12);
            check("exp2 trans // trfv - 1   ", S3.get_T()(0), tau3.get_T()(0), 1e-12);
            check("exp2 trans // trfv - 2   ", S3.get_T()(1), tau3.get_T()(1), 1e-12);
            check("exp2 trans // trfv - 3   ", S3.get_T()(2), tau3.get_T()(2), 1e-12);
            check("exp2 trans // dual - 1   ", S3.get_T()(0), Z3.get_T()(0), 1e-12);
            check("exp2 trans // dual - 2   ", S3.get_T()(1), Z3.get_T()(1), 1e-12);
            check("exp2 trans // dual - 3   ", S3.get_T()(2), Z3.get_T()(2), 1e-12);

            ang::speu_rodrigues gq4 = gq.pow(-0.3);
            ang::speu_dcm gR4 = gR.pow(-0.3);
            ang::homogeneous M4 = M.pow(-0.3);
            ang::trfv tau4 = tau.pow(-0.3);
            ang::dual Z4 = Z.pow(-0.3);
            ang::screw S4 = S.pow(-0.3);
            std::cout << S.get_rotv()() << std::endl;
            std::cout << tau4.get_rotv()() << std::endl;
            std::cout << S4.get_rotv()() << std::endl;

            check("exp3 rotv // gq - 1     ", S4.get_rotv()()(0), gq4.get_rotv()()(0), 1e-12);
            check("exp3 rotv // gq - 2     ", S4.get_rotv()()(1), gq4.get_rotv()()(1), 1e-12);
            check("exp3 rotv // gq - 3     ", S4.get_rotv()()(2), gq4.get_rotv()()(2), 1e-12);
            check("exp3 rotv // gR - 1     ", S4.get_rotv()()(0), gR4.get_rotv()()(0), 1e-12);
            check("exp3 rotv // gR - 2     ", S4.get_rotv()()(1), gR4.get_rotv()()(1), 1e-12);
            check("exp3 rotv // gR - 3     ", S4.get_rotv()()(2), gR4.get_rotv()()(2), 1e-12);
            check("exp3 rotv // homo - 1   ", S4.get_rotv()()(0), M4.get_rotv()()(0), 1e-12);
            check("exp3 rotv // homo - 2   ", S4.get_rotv()()(1), M4.get_rotv()()(1), 1e-12);
            check("exp3 rotv // homo - 3   ", S4.get_rotv()()(2), M4.get_rotv()()(2), 1e-12);
            check("exp3 rotv // trfv - 1   ", S4.get_rotv()()(0), tau4.get_rotv()()(0), 1e-12);
            check("exp3 rotv // trfv - 2   ", S4.get_rotv()()(1), tau4.get_rotv()()(1), 1e-12);
            check("exp3 rotv // trfv - 3   ", S4.get_rotv()()(2), tau4.get_rotv()()(2), 1e-12);
            check("exp3 rotv // dual - 1   ", S4.get_rotv()()(0), Z4.get_rotv()()(0), 1e-12);
            check("exp3 rotv // dual - 2   ", S4.get_rotv()()(1), Z4.get_rotv()()(1), 1e-12);
            check("exp3 rotv // dual - 3   ", S4.get_rotv()()(2), Z4.get_rotv()()(2), 1e-12);

            check("exp3 trans // gq - 1     ", S4.get_T()(0), gq4.get_T()(0), 1e-12);
            check("exp3 trans // gq - 2     ", S4.get_T()(1), gq4.get_T()(1), 1e-12);
            check("exp3 trans // gq - 3     ", S4.get_T()(2), gq4.get_T()(2), 1e-12);
            check("exp3 trans // gR - 1     ", S4.get_T()(0), gR4.get_T()(0), 1e-12);
            check("exp3 trans // gR - 2     ", S4.get_T()(1), gR4.get_T()(1), 1e-12);
            check("exp3 trans // gR - 3     ", S4.get_T()(2), gR4.get_T()(2), 1e-12);
            check("exp3 trans // homo - 1   ", S4.get_T()(0), M4.get_T()(0), 1e-12);
            check("exp3 trans // homo - 2   ", S4.get_T()(1), M4.get_T()(1), 1e-12);
            check("exp3 trans // homo - 3   ", S4.get_T()(2), M4.get_T()(2), 1e-12);
            check("exp3 trans // trfv - 1   ", S4.get_T()(0), tau4.get_T()(0), 1e-12);
            check("exp3 trans // trfv - 2   ", S4.get_T()(1), tau4.get_T()(1), 1e-12);
            check("exp3 trans // trfv - 3   ", S4.get_T()(2), tau4.get_T()(2), 1e-12);
            check("exp3 trans // dual - 1   ", S4.get_T()(0), Z4.get_T()(0), 1e-12);
            check("exp3 trans // dual - 2   ", S4.get_T()(1), Z4.get_T()(1), 1e-12);
            check("exp3 trans // dual - 3   ", S4.get_T()(2), Z4.get_T()(2), 1e-12);

            ang::speu_rodrigues gq5 = gq.pow(-3.3);
            ang::speu_dcm gR5 = gR.pow(-3.3);
            ang::homogeneous M5 = M.pow(-3.3);
            ang::trfv tau5 = tau.pow(-3.3);
            ang::dual Z5 = Z.pow(-3.3);
            ang::screw S5 = S.pow(-3.3);

            check("exp4 rotv // gq - 1     ", S5.get_rotv()()(0), gq5.get_rotv()()(0), 1e-12);
            check("exp4 rotv // gq - 2     ", S5.get_rotv()()(1), gq5.get_rotv()()(1), 1e-12);
            check("exp4 rotv // gq - 3     ", S5.get_rotv()()(2), gq5.get_rotv()()(2), 1e-12);
            check("exp4 rotv // gR - 1     ", S5.get_rotv()()(0), gR5.get_rotv()()(0), 1e-12);
            check("exp4 rotv // gR - 2     ", S5.get_rotv()()(1), gR5.get_rotv()()(1), 1e-12);
            check("exp4 rotv // gR - 3     ", S5.get_rotv()()(2), gR5.get_rotv()()(2), 1e-12);
            check("exp4 rotv // homo - 1   ", S5.get_rotv()()(0), M5.get_rotv()()(0), 1e-12);
            check("exp4 rotv // homo - 2   ", S5.get_rotv()()(1), M5.get_rotv()()(1), 1e-12);
            check("exp4 rotv // homo - 3   ", S5.get_rotv()()(2), M5.get_rotv()()(2), 1e-12);
            check("exp4 rotv // trfv - 1   ", S5.get_rotv()()(0), tau5.get_rotv()()(0), 1e-12);
            check("exp4 rotv // trfv - 2   ", S5.get_rotv()()(1), tau5.get_rotv()()(1), 1e-12);
            check("exp4 rotv // trfv - 3   ", S5.get_rotv()()(2), tau5.get_rotv()()(2), 1e-12);
            check("exp4 rotv // dual - 1   ", S5.get_rotv()()(0), Z5.get_rotv()()(0), 1e-12);
            check("exp4 rotv // dual - 2   ", S5.get_rotv()()(1), Z5.get_rotv()()(1), 1e-12);
            check("exp4 rotv // dual - 3   ", S5.get_rotv()()(2), Z5.get_rotv()()(2), 1e-12);

            check("exp4 trans // gq - 1     ", S5.get_T()(0), gq5.get_T()(0), 1e-12);
            check("exp4 trans // gq - 2     ", S5.get_T()(1), gq5.get_T()(1), 1e-12);
            check("exp4 trans // gq - 3     ", S5.get_T()(2), gq5.get_T()(2), 1e-12);
            check("exp4 trans // gR - 1     ", S5.get_T()(0), gR5.get_T()(0), 1e-12);
            check("exp4 trans // gR - 2     ", S5.get_T()(1), gR5.get_T()(1), 1e-12);
            check("exp4 trans // gR - 3     ", S5.get_T()(2), gR5.get_T()(2), 1e-12);
            check("exp4 trans // homo - 1   ", S5.get_T()(0), M5.get_T()(0), 1e-12);
            check("exp4 trans // homo - 2   ", S5.get_T()(1), M5.get_T()(1), 1e-12);
            check("exp4 trans // homo - 3   ", S5.get_T()(2), M5.get_T()(2), 1e-12);
            check("exp4 trans // trfv - 1   ", S5.get_T()(0), tau5.get_T()(0), 1e-12);
            check("exp4 trans // trfv - 2   ", S5.get_T()(1), tau5.get_T()(1), 1e-12);
            check("exp4 trans // trfv - 3   ", S5.get_T()(2), tau5.get_T()(2), 1e-12);
            check("exp4 trans // dual - 1   ", S5.get_T()(0), Z5.get_T()(0), 1e-12);
            check("exp4 trans // dual - 2   ", S5.get_T()(1), Z5.get_T()(1), 1e-12);
            check("exp4 trans // dual - 3   ", S5.get_T()(2), Z5.get_T()(2), 1e-12);

            ang::speu_rodrigues gq6 = gq.pow(0.);
            ang::speu_dcm gR6 = gR.pow(0.);
            ang::homogeneous M6 = M.pow(0.);
            ang::trfv tau6 = tau.pow(0.);
            ang::dual Z6 = Z.pow(0.);
            ang::screw S6 = S.pow(0.);

            check("exp5 rotv // gq - 1     ", S6.get_rotv()()(0), gq6.get_rotv()()(0), 1e-12);
            check("exp5 rotv // gq - 2     ", S6.get_rotv()()(1), gq6.get_rotv()()(1), 1e-12);
            check("exp5 rotv // gq - 3     ", S6.get_rotv()()(2), gq6.get_rotv()()(2), 1e-12);
            check("exp5 rotv // gR - 1     ", S6.get_rotv()()(0), gR6.get_rotv()()(0), 1e-12);
            check("exp5 rotv // gR - 2     ", S6.get_rotv()()(1), gR6.get_rotv()()(1), 1e-12);
            check("exp5 rotv // gR - 3     ", S6.get_rotv()()(2), gR6.get_rotv()()(2), 1e-12);
            check("exp5 rotv // homo - 1   ", S6.get_rotv()()(0), M6.get_rotv()()(0), 1e-12);
            check("exp5 rotv // homo - 2   ", S6.get_rotv()()(1), M6.get_rotv()()(1), 1e-12);
            check("exp5 rotv // homo - 3   ", S6.get_rotv()()(2), M6.get_rotv()()(2), 1e-12);
            check("exp5 rotv // trfv - 1   ", S6.get_rotv()()(0), tau6.get_rotv()()(0), 1e-12);
            check("exp5 rotv // trfv - 2   ", S6.get_rotv()()(1), tau6.get_rotv()()(1), 1e-12);
            check("exp5 rotv // trfv - 3   ", S6.get_rotv()()(2), tau6.get_rotv()()(2), 1e-12);
            check("exp5 rotv // dual - 1   ", S6.get_rotv()()(0), Z6.get_rotv()()(0), 1e-12);
            check("exp5 rotv // dual - 2   ", S6.get_rotv()()(1), Z6.get_rotv()()(1), 1e-12);
            check("exp5 rotv // dual - 3   ", S6.get_rotv()()(2), Z6.get_rotv()()(2), 1e-12);

            check("exp5 trans // gq - 1     ", S6.get_T()(0), gq6.get_T()(0), 1e-12);
            check("exp5 trans // gq - 2     ", S6.get_T()(1), gq6.get_T()(1), 1e-12);
            check("exp5 trans // gq - 3     ", S6.get_T()(2), gq6.get_T()(2), 1e-12);
            check("exp5 trans // gR - 1     ", S6.get_T()(0), gR6.get_T()(0), 1e-12);
            check("exp5 trans // gR - 2     ", S6.get_T()(1), gR6.get_T()(1), 1e-12);
            check("exp5 trans // gR - 3     ", S6.get_T()(2), gR6.get_T()(2), 1e-12);
            check("exp5 trans // homo - 1   ", S6.get_T()(0), M6.get_T()(0), 1e-12);
            check("exp5 trans // homo - 2   ", S6.get_T()(1), M6.get_T()(1), 1e-12);
            check("exp5 trans // homo - 3   ", S6.get_T()(2), M6.get_T()(2), 1e-12);
            check("exp5 trans // trfv - 1   ", S6.get_T()(0), tau6.get_T()(0), 1e-12);
            check("exp5 trans // trfv - 2   ", S6.get_T()(1), tau6.get_T()(1), 1e-12);
            check("exp5 trans // trfv - 3   ", S6.get_T()(2), tau6.get_T()(2), 1e-12);
            check("exp5 trans // dual - 1   ", S6.get_T()(0), Z6.get_T()(0), 1e-12);
            check("exp5 trans // dual - 2   ", S6.get_T()(1), Z6.get_T()(1), 1e-12);
            check("exp5 trans // dual - 3   ", S6.get_T()(2), Z6.get_T()(2), 1e-12);

            check("exp5 speu rodrigues // gq - 1     ", gq6.get_rodrigues()()(0), 1., 1e-12);
            check("exp5 speu rodrigues // gq - 2     ", gq6.get_rodrigues()()(1), 0., 1e-12);
            check("exp5 speu rodrigues // gq - 3     ", gq6.get_rodrigues()()(2), 0., 1e-12);
            check("exp5 speu rodrigues // gq - 4     ", gq6.get_rodrigues()()(3), 0., 1e-12);
            check("exp5 speu rodrigues // gq - 1     ", gq6.get_T()(0), 0., 1e-12);
            check("exp5 speu rodrigues // gq - 2     ", gq6.get_T()(1), 0., 1e-12);
            check("exp5 speu rodrigues // gq - 3     ", gq6.get_T()(2), 0., 1e-12);

            check("exp5 speu dcm // gR - 11     ", gR6.get_dcm()()(0,0), 1., 1e-12);
            check("exp5 speu dcm // gR - 12     ", gR6.get_dcm()()(0,1), 0., 1e-12);
            check("exp5 speu dcm // gR - 13     ", gR6.get_dcm()()(0,2), 0., 1e-12);
            check("exp5 speu dcm // gR - 21     ", gR6.get_dcm()()(1,0), 0., 1e-12);
            check("exp5 speu dcm // gR - 22     ", gR6.get_dcm()()(1,1), 1., 1e-12);
            check("exp5 speu dcm // gR - 23     ", gR6.get_dcm()()(1,2), 0., 1e-12);
            check("exp5 speu dcm // gR - 31     ", gR6.get_dcm()()(2,0), 0., 1e-12);
            check("exp5 speu dcm // gR - 32     ", gR6.get_dcm()()(2,1), 0., 1e-12);
            check("exp5 speu dcm // gR - 33     ", gR6.get_dcm()()(2,2), 1., 1e-12);
            check("exp5 speu dcm // gR - 1      ", gq6.get_T()(0), 0., 1e-12);
            check("exp5 speu dcm // gR - 2      ", gq6.get_T()(1), 0., 1e-12);
            check("exp5 speu dcm // gR - 3      ", gq6.get_T()(2), 0., 1e-12);

            check("exp5 homogeneous // M - 11   ", M6()(0,0), 1., 1e-12);
            check("exp5 homogeneous // M - 12   ", M6()(0,1), 0., 1e-12);
            check("exp5 homogeneous // M - 13   ", M6()(0,2), 0., 1e-12);
            check("exp5 homogeneous // M - 14   ", M6()(0,3), 0., 1e-12);
            check("exp5 homogeneous // M - 21   ", M6()(1,0), 0., 1e-12);
            check("exp5 homogeneous // M - 22   ", M6()(1,1), 1., 1e-12);
            check("exp5 homogeneous // M - 23   ", M6()(1,2), 0., 1e-12);
            check("exp5 homogeneous // M - 24   ", M6()(1,3), 0., 1e-12);
            check("exp5 homogeneous // M - 31   ", M6()(2,0), 0., 1e-12);
            check("exp5 homogeneous // M - 32   ", M6()(2,1), 0., 1e-12);
            check("exp5 homogeneous // M - 33   ", M6()(2,2), 1., 1e-12);
            check("exp5 homogeneous // M - 34   ", M6()(2,3), 0., 1e-12);
            check("exp5 homogeneous // M - 41   ", M6()(3,0), 0., 1e-12);
            check("exp5 homogeneous // M - 42   ", M6()(3,1), 0., 1e-12);
            check("exp5 homogeneous // M - 43   ", M6()(3,2), 0., 1e-12);
            check("exp5 homogeneous // M - 44   ", M6()(3,3), 1., 1e-12);

            check("exp5 trfv // tau - 1         ", tau6.get_s()(0), 0., 1e-12);
            check("exp5 trfv // tau - 2         ", tau6.get_s()(1), 0., 1e-12);
            check("exp5 trfv // tau - 3         ", tau6.get_s()(2), 0., 1e-12);
            check("exp5 trfv // tau - 1         ", tau6.get_rotv()()(0), 0., 1e-12);
            check("exp5 trfv // tau - 2         ", tau6.get_rotv()()(1), 0., 1e-12);
            check("exp5 trfv // tau - 3         ", tau6.get_rotv()()(2), 0., 1e-12);

            check("exp5 dual // zeta - 1        ", Z6.get_qr()()(0), 1., 1e-12);
            check("exp5 dual // zeta - 2        ", Z6.get_qr()()(1), 0., 1e-12);
            check("exp5 dual // zeta - 3        ", Z6.get_qr()()(2), 0., 1e-12);
            check("exp5 dual // zeta - 4        ", Z6.get_qr()()(3), 0., 1e-12);
            check("exp5 dual // zeta - 5        ", Z6.get_qd()(0), 0., 1e-12);
            check("exp5 dual // zeta - 6        ", Z6.get_qd()(1), 0., 1e-12);
            check("exp5 dual // zeta - 7        ", Z6.get_qd()(2), 0., 1e-12);
            check("exp5 dual // zeta - 8        ", Z6.get_qd()(3), 0., 1e-12);

            check("exp5 screw // S - 1          ", std::isnan(S6.get_h()), true, 1e-12);
            check("exp5 screw // S - 2          ", S6.get_phi(), 0., 1e-12);

            ang::speu_rodrigues gq7 = gq.pow(1.);
            ang::speu_dcm gR7 = gR.pow(1.);
            ang::homogeneous M7 = M.pow(1.);
            ang::trfv tau7 = tau.pow(1.);
            ang::dual Z7 = Z.pow(1.);
            ang::screw S7 = S.pow(1.);

            check("exp6 rotv // gq - 1     ", S7.get_rotv()()(0), gq7.get_rotv()()(0), 1e-12);
            check("exp6 rotv // gq - 2     ", S7.get_rotv()()(1), gq7.get_rotv()()(1), 1e-12);
            check("exp6 rotv // gq - 3     ", S7.get_rotv()()(2), gq7.get_rotv()()(2), 1e-12);
            check("exp6 rotv // gR - 1     ", S7.get_rotv()()(0), gR7.get_rotv()()(0), 1e-12);
            check("exp6 rotv // gR - 2     ", S7.get_rotv()()(1), gR7.get_rotv()()(1), 1e-12);
            check("exp6 rotv // gR - 3     ", S7.get_rotv()()(2), gR7.get_rotv()()(2), 1e-12);
            check("exp6 rotv // homo - 1   ", S7.get_rotv()()(0), M7.get_rotv()()(0), 1e-12);
            check("exp6 rotv // homo - 2   ", S7.get_rotv()()(1), M7.get_rotv()()(1), 1e-12);
            check("exp6 rotv // homo - 3   ", S7.get_rotv()()(2), M7.get_rotv()()(2), 1e-12);
            check("exp6 rotv // trfv - 1   ", S7.get_rotv()()(0), tau7.get_rotv()()(0), 1e-12);
            check("exp6 rotv // trfv - 2   ", S7.get_rotv()()(1), tau7.get_rotv()()(1), 1e-12);
            check("exp6 rotv // trfv - 3   ", S7.get_rotv()()(2), tau7.get_rotv()()(2), 1e-12);
            check("exp6 rotv // dual - 1   ", S7.get_rotv()()(0), Z7.get_rotv()()(0), 1e-12);
            check("exp6 rotv // dual - 2   ", S7.get_rotv()()(1), Z7.get_rotv()()(1), 1e-12);
            check("exp6 rotv // dual - 3   ", S7.get_rotv()()(2), Z7.get_rotv()()(2), 1e-12);

            check("exp6 trans // gq - 1     ", S7.get_T()(0), gq7.get_T()(0), 1e-12);
            check("exp6 trans // gq - 2     ", S7.get_T()(1), gq7.get_T()(1), 1e-12);
            check("exp6 trans // gq - 3     ", S7.get_T()(2), gq7.get_T()(2), 1e-12);
            check("exp6 trans // gR - 1     ", S7.get_T()(0), gR7.get_T()(0), 1e-12);
            check("exp6 trans // gR - 2     ", S7.get_T()(1), gR7.get_T()(1), 1e-12);
            check("exp6 trans // gR - 3     ", S7.get_T()(2), gR7.get_T()(2), 1e-12);
            check("exp6 trans // homo - 1   ", S7.get_T()(0), M7.get_T()(0), 1e-12);
            check("exp6 trans // homo - 2   ", S7.get_T()(1), M7.get_T()(1), 1e-12);
            check("exp6 trans // homo - 3   ", S7.get_T()(2), M7.get_T()(2), 1e-12);
            check("exp6 trans // trfv - 1   ", S7.get_T()(0), tau7.get_T()(0), 1e-12);
            check("exp6 trans // trfv - 2   ", S7.get_T()(1), tau7.get_T()(1), 1e-12);
            check("exp6 trans // trfv - 3   ", S7.get_T()(2), tau7.get_T()(2), 1e-12);
            check("exp6 trans // dual - 1   ", S7.get_T()(0), Z7.get_T()(0), 1e-12);
            check("exp6 trans // dual - 2   ", S7.get_T()(1), Z7.get_T()(1), 1e-12);
            check("exp6 trans // dual - 3   ", S7.get_T()(2), Z7.get_T()(2), 1e-12);

            check("exp6 speu rodrigues // gq - 1     ", gq7.get_rodrigues()()(0), gq.get_rodrigues()()(0), 1e-12);
            check("exp6 speu rodrigues // gq - 2     ", gq7.get_rodrigues()()(1), gq.get_rodrigues()()(1), 1e-12);
            check("exp6 speu rodrigues // gq - 3     ", gq7.get_rodrigues()()(2), gq.get_rodrigues()()(2), 1e-12);
            check("exp6 speu rodrigues // gq - 4     ", gq7.get_rodrigues()()(3), gq.get_rodrigues()()(3), 1e-12);
            check("exp6 speu rodrigues // gq - 1     ", gq7.get_T()(0), gq.get_T()(0), 1e-12);
            check("exp6 speu rodrigues // gq - 2     ", gq7.get_T()(1), gq.get_T()(1), 1e-12);
            check("exp6 speu rodrigues // gq - 3     ", gq7.get_T()(2), gq.get_T()(2), 1e-12);

            check("exp6 speu dcm // gR - 11     ", gR7.get_dcm()()(0,0), gR.get_dcm()()(0,0), 1e-12);
            check("exp6 speu dcm // gR - 12     ", gR7.get_dcm()()(0,1), gR.get_dcm()()(0,1), 1e-12);
            check("exp6 speu dcm // gR - 13     ", gR7.get_dcm()()(0,2), gR.get_dcm()()(0,2), 1e-12);
            check("exp6 speu dcm // gR - 21     ", gR7.get_dcm()()(1,0), gR.get_dcm()()(1,0), 1e-12);
            check("exp6 speu dcm // gR - 22     ", gR7.get_dcm()()(1,1), gR.get_dcm()()(1,1), 1e-12);
            check("exp6 speu dcm // gR - 23     ", gR7.get_dcm()()(1,2), gR.get_dcm()()(1,2), 1e-12);
            check("exp6 speu dcm // gR - 31     ", gR7.get_dcm()()(2,0), gR.get_dcm()()(2,0), 1e-12);
            check("exp6 speu dcm // gR - 32     ", gR7.get_dcm()()(2,1), gR.get_dcm()()(2,1), 1e-12);
            check("exp6 speu dcm // gR - 33     ", gR7.get_dcm()()(2,2), gR.get_dcm()()(2,2), 1e-12);
            check("exp6 speu dcm // gR - 1      ", gR7.get_T()(0), gR.get_T()(0), 1e-12);
            check("exp6 speu dcm // gR - 2      ", gR7.get_T()(1), gR.get_T()(1), 1e-12);
            check("exp6 speu dcm // gR - 3      ", gR7.get_T()(2), gR.get_T()(2), 1e-12);

            check("exp6 homogeneous // M - 11   ", M7()(0,0), M()(0,0), 1e-12);
            check("exp6 homogeneous // M - 12   ", M7()(0,1), M()(0,1), 1e-12);
            check("exp6 homogeneous // M - 13   ", M7()(0,2), M()(0,2), 1e-12);
            check("exp6 homogeneous // M - 14   ", M7()(0,3), M()(0,3), 1e-12);
            check("exp6 homogeneous // M - 21   ", M7()(1,0), M()(1,0), 1e-12);
            check("exp6 homogeneous // M - 22   ", M7()(1,1), M()(1,1), 1e-12);
            check("exp6 homogeneous // M - 23   ", M7()(1,2), M()(1,2), 1e-12);
            check("exp6 homogeneous // M - 24   ", M7()(1,3), M()(1,3), 1e-12);
            check("exp6 homogeneous // M - 31   ", M7()(2,0), M()(2,0), 1e-12);
            check("exp6 homogeneous // M - 32   ", M7()(2,1), M()(2,1), 1e-12);
            check("exp6 homogeneous // M - 33   ", M7()(2,2), M()(2,2), 1e-12);
            check("exp6 homogeneous // M - 34   ", M7()(2,3), M()(2,3), 1e-12);
            check("exp6 homogeneous // M - 41   ", M7()(3,0), M()(3,0), 1e-12);
            check("exp6 homogeneous // M - 42   ", M7()(3,1), M()(3,1), 1e-12);
            check("exp6 homogeneous // M - 43   ", M7()(3,2), M()(3,2), 1e-12);
            check("exp6 homogeneous // M - 44   ", M7()(3,3), M()(3,3), 1e-12);

            check("exp6 trfv // tau - 1         ", tau7.get_s()(0), tau.get_s()(0), 1e-12);
            check("exp6 trfv // tau - 2         ", tau7.get_s()(1), tau.get_s()(1), 1e-12);
            check("exp6 trfv // tau - 3         ", tau7.get_s()(2), tau.get_s()(2), 1e-12);
            check("exp6 trfv // tau - 1         ", tau7.get_rotv()()(0), tau.get_rotv()()(0), 1e-12);
            check("exp6 trfv // tau - 2         ", tau7.get_rotv()()(1), tau.get_rotv()()(1), 1e-12);
            check("exp6 trfv // tau - 3         ", tau7.get_rotv()()(2), tau.get_rotv()()(2), 1e-12);

            check("exp6 dual // zeta - 1        ", Z7.get_qr()()(0), Z.get_qr()()(0), 1e-12);
            check("exp6 dual // zeta - 2        ", Z7.get_qr()()(1), Z.get_qr()()(1), 1e-12);
            check("exp6 dual // zeta - 3        ", Z7.get_qr()()(2), Z.get_qr()()(2), 1e-12);
            check("exp6 dual // zeta - 4        ", Z7.get_qr()()(3), Z.get_qr()()(3), 1e-12);
            check("exp6 dual // zeta - 5        ", Z7.get_qd()(0), Z.get_qd()(0), 1e-12);
            check("exp6 dual // zeta - 6        ", Z7.get_qd()(1), Z.get_qd()(1), 1e-12);
            check("exp6 dual // zeta - 7        ", Z7.get_qd()(2), Z.get_qd()(2), 1e-12);
            check("exp6 dual // zeta - 8        ", Z7.get_qd()(3), Z.get_qd()(3), 1e-12);

            check("exp6 screw // S - 1          ", S7.get_n()(0), S.get_n()(0), 1e-12);
            check("exp6 screw // S - 2          ", S7.get_n()(1), S.get_n()(1), 1e-12);
            check("exp6 screw // S - 3          ", S7.get_n()(2), S.get_n()(2), 1e-12);
            check("exp6 screw // S - phi        ", S7.get_phi(), S.get_phi(), 1e-12);
            check("exp6 screw // S - d          ", S7.get_d(), S.get_d(), 1e-12);
            // tests with m do not work because values are huge for small angles, but they are the same as inputs

        }
    }
} // closes test_power

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_sclerp() {

    Eigen::Vector3d l_big(3., 8., 5.);
    Eigen::Vector3d l = l_big / l_big.norm();
    Eigen::Vector3d p(2.0, -1.7, 0.3);
    Eigen::Vector3d m = p.cross(l);
    double phi = 0.2;
    double d = 3.3;
    double h = d / phi;

    // I make them multiples so I can manually compute the interpolated values
    ang::screw S00(l, m, phi, d, h);
    ang::screw S22(l, m, phi * 2., d * 2., h);
    ang::screw S11(l, m, phi * 1.2, d * 1.2, h);

    ang::screw S0 = ang::screw::sclerp(S00, S22, 0.);
    ang::screw S1 = ang::screw::sclerp(S00, S22, 0.2);
    ang::screw S2 = ang::screw::sclerp(S00, S22, 1.);

    check("sclerp screw l 00   ", S00.get_n()(0), S0.get_n()(0), 1e-12);
    check("sclerp screw l 01   ", S00.get_n()(1), S0.get_n()(1), 1e-12);
    check("sclerp screw l 02   ", S00.get_n()(2), S0.get_n()(2), 1e-12);
    check("sclerp screw m 00   ", S00.get_m()(0), S0.get_m()(0), 1e-12);
    check("sclerp screw m 01   ", S00.get_m()(1), S0.get_m()(1), 1e-12);
    check("sclerp screw m 02   ", S00.get_m()(2), S0.get_m()(2), 1e-12);
    check("sclerp screw phi 0  ", S00.get_phi(),  S0.get_phi(),  1e-12);
    check("sclerp screw d 0    ", S00.get_d(),    S0.get_d(),    1e-12);
    check("sclerp screw h 0    ", S00.get_h(),    S0.get_h(),    1e-12);

    check("sclerp screw l 10   ", S11.get_n()(0), S1.get_n()(0), 1e-12);
    check("sclerp screw l 11   ", S11.get_n()(1), S1.get_n()(1), 1e-12);
    check("sclerp screw l 12   ", S11.get_n()(2), S1.get_n()(2), 1e-12);
    check("sclerp screw m 10   ", S11.get_m()(0), S1.get_m()(0), 1e-12);
    check("sclerp screw m 11   ", S11.get_m()(1), S1.get_m()(1), 1e-12);
    check("sclerp screw m 12   ", S11.get_m()(2), S1.get_m()(2), 1e-12);
    check("sclerp screw phi 1  ", S11.get_phi(),  S1.get_phi(),  1e-12);
    check("sclerp screw d 1    ", S11.get_d(),    S1.get_d(),    1e-12);
    check("sclerp screw h 1    ", S11.get_h(),    S1.get_h(),    1e-12);

    check("sclerp screw l 20   ", S22.get_n()(0), S2.get_n()(0), 1e-12);
    check("sclerp screw l 21   ", S22.get_n()(1), S2.get_n()(1), 1e-12);
    check("sclerp screw l 22   ", S22.get_n()(2), S2.get_n()(2), 1e-12);
    check("sclerp screw m 20   ", S22.get_m()(0), S2.get_m()(0), 1e-12);
    check("sclerp screw m 21   ", S22.get_m()(1), S2.get_m()(1), 1e-12);
    check("sclerp screw m 22   ", S22.get_m()(2), S2.get_m()(2), 1e-12);
    check("sclerp screw phi 2  ", S22.get_phi(),  S2.get_phi(),  1e-12);
    check("sclerp screw d 2    ", S22.get_d(),    S2.get_d(),    1e-12);
    check("sclerp screw h 2    ", S22.get_h(),    S2.get_h(),    1e-12);

    ang::speu_rodrigues Ogq00(S00);
    ang::speu_rodrigues Ogq22(S22);
    ang::speu_rodrigues Ogq11(S11);

    ang::speu_rodrigues Ogq0 = ang::speu_rodrigues::sclerp(Ogq00, Ogq22, 0.);
    ang::speu_rodrigues Ogq1 = ang::speu_rodrigues::sclerp(Ogq00, Ogq22, 0.2);
    ang::speu_rodrigues Ogq2 = ang::speu_rodrigues::sclerp(Ogq00, Ogq22, 1.);

    check("sclerp speu_rodrig 00 ", Ogq00.get_rotv()()(0), Ogq0.get_rotv()()(0), 1e-12);
    check("sclerp speu_rodrig 01 ", Ogq00.get_rotv()()(1), Ogq0.get_rotv()()(1), 1e-12);
    check("sclerp speu_rodrig 02 ", Ogq00.get_rotv()()(2), Ogq0.get_rotv()()(2), 1e-12);
    check("sclerp speu_rodrig 00 ", Ogq00.get_T()(0), Ogq0.get_T()(0), 1e-12);
    check("sclerp speu_rodrig 01 ", Ogq00.get_T()(1), Ogq0.get_T()(1), 1e-12);
    check("sclerp speu_rodrig 02 ", Ogq00.get_T()(2), Ogq0.get_T()(2), 1e-12);

    check("sclerp speu_rodrig 10 ", Ogq11.get_rotv()()(0), Ogq1.get_rotv()()(0), 1e-12);
    check("sclerp speu_rodrig 11 ", Ogq11.get_rotv()()(1), Ogq1.get_rotv()()(1), 1e-12);
    check("sclerp speu_rodrig 12 ", Ogq11.get_rotv()()(2), Ogq1.get_rotv()()(2), 1e-12);
    check("sclerp speu_rodrig 10 ", Ogq11.get_T()(0), Ogq1.get_T()(0), 1e-12);
    check("sclerp speu_rodrig 11 ", Ogq11.get_T()(1), Ogq1.get_T()(1), 1e-12);
    check("sclerp speu_rodrig 12 ", Ogq11.get_T()(2), Ogq1.get_T()(2), 1e-12);

    check("sclerp speu_rodrig 20 ", Ogq22.get_rotv()()(0), Ogq2.get_rotv()()(0), 1e-12);
    check("sclerp speu_rodrig 21 ", Ogq22.get_rotv()()(1), Ogq2.get_rotv()()(1), 1e-12);
    check("sclerp speu_rodrig 22 ", Ogq22.get_rotv()()(2), Ogq2.get_rotv()()(2), 1e-12);
    check("sclerp speu_rodrig 20 ", Ogq22.get_T()(0), Ogq2.get_T()(0), 1e-12);
    check("sclerp speu_rodrig 21 ", Ogq22.get_T()(1), Ogq2.get_T()(1), 1e-12);
    check("sclerp speu_rodrig 22 ", Ogq22.get_T()(2), Ogq2.get_T()(2), 1e-12);

    ang::speu_dcm OgR00(S00);
    ang::speu_dcm OgR22(S22);
    ang::speu_dcm OgR11(S11);

    ang::speu_dcm OgR0 = ang::speu_dcm::sclerp(OgR00, OgR22, 0.);
    ang::speu_dcm OgR1 = ang::speu_dcm::sclerp(OgR00, OgR22, 0.2);
    ang::speu_dcm OgR2 = ang::speu_dcm::sclerp(OgR00, OgR22, 1.);

    check("sclerp speu_dcm 00 ", OgR00.get_rotv()()(0), OgR0.get_rotv()()(0), 1e-12);
    check("sclerp speu_dcm 01 ", OgR00.get_rotv()()(1), OgR0.get_rotv()()(1), 1e-12);
    check("sclerp speu_dcm 02 ", OgR00.get_rotv()()(2), OgR0.get_rotv()()(2), 1e-12);
    check("sclerp speu_dcm 00 ", OgR00.get_T()(0), OgR0.get_T()(0), 1e-12);
    check("sclerp speu_dcm 01 ", OgR00.get_T()(1), OgR0.get_T()(1), 1e-12);
    check("sclerp speu_dcm 02 ", OgR00.get_T()(2), OgR0.get_T()(2), 1e-12);

    check("sclerp speu_dcm 10 ", OgR11.get_rotv()()(0), OgR1.get_rotv()()(0), 1e-12);
    check("sclerp speu_dcm 11 ", OgR11.get_rotv()()(1), OgR1.get_rotv()()(1), 1e-12);
    check("sclerp speu_dcm 12 ", OgR11.get_rotv()()(2), OgR1.get_rotv()()(2), 1e-12);
    check("sclerp speu_dcm 10 ", OgR11.get_T()(0), OgR1.get_T()(0), 1e-12);
    check("sclerp speu_dcm 11 ", OgR11.get_T()(1), OgR1.get_T()(1), 1e-12);
    check("sclerp speu_dcm 12 ", OgR11.get_T()(2), OgR1.get_T()(2), 1e-12);

    check("sclerp speu_dcm 20 ", OgR22.get_rotv()()(0), OgR2.get_rotv()()(0), 1e-12);
    check("sclerp speu_dcm 21 ", OgR22.get_rotv()()(1), OgR2.get_rotv()()(1), 1e-12);
    check("sclerp speu_dcm 22 ", OgR22.get_rotv()()(2), OgR2.get_rotv()()(2), 1e-12);
    check("sclerp speu_dcm 20 ", OgR22.get_T()(0), OgR2.get_T()(0), 1e-12);
    check("sclerp speu_dcm 21 ", OgR22.get_T()(1), OgR2.get_T()(1), 1e-12);
    check("sclerp speu_dcm 22 ", OgR22.get_T()(2), OgR2.get_T()(2), 1e-12);

    ang::homogeneous M00(S00);
    ang::homogeneous M22(S22);
    ang::homogeneous M11(S11);

    ang::homogeneous M0 = ang::homogeneous::sclerp(M00, M22, 0.);
    ang::homogeneous M1 = ang::homogeneous::sclerp(M00, M22, 0.2);
    ang::homogeneous M2 = ang::homogeneous::sclerp(M00, M22, 1.);

    check("sclerp homogeneous 00 ", M00.get_rotv()()(0), M0.get_rotv()()(0), 1e-12);
    check("sclerp homogeneous 01 ", M00.get_rotv()()(1), M0.get_rotv()()(1), 1e-12);
    check("sclerp homogeneous 02 ", M00.get_rotv()()(2), M0.get_rotv()()(2), 1e-12);
    check("sclerp homogeneous 00 ", M00.get_T()(0), M0.get_T()(0), 1e-12);
    check("sclerp homogeneous 01 ", M00.get_T()(1), M0.get_T()(1), 1e-12);
    check("sclerp homogeneous 02 ", M00.get_T()(2), M0.get_T()(2), 1e-12);

    check("sclerp homogeneous 10 ", M11.get_rotv()()(0), M1.get_rotv()()(0), 1e-12);
    check("sclerp homogeneous 11 ", M11.get_rotv()()(1), M1.get_rotv()()(1), 1e-12);
    check("sclerp homogeneous 12 ", M11.get_rotv()()(2), M1.get_rotv()()(2), 1e-12);
    check("sclerp homogeneous 10 ", M11.get_T()(0), M1.get_T()(0), 1e-12);
    check("sclerp homogeneous 11 ", M11.get_T()(1), M1.get_T()(1), 1e-12);
    check("sclerp homogeneous 12 ", M11.get_T()(2), M1.get_T()(2), 1e-12);

    check("sclerp homogeneous 20 ", M22.get_rotv()()(0), M2.get_rotv()()(0), 1e-12);
    check("sclerp homogeneous 21 ", M22.get_rotv()()(1), M2.get_rotv()()(1), 1e-12);
    check("sclerp homogeneous 22 ", M22.get_rotv()()(2), M2.get_rotv()()(2), 1e-12);
    check("sclerp homogeneous 20 ", M22.get_T()(0), M2.get_T()(0), 1e-12);
    check("sclerp homogeneous 21 ", M22.get_T()(1), M2.get_T()(1), 1e-12);
    check("sclerp homogeneous 22 ", M22.get_T()(2), M2.get_T()(2), 1e-12);

    ang::trfv tau00(S00);
    ang::trfv tau22(S22);
    ang::trfv tau11(S11);

    ang::trfv tau0 = ang::trfv::sclerp(tau00, tau22, 0.);
    ang::trfv tau1 = ang::trfv::sclerp(tau00, tau22, 0.2);
    ang::trfv tau2 = ang::trfv::sclerp(tau00, tau22, 1.);

    check("sclerp trfv 00 ", tau00.get_rotv()()(0), tau0.get_rotv()()(0), 1e-12);
    check("sclerp trfv 01 ", tau00.get_rotv()()(1), tau0.get_rotv()()(1), 1e-12);
    check("sclerp trfv 02 ", tau00.get_rotv()()(2), tau0.get_rotv()()(2), 1e-12);
    check("sclerp trfv 00 ", tau00.get_T()(0), tau0.get_T()(0), 1e-12);
    check("sclerp trfv 01 ", tau00.get_T()(1), tau0.get_T()(1), 1e-12);
    check("sclerp trfv 02 ", tau00.get_T()(2), tau0.get_T()(2), 1e-12);

    check("sclerp trfv 10 ", tau11.get_rotv()()(0), tau1.get_rotv()()(0), 1e-12);
    check("sclerp trfv 11 ", tau11.get_rotv()()(1), tau1.get_rotv()()(1), 1e-12);
    check("sclerp trfv 12 ", tau11.get_rotv()()(2), tau1.get_rotv()()(2), 1e-12);
    check("sclerp trfv 10 ", tau11.get_T()(0), tau1.get_T()(0), 1e-12);
    check("sclerp trfv 11 ", tau11.get_T()(1), tau1.get_T()(1), 1e-12);
    check("sclerp trfv 12 ", tau11.get_T()(2), tau1.get_T()(2), 1e-12);

    check("sclerp trfv 20 ", tau22.get_rotv()()(0), tau2.get_rotv()()(0), 1e-12);
    check("sclerp trfv 21 ", tau22.get_rotv()()(1), tau2.get_rotv()()(1), 1e-12);
    check("sclerp trfv 22 ", tau22.get_rotv()()(2), tau2.get_rotv()()(2), 1e-12);
    check("sclerp trfv 20 ", tau22.get_T()(0), tau2.get_T()(0), 1e-12);
    check("sclerp trfv 21 ", tau22.get_T()(1), tau2.get_T()(1), 1e-12);
    check("sclerp trfv 22 ", tau22.get_T()(2), tau2.get_T()(2), 1e-12);

    ang::dual Z00(S00);
    ang::dual Z22(S22);
    ang::dual Z11(S11);

    ang::dual Z0 = ang::dual::sclerp(Z00, Z22, 0.);
    ang::dual Z1 = ang::dual::sclerp(Z00, Z22, 0.2);
    ang::dual Z2 = ang::dual::sclerp(Z00, Z22, 1.);

    check("sclerp dual 00 ", Z00.get_rotv()()(0), Z0.get_rotv()()(0), 1e-12);
    check("sclerp dual 01 ", Z00.get_rotv()()(1), Z0.get_rotv()()(1), 1e-12);
    check("sclerp dual 02 ", Z00.get_rotv()()(2), Z0.get_rotv()()(2), 1e-12);
    check("sclerp dual 00 ", Z00.get_T()(0), Z0.get_T()(0), 1e-12);
    check("sclerp dual 01 ", Z00.get_T()(1), Z0.get_T()(1), 1e-12);
    check("sclerp dual 02 ", Z00.get_T()(2), Z0.get_T()(2), 1e-12);

    check("sclerp dual 10 ", Z11.get_rotv()()(0), Z1.get_rotv()()(0), 1e-12);
    check("sclerp dual 11 ", Z11.get_rotv()()(1), Z1.get_rotv()()(1), 1e-12);
    check("sclerp dual 12 ", Z11.get_rotv()()(2), Z1.get_rotv()()(2), 1e-12);
    check("sclerp dual 10 ", Z11.get_T()(0), Z1.get_T()(0), 1e-12);
    check("sclerp dual 11 ", Z11.get_T()(1), Z1.get_T()(1), 1e-12);
    check("sclerp dual 12 ", Z11.get_T()(2), Z1.get_T()(2), 1e-12);

    check("sclerp dual 20 ", Z22.get_rotv()()(0), Z2.get_rotv()()(0), 1e-12);
    check("sclerp dual 21 ", Z22.get_rotv()()(1), Z2.get_rotv()()(1), 1e-12);
    check("sclerp dual 22 ", Z22.get_rotv()()(2), Z2.get_rotv()()(2), 1e-12);
    check("sclerp dual 20 ", Z22.get_T()(0), Z2.get_T()(0), 1e-12);
    check("sclerp dual 21 ", Z22.get_T()(1), Z2.get_T()(1), 1e-12);
    check("sclerp dual 22 ", Z22.get_T()(2), Z2.get_T()(2), 1e-12);

} // closes test_sclerp

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_plus_minus() {

    // the purpose of this test is to verify that the plus and minus functions of a rotation
    ang::rotv r_a(0.2, -0.3, 0.12);
    Eigen::Vector3d x_a(13.4, -7.6, 2.8);
    ang::trfv tau_a(r_a, x_a);
    ang::speu_dcm gR_a(tau_a);
    ang::speu_rodrigues gq_a(tau_a);
    ang::homogeneous M_a(tau_a);
    ang::screw S_a(tau_a);
    ang::dual Z_a(tau_a);

    ang::rotv Delta_r(-0.01, 0.03, 0.05);
    Eigen::Vector3d Delta_x(3.1, 1.6, -12.3);
    ang::trfv Delta_tau(Delta_r, Delta_x);
    ang::screw Delta_S(Delta_tau);

    ang::speu_rodrigues gq_b1  = gq_a.plus_right(Delta_tau);
    ang::speu_dcm       gR_b1  = gR_a.plus_right(Delta_tau);
    ang::homogeneous    M_b1   = M_a.plus_right(Delta_tau);
    ang::trfv           tau_b1 = tau_a.plus_right(Delta_tau);
    ang::screw          S_b1   = S_a.plus_right(Delta_tau);
    ang::dual           Z_b1   = Z_a.plus_right(Delta_tau);

    ang::trfv tau_bq1(gq_b1);
    ang::trfv tau_bR1(gR_b1);
    ang::trfv tau_bM1(M_b1);
    ang::trfv tau_bS1(S_b1);
    ang::trfv tau_bZ1(Z_b1);

    ang::speu_rodrigues gq_b2  = gq_a.plus_left(Delta_tau);
    ang::speu_dcm       gR_b2  = gR_a.plus_left(Delta_tau);
    ang::homogeneous    M_b2   = M_a.plus_left(Delta_tau);
    ang::trfv           tau_b2 = tau_a.plus_left(Delta_tau);
    ang::screw          S_b2   = S_a.plus_left(Delta_tau);
    ang::dual           Z_b2   = Z_a.plus_left(Delta_tau);

    ang::trfv tau_bq2(gq_b2);
    ang::trfv tau_bR2(gR_b2);
    ang::trfv tau_bM2(M_b2);
    ang::trfv tau_bS2(S_b2);
    ang::trfv tau_bZ2(Z_b2);

    check("plus right sp rod - sp dcm 1      ", tau_bq1()(0), tau_bR1()(0), 1e-12);
    check("plus right sp rod - sp dcm 2      ", tau_bq1()(1), tau_bR1()(1), 1e-12);
    check("plus right sp rod - sp dcm 3      ", tau_bq1()(2), tau_bR1()(2), 1e-12);
    check("plus right sp rod - sp dcm 4      ", tau_bq1()(3), tau_bR1()(3), 1e-12);
    check("plus right sp rod - sp dcm 5      ", tau_bq1()(4), tau_bR1()(4), 1e-12);
    check("plus right sp rod - sp dcm 6      ", tau_bq1()(5), tau_bR1()(5), 1e-12);
    check("plus right sp rod - homog 1       ", tau_bq1()(0), tau_bM1()(0), 1e-12);
    check("plus right sp rod - homog 2       ", tau_bq1()(1), tau_bM1()(1), 1e-12);
    check("plus right sp rod - homog 3       ", tau_bq1()(2), tau_bM1()(2), 1e-12);
    check("plus right sp rod - homog 4       ", tau_bq1()(3), tau_bM1()(3), 1e-12);
    check("plus right sp rod - homog 5       ", tau_bq1()(4), tau_bM1()(4), 1e-12);
    check("plus right sp rod - homog 6       ", tau_bq1()(5), tau_bM1()(5), 1e-12);
    check("plus right sp rod - trfv 1        ", tau_bq1()(0), tau_b1()(0), 1e-12);
    check("plus right sp rod - trfv 2        ", tau_bq1()(1), tau_b1()(1), 1e-12);
    check("plus right sp rod - trfv 3        ", tau_bq1()(2), tau_b1()(2), 1e-12);
    check("plus right sp rod - trfv 4        ", tau_bq1()(3), tau_b1()(3), 1e-12);
    check("plus right sp rod - trfv 5        ", tau_bq1()(4), tau_b1()(4), 1e-12);
    check("plus right sp rod - trfv 6        ", tau_bq1()(5), tau_b1()(5), 1e-12);
    check("plus right sp rod - screw 1       ", tau_bq1()(0), tau_bS1()(0), 1e-12);
    check("plus right sp rod - screw 2       ", tau_bq1()(1), tau_bS1()(1), 1e-12);
    check("plus right sp rod - screw 3       ", tau_bq1()(2), tau_bS1()(2), 1e-12);
    check("plus right sp rod - screw 4       ", tau_bq1()(3), tau_bS1()(3), 1e-12);
    check("plus right sp rod - screw 5       ", tau_bq1()(4), tau_bS1()(4), 1e-12);
    check("plus right sp rod - screw 6       ", tau_bq1()(5), tau_bS1()(5), 1e-12);
    check("plus right sp rod - dual 1        ", tau_bq1()(0), tau_bZ1()(0), 1e-12);
    check("plus right sp rod - dual 2        ", tau_bq1()(1), tau_bZ1()(1), 1e-12);
    check("plus right sp rod - dual 3        ", tau_bq1()(2), tau_bZ1()(2), 1e-12);
    check("plus right sp rod - dual 4        ", tau_bq1()(3), tau_bZ1()(3), 1e-12);
    check("plus right sp rod - dual 5        ", tau_bq1()(4), tau_bZ1()(4), 1e-12);
    check("plus right sp rod - dual 6        ", tau_bq1()(5), tau_bZ1()(5), 1e-12);

    check("plus left sp rod - sp dcm 1      ", tau_bq2()(0), tau_bR2()(0), 1e-12);
    check("plus left sp rod - sp dcm 2      ", tau_bq2()(1), tau_bR2()(1), 1e-12);
    check("plus left sp rod - sp dcm 3      ", tau_bq2()(2), tau_bR2()(2), 1e-12);
    check("plus left sp rod - sp dcm 4      ", tau_bq2()(3), tau_bR2()(3), 1e-12);
    check("plus left sp rod - sp dcm 5      ", tau_bq2()(4), tau_bR2()(4), 1e-12);
    check("plus left sp rod - sp dcm 6      ", tau_bq2()(5), tau_bR2()(5), 1e-12);
    check("plus left sp rod - homog 1       ", tau_bq2()(0), tau_bM2()(0), 1e-12);
    check("plus left sp rod - homog 2       ", tau_bq2()(1), tau_bM2()(1), 1e-12);
    check("plus left sp rod - homog 3       ", tau_bq2()(2), tau_bM2()(2), 1e-12);
    check("plus left sp rod - homog 4       ", tau_bq2()(3), tau_bM2()(3), 1e-12);
    check("plus left sp rod - homog 5       ", tau_bq2()(4), tau_bM2()(4), 1e-12);
    check("plus left sp rod - homog 6       ", tau_bq2()(5), tau_bM2()(5), 1e-12);
    check("plus left sp rod - trfv 1        ", tau_bq2()(0), tau_b2()(0), 1e-12);
    check("plus left sp rod - trfv 2        ", tau_bq2()(1), tau_b2()(1), 1e-12);
    check("plus left sp rod - trfv 3        ", tau_bq2()(2), tau_b2()(2), 1e-12);
    check("plus left sp rod - trfv 4        ", tau_bq2()(3), tau_b2()(3), 1e-12);
    check("plus left sp rod - trfv 5        ", tau_bq2()(4), tau_b2()(4), 1e-12);
    check("plus left sp rod - trfv 6        ", tau_bq2()(5), tau_b2()(5), 1e-12);
    check("plus left sp rod - screw 1       ", tau_bq2()(0), tau_bS2()(0), 1e-12);
    check("plus left sp rod - screw 2       ", tau_bq2()(1), tau_bS2()(1), 1e-12);
    check("plus left sp rod - screw 3       ", tau_bq2()(2), tau_bS2()(2), 1e-12);
    check("plus left sp rod - screw 4       ", tau_bq2()(3), tau_bS2()(3), 1e-12);
    check("plus left sp rod - screw 5       ", tau_bq2()(4), tau_bS2()(4), 1e-12);
    check("plus left sp rod - screw 6       ", tau_bq2()(5), tau_bS2()(5), 1e-12);
    check("plus left sp rod - dual 1        ", tau_bq2()(0), tau_bZ2()(0), 1e-12);
    check("plus left sp rod - dual 2        ", tau_bq2()(1), tau_bZ2()(1), 1e-12);
    check("plus left sp rod - dual 3        ", tau_bq2()(2), tau_bZ2()(2), 1e-12);
    check("plus left sp rod - dual 4        ", tau_bq2()(3), tau_bZ2()(3), 1e-12);
    check("plus left sp rod - dual 5        ", tau_bq2()(4), tau_bZ2()(4), 1e-12);
    check("plus left sp rod - dual 6        ", tau_bq2()(5), tau_bZ2()(5), 1e-12);

    ang::speu_rodrigues Xgq_b1  = gq_a.plus_right(Delta_S);
    ang::speu_dcm       XgR_b1  = gR_a.plus_right(Delta_S);
    ang::homogeneous    XM_b1   = M_a.plus_right(Delta_S);
    ang::trfv           Xtau_b1 = tau_a.plus_right(Delta_S);
    ang::screw          XS_b1   = S_a.plus_right(Delta_S);
    ang::dual           XZ_b1   = Z_a.plus_right(Delta_S);

    ang::trfv Xtau_bq1(Xgq_b1);
    ang::trfv Xtau_bR1(XgR_b1);
    ang::trfv Xtau_bM1(XM_b1);
    ang::trfv Xtau_bS1(XS_b1);
    ang::trfv Xtau_bZ1(XZ_b1);

    ang::speu_rodrigues Xgq_b2  = gq_a.plus_left(Delta_S);
    ang::speu_dcm       XgR_b2  = gR_a.plus_left(Delta_S);
    ang::homogeneous    XM_b2   = M_a.plus_left(Delta_S);
    ang::trfv           Xtau_b2 = tau_a.plus_left(Delta_S);
    ang::screw          XS_b2   = S_a.plus_left(Delta_S);
    ang::dual           XZ_b2   = Z_a.plus_left(Delta_S);

    ang::trfv Xtau_bq2(Xgq_b2);
    ang::trfv Xtau_bR2(XgR_b2);
    ang::trfv Xtau_bM2(XM_b2);
    ang::trfv Xtau_bS2(XS_b2);
    ang::trfv Xtau_bZ2(XZ_b2);

    check("plus right sp rod - sp dcm 1      ", Xtau_bq1()(0), Xtau_bR1()(0), 1e-12);
    check("plus right sp rod - sp dcm 2      ", Xtau_bq1()(1), Xtau_bR1()(1), 1e-12);
    check("plus right sp rod - sp dcm 3      ", Xtau_bq1()(2), Xtau_bR1()(2), 1e-12);
    check("plus right sp rod - sp dcm 4      ", Xtau_bq1()(3), Xtau_bR1()(3), 1e-12);
    check("plus right sp rod - sp dcm 5      ", Xtau_bq1()(4), Xtau_bR1()(4), 1e-12);
    check("plus right sp rod - sp dcm 6      ", Xtau_bq1()(5), Xtau_bR1()(5), 1e-12);
    check("plus right sp rod - homog 1       ", Xtau_bq1()(0), Xtau_bM1()(0), 1e-12);
    check("plus right sp rod - homog 2       ", Xtau_bq1()(1), Xtau_bM1()(1), 1e-12);
    check("plus right sp rod - homog 3       ", Xtau_bq1()(2), Xtau_bM1()(2), 1e-12);
    check("plus right sp rod - homog 4       ", Xtau_bq1()(3), Xtau_bM1()(3), 1e-12);
    check("plus right sp rod - homog 5       ", Xtau_bq1()(4), Xtau_bM1()(4), 1e-12);
    check("plus right sp rod - homog 6       ", Xtau_bq1()(5), Xtau_bM1()(5), 1e-12);
    check("plus right sp rod - trfv 1        ", Xtau_bq1()(0), Xtau_b1()(0), 1e-12);
    check("plus right sp rod - trfv 2        ", Xtau_bq1()(1), Xtau_b1()(1), 1e-12);
    check("plus right sp rod - trfv 3        ", Xtau_bq1()(2), Xtau_b1()(2), 1e-12);
    check("plus right sp rod - trfv 4        ", Xtau_bq1()(3), Xtau_b1()(3), 1e-12);
    check("plus right sp rod - trfv 5        ", Xtau_bq1()(4), Xtau_b1()(4), 1e-12);
    check("plus right sp rod - trfv 6        ", Xtau_bq1()(5), Xtau_b1()(5), 1e-12);
    check("plus right sp rod - screw 1       ", Xtau_bq1()(0), Xtau_bS1()(0), 1e-12);
    check("plus right sp rod - screw 2       ", Xtau_bq1()(1), Xtau_bS1()(1), 1e-12);
    check("plus right sp rod - screw 3       ", Xtau_bq1()(2), Xtau_bS1()(2), 1e-12);
    check("plus right sp rod - screw 4       ", Xtau_bq1()(3), Xtau_bS1()(3), 1e-12);
    check("plus right sp rod - screw 5       ", Xtau_bq1()(4), Xtau_bS1()(4), 1e-12);
    check("plus right sp rod - screw 6       ", Xtau_bq1()(5), Xtau_bS1()(5), 1e-12);
    check("plus right sp rod - dual 1        ", Xtau_bq1()(0), Xtau_bZ1()(0), 1e-12);
    check("plus right sp rod - dual 2        ", Xtau_bq1()(1), Xtau_bZ1()(1), 1e-12);
    check("plus right sp rod - dual 3        ", Xtau_bq1()(2), Xtau_bZ1()(2), 1e-12);
    check("plus right sp rod - dual 4        ", Xtau_bq1()(3), Xtau_bZ1()(3), 1e-12);
    check("plus right sp rod - dual 5        ", Xtau_bq1()(4), Xtau_bZ1()(4), 1e-12);
    check("plus right sp rod - dual 6        ", Xtau_bq1()(5), Xtau_bZ1()(5), 1e-12);

    check("plus left sp rod - sp dcm 1      ", Xtau_bq2()(0), Xtau_bR2()(0), 1e-12);
    check("plus left sp rod - sp dcm 2      ", Xtau_bq2()(1), Xtau_bR2()(1), 1e-12);
    check("plus left sp rod - sp dcm 3      ", Xtau_bq2()(2), Xtau_bR2()(2), 1e-12);
    check("plus left sp rod - sp dcm 4      ", Xtau_bq2()(3), Xtau_bR2()(3), 1e-12);
    check("plus left sp rod - sp dcm 5      ", Xtau_bq2()(4), Xtau_bR2()(4), 1e-12);
    check("plus left sp rod - sp dcm 6      ", Xtau_bq2()(5), Xtau_bR2()(5), 1e-12);
    check("plus left sp rod - homog 1       ", Xtau_bq2()(0), Xtau_bM2()(0), 1e-12);
    check("plus left sp rod - homog 2       ", Xtau_bq2()(1), Xtau_bM2()(1), 1e-12);
    check("plus left sp rod - homog 3       ", Xtau_bq2()(2), Xtau_bM2()(2), 1e-12);
    check("plus left sp rod - homog 4       ", Xtau_bq2()(3), Xtau_bM2()(3), 1e-12);
    check("plus left sp rod - homog 5       ", Xtau_bq2()(4), Xtau_bM2()(4), 1e-12);
    check("plus left sp rod - homog 6       ", Xtau_bq2()(5), Xtau_bM2()(5), 1e-12);
    check("plus left sp rod - trfv 1        ", Xtau_bq2()(0), Xtau_b2()(0), 1e-12);
    check("plus left sp rod - trfv 2        ", Xtau_bq2()(1), Xtau_b2()(1), 1e-12);
    check("plus left sp rod - trfv 3        ", Xtau_bq2()(2), Xtau_b2()(2), 1e-12);
    check("plus left sp rod - trfv 4        ", Xtau_bq2()(3), Xtau_b2()(3), 1e-12);
    check("plus left sp rod - trfv 5        ", Xtau_bq2()(4), Xtau_b2()(4), 1e-12);
    check("plus left sp rod - trfv 6        ", Xtau_bq2()(5), Xtau_b2()(5), 1e-12);
    check("plus left sp rod - screw 1       ", Xtau_bq2()(0), Xtau_bS2()(0), 1e-12);
    check("plus left sp rod - screw 2       ", Xtau_bq2()(1), Xtau_bS2()(1), 1e-12);
    check("plus left sp rod - screw 3       ", Xtau_bq2()(2), Xtau_bS2()(2), 1e-12);
    check("plus left sp rod - screw 4       ", Xtau_bq2()(3), Xtau_bS2()(3), 1e-12);
    check("plus left sp rod - screw 5       ", Xtau_bq2()(4), Xtau_bS2()(4), 1e-12);
    check("plus left sp rod - screw 6       ", Xtau_bq2()(5), Xtau_bS2()(5), 1e-12);
    check("plus left sp rod - dual 1        ", Xtau_bq2()(0), Xtau_bZ2()(0), 1e-12);
    check("plus left sp rod - dual 2        ", Xtau_bq2()(1), Xtau_bZ2()(1), 1e-12);
    check("plus left sp rod - dual 3        ", Xtau_bq2()(2), Xtau_bZ2()(2), 1e-12);
    check("plus left sp rod - dual 4        ", Xtau_bq2()(3), Xtau_bZ2()(3), 1e-12);
    check("plus left sp rod - dual 5        ", Xtau_bq2()(4), Xtau_bZ2()(4), 1e-12);
    check("plus left sp rod - dual 6        ", Xtau_bq2()(5), Xtau_bZ2()(5), 1e-12);

    ang::trfv tau_cq1   = gq_b1.minus_right_trfv(gq_a);
    ang::trfv tau_cR1   = gR_b1.minus_right_trfv(gR_a);
    ang::trfv tau_cM1   = M_b1.minus_right_trfv(M_a);
    ang::trfv tau_ctau1 = tau_b1.minus_right_trfv(tau_a);
    ang::trfv tau_cS1   = S_b1.minus_right_trfv(S_a);
    ang::trfv tau_cZ1   = Z_b1.minus_right_trfv(Z_a);

    ang::trfv tau_cq2   = gq_b2.minus_left_trfv(gq_a);
    ang::trfv tau_cR2   = gR_b2.minus_left_trfv(gR_a);
    ang::trfv tau_cM2   = M_b2.minus_left_trfv(M_a);
    ang::trfv tau_ctau2 = tau_b2.minus_left_trfv(tau_a);
    ang::trfv tau_cS2   = S_b2.minus_left_trfv(S_a);
    ang::trfv tau_cZ2   = Z_b2.minus_left_trfv(Z_a);

    check("minus right sp rod 1              ", tau_cq1()(0), Delta_tau()(0), 1e-12);
    check("minus right sp rod 2              ", tau_cq1()(1), Delta_tau()(1), 1e-12);
    check("minus right sp rod 3              ", tau_cq1()(2), Delta_tau()(2), 1e-12);
    check("minus right sp rod 4              ", tau_cq1()(3), Delta_tau()(3), 1e-12);
    check("minus right sp rod 5              ", tau_cq1()(4), Delta_tau()(4), 1e-12);
    check("minus right sp rod 6              ", tau_cq1()(5), Delta_tau()(5), 1e-12);
    check("minus right sp dcm 1              ", tau_cR1()(0), Delta_tau()(0), 1e-12);
    check("minus right sp dcm 2              ", tau_cR1()(1), Delta_tau()(1), 1e-12);
    check("minus right sp dcm 3              ", tau_cR1()(2), Delta_tau()(2), 1e-12);
    check("minus right sp dcm 4              ", tau_cR1()(3), Delta_tau()(3), 1e-12);
    check("minus right sp dcm 5              ", tau_cR1()(4), Delta_tau()(4), 1e-12);
    check("minus right sp dcm 6              ", tau_cR1()(5), Delta_tau()(5), 1e-12);
    check("minus right homog 1               ", tau_cM1()(0), Delta_tau()(0), 1e-12);
    check("minus right homog 2               ", tau_cM1()(1), Delta_tau()(1), 1e-12);
    check("minus right homog 3               ", tau_cM1()(2), Delta_tau()(2), 1e-12);
    check("minus right homog 4               ", tau_cM1()(3), Delta_tau()(3), 1e-12);
    check("minus right homog 5               ", tau_cM1()(4), Delta_tau()(4), 1e-12);
    check("minus right homog 6               ", tau_cM1()(5), Delta_tau()(5), 1e-12);
    check("minus right trfv 1                ", tau_ctau1()(0), Delta_tau()(0), 1e-12);
    check("minus right trfv 2                ", tau_ctau1()(1), Delta_tau()(1), 1e-12);
    check("minus right trfv 3                ", tau_ctau1()(2), Delta_tau()(2), 1e-12);
    check("minus right trfv 4                ", tau_ctau1()(3), Delta_tau()(3), 1e-12);
    check("minus right trfv 5                ", tau_ctau1()(4), Delta_tau()(4), 1e-12);
    check("minus right trfv 6                ", tau_ctau1()(5), Delta_tau()(5), 1e-12);
    check("minus right screw 1               ", tau_cS1()(0), Delta_tau()(0), 1e-12);
    check("minus right screw 2               ", tau_cS1()(1), Delta_tau()(1), 1e-12);
    check("minus right screw 3               ", tau_cS1()(2), Delta_tau()(2), 1e-12);
    check("minus right screw 4               ", tau_cS1()(3), Delta_tau()(3), 1e-12);
    check("minus right screw 5               ", tau_cS1()(4), Delta_tau()(4), 1e-12);
    check("minus right screw 6               ", tau_cS1()(5), Delta_tau()(5), 1e-12);
    check("minus right dual 1                ", tau_cZ1()(0), Delta_tau()(0), 1e-12);
    check("minus right dual 2                ", tau_cZ1()(1), Delta_tau()(1), 1e-12);
    check("minus right dual 3                ", tau_cZ1()(2), Delta_tau()(2), 1e-12);
    check("minus right dual 4                ", tau_cZ1()(3), Delta_tau()(3), 1e-12);
    check("minus right dual 5                ", tau_cZ1()(4), Delta_tau()(4), 1e-12);
    check("minus right dual 6                ", tau_cZ1()(5), Delta_tau()(5), 1e-12);

    check("minus left sp rod 1              ", tau_cq2()(0), Delta_tau()(0), 1e-12);
    check("minus left sp rod 2              ", tau_cq2()(1), Delta_tau()(1), 1e-12);
    check("minus left sp rod 3              ", tau_cq2()(2), Delta_tau()(2), 1e-12);
    check("minus left sp rod 4              ", tau_cq2()(3), Delta_tau()(3), 1e-12);
    check("minus left sp rod 5              ", tau_cq2()(4), Delta_tau()(4), 1e-12);
    check("minus left sp rod 6              ", tau_cq2()(5), Delta_tau()(5), 1e-12);
    check("minus left sp dcm 1              ", tau_cR2()(0), Delta_tau()(0), 1e-12);
    check("minus left sp dcm 2              ", tau_cR2()(1), Delta_tau()(1), 1e-12);
    check("minus left sp dcm 3              ", tau_cR2()(2), Delta_tau()(2), 1e-12);
    check("minus left sp dcm 4              ", tau_cR2()(3), Delta_tau()(3), 1e-12);
    check("minus left sp dcm 5              ", tau_cR2()(4), Delta_tau()(4), 1e-12);
    check("minus left sp dcm 6              ", tau_cR2()(5), Delta_tau()(5), 1e-12);
    check("minus left homog 1               ", tau_cM2()(0), Delta_tau()(0), 1e-12);
    check("minus left homog 2               ", tau_cM2()(1), Delta_tau()(1), 1e-12);
    check("minus left homog 3               ", tau_cM2()(2), Delta_tau()(2), 1e-12);
    check("minus left homog 4               ", tau_cM2()(3), Delta_tau()(3), 1e-12);
    check("minus left homog 5               ", tau_cM2()(4), Delta_tau()(4), 1e-12);
    check("minus left homog 6               ", tau_cM2()(5), Delta_tau()(5), 1e-12);
    check("minus left trfv 1                ", tau_ctau2()(0), Delta_tau()(0), 1e-12);
    check("minus left trfv 2                ", tau_ctau2()(1), Delta_tau()(1), 1e-12);
    check("minus left trfv 3                ", tau_ctau2()(2), Delta_tau()(2), 1e-12);
    check("minus left trfv 4                ", tau_ctau2()(3), Delta_tau()(3), 1e-12);
    check("minus left trfv 5                ", tau_ctau2()(4), Delta_tau()(4), 1e-12);
    check("minus left trfv 6                ", tau_ctau2()(5), Delta_tau()(5), 1e-12);
    check("minus left screw 1               ", tau_cS2()(0), Delta_tau()(0), 1e-12);
    check("minus left screw 2               ", tau_cS2()(1), Delta_tau()(1), 1e-12);
    check("minus left screw 3               ", tau_cS2()(2), Delta_tau()(2), 1e-12);
    check("minus left screw 4               ", tau_cS2()(3), Delta_tau()(3), 1e-12);
    check("minus left screw 5               ", tau_cS2()(4), Delta_tau()(4), 1e-12);
    check("minus left screw 6               ", tau_cS2()(5), Delta_tau()(5), 1e-12);
    check("minus left dual 1                ", tau_cZ2()(0), Delta_tau()(0), 1e-12);
    check("minus left dual 2                ", tau_cZ2()(1), Delta_tau()(1), 1e-12);
    check("minus left dual 3                ", tau_cZ2()(2), Delta_tau()(2), 1e-12);
    check("minus left dual 4                ", tau_cZ2()(3), Delta_tau()(3), 1e-12);
    check("minus left dual 5                ", tau_cZ2()(4), Delta_tau()(4), 1e-12);
    check("minus left dual 6                ", tau_cZ2()(5), Delta_tau()(5), 1e-12);

    ang::screw XS_cq1   = Xgq_b1.minus_right_screw(gq_a);
    ang::screw XS_cR1   = XgR_b1.minus_right_screw(gR_a);
    ang::screw XS_cM1   = XM_b1.minus_right_screw(M_a);
    ang::screw XS_ctau1 = Xtau_b1.minus_right_screw(tau_a);
    ang::screw XS_cS1   = XS_b1.minus_right_screw(S_a);
    ang::screw XS_cZ1   = XZ_b1.minus_right_screw(Z_a);

    ang::trfv Xtau_cq1(XS_cq1);
    ang::trfv Xtau_cR1(XS_cR1);
    ang::trfv Xtau_cM1(XS_cM1);
    ang::trfv Xtau_ctau1(XS_ctau1);
    ang::trfv Xtau_cS1(XS_cS1);
    ang::trfv Xtau_cZ1(XS_cZ1);

    ang::screw XS_cq2   = Xgq_b2.minus_left_screw(gq_a);
    ang::screw XS_cR2   = XgR_b2.minus_left_screw(gR_a);
    ang::screw XS_cM2   = XM_b2.minus_left_screw(M_a);
    ang::screw XS_ctau2 = tau_b2.minus_left_screw(tau_a);
    ang::screw XS_cS2   = S_b2.minus_left_screw(S_a);
    ang::screw XS_cZ2   = Z_b2.minus_left_screw(Z_a);

    ang::trfv Xtau_cq2(XS_cq2);
    ang::trfv Xtau_cR2(XS_cR2);
    ang::trfv Xtau_cM2(XS_cM2);
    ang::trfv Xtau_ctau2(XS_ctau2);
    ang::trfv Xtau_cS2(XS_cS2);
    ang::trfv Xtau_cZ2(XS_cZ2);

    check("minus right sp rod 1              ", Xtau_cq1()(0), Delta_tau()(0), 1e-12);
    check("minus right sp rod 2              ", Xtau_cq1()(1), Delta_tau()(1), 1e-12);
    check("minus right sp rod 3              ", Xtau_cq1()(2), Delta_tau()(2), 1e-12);
    check("minus right sp rod 4              ", Xtau_cq1()(3), Delta_tau()(3), 1e-12);
    check("minus right sp rod 5              ", Xtau_cq1()(4), Delta_tau()(4), 1e-12);
    check("minus right sp rod 6              ", Xtau_cq1()(5), Delta_tau()(5), 1e-12);
    check("minus right sp dcm 1              ", Xtau_cR1()(0), Delta_tau()(0), 1e-12);
    check("minus right sp dcm 2              ", Xtau_cR1()(1), Delta_tau()(1), 1e-12);
    check("minus right sp dcm 3              ", Xtau_cR1()(2), Delta_tau()(2), 1e-12);
    check("minus right sp dcm 4              ", Xtau_cR1()(3), Delta_tau()(3), 1e-12);
    check("minus right sp dcm 5              ", Xtau_cR1()(4), Delta_tau()(4), 1e-12);
    check("minus right sp dcm 6              ", Xtau_cR1()(5), Delta_tau()(5), 1e-12);
    check("minus right homog 1               ", Xtau_cM1()(0), Delta_tau()(0), 1e-12);
    check("minus right homog 2               ", Xtau_cM1()(1), Delta_tau()(1), 1e-12);
    check("minus right homog 3               ", Xtau_cM1()(2), Delta_tau()(2), 1e-12);
    check("minus right homog 4               ", Xtau_cM1()(3), Delta_tau()(3), 1e-12);
    check("minus right homog 5               ", Xtau_cM1()(4), Delta_tau()(4), 1e-12);
    check("minus right homog 6               ", Xtau_cM1()(5), Delta_tau()(5), 1e-12);
    check("minus right trfv 1                ", Xtau_ctau1()(0), Delta_tau()(0), 1e-12);
    check("minus right trfv 2                ", Xtau_ctau1()(1), Delta_tau()(1), 1e-12);
    check("minus right trfv 3                ", Xtau_ctau1()(2), Delta_tau()(2), 1e-12);
    check("minus right trfv 4                ", Xtau_ctau1()(3), Delta_tau()(3), 1e-12);
    check("minus right trfv 5                ", Xtau_ctau1()(4), Delta_tau()(4), 1e-12);
    check("minus right trfv 6                ", Xtau_ctau1()(5), Delta_tau()(5), 1e-12);
    check("minus right screw 1               ", Xtau_cS1()(0), Delta_tau()(0), 1e-12);
    check("minus right screw 2               ", Xtau_cS1()(1), Delta_tau()(1), 1e-12);
    check("minus right screw 3               ", Xtau_cS1()(2), Delta_tau()(2), 1e-12);
    check("minus right screw 4               ", Xtau_cS1()(3), Delta_tau()(3), 1e-12);
    check("minus right screw 5               ", Xtau_cS1()(4), Delta_tau()(4), 1e-12);
    check("minus right screw 6               ", Xtau_cS1()(5), Delta_tau()(5), 1e-12);
    check("minus right dual 1                ", Xtau_cZ1()(0), Delta_tau()(0), 1e-12);
    check("minus right dual 2                ", Xtau_cZ1()(1), Delta_tau()(1), 1e-12);
    check("minus right dual 3                ", Xtau_cZ1()(2), Delta_tau()(2), 1e-12);
    check("minus right dual 4                ", Xtau_cZ1()(3), Delta_tau()(3), 1e-12);
    check("minus right dual 5                ", Xtau_cZ1()(4), Delta_tau()(4), 1e-12);
    check("minus right dual 6                ", Xtau_cZ1()(5), Delta_tau()(5), 1e-12);

    check("minus left sp rod 1              ", Xtau_cq2()(0), Delta_tau()(0), 1e-12);
    check("minus left sp rod 2              ", Xtau_cq2()(1), Delta_tau()(1), 1e-12);
    check("minus left sp rod 3              ", Xtau_cq2()(2), Delta_tau()(2), 1e-12);
    check("minus left sp rod 4              ", Xtau_cq2()(3), Delta_tau()(3), 1e-12);
    check("minus left sp rod 5              ", Xtau_cq2()(4), Delta_tau()(4), 1e-12);
    check("minus left sp rod 6              ", Xtau_cq2()(5), Delta_tau()(5), 1e-12);
    check("minus left sp dcm 1              ", Xtau_cR2()(0), Delta_tau()(0), 1e-12);
    check("minus left sp dcm 2              ", Xtau_cR2()(1), Delta_tau()(1), 1e-12);
    check("minus left sp dcm 3              ", Xtau_cR2()(2), Delta_tau()(2), 1e-12);
    check("minus left sp dcm 4              ", Xtau_cR2()(3), Delta_tau()(3), 1e-12);
    check("minus left sp dcm 5              ", Xtau_cR2()(4), Delta_tau()(4), 1e-12);
    check("minus left sp dcm 6              ", Xtau_cR2()(5), Delta_tau()(5), 1e-12);
    check("minus left homog 1               ", Xtau_cM2()(0), Delta_tau()(0), 1e-12);
    check("minus left homog 2               ", Xtau_cM2()(1), Delta_tau()(1), 1e-12);
    check("minus left homog 3               ", Xtau_cM2()(2), Delta_tau()(2), 1e-12);
    check("minus left homog 4               ", Xtau_cM2()(3), Delta_tau()(3), 1e-12);
    check("minus left homog 5               ", Xtau_cM2()(4), Delta_tau()(4), 1e-12);
    check("minus left homog 6               ", Xtau_cM2()(5), Delta_tau()(5), 1e-12);
    check("minus left trfv 1                ", Xtau_ctau2()(0), Delta_tau()(0), 1e-12);
    check("minus left trfv 2                ", Xtau_ctau2()(1), Delta_tau()(1), 1e-12);
    check("minus left trfv 3                ", Xtau_ctau2()(2), Delta_tau()(2), 1e-12);
    check("minus left trfv 4                ", Xtau_ctau2()(3), Delta_tau()(3), 1e-12);
    check("minus left trfv 5                ", Xtau_ctau2()(4), Delta_tau()(4), 1e-12);
    check("minus left trfv 6                ", Xtau_ctau2()(5), Delta_tau()(5), 1e-12);
    check("minus left screw 1               ", Xtau_cS2()(0), Delta_tau()(0), 1e-12);
    check("minus left screw 2               ", Xtau_cS2()(1), Delta_tau()(1), 1e-12);
    check("minus left screw 3               ", Xtau_cS2()(2), Delta_tau()(2), 1e-12);
    check("minus left screw 4               ", Xtau_cS2()(3), Delta_tau()(3), 1e-12);
    check("minus left screw 5               ", Xtau_cS2()(4), Delta_tau()(4), 1e-12);
    check("minus left screw 6               ", Xtau_cS2()(5), Delta_tau()(5), 1e-12);
    check("minus left dual 1                ", Xtau_cZ2()(0), Delta_tau()(0), 1e-12);
    check("minus left dual 2                ", Xtau_cZ2()(1), Delta_tau()(1), 1e-12);
    check("minus left dual 3                ", Xtau_cZ2()(2), Delta_tau()(2), 1e-12);
    check("minus left dual 4                ", Xtau_cZ2()(3), Delta_tau()(3), 1e-12);
    check("minus left dual 5                ", Xtau_cZ2()(4), Delta_tau()(4), 1e-12);
    check("minus left dual 6                ", Xtau_cZ2()(5), Delta_tau()(5), 1e-12);

} // closes test_plus_minus

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_adjoint() {
    ang::se3_tangent xi_ebb_0(0.24, 2.37, -5.02, 0.51, -0.38, -0.14);
    ang::se3_tangent_homo xi_ebb_0_homo(xi_ebb_0);
    ang::se3_tangent_dual xi_ebb_0_dual(xi_ebb_0);

    ang::rotv r_xb(-0.11, -0.42, 0.69);
    Eigen::Vector3d x_ebe(0.84, -0.09, 0.39);
    ang::trfv           tau_xb(r_xb, x_ebe);
    ang::speu_rodrigues Gq_xb(tau_xb);
    ang::speu_dcm       Gr_xb(tau_xb);
    ang::homogeneous    M_xb(tau_xb);
    ang::dual           Z_xb(tau_xb);
    ang::screw          S_xb(tau_xb);

    ang::se3_tangent xi_ebe_01;
    xi_ebe_01.set_w(ang::so3_tangent(tau_xb ^ xi_ebb_0.get_w()()));
    xi_ebe_01.set_vi((tau_xb ^ xi_ebb_0.get_vi()) + ang::tools::skew3(tau_xb.get_T()) * xi_ebe_01.get_w()()); // do not remove parenthesis
    ang::se3_tangent xi_ebe_02 = M_xb | xi_ebb_0;
    ang::se3_tangent xi_ebe_03 = Gr_xb | xi_ebb_0;
    ang::se3_tangent xi_ebe_04 = Gq_xb | xi_ebb_0;
    ang::se3_tangent xi_ebe_05 = tau_xb | xi_ebb_0;
    ang::se3_tangent xi_ebe_06 = Z_xb | xi_ebb_0;
    ang::se3_tangent xi_ebe_07 = ang::se3_tangent::wedge(M_xb | xi_ebb_0_homo);
    ang::se3_tangent xi_ebe_08 = ang::se3_tangent::wedge(Gr_xb | xi_ebb_0_homo);
    ang::se3_tangent xi_ebe_09 = ang::se3_tangent::wedge(Gq_xb | xi_ebb_0_homo);
    ang::se3_tangent xi_ebe_11 = ang::se3_tangent::wedge(Z_xb | xi_ebb_0_dual);

    check("adjoint 21        ", xi_ebe_01()(0), xi_ebe_02()(0), 1e-12);
    check("adjoint 22        ", xi_ebe_01()(1), xi_ebe_02()(1), 1e-12);
    check("adjoint 23        ", xi_ebe_01()(2), xi_ebe_02()(2), 1e-12);
    check("adjoint 24        ", xi_ebe_01()(3), xi_ebe_02()(3), 1e-12);
    check("adjoint 25        ", xi_ebe_01()(4), xi_ebe_02()(4), 1e-12);
    check("adjoint 26        ", xi_ebe_01()(5), xi_ebe_02()(5), 1e-12);

    check("adjoint 31        ", xi_ebe_01()(0), xi_ebe_03()(0), 1e-12);
    check("adjoint 32        ", xi_ebe_01()(1), xi_ebe_03()(1), 1e-12);
    check("adjoint 33        ", xi_ebe_01()(2), xi_ebe_03()(2), 1e-12);
    check("adjoint 34        ", xi_ebe_01()(3), xi_ebe_03()(3), 1e-12);
    check("adjoint 35        ", xi_ebe_01()(4), xi_ebe_03()(4), 1e-12);
    check("adjoint 36        ", xi_ebe_01()(5), xi_ebe_03()(5), 1e-12);

    check("adjoint 41        ", xi_ebe_01()(0), xi_ebe_04()(0), 1e-12);
    check("adjoint 42        ", xi_ebe_01()(1), xi_ebe_04()(1), 1e-12);
    check("adjoint 43        ", xi_ebe_01()(2), xi_ebe_04()(2), 1e-12);
    check("adjoint 44        ", xi_ebe_01()(3), xi_ebe_04()(3), 1e-12);
    check("adjoint 45        ", xi_ebe_01()(4), xi_ebe_04()(4), 1e-12);
    check("adjoint 46        ", xi_ebe_01()(5), xi_ebe_04()(5), 1e-12);

    check("adjoint 51        ", xi_ebe_01()(0), xi_ebe_05()(0), 1e-12);
    check("adjoint 52        ", xi_ebe_01()(1), xi_ebe_05()(1), 1e-12);
    check("adjoint 53        ", xi_ebe_01()(2), xi_ebe_05()(2), 1e-12);
    check("adjoint 54        ", xi_ebe_01()(3), xi_ebe_05()(3), 1e-12);
    check("adjoint 55        ", xi_ebe_01()(4), xi_ebe_05()(4), 1e-12);
    check("adjoint 56        ", xi_ebe_01()(5), xi_ebe_05()(5), 1e-12);

    check("adjoint 61        ", xi_ebe_01()(0), xi_ebe_06()(0), 1e-12);
    check("adjoint 62        ", xi_ebe_01()(1), xi_ebe_06()(1), 1e-12);
    check("adjoint 63        ", xi_ebe_01()(2), xi_ebe_06()(2), 1e-12);
    check("adjoint 64        ", xi_ebe_01()(3), xi_ebe_06()(3), 1e-12);
    check("adjoint 65        ", xi_ebe_01()(4), xi_ebe_06()(4), 1e-12);
    check("adjoint 66        ", xi_ebe_01()(5), xi_ebe_06()(5), 1e-12);

    check("adjoint 71        ", xi_ebe_01()(0), xi_ebe_07()(0), 1e-12);
    check("adjoint 72        ", xi_ebe_01()(1), xi_ebe_07()(1), 1e-12);
    check("adjoint 73        ", xi_ebe_01()(2), xi_ebe_07()(2), 1e-12);
    check("adjoint 74        ", xi_ebe_01()(3), xi_ebe_07()(3), 1e-12);
    check("adjoint 75        ", xi_ebe_01()(4), xi_ebe_07()(4), 1e-12);
    check("adjoint 76        ", xi_ebe_01()(5), xi_ebe_07()(5), 1e-12);

    check("adjoint 81        ", xi_ebe_01()(0), xi_ebe_08()(0), 1e-12);
    check("adjoint 82        ", xi_ebe_01()(1), xi_ebe_08()(1), 1e-12);
    check("adjoint 83        ", xi_ebe_01()(2), xi_ebe_08()(2), 1e-12);
    check("adjoint 84        ", xi_ebe_01()(3), xi_ebe_08()(3), 1e-12);
    check("adjoint 85        ", xi_ebe_01()(4), xi_ebe_08()(4), 1e-12);
    check("adjoint 86        ", xi_ebe_01()(5), xi_ebe_08()(5), 1e-12);

    check("adjoint 91        ", xi_ebe_01()(0), xi_ebe_09()(0), 1e-12);
    check("adjoint 92        ", xi_ebe_01()(1), xi_ebe_09()(1), 1e-12);
    check("adjoint 93        ", xi_ebe_01()(2), xi_ebe_09()(2), 1e-12);
    check("adjoint 94        ", xi_ebe_01()(3), xi_ebe_09()(3), 1e-12);
    check("adjoint 95        ", xi_ebe_01()(4), xi_ebe_09()(4), 1e-12);
    check("adjoint 96        ", xi_ebe_01()(5), xi_ebe_09()(5), 1e-12);

    check("adjoint 111       ", xi_ebe_01()(0), xi_ebe_11()(0), 1e-12);
    check("adjoint 112       ", xi_ebe_01()(1), xi_ebe_11()(1), 1e-12);
    check("adjoint 113       ", xi_ebe_01()(2), xi_ebe_11()(2), 1e-12);
    check("adjoint 114       ", xi_ebe_01()(3), xi_ebe_11()(3), 1e-12);
    check("adjoint 115       ", xi_ebe_01()(4), xi_ebe_11()(4), 1e-12);
    check("adjoint 116       ", xi_ebe_01()(5), xi_ebe_11()(5), 1e-12);

    ang::se3_tangent xi_ebb_01;
    xi_ebb_01.set_w(ang::so3_tangent(tau_xb & xi_ebe_01.get_w()()));
    xi_ebb_01.set_vi((tau_xb & xi_ebe_01.get_vi()) - (tau_xb & (ang::tools::skew3(tau_xb.get_T()) * xi_ebe_01.get_w()()))); // do not remove parenthesis
    ang::se3_tangent xi_ebb_02 = M_xb % xi_ebe_02;
    ang::se3_tangent xi_ebb_03 = Gr_xb % xi_ebe_03;
    ang::se3_tangent xi_ebb_04 = Gq_xb % xi_ebe_04;
    ang::se3_tangent xi_ebb_05 = tau_xb % xi_ebe_05;
    ang::se3_tangent xi_ebb_06 = Z_xb % xi_ebe_06;
    ang::se3_tangent xi_ebb_07 = ang::se3_tangent::wedge(M_xb % ang::se3_tangent::hat_homo(xi_ebe_07));
    ang::se3_tangent xi_ebb_08 = ang::se3_tangent::wedge(Gr_xb % ang::se3_tangent::hat_homo(xi_ebe_08));
    ang::se3_tangent xi_ebb_09 = ang::se3_tangent::wedge(Gq_xb % ang::se3_tangent::hat_homo(xi_ebe_09));
    ang::se3_tangent xi_ebb_11 = ang::se3_tangent::wedge(Z_xb % ang::se3_tangent::hat_dual(xi_ebe_11));

    check("adjoint inv 11    ", xi_ebb_0()(0), xi_ebb_01()(0), 1e-12);
    check("adjoint inv 12    ", xi_ebb_0()(1), xi_ebb_01()(1), 1e-12);
    check("adjoint inv 13    ", xi_ebb_0()(2), xi_ebb_01()(2), 1e-12);
    check("adjoint inv 14    ", xi_ebb_0()(3), xi_ebb_01()(3), 1e-12);
    check("adjoint inv 15    ", xi_ebb_0()(4), xi_ebb_01()(4), 1e-12);
    check("adjoint inv 16    ", xi_ebb_0()(5), xi_ebb_01()(5), 1e-12);

    check("adjoint inv 21    ", xi_ebb_0()(0), xi_ebb_02()(0), 1e-12);
    check("adjoint inv 22    ", xi_ebb_0()(1), xi_ebb_02()(1), 1e-12);
    check("adjoint inv 23    ", xi_ebb_0()(2), xi_ebb_02()(2), 1e-12);
    check("adjoint inv 24    ", xi_ebb_0()(3), xi_ebb_02()(3), 1e-12);
    check("adjoint inv 25    ", xi_ebb_0()(4), xi_ebb_02()(4), 1e-12);
    check("adjoint inv 26    ", xi_ebb_0()(5), xi_ebb_02()(5), 1e-12);

    check("adjoint inv 31    ", xi_ebb_0()(0), xi_ebb_03()(0), 1e-12);
    check("adjoint inv 32    ", xi_ebb_0()(1), xi_ebb_03()(1), 1e-12);
    check("adjoint inv 33    ", xi_ebb_0()(2), xi_ebb_03()(2), 1e-12);
    check("adjoint inv 34    ", xi_ebb_0()(3), xi_ebb_03()(3), 1e-12);
    check("adjoint inv 35    ", xi_ebb_0()(4), xi_ebb_03()(4), 1e-12);
    check("adjoint inv 36    ", xi_ebb_0()(5), xi_ebb_03()(5), 1e-12);

    check("adjoint inv 41    ", xi_ebb_0()(0), xi_ebb_04()(0), 1e-12);
    check("adjoint inv 42    ", xi_ebb_0()(1), xi_ebb_04()(1), 1e-12);
    check("adjoint inv 43    ", xi_ebb_0()(2), xi_ebb_04()(2), 1e-12);
    check("adjoint inv 44    ", xi_ebb_0()(3), xi_ebb_04()(3), 1e-12);
    check("adjoint inv 45    ", xi_ebb_0()(4), xi_ebb_04()(4), 1e-12);
    check("adjoint inv 46    ", xi_ebb_0()(5), xi_ebb_04()(5), 1e-12);

    check("adjoint inv 51    ", xi_ebb_0()(0), xi_ebb_05()(0), 1e-12);
    check("adjoint inv 52    ", xi_ebb_0()(1), xi_ebb_05()(1), 1e-12);
    check("adjoint inv 53    ", xi_ebb_0()(2), xi_ebb_05()(2), 1e-12);
    check("adjoint inv 54    ", xi_ebb_0()(3), xi_ebb_05()(3), 1e-12);
    check("adjoint inv 55    ", xi_ebb_0()(4), xi_ebb_05()(4), 1e-12);
    check("adjoint inv 56    ", xi_ebb_0()(5), xi_ebb_05()(5), 1e-12);

    check("adjoint inv 61    ", xi_ebb_0()(0), xi_ebb_06()(0), 1e-12);
    check("adjoint inv 62    ", xi_ebb_0()(1), xi_ebb_06()(1), 1e-12);
    check("adjoint inv 63    ", xi_ebb_0()(2), xi_ebb_06()(2), 1e-12);
    check("adjoint inv 64    ", xi_ebb_0()(3), xi_ebb_06()(3), 1e-12);
    check("adjoint inv 65    ", xi_ebb_0()(4), xi_ebb_06()(4), 1e-12);
    check("adjoint inv 66    ", xi_ebb_0()(5), xi_ebb_06()(5), 1e-12);

    check("adjoint inv 71    ", xi_ebb_0()(0), xi_ebb_07()(0), 1e-12);
    check("adjoint inv 72    ", xi_ebb_0()(1), xi_ebb_07()(1), 1e-12);
    check("adjoint inv 73    ", xi_ebb_0()(2), xi_ebb_07()(2), 1e-12);
    check("adjoint inv 74    ", xi_ebb_0()(3), xi_ebb_07()(3), 1e-12);
    check("adjoint inv 75    ", xi_ebb_0()(4), xi_ebb_07()(4), 1e-12);
    check("adjoint inv 76    ", xi_ebb_0()(5), xi_ebb_07()(5), 1e-12);

    check("adjoint inv 81    ", xi_ebb_0()(0), xi_ebb_08()(0), 1e-12);
    check("adjoint inv 82    ", xi_ebb_0()(1), xi_ebb_08()(1), 1e-12);
    check("adjoint inv 83    ", xi_ebb_0()(2), xi_ebb_08()(2), 1e-12);
    check("adjoint inv 84    ", xi_ebb_0()(3), xi_ebb_08()(3), 1e-12);
    check("adjoint inv 85    ", xi_ebb_0()(4), xi_ebb_08()(4), 1e-12);
    check("adjoint inv 86    ", xi_ebb_0()(5), xi_ebb_08()(5), 1e-12);

    check("adjoint inv 91    ", xi_ebb_0()(0), xi_ebb_09()(0), 1e-12);
    check("adjoint inv 92    ", xi_ebb_0()(1), xi_ebb_09()(1), 1e-12);
    check("adjoint inv 93    ", xi_ebb_0()(2), xi_ebb_09()(2), 1e-12);
    check("adjoint inv 94    ", xi_ebb_0()(3), xi_ebb_09()(3), 1e-12);
    check("adjoint inv 95    ", xi_ebb_0()(4), xi_ebb_09()(4), 1e-12);
    check("adjoint inv 96    ", xi_ebb_0()(5), xi_ebb_09()(5), 1e-12);

    check("adjoint inv 111   ", xi_ebb_0()(0), xi_ebb_11()(0), 1e-12);
    check("adjoint inv 112   ", xi_ebb_0()(1), xi_ebb_11()(1), 1e-12);
    check("adjoint inv 113   ", xi_ebb_0()(2), xi_ebb_11()(2), 1e-12);
    check("adjoint inv 114   ", xi_ebb_0()(3), xi_ebb_11()(3), 1e-12);
    check("adjoint inv 115   ", xi_ebb_0()(4), xi_ebb_11()(4), 1e-12);
    check("adjoint inv 116   ", xi_ebb_0()(5), xi_ebb_11()(5), 1e-12);



    Eigen::Matrix6d Adf1 = Gr_xb.adjoint_matrix_forward();
    Eigen::Matrix6d Adf2 = Gq_xb.adjoint_matrix_forward();
    Eigen::Matrix6d Adf3 = M_xb.adjoint_matrix_forward();
    Eigen::Matrix6d Adf4 = Z_xb.adjoint_matrix_forward();
    Eigen::Matrix6d Adf5 = tau_xb.adjoint_matrix_forward();
    Eigen::Matrix6d Adf6 = S_xb.adjoint_matrix_forward();
    Eigen::Vector6d xi_ebe_12 = Adf1 * xi_ebb_0();

    check("adjoint for 11    ", Adf1(0,0), Adf2(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf2(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf2(0,2), 1e-12);
    check("adjoint for 14    ", Adf1(0,3), Adf2(0,3), 1e-12);
    check("adjoint for 15    ", Adf1(0,4), Adf2(0,4), 1e-12);
    check("adjoint for 16    ", Adf1(0,5), Adf2(0,5), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf2(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf2(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf2(1,2), 1e-12);
    check("adjoint for 24    ", Adf1(1,3), Adf2(1,3), 1e-12);
    check("adjoint for 25    ", Adf1(1,4), Adf2(1,4), 1e-12);
    check("adjoint for 26    ", Adf1(1,5), Adf2(1,5), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf2(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf2(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf2(2,2), 1e-12);
    check("adjoint for 34    ", Adf1(2,3), Adf2(2,3), 1e-12);
    check("adjoint for 35    ", Adf1(2,4), Adf2(2,4), 1e-12);
    check("adjoint for 36    ", Adf1(2,5), Adf2(2,5), 1e-12);
    check("adjoint for 41    ", Adf1(3,0), Adf2(3,0), 1e-12);
    check("adjoint for 42    ", Adf1(3,1), Adf2(3,1), 1e-12);
    check("adjoint for 43    ", Adf1(3,2), Adf2(3,2), 1e-12);
    check("adjoint for 44    ", Adf1(3,3), Adf2(3,3), 1e-12);
    check("adjoint for 45    ", Adf1(3,4), Adf2(3,4), 1e-12);
    check("adjoint for 46    ", Adf1(3,5), Adf2(3,5), 1e-12);
    check("adjoint for 51    ", Adf1(4,0), Adf2(4,0), 1e-12);
    check("adjoint for 52    ", Adf1(4,1), Adf2(4,1), 1e-12);
    check("adjoint for 53    ", Adf1(4,2), Adf2(4,2), 1e-12);
    check("adjoint for 54    ", Adf1(4,3), Adf2(4,3), 1e-12);
    check("adjoint for 55    ", Adf1(4,4), Adf2(4,4), 1e-12);
    check("adjoint for 56    ", Adf1(4,5), Adf2(4,5), 1e-12);
    check("adjoint for 61    ", Adf1(5,0), Adf2(5,0), 1e-12);
    check("adjoint for 62    ", Adf1(5,1), Adf2(5,1), 1e-12);
    check("adjoint for 63    ", Adf1(5,2), Adf2(5,2), 1e-12);
    check("adjoint for 64    ", Adf1(5,3), Adf2(5,3), 1e-12);
    check("adjoint for 65    ", Adf1(5,4), Adf2(5,4), 1e-12);
    check("adjoint for 66    ", Adf1(5,5), Adf2(5,5), 1e-12);

    check("adjoint for 11    ", Adf1(0,0), Adf3(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf3(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf3(0,2), 1e-12);
    check("adjoint for 14    ", Adf1(0,3), Adf3(0,3), 1e-12);
    check("adjoint for 15    ", Adf1(0,4), Adf3(0,4), 1e-12);
    check("adjoint for 16    ", Adf1(0,5), Adf3(0,5), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf3(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf3(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf3(1,2), 1e-12);
    check("adjoint for 24    ", Adf1(1,3), Adf3(1,3), 1e-12);
    check("adjoint for 25    ", Adf1(1,4), Adf3(1,4), 1e-12);
    check("adjoint for 26    ", Adf1(1,5), Adf3(1,5), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf3(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf3(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf3(2,2), 1e-12);
    check("adjoint for 34    ", Adf1(2,3), Adf3(2,3), 1e-12);
    check("adjoint for 35    ", Adf1(2,4), Adf3(2,4), 1e-12);
    check("adjoint for 36    ", Adf1(2,5), Adf3(2,5), 1e-12);
    check("adjoint for 41    ", Adf1(3,0), Adf3(3,0), 1e-12);
    check("adjoint for 42    ", Adf1(3,1), Adf3(3,1), 1e-12);
    check("adjoint for 43    ", Adf1(3,2), Adf3(3,2), 1e-12);
    check("adjoint for 44    ", Adf1(3,3), Adf3(3,3), 1e-12);
    check("adjoint for 45    ", Adf1(3,4), Adf3(3,4), 1e-12);
    check("adjoint for 46    ", Adf1(3,5), Adf3(3,5), 1e-12);
    check("adjoint for 51    ", Adf1(4,0), Adf3(4,0), 1e-12);
    check("adjoint for 52    ", Adf1(4,1), Adf3(4,1), 1e-12);
    check("adjoint for 53    ", Adf1(4,2), Adf3(4,2), 1e-12);
    check("adjoint for 54    ", Adf1(4,3), Adf3(4,3), 1e-12);
    check("adjoint for 55    ", Adf1(4,4), Adf3(4,4), 1e-12);
    check("adjoint for 56    ", Adf1(4,5), Adf3(4,5), 1e-12);
    check("adjoint for 61    ", Adf1(5,0), Adf3(5,0), 1e-12);
    check("adjoint for 62    ", Adf1(5,1), Adf3(5,1), 1e-12);
    check("adjoint for 63    ", Adf1(5,2), Adf3(5,2), 1e-12);
    check("adjoint for 64    ", Adf1(5,3), Adf3(5,3), 1e-12);
    check("adjoint for 65    ", Adf1(5,4), Adf3(5,4), 1e-12);
    check("adjoint for 66    ", Adf1(5,5), Adf3(5,5), 1e-12);

    check("adjoint for 11    ", Adf1(0,0), Adf4(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf4(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf4(0,2), 1e-12);
    check("adjoint for 14    ", Adf1(0,3), Adf4(0,3), 1e-12);
    check("adjoint for 15    ", Adf1(0,4), Adf4(0,4), 1e-12);
    check("adjoint for 16    ", Adf1(0,5), Adf4(0,5), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf4(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf4(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf4(1,2), 1e-12);
    check("adjoint for 24    ", Adf1(1,3), Adf4(1,3), 1e-12);
    check("adjoint for 25    ", Adf1(1,4), Adf4(1,4), 1e-12);
    check("adjoint for 26    ", Adf1(1,5), Adf4(1,5), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf4(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf4(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf4(2,2), 1e-12);
    check("adjoint for 34    ", Adf1(2,3), Adf4(2,3), 1e-12);
    check("adjoint for 35    ", Adf1(2,4), Adf4(2,4), 1e-12);
    check("adjoint for 36    ", Adf1(2,5), Adf4(2,5), 1e-12);
    check("adjoint for 41    ", Adf1(3,0), Adf4(3,0), 1e-12);
    check("adjoint for 42    ", Adf1(3,1), Adf4(3,1), 1e-12);
    check("adjoint for 43    ", Adf1(3,2), Adf4(3,2), 1e-12);
    check("adjoint for 44    ", Adf1(3,3), Adf4(3,3), 1e-12);
    check("adjoint for 45    ", Adf1(3,4), Adf4(3,4), 1e-12);
    check("adjoint for 46    ", Adf1(3,5), Adf4(3,5), 1e-12);
    check("adjoint for 51    ", Adf1(4,0), Adf4(4,0), 1e-12);
    check("adjoint for 52    ", Adf1(4,1), Adf4(4,1), 1e-12);
    check("adjoint for 53    ", Adf1(4,2), Adf4(4,2), 1e-12);
    check("adjoint for 54    ", Adf1(4,3), Adf4(4,3), 1e-12);
    check("adjoint for 55    ", Adf1(4,4), Adf4(4,4), 1e-12);
    check("adjoint for 56    ", Adf1(4,5), Adf4(4,5), 1e-12);
    check("adjoint for 61    ", Adf1(5,0), Adf4(5,0), 1e-12);
    check("adjoint for 62    ", Adf1(5,1), Adf4(5,1), 1e-12);
    check("adjoint for 63    ", Adf1(5,2), Adf4(5,2), 1e-12);
    check("adjoint for 64    ", Adf1(5,3), Adf4(5,3), 1e-12);
    check("adjoint for 65    ", Adf1(5,4), Adf4(5,4), 1e-12);
    check("adjoint for 66    ", Adf1(5,5), Adf4(5,5), 1e-12);

    check("adjoint for 11    ", Adf1(0,0), Adf5(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf5(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf5(0,2), 1e-12);
    check("adjoint for 14    ", Adf1(0,3), Adf5(0,3), 1e-12);
    check("adjoint for 15    ", Adf1(0,4), Adf5(0,4), 1e-12);
    check("adjoint for 16    ", Adf1(0,5), Adf5(0,5), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf5(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf5(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf5(1,2), 1e-12);
    check("adjoint for 24    ", Adf1(1,3), Adf5(1,3), 1e-12);
    check("adjoint for 25    ", Adf1(1,4), Adf5(1,4), 1e-12);
    check("adjoint for 26    ", Adf1(1,5), Adf5(1,5), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf5(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf5(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf5(2,2), 1e-12);
    check("adjoint for 34    ", Adf1(2,3), Adf5(2,3), 1e-12);
    check("adjoint for 35    ", Adf1(2,4), Adf5(2,4), 1e-12);
    check("adjoint for 36    ", Adf1(2,5), Adf5(2,5), 1e-12);
    check("adjoint for 41    ", Adf1(3,0), Adf5(3,0), 1e-12);
    check("adjoint for 42    ", Adf1(3,1), Adf5(3,1), 1e-12);
    check("adjoint for 43    ", Adf1(3,2), Adf5(3,2), 1e-12);
    check("adjoint for 44    ", Adf1(3,3), Adf5(3,3), 1e-12);
    check("adjoint for 45    ", Adf1(3,4), Adf5(3,4), 1e-12);
    check("adjoint for 46    ", Adf1(3,5), Adf5(3,5), 1e-12);
    check("adjoint for 51    ", Adf1(4,0), Adf5(4,0), 1e-12);
    check("adjoint for 52    ", Adf1(4,1), Adf5(4,1), 1e-12);
    check("adjoint for 53    ", Adf1(4,2), Adf5(4,2), 1e-12);
    check("adjoint for 54    ", Adf1(4,3), Adf5(4,3), 1e-12);
    check("adjoint for 55    ", Adf1(4,4), Adf5(4,4), 1e-12);
    check("adjoint for 56    ", Adf1(4,5), Adf5(4,5), 1e-12);
    check("adjoint for 61    ", Adf1(5,0), Adf5(5,0), 1e-12);
    check("adjoint for 62    ", Adf1(5,1), Adf5(5,1), 1e-12);
    check("adjoint for 63    ", Adf1(5,2), Adf5(5,2), 1e-12);
    check("adjoint for 64    ", Adf1(5,3), Adf5(5,3), 1e-12);
    check("adjoint for 65    ", Adf1(5,4), Adf5(5,4), 1e-12);
    check("adjoint for 66    ", Adf1(5,5), Adf5(5,5), 1e-12);

    check("adjoint for 11    ", Adf1(0,0), Adf6(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf6(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf6(0,2), 1e-12);
    check("adjoint for 14    ", Adf1(0,3), Adf6(0,3), 1e-12);
    check("adjoint for 15    ", Adf1(0,4), Adf6(0,4), 1e-12);
    check("adjoint for 16    ", Adf1(0,5), Adf6(0,5), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf6(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf6(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf6(1,2), 1e-12);
    check("adjoint for 24    ", Adf1(1,3), Adf6(1,3), 1e-12);
    check("adjoint for 25    ", Adf1(1,4), Adf6(1,4), 1e-12);
    check("adjoint for 26    ", Adf1(1,5), Adf6(1,5), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf6(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf6(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf6(2,2), 1e-12);
    check("adjoint for 34    ", Adf1(2,3), Adf6(2,3), 1e-12);
    check("adjoint for 35    ", Adf1(2,4), Adf6(2,4), 1e-12);
    check("adjoint for 36    ", Adf1(2,5), Adf6(2,5), 1e-12);
    check("adjoint for 41    ", Adf1(3,0), Adf6(3,0), 1e-12);
    check("adjoint for 42    ", Adf1(3,1), Adf6(3,1), 1e-12);
    check("adjoint for 43    ", Adf1(3,2), Adf6(3,2), 1e-12);
    check("adjoint for 44    ", Adf1(3,3), Adf6(3,3), 1e-12);
    check("adjoint for 45    ", Adf1(3,4), Adf6(3,4), 1e-12);
    check("adjoint for 46    ", Adf1(3,5), Adf6(3,5), 1e-12);
    check("adjoint for 51    ", Adf1(4,0), Adf6(4,0), 1e-12);
    check("adjoint for 52    ", Adf1(4,1), Adf6(4,1), 1e-12);
    check("adjoint for 53    ", Adf1(4,2), Adf6(4,2), 1e-12);
    check("adjoint for 54    ", Adf1(4,3), Adf6(4,3), 1e-12);
    check("adjoint for 55    ", Adf1(4,4), Adf6(4,4), 1e-12);
    check("adjoint for 56    ", Adf1(4,5), Adf6(4,5), 1e-12);
    check("adjoint for 61    ", Adf1(5,0), Adf6(5,0), 1e-12);
    check("adjoint for 62    ", Adf1(5,1), Adf6(5,1), 1e-12);
    check("adjoint for 63    ", Adf1(5,2), Adf6(5,2), 1e-12);
    check("adjoint for 64    ", Adf1(5,3), Adf6(5,3), 1e-12);
    check("adjoint for 65    ", Adf1(5,4), Adf6(5,4), 1e-12);
    check("adjoint for 66    ", Adf1(5,5), Adf6(5,5), 1e-12);

    check("adjoint for 1     ", xi_ebe_12(0), xi_ebe_01()(0), 1e-12);
    check("adjoint for 2     ", xi_ebe_12(1), xi_ebe_01()(1), 1e-12);
    check("adjoint for 3     ", xi_ebe_12(2), xi_ebe_01()(2), 1e-12);
    check("adjoint for 4     ", xi_ebe_12(3), xi_ebe_01()(3), 1e-12);
    check("adjoint for 5     ", xi_ebe_12(4), xi_ebe_01()(4), 1e-12);
    check("adjoint for 6     ", xi_ebe_12(5), xi_ebe_01()(5), 1e-12);

    Eigen::Matrix6d Adb1 = Gr_xb.adjoint_matrix_backward();
    Eigen::Matrix6d Adb2 = Gq_xb.adjoint_matrix_backward();
    Eigen::Matrix6d Adb3 = M_xb.adjoint_matrix_backward();
    Eigen::Matrix6d Adb4 = Z_xb.adjoint_matrix_backward();
    Eigen::Matrix6d Adb5 = tau_xb.adjoint_matrix_backward();
    Eigen::Matrix6d Adb6 = S_xb.adjoint_matrix_backward();
    Eigen::Vector6d xi_ebb_12 = Adb1 * xi_ebe_12;

    check("adjoint bac 11    ", Adb1(0,0), Adb2(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb2(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb2(0,2), 1e-12);
    check("adjoint bac 14    ", Adb1(0,3), Adb2(0,3), 1e-12);
    check("adjoint bac 15    ", Adb1(0,4), Adb2(0,4), 1e-12);
    check("adjoint bac 16    ", Adb1(0,5), Adb2(0,5), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb2(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb2(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb2(1,2), 1e-12);
    check("adjoint bac 24    ", Adb1(1,3), Adb2(1,3), 1e-12);
    check("adjoint bac 25    ", Adb1(1,4), Adb2(1,4), 1e-12);
    check("adjoint bac 26    ", Adb1(1,5), Adb2(1,5), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb2(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb2(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb2(2,2), 1e-12);
    check("adjoint bac 34    ", Adb1(2,3), Adb2(2,3), 1e-12);
    check("adjoint bac 35    ", Adb1(2,4), Adb2(2,4), 1e-12);
    check("adjoint bac 36    ", Adb1(2,5), Adb2(2,5), 1e-12);
    check("adjoint bac 41    ", Adb1(3,0), Adb2(3,0), 1e-12);
    check("adjoint bac 42    ", Adb1(3,1), Adb2(3,1), 1e-12);
    check("adjoint bac 43    ", Adb1(3,2), Adb2(3,2), 1e-12);
    check("adjoint bac 44    ", Adb1(3,3), Adb2(3,3), 1e-12);
    check("adjoint bac 45    ", Adb1(3,4), Adb2(3,4), 1e-12);
    check("adjoint bac 46    ", Adb1(3,5), Adb2(3,5), 1e-12);
    check("adjoint bac 51    ", Adb1(4,0), Adb2(4,0), 1e-12);
    check("adjoint bac 52    ", Adb1(4,1), Adb2(4,1), 1e-12);
    check("adjoint bac 53    ", Adb1(4,2), Adb2(4,2), 1e-12);
    check("adjoint bac 54    ", Adb1(4,3), Adb2(4,3), 1e-12);
    check("adjoint bac 55    ", Adb1(4,4), Adb2(4,4), 1e-12);
    check("adjoint bac 56    ", Adb1(4,5), Adb2(4,5), 1e-12);
    check("adjoint bac 61    ", Adb1(5,0), Adb2(5,0), 1e-12);
    check("adjoint bac 62    ", Adb1(5,1), Adb2(5,1), 1e-12);
    check("adjoint bac 63    ", Adb1(5,2), Adb2(5,2), 1e-12);
    check("adjoint bac 64    ", Adb1(5,3), Adb2(5,3), 1e-12);
    check("adjoint bac 65    ", Adb1(5,4), Adb2(5,4), 1e-12);
    check("adjoint bac 66    ", Adb1(5,5), Adb2(5,5), 1e-12);

    check("adjoint bac 11    ", Adb1(0,0), Adb3(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb3(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb3(0,2), 1e-12);
    check("adjoint bac 14    ", Adb1(0,3), Adb3(0,3), 1e-12);
    check("adjoint bac 15    ", Adb1(0,4), Adb3(0,4), 1e-12);
    check("adjoint bac 16    ", Adb1(0,5), Adb3(0,5), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb3(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb3(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb3(1,2), 1e-12);
    check("adjoint bac 24    ", Adb1(1,3), Adb3(1,3), 1e-12);
    check("adjoint bac 25    ", Adb1(1,4), Adb3(1,4), 1e-12);
    check("adjoint bac 26    ", Adb1(1,5), Adb3(1,5), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb3(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb3(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb3(2,2), 1e-12);
    check("adjoint bac 34    ", Adb1(2,3), Adb3(2,3), 1e-12);
    check("adjoint bac 35    ", Adb1(2,4), Adb3(2,4), 1e-12);
    check("adjoint bac 36    ", Adb1(2,5), Adb3(2,5), 1e-12);
    check("adjoint bac 41    ", Adb1(3,0), Adb3(3,0), 1e-12);
    check("adjoint bac 42    ", Adb1(3,1), Adb3(3,1), 1e-12);
    check("adjoint bac 43    ", Adb1(3,2), Adb3(3,2), 1e-12);
    check("adjoint bac 44    ", Adb1(3,3), Adb3(3,3), 1e-12);
    check("adjoint bac 45    ", Adb1(3,4), Adb3(3,4), 1e-12);
    check("adjoint bac 46    ", Adb1(3,5), Adb3(3,5), 1e-12);
    check("adjoint bac 51    ", Adb1(4,0), Adb3(4,0), 1e-12);
    check("adjoint bac 52    ", Adb1(4,1), Adb3(4,1), 1e-12);
    check("adjoint bac 53    ", Adb1(4,2), Adb3(4,2), 1e-12);
    check("adjoint bac 54    ", Adb1(4,3), Adb3(4,3), 1e-12);
    check("adjoint bac 55    ", Adb1(4,4), Adb3(4,4), 1e-12);
    check("adjoint bac 56    ", Adb1(4,5), Adb3(4,5), 1e-12);
    check("adjoint bac 61    ", Adb1(5,0), Adb3(5,0), 1e-12);
    check("adjoint bac 62    ", Adb1(5,1), Adb3(5,1), 1e-12);
    check("adjoint bac 63    ", Adb1(5,2), Adb3(5,2), 1e-12);
    check("adjoint bac 64    ", Adb1(5,3), Adb3(5,3), 1e-12);
    check("adjoint bac 65    ", Adb1(5,4), Adb3(5,4), 1e-12);
    check("adjoint bac 66    ", Adb1(5,5), Adb3(5,5), 1e-12);

    check("adjoint bac 11    ", Adb1(0,0), Adb4(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb4(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb4(0,2), 1e-12);
    check("adjoint bac 14    ", Adb1(0,3), Adb4(0,3), 1e-12);
    check("adjoint bac 15    ", Adb1(0,4), Adb4(0,4), 1e-12);
    check("adjoint bac 16    ", Adb1(0,5), Adb4(0,5), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb4(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb4(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb4(1,2), 1e-12);
    check("adjoint bac 24    ", Adb1(1,3), Adb4(1,3), 1e-12);
    check("adjoint bac 25    ", Adb1(1,4), Adb4(1,4), 1e-12);
    check("adjoint bac 26    ", Adb1(1,5), Adb4(1,5), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb4(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb4(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb4(2,2), 1e-12);
    check("adjoint bac 34    ", Adb1(2,3), Adb4(2,3), 1e-12);
    check("adjoint bac 35    ", Adb1(2,4), Adb4(2,4), 1e-12);
    check("adjoint bac 36    ", Adb1(2,5), Adb4(2,5), 1e-12);
    check("adjoint bac 41    ", Adb1(3,0), Adb4(3,0), 1e-12);
    check("adjoint bac 42    ", Adb1(3,1), Adb4(3,1), 1e-12);
    check("adjoint bac 43    ", Adb1(3,2), Adb4(3,2), 1e-12);
    check("adjoint bac 44    ", Adb1(3,3), Adb4(3,3), 1e-12);
    check("adjoint bac 45    ", Adb1(3,4), Adb4(3,4), 1e-12);
    check("adjoint bac 46    ", Adb1(3,5), Adb4(3,5), 1e-12);
    check("adjoint bac 51    ", Adb1(4,0), Adb4(4,0), 1e-12);
    check("adjoint bac 52    ", Adb1(4,1), Adb4(4,1), 1e-12);
    check("adjoint bac 53    ", Adb1(4,2), Adb4(4,2), 1e-12);
    check("adjoint bac 54    ", Adb1(4,3), Adb4(4,3), 1e-12);
    check("adjoint bac 55    ", Adb1(4,4), Adb4(4,4), 1e-12);
    check("adjoint bac 56    ", Adb1(4,5), Adb4(4,5), 1e-12);
    check("adjoint bac 61    ", Adb1(5,0), Adb4(5,0), 1e-12);
    check("adjoint bac 62    ", Adb1(5,1), Adb4(5,1), 1e-12);
    check("adjoint bac 63    ", Adb1(5,2), Adb4(5,2), 1e-12);
    check("adjoint bac 64    ", Adb1(5,3), Adb4(5,3), 1e-12);
    check("adjoint bac 65    ", Adb1(5,4), Adb4(5,4), 1e-12);
    check("adjoint bac 66    ", Adb1(5,5), Adb4(5,5), 1e-12);

    check("adjoint bac 11    ", Adb1(0,0), Adb5(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb5(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb5(0,2), 1e-12);
    check("adjoint bac 14    ", Adb1(0,3), Adb5(0,3), 1e-12);
    check("adjoint bac 15    ", Adb1(0,4), Adb5(0,4), 1e-12);
    check("adjoint bac 16    ", Adb1(0,5), Adb5(0,5), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb5(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb5(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb5(1,2), 1e-12);
    check("adjoint bac 24    ", Adb1(1,3), Adb5(1,3), 1e-12);
    check("adjoint bac 25    ", Adb1(1,4), Adb5(1,4), 1e-12);
    check("adjoint bac 26    ", Adb1(1,5), Adb5(1,5), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb5(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb5(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb5(2,2), 1e-12);
    check("adjoint bac 34    ", Adb1(2,3), Adb5(2,3), 1e-12);
    check("adjoint bac 35    ", Adb1(2,4), Adb5(2,4), 1e-12);
    check("adjoint bac 36    ", Adb1(2,5), Adb5(2,5), 1e-12);
    check("adjoint bac 41    ", Adb1(3,0), Adb5(3,0), 1e-12);
    check("adjoint bac 42    ", Adb1(3,1), Adb5(3,1), 1e-12);
    check("adjoint bac 43    ", Adb1(3,2), Adb5(3,2), 1e-12);
    check("adjoint bac 44    ", Adb1(3,3), Adb5(3,3), 1e-12);
    check("adjoint bac 45    ", Adb1(3,4), Adb5(3,4), 1e-12);
    check("adjoint bac 46    ", Adb1(3,5), Adb5(3,5), 1e-12);
    check("adjoint bac 51    ", Adb1(4,0), Adb5(4,0), 1e-12);
    check("adjoint bac 52    ", Adb1(4,1), Adb5(4,1), 1e-12);
    check("adjoint bac 53    ", Adb1(4,2), Adb5(4,2), 1e-12);
    check("adjoint bac 54    ", Adb1(4,3), Adb5(4,3), 1e-12);
    check("adjoint bac 55    ", Adb1(4,4), Adb5(4,4), 1e-12);
    check("adjoint bac 56    ", Adb1(4,5), Adb5(4,5), 1e-12);
    check("adjoint bac 61    ", Adb1(5,0), Adb5(5,0), 1e-12);
    check("adjoint bac 62    ", Adb1(5,1), Adb5(5,1), 1e-12);
    check("adjoint bac 63    ", Adb1(5,2), Adb5(5,2), 1e-12);
    check("adjoint bac 64    ", Adb1(5,3), Adb5(5,3), 1e-12);
    check("adjoint bac 65    ", Adb1(5,4), Adb5(5,4), 1e-12);
    check("adjoint bac 66    ", Adb1(5,5), Adb5(5,5), 1e-12);

    check("adjoint bac 11    ", Adb1(0,0), Adb6(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb6(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb6(0,2), 1e-12);
    check("adjoint bac 14    ", Adb1(0,3), Adb6(0,3), 1e-12);
    check("adjoint bac 15    ", Adb1(0,4), Adb6(0,4), 1e-12);
    check("adjoint bac 16    ", Adb1(0,5), Adb6(0,5), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb6(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb6(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb6(1,2), 1e-12);
    check("adjoint bac 24    ", Adb1(1,3), Adb6(1,3), 1e-12);
    check("adjoint bac 25    ", Adb1(1,4), Adb6(1,4), 1e-12);
    check("adjoint bac 26    ", Adb1(1,5), Adb6(1,5), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb6(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb6(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb6(2,2), 1e-12);
    check("adjoint bac 34    ", Adb1(2,3), Adb6(2,3), 1e-12);
    check("adjoint bac 35    ", Adb1(2,4), Adb6(2,4), 1e-12);
    check("adjoint bac 36    ", Adb1(2,5), Adb6(2,5), 1e-12);
    check("adjoint bac 41    ", Adb1(3,0), Adb6(3,0), 1e-12);
    check("adjoint bac 42    ", Adb1(3,1), Adb6(3,1), 1e-12);
    check("adjoint bac 43    ", Adb1(3,2), Adb6(3,2), 1e-12);
    check("adjoint bac 44    ", Adb1(3,3), Adb6(3,3), 1e-12);
    check("adjoint bac 45    ", Adb1(3,4), Adb6(3,4), 1e-12);
    check("adjoint bac 46    ", Adb1(3,5), Adb6(3,5), 1e-12);
    check("adjoint bac 51    ", Adb1(4,0), Adb6(4,0), 1e-12);
    check("adjoint bac 52    ", Adb1(4,1), Adb6(4,1), 1e-12);
    check("adjoint bac 53    ", Adb1(4,2), Adb6(4,2), 1e-12);
    check("adjoint bac 54    ", Adb1(4,3), Adb6(4,3), 1e-12);
    check("adjoint bac 55    ", Adb1(4,4), Adb6(4,4), 1e-12);
    check("adjoint bac 56    ", Adb1(4,5), Adb6(4,5), 1e-12);
    check("adjoint bac 61    ", Adb1(5,0), Adb6(5,0), 1e-12);
    check("adjoint bac 62    ", Adb1(5,1), Adb6(5,1), 1e-12);
    check("adjoint bac 63    ", Adb1(5,2), Adb6(5,2), 1e-12);
    check("adjoint bac 64    ", Adb1(5,3), Adb6(5,3), 1e-12);
    check("adjoint bac 65    ", Adb1(5,4), Adb6(5,4), 1e-12);
    check("adjoint bac 66    ", Adb1(5,5), Adb6(5,5), 1e-12);

    check("adjoint bac 1     ", xi_ebb_12(0), xi_ebb_01()(0), 1e-12);
    check("adjoint bac 2     ", xi_ebb_12(1), xi_ebb_01()(1), 1e-12);
    check("adjoint bac 3     ", xi_ebb_12(2), xi_ebb_01()(2), 1e-12);
    check("adjoint bac 4     ", xi_ebb_12(3), xi_ebb_01()(3), 1e-12);
    check("adjoint bac 5     ", xi_ebb_12(4), xi_ebb_01()(4), 1e-12);
    check("adjoint bac 6     ", xi_ebb_12(5), xi_ebb_01()(5), 1e-12);

    ang::se3_tangent xi_ebb_1(+0.11, +3.87, -5.02, +2.08, -0.68, -0.14);
    ang::se3_tangent xi_ebb_2(-0.27, -0.55, +1.11, +2.42, +1.34, -1.58);
    ang::se3_tangent xi_ebb_3 = xi_ebb_1 + xi_ebb_2;
    ang::se3_tangent xi_ebb_4 = xi_ebb_1 - xi_ebb_2;

    ang::se3_tangent xi_ebe_1 = Z_xb | xi_ebb_1;
    ang::se3_tangent xi_ebe_2 = Z_xb | xi_ebb_2;
    ang::se3_tangent xi_ebe_3 = Z_xb | xi_ebb_3;
    ang::se3_tangent xi_ebe_4 = Z_xb | xi_ebb_4;

    check("adjoint + 111   ", xi_ebe_3()(0), xi_ebe_1()(0) + xi_ebe_2()(0), 1e-12);
    check("adjoint + 112   ", xi_ebe_3()(1), xi_ebe_1()(1) + xi_ebe_2()(1), 1e-12);
    check("adjoint + 113   ", xi_ebe_3()(2), xi_ebe_1()(2) + xi_ebe_2()(2), 1e-12);
    check("adjoint + 114   ", xi_ebe_3()(3), xi_ebe_1()(3) + xi_ebe_2()(3), 1e-12);
    check("adjoint + 115   ", xi_ebe_3()(4), xi_ebe_1()(4) + xi_ebe_2()(4), 1e-12);
    check("adjoint + 116   ", xi_ebe_3()(5), xi_ebe_1()(5) + xi_ebe_2()(5), 1e-12);

    check("adjoint - 111   ", xi_ebe_4()(0), xi_ebe_1()(0) - xi_ebe_2()(0), 1e-12);
    check("adjoint - 112   ", xi_ebe_4()(1), xi_ebe_1()(1) - xi_ebe_2()(1), 1e-12);
    check("adjoint - 113   ", xi_ebe_4()(2), xi_ebe_1()(2) - xi_ebe_2()(2), 1e-12);
    check("adjoint - 114   ", xi_ebe_4()(3), xi_ebe_1()(3) - xi_ebe_2()(3), 1e-12);
    check("adjoint - 115   ", xi_ebe_4()(4), xi_ebe_1()(4) - xi_ebe_2()(4), 1e-12);
    check("adjoint - 116   ", xi_ebe_4()(5), xi_ebe_1()(5) - xi_ebe_2()(5), 1e-12);

} // closes test_adjoint

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_velocity() {

    ang::se3_tangent xi_xbb(-0.89, 0.11, 1.34, 0.51, -0.38, -0.14);

    Eigen::Vector3d T_xbx(3.2, 4.1, -1.8);
    ang::euler euler_xb(0.12, 0.17, -0.03);
    ang::dcm R_xb(euler_xb);
    ang::speu_dcm GR_xb(R_xb, T_xbx);

    ang::se3_tangent xi_xbx = GR_xb | xi_xbb;

    Eigen::Vector3d p_b(7.1, -11.8, -9.4);
    Eigen::Vector3d v_bpb_0 = xi_xbb.get_vi() + xi_xbb.get_w()().cross(p_b);
    Eigen::Vector3d v_bpb_1 = xi_xbb.point_velocity(p_b);

    check("velocity 11    ", v_bpb_0(0), v_bpb_1(0), 1e-12);
    check("velocity 12    ", v_bpb_0(1), v_bpb_1(1), 1e-12);
    check("velocity 13    ", v_bpb_0(2), v_bpb_1(2), 1e-12);

    Eigen::Vector3d p_x = GR_xb * p_b;
    Eigen::Vector3d v_xpx_0 = GR_xb ^ v_bpb_0; // origin of x does NOT coincide with origin of b
    Eigen::Vector3d v_xpx_1 = xi_xbx.point_velocity(p_x);

    check("velocity 21    ", v_xpx_0(0), v_xpx_1(0), 1e-12);
    check("velocity 22    ", v_xpx_0(1), v_xpx_1(1), 1e-12);
    check("velocity 23    ", v_xpx_0(2), v_xpx_1(2), 1e-12);

} // closes test_velocity

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void ang::test::Tse3::test_dual_set() {

    Eigen::Vector3d T_ebe_m(3.2, 4.1, -1.8);
    ang::euler euler_eb(0.12, 0.17, -0.03);
    ang::rodrigues q_eb(euler_eb);
    ang::rotv r_eb(euler_eb);
    ang::dual z_eb1(r_eb, T_ebe_m);

    Eigen::Vector3d T_ebe_m_dummy(-0.7, 4.4, 2.8);
    ang::rotv r_eb_dummy(0.33, 0.14, -0.68);
    ang::dual z_eb2(r_eb_dummy, T_ebe_m_dummy);
    ang::dual z_eb3(r_eb_dummy, T_ebe_m_dummy);
    ang::dual z_eb4(r_eb_dummy, T_ebe_m_dummy);
    ang::dual z_eb5(r_eb_dummy, T_ebe_m_dummy);

    z_eb2.set_T(T_ebe_m);
    z_eb2.set_rotv(r_eb);

    z_eb3.set_T(T_ebe_m);
    z_eb3.set_rodrigues(q_eb);

    z_eb4.set_rotv(r_eb);
    z_eb4.set_T(T_ebe_m);

    z_eb5.set_rodrigues(q_eb);
    z_eb5.set_T(T_ebe_m);


    check("1st T 2nd rot 0    ", z_eb1.get_qr()()(0), z_eb2.get_qr()()(0),  1e-12);
    check("1st T 2nd rot 1    ", z_eb1.get_qr()()(1), z_eb2.get_qr()()(1),  1e-12);
    check("1st T 2nd rot 2    ", z_eb1.get_qr()()(2), z_eb2.get_qr()()(2),  1e-12);
    check("1st T 2nd rot 3    ", z_eb1.get_qr()()(3), z_eb2.get_qr()()(3),  1e-12);
    check("1st T 2nd rot 4    ", z_eb1.get_qd()(0),   z_eb2.get_qd()(0),  1e-12);
    check("1st T 2nd rot 5    ", z_eb1.get_qd()(1),   z_eb2.get_qd()(1),  1e-12);
    check("1st T 2nd rot 6    ", z_eb1.get_qd()(2),   z_eb2.get_qd()(2),  1e-12);
    check("1st T 2nd rot 7    ", z_eb1.get_qd()(3),   z_eb2.get_qd()(3),  1e-12);

    check("1st T 2nd rod 0    ", z_eb1.get_qr()()(0), z_eb3.get_qr()()(0),  1e-12);
    check("1st T 2nd rod 1    ", z_eb1.get_qr()()(1), z_eb3.get_qr()()(1),  1e-12);
    check("1st T 2nd rod 2    ", z_eb1.get_qr()()(2), z_eb3.get_qr()()(2),  1e-12);
    check("1st T 2nd rod 3    ", z_eb1.get_qr()()(3), z_eb3.get_qr()()(3),  1e-12);
    check("1st T 2nd rod 4    ", z_eb1.get_qd()(0),   z_eb3.get_qd()(0),  1e-12);
    check("1st T 2nd rod 5    ", z_eb1.get_qd()(1),   z_eb3.get_qd()(1),  1e-12);
    check("1st T 2nd rod 6    ", z_eb1.get_qd()(2),   z_eb3.get_qd()(2),  1e-12);
    check("1st T 2nd rod 7    ", z_eb1.get_qd()(3),   z_eb3.get_qd()(3),  1e-12);

    check("1st rot 2nd T 0    ", z_eb1.get_qr()()(0), z_eb4.get_qr()()(0),  1e-12);
    check("1st rot 2nd T 1    ", z_eb1.get_qr()()(1), z_eb4.get_qr()()(1),  1e-12);
    check("1st rot 2nd T 2    ", z_eb1.get_qr()()(2), z_eb4.get_qr()()(2),  1e-12);
    check("1st rot 2nd T 3    ", z_eb1.get_qr()()(3), z_eb4.get_qr()()(3),  1e-12);
    check("1st rot 2nd T 4    ", z_eb1.get_qd()(0),   z_eb4.get_qd()(0),  1e-12);
    check("1st rot 2nd T 5    ", z_eb1.get_qd()(1),   z_eb4.get_qd()(1),  1e-12);
    check("1st rot 2nd T 6    ", z_eb1.get_qd()(2),   z_eb4.get_qd()(2),  1e-12);
    check("1st rot 2nd T 7    ", z_eb1.get_qd()(3),   z_eb4.get_qd()(3),  1e-12);

    check("1st rod 2nd T 0    ", z_eb1.get_qr()()(0), z_eb5.get_qr()()(0),  1e-12);
    check("1st rod 2nd T 1    ", z_eb1.get_qr()()(1), z_eb5.get_qr()()(1),  1e-12);
    check("1st rod 2nd T 2    ", z_eb1.get_qr()()(2), z_eb5.get_qr()()(2),  1e-12);
    check("1st rod 2nd T 3    ", z_eb1.get_qr()()(3), z_eb5.get_qr()()(3),  1e-12);
    check("1st rod 2nd T 4    ", z_eb1.get_qd()(0),   z_eb5.get_qd()(0),  1e-12);
    check("1st rod 2nd T 5    ", z_eb1.get_qd()(1),   z_eb5.get_qd()(1),  1e-12);
    check("1st rod 2nd T 6    ", z_eb1.get_qd()(2),   z_eb5.get_qd()(2),  1e-12);
    check("1st rod 2nd T 7    ", z_eb1.get_qd()(3),   z_eb5.get_qd()(3),  1e-12);



} // closes test_dual_set

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




























