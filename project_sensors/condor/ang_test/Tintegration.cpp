#include "Tintegration.h"

#include <iostream>
#include "../ang/rotate/dcm.h"
#include "../ang/rotate/euler.h"
#include "../ang/rotate/rodrigues.h"
#include "../ang/rotate/rotv.h"
#include "../ang/transform/speu_rodrigues.h"
#include "../ang/transform/speu_dcm.h"
#include "../ang/transform/homogeneous.h"
#include "../ang/transform/trfv.h"
#include "../ang/transform/dual.h"
#include "../ang/transform/se3_tangent.h"
#include "../ang/auxiliary.h"
#include "../ang/tools.h"

using namespace std;

ang::test::Tintegration::Tintegration(jail::counter& Ocounter)
        : ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tintegration::run() {
    ::jail::unit_test::run();

    // w --> world (inertial)
    // b --> body
    // p --> point

    // Initial conditions - position
    Eigen::Vector3d x0_wbw_m(0.0, 0.0, 0.0);

    // Fixed position of point in body
    Eigen::Vector3d x_bpb_m(3.5, -1.2, 4.6);

    // Initial conditions - attitude
    ang::rotv rv0_wb(0.15, -0.04, 0.07);
    ang::rodrigues q0_wb(rv0_wb);
    ang::dcm R0_wb(rv0_wb);
    ang::euler euler0_wb(rv0_wb);

    // Initial conditions - transformation
    ang::speu_rodrigues Gq0_wb(q0_wb, x0_wbw_m);
    ang::speu_dcm GR0_wb(R0_wb, x0_wbw_m);
    ang::homogeneous M0_wb(R0_wb, x0_wbw_m);
    ang::trfv tau0_wb(rv0_wb, x0_wbw_m);
    ang::dual Z0_wb(q0_wb, x0_wbw_m);

    // Initial conditions - speed
    Eigen::Vector3d v0_w_mps(1., 2., 3.);
    Eigen::Vector3d v0_b_mps = q0_wb / v0_w_mps;

    // Initial rotation speed
    ang::so3_tangent w0_wbw_rps(0.015, 0.02, 0.039);
    ang::so3_tangent w0_wbb_rps = q0_wb % w0_wbw_rps;

    // Fixed body linear acceleration
    Eigen::Vector3d a_b_mps2(0.02, -0.03, 0.04);

    // Initial body angular acceleration
    Eigen::Vector3d alpha0_b_rps2(3.1e-3, -1.6e-3, 0.8e-3);

    // 1st order integration to simplify tests --> go to very high frequency to compensate
    double Deltat_sec = 0.001;
    double tend_sec = 100.0;

    // specific tests to analyze the results
    test_differences(GR0_wb, v0_w_mps, w0_wbw_rps, a_b_mps2, alpha0_b_rps2);
    test_analysis(GR0_wb, v0_w_mps, w0_wbw_rps, a_b_mps2, alpha0_b_rps2);

    std::vector<double> res_space_rodrigues(25), res_space_dcm(25), res_space_rotv(25), res_body_rodrigues(25), res_body_dcm(25), res_body_rotv(25);
    std::vector<double> res_space_speu_rodrigues(25), res_space_speu_dcm(25), res_space_homogeneous(25), res_space_trfv(25), res_space_dual(25), res_body_speu_rodrigues(25), res_body_speu_dcm(25), res_body_homogeneous(25), res_body_trfv(25), res_body_dual(25);

    test_space_rodrigues     (res_space_rodrigues,      Deltat_sec, tend_sec, x0_wbw_m, q0_wb,   v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_dcm           (res_space_dcm,            Deltat_sec, tend_sec, x0_wbw_m, R0_wb,   v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_rotv          (res_space_rotv,           Deltat_sec, tend_sec, x0_wbw_m, rv0_wb,  v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_rodrigues      (res_body_rodrigues,       Deltat_sec, tend_sec, x0_wbw_m, q0_wb,   v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_dcm            (res_body_dcm,             Deltat_sec, tend_sec, x0_wbw_m, R0_wb,   v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_rotv           (res_body_rotv,            Deltat_sec, tend_sec, x0_wbw_m, rv0_wb,  v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_speu_rodrigues(res_space_speu_rodrigues, Deltat_sec, tend_sec,           Gq0_wb,  v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_speu_dcm      (res_space_speu_dcm,       Deltat_sec, tend_sec,           GR0_wb,  v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_homogeneous   (res_space_homogeneous,    Deltat_sec, tend_sec,           M0_wb,   v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_trfv          (res_space_trfv,           Deltat_sec, tend_sec,           tau0_wb, v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_space_dual          (res_space_dual,           Deltat_sec, tend_sec,           Z0_wb,   v0_w_mps, w0_wbw_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_speu_rodrigues (res_body_speu_rodrigues,  Deltat_sec, tend_sec,           Gq0_wb,  v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_speu_dcm       (res_body_speu_dcm,        Deltat_sec, tend_sec,           GR0_wb,  v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_homogeneous    (res_body_homogeneous,     Deltat_sec, tend_sec,           M0_wb,   v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_trfv           (res_body_trfv,            Deltat_sec, tend_sec,           tau0_wb, v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);
    test_body_dual           (res_body_dual,            Deltat_sec, tend_sec,           Z0_wb,   v0_b_mps, w0_wbb_rps, x_bpb_m, a_b_mps2, alpha0_b_rps2);

    check("final 0102 t [sec]          ", res_space_rodrigues[0],  res_space_dcm[0],  1e-8);
    check("final 0102 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_dcm[1],  1e-8);
    check("final 0102 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_dcm[2],  1e-8);
    check("final 0102 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_dcm[3],  1e-8);
    check("final 0102 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_dcm[4],  1e-8);
    check("final 0102 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_dcm[5],  1e-8);
    check("final 0102 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_dcm[6],  1e-8);
    check("final 0102 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_dcm[7],  1e-8);
    check("final 0102 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_dcm[8],  1e-8);
    check("final 0102 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_dcm[9],  1e-8);
    check("final 0102 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_dcm[10], 1e-8);
    check("final 0102 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_dcm[11], 1e-8);
    check("final 0102 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_dcm[12], 1e-8);
    check("final 0102 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_dcm[13], 1e-6);
    check("final 0102 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_dcm[14], 1e-6);
    check("final 0102 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_dcm[15], 1e-6);
    check("final 0102 body yaw [deg]   ", res_space_rodrigues[16], res_space_dcm[16], 1e-6);
    check("final 0102 body pitch [deg] ", res_space_rodrigues[17], res_space_dcm[17], 1e-6);
    check("final 0102 body bank [deg]  ", res_space_rodrigues[18], res_space_dcm[18], 1e-6);
    check("final 0102 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_dcm[19], 1e-6);
    check("final 0102 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_dcm[20], 1e-6);
    check("final 0102 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_dcm[21], 1e-6);
    check("final 0102 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_dcm[22], 1e-8);
    check("final 0102 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_dcm[23], 1e-8);
    check("final 0102 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_dcm[24], 1e-8);

    check("final 0103 t [sec]          ", res_space_rodrigues[0],  res_space_rotv[0],  1e-8);
    check("final 0103 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_rotv[1],  1e-8);
    check("final 0103 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_rotv[2],  1e-8);
    check("final 0103 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_rotv[3],  1e-8);
    check("final 0103 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_rotv[4],  1e-8);
    check("final 0103 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_rotv[5],  1e-8);
    check("final 0103 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_rotv[6],  1e-8);
    check("final 0103 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_rotv[7],  1e-8);
    check("final 0103 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_rotv[8],  1e-8);
    check("final 0103 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_rotv[9],  1e-8);
    check("final 0103 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_rotv[10], 1e-8);
    check("final 0103 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_rotv[11], 1e-8);
    check("final 0103 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_rotv[12], 1e-8);
    check("final 0103 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_rotv[13], 1e-6);
    check("final 0103 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_rotv[14], 1e-6);
    check("final 0103 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_rotv[15], 1e-6);
    check("final 0103 body yaw [deg]   ", res_space_rodrigues[16], res_space_rotv[16], 1e-6);
    check("final 0103 body pitch [deg] ", res_space_rodrigues[17], res_space_rotv[17], 1e-6);
    check("final 0103 body bank [deg]  ", res_space_rodrigues[18], res_space_rotv[18], 1e-6);
    check("final 0103 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_rotv[19], 1e-6);
    check("final 0103 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_rotv[20], 1e-6);
    check("final 0103 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_rotv[21], 1e-6);
    check("final 0103 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_rotv[22], 1e-8);
    check("final 0103 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_rotv[23], 1e-8);
    check("final 0103 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_rotv[24], 1e-8);

    check("final 0111 t [sec]          ", res_space_rodrigues[0],  res_body_rodrigues[0],  1e-8);
    check("final 0111 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_rodrigues[1],  1e-8);
    check("final 0111 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_rodrigues[2],  1e-8);
    check("final 0111 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_rodrigues[3],  1e-8);
    check("final 0111 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_rodrigues[4],  1e-8);
    check("final 0111 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_rodrigues[5],  1e-8);
    check("final 0111 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_rodrigues[6],  1e-8);
    check("final 0111 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_rodrigues[7],  1e-8);
    check("final 0111 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_rodrigues[8],  1e-8);
    check("final 0111 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_rodrigues[9],  1e-8);
    check("final 0111 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_rodrigues[10], 1e-8);
    check("final 0111 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_rodrigues[11], 1e-8);
    check("final 0111 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_rodrigues[12], 1e-8);
    check("final 0111 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_rodrigues[13], 1e-6);
    check("final 0111 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_rodrigues[14], 1e-6);
    check("final 0111 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_rodrigues[15], 1e-6);
    check("final 0111 body yaw [deg]   ", res_space_rodrigues[16], res_body_rodrigues[16], 1e-6);
    check("final 0111 body pitch [deg] ", res_space_rodrigues[17], res_body_rodrigues[17], 1e-6);
    check("final 0111 body bank [deg]  ", res_space_rodrigues[18], res_body_rodrigues[18], 1e-6);
    check("final 0111 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_rodrigues[19], 1e-6);
    check("final 0111 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_rodrigues[20], 1e-6);
    check("final 0111 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_rodrigues[21], 1e-6);
    check("final 0111 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_rodrigues[22], 1e-8);
    check("final 0111 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_rodrigues[23], 1e-8);
    check("final 0111 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_rodrigues[24], 1e-8);

    check("final 0112 t [sec]          ", res_space_rodrigues[0],  res_body_dcm[0],  1e-8);
    check("final 0112 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_dcm[1],  1e-8);
    check("final 0112 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_dcm[2],  1e-8);
    check("final 0112 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_dcm[3],  1e-8);
    check("final 0112 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_dcm[4],  1e-8);
    check("final 0112 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_dcm[5],  1e-8);
    check("final 0112 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_dcm[6],  1e-8);
    check("final 0112 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_dcm[7],  1e-8);
    check("final 0112 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_dcm[8],  1e-8);
    check("final 0112 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_dcm[9],  1e-8);
    check("final 0112 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_dcm[10], 1e-8);
    check("final 0112 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_dcm[11], 1e-8);
    check("final 0112 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_dcm[12], 1e-8);
    check("final 0112 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_dcm[13], 1e-6);
    check("final 0112 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_dcm[14], 1e-6);
    check("final 0112 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_dcm[15], 1e-6);
    check("final 0112 body yaw [deg]   ", res_space_rodrigues[16], res_body_dcm[16], 1e-6);
    check("final 0112 body pitch [deg] ", res_space_rodrigues[17], res_body_dcm[17], 1e-6);
    check("final 0112 body bank [deg]  ", res_space_rodrigues[18], res_body_dcm[18], 1e-6);
    check("final 0112 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_dcm[19], 1e-6);
    check("final 0112 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_dcm[20], 1e-6);
    check("final 0112 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_dcm[21], 1e-6);
    check("final 0112 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_dcm[22], 1e-8);
    check("final 0112 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_dcm[23], 1e-8);
    check("final 0112 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_dcm[24], 1e-8);

    check("final 0113 t [sec]          ", res_space_rodrigues[0],  res_body_rotv[0],  1e-8);
    check("final 0113 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_rotv[1],  1e-8);
    check("final 0113 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_rotv[2],  1e-8);
    check("final 0113 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_rotv[3],  1e-8);
    check("final 0113 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_rotv[4],  1e-8);
    check("final 0113 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_rotv[5],  1e-8);
    check("final 0113 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_rotv[6],  1e-8);
    check("final 0113 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_rotv[7],  1e-8);
    check("final 0113 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_rotv[8],  1e-8);
    check("final 0113 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_rotv[9],  1e-8);
    check("final 0113 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_rotv[10], 1e-8);
    check("final 0113 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_rotv[11], 1e-8);
    check("final 0113 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_rotv[12], 1e-8);
    check("final 0113 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_rotv[13], 1e-6);
    check("final 0113 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_rotv[14], 1e-6);
    check("final 0113 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_rotv[15], 1e-6);
    check("final 0113 body yaw [deg]   ", res_space_rodrigues[16], res_body_rotv[16], 1e-6);
    check("final 0113 body pitch [deg] ", res_space_rodrigues[17], res_body_rotv[17], 1e-6);
    check("final 0113 body bank [deg]  ", res_space_rodrigues[18], res_body_rotv[18], 1e-6);
    check("final 0113 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_rotv[19], 1e-6);
    check("final 0113 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_rotv[20], 1e-6);
    check("final 0113 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_rotv[21], 1e-6);
    check("final 0113 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_rotv[22], 1e-8);
    check("final 0113 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_rotv[23], 1e-7);
    check("final 0113 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_rotv[24], 1e-7);

    check("final 0121 t [sec]          ", res_space_rodrigues[0],  res_space_speu_rodrigues[0],  1e-8);
    check("final 0121 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_speu_rodrigues[1],  1e-8);
    check("final 0121 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_speu_rodrigues[2],  1e-8);
    check("final 0121 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_speu_rodrigues[3],  1e-8);
    check("final 0121 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_speu_rodrigues[4],  1e-8);
    check("final 0121 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_speu_rodrigues[5],  1e-8);
    check("final 0121 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_speu_rodrigues[6],  1e-8);
    check("final 0121 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_speu_rodrigues[7],  1e-8);
    check("final 0121 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_speu_rodrigues[8],  1e-8);
    check("final 0121 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_speu_rodrigues[9],  1e-8);
    check("final 0121 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_speu_rodrigues[10], 1e-8);
    check("final 0121 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_speu_rodrigues[11], 1e-8);
    check("final 0121 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_speu_rodrigues[12], 1e-8);
    check("final 0121 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_speu_rodrigues[13], 1e-6);
    check("final 0121 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_speu_rodrigues[14], 1e-6);
    check("final 0121 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_speu_rodrigues[15], 1e-6);
    check("final 0121 body yaw [deg]   ", res_space_rodrigues[16], res_space_speu_rodrigues[16], 1e-6);
    check("final 0121 body pitch [deg] ", res_space_rodrigues[17], res_space_speu_rodrigues[17], 1e-6);
    check("final 0121 body bank [deg]  ", res_space_rodrigues[18], res_space_speu_rodrigues[18], 1e-6);
    check("final 0121 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_speu_rodrigues[19], 1e-6);
    check("final 0121 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_speu_rodrigues[20], 1e-6);
    check("final 0121 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_speu_rodrigues[21], 1e-6);
    check("final 0121 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_speu_rodrigues[22], 1e-8);
    check("final 0121 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_speu_rodrigues[23], 1e-8);
    check("final 0121 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_speu_rodrigues[24], 1e-8);

    check("final 0122 t [sec]          ", res_space_rodrigues[0],  res_space_speu_dcm[0],  1e-8);
    check("final 0122 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_speu_dcm[1],  1e-8);
    check("final 0122 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_speu_dcm[2],  1e-8);
    check("final 0122 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_speu_dcm[3],  1e-8);
    check("final 0122 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_speu_dcm[4],  1e-8);
    check("final 0122 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_speu_dcm[5],  1e-8);
    check("final 0122 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_speu_dcm[6],  1e-8);
    check("final 0122 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_speu_dcm[7],  1e-8);
    check("final 0122 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_speu_dcm[8],  1e-8);
    check("final 0122 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_speu_dcm[9],  1e-8);
    check("final 0122 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_speu_dcm[10], 1e-8);
    check("final 0122 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_speu_dcm[11], 1e-8);
    check("final 0122 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_speu_dcm[12], 1e-8);
    check("final 0122 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_speu_dcm[13], 1e-6);
    check("final 0122 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_speu_dcm[14], 1e-6);
    check("final 0122 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_speu_dcm[15], 1e-6);
    check("final 0122 body yaw [deg]   ", res_space_rodrigues[16], res_space_speu_dcm[16], 1e-6);
    check("final 0122 body pitch [deg] ", res_space_rodrigues[17], res_space_speu_dcm[17], 1e-6);
    check("final 0122 body bank [deg]  ", res_space_rodrigues[18], res_space_speu_dcm[18], 1e-6);
    check("final 0122 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_speu_dcm[19], 1e-6);
    check("final 0122 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_speu_dcm[20], 1e-6);
    check("final 0122 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_speu_dcm[21], 1e-6);
    check("final 0122 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_speu_dcm[22], 1e-8);
    check("final 0122 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_speu_dcm[23], 1e-8);
    check("final 0122 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_speu_dcm[24], 1e-8);

    check("final 0123 t [sec]          ", res_space_rodrigues[0],  res_space_homogeneous[0],  1e-8);
    check("final 0123 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_homogeneous[1],  1e-8);
    check("final 0123 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_homogeneous[2],  1e-8);
    check("final 0123 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_homogeneous[3],  1e-8);
    check("final 0123 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_homogeneous[4],  1e-8);
    check("final 0123 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_homogeneous[5],  1e-8);
    check("final 0123 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_homogeneous[6],  1e-8);
    check("final 0123 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_homogeneous[7],  1e-8);
    check("final 0123 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_homogeneous[8],  1e-8);
    check("final 0123 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_homogeneous[9],  1e-8);
    check("final 0123 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_homogeneous[10], 1e-8);
    check("final 0123 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_homogeneous[11], 1e-8);
    check("final 0123 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_homogeneous[12], 1e-8);
    check("final 0123 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_homogeneous[13], 1e-6);
    check("final 0123 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_homogeneous[14], 1e-6);
    check("final 0123 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_homogeneous[15], 1e-6);
    check("final 0123 body yaw [deg]   ", res_space_rodrigues[16], res_space_homogeneous[16], 1e-6);
    check("final 0123 body pitch [deg] ", res_space_rodrigues[17], res_space_homogeneous[17], 1e-6);
    check("final 0123 body bank [deg]  ", res_space_rodrigues[18], res_space_homogeneous[18], 1e-6);
    check("final 0123 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_homogeneous[19], 1e-6);
    check("final 0123 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_homogeneous[20], 1e-6);
    check("final 0123 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_homogeneous[21], 1e-6);
    check("final 0123 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_homogeneous[22], 1e-8);
    check("final 0123 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_homogeneous[23], 1e-8);
    check("final 0123 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_homogeneous[24], 1e-8);

    check("final 0124 t [sec]          ", res_space_rodrigues[0],  res_space_trfv[0],  1e-8);
    check("final 0124 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_trfv[1],  1e-8);
    check("final 0124 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_trfv[2],  1e-8);
    check("final 0124 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_trfv[3],  1e-8);
    check("final 0124 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_trfv[4],  1e-8);
    check("final 0124 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_trfv[5],  1e-8);
    check("final 0124 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_trfv[6],  1e-8);
    check("final 0124 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_trfv[7],  1e-8);
    check("final 0124 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_trfv[8],  1e-8);
    check("final 0124 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_trfv[9],  1e-8);
    check("final 0124 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_trfv[10], 1e-8);
    check("final 0124 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_trfv[11], 1e-8);
    check("final 0124 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_trfv[12], 1e-8);
    check("final 0124 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_trfv[13], 1e-6);
    check("final 0124 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_trfv[14], 1e-6);
    check("final 0124 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_trfv[15], 1e-6);
    check("final 0124 body yaw [deg]   ", res_space_rodrigues[16], res_space_trfv[16], 1e-6);
    check("final 0124 body pitch [deg] ", res_space_rodrigues[17], res_space_trfv[17], 1e-6);
    check("final 0124 body bank [deg]  ", res_space_rodrigues[18], res_space_trfv[18], 1e-6);
    check("final 0124 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_trfv[19], 1e-6);
    check("final 0124 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_trfv[20], 1e-6);
    check("final 0124 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_trfv[21], 1e-6);
    check("final 0124 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_trfv[22], 1e-8);
    check("final 0124 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_trfv[23], 1e-8);
    check("final 0124 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_trfv[24], 1e-8);

    check("final 0125 t [sec]          ", res_space_rodrigues[0],  res_space_dual[0],  1e-8);
    check("final 0125 v_b1 [mps]       ", res_space_rodrigues[1],  res_space_dual[1],  1e-8);
    check("final 0125 v_b2 [mps]       ", res_space_rodrigues[2],  res_space_dual[2],  1e-8);
    check("final 0125 v_b3 [mps]       ", res_space_rodrigues[3],  res_space_dual[3],  1e-8);
    check("final 0125 v_w1 [mps]       ", res_space_rodrigues[4],  res_space_dual[4],  1e-8);
    check("final 0125 v_w2 [mps]       ", res_space_rodrigues[5],  res_space_dual[5],  1e-8);
    check("final 0125 v_w3 [mps]       ", res_space_rodrigues[6],  res_space_dual[6],  1e-8);
    check("final 0125 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_space_dual[7],  1e-8);
    check("final 0125 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_space_dual[8],  1e-8);
    check("final 0125 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_space_dual[9],  1e-8);
    check("final 0125 w_wbw1 [rps]     ", res_space_rodrigues[10], res_space_dual[10], 1e-8);
    check("final 0125 w_wbw2 [rps]     ", res_space_rodrigues[11], res_space_dual[11], 1e-8);
    check("final 0125 w_wbw3 [rps]     ", res_space_rodrigues[12], res_space_dual[12], 1e-8);
    check("final 0125 x_wbw1 [m]       ", res_space_rodrigues[13], res_space_dual[13], 1e-6);
    check("final 0125 x_wbw2 [m]       ", res_space_rodrigues[14], res_space_dual[14], 1e-6);
    check("final 0125 x_wbw3 [m]       ", res_space_rodrigues[15], res_space_dual[15], 1e-6);
    check("final 0125 body yaw [deg]   ", res_space_rodrigues[16], res_space_dual[16], 1e-6);
    check("final 0125 body pitch [deg] ", res_space_rodrigues[17], res_space_dual[17], 1e-6);
    check("final 0125 body bank [deg]  ", res_space_rodrigues[18], res_space_dual[18], 1e-6);
    check("final 0125 x_wpw1 [m]       ", res_space_rodrigues[19], res_space_dual[19], 1e-6);
    check("final 0125 x_wpw2 [m]       ", res_space_rodrigues[20], res_space_dual[20], 1e-6);
    check("final 0125 x_wpw3 [m]       ", res_space_rodrigues[21], res_space_dual[21], 1e-6);
    check("final 0125 v_wpw1 [mps]     ", res_space_rodrigues[22], res_space_dual[22], 1e-8);
    check("final 0125 v_wpw2 [mps]     ", res_space_rodrigues[23], res_space_dual[23], 1e-8);
    check("final 0125 v_wpw3 [mps]     ", res_space_rodrigues[24], res_space_dual[24], 1e-8);

    check("final 0131 t [sec]          ", res_space_rodrigues[0],  res_body_speu_rodrigues[0],  1e-8);
    check("final 0131 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_speu_rodrigues[1],  1e-8);
    check("final 0131 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_speu_rodrigues[2],  1e-8);
    check("final 0131 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_speu_rodrigues[3],  1e-8);
    check("final 0131 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_speu_rodrigues[4],  1e-8);
    check("final 0131 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_speu_rodrigues[5],  1e-8);
    check("final 0131 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_speu_rodrigues[6],  1e-8);
    check("final 0131 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_speu_rodrigues[7],  1e-8);
    check("final 0131 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_speu_rodrigues[8],  1e-8);
    check("final 0131 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_speu_rodrigues[9],  1e-8);
    check("final 0131 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_speu_rodrigues[10], 1e-8);
    check("final 0131 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_speu_rodrigues[11], 1e-8);
    check("final 0131 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_speu_rodrigues[12], 1e-8);
    check("final 0131 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_speu_rodrigues[13], 1e-6);
    check("final 0131 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_speu_rodrigues[14], 1e-6);
    check("final 0131 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_speu_rodrigues[15], 1e-6);
    check("final 0131 body yaw [deg]   ", res_space_rodrigues[16], res_body_speu_rodrigues[16], 1e-6);
    check("final 0131 body pitch [deg] ", res_space_rodrigues[17], res_body_speu_rodrigues[17], 1e-6);
    check("final 0131 body bank [deg]  ", res_space_rodrigues[18], res_body_speu_rodrigues[18], 1e-6);
    check("final 0131 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_speu_rodrigues[19], 1e-6);
    check("final 0131 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_speu_rodrigues[20], 1e-6);
    check("final 0131 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_speu_rodrigues[21], 1e-6);
    check("final 0131 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_speu_rodrigues[22], 1e-8);
    check("final 0131 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_speu_rodrigues[23], 1e-8);
    check("final 0131 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_speu_rodrigues[24], 1e-8);

    check("final 0132 t [sec]          ", res_space_rodrigues[0],  res_body_speu_dcm[0],  1e-8);
    check("final 0132 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_speu_dcm[1],  1e-8);
    check("final 0132 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_speu_dcm[2],  1e-8);
    check("final 0132 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_speu_dcm[3],  1e-8);
    check("final 0132 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_speu_dcm[4],  1e-8);
    check("final 0132 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_speu_dcm[5],  1e-8);
    check("final 0132 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_speu_dcm[6],  1e-8);
    check("final 0132 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_speu_dcm[7],  1e-8);
    check("final 0132 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_speu_dcm[8],  1e-8);
    check("final 0132 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_speu_dcm[9],  1e-8);
    check("final 0132 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_speu_dcm[10], 1e-8);
    check("final 0132 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_speu_dcm[11], 1e-8);
    check("final 0132 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_speu_dcm[12], 1e-8);
    check("final 0132 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_speu_dcm[13], 1e-6);
    check("final 0132 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_speu_dcm[14], 1e-6);
    check("final 0132 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_speu_dcm[15], 1e-6);
    check("final 0132 body yaw [deg]   ", res_space_rodrigues[16], res_body_speu_dcm[16], 1e-6);
    check("final 0132 body pitch [deg] ", res_space_rodrigues[17], res_body_speu_dcm[17], 1e-6);
    check("final 0132 body bank [deg]  ", res_space_rodrigues[18], res_body_speu_dcm[18], 1e-6);
    check("final 0132 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_speu_dcm[19], 1e-6);
    check("final 0132 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_speu_dcm[20], 1e-6);
    check("final 0132 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_speu_dcm[21], 1e-6);
    check("final 0132 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_speu_dcm[22], 1e-8);
    check("final 0132 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_speu_dcm[23], 1e-8);
    check("final 0132 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_speu_dcm[24], 1e-8);

    check("final 0133 t [sec]          ", res_space_rodrigues[0],  res_body_homogeneous[0],  1e-8);
    check("final 0133 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_homogeneous[1],  1e-8);
    check("final 0133 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_homogeneous[2],  1e-8);
    check("final 0133 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_homogeneous[3],  1e-8);
    check("final 0133 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_homogeneous[4],  1e-8);
    check("final 0133 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_homogeneous[5],  1e-8);
    check("final 0133 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_homogeneous[6],  1e-8);
    check("final 0133 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_homogeneous[7],  1e-8);
    check("final 0133 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_homogeneous[8],  1e-8);
    check("final 0133 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_homogeneous[9],  1e-8);
    check("final 0133 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_homogeneous[10], 1e-8);
    check("final 0133 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_homogeneous[11], 1e-8);
    check("final 0133 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_homogeneous[12], 1e-8);
    check("final 0133 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_homogeneous[13], 1e-6);
    check("final 0133 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_homogeneous[14], 1e-6);
    check("final 0133 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_homogeneous[15], 1e-6);
    check("final 0133 body yaw [deg]   ", res_space_rodrigues[16], res_body_homogeneous[16], 1e-6);
    check("final 0133 body pitch [deg] ", res_space_rodrigues[17], res_body_homogeneous[17], 1e-6);
    check("final 0133 body bank [deg]  ", res_space_rodrigues[18], res_body_homogeneous[18], 1e-6);
    check("final 0133 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_homogeneous[19], 1e-6);
    check("final 0133 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_homogeneous[20], 1e-6);
    check("final 0133 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_homogeneous[21], 1e-6);
    check("final 0133 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_homogeneous[22], 1e-8);
    check("final 0133 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_homogeneous[23], 1e-8);
    check("final 0133 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_homogeneous[24], 1e-8);

    check("final 0134 t [sec]          ", res_space_rodrigues[0],  res_body_trfv[0],  1e-8);
    check("final 0134 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_trfv[1],  1e-8);
    check("final 0134 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_trfv[2],  1e-8);
    check("final 0134 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_trfv[3],  1e-8);
    check("final 0134 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_trfv[4],  1e-8);
    check("final 0134 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_trfv[5],  1e-8);
    check("final 0134 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_trfv[6],  1e-8);
    check("final 0134 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_trfv[7],  1e-8);
    check("final 0134 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_trfv[8],  1e-8);
    check("final 0134 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_trfv[9],  1e-8);
    check("final 0134 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_trfv[10], 1e-8);
    check("final 0134 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_trfv[11], 1e-8);
    check("final 0134 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_trfv[12], 1e-8);
    check("final 0134 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_trfv[13], 1e-6);
    check("final 0134 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_trfv[14], 1e-6);
    check("final 0134 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_trfv[15], 1e-6);
    check("final 0134 body yaw [deg]   ", res_space_rodrigues[16], res_body_trfv[16], 1e-6);
    check("final 0134 body pitch [deg] ", res_space_rodrigues[17], res_body_trfv[17], 1e-6);
    check("final 0134 body bank [deg]  ", res_space_rodrigues[18], res_body_trfv[18], 1e-6);
    check("final 0134 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_trfv[19], 1e-6);
    check("final 0134 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_trfv[20], 1e-6);
    check("final 0134 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_trfv[21], 1e-6);
    check("final 0134 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_trfv[22], 1e-8);
    check("final 0134 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_trfv[23], 1e-8);
    check("final 0134 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_trfv[24], 1e-8);

    check("final 0135 t [sec]          ", res_space_rodrigues[0],  res_body_dual[0],  1e-8);
    check("final 0135 v_b1 [mps]       ", res_space_rodrigues[1],  res_body_dual[1],  1e-8);
    check("final 0135 v_b2 [mps]       ", res_space_rodrigues[2],  res_body_dual[2],  1e-8);
    check("final 0135 v_b3 [mps]       ", res_space_rodrigues[3],  res_body_dual[3],  1e-8);
    check("final 0135 v_w1 [mps]       ", res_space_rodrigues[4],  res_body_dual[4],  1e-8);
    check("final 0135 v_w2 [mps]       ", res_space_rodrigues[5],  res_body_dual[5],  1e-8);
    check("final 0135 v_w3 [mps]       ", res_space_rodrigues[6],  res_body_dual[6],  1e-8);
    check("final 0135 w_wbb1 [rps]     ", res_space_rodrigues[7],  res_body_dual[7],  1e-8);
    check("final 0135 w_wbb2 [rps]     ", res_space_rodrigues[8],  res_body_dual[8],  1e-8);
    check("final 0135 w_wbb3 [rps]     ", res_space_rodrigues[9],  res_body_dual[9],  1e-8);
    check("final 0135 w_wbw1 [rps]     ", res_space_rodrigues[10], res_body_dual[10], 1e-8);
    check("final 0135 w_wbw2 [rps]     ", res_space_rodrigues[11], res_body_dual[11], 1e-8);
    check("final 0135 w_wbw3 [rps]     ", res_space_rodrigues[12], res_body_dual[12], 1e-8);
    check("final 0135 x_wbw1 [m]       ", res_space_rodrigues[13], res_body_dual[13], 1e-6);
    check("final 0135 x_wbw2 [m]       ", res_space_rodrigues[14], res_body_dual[14], 1e-6);
    check("final 0135 x_wbw3 [m]       ", res_space_rodrigues[15], res_body_dual[15], 1e-6);
    check("final 0135 body yaw [deg]   ", res_space_rodrigues[16], res_body_dual[16], 1e-6);
    check("final 0135 body pitch [deg] ", res_space_rodrigues[17], res_body_dual[17], 1e-6);
    check("final 0135 body bank [deg]  ", res_space_rodrigues[18], res_body_dual[18], 1e-6);
    check("final 0135 x_wpw1 [m]       ", res_space_rodrigues[19], res_body_dual[19], 1e-6);
    check("final 0135 x_wpw2 [m]       ", res_space_rodrigues[20], res_body_dual[20], 1e-6);
    check("final 0135 x_wpw3 [m]       ", res_space_rodrigues[21], res_body_dual[21], 1e-6);
    check("final 0135 v_wpw1 [mps]     ", res_space_rodrigues[22], res_body_dual[22], 1e-8);
    check("final 0135 v_wpw2 [mps]     ", res_space_rodrigues[23], res_body_dual[23], 1e-8);
    check("final 0135 v_wpw3 [mps]     ", res_space_rodrigues[24], res_body_dual[24], 1e-8);

    finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_rodrigues(vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::rodrigues& q0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    vector<double>          t_sec(nel);
    vector<Eigen::Vector3d> v_w_mps(nel);
    vector<ang::so3_tangent> w_wbw_rps(nel);
    vector<Eigen::Vector3d> x_wbw_m(nel);
    vector<ang::rodrigues>  q_wb(nel);

    t_sec[0]     = 0.;
    v_w_mps[0]   = v0_w_mps;
    w_wbw_rps[0] = w0_wbw_rps;
    x_wbw_m[0]   = x0_wbw_m;
    q_wb[0]      = q0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tw_wbw_rps, Tx_wbw_m;
    ang::rodrigues  Tq_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_w_mps_dt   = q_wb[i-1] * a_b_mps2;
        Tdw_wbw_rps_dt = q_wb[i-1] * alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps[i-1];

        Tv_w_mps       = v_w_mps[i-1]   + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps[i-1]() + Tdw_wbw_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        Tq_wb          = ang::rotv(w_wbw_rps[i-1]() * Deltat_sec).exp_map_rodrigues() * q_wb[i-1];

        Rdv_w_mps_dt   = Tq_wb * a_b_mps2;
        Rdw_wbw_rps_dt = Tq_wb * alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_w_mps[i]     = v_w_mps[i-1]   + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps[i]   = w_wbw_rps[i-1]() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        q_wb[i]        = ang::rotv((w_wbw_rps[i-1]() + Tw_wbw_rps) * 0.5 * Deltat_sec).exp_map_rodrigues() * q_wb[i-1];
    }

    ang::euler euler_wb_rad_back(q_wb.back());
    Eigen::Vector3d v_b_mps_back   = q_wb.back() / v_w_mps.back();
    ang::so3_tangent w_wbb_rps_back = q_wb.back() % w_wbw_rps.back();
    Eigen::Vector3d x_wpw_m_back   = x_wbw_m.back() + q_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back = v_w_mps.back() + w_wbw_rps.back().point_velocity(q_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps.back()(0);
    res[5]  = v_w_mps.back()(1);
    res[6]  = v_w_mps.back()(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps.back()()(0);
    res[11] = w_wbw_rps.back()()(1);
    res[12] = w_wbw_rps.back()()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_rodrigues integration." << std::endl;
} // closes test_space_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_dcm(vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::dcm& R0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    vector<double>          t_sec(nel);
    vector<Eigen::Vector3d> v_w_mps(nel);
    vector<ang::so3_tangent> w_wbw_rps(nel);
    vector<Eigen::Vector3d> x_wbw_m(nel);
    vector<ang::dcm>        R_wb(nel);

    t_sec[0]     = 0.;
    v_w_mps[0]   = v0_w_mps;
    w_wbw_rps[0] = w0_wbw_rps;
    x_wbw_m[0]   = x0_wbw_m;
    R_wb[0]      = R0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tw_wbw_rps, Tx_wbw_m;
    ang::dcm TR_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_w_mps_dt   = R_wb[i-1] * a_b_mps2;
        Tdw_wbw_rps_dt = R_wb[i-1] * alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps[i-1];

        Tv_w_mps       = v_w_mps[i-1]   + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps[i-1]() + Tdw_wbw_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        TR_wb          = ang::rotv(w_wbw_rps[i-1]() * Deltat_sec).exp_map_dcm() * R_wb[i-1];

        Rdv_w_mps_dt   = TR_wb * a_b_mps2;
        Rdw_wbw_rps_dt = TR_wb * alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_w_mps[i]     = v_w_mps[i-1]   + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps[i]   = w_wbw_rps[i-1]() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        R_wb[i]        = ang::rotv((w_wbw_rps[i-1]() + Tw_wbw_rps) * 0.5 * Deltat_sec).exp_map_dcm() * R_wb[i-1];
    }

    ang::euler euler_wb_rad_back(R_wb.back());
    Eigen::Vector3d v_b_mps_back    = R_wb.back() / v_w_mps.back();
    ang::so3_tangent w_wbb_rps_back = R_wb.back() % w_wbw_rps.back();
    Eigen::Vector3d x_wpw_m_back    = x_wbw_m.back() + R_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back  = v_w_mps.back() + w_wbw_rps.back().point_velocity(R_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps.back()(0);
    res[5]  = v_w_mps.back()(1);
    res[6]  = v_w_mps.back()(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps.back()()(0);
    res[11] = w_wbw_rps.back()()(1);
    res[12] = w_wbw_rps.back()()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_dcm integration." << std::endl;
} // closes test_space_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_rotv(vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::rotv& rv0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    vector<double>          t_sec(nel);
    vector<Eigen::Vector3d> v_w_mps(nel);
    vector<ang::so3_tangent> w_wbw_rps(nel);
    vector<Eigen::Vector3d> x_wbw_m(nel);
    vector<ang::rotv>       rv_wb(nel);

    t_sec[0]     = 0.;
    v_w_mps[0]   = v0_w_mps;
    w_wbw_rps[0] = w0_wbw_rps;
    x_wbw_m[0]   = x0_wbw_m;
    rv_wb[0]     = rv0_wb;

    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tw_wbw_rps, Tx_wbw_m;
    ang::rotv Trv_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d alpha_b_rps2;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_w_mps_dt   = rv_wb[i-1] * a_b_mps2;
        Tdw_wbw_rps_dt = rv_wb[i-1] * alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps[i-1];

        Tv_w_mps       = v_w_mps[i-1]   + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps[i-1]() + Tdw_wbw_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        Trv_wb         = ang::rotv(w_wbw_rps[i-1]() * Deltat_sec) * rv_wb[i-1];

        Rdv_w_mps_dt   = Trv_wb * a_b_mps2;
        Rdw_wbw_rps_dt = Trv_wb * alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_w_mps[i]     = v_w_mps[i-1]   + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps[i]   = w_wbw_rps[i-1]() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        rv_wb[i]       = ang::rotv((w_wbw_rps[i-1]() + Tw_wbw_rps) * 0.5 * Deltat_sec) * rv_wb[i-1];
    }

    ang::euler euler_wb_rad_back(rv_wb.back());
    Eigen::Vector3d v_b_mps_back    = rv_wb.back() / v_w_mps.back();
    ang::so3_tangent w_wbb_rps_back = rv_wb.back() % w_wbw_rps.back();
    Eigen::Vector3d x_wpw_m_back    = x_wbw_m.back() + rv_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back  = v_w_mps.back() + w_wbw_rps.back().point_velocity(rv_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps.back()(0);
    res[5]  = v_w_mps.back()(1);
    res[6]  = v_w_mps.back()(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps.back()()(0);
    res[11] = w_wbw_rps.back()()(1);
    res[12] = w_wbw_rps.back()()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_rotv integration." << std::endl;
} // closes test_space_rotv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_rodrigues(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::rodrigues& q0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>          t_sec(nel);
    std::vector<Eigen::Vector3d> v_b_mps(nel);
    std::vector<ang::so3_tangent> w_wbb_rps(nel);
    std::vector<Eigen::Vector3d> x_wbw_m(nel);
    std::vector<ang::rodrigues>  q_wb(nel);

    t_sec[0]     = 0.;
    v_b_mps[0]   = v0_b_mps;
    w_wbb_rps[0] = w0_wbb_rps;
    x_wbw_m[0]   = x0_wbw_m;
    q_wb[0]      = q0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_b_mps_dt, Tdw_wbb_rps_dt;
    Eigen::Vector3d Tv_b_mps, Tw_wbb_rps, Tx_wbw_m;
    ang::rodrigues  Tq_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_b_mps_dt, Rdw_wbb_rps_dt;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_b_mps_dt   = - w_wbb_rps[i-1]().cross(v_b_mps[i-1]) + a_b_mps2;
        Tdw_wbb_rps_dt = alpha_b_rps2;
        Tdx_wbw_m_dt   = q_wb[i-1] * v_b_mps[i-1];

        Tv_b_mps       = v_b_mps[i-1]   + Tdv_b_mps_dt   * Deltat_sec;
        Tw_wbb_rps     = w_wbb_rps[i-1]() + Tdw_wbb_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        Tq_wb          = q_wb[i-1] * ang::rotv(w_wbb_rps[i-1]() * Deltat_sec).exp_map_rodrigues();

        Rdv_b_mps_dt   = - Tw_wbb_rps.cross(Tv_b_mps) + a_b_mps2;
        Rdw_wbb_rps_dt = alpha_b_rps2;
        Rdx_wbw_m_dt   = Tq_wb * Tv_b_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_b_mps[i]     = v_b_mps[i-1]   + (Tdv_b_mps_dt   + Rdv_b_mps_dt)   * 0.5 * Deltat_sec;
        w_wbb_rps[i]   = w_wbb_rps[i-1]() + (Tdw_wbb_rps_dt + Rdw_wbb_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        q_wb[i]        = q_wb[i-1] * ang::rotv((w_wbb_rps[i-1]() + Tw_wbb_rps) * 0.5 * Deltat_sec).exp_map_rodrigues();
    }

    ang::euler euler_wb_rad_back(q_wb.back());
    Eigen::Vector3d v_w_mps_back    = q_wb.back() * v_b_mps.back();
    ang::so3_tangent w_wbw_rps_back = q_wb.back() | w_wbb_rps.back();
    Eigen::Vector3d x_wpw_m_back    = x_wbw_m.back() + q_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back  = v_w_mps_back + w_wbw_rps_back.point_velocity(q_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps.back()(0);
    res[2]  = v_b_mps.back()(1);
    res[3]  = v_b_mps.back()(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps.back()()(0);
    res[8]  = w_wbb_rps.back()()(1);
    res[9]  = w_wbb_rps.back()()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_rodrigues integration." << std::endl;
} // closes test_body_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_dcm(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::dcm& R0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>          t_sec(nel);
    std::vector<Eigen::Vector3d> v_b_mps(nel);
    std::vector<ang::so3_tangent> w_wbb_rps(nel);
    std::vector<Eigen::Vector3d> x_wbw_m(nel);
    std::vector<ang::dcm>        R_wb(nel);

    t_sec[0]     = 0.;
    v_b_mps[0]   = v0_b_mps;
    w_wbb_rps[0] = w0_wbb_rps;
    x_wbw_m[0]   = x0_wbw_m;
    R_wb[0]      = R0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_b_mps_dt, Tdw_wbb_rps_dt;
    Eigen::Vector3d Tv_b_mps, Tw_wbb_rps, Tx_wbw_m;
    ang::dcm        TR_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_b_mps_dt, Rdw_wbb_rps_dt;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_b_mps_dt   = - w_wbb_rps[i-1]().cross(v_b_mps[i-1]) + a_b_mps2;
        Tdw_wbb_rps_dt = alpha_b_rps2;
        Tdx_wbw_m_dt   = R_wb[i-1] * v_b_mps[i-1];

        Tv_b_mps       = v_b_mps[i-1]   + Tdv_b_mps_dt   * Deltat_sec;
        Tw_wbb_rps     = w_wbb_rps[i-1]() + Tdw_wbb_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        TR_wb          = R_wb[i-1] * ang::rotv(w_wbb_rps[i-1]() * Deltat_sec).exp_map_dcm();

        Rdv_b_mps_dt   = - Tw_wbb_rps.cross(Tv_b_mps) + a_b_mps2;
        Rdw_wbb_rps_dt = alpha_b_rps2;
        Rdx_wbw_m_dt   = TR_wb * Tv_b_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_b_mps[i]     = v_b_mps[i-1]   + (Tdv_b_mps_dt   + Rdv_b_mps_dt)   * 0.5 * Deltat_sec;
        w_wbb_rps[i]   = w_wbb_rps[i-1]() + (Tdw_wbb_rps_dt + Rdw_wbb_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        R_wb[i]        = R_wb[i-1] * ang::rotv((w_wbb_rps[i-1]() + Tw_wbb_rps) * 0.5 * Deltat_sec).exp_map_dcm();
    }

    ang::euler euler_wb_rad_back(R_wb.back());
    Eigen::Vector3d v_w_mps_back    = R_wb.back() * v_b_mps.back();
    ang::so3_tangent w_wbw_rps_back = R_wb.back() | w_wbb_rps.back();
    Eigen::Vector3d x_wpw_m_back    = x_wbw_m.back() + R_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back  = v_w_mps_back + w_wbw_rps_back.point_velocity(R_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps.back()(0);
    res[2]  = v_b_mps.back()(1);
    res[3]  = v_b_mps.back()(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps.back()()(0);
    res[8]  = w_wbb_rps.back()()(1);
    res[9]  = w_wbb_rps.back()()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_dcm integration." << std::endl;
} // closes test_body_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_rotv(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const Eigen::Vector3d& x0_wbw_m,
                                     const ang::rotv& rv0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>          t_sec(nel);
    std::vector<Eigen::Vector3d> v_b_mps(nel);
    std::vector<ang::so3_tangent> w_wbb_rps(nel);
    std::vector<Eigen::Vector3d> x_wbw_m(nel);
    std::vector<ang::rotv>       rv_wb(nel);

    t_sec[0]     = 0.;
    v_b_mps[0]   = v0_b_mps;
    w_wbb_rps[0] = w0_wbb_rps;
    x_wbw_m[0]   = x0_wbw_m;
    rv_wb[0]     = rv0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_b_mps_dt, Tdw_wbb_rps_dt;
    Eigen::Vector3d Tv_b_mps, Tw_wbb_rps, Tx_wbw_m;
    ang::rotv Trv_wb;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_b_mps_dt, Rdw_wbb_rps_dt;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdv_b_mps_dt   = - w_wbb_rps[i-1]().cross(v_b_mps[i-1]) + a_b_mps2;
        Tdw_wbb_rps_dt = alpha_b_rps2;
        Tdx_wbw_m_dt   = rv_wb[i-1] * v_b_mps[i-1];

        Tv_b_mps       = v_b_mps[i-1]   + Tdv_b_mps_dt   * Deltat_sec;
        Tw_wbb_rps     = w_wbb_rps[i-1]() + Tdw_wbb_rps_dt * Deltat_sec;
        Tx_wbw_m       = x_wbw_m[i-1]   + Tdx_wbw_m_dt   * Deltat_sec;
        Trv_wb         = rv_wb[i-1] * ang::rotv(w_wbb_rps[i-1]() * Deltat_sec);

        Rdv_b_mps_dt   = - Tw_wbb_rps.cross(Tv_b_mps) + a_b_mps2;
        Rdw_wbb_rps_dt = alpha_b_rps2;
        Rdx_wbw_m_dt   = Trv_wb * Tv_b_mps;

        t_sec[i]       = t_sec[i-1]     + Deltat_sec;
        v_b_mps[i]     = v_b_mps[i-1]   + (Tdv_b_mps_dt   + Rdv_b_mps_dt)   * 0.5 * Deltat_sec;
        w_wbb_rps[i]   = w_wbb_rps[i-1]() + (Tdw_wbb_rps_dt + Rdw_wbb_rps_dt) * 0.5 * Deltat_sec;
        x_wbw_m[i]     = x_wbw_m[i-1]   + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        rv_wb[i]       = rv_wb[i-1] * ang::rotv((w_wbb_rps[i-1]() + Tw_wbb_rps) * 0.5 * Deltat_sec);
    }

    ang::euler euler_wb_rad_back(rv_wb.back());
    Eigen::Vector3d v_w_mps_back    = rv_wb.back() * v_b_mps.back();
    ang::so3_tangent w_wbw_rps_back = rv_wb.back() | w_wbb_rps.back();
    Eigen::Vector3d x_wpw_m_back    = x_wbw_m.back() + rv_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back  = v_w_mps_back + w_wbw_rps_back.point_velocity(rv_wb.back() * x_bpb_m);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps.back()(0);
    res[2]  = v_b_mps.back()(1);
    res[3]  = v_b_mps.back()(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps.back()()(0);
    res[8]  = w_wbb_rps.back()()(1);
    res[9]  = w_wbb_rps.back()()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m.back()(0);
    res[14] = x_wbw_m.back()(1);
    res[15] = x_wbw_m.back()(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_rotv integration." << std::endl;
} // closes test_body_rotv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_speu_rodrigues(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::speu_rodrigues& Gq0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>               t_sec(nel);
    std::vector<ang::se3_tangent>     xi_wbw_mrps(nel); // contains [vi_wbw_mps = v_w_mps + x_wbw_m x w_wbw_rps, w_wbw_rps]
    std::vector<ang::speu_rodrigues>  Gq_wb(nel);

    t_sec[0] = 0.;
    xi_wbw_mrps[0].set_vi(v0_w_mps + Gq0_wb.get_T().cross(w0_wbw_rps()));
    xi_wbw_mrps[0].set_w(w0_wbw_rps);
    Gq_wb[0] = Gq0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tx_wbw_m, Tvi_wbw_mps;
    ang::so3_tangent Tw_wbw_rps;
    ang::se3_tangent Txi_wbw_mrps;
    ang::speu_rodrigues TGq_wb;
    Eigen::Vector3d x_wbw_m_prev, vi_wbw_mps_prev, v_w_mps_prev;
    Eigen::Vector3d x_wbw_m_post, vi_wbw_mps_post, v_w_mps_post;
    ang::so3_tangent w_wbw_rps_prev;
    ang::so3_tangent w_wbw_rps_post;

    // WATCH OUT --> Do it this way to obtain results similar to space so3, body so3, and body se3
    // Refer to test_differences and test_analyze for the reasons why this should not be done the obvious way
    // Instead, integrate the three members of vi_wbw_mps independently
    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        x_wbw_m_prev    = Gq_wb[i-1].get_T();
        vi_wbw_mps_prev = xi_wbw_mrps[i-1].get_vi();
        w_wbw_rps_prev  = xi_wbw_mrps[i-1].get_w();
        v_w_mps_prev    = vi_wbw_mps_prev + w_wbw_rps_prev().cross(x_wbw_m_prev);

        Tdv_w_mps_dt   = Gq_wb[i-1] ^ a_b_mps2;
        Tdw_wbw_rps_dt = Gq_wb[i-1] ^ alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps_prev;

        Tx_wbw_m       = x_wbw_m_prev     + Tdx_wbw_m_dt   * Deltat_sec;
        Tv_w_mps       = v_w_mps_prev     + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps_prev() + Tdw_wbw_rps_dt * Deltat_sec;
        Tvi_wbw_mps    = Tv_w_mps + Tx_wbw_m.cross(Tw_wbw_rps());

        Txi_wbw_mrps.set_vi(Tvi_wbw_mps);
        Txi_wbw_mrps.set_w(Tw_wbw_rps);
        TGq_wb         = ang::trfv(xi_wbw_mrps[i-1]() * Deltat_sec).exp_map_speu_rodrigues() * Gq_wb[i-1];

        Rdv_w_mps_dt   = TGq_wb ^ a_b_mps2;
        Rdw_wbw_rps_dt = TGq_wb ^ alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]        = t_sec[i-1]       + Deltat_sec;
        x_wbw_m_post    = x_wbw_m_prev     + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        v_w_mps_post    = v_w_mps_prev     + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps_post  = w_wbw_rps_prev() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        vi_wbw_mps_post = v_w_mps_post + x_wbw_m_post.cross(w_wbw_rps_post());

        xi_wbw_mrps[i].set_vi(vi_wbw_mps_post);
        xi_wbw_mrps[i].set_w(w_wbw_rps_post);
        Gq_wb[i] = ang::trfv((xi_wbw_mrps[i-1]() + Txi_wbw_mrps()) * 0.5 * Deltat_sec).exp_map_speu_rodrigues() * Gq_wb[i-1];
    }

    ang::euler euler_wb_rad_back(Gq_wb.back().get_rodrigues());
    Eigen::Vector3d x_wbw_m_back      = Gq_wb.back().get_T();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps.back().get_w();
    Eigen::Vector3d vi_wbw_mps_back   = xi_wbw_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = vi_wbw_mps_back + w_wbw_rps_back().cross(x_wbw_m_back);
    Eigen::Vector3d v_b_mps_back      = Gq_wb.back() & v_w_mps_back;
    Eigen::Vector3d x_wpw_m_back      = Gq_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps.back().point_velocity(x_wpw_m_back);

    ang::se3_tangent xi_wbb_mrps_back = Gq_wb.back() % xi_wbw_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps_back.get_w();
    // IMPORTANT NOTE: Do not obtain v_b_mps_back from xi_wbb_mrps_back --> see reason in test_differences and test_analyze

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_speu_rodrigues integration." << std::endl;
} // closes test_space_speu_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_speu_dcm(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::speu_dcm& GR0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbw_mrps(nel); // contains [vi_wbw_mps = v_w_mps + x_wbw_m x w_wbw_rps, w_wbw_rps]
    std::vector<ang::speu_dcm>    GR_wb(nel);

    t_sec[0] = 0.;
    xi_wbw_mrps[0].set_vi(v0_w_mps + GR0_wb.get_T().cross(w0_wbw_rps()));
    xi_wbw_mrps[0].set_w(w0_wbw_rps);
    GR_wb[0] = GR0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tx_wbw_m, Tvi_wbw_mps;
    ang::so3_tangent Tw_wbw_rps;
    ang::se3_tangent Txi_wbw_mrps;
    ang::speu_dcm TGR_wb;
    Eigen::Vector3d x_wbw_m_prev, vi_wbw_mps_prev, v_w_mps_prev;
    Eigen::Vector3d x_wbw_m_post, vi_wbw_mps_post, v_w_mps_post;
    ang::so3_tangent w_wbw_rps_prev;
    ang::so3_tangent w_wbw_rps_post;

    // WATCH OUT --> Do it this way to obtain results similar to space so3, body so3, and body se3
    // Refer to test_differences and test_analyze for the reasons why this should not be done the obvious way
    // Instead, integrate the three members of vi_wbw_mps independently
    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        x_wbw_m_prev    = GR_wb[i-1].get_T();
        vi_wbw_mps_prev = xi_wbw_mrps[i-1].get_vi();
        w_wbw_rps_prev  = xi_wbw_mrps[i-1].get_w();
        v_w_mps_prev    = vi_wbw_mps_prev + w_wbw_rps_prev().cross(x_wbw_m_prev);

        Tdv_w_mps_dt   = GR_wb[i-1] ^ a_b_mps2;
        Tdw_wbw_rps_dt = GR_wb[i-1] ^ alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps_prev;

        Tx_wbw_m       = x_wbw_m_prev     + Tdx_wbw_m_dt   * Deltat_sec;
        Tv_w_mps       = v_w_mps_prev     + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps_prev() + Tdw_wbw_rps_dt * Deltat_sec;
        Tvi_wbw_mps    = Tv_w_mps + Tx_wbw_m.cross(Tw_wbw_rps());

        Txi_wbw_mrps.set_vi(Tvi_wbw_mps);
        Txi_wbw_mrps.set_w(Tw_wbw_rps);
        TGR_wb         = ang::trfv(xi_wbw_mrps[i-1]() * Deltat_sec).exp_map_speu_dcm() * GR_wb[i-1];

        Rdv_w_mps_dt   = TGR_wb ^ a_b_mps2;
        Rdw_wbw_rps_dt = TGR_wb ^ alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]        = t_sec[i-1]       + Deltat_sec;
        x_wbw_m_post    = x_wbw_m_prev     + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        v_w_mps_post    = v_w_mps_prev     + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps_post  = w_wbw_rps_prev() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        vi_wbw_mps_post = v_w_mps_post + x_wbw_m_post.cross(w_wbw_rps_post());

        xi_wbw_mrps[i].set_vi(vi_wbw_mps_post);
        xi_wbw_mrps[i].set_w(w_wbw_rps_post);
        GR_wb[i]        = ang::trfv((xi_wbw_mrps[i-1]() + Txi_wbw_mrps()) * 0.5 * Deltat_sec).exp_map_speu_dcm() * GR_wb[i-1];
    }

    ang::euler euler_wb_rad_back(GR_wb.back().get_dcm());
    Eigen::Vector3d x_wbw_m_back      = GR_wb.back().get_T();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps.back().get_w();
    Eigen::Vector3d vi_wbw_mps_back   = xi_wbw_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = vi_wbw_mps_back + w_wbw_rps_back().cross(x_wbw_m_back);
    Eigen::Vector3d v_b_mps_back      = GR_wb.back() & v_w_mps_back;
    Eigen::Vector3d x_wpw_m_back      = GR_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps.back().point_velocity(x_wpw_m_back);

    ang::se3_tangent xi_wbb_mrps_back = GR_wb.back() % xi_wbw_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps_back.get_w();
    // IMPORTANT NOTE: Do not obtain v_b_mps_back from xi_wbb_mrps_back --> see reason in test_differences and test_analyze

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_speu_dcm integration." << std::endl;
} // closes test_space_speu_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_homogeneous(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::homogeneous& M0_wb,
                                     const Eigen::Vector3d& v0_w_mps,
                                     const ang::so3_tangent& w0_wbw_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbw_mrps(nel); // contains [vi_wbw_mps = v_w_mps + x_wbw_m x w_wbw_rps, w_wbw_rps]
    std::vector<ang::homogeneous> M_wb(nel);

    t_sec[0] = 0.;
    xi_wbw_mrps[0].set_vi(v0_w_mps + M0_wb.get_T().cross(w0_wbw_rps()));
    xi_wbw_mrps[0].set_w(w0_wbw_rps);
    M_wb[0] = M0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tx_wbw_m, Tvi_wbw_mps;
    ang::so3_tangent Tw_wbw_rps;
    ang::se3_tangent Txi_wbw_mrps;
    ang::homogeneous TM_wb;
    Eigen::Vector3d x_wbw_m_prev, vi_wbw_mps_prev, v_w_mps_prev;
    Eigen::Vector3d x_wbw_m_post, vi_wbw_mps_post, v_w_mps_post;
    ang::so3_tangent w_wbw_rps_prev;
    ang::so3_tangent w_wbw_rps_post;

    // WATCH OUT --> Do it this way to obtain results similar to space so3, body so3, and body se3
    // Refer to test_differences and test_analyze for the reasons why this should not be done the obvious way
    // Instead, integrate the three members of vi_wbw_mps independently
    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        x_wbw_m_prev    = M_wb[i-1].get_T();
        vi_wbw_mps_prev = xi_wbw_mrps[i-1].get_vi();
        w_wbw_rps_prev  = xi_wbw_mrps[i-1].get_w();
        v_w_mps_prev    = vi_wbw_mps_prev + w_wbw_rps_prev().cross(x_wbw_m_prev);

        Tdv_w_mps_dt   = M_wb[i-1] ^ a_b_mps2;
        Tdw_wbw_rps_dt = M_wb[i-1] ^ alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps_prev;

        Tx_wbw_m       = x_wbw_m_prev     + Tdx_wbw_m_dt   * Deltat_sec;
        Tv_w_mps       = v_w_mps_prev     + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps_prev() + Tdw_wbw_rps_dt * Deltat_sec;
        Tvi_wbw_mps    = Tv_w_mps + Tx_wbw_m.cross(Tw_wbw_rps());

        Txi_wbw_mrps.set_vi(Tvi_wbw_mps);
        Txi_wbw_mrps.set_w(Tw_wbw_rps);
        TM_wb          = ang::trfv(xi_wbw_mrps[i-1]() * Deltat_sec).exp_map_homogeneous() * M_wb[i-1];

        Rdv_w_mps_dt   = TM_wb ^ a_b_mps2;
        Rdw_wbw_rps_dt = TM_wb ^ alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]        = t_sec[i-1]       + Deltat_sec;
        x_wbw_m_post    = x_wbw_m_prev     + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        v_w_mps_post    = v_w_mps_prev     + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps_post  = w_wbw_rps_prev() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        vi_wbw_mps_post = v_w_mps_post + x_wbw_m_post.cross(w_wbw_rps_post());

        xi_wbw_mrps[i].set_vi(vi_wbw_mps_post);
        xi_wbw_mrps[i].set_w(w_wbw_rps_post);
        M_wb[i]        = ang::trfv((xi_wbw_mrps[i-1]() + Txi_wbw_mrps()) * 0.5 * Deltat_sec).exp_map_homogeneous() * M_wb[i-1];
    }

    ang::euler euler_wb_rad_back(M_wb.back().get_dcm());
    Eigen::Vector3d x_wbw_m_back      = M_wb.back().get_T();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps.back().get_w();
    Eigen::Vector3d vi_wbw_mps_back   = xi_wbw_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = vi_wbw_mps_back + w_wbw_rps_back().cross(x_wbw_m_back);
    Eigen::Vector3d v_b_mps_back      = M_wb.back() & v_w_mps_back;
    Eigen::Vector3d x_wpw_m_back      = M_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps.back().point_velocity( x_wpw_m_back);

    ang::se3_tangent xi_wbb_mrps_back = M_wb.back() % xi_wbw_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps_back.get_w();
    // IMPORTANT NOTE: Do not obtain v_b_mps_back from xi_wbb_mrps_back --> see reason in test_differences and test_analyze

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_homogeneous integration." << std::endl;
} // closes test_space_homogeneous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_trfv(std::vector<double>& res,
                                            const double& Deltat_sec,
                                            const double& t_sec_end,
                                            const ang::trfv& tau0_wb,
                                            const Eigen::Vector3d& v0_w_mps,
                                            const ang::so3_tangent& w0_wbw_rps,
                                            const Eigen::Vector3d& x_bpb_m,
                                            const Eigen::Vector3d& a_b_mps2,
                                            const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbw_mrps(nel); // contains [vi_wbw_mps = v_w_mps + x_wbw_m x w_wbw_rps, w_wbw_rps]
    std::vector<ang::trfv>        tau_wb(nel);

    t_sec[0] = 0.;
    xi_wbw_mrps[0].set_vi(v0_w_mps + tau0_wb.get_T().cross(w0_wbw_rps()));
    xi_wbw_mrps[0].set_w(w0_wbw_rps);
    tau_wb[0] = tau0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tx_wbw_m, Tvi_wbw_mps;
    ang::so3_tangent Tw_wbw_rps;
    ang::se3_tangent Txi_wbw_mrps;
    ang::trfv Ttau_wb;
    Eigen::Vector3d x_wbw_m_prev, vi_wbw_mps_prev, v_w_mps_prev;
    Eigen::Vector3d x_wbw_m_post, vi_wbw_mps_post, v_w_mps_post;
    ang::so3_tangent w_wbw_rps_prev;
    ang::so3_tangent w_wbw_rps_post;

    // WATCH OUT --> Do it this way to obtain results similar to space so3, body so3, and body se3
    // Refer to test_differences and test_analyze for the reasons why this should not be done the obvious way
    // Instead, integrate the three members of vi_wbw_mps independently
    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        x_wbw_m_prev    = tau_wb[i-1].get_T();
        vi_wbw_mps_prev = xi_wbw_mrps[i-1].get_vi();
        w_wbw_rps_prev  = xi_wbw_mrps[i-1].get_w();
        v_w_mps_prev    = vi_wbw_mps_prev + w_wbw_rps_prev().cross(x_wbw_m_prev);

        Tdv_w_mps_dt   = tau_wb[i-1] ^ a_b_mps2;
        Tdw_wbw_rps_dt = tau_wb[i-1] ^ alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps_prev;

        Tx_wbw_m       = x_wbw_m_prev     + Tdx_wbw_m_dt   * Deltat_sec;
        Tv_w_mps       = v_w_mps_prev     + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps_prev() + Tdw_wbw_rps_dt * Deltat_sec;
        Tvi_wbw_mps    = Tv_w_mps + Tx_wbw_m.cross(Tw_wbw_rps());

        Txi_wbw_mrps.set_vi(Tvi_wbw_mps);
        Txi_wbw_mrps.set_w(Tw_wbw_rps);
        Ttau_wb          = ang::trfv(xi_wbw_mrps[i-1]() * Deltat_sec) * tau_wb[i-1];

        Rdv_w_mps_dt   = Ttau_wb ^ a_b_mps2;
        Rdw_wbw_rps_dt = Ttau_wb ^ alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]        = t_sec[i-1]       + Deltat_sec;
        x_wbw_m_post    = x_wbw_m_prev     + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        v_w_mps_post    = v_w_mps_prev     + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps_post  = w_wbw_rps_prev() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        vi_wbw_mps_post = v_w_mps_post + x_wbw_m_post.cross(w_wbw_rps_post());

        xi_wbw_mrps[i].set_vi(vi_wbw_mps_post);
        xi_wbw_mrps[i].set_w(w_wbw_rps_post);
        tau_wb[i]       = ang::trfv((xi_wbw_mrps[i-1]() + Txi_wbw_mrps()) * 0.5 * Deltat_sec) * tau_wb[i-1];
    }

    ang::euler euler_wb_rad_back(tau_wb.back().get_rotv());
    Eigen::Vector3d x_wbw_m_back      = tau_wb.back().get_T();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps.back().get_w();
    Eigen::Vector3d vi_wbw_mps_back   = xi_wbw_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = vi_wbw_mps_back + w_wbw_rps_back().cross(x_wbw_m_back);
    Eigen::Vector3d v_b_mps_back      = tau_wb.back() & v_w_mps_back;
    Eigen::Vector3d x_wpw_m_back      = tau_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps.back().point_velocity(x_wpw_m_back);

    ang::se3_tangent xi_wbb_mrps_back = tau_wb.back() % xi_wbw_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps_back.get_w();
    // IMPORTANT NOTE: Do not obtain v_b_mps_back from xi_wbb_mrps_back --> see reason in test_differences and test_analyze

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_trfv integration." << std::endl;
} // closes test_space_trfv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_space_dual(std::vector<double>& res,
                                              const double& Deltat_sec,
                                              const double& t_sec_end,
                                              const ang::dual& Z0_wb,
                                              const Eigen::Vector3d& v0_w_mps,
                                              const ang::so3_tangent& w0_wbw_rps,
                                              const Eigen::Vector3d& x_bpb_m,
                                              const Eigen::Vector3d& a_b_mps2,
                                              const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbw_mrps(nel); // contains [vi_wbw_mps = v_w_mps + x_wbw_m x w_wbw_rps, w_wbw_rps]
    std::vector<ang::dual> Z_wb(nel);

    t_sec[0] = 0.;
    xi_wbw_mrps[0].set_vi(v0_w_mps + Z0_wb.get_T().cross(w0_wbw_rps()));
    xi_wbw_mrps[0].set_w(w0_wbw_rps);
    Z_wb[0] = Z0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector3d Tdx_wbw_m_dt, Tdv_w_mps_dt, Tdw_wbw_rps_dt;
    Eigen::Vector3d Rdx_wbw_m_dt, Rdv_w_mps_dt, Rdw_wbw_rps_dt;
    Eigen::Vector3d Tv_w_mps, Tx_wbw_m, Tvi_wbw_mps;
    ang::so3_tangent Tw_wbw_rps;
    ang::se3_tangent Txi_wbw_mrps;
    ang::dual TZ_wb;
    Eigen::Vector3d x_wbw_m_prev, vi_wbw_mps_prev, v_w_mps_prev;
    Eigen::Vector3d x_wbw_m_post, vi_wbw_mps_post, v_w_mps_post;
    ang::so3_tangent w_wbw_rps_prev;
    ang::so3_tangent w_wbw_rps_post;

    // WATCH OUT --> Do it this way to obtain results similar to space so3, body so3, and body se3
    // Refer to test_differences and test_analyze for the reasons why this should not be done the obvious way
    // Instead, integrate the three members of vi_wbw_mps independently
    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        x_wbw_m_prev    = Z_wb[i-1].get_T();
        vi_wbw_mps_prev = xi_wbw_mrps[i-1].get_vi();
        w_wbw_rps_prev  = xi_wbw_mrps[i-1].get_w();
        v_w_mps_prev    = vi_wbw_mps_prev + w_wbw_rps_prev().cross(x_wbw_m_prev);

        Tdv_w_mps_dt   = Z_wb[i-1] ^ a_b_mps2;
        Tdw_wbw_rps_dt = Z_wb[i-1] ^ alpha_b_rps2;
        Tdx_wbw_m_dt   = v_w_mps_prev;

        Tx_wbw_m       = x_wbw_m_prev     + Tdx_wbw_m_dt   * Deltat_sec;
        Tv_w_mps       = v_w_mps_prev     + Tdv_w_mps_dt   * Deltat_sec;
        Tw_wbw_rps     = w_wbw_rps_prev() + Tdw_wbw_rps_dt * Deltat_sec;
        Tvi_wbw_mps    = Tv_w_mps + Tx_wbw_m.cross(Tw_wbw_rps());

        Txi_wbw_mrps.set_vi(Tvi_wbw_mps);
        Txi_wbw_mrps.set_w(Tw_wbw_rps);
        TZ_wb          = ang::trfv(xi_wbw_mrps[i-1]() * Deltat_sec).exp_map_dual() * Z_wb[i-1];

        Rdv_w_mps_dt   = TZ_wb ^ a_b_mps2;
        Rdw_wbw_rps_dt = TZ_wb ^ alpha_b_rps2;
        Rdx_wbw_m_dt   = Tv_w_mps;

        t_sec[i]        = t_sec[i-1]       + Deltat_sec;
        x_wbw_m_post    = x_wbw_m_prev     + (Tdx_wbw_m_dt   + Rdx_wbw_m_dt)   * 0.5 * Deltat_sec;
        v_w_mps_post    = v_w_mps_prev     + (Tdv_w_mps_dt   + Rdv_w_mps_dt)   * 0.5 * Deltat_sec;
        w_wbw_rps_post  = w_wbw_rps_prev() + (Tdw_wbw_rps_dt + Rdw_wbw_rps_dt) * 0.5 * Deltat_sec;
        vi_wbw_mps_post = v_w_mps_post + x_wbw_m_post.cross(w_wbw_rps_post());

        xi_wbw_mrps[i].set_vi(vi_wbw_mps_post);
        xi_wbw_mrps[i].set_w(w_wbw_rps_post);
        Z_wb[i]        = ang::trfv((xi_wbw_mrps[i-1]() + Txi_wbw_mrps()) * 0.5 * Deltat_sec).exp_map_dual() * Z_wb[i-1];
    }

    ang::euler euler_wb_rad_back(Z_wb.back().get_rodrigues());
    Eigen::Vector3d x_wbw_m_back      = Z_wb.back().get_T();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps.back().get_w();
    Eigen::Vector3d vi_wbw_mps_back   = xi_wbw_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = vi_wbw_mps_back + w_wbw_rps_back().cross(x_wbw_m_back);
    Eigen::Vector3d v_b_mps_back      = Z_wb.back() & v_w_mps_back;
    Eigen::Vector3d x_wpw_m_back      = Z_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps.back().point_velocity( x_wpw_m_back);

    ang::se3_tangent xi_wbb_mrps_back = Z_wb.back() % xi_wbw_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps_back.get_w();
    // IMPORTANT NOTE: Do not obtain v_b_mps_back from xi_wbb_mrps_back --> see reason in test_differences and test_analyze

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded space_dual integration." << std::endl;
} // closes test_space_dual

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_speu_rodrigues(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::speu_rodrigues& Gq0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>              t_sec(nel);
    std::vector<ang::se3_tangent>    xi_wbb_mrps(nel); // contains [vi_wbbb_mps, w_wbb_rps]
    std::vector<ang::speu_rodrigues> Gq_wb(nel);

    t_sec[0] = 0.;
    xi_wbb_mrps[0].set_vi(v0_b_mps);
    xi_wbb_mrps[0].set_w(w0_wbb_rps);
    Gq_wb[0] = Gq0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector6d Tdxi_wbb_mrps_dt, Rdxi_wbb_mrps_dt;
    ang::se3_tangent Txi_wbb_mrps;
    ang::speu_rodrigues TGq_wb;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdxi_wbb_mrps_dt.head<3>() = - xi_wbb_mrps[i-1].get_w()().cross(xi_wbb_mrps[i-1].get_vi()) + a_b_mps2;
        Tdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        Txi_wbb_mrps() = xi_wbb_mrps[i-1]() + Tdxi_wbb_mrps_dt * Deltat_sec;
        TGq_wb         = Gq_wb[i-1] * ang::trfv(xi_wbb_mrps[i-1]() * Deltat_sec).exp_map_speu_rodrigues();

        Rdxi_wbb_mrps_dt.head<3>() = - Txi_wbb_mrps.get_w()().cross(Txi_wbb_mrps.get_vi()) + a_b_mps2;
        Rdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        t_sec[i]         = t_sec[i-1] + Deltat_sec;
        xi_wbb_mrps[i]() = xi_wbb_mrps[i-1]() + (Tdxi_wbb_mrps_dt + Rdxi_wbb_mrps_dt) * 0.5 * Deltat_sec;
        Gq_wb[i]         = Gq_wb[i-1] * ang::trfv((xi_wbb_mrps[i-1]() + Txi_wbb_mrps()) * 0.5 * Deltat_sec).exp_map_speu_rodrigues();
    }

    ang::euler euler_wb_rad_back(Gq_wb.back().get_rodrigues());
    Eigen::Vector3d x_wbw_m_back      = Gq_wb.back().get_T();
    ang::se3_tangent xi_wbw_mrps_back = Gq_wb.back() | xi_wbb_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps.back().get_w();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps_back.get_w();
    Eigen::Vector3d v_b_mps_back      = xi_wbb_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = Gq_wb.back() ^ v_b_mps_back;
    Eigen::Vector3d x_wpw_m_back      = Gq_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps_back.point_velocity(x_wpw_m_back);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_speu_rodrigues integration." << std::endl;
} // closes test_body_speu_rodrigues

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_speu_dcm(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::speu_dcm& GR0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbb_mrps(nel); // contains [vi_wbb_mps, w_wbb_rps]
    std::vector<ang::speu_dcm>    GR_wb(nel);

    t_sec[0] = 0.;
    xi_wbb_mrps[0].set_vi(v0_b_mps);
    xi_wbb_mrps[0].set_w(w0_wbb_rps);
    GR_wb[0] = GR0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector6d Tdxi_wbb_mrps_dt, Rdxi_wbb_mrps_dt;
    ang::se3_tangent Txi_wbb_mrps;
    ang::speu_dcm TGR_wb;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdxi_wbb_mrps_dt.head<3>() = - xi_wbb_mrps[i-1].get_w()().cross(xi_wbb_mrps[i-1].get_vi()) + a_b_mps2;
        Tdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        Txi_wbb_mrps() = xi_wbb_mrps[i-1]() + Tdxi_wbb_mrps_dt * Deltat_sec;
        TGR_wb         = GR_wb[i-1] * ang::trfv(xi_wbb_mrps[i-1]() * Deltat_sec).exp_map_speu_dcm();

        Rdxi_wbb_mrps_dt.head<3>() = - Txi_wbb_mrps.get_w()().cross(Txi_wbb_mrps.get_vi()) + a_b_mps2;
        Rdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        t_sec[i]         = t_sec[i-1] + Deltat_sec;
        xi_wbb_mrps[i]() = xi_wbb_mrps[i-1]() + (Tdxi_wbb_mrps_dt + Rdxi_wbb_mrps_dt) * 0.5 * Deltat_sec;
        GR_wb[i]         = GR_wb[i-1] * ang::trfv((xi_wbb_mrps[i-1]() + Txi_wbb_mrps()) * 0.5 * Deltat_sec).exp_map_speu_dcm();
    }

    ang::euler euler_wb_rad_back(GR_wb.back().get_dcm());
    Eigen::Vector3d x_wbw_m_back      = GR_wb.back().get_T();
    ang::se3_tangent xi_wbw_mrps_back = GR_wb.back() | xi_wbb_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps.back().get_w();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps_back.get_w();
    Eigen::Vector3d v_b_mps_back      = xi_wbb_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = GR_wb.back() ^ v_b_mps_back;
    Eigen::Vector3d x_wpw_m_back      = GR_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps_back.point_velocity(x_wpw_m_back);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_speu_dcm integration." << std::endl;
} // closes test_body_speu_dcm

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_homogeneous(std::vector<double>& res,
                                     const double& Deltat_sec,
                                     const double& t_sec_end,
                                     const ang::homogeneous& M0_wb,
                                     const Eigen::Vector3d& v0_b_mps,
                                     const ang::so3_tangent& w0_wbb_rps,
                                     const Eigen::Vector3d& x_bpb_m,
                                     const Eigen::Vector3d& a_b_mps2,
                                     const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbb_mrps(nel); // contains [vi_wbb_mps, w_wbb_rps]
    std::vector<ang::homogeneous> M_wb(nel);

    t_sec[0] = 0.;
    xi_wbb_mrps[0].set_vi(v0_b_mps);
    xi_wbb_mrps[0].set_w(w0_wbb_rps);
    M_wb[0] = M0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector6d Tdxi_wbb_mrps_dt, Rdxi_wbb_mrps_dt;
    ang::se3_tangent Txi_wbb_mrps;
    ang::homogeneous TM_wb;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdxi_wbb_mrps_dt.head<3>() = - xi_wbb_mrps[i-1].get_w()().cross(xi_wbb_mrps[i-1].get_vi()) + a_b_mps2;
        Tdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        Txi_wbb_mrps() = xi_wbb_mrps[i-1]() + Tdxi_wbb_mrps_dt * Deltat_sec;
        TM_wb          = M_wb[i-1] * ang::trfv(xi_wbb_mrps[i-1]() * Deltat_sec).exp_map_homogeneous();

        Rdxi_wbb_mrps_dt.head<3>() = - Txi_wbb_mrps.get_w()().cross(Txi_wbb_mrps.get_vi()) + a_b_mps2;
        Rdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        t_sec[i]         = t_sec[i-1] + Deltat_sec;
        xi_wbb_mrps[i]() = xi_wbb_mrps[i-1]() + (Tdxi_wbb_mrps_dt + Rdxi_wbb_mrps_dt) * 0.5 * Deltat_sec;
        M_wb[i]          = M_wb[i-1] * ang::trfv((xi_wbb_mrps[i-1]() + Txi_wbb_mrps()) * 0.5 * Deltat_sec).exp_map_homogeneous();
    }

    ang::euler euler_wb_rad_back(M_wb.back().get_dcm());
    Eigen::Vector3d x_wbw_m_back      = M_wb.back().get_T();
    ang::se3_tangent xi_wbw_mrps_back = M_wb.back() | xi_wbb_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps.back().get_w();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps_back.get_w();
    Eigen::Vector3d v_b_mps_back      = xi_wbb_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = M_wb.back() ^ v_b_mps_back;
    Eigen::Vector3d x_wpw_m_back      = M_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps_back.point_velocity(x_wpw_m_back);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_homogeneous integration." << std::endl;
} // closes test_body_homogeneous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_trfv(std::vector<double>& res,
                                                    const double& Deltat_sec,
                                                    const double& t_sec_end,
                                                    const ang::trfv& tau0_wb,
                                                    const Eigen::Vector3d& v0_b_mps,
                                                    const ang::so3_tangent& w0_wbb_rps,
                                                    const Eigen::Vector3d& x_bpb_m,
                                                    const Eigen::Vector3d& a_b_mps2,
                                                    const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbb_mrps(nel); // contains [vi_wbb_mps, w_wbb_rps]
    std::vector<ang::trfv>        tau_wb(nel);

    t_sec[0] = 0.;
    xi_wbb_mrps[0].set_vi(v0_b_mps);
    xi_wbb_mrps[0].set_w(w0_wbb_rps);
    tau_wb[0] = tau0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector6d Tdxi_wbb_mrps_dt, Rdxi_wbb_mrps_dt;
    ang::se3_tangent Txi_wbb_mrps;
    ang::trfv Ttau_wb;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdxi_wbb_mrps_dt.head<3>() = - xi_wbb_mrps[i-1].get_w()().cross(xi_wbb_mrps[i-1].get_vi()) + a_b_mps2;
        Tdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        Txi_wbb_mrps() = xi_wbb_mrps[i-1]() + Tdxi_wbb_mrps_dt * Deltat_sec;
        Ttau_wb        = tau_wb[i-1] * ang::trfv(xi_wbb_mrps[i-1]() * Deltat_sec);

        Rdxi_wbb_mrps_dt.head<3>() = - Txi_wbb_mrps.get_w()().cross(Txi_wbb_mrps.get_vi()) + a_b_mps2;
        Rdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        t_sec[i]         = t_sec[i-1]     + Deltat_sec;
        xi_wbb_mrps[i]() = xi_wbb_mrps[i-1]() + (Tdxi_wbb_mrps_dt + Rdxi_wbb_mrps_dt) * 0.5 * Deltat_sec;
        tau_wb[i]        = tau_wb[i-1] * ang::trfv((xi_wbb_mrps[i-1]() + Txi_wbb_mrps()) * 0.5 * Deltat_sec);
    }

    ang::euler euler_wb_rad_back(tau_wb.back().get_rotv());
    Eigen::Vector3d x_wbw_m_back      = tau_wb.back().get_T();
    ang::se3_tangent xi_wbw_mrps_back = tau_wb.back() | xi_wbb_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps.back().get_w();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps_back.get_w();
    Eigen::Vector3d v_b_mps_back      = xi_wbb_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = tau_wb.back() ^ v_b_mps_back;
    Eigen::Vector3d x_wpw_m_back      = tau_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps_back.point_velocity(x_wpw_m_back);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_trfv integration." << std::endl;
} // closes test_body_trfv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_body_dual(std::vector<double>& res,
                                             const double& Deltat_sec,
                                             const double& t_sec_end,
                                             const ang::dual& Z0_wb,
                                             const Eigen::Vector3d& v0_b_mps,
                                             const ang::so3_tangent& w0_wbb_rps,
                                             const Eigen::Vector3d& x_bpb_m,
                                             const Eigen::Vector3d& a_b_mps2,
                                             const Eigen::Vector3d& alpha0_b_rps2) {
    unsigned nel = (unsigned)floor(t_sec_end / Deltat_sec) + 1;
    double pi = math::constant::PI();

    std::vector<double>           t_sec(nel);
    std::vector<ang::se3_tangent> xi_wbb_mrps(nel); // contains [vi_wbb_mps, w_wbb_rps]
    std::vector<ang::dual> Z_wb(nel);

    t_sec[0] = 0.;
    xi_wbb_mrps[0].set_vi(v0_b_mps);
    xi_wbb_mrps[0].set_w(w0_wbb_rps);
    Z_wb[0] = Z0_wb;

    Eigen::Vector3d alpha_b_rps2;
    Eigen::Vector6d Tdxi_wbb_mrps_dt, Rdxi_wbb_mrps_dt;
    ang::se3_tangent Txi_wbb_mrps;
    ang::dual TZ_wb;

    for (int i = 1; i != nel; ++i) {
        alpha_b_rps2(0) = alpha0_b_rps2(0) * std::cos(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(1) = alpha0_b_rps2(1) * std::sin(i * Deltat_sec * 2 * pi / 10.0);
        alpha_b_rps2(2) = alpha0_b_rps2(2) * std::cos(i * Deltat_sec * 2 * pi /  5.0);

        Tdxi_wbb_mrps_dt.head<3>() = - xi_wbb_mrps[i-1].get_w()().cross(xi_wbb_mrps[i-1].get_vi()) + a_b_mps2;
        Tdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        Txi_wbb_mrps() = xi_wbb_mrps[i-1]() + Tdxi_wbb_mrps_dt * Deltat_sec;
        TZ_wb          = Z_wb[i-1] * ang::trfv(xi_wbb_mrps[i-1]() * Deltat_sec).exp_map_dual();

        Rdxi_wbb_mrps_dt.head<3>() = - Txi_wbb_mrps.get_w()().cross(Txi_wbb_mrps.get_vi()) + a_b_mps2;
        Rdxi_wbb_mrps_dt.tail<3>() = alpha_b_rps2;

        t_sec[i]         = t_sec[i-1] + Deltat_sec;
        xi_wbb_mrps[i]() = xi_wbb_mrps[i-1]() + (Tdxi_wbb_mrps_dt + Rdxi_wbb_mrps_dt) * 0.5 * Deltat_sec;
        Z_wb[i]          = Z_wb[i-1] * ang::trfv((xi_wbb_mrps[i-1]() + Txi_wbb_mrps()) * 0.5 * Deltat_sec).exp_map_dual();
    }

    ang::euler euler_wb_rad_back(Z_wb.back().get_rodrigues());
    Eigen::Vector3d x_wbw_m_back      = Z_wb.back().get_T();
    ang::se3_tangent xi_wbw_mrps_back = Z_wb.back() | xi_wbb_mrps.back();
    ang::so3_tangent w_wbb_rps_back   = xi_wbb_mrps.back().get_w();
    ang::so3_tangent w_wbw_rps_back   = xi_wbw_mrps_back.get_w();
    Eigen::Vector3d v_b_mps_back      = xi_wbb_mrps.back().get_vi();
    Eigen::Vector3d v_w_mps_back      = Z_wb.back() ^ v_b_mps_back;
    Eigen::Vector3d x_wpw_m_back      = Z_wb.back() * x_bpb_m;
    Eigen::Vector3d v_wpw_mps_back    = xi_wbw_mrps_back.point_velocity(x_wpw_m_back);

    res[0]  = t_sec.back();
    res[1]  = v_b_mps_back(0);
    res[2]  = v_b_mps_back(1);
    res[3]  = v_b_mps_back(2);
    res[4]  = v_w_mps_back(0);
    res[5]  = v_w_mps_back(1);
    res[6]  = v_w_mps_back(2);
    res[7]  = w_wbb_rps_back()(0);
    res[8]  = w_wbb_rps_back()(1);
    res[9]  = w_wbb_rps_back()(2);
    res[10] = w_wbw_rps_back()(0);
    res[11] = w_wbw_rps_back()(1);
    res[12] = w_wbw_rps_back()(2);
    res[13] = x_wbw_m_back(0);
    res[14] = x_wbw_m_back(1);
    res[15] = x_wbw_m_back(2);
    res[16] = euler_wb_rad_back.get_yaw_rad() * math::constant::R2D();
    res[17] = euler_wb_rad_back.get_pitch_rad() * math::constant::R2D();
    res[18] = euler_wb_rad_back.get_bank_rad() * math::constant::R2D();
    res[19] = x_wpw_m_back(0);
    res[20] = x_wpw_m_back(1);
    res[21] = x_wpw_m_back(2);
    res[22] = v_wpw_mps_back(0);
    res[23] = v_wpw_mps_back(1);
    res[24] = v_wpw_mps_back(2);

    std::cout << "Concluded body_dual integration." << std::endl;
} // closes test_body_dual

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_differences(const ang::speu_dcm& GR0_wb,
                                               const Eigen::Vector3d& v0_w_mps,
                                               const ang::so3_tangent& w0_wbw_rps,
                                               const Eigen::Vector3d& a_b_mps2,
                                               const Eigen::Vector3d& alpha_b_rps2) {
    // With a single integration step using Euler's method, this test determines the differences between the four ways of
    // integrating (SO3 and SE3, space and body). Useful to find out the right way of implementing space SE3 integration
    // if we want the results to coincide with the others). No problems with the other three (space SO3, body SO3, body SE3).

    double Deltat_sec = 0.001;

    ang::dcm         R0_wb        = GR0_wb.get_dcm();
    Eigen::Vector3d  T0_wbw_m     = GR0_wb.get_T();
    Eigen::Vector3d  v0_b_mps     = GR0_wb & v0_w_mps;
    ang::so3_tangent w0_wbb_rps   = GR0_wb.get_dcm() % w0_wbw_rps;
    Eigen::Vector3d  vi0_wbb_mps  = v0_b_mps;
    Eigen::Vector3d  T0_wbb_m     = GR0_wb / T0_wbw_m;
    ang::se3_tangent xi0_wbb_mrps(vi0_wbb_mps, w0_wbb_rps);
    ang::se3_tangent xi0_wbw_mrps = GR0_wb | xi0_wbb_mrps;
    Eigen::Vector3d  vi0_wbw_mps  = xi0_wbw_mrps.get_vi();

    ////////////////////////////////////////////////////////////////////////////////////////
    // SO3 space
    ang::so3_tangent Uw_wbw_rps       (w0_wbw_rps() + (GR0_wb ^ alpha_b_rps2) * Deltat_sec);
    Eigen::Vector3d  Uv_w_mps       = v0_w_mps + (GR0_wb ^ a_b_mps2) * Deltat_sec;
    Eigen::Vector3d  UT_wbw_m       = T0_wbw_m + v0_w_mps * Deltat_sec;
    ang::dcm         UR_wb          = ang::rotv(w0_wbw_rps() * Deltat_sec).exp_map_dcm() * R0_wb;

    ang::so3_tangent Uw_wbb_rps     = UR_wb % Uw_wbw_rps;
    Eigen::Vector3d  Uv_b_mps       = UR_wb / Uv_w_mps;
    Eigen::Vector3d  Uvi_wbw_mps    = Uv_w_mps - Uw_wbw_rps().cross(UT_wbw_m);
    ang::rotv        Ur_wb            (UR_wb);

    ////////////////////////////////////////////////////////////////////////////////////////
    // SO3 body
    ang::so3_tangent Rw_wbb_rps       (w0_wbb_rps() + alpha_b_rps2 * Deltat_sec);
    Eigen::Vector3d  Rv_b_mps       = v0_b_mps + (a_b_mps2 - w0_wbb_rps().cross(v0_b_mps)) * Deltat_sec;
    Eigen::Vector3d  RT_wbw_m       = T0_wbw_m + (R0_wb * v0_b_mps) * Deltat_sec;
    ang::dcm         RR_wb          = R0_wb * ang::rotv(w0_wbb_rps() * Deltat_sec).exp_map_dcm();

    ang::so3_tangent Rw_wbw_rps     = RR_wb | Rw_wbb_rps;
    Eigen::Vector3d  Rv_w_mps       = RR_wb * Rv_b_mps;
    Eigen::Vector3d  Rvi_wbw_mps    = Rv_w_mps - Rw_wbw_rps().cross(RT_wbw_m);
    ang::rotv        Rr_wb            (RR_wb);

    ////////////////////////////////////////////////////////////////////////////////////////
    // SE3 body
    Eigen::Vector6d  Sdxi_wbb_mrps_dt;
    Sdxi_wbb_mrps_dt.tail<3>()      = alpha_b_rps2;
    Sdxi_wbb_mrps_dt.head<3>()      = a_b_mps2 - xi0_wbb_mrps.get_w()().cross(xi0_wbb_mrps.get_vi());

    ang::se3_tangent Sxi_wbb_mrps     (xi0_wbb_mrps() + Sdxi_wbb_mrps_dt * Deltat_sec);
    ang::speu_dcm    SGR_wb         = GR0_wb * ang::trfv(xi0_wbb_mrps() * Deltat_sec).exp_map_speu_dcm();

    Eigen::Vector3d  ST_wbw_m       = SGR_wb.get_T();
    ang::rotv        Sr_wb          = SGR_wb.get_rotv();
    ang::se3_tangent Sxi_wbw_mrps   = SGR_wb | Sxi_wbb_mrps;
    ang::so3_tangent Sw_wbb_rps     = Sxi_wbb_mrps.get_w();
    ang::so3_tangent Sw_wbw_rps     = Sxi_wbw_mrps.get_w();
    Eigen::Vector3d  Sv_b_mps       = Sxi_wbb_mrps.get_vi();
    Eigen::Vector3d  Svi_wbw_mps    = Sxi_wbw_mrps.get_vi();
    Eigen::Vector3d  SAv_w_mps      = SGR_wb ^ Sv_b_mps;
    Eigen::Vector3d  SBv_w_mps      = Svi_wbw_mps + Sw_wbw_rps().cross(ST_wbw_m);

    ////////////////////////////////////////////////////////////////////////////////////////
    // SE3 space - GOOD WAY
    Eigen::Vector3d  Mv_w_mps       = v0_w_mps + (GR0_wb ^ a_b_mps2) * Deltat_sec;
    Eigen::Vector3d  MT_wbw_m       = T0_wbw_m + v0_w_mps * Deltat_sec;
    ang::so3_tangent Mw_wbw_rps      (w0_wbw_rps() + (GR0_wb ^ alpha_b_rps2) * Deltat_sec);
    Eigen::Vector3d  Mvi_wbw_mps    = Mv_w_mps - Mw_wbw_rps().cross(MT_wbw_m);

    ang::se3_tangent Mxi_wbw_mrps    (Mvi_wbw_mps, Mw_wbw_rps);
    ang::speu_dcm    MGR_wb         = ang::trfv(xi0_wbw_mrps() * Deltat_sec).exp_map_speu_dcm() * GR0_wb;

    ang::rotv        Mr_wb          = MGR_wb.get_rotv();
    ang::se3_tangent Mxi_wbb_mrps   = MGR_wb % Mxi_wbw_mrps;
    ang::so3_tangent Mw_wbb_rps     = Mxi_wbb_mrps.get_w();

    Eigen::Vector3d  MAv_w_mps      = Mvi_wbw_mps + Mw_wbw_rps().cross(MT_wbw_m); // same a Mv_w_mps
    Eigen::Vector3d  MAv_b_mps      = MGR_wb & Mv_w_mps;
    Eigen::Vector3d  MBv_b_mps      = Mxi_wbb_mrps.get_vi(); // not bad but different than space SO3 (introduces SO3 - SE3 diffs). Avoid it for commonality
    Eigen::Vector3d  MBv_w_mps      = MGR_wb ^ MBv_b_mps;  // not bad but different than space SO3 (introduces SO3 - SE3 diffs). Avoid it for commonality

    ////////////////////////////////////////////////////////////////////////////////////////
    // SE3 space - WRONG WAY or at least DIFFERENT THAN OTHERS
    // Note that in the Zxi_wbw_mrps, there is a term with a square of Deltat_sec that is lacking
    // This term is present in the other executions, and that is the reason for the differences
    // It does not mean that is better or worse, just that by doing things this way the Euler
    // and Heun's integrations work slightly differently
    // The missing term is Mdx_wbw_m_dt.cross(Mdw_wbw_rps_dt) * Deltat_sec * Deltat_sec
    Eigen::Vector6d  Zdxi_wbw_mrps_dt;
    Zdxi_wbw_mrps_dt.tail<3>()      = GR0_wb ^ alpha_b_rps2;
    Zdxi_wbw_mrps_dt.head<3>()      = (GR0_wb ^ a_b_mps2) + T0_wbw_m.cross(GR0_wb ^ alpha_b_rps2) + v0_w_mps.cross(w0_wbw_rps());

    ang::se3_tangent Zxi_wbw_mrps     (xi0_wbw_mrps() + Zdxi_wbw_mrps_dt * Deltat_sec);
    ang::speu_dcm    ZGR_wb         = ang::trfv(xi0_wbw_mrps() * Deltat_sec).exp_map_speu_dcm() * GR0_wb;

    Eigen::Vector3d  ZT_wbw_m       = ZGR_wb.get_T();
    ang::rotv        Zr_wb          = ZGR_wb.get_rotv();
    ang::so3_tangent Zw_wbw_rps     = Zxi_wbw_mrps.get_w();
    ang::se3_tangent Zxi_wbb_mrps   = ZGR_wb % Zxi_wbw_mrps;
    ang::so3_tangent Zw_wbb_rps     = Zxi_wbb_mrps.get_w();
    Eigen::Vector3d  Zvi_wbw_mps    = Zxi_wbw_mrps.get_vi();

    Eigen::Vector3d  ZAv_w_mps      = Zvi_wbw_mps + Zw_wbw_rps().cross(ZT_wbw_m);
    Eigen::Vector3d  ZAv_b_mps      = ZGR_wb & ZAv_w_mps;
    Eigen::Vector3d  ZBv_b_mps      = Zxi_wbb_mrps.get_vi();
    Eigen::Vector3d  ZBv_w_mps      = ZGR_wb ^ ZAv_b_mps;

    ////////////////////////////////////////////////////////////////////////////////////////

    cout << endl << "x_wbw_m" << endl;
    cout << "Slight difference between SO3 and SE3" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << UT_wbw_m(0) << setw(22) << UT_wbw_m(1) << setw(22) << UT_wbw_m(2) << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << RT_wbw_m(0) << setw(22) << RT_wbw_m(1) << setw(22) << RT_wbw_m(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << ST_wbw_m(0) << setw(22) << ST_wbw_m(1) << setw(22) << ST_wbw_m(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << MT_wbw_m(0) << setw(22) << MT_wbw_m(1) << setw(22) << MT_wbw_m(2) << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << ZT_wbw_m(0) << setw(22) << ZT_wbw_m(1) << setw(22) << ZT_wbw_m(2) << endl;

    cout << endl << "r_wb" << endl;
    cout << "No differences" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Ur_wb()(0) << setw(22) << Ur_wb()(1) << setw(22) << Ur_wb()(2) << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rr_wb()(0) << setw(22) << Rr_wb()(1) << setw(22) << Rr_wb()(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << Sr_wb()(0) << setw(22) << Sr_wb()(1) << setw(22) << Sr_wb()(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << Mr_wb()(0) << setw(22) << Mr_wb()(1) << setw(22) << Mr_wb()(2) << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << Zr_wb()(0) << setw(22) << Zr_wb()(1) << setw(22) << Zr_wb()(2) << endl;

    cout << endl << "w_wbw_rps" << endl;
    cout << "Slight difference between space and body" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Uw_wbw_rps()(0) << setw(22) << Uw_wbw_rps()(1) << setw(22) << Uw_wbw_rps()(2) << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rw_wbw_rps()(0) << setw(22) << Rw_wbw_rps()(1) << setw(22) << Rw_wbw_rps()(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << Sw_wbw_rps()(0) << setw(22) << Sw_wbw_rps()(1) << setw(22) << Sw_wbw_rps()(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << Mw_wbw_rps()(0) << setw(22) << Mw_wbw_rps()(1) << setw(22) << Mw_wbw_rps()(2) << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << Zw_wbw_rps()(0) << setw(22) << Zw_wbw_rps()(1) << setw(22) << Zw_wbw_rps()(2) << endl;

    cout << endl << "w_wbb_rps" << endl;
    cout << "Slight difference between space and body" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Uw_wbb_rps()(0) << setw(22) << Uw_wbb_rps()(1) << setw(22) << Uw_wbb_rps()(2) << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rw_wbb_rps()(0) << setw(22) << Rw_wbb_rps()(1) << setw(22) << Rw_wbb_rps()(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << Sw_wbb_rps()(0) << setw(22) << Sw_wbb_rps()(1) << setw(22) << Sw_wbb_rps()(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << Mw_wbb_rps()(0) << setw(22) << Mw_wbb_rps()(1) << setw(22) << Mw_wbb_rps()(2) << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << Zw_wbb_rps()(0) << setw(22) << Zw_wbb_rps()(1) << setw(22) << Zw_wbb_rps()(2) << endl;

    cout << endl << "vi_wbw_mps" << endl;
    cout << "Slight difference between space and body" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Uvi_wbw_mps(0) << setw(22) << Uvi_wbw_mps(1) << setw(22) << Uvi_wbw_mps(2) << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rvi_wbw_mps(0) << setw(22) << Rvi_wbw_mps(1) << setw(22) << Rvi_wbw_mps(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << Svi_wbw_mps(0) << setw(22) << Svi_wbw_mps(1) << setw(22) << Svi_wbw_mps(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << Mvi_wbw_mps(0) << setw(22) << Mvi_wbw_mps(1) << setw(22) << Mvi_wbw_mps(2) << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << Zvi_wbw_mps(0) << setw(22) << Zvi_wbw_mps(1) << setw(22) << Zvi_wbw_mps(2) << "  WRONG" << endl;

    cout << endl << "v_b_mps = vi_wbb_mps" << endl;
    cout << "Slight difference between space and body" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Uv_b_mps(0)  << setw(22) << Uv_b_mps(1)  << setw(22) << Uv_b_mps(2)  << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rv_b_mps(0)  << setw(22) << Rv_b_mps(1)  << setw(22) << Rv_b_mps(2)  << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << Sv_b_mps(0)  << setw(22) << Sv_b_mps(1)  << setw(22) << Sv_b_mps(2)  << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << MAv_b_mps(0) << setw(22) << MAv_b_mps(1) << setw(22) << MAv_b_mps(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << MBv_b_mps(0) << setw(22) << MBv_b_mps(1) << setw(22) << MBv_b_mps(2) << "  DIFF" << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << ZAv_b_mps(0) << setw(22) << ZAv_b_mps(1) << setw(22) << ZAv_b_mps(2) << "  WRONG" << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << ZBv_b_mps(0) << setw(22) << ZBv_b_mps(1) << setw(22) << ZBv_b_mps(2) << "  WRONG" << endl;

    cout << endl << "v_w_mps" << endl;
    cout << "Slight difference between space and body" << endl;
    cout << "Space SO3  " << scientific << setw(22) << setprecision(14) << showpos << Uv_w_mps(0)  << setw(22) << Uv_w_mps(1)  << setw(22) << Uv_w_mps(2)  << endl;
    cout << "Body  SO3  " << scientific << setw(22) << setprecision(14) << showpos << Rv_w_mps(0)  << setw(22) << Rv_w_mps(1)  << setw(22) << Rv_w_mps(2)  << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << SAv_w_mps(0) << setw(22) << SAv_w_mps(1) << setw(22) << SAv_w_mps(2) << endl;
    cout << "Body  SE3  " << scientific << setw(22) << setprecision(14) << showpos << SBv_w_mps(0) << setw(22) << SBv_w_mps(1) << setw(22) << SBv_w_mps(2) << endl;
    cout << "Space SE3  " << scientific << setw(22) << setprecision(14) << showpos << Mv_w_mps(0)  << setw(22) << Mv_w_mps(1) << setw(22)  << Mv_w_mps(2)  << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << MAv_w_mps(0) << setw(22) << MAv_w_mps(1) << setw(22) << MAv_w_mps(2) << endl;
    cout << "Space SE3 G" << scientific << setw(22) << setprecision(14) << showpos << MBv_w_mps(0) << setw(22) << MBv_w_mps(1) << setw(22) << MBv_w_mps(2) << "  DIFF" << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << ZAv_w_mps(0) << setw(22) << ZAv_w_mps(1) << setw(22) << ZAv_w_mps(2) << "  WRONG" << endl;
    cout << "Space SE3 B" << scientific << setw(22) << setprecision(14) << showpos << ZBv_w_mps(0) << setw(22) << ZBv_w_mps(1) << setw(22) << ZBv_w_mps(2) << "  WRONG" << endl;
    cout << endl << endl << endl;

} // closes test_differences

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegration::test_analysis(const ang::speu_dcm& GR0_wb,
                                          const Eigen::Vector3d& v0_w_mps,
                                          const ang::so3_tangent& w0_wbw_rps,
                                          const Eigen::Vector3d& a_b_mps2,
                                          const Eigen::Vector3d& alpha_b_rps2) {
    // With a single integration step using Euler's method, this test analyzes the differences between the two ways of
    // executing the space SE3 integration. One of them coincides with the other integration methods, while the other
    // (the obvious way of doing things) does not. The reason is a second order term (Diff below) which is missing.
    // Nothing intrinsically wrong, it just comes from using Euler integration
    // It would be the same for Heun's, as error accumulates.

    double Deltat_sec = 0.001;

    ang::speu_dcm GRX_wb(GR0_wb.get_dcm(), {1.2, -1.3, 1.7}); // for this test only, so it is not all zeros

    ang::dcm         R0_wb        = GRX_wb.get_dcm();
    Eigen::Vector3d  T0_wbw_m     = GRX_wb.get_T();
    Eigen::Vector3d  v0_b_mps     = GRX_wb & v0_w_mps;
    ang::so3_tangent w0_wbb_rps   = R0_wb % w0_wbw_rps;
    Eigen::Vector3d  vi0_wbb_mps  = v0_b_mps;
    Eigen::Vector3d  T0_wbb_m     = GRX_wb / T0_wbw_m;

    ang::se3_tangent xi0_wbb_mrps(vi0_wbb_mps, w0_wbb_rps);
    ang::se3_tangent xi0_wbw_mrps = GRX_wb | xi0_wbb_mrps;
    Eigen::Vector3d  vi0_wbw_mps  = xi0_wbw_mrps.get_vi();

    Eigen::Vector3d  dv_w_mps_dt     = R0_wb * a_b_mps2;
    Eigen::Vector3d  dT_wbw_m_dt     = v0_w_mps;
    Eigen::Vector3d  dw_wbw_rps_dt   = R0_wb * alpha_b_rps2;
    Eigen::Vector3d  dvi_wbw_mps_dt  = dv_w_mps_dt + T0_wbw_m.cross(dw_wbw_rps_dt) + dT_wbw_m_dt.cross(w0_wbw_rps());

    Eigen::Vector3d  v1_w_mps        = v0_w_mps     + dv_w_mps_dt * Deltat_sec;
    Eigen::Vector3d  T1_wbw_m        = T0_wbw_m     + dT_wbw_m_dt * Deltat_sec;
    Eigen::Vector3d  w1_wbw_rps      = w0_wbw_rps() + dw_wbw_rps_dt * Deltat_sec;
    Eigen::Vector3d  vi1_wbw_mps     = v1_w_mps + T1_wbw_m.cross(w1_wbw_rps);

    // three components
    Eigen::Vector3d  Delta_v         = dv_w_mps_dt * Deltat_sec;
    Eigen::Vector3d  Delta_T         = dT_wbw_m_dt.cross(w0_wbw_rps()) * Deltat_sec;
    Eigen::Vector3d  Delta_w         = T0_wbw_m.cross(dw_wbw_rps_dt) * Deltat_sec;
    Eigen::Vector3d  Delta_2         = dT_wbw_m_dt.cross(dw_wbw_rps_dt) * Deltat_sec * Deltat_sec;
    Eigen::Vector3d  vi2_wbw_mps     = vi0_wbw_mps  + dvi_wbw_mps_dt * Deltat_sec;
    Eigen::Vector3d  Diff            = vi2_wbw_mps - vi1_wbw_mps;

    cout << "vi0       " << scientific << setw(22) << setprecision(14) << showpos << vi0_wbw_mps(0) << setw(22) << vi0_wbw_mps(1) << setw(22) << vi0_wbw_mps(2) << endl;
    cout << "Delta_v   " << scientific << setw(22) << setprecision(14) << showpos << Delta_v(0) << setw(22) << Delta_v(1) << setw(22) << Delta_v(2) << endl;
    cout << "Delta_T   " << scientific << setw(22) << setprecision(14) << showpos << Delta_T(0) << setw(22) << Delta_T(1) << setw(22) << Delta_T(2) << endl;
    cout << "Delta_w   " << scientific << setw(22) << setprecision(14) << showpos << Delta_w(0) << setw(22) << Delta_w(1) << setw(22) << Delta_w(2) << endl;
    cout << "Delta_2   " << scientific << setw(22) << setprecision(14) << showpos << Delta_2(0) << setw(22) << Delta_2(1) << setw(22) << Delta_2(2) << endl;
    cout << "vi1       " << scientific << setw(22) << setprecision(14) << showpos << vi1_wbw_mps(0) << setw(22) << vi1_wbw_mps(1) << setw(22) << vi1_wbw_mps(2) << endl;
    cout << "vi2       " << scientific << setw(22) << setprecision(14) << showpos << vi2_wbw_mps(0) << setw(22) << vi2_wbw_mps(1) << setw(22) << vi2_wbw_mps(2) << endl;
    cout << "Diff      " << scientific << setw(22) << setprecision(14) << showpos << Diff(0) << setw(22) << Diff(1) << setw(22) << Diff(2) << endl;
    cout << endl << endl << endl;

} // closes test_analysis

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////















