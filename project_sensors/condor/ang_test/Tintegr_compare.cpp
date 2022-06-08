#include "Tintegr_compare.h"

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

ang::test::Tintegr_compare::Tintegr_compare(jail::counter& Ocounter)
        : ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tintegr_compare::run() {
    ::jail::unit_test::run();

    test_so3_local();	                   cout << endl << endl;
    test_so3_global();	                   cout << endl << endl;
    test_so3_local_vs_global();            cout << endl << endl;
    test_se3_local();		               cout << endl << endl;
    test_se3_global();		               cout << endl << endl;
    test_se3_local_vs_global();            cout << endl << endl;

    finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_so3_local() {
    double d2r = math::constant::D2R();
    int nel = 100001;
    double Deltat_sec = 0.01;

    so3_tangent w_nbb_rps(1e-3, 2e-3, 3e-3);
    rotv Deltar_nbb(w_nbb_rps() * Deltat_sec);

    vector<dcm> VR_nb1(nel), VR_nb2(nel);
    vector<rodrigues> Vq_nb1(nel), Vq_nb2(nel);
    vector<rotv> Vr_nb1(nel), Vr_nb2(nel);

    euler euler_nb_init(35 * d2r, 7 * d2r, 1 * d2r);
    VR_nb1[0] = euler_nb_init;
    VR_nb2[0] = euler_nb_init;
    Vq_nb1[0] = euler_nb_init;
    Vq_nb2[0] = euler_nb_init;
    Vr_nb1[0] = euler_nb_init;
    Vr_nb2[0] = euler_nb_init;

    for (int i = 1; i != nel; ++i) {

        Eigen::Matrix3d R_nb_dot  = VR_nb1[i-1].omegabody2dot(w_nbb_rps);
        Eigen::Vector4d q_nb_dot  = Vq_nb1[i-1].omegabody2dot(w_nbb_rps).get();
        Eigen::Vector3d r_nb_dot  = Vr_nb1[i-1].omegabody2dot(w_nbb_rps);

        Eigen::Vector9d R_nb_new  = dcm::wedge(VR_nb1[i-1])       + dcm::wedge(R_nb_dot * Deltat_sec);
        Eigen::Vector4d q_nb_new  = rodrigues::wedge(Vq_nb1[i-1]) + q_nb_dot * Deltat_sec;
        Eigen::Vector3d r_nb_new  = rotv::wedge(Vr_nb1[i-1])      + r_nb_dot * Deltat_sec;

        VR_nb1[i]                 = dcm::hat(R_nb_new);
        Vq_nb1[i]                 = rodrigues::hat(q_nb_new);
        Vr_nb1[i]                 = rotv::hat(r_nb_new);

        VR_nb2[i]                  = VR_nb2[i-1].plus_right(Deltar_nbb);
        Vq_nb2[i]                  = Vq_nb2[i-1].plus_right(Deltar_nbb);
        Vr_nb2[i]                  = Vr_nb2[i-1].plus_right(Deltar_nbb);
    }

    euler euler1_R(VR_nb1.back());
    euler euler1_q(Vq_nb1.back());
    euler euler1_r(Vr_nb1.back());
    euler euler2_R(VR_nb2.back());
    euler euler2_q(Vq_nb2.back());
    euler euler2_r(Vr_nb2.back());

    std::cout << "R1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_bank_rad() / d2r
              << std::endl;
    std::cout << "R2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_bank_rad() / d2r
              << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////
    // ============================ RESULTS ============================
    // Integrating in so(3) the results are identical for all three representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for rotation matrix and unit quaternion
    // Bigger difference with integrating in euclidean with rotation vector
    ///////////////////////////////////////////////////////////////////

    // Differences between both methods are way higher for rotation vector than rotation matrix or unit quaternion
    check("R1 - R2 1    ", euler1_R.get_yaw_rad(),   euler2_R.get_yaw_rad(),   1e-8);
    check("R1 - R2 2    ", euler1_R.get_pitch_rad(), euler2_R.get_pitch_rad(), 1e-8);
    check("R1 - R2 3    ", euler1_R.get_bank_rad(),  euler2_R.get_bank_rad(),  1e-8);

    check("q1 - q2 1    ", euler1_q.get_yaw_rad(),   euler2_q.get_yaw_rad(),   1e-9);
    check("q1 - q2 2    ", euler1_q.get_pitch_rad(), euler2_q.get_pitch_rad(), 1e-9);
    check("q1 - q2 3    ", euler1_q.get_bank_rad(),  euler2_q.get_bank_rad(),  1e-9);

    check("r1 - r2 1    ", euler1_r.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-5);
    check("r1 - r2 2    ", euler1_r.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-5);
    check("r1 - r2 3    ", euler1_r.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-5);

    // Euclidean integration
    check("R1 - q1 1    ", euler1_R.get_yaw_rad(),   euler1_q.get_yaw_rad(),   1e-8);
    check("R1 - q1 2    ", euler1_R.get_pitch_rad(), euler1_q.get_pitch_rad(), 1e-8);
    check("R1 - q1 3    ", euler1_R.get_bank_rad(),  euler1_q.get_bank_rad(),  1e-8);

    check("R1 - r1 1    ", euler1_R.get_yaw_rad(),   euler1_r.get_yaw_rad(),   1e-5);
    check("R1 - r1 2    ", euler1_R.get_pitch_rad(), euler1_r.get_pitch_rad(), 1e-5);
    check("R1 - r1 3    ", euler1_R.get_bank_rad(),  euler1_r.get_bank_rad(),  1e-5);

    check("q1 - r1 1    ", euler1_q.get_yaw_rad(),   euler1_r.get_yaw_rad(),   1e-5);
    check("q1 - r1 2    ", euler1_q.get_pitch_rad(), euler1_r.get_pitch_rad(), 1e-5);
    check("q1 - r1 3    ", euler1_q.get_bank_rad(),  euler1_r.get_bank_rad(),  1e-5);

    // so(3) integration
    check("R2 - q2 1    ", euler2_R.get_yaw_rad(),   euler2_q.get_yaw_rad(),   1e-12);
    check("R2 - q2 2    ", euler2_R.get_pitch_rad(), euler2_q.get_pitch_rad(), 1e-12);
    check("R2 - q2 3    ", euler2_R.get_bank_rad(),  euler2_q.get_bank_rad(),  1e-12);

    check("R2 - r2 1    ", euler2_R.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-12);
    check("R2 - r2 2    ", euler2_R.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-12);
    check("R2 - r2 3    ", euler2_R.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-12);

    check("q2 - r2 1    ", euler2_q.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-12);
    check("q2 - r2 2    ", euler2_q.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-12);
    check("q2 - r2 3    ", euler2_q.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-12);

} // closes test_so3_local

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_so3_global() {
    double d2r = math::constant::D2R();
    int nel = 100001;
    double Deltat_sec = 0.01;

    so3_tangent w_nbn_rps(1e-3, 2e-3, 3e-3);
    rotv Deltar_nbn(w_nbn_rps() * Deltat_sec);

    vector<dcm> VR_nb1(nel), VR_nb2(nel);
    vector<rodrigues> Vq_nb1(nel), Vq_nb2(nel);
    vector<rotv> Vr_nb1(nel), Vr_nb2(nel);

    euler euler_nb_init(35 * d2r, 7 * d2r, 1 * d2r);
    VR_nb1[0] = euler_nb_init;
    VR_nb2[0] = euler_nb_init;
    Vq_nb1[0] = euler_nb_init;
    Vq_nb2[0] = euler_nb_init;
    Vr_nb1[0] = euler_nb_init;
    Vr_nb2[0] = euler_nb_init;

    for (int i = 1; i != nel; ++i) {

        Eigen::Matrix3d R_nb_dot  = VR_nb1[i-1].omegaspace2dot(w_nbn_rps);
        Eigen::Vector4d q_nb_dot  = Vq_nb1[i-1].omegaspace2dot(w_nbn_rps).get();
        Eigen::Vector3d r_nb_dot  = Vr_nb1[i-1].omegaspace2dot(w_nbn_rps);

        Eigen::Vector9d R_nb_new  = dcm::wedge(VR_nb1[i-1])       + dcm::wedge(R_nb_dot * Deltat_sec);
        Eigen::Vector4d q_nb_new  = rodrigues::wedge(Vq_nb1[i-1]) + q_nb_dot * Deltat_sec;
        Eigen::Vector3d r_nb_new  = rotv::wedge(Vr_nb1[i-1])      + r_nb_dot * Deltat_sec;

        VR_nb1[i]                 = dcm::hat(R_nb_new);
        Vq_nb1[i]                 = rodrigues::hat(q_nb_new);
        Vr_nb1[i]                 = rotv::hat(r_nb_new);

        VR_nb2[i]                  = VR_nb2[i-1].plus_left(Deltar_nbn);
        Vq_nb2[i]                  = Vq_nb2[i-1].plus_left(Deltar_nbn);
        Vr_nb2[i]                  = Vr_nb2[i-1].plus_left(Deltar_nbn);
    }

    euler euler1_R(VR_nb1.back());
    euler euler1_q(Vq_nb1.back());
    euler euler1_r(Vr_nb1.back());
    euler euler2_R(VR_nb2.back());
    euler euler2_q(Vq_nb2.back());
    euler euler2_r(Vr_nb2.back());

    std::cout << "R1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_R.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_q.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_r.get_bank_rad() / d2r
              << std::endl;
    std::cout << "R2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_R.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_q.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_r.get_bank_rad() / d2r
              << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////
    // ============================ RESULTS ============================
    // Integrating in so(3) the results are identical for all three representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for rotation matrix and unit quaternion
    // Bigger difference with integrating in euclidean with rotation vector
    ///////////////////////////////////////////////////////////////////

    // Differences between both methods are way higher for rotation vector than rotation matrix or unit quaternion
    check("R1 - R2 1    ", euler1_R.get_yaw_rad(),   euler2_R.get_yaw_rad(),   1e-8);
    check("R1 - R2 2    ", euler1_R.get_pitch_rad(), euler2_R.get_pitch_rad(), 1e-8);
    check("R1 - R2 3    ", euler1_R.get_bank_rad(),  euler2_R.get_bank_rad(),  1e-8);

    check("q1 - q2 1    ", euler1_q.get_yaw_rad(),   euler2_q.get_yaw_rad(),   1e-9);
    check("q1 - q2 2    ", euler1_q.get_pitch_rad(), euler2_q.get_pitch_rad(), 1e-9);
    check("q1 - q2 3    ", euler1_q.get_bank_rad(),  euler2_q.get_bank_rad(),  1e-9);

    check("r1 - r2 1    ", euler1_r.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-5);
    check("r1 - r2 2    ", euler1_r.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-5);
    check("r1 - r2 3    ", euler1_r.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-5);

    // Euclidean integration
    check("R1 - q1 1    ", euler1_R.get_yaw_rad(),   euler1_q.get_yaw_rad(),   1e-8);
    check("R1 - q1 2    ", euler1_R.get_pitch_rad(), euler1_q.get_pitch_rad(), 1e-8);
    check("R1 - q1 3    ", euler1_R.get_bank_rad(),  euler1_q.get_bank_rad(),  1e-8);

    check("R1 - r1 1    ", euler1_R.get_yaw_rad(),   euler1_r.get_yaw_rad(),   1e-5);
    check("R1 - r1 2    ", euler1_R.get_pitch_rad(), euler1_r.get_pitch_rad(), 1e-5);
    check("R1 - r1 3    ", euler1_R.get_bank_rad(),  euler1_r.get_bank_rad(),  1e-5);

    check("q1 - r1 1    ", euler1_q.get_yaw_rad(),   euler1_r.get_yaw_rad(),   1e-5);
    check("q1 - r1 2    ", euler1_q.get_pitch_rad(), euler1_r.get_pitch_rad(), 1e-5);
    check("q1 - r1 3    ", euler1_q.get_bank_rad(),  euler1_r.get_bank_rad(),  1e-5);

    // so(3) integration
    check("R2 - q2 1    ", euler2_R.get_yaw_rad(),   euler2_q.get_yaw_rad(),   1e-11);
    check("R2 - q2 2    ", euler2_R.get_pitch_rad(), euler2_q.get_pitch_rad(), 1e-11);
    check("R2 - q2 3    ", euler2_R.get_bank_rad(),  euler2_q.get_bank_rad(),  1e-11);

    check("R2 - r2 1    ", euler2_R.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-11);
    check("R2 - r2 2    ", euler2_R.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-11);
    check("R2 - r2 3    ", euler2_R.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-11);

    check("q2 - r2 1    ", euler2_q.get_yaw_rad(),   euler2_r.get_yaw_rad(),   1e-11);
    check("q2 - r2 2    ", euler2_q.get_pitch_rad(), euler2_r.get_pitch_rad(), 1e-11);
    check("q2 - r2 3    ", euler2_q.get_bank_rad(),  euler2_r.get_bank_rad(),  1e-11);

} // closes test_so3_global

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_so3_local_vs_global() {
    double d2r = math::constant::D2R();
    int nel = 100001;
    double Deltat_sec = 0.01;

    so3_tangent w_nbb_rps(1e-3, 2e-3, 3e-3);
    rotv Deltar_nbb(w_nbb_rps() * Deltat_sec);
    so3_tangent w_nbn_rps_R, w_nbn_rps_q, w_nbn_rps_r;
    rotv Deltar_nbn_R, Deltar_nbn_q, Deltar_nbn_r;

    vector<dcm> VR_nb_l(nel), VR_nb_g(nel);
    vector<rodrigues> Vq_nb_l(nel), Vq_nb_g(nel);
    vector<rotv> Vr_nb_l(nel), Vr_nb_g(nel);

    euler euler_nb_init(35 * d2r, 7 * d2r, 1 * d2r);
    VR_nb_l[0] = euler_nb_init;
    VR_nb_g[0] = euler_nb_init;
    Vq_nb_l[0] = euler_nb_init;
    Vq_nb_g[0] = euler_nb_init;
    Vr_nb_l[0] = euler_nb_init;
    Vr_nb_g[0] = euler_nb_init;

    for (int i = 1; i != nel; ++i) {
        VR_nb_l[i]    = VR_nb_l[i-1].plus_right(Deltar_nbb);
        Vq_nb_l[i]    = Vq_nb_l[i-1].plus_right(Deltar_nbb);
        Vr_nb_l[i]    = Vr_nb_l[i-1].plus_right(Deltar_nbb);

        w_nbn_rps_R() = VR_nb_l[i-1].adjoint_matrix_forward() * w_nbb_rps();
        w_nbn_rps_q() = Vq_nb_l[i-1].adjoint_matrix_forward() * w_nbb_rps();
        w_nbn_rps_r() = Vr_nb_l[i-1].adjoint_matrix_forward() * w_nbb_rps();

        Deltar_nbn_R() = w_nbn_rps_R() * Deltat_sec;
        Deltar_nbn_q() = w_nbn_rps_q() * Deltat_sec;
        Deltar_nbn_r() = w_nbn_rps_r() * Deltat_sec;

        VR_nb_g[i]    = VR_nb_g[i-1].plus_left(Deltar_nbn_R);
        Vq_nb_g[i]    = Vq_nb_g[i-1].plus_left(Deltar_nbn_q);
        Vr_nb_g[i]    = Vr_nb_g[i-1].plus_left(Deltar_nbn_r);
    }

    euler euler_R_l(VR_nb_l.back());
    euler euler_q_l(Vq_nb_l.back());
    euler euler_r_l(Vr_nb_l.back());
    euler euler_R_g(VR_nb_g.back());
    euler euler_q_g(Vq_nb_g.back());
    euler euler_r_g(Vr_nb_g.back());

    std::cout << "R_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "R_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_R_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "q_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_q_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "r_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_r_g.get_bank_rad() / d2r
              << std::endl << std::endl;

    check("R local - q local 1    ", euler_R_l.get_yaw_rad(),   euler_q_l.get_yaw_rad(),   1e-12);
    check("R local - q local 2    ", euler_R_l.get_pitch_rad(), euler_q_l.get_pitch_rad(), 1e-12);
    check("R local - q local 3    ", euler_R_l.get_bank_rad(),  euler_q_l.get_bank_rad(),  1e-12);

    check("R local - r local 1    ", euler_R_l.get_yaw_rad(),   euler_r_l.get_yaw_rad(),   1e-12);
    check("R local - r local 2    ", euler_R_l.get_pitch_rad(), euler_r_l.get_pitch_rad(), 1e-12);
    check("R local - r local 3    ", euler_R_l.get_bank_rad(),  euler_r_l.get_bank_rad(),  1e-12);

    check("R local - R global 1   ", euler_R_l.get_yaw_rad(),   euler_R_g.get_yaw_rad(),   1e-12);
    check("R local - R global 2   ", euler_R_l.get_pitch_rad(), euler_R_g.get_pitch_rad(), 1e-12);
    check("R local - R global 3   ", euler_R_l.get_bank_rad(),  euler_R_g.get_bank_rad(),  1e-12);

    check("R local - q global 1   ", euler_R_l.get_yaw_rad(),   euler_q_g.get_yaw_rad(),   1e-12);
    check("R local - q global 2   ", euler_R_l.get_pitch_rad(), euler_q_g.get_pitch_rad(), 1e-12);
    check("R local - q global 3   ", euler_R_l.get_bank_rad(),  euler_q_g.get_bank_rad(),  1e-12);

    check("R local - r global 1   ", euler_R_l.get_yaw_rad(),   euler_r_g.get_yaw_rad(),   1e-12);
    check("R local - r global 2   ", euler_R_l.get_pitch_rad(), euler_r_g.get_pitch_rad(), 1e-12);
    check("R local - r global 3   ", euler_R_l.get_bank_rad(),  euler_r_g.get_bank_rad(),  1e-12);

} // closes test_so3_local_vs_global

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_se3_local() {
    double d2r = math::constant::D2R();
    int nel = 100001;
    double Deltat_sec = 0.01;

    so3_tangent w_ebb_rps(1e-3, 2e-3, 3e-3);
    Eigen::Vector3d v_b_mps(25.0, 1.2, -0.7);
    se3_tangent xi_ebb_mrps(v_b_mps, w_ebb_rps);
    trfv Deltatau_ebb(xi_ebb_mrps() * Deltat_sec);

    vector<speu_dcm> VgR_eb1(nel), VgR_eb2(nel);
    vector<speu_rodrigues> Vgq_eb1(nel), Vgq_eb2(nel);
    vector<homogeneous> VM_eb1(nel), VM_eb2(nel);
    vector<dual> Vz_eb1(nel), Vz_eb2(nel);
    vector<trfv> Vtau_eb2(nel); // euclidean integration not available

    euler euler_eb_init(35 * d2r, 7 * d2r, 1 * d2r);
    dcm R_eb_init(euler_eb_init);
    rodrigues q_eb_init(euler_eb_init);
    rotv r_eb_init(euler_eb_init);
    Eigen::Vector3d T_ebe_m_init = {300.0, -650.0, 930.0};
    VgR_eb1[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    VgR_eb2[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    Vgq_eb1[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    Vgq_eb2[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    VM_eb1[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    VM_eb2[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    Vz_eb1[0]   = dual(q_eb_init, T_ebe_m_init);
    Vz_eb2[0]   = dual(q_eb_init, T_ebe_m_init);
    Vtau_eb2[0] = trfv(r_eb_init, T_ebe_m_init);

    for (int i = 1; i != nel; ++i) {
        Eigen::Matrix34d gR_eb_dot  = VgR_eb1[i-1].xibody2dot(xi_ebb_mrps);
        Eigen::Vector7d  gq_eb_dot  = Vgq_eb1[i-1].xibody2dot(xi_ebb_mrps);
        Eigen::Matrix4d  M_eb_dot   = VM_eb1[i-1].xibody2dot(xi_ebb_mrps);
        Eigen::Vector8d  z_eb_dot   = Vz_eb1[i-1].xibody2dot(xi_ebb_mrps).get();

        Eigen::Vector12d gR_eb_new  = speu_dcm::wedge(VgR_eb1[i-1])       + speu_dcm::wedge(gR_eb_dot * Deltat_sec);
        Eigen::Vector7d  gq_eb_new  = speu_rodrigues::wedge(Vgq_eb1[i-1]) + gq_eb_dot * Deltat_sec;
        Eigen::Vector12d M_eb_new   = homogeneous::wedge(VM_eb1[i-1])     + homogeneous::wedge(M_eb_dot * Deltat_sec);
        Eigen::Vector8d  z_eb_new   = dual::wedge(Vz_eb1[i-1])            + z_eb_dot * Deltat_sec;

        VgR_eb1[i]                  = speu_dcm::hat(gR_eb_new);
        Vgq_eb1[i]                  = speu_rodrigues::hat(gq_eb_new);
        VM_eb1[i]                   = homogeneous::hat(M_eb_new);
        Vz_eb1[i]                   = dual::hat(z_eb_new);

        VgR_eb2[i]                  = VgR_eb2[i-1].plus_right(Deltatau_ebb);
        Vgq_eb2[i]                  = Vgq_eb2[i-1].plus_right(Deltatau_ebb);
        VM_eb2[i]                   = VM_eb2[i-1].plus_right(Deltatau_ebb);
        Vz_eb2[i]                   = Vz_eb2[i-1].plus_right(Deltatau_ebb);
        Vtau_eb2[i]                 = Vtau_eb2[i-1].plus_right(Deltatau_ebb);
    }

    euler euler1_gR(VgR_eb1.back().get_dcm());
    euler euler1_gq(Vgq_eb1.back().get_rodrigues());
    euler euler1_M(VM_eb1.back().get_dcm());
    euler euler1_z(Vz_eb1.back().get_rodrigues());
    euler euler2_gR(VgR_eb2.back().get_dcm());
    euler euler2_gq(Vgq_eb2.back().get_rodrigues());
    euler euler2_M(VM_eb2.back().get_dcm());
    euler euler2_z(Vz_eb2.back().get_rodrigues());
    euler euler2_tau(Vtau_eb2.back().get_rotv());

    std::cout << "qR1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qq1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qR2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qq2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_bank_rad() / d2r
              << std::endl;
    std::cout << "tau2 "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_bank_rad() / d2r
              << std::endl << std::endl;

    Eigen::Vector3d T1_ebe_m_gR  = VgR_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_gq  = Vgq_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_M   = VM_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_z   = Vz_eb1.back().get_T();
    Eigen::Vector3d T2_ebe_m_gR  = VgR_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_gq  = Vgq_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_M   = VM_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_z   = Vz_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_tau = Vtau_eb2.back().get_T();

    std::cout << "qR1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(2)
              << std::endl;
    std::cout << "qq1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(2)
              << std::endl;
    std::cout << "M1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(2)
              << std:: endl;
    std::cout << "z1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(2)
              << std::endl;
    std::cout << "qR2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(2)
              << std::endl;
    std::cout << "qq2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(2)
              << std::endl;
    std::cout << "M2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(2)
              << std::endl;
    std::cout << "z2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(2)
              << std::endl;
    std::cout << "tau2 "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(2)
              << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////
    // ======================== RESULTS ROTATION =======================
    // Integrating in se(3) the results are identical for all five representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for all representations but
    // - euclidean dcm equal to homogeneous (and to so3 dcm)
    // - euclidean rodrigues equal to dual (and to so3 rodrigues)
    //
    // ======================== RESULTS TRASLATION =======================
    // Integrating in se(3) the results are identical for all five representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for all representations but
    // - euclidean dcm equal to homogeneous
    ///////////////////////////////////////////////////////////////////

    // euclidean integration
    check("gR1 - gR2 1    ", euler1_gR.get_yaw_rad(),   euler2_gR.get_yaw_rad(),   1e-8);
    check("gR1 - gR2 2    ", euler1_gR.get_pitch_rad(), euler2_gR.get_pitch_rad(), 1e-8);
    check("gR1 - gR2 3    ", euler1_gR.get_bank_rad(),  euler2_gR.get_bank_rad(),  1e-8);

    check("gq1 - gq2 1    ", euler1_gq.get_yaw_rad(),   euler2_gq.get_yaw_rad(),   1e-9);
    check("gq1 - gq2 2    ", euler1_gq.get_pitch_rad(), euler2_gq.get_pitch_rad(), 1e-9);
    check("gq1 - gq2 3    ", euler1_gq.get_bank_rad(),  euler2_gq.get_bank_rad(),  1e-9);

    check("M1 - M2 1      ", euler1_M.get_yaw_rad(),    euler2_M.get_yaw_rad(),    1e-8);
    check("M1 - M2 2      ", euler1_M.get_pitch_rad(),  euler2_M.get_pitch_rad(),  1e-8);
    check("M1 - M2 3      ", euler1_M.get_bank_rad(),   euler2_M.get_bank_rad(),   1e-8);

    check("z1 - z2 1      ", euler1_z.get_yaw_rad(),    euler2_z.get_yaw_rad(),    1e-9);
    check("z1 - z2 2      ", euler1_z.get_pitch_rad(),  euler2_z.get_pitch_rad(),  1e-9);
    check("z1 - z2 3      ", euler1_z.get_bank_rad(),   euler2_z.get_bank_rad(),   1e-9);

    // se(3) integration
    check("gR2 - gq2 1    ", euler2_gR.get_yaw_rad(),   euler2_gq.get_yaw_rad(),   1e-12);
    check("gR2 - gq2 2    ", euler2_gR.get_pitch_rad(), euler2_gq.get_pitch_rad(), 1e-12);
    check("gR2 - gq2 3    ", euler2_gR.get_bank_rad(),  euler2_gq.get_bank_rad(),  1e-12);

    check("gR2 - M2 1     ", euler2_gR.get_yaw_rad(),   euler2_M.get_yaw_rad(),   1e-12);
    check("gR2 - M2 2     ", euler2_gR.get_pitch_rad(), euler2_M.get_pitch_rad(), 1e-12);
    check("gR2 - M2 3     ", euler2_gR.get_bank_rad(),  euler2_M.get_bank_rad(),  1e-12);

    check("gR2 - z2 1     ", euler2_gR.get_yaw_rad(),   euler2_z.get_yaw_rad(),   1e-12);
    check("gR2 - z2 2     ", euler2_gR.get_pitch_rad(), euler2_z.get_pitch_rad(), 1e-12);
    check("gR2 - z2 3     ", euler2_gR.get_bank_rad(),  euler2_z.get_bank_rad(),  1e-12);

    check("gR2 - tau2 1   ", euler2_gR.get_yaw_rad(),   euler2_tau.get_yaw_rad(),   1e-12);
    check("gR2 - tau2 2   ", euler2_gR.get_pitch_rad(), euler2_tau.get_pitch_rad(), 1e-12);
    check("gR2 - tau2 3   ", euler2_gR.get_bank_rad(),  euler2_tau.get_bank_rad(),  1e-12);

    check("gR1 - gR2 1    ", euler1_gR.get_yaw_rad(),   euler2_gR.get_yaw_rad(),   1e-8);
    check("gR1 - gR2 2    ", euler1_gR.get_pitch_rad(), euler2_gR.get_pitch_rad(), 1e-8);
    check("gR1 - gR2 3    ", euler1_gR.get_bank_rad(),  euler2_gR.get_bank_rad(),  1e-8);

    // euclidean integration
    check("gR1 - gR2 1    ", T1_ebe_m_gR(0), T2_ebe_m_gR(0), 1e-0);
    check("gR1 - gR2 2    ", T1_ebe_m_gR(1), T2_ebe_m_gR(1), 1e-0);
    check("gR1 - gR2 3    ", T1_ebe_m_gR(2), T2_ebe_m_gR(2), 1e-0);

    check("gq1 - gq2 1    ", T1_ebe_m_gq(0), T2_ebe_m_gq(0), 1e-0);
    check("gq1 - gq2 2    ", T1_ebe_m_gq(1), T2_ebe_m_gq(1), 1e-0);
    check("gq1 - gq2 3    ", T1_ebe_m_gq(2), T2_ebe_m_gq(2), 1e-0);

    check("M1 - M2 1      ", T1_ebe_m_M(0), T2_ebe_m_M(0), 1e-0);
    check("M1 - M2 2      ", T1_ebe_m_M(1), T2_ebe_m_M(1), 1e-0);
    check("M1 - M2 3      ", T1_ebe_m_M(2), T2_ebe_m_M(2), 1e-0);

    check("z1 - z2 1      ", T1_ebe_m_z(0), T2_ebe_m_z(0), 1e-0);
    check("z1 - z2 2      ", T1_ebe_m_z(1), T2_ebe_m_z(1), 1e-0);
    check("z1 - z2 3      ", T1_ebe_m_z(2), T2_ebe_m_z(2), 1e-0);

    // se(3) integration
    check("gR2 - gq2 1    ", T2_ebe_m_gR(0), T2_ebe_m_gq(0), 1e-8);
    check("gR2 - gq2 2    ", T2_ebe_m_gR(1), T2_ebe_m_gq(1), 1e-8);
    check("gR2 - gq2 3    ", T2_ebe_m_gR(2), T2_ebe_m_gq(2), 1e-8);

    check("gR2 - M2 1     ", T2_ebe_m_gR(0), T2_ebe_m_M(0), 1e-7);
    check("gR2 - M2 2     ", T2_ebe_m_gR(1), T2_ebe_m_M(1), 1e-7);
    check("gR2 - M2 3     ", T2_ebe_m_gR(2), T2_ebe_m_M(2), 1e-7);

    check("gR2 - z2 1     ", T2_ebe_m_gR(0), T2_ebe_m_z(0), 1e-7);
    check("gR2 - z2 2     ", T2_ebe_m_gR(1), T2_ebe_m_z(1), 1e-7);
    check("gR2 - z2 3     ", T2_ebe_m_gR(2), T2_ebe_m_z(2), 1e-7);

    check("gR2 - tau2 1   ", T2_ebe_m_gR(0), T2_ebe_m_tau(0), 1e-8);
    check("gR2 - tau2 2   ", T2_ebe_m_gR(1), T2_ebe_m_tau(1), 1e-8);
    check("gR2 - tau2 3   ", T2_ebe_m_gR(2), T2_ebe_m_tau(2), 1e-8);

} // closes test_se3_local

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_se3_global() {
    double d2r = math::constant::D2R();
    int nel = 3; ///100001;
    double Deltat_sec = 0.01;

    so3_tangent w_ebe_rps(1e-3, 2e-3, 3e-3);
    Eigen::Vector3d vi_e_mps(25.0, 1.2, -0.7);
    se3_tangent xi_ebe_mrps(vi_e_mps, w_ebe_rps);
    trfv Deltatau_ebe(xi_ebe_mrps() * Deltat_sec);

    vector<speu_dcm> VgR_eb1(nel), VgR_eb2(nel);
    vector<speu_rodrigues> Vgq_eb1(nel), Vgq_eb2(nel);
    vector<homogeneous> VM_eb1(nel), VM_eb2(nel);
    vector<dual> Vz_eb1(nel), Vz_eb2(nel);
    vector<trfv> Vtau_eb2(nel); // euclidean integration not available

    euler euler_eb_init(35 * d2r, 7 * d2r, 1 * d2r);
    dcm R_eb_init(euler_eb_init);
    rodrigues q_eb_init(euler_eb_init);
    rotv r_eb_init(euler_eb_init);
    Eigen::Vector3d T_ebe_m_init = {300.0, -650.0, 930.0};
    VgR_eb1[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    VgR_eb2[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    Vgq_eb1[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    Vgq_eb2[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    VM_eb1[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    VM_eb2[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    Vz_eb1[0]   = dual(q_eb_init, T_ebe_m_init);
    Vz_eb2[0]   = dual(q_eb_init, T_ebe_m_init);
    Vtau_eb2[0] = trfv(r_eb_init, T_ebe_m_init);

    for (int i = 1; i != nel; ++i) {
        Eigen::Matrix34d gR_eb_dot  = VgR_eb1[i-1].xispace2dot(xi_ebe_mrps);
        Eigen::Vector7d  gq_eb_dot  = Vgq_eb1[i-1].xispace2dot(xi_ebe_mrps);
        Eigen::Matrix4d  M_eb_dot   = VM_eb1[i-1].xispace2dot(xi_ebe_mrps);
        Eigen::Vector8d  z_eb_dot   = Vz_eb1[i-1].xispace2dot(xi_ebe_mrps).get();

        Eigen::Vector12d gR_eb_new  = speu_dcm::wedge(VgR_eb1[i-1])       + speu_dcm::wedge(gR_eb_dot * Deltat_sec);
        Eigen::Vector7d  gq_eb_new  = speu_rodrigues::wedge(Vgq_eb1[i-1]) + gq_eb_dot * Deltat_sec;
        Eigen::Vector12d M_eb_new   = homogeneous::wedge(VM_eb1[i-1])     + homogeneous::wedge(M_eb_dot * Deltat_sec);
        Eigen::Vector8d  z_eb_new   = dual::wedge(Vz_eb1[i-1])            + z_eb_dot * Deltat_sec;

        VgR_eb1[i]                  = speu_dcm::hat(gR_eb_new);
        Vgq_eb1[i]                  = speu_rodrigues::hat(gq_eb_new);
        VM_eb1[i]                   = homogeneous::hat(M_eb_new);
        Vz_eb1[i]                   = dual::hat(z_eb_new);

        VgR_eb2[i]                  = VgR_eb2[i-1].plus_left(Deltatau_ebe);
        Vgq_eb2[i]                  = Vgq_eb2[i-1].plus_left(Deltatau_ebe);
        VM_eb2[i]                   = VM_eb2[i-1].plus_left(Deltatau_ebe);
        Vz_eb2[i]                   = Vz_eb2[i-1].plus_left(Deltatau_ebe);
        Vtau_eb2[i]                 = Vtau_eb2[i-1].plus_left(Deltatau_ebe);
    }

    euler euler1_gR(VgR_eb1.back().get_dcm());
    euler euler1_gq(Vgq_eb1.back().get_rodrigues());
    euler euler1_M(VM_eb1.back().get_dcm());
    euler euler1_z(Vz_eb1.back().get_rodrigues());
    euler euler2_gR(VgR_eb2.back().get_dcm());
    euler euler2_gq(Vgq_eb2.back().get_rodrigues());
    euler euler2_M(VM_eb2.back().get_dcm());
    euler euler2_z(Vz_eb2.back().get_rodrigues());
    euler euler2_tau(Vtau_eb2.back().get_rotv());

    std::cout << "qR1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gR.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qq1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_gq.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_M.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler1_z.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qR2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gR.get_bank_rad() / d2r
              << std::endl;
    std::cout << "qq2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_gq.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_M.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_z.get_bank_rad() / d2r
              << std::endl;
    std::cout << "tau2 "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler2_tau.get_bank_rad() / d2r
              << std::endl << std::endl;

    Eigen::Vector3d T1_ebe_m_gR  = VgR_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_gq  = Vgq_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_M   = VM_eb1.back().get_T();
    Eigen::Vector3d T1_ebe_m_z   = Vz_eb1.back().get_T();
    Eigen::Vector3d T2_ebe_m_gR  = VgR_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_gq  = Vgq_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_M   = VM_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_z   = Vz_eb2.back().get_T();
    Eigen::Vector3d T2_ebe_m_tau = Vtau_eb2.back().get_T();

    std::cout << "qR1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gR(2)
              << std::endl;
    std::cout << "qq1  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_gq(2)
              << std::endl;
    std::cout << "M1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_M(2)
              << std:: endl;
    std::cout << "z1   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T1_ebe_m_z(2)
              << std::endl;
    std::cout << "qR2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gR(2)
              << std::endl;
    std::cout << "qq2  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_gq(2)
              << std::endl;
    std::cout << "M2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_M(2)
              << std::endl;
    std::cout << "z2   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_z(2)
              << std::endl;
    std::cout << "tau2 "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T2_ebe_m_tau(2)
              << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////
    // ======================== RESULTS ROTATION =======================
    // Integrating in se(3) the results are identical for all five representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for all representations but
    // - euclidean dcm equal to homogeneous (and to so3 dcm)
    // - euclidean rodrigues equal to dual (and to so3 rodrigues)
    //
    // ======================== RESULTS TRASLATION =======================
    // Integrating in se(3) the results are identical for all five representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Small differences with integrating in euclidean for all representations but
    // - euclidean dcm equal to homogeneous
    ///////////////////////////////////////////////////////////////////

    // euclidean integration
    check("gR1 - gR2 1    ", euler1_gR.get_yaw_rad(),   euler2_gR.get_yaw_rad(),   1e-8);
    check("gR1 - gR2 2    ", euler1_gR.get_pitch_rad(), euler2_gR.get_pitch_rad(), 1e-8);
    check("gR1 - gR2 3    ", euler1_gR.get_bank_rad(),  euler2_gR.get_bank_rad(),  1e-8);

    check("gq1 - gq2 1    ", euler1_gq.get_yaw_rad(),   euler2_gq.get_yaw_rad(),   1e-9);
    check("gq1 - gq2 2    ", euler1_gq.get_pitch_rad(), euler2_gq.get_pitch_rad(), 1e-9);
    check("gq1 - gq2 3    ", euler1_gq.get_bank_rad(),  euler2_gq.get_bank_rad(),  1e-9);

    check("M1 - M2 1      ", euler1_M.get_yaw_rad(),    euler2_M.get_yaw_rad(),    1e-8);
    check("M1 - M2 2      ", euler1_M.get_pitch_rad(),  euler2_M.get_pitch_rad(),  1e-8);
    check("M1 - M2 3      ", euler1_M.get_bank_rad(),   euler2_M.get_bank_rad(),   1e-8);

    check("z1 - z2 1      ", euler1_z.get_yaw_rad(),    euler2_z.get_yaw_rad(),    1e-9);
    check("z1 - z2 2      ", euler1_z.get_pitch_rad(),  euler2_z.get_pitch_rad(),  1e-9);
    check("z1 - z2 3      ", euler1_z.get_bank_rad(),   euler2_z.get_bank_rad(),   1e-9);

    // se(3) integration
    check("gR2 - gq2 1    ", euler2_gR.get_yaw_rad(),   euler2_gq.get_yaw_rad(),   1e-12);
    check("gR2 - gq2 2    ", euler2_gR.get_pitch_rad(), euler2_gq.get_pitch_rad(), 1e-12);
    check("gR2 - gq2 3    ", euler2_gR.get_bank_rad(),  euler2_gq.get_bank_rad(),  1e-12);

    check("gR2 - M2 1     ", euler2_gR.get_yaw_rad(),   euler2_M.get_yaw_rad(),   1e-12);
    check("gR2 - M2 2     ", euler2_gR.get_pitch_rad(), euler2_M.get_pitch_rad(), 1e-12);
    check("gR2 - M2 3     ", euler2_gR.get_bank_rad(),  euler2_M.get_bank_rad(),  1e-12);

    check("gR2 - z2 1     ", euler2_gR.get_yaw_rad(),   euler2_z.get_yaw_rad(),   1e-12);
    check("gR2 - z2 2     ", euler2_gR.get_pitch_rad(), euler2_z.get_pitch_rad(), 1e-12);
    check("gR2 - z2 3     ", euler2_gR.get_bank_rad(),  euler2_z.get_bank_rad(),  1e-12);

    check("gR2 - tau2 1   ", euler2_gR.get_yaw_rad(),   euler2_tau.get_yaw_rad(),   1e-11);
    check("gR2 - tau2 2   ", euler2_gR.get_pitch_rad(), euler2_tau.get_pitch_rad(), 1e-11);
    check("gR2 - tau2 3   ", euler2_gR.get_bank_rad(),  euler2_tau.get_bank_rad(),  1e-11);

    check("gR1 - gR2 1    ", euler1_gR.get_yaw_rad(),   euler2_gR.get_yaw_rad(),   1e-8);
    check("gR1 - gR2 2    ", euler1_gR.get_pitch_rad(), euler2_gR.get_pitch_rad(), 1e-8);
    check("gR1 - gR2 3    ", euler1_gR.get_bank_rad(),  euler2_gR.get_bank_rad(),  1e-8);

    // euclidean integration
    check("gR1 - gR2 1    ", T1_ebe_m_gR(0), T2_ebe_m_gR(0), 1e-0);
    check("gR1 - gR2 2    ", T1_ebe_m_gR(1), T2_ebe_m_gR(1), 1e-0);
    check("gR1 - gR2 3    ", T1_ebe_m_gR(2), T2_ebe_m_gR(2), 1e-0);

    check("gq1 - gq2 1    ", T1_ebe_m_gq(0), T2_ebe_m_gq(0), 1e-0);
    check("gq1 - gq2 2    ", T1_ebe_m_gq(1), T2_ebe_m_gq(1), 1e-0);
    check("gq1 - gq2 3    ", T1_ebe_m_gq(2), T2_ebe_m_gq(2), 1e-0);

    check("M1 - M2 1      ", T1_ebe_m_M(0), T2_ebe_m_M(0), 1e-0);
    check("M1 - M2 2      ", T1_ebe_m_M(1), T2_ebe_m_M(1), 1e-0);
    check("M1 - M2 3      ", T1_ebe_m_M(2), T2_ebe_m_M(2), 1e-0);

    check("z1 - z2 1      ", T1_ebe_m_z(0), T2_ebe_m_z(0), 1e-0);
    check("z1 - z2 2      ", T1_ebe_m_z(1), T2_ebe_m_z(1), 1e-0);
    check("z1 - z2 3      ", T1_ebe_m_z(2), T2_ebe_m_z(2), 1e-0);

    // se(3) integration
    check("gR2 - gq2 1    ", T2_ebe_m_gR(0), T2_ebe_m_gq(0), 1e-7);
    check("gR2 - gq2 2    ", T2_ebe_m_gR(1), T2_ebe_m_gq(1), 1e-7);
    check("gR2 - gq2 3    ", T2_ebe_m_gR(2), T2_ebe_m_gq(2), 1e-7);

    check("gR2 - M2 1     ", T2_ebe_m_gR(0), T2_ebe_m_M(0), 1e-7);
    check("gR2 - M2 2     ", T2_ebe_m_gR(1), T2_ebe_m_M(1), 1e-7);
    check("gR2 - M2 3     ", T2_ebe_m_gR(2), T2_ebe_m_M(2), 1e-7);

    check("gR2 - z2 1     ", T2_ebe_m_gR(0), T2_ebe_m_z(0), 1e-7);
    check("gR2 - z2 2     ", T2_ebe_m_gR(1), T2_ebe_m_z(1), 1e-7);
    check("gR2 - z2 3     ", T2_ebe_m_gR(2), T2_ebe_m_z(2), 1e-7);

    check("gR2 - tau2 1   ", T2_ebe_m_gR(0), T2_ebe_m_tau(0), 1e-7);
    check("gR2 - tau2 2   ", T2_ebe_m_gR(1), T2_ebe_m_tau(1), 1e-7);
    check("gR2 - tau2 3   ", T2_ebe_m_gR(2), T2_ebe_m_tau(2), 1e-7);

} // closes test_se3_global

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tintegr_compare::test_se3_local_vs_global() {
    double d2r = math::constant::D2R();
    int nel = 100001;
    double Deltat_sec = 0.01;

    so3_tangent w_ebb_rps(1e-3, 2e-3, 3e-3);
    Eigen::Vector3d v_b_mps(25.0, 1.2, -0.7);
    se3_tangent xi_ebb_mrps(v_b_mps, w_ebb_rps);
    trfv Deltatau_ebb(xi_ebb_mrps() * Deltat_sec);
    se3_tangent xi_ebe_mrps_gR, xi_ebe_mrps_gq, xi_ebe_mrps_M, xi_ebe_mrps_z, xi_ebe_mrps_tau;
    trfv Deltatau_ebe_gR, Deltatau_ebe_gq, Deltatau_ebe_M, Deltatau_ebe_z, Deltatau_ebe_tau;

    vector<speu_dcm> VgR_eb_l(nel), VgR_eb_g(nel);
    vector<speu_rodrigues> Vgq_eb_l(nel), Vgq_eb_g(nel);
    vector<homogeneous> VM_eb_l(nel), VM_eb_g(nel);
    vector<dual> Vz_eb_l(nel), Vz_eb_g(nel);
    vector<trfv> Vtau_eb_l(nel), Vtau_eb_g(nel);

    euler euler_eb_init(35 * d2r, 7 * d2r, 1 * d2r);
    dcm R_eb_init(euler_eb_init);
    rodrigues q_eb_init(euler_eb_init);
    rotv r_eb_init(euler_eb_init);
    Eigen::Vector3d T_ebe_m_init = {300.0, -650.0, 930.0};
    VgR_eb_l[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    VgR_eb_g[0]  = speu_dcm(R_eb_init, T_ebe_m_init);
    Vgq_eb_l[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    Vgq_eb_g[0]  = speu_rodrigues(q_eb_init, T_ebe_m_init);
    VM_eb_l[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    VM_eb_g[0]   = homogeneous(R_eb_init, T_ebe_m_init);
    Vz_eb_l[0]   = dual(q_eb_init, T_ebe_m_init);
    Vz_eb_g[0]   = dual(q_eb_init, T_ebe_m_init);
    Vtau_eb_l[0] = trfv(r_eb_init, T_ebe_m_init);
    Vtau_eb_g[0] = trfv(r_eb_init, T_ebe_m_init);

    for (int i = 1; i != nel; ++i) {
        VgR_eb_l[i]       = VgR_eb_l[i-1].plus_right(Deltatau_ebb);
        Vgq_eb_l[i]       = Vgq_eb_l[i-1].plus_right(Deltatau_ebb);
        VM_eb_l[i]        = VM_eb_l[i-1].plus_right(Deltatau_ebb);
        Vz_eb_l[i]        = Vz_eb_l[i-1].plus_right(Deltatau_ebb);
        Vtau_eb_l[i]      = Vtau_eb_l[i-1].plus_right(Deltatau_ebb);

        xi_ebe_mrps_gR()  = VgR_eb_l[i-1].adjoint_matrix_forward() * xi_ebb_mrps();
        xi_ebe_mrps_gq()  = Vgq_eb_l[i-1].adjoint_matrix_forward() * xi_ebb_mrps();
        xi_ebe_mrps_M()   = VM_eb_l[i-1].adjoint_matrix_forward() * xi_ebb_mrps();
        xi_ebe_mrps_z()   = Vz_eb_l[i-1].adjoint_matrix_forward() * xi_ebb_mrps();
        xi_ebe_mrps_tau() = Vtau_eb_l[i-1].adjoint_matrix_forward() * xi_ebb_mrps();

        Deltatau_ebe_gR()  = xi_ebe_mrps_gR()  * Deltat_sec;
        Deltatau_ebe_gq()  = xi_ebe_mrps_gq()  * Deltat_sec;
        Deltatau_ebe_M()   = xi_ebe_mrps_M()   * Deltat_sec;
        Deltatau_ebe_z()   = xi_ebe_mrps_z()   * Deltat_sec;
        Deltatau_ebe_tau() = xi_ebe_mrps_tau() * Deltat_sec;

        VgR_eb_g[i]       = VgR_eb_g[i-1].plus_left(Deltatau_ebe_gR);
        Vgq_eb_g[i]       = Vgq_eb_g[i-1].plus_left(Deltatau_ebe_gq);
        VM_eb_g[i]        = VM_eb_g[i-1].plus_left(Deltatau_ebe_M);
        Vz_eb_g[i]        = Vz_eb_g[i-1].plus_left(Deltatau_ebe_z);
        Vtau_eb_g[i]      = Vtau_eb_g[i-1].plus_left(Deltatau_ebe_tau);
    }

    euler euler_gR_l(VgR_eb_l.back().get_dcm());
    euler euler_gq_l(Vgq_eb_l.back().get_rodrigues());
    euler euler_M_l(VM_eb_l.back().get_dcm());
    euler euler_z_l(Vz_eb_l.back().get_rodrigues());
    euler euler_tau_l(Vtau_eb_l.back().get_rotv());

    euler euler_gR_g(VgR_eb_g.back().get_dcm());
    euler euler_gq_g(Vgq_eb_g.back().get_rodrigues());
    euler euler_M_g(VM_eb_g.back().get_dcm());
    euler euler_z_g(Vz_eb_g.back().get_rodrigues());
    euler euler_tau_g(Vtau_eb_g.back().get_rotv());

    std::cout << "gR_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "gq_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M_local   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z_local   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "tau_local "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_l.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_l.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_l.get_bank_rad() / d2r
              << std::endl;
    std::cout << "gR_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gR_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "gq_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_gq_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "M_global  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_M_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "z_global  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_z_g.get_bank_rad() / d2r
              << std::endl;
    std::cout << "tau_global"
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_g.get_yaw_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_g.get_pitch_rad() / d2r
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << euler_tau_g.get_bank_rad() / d2r
              << std::endl << std::endl;

    Eigen::Vector3d T_ebe_m_gR_l  = VgR_eb_l.back().get_T();
    Eigen::Vector3d T_ebe_m_gq_l  = Vgq_eb_l.back().get_T();
    Eigen::Vector3d T_ebe_m_M_l   = VM_eb_l.back().get_T();
    Eigen::Vector3d T_ebe_m_z_l   = Vz_eb_l.back().get_T();
    Eigen::Vector3d T_ebe_m_tau_l = Vtau_eb_l.back().get_T();
    Eigen::Vector3d T_ebe_m_gR_g  = VgR_eb_g.back().get_T();
    Eigen::Vector3d T_ebe_m_gq_g  = Vgq_eb_g.back().get_T();
    Eigen::Vector3d T_ebe_m_M_g   = VM_eb_g.back().get_T();
    Eigen::Vector3d T_ebe_m_z_g   = Vz_eb_g.back().get_T();
    Eigen::Vector3d T_ebe_m_tau_g = Vtau_eb_g.back().get_T();


    std::cout << "qR_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_l(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_l(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_l(2)
              << std::endl;
    std::cout << "qq_local  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_l(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_l(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_l(2)
              << std::endl;
    std::cout << "M_local   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_l(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_l(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_l(2)
              << std:: endl;
    std::cout << "z_local   "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_l(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_l(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_l(2)
              << std::endl;
    std::cout << "tau_local "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_l(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_l(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_l(2)
              << std::endl;
    std::cout << "qR_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_g(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_g(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gR_g(2)
              << std::endl;
    std::cout << "qq_global "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_g(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_g(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_gq_g(2)
              << std::endl;
    std::cout << "M_global  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_g(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_g(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_M_g(2)
              << std::endl;
    std::cout << "z_global  "
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_g(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_g(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_z_g(2)
              << std::endl;
    std::cout << "tau_global"
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_g(0)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_g(1)
              << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos << T_ebe_m_tau_g(2)
              << std::endl << std::endl;


    check("gR local - gq local 1    ", euler_gR_l.get_yaw_rad(),   euler_gq_l.get_yaw_rad(),   1e-12);
    check("gR local - gq local 2    ", euler_gR_l.get_pitch_rad(), euler_gq_l.get_pitch_rad(), 1e-12);
    check("gR local - gq local 3    ", euler_gR_l.get_bank_rad(),  euler_gq_l.get_bank_rad(),  1e-12);

    check("gR local - M local 1     ", euler_gR_l.get_yaw_rad(),   euler_M_l.get_yaw_rad(),   1e-12);
    check("gR local - M local 2     ", euler_gR_l.get_pitch_rad(), euler_M_l.get_pitch_rad(), 1e-12);
    check("gR local - M local 3     ", euler_gR_l.get_bank_rad(),  euler_M_l.get_bank_rad(),  1e-12);

    check("gR local - z local 1     ", euler_gR_l.get_yaw_rad(),   euler_z_l.get_yaw_rad(),   1e-12);
    check("gR local - z local 2     ", euler_gR_l.get_pitch_rad(), euler_z_l.get_pitch_rad(), 1e-12);
    check("gR local - z local 3     ", euler_gR_l.get_bank_rad(),  euler_z_l.get_bank_rad(),  1e-12);

    check("gR local - tau local 1   ", euler_gR_l.get_yaw_rad(),   euler_tau_l.get_yaw_rad(),   1e-12);
    check("gR local - tau local 2   ", euler_gR_l.get_pitch_rad(), euler_tau_l.get_pitch_rad(), 1e-12);
    check("gR local - tau local 3   ", euler_gR_l.get_bank_rad(),  euler_tau_l.get_bank_rad(),  1e-12);

    check("gR local - gR global 1   ", euler_gR_l.get_yaw_rad(),   euler_gR_g.get_yaw_rad(),   1e-12);
    check("gR local - gR global 2   ", euler_gR_l.get_pitch_rad(), euler_gR_g.get_pitch_rad(), 1e-12);
    check("gR local - gR global 3   ", euler_gR_l.get_bank_rad(),  euler_gR_g.get_bank_rad(),  1e-12);

    check("gR local - gq global 1   ", euler_gR_l.get_yaw_rad(),   euler_gq_g.get_yaw_rad(),   1e-12);
    check("gR local - gq global 2   ", euler_gR_l.get_pitch_rad(), euler_gq_g.get_pitch_rad(), 1e-12);
    check("gR local - gq global 3   ", euler_gR_l.get_bank_rad(),  euler_gq_g.get_bank_rad(),  1e-12);

    check("gR local - M global 1    ", euler_gR_l.get_yaw_rad(),   euler_M_g.get_yaw_rad(),   1e-10);
    check("gR local - M global 2    ", euler_gR_l.get_pitch_rad(), euler_M_g.get_pitch_rad(), 1e-10);
    check("gR local - M global 3    ", euler_gR_l.get_bank_rad(),  euler_M_g.get_bank_rad(),  1e-10);

    check("gR local - z global 1    ", euler_gR_l.get_yaw_rad(),   euler_z_g.get_yaw_rad(),   1e-12);
    check("gR local - z global 2    ", euler_gR_l.get_pitch_rad(), euler_z_g.get_pitch_rad(), 1e-12);
    check("gR local - z global 3    ", euler_gR_l.get_bank_rad(),  euler_z_g.get_bank_rad(),  1e-12);

    check("gR local - tau global 1  ", euler_gR_l.get_yaw_rad(),   euler_tau_g.get_yaw_rad(),   1e-12);
    check("gR local - tau global 2  ", euler_gR_l.get_pitch_rad(), euler_tau_g.get_pitch_rad(), 1e-12);
    check("gR local - tau global 3  ", euler_gR_l.get_bank_rad(),  euler_tau_g.get_bank_rad(),  1e-12);


    check("gR local - gq local 1    ", T_ebe_m_gR_l(0), T_ebe_m_gq_l(0),  1e-8);
    check("gR local - gq local 2    ", T_ebe_m_gR_l(1), T_ebe_m_gq_l(1),  1e-8);
    check("gR local - gq local 3    ", T_ebe_m_gR_l(2), T_ebe_m_gq_l(2),  1e-8);

    check("gR local - M local 1     ", T_ebe_m_gR_l(0), T_ebe_m_M_l(0),  1e-7);
    check("gR local - M local 2     ", T_ebe_m_gR_l(1), T_ebe_m_M_l(1),  1e-7);
    check("gR local - M local 3     ", T_ebe_m_gR_l(2), T_ebe_m_M_l(2),  1e-7);

    check("gR local - z local 1     ", T_ebe_m_gR_l(0), T_ebe_m_z_l(0),  1e-7);
    check("gR local - z local 2     ", T_ebe_m_gR_l(1), T_ebe_m_z_l(1),  1e-7);
    check("gR local - z local 3     ", T_ebe_m_gR_l(2), T_ebe_m_z_l(2),  1e-7);

    check("gR local - tau local 1   ", T_ebe_m_gR_l(0), T_ebe_m_tau_l(0),  1e-8);
    check("gR local - tau local 2   ", T_ebe_m_gR_l(1), T_ebe_m_tau_l(1),  1e-8);
    check("gR local - tau local 3   ", T_ebe_m_gR_l(2), T_ebe_m_tau_l(2),  1e-8);

    check("gR local - gR global 1   ", T_ebe_m_gR_l(0), T_ebe_m_gR_g(0),  1e-7);
    check("gR local - gR global 2   ", T_ebe_m_gR_l(1), T_ebe_m_gR_g(1),  1e-7);
    check("gR local - gR global 3   ", T_ebe_m_gR_l(2), T_ebe_m_gR_g(2),  1e-7);

    check("gR local - gq global 1   ", T_ebe_m_gR_l(0), T_ebe_m_gq_g(0),  1e-7);
    check("gR local - gq global 2   ", T_ebe_m_gR_l(1), T_ebe_m_gq_g(1),  1e-7);
    check("gR local - gq global 3   ", T_ebe_m_gR_l(2), T_ebe_m_gq_g(2),  1e-7);

    check("gR local - M global 1    ", T_ebe_m_gR_l(0), T_ebe_m_M_g(0),  1e-7);
    check("gR local - M global 2    ", T_ebe_m_gR_l(1), T_ebe_m_M_g(1),  1e-7);
    check("gR local - M global 3    ", T_ebe_m_gR_l(2), T_ebe_m_M_g(2),  1e-7);

    check("gR local - z global 1    ", T_ebe_m_gR_l(0), T_ebe_m_z_g(0),  1e-7);
    check("gR local - z global 2    ", T_ebe_m_gR_l(1), T_ebe_m_z_g(1),  1e-7);
    check("gR local - z global 3    ", T_ebe_m_gR_l(2), T_ebe_m_z_g(2),  1e-7);

    check("gR local - tau global 1  ", T_ebe_m_gR_l(0), T_ebe_m_tau_g(0), 1e-7);
    check("gR local - tau global 2  ", T_ebe_m_gR_l(1), T_ebe_m_tau_g(1), 1e-7);
    check("gR local - tau global 3  ", T_ebe_m_gR_l(2), T_ebe_m_tau_g(2), 1e-7);


} // closes test_se3_local_vs_global

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////











