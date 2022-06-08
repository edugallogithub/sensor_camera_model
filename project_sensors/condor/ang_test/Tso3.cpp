#include "Tso3.h"

#include "../ang/rotate/dcm.h"
#include "../ang/rotate/rodrigues.h"
#include "../ang/rotate/euler.h"
#include "../ang/rotate/rotv.h"
#include "../ang/rotate/so3_tangent.h"
#include "../ang/tools.h"
#include <iostream>

using namespace std;

ang::test::Tso3::Tso3(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Tso3::run() {
	::jail::unit_test::run();

    test_so3();                            cout << endl << endl;
    test_euler();                          cout << endl << endl;
    test_exp_log_maps();                   cout << endl << endl;
    test_exp_log_maps_small();             cout << endl << endl;
    test_power();                          cout << endl << endl;
    test_slerp();                          cout << endl << endl;
    test_plus_minus();                     cout << endl << endl;
    test_quat_4_solutions();               cout << endl << endl;
    test_rotv_euclidean_diff();            cout << endl << endl;
    test_adjoint();                        cout << endl << endl;
    test_point_velocity();                 cout << endl << endl;

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_so3() {
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
    vector<double> Vangle_rad(16);
    Vangle_rad[0]  = 0.7 * math::constant::PI(); // between 90 and 180 deg
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
    Vangle_rad[11] = 1.4 * math::constant::PI(); // between 180 and 270 deg
    Vangle_rad[12] = 1.9 * math::constant::PI(); // between 270 and 360 deg
    Vangle_rad[13] = 2.0 * math::constant::PI() - 1e-9;
    Vangle_rad[14] = 2.0 * math::constant::PI() - 1e-15;
    Vangle_rad[14] = 2.0 * math::constant::PI() - 1e-16;

    Eigen::Vector3d  dir_unit;
    for (unsigned int i = 0; i != Vdir_full.size(); ++i) {
        dir_unit = Vdir_full[i] / Vdir_full[i].norm();
        for (unsigned int j = 0; j != Vangle_rad.size(); ++j) {

            cout << endl << "Iteration i = " << i << "  j = " << j << endl;

            // back and forth from rotation vector to all other representations
            ang::rotv rv1_nb(dir_unit * Vangle_rad[j]);
            ang::dcm R1_nb(rv1_nb);
            ang::rodrigues q1_nb(rv1_nb);
            ang::euler euler1_nb(rv1_nb);

            ang::rotv rv2_nb(R1_nb);
            ang::rotv rv3_nb(q1_nb);
            ang::rotv rv4_nb(euler1_nb);

            // back and forth from rotation matrix to all other representations
            ang::rodrigues q2_nb(R1_nb);
            ang::euler euler2_nb(R1_nb);

            ang::dcm R2_nb(q2_nb);
            ang::dcm R3_nb(euler2_nb);

            // back and forth from quaternion to all other representations
            ang::euler euler3_nb(q1_nb);

            ang::rodrigues q3_nb(euler3_nb);

            check("rotv - rotv - 1a         ", rv1_nb()(0), rv2_nb()(0), 1e-12);
            check("rotv - rotv - 2a         ", rv1_nb()(1), rv2_nb()(1), 1e-12);
            check("rotv - rotv - 3a         ", rv1_nb()(2), rv2_nb()(2), 1e-12);

            check("rotv - rotv - 1b         ", rv1_nb()(0), rv3_nb()(0), 1e-12);
            check("rotv - rotv - 2b         ", rv1_nb()(1), rv3_nb()(1), 1e-12);
            check("rotv - rotv - 3b         ", rv1_nb()(2), rv3_nb()(2), 1e-12);

            check("rotv - rotv - 1c         ", rv1_nb()(0), rv4_nb()(0), 1e-12);
            check("rotv - rotv - 2c         ", rv1_nb()(1), rv4_nb()(1), 1e-12);
            check("rotv - rotv - 3c         ", rv1_nb()(2), rv4_nb()(2), 1e-12);

            check("dcm - dcm - 11a          ", R1_nb()(0,0), R2_nb()(0,0), 1e-12);
            check("dcm - dcm - 12a          ", R1_nb()(0,1), R2_nb()(0,1), 1e-12);
            check("dcm - dcm - 13a          ", R1_nb()(0,2), R2_nb()(0,2), 1e-12);
            check("dcm - dcm - 21a          ", R1_nb()(1,0), R2_nb()(1,0), 1e-12);
            check("dcm - dcm - 22a          ", R1_nb()(1,1), R2_nb()(1,1), 1e-12);
            check("dcm - dcm - 23a          ", R1_nb()(1,2), R2_nb()(1,2), 1e-12);
            check("dcm - dcm - 31a          ", R1_nb()(2,0), R2_nb()(2,0), 1e-12);
            check("dcm - dcm - 32a          ", R1_nb()(2,1), R2_nb()(2,1), 1e-12);
            check("dcm - dcm - 33a          ", R1_nb()(2,2), R2_nb()(2,2), 1e-12);

            check("dcm - dcm - 11b          ", R1_nb()(0,0), R3_nb()(0,0), 1e-12);
            check("dcm - dcm - 12b          ", R1_nb()(0,1), R3_nb()(0,1), 1e-12);
            check("dcm - dcm - 13b          ", R1_nb()(0,2), R3_nb()(0,2), 1e-12);
            check("dcm - dcm - 21b          ", R1_nb()(1,0), R3_nb()(1,0), 1e-12);
            check("dcm - dcm - 22b          ", R1_nb()(1,1), R3_nb()(1,1), 1e-12);
            check("dcm - dcm - 23b          ", R1_nb()(1,2), R3_nb()(1,2), 1e-12);
            check("dcm - dcm - 31b          ", R1_nb()(2,0), R3_nb()(2,0), 1e-12);
            check("dcm - dcm - 32b          ", R1_nb()(2,1), R3_nb()(2,1), 1e-12);
            check("dcm - dcm - 33b          ", R1_nb()(2,2), R3_nb()(2,2), 1e-12);

            check("quat - quat - 0a         ", q1_nb()(0), q2_nb()(0), 1e-12);
            check("quat - quat - 1a         ", q1_nb()(1), q2_nb()(1), 1e-12);
            check("quat - quat - 2a         ", q1_nb()(2), q2_nb()(2), 1e-12);
            check("quat - quat - 3a         ", q1_nb()(3), q2_nb()(3), 1e-12);

            check("quat - quat - 0b         ", q1_nb()(0), q3_nb()(0), 1e-12);
            check("quat - quat - 1b         ", q1_nb()(1), q3_nb()(1), 1e-12);
            check("quat - quat - 2b         ", q1_nb()(2), q3_nb()(2), 1e-12);
            check("quat - quat - 3b         ", q1_nb()(3), q3_nb()(3), 1e-12);

            check("euler - euler - 1a       ", euler1_nb.get_yaw_rad(),   euler2_nb.get_yaw_rad(),   1e-12);
            check("euler - euler - 2a       ", euler1_nb.get_pitch_rad(), euler2_nb.get_pitch_rad(), 1e-12);
            check("euler - euler - 3a       ", euler1_nb.get_bank_rad(),  euler2_nb.get_bank_rad(),  1e-12);

            check("euler - euler - 1b       ", euler1_nb.get_yaw_rad(),   euler3_nb.get_yaw_rad(),   1e-12);
            check("euler - euler - 2b       ", euler1_nb.get_pitch_rad(), euler3_nb.get_pitch_rad(), 1e-12);
            check("euler - euler - 3b       ", euler1_nb.get_bank_rad(),  euler3_nb.get_bank_rad(), 1e-12);

            Eigen::Vector3d vec_n(-4.0, 3.0, 7.0);
            Eigen::Vector3d vec_b_rotv, vec_b_dcm, vec_b_rodrigues, vec_b_euler;
            Eigen::Vector3d vec_b_rotvP, vec_b_dcmP, vec_b_rodriguesP, vec_b_eulerP;

            vec_b_rotv         = rv1_nb / vec_n;
            vec_b_dcm          = R1_nb / vec_n;
            vec_b_rodrigues    = q1_nb / vec_n;
            vec_b_euler        = euler1_nb / vec_n;

            vec_b_rotvP         = rv1_nb.inverse() * vec_n;
            vec_b_dcmP          = R1_nb.inverse() * vec_n;
            vec_b_rodriguesP    = q1_nb.inverse() * vec_n;
            //vec_b_eulerP        = euler1_nb.inverse() * vec_n; // not implemented

            check("method / - 1a             ", vec_b_rotv(0), vec_b_dcm(0),       1e-10);
            check("method / - 2a             ", vec_b_rotv(1), vec_b_dcm(1),       1e-10);
            check("method / - 3a             ", vec_b_rotv(2), vec_b_dcm(2),       1e-10);
            check("method / - 1b             ", vec_b_rotv(0), vec_b_rodrigues(0), 1e-10);
            check("method / - 2b             ", vec_b_rotv(1), vec_b_rodrigues(1), 1e-10);
            check("method / - 3b             ", vec_b_rotv(2), vec_b_rodrigues(2), 1e-10);
            check("method / - 1c             ", vec_b_rotv(0), vec_b_euler(0),     1e-10);
            check("method / - 2c             ", vec_b_rotv(1), vec_b_euler(1),     1e-10);
            check("method / - 3c             ", vec_b_rotv(2), vec_b_euler(2),     1e-10);

            check("method inverse and * - 1a ", vec_b_rotv(0), vec_b_rotvP(0),      1e-10);
            check("method inverse and * - 2a ", vec_b_rotv(1), vec_b_rotvP(1),      1e-10);
            check("method inverse and * - 3a ", vec_b_rotv(2), vec_b_rotvP(2),      1e-10);
            check("method inverse and * - 1b ", vec_b_rotv(0), vec_b_dcmP(0),       1e-10);
            check("method inverse and * - 2b ", vec_b_rotv(1), vec_b_dcmP(1),       1e-10);
            check("method inverse and * - 3b ", vec_b_rotv(2), vec_b_dcmP(2),       1e-10);
            check("method inverse and * - 1c ", vec_b_rotv(0), vec_b_rodriguesP(0), 1e-10);
            check("method inverse and * - 2c ", vec_b_rotv(1), vec_b_rodriguesP(1), 1e-10);
            check("method inverse and * - 3c ", vec_b_rotv(2), vec_b_rodriguesP(2), 1e-10);

            Eigen::Vector3d vec_n_rotv, vec_n_dcm, vec_n_rodrigues, vec_n_euler;
            Eigen::Vector3d vec_n_rotvP, vec_n_dcmP, vec_n_rodriguesP, vec_n_eulerP;

            vec_n_rotv         = rv1_nb * vec_b_rotv;
            vec_n_dcm          = R1_nb * vec_b_dcm;
            vec_n_rodrigues    = q1_nb * vec_b_rodrigues;
            vec_n_euler        = euler1_nb * vec_b_euler;

            vec_n_rotvP         = rv1_nb.inverse() / vec_b_rotvP;
            vec_n_dcmP          = R1_nb.inverse() / vec_b_dcmP;
            vec_n_rodriguesP    = q1_nb.inverse() / vec_b_rodriguesP;
            //vec_n_eulerP        = euler_nb.inverse() / vec_b_eulerP; // not implemented

            check("method * - 1a            ", vec_n(0), vec_n_rotv(0), 1e-12);
            check("method * - 2a            ", vec_n(1), vec_n_rotv(1), 1e-12);
            check("method * - 3a            ", vec_n(2), vec_n_rotv(2), 1e-12);
            check("method * - 1b            ", vec_n(0), vec_n_dcm(0), 1e-11);
            check("method * - 2b            ", vec_n(1), vec_n_dcm(1), 1e-11);
            check("method * - 3b            ", vec_n(2), vec_n_dcm(2), 1e-11);
            check("method * - 1c            ", vec_n(0), vec_n_rodrigues(0), 1e-12);
            check("method * - 2c            ", vec_n(1), vec_n_rodrigues(1), 1e-12);
            check("method * - 3c            ", vec_n(2), vec_n_rodrigues(2), 1e-12);
            check("method * - 1d            ", vec_n(0), vec_n_euler(0), 1e-12);
            check("method * - 2d            ", vec_n(1), vec_n_euler(1), 1e-12);
            check("method * - 3d            ", vec_n(2), vec_n_euler(2), 1e-12);

            check("method inverse and / 1a  ", vec_n(0), vec_n_rotvP(0), 1e-12);
            check("method inverse and / 2a  ", vec_n(1), vec_n_rotvP(1), 1e-12);
            check("method inverse and / 3a  ", vec_n(2), vec_n_rotvP(2), 1e-12);
            check("method inverse and / 1b  ", vec_n(0), vec_n_dcmP(0), 1e-11);
            check("method inverse and / 2b  ", vec_n(1), vec_n_dcmP(1), 1e-11);
            check("method inverse and / 3b  ", vec_n(2), vec_n_dcmP(2), 1e-11);
            check("method inverse and / 1c  ", vec_n(0), vec_n_rodriguesP(0), 1e-12);
            check("method inverse and / 2c  ", vec_n(1), vec_n_rodriguesP(1), 1e-12);
            check("method inverse and / 3c  ", vec_n(2), vec_n_rodriguesP(2), 1e-12);

            ang::so3_tangent w_nbb_rps(0.2*d2r, -0.1*d2r, 0.3*d2r);

            Eigen::Vector3d rvdot1_nb, eulerdot1_nb;
            Eigen::Matrix3d Rdot1_nb;
            ang::quat qdot1_nb;
            ang::so3_tangent w_nbb_rps_rotv, w_nbb_rps_dcm, w_nbb_rps_rodrigues, w_nbb_rps_euler;

            rvdot1_nb        = rv1_nb.omegabody2dot(w_nbb_rps);
            Rdot1_nb         = R1_nb.omegabody2dot(w_nbb_rps);
            qdot1_nb         = q1_nb.omegabody2dot(w_nbb_rps);
            eulerdot1_nb     = euler1_nb.omegabody2dot(w_nbb_rps);

            w_nbb_rps_rotv      = rv1_nb.dot2omegabody(rvdot1_nb);
            w_nbb_rps_dcm       = R1_nb.dot2omegabody(Rdot1_nb);
            w_nbb_rps_rodrigues = q1_nb.dot2omegabody(qdot1_nb);
            w_nbb_rps_euler     = euler1_nb.dot2omegabody(eulerdot1_nb);

            check("omega body - 1a    ", w_nbb_rps()(0), w_nbb_rps_rotv()(0), 1e-10);
            check("omega body - 2a    ", w_nbb_rps()(1), w_nbb_rps_rotv()(1), 1e-10);
            check("omega body - 3a    ", w_nbb_rps()(2), w_nbb_rps_rotv()(2), 1e-10);
            check("omega body - 1b    ", w_nbb_rps()(0), w_nbb_rps_dcm()(0), 1e-12);
            check("omega body - 2b    ", w_nbb_rps()(1), w_nbb_rps_dcm()(1), 1e-12);
            check("omega body - 3b    ", w_nbb_rps()(2), w_nbb_rps_dcm()(2), 1e-12);
            check("omega body - 1c    ", w_nbb_rps()(0), w_nbb_rps_rodrigues()(0), 1e-12);
            check("omega body - 2c    ", w_nbb_rps()(1), w_nbb_rps_rodrigues()(1), 1e-12);
            check("omega body - 3c    ", w_nbb_rps()(2), w_nbb_rps_rodrigues()(2), 1e-12);
            check("omega body - 1d    ", w_nbb_rps()(0), w_nbb_rps_euler()(0), 1e-12);
            check("omega body - 2d    ", w_nbb_rps()(1), w_nbb_rps_euler()(1), 1e-12);
            check("omega body - 3d    ", w_nbb_rps()(2), w_nbb_rps_euler()(2), 1e-12);

            ang::so3_tangent w_nbn_rps = q1_nb | w_nbb_rps;

            Eigen::Vector3d rvdot2_nb, eulerdot2_nb;
            Eigen::Matrix3d Rdot2_nb;
            ang::quat qdot2_nb;
            ang::so3_tangent w_nbn_rps_rotv, w_nbn_rps_dcm, w_nbn_rps_rodrigues, w_nbn_rps_euler;

            rvdot2_nb       = rv1_nb.omegaspace2dot(w_nbn_rps);
            Rdot2_nb        = R1_nb.omegaspace2dot(w_nbn_rps);
            qdot2_nb        = q1_nb.omegaspace2dot(w_nbn_rps);
            eulerdot2_nb    = euler1_nb.omegaspace2dot(w_nbn_rps);

            w_nbn_rps_rotv      = rv1_nb.dot2omegaspace(rvdot2_nb);
            w_nbn_rps_dcm       = R1_nb.dot2omegaspace(Rdot2_nb);
            w_nbn_rps_rodrigues = q1_nb.dot2omegaspace(qdot2_nb);
            w_nbn_rps_euler     = euler1_nb.dot2omegaspace(eulerdot2_nb);

            check("omega body-space rotv1   ", rvdot1_nb(0), rvdot2_nb(0), 1e-10);
            check("omega body-space rotv2   ", rvdot1_nb(1), rvdot2_nb(1), 1e-10);
            check("omega body-space rotv3   ", rvdot1_nb(2), rvdot2_nb(2), 1e-10);
            check("omega body-space dcm11   ", Rdot1_nb(0,0), Rdot2_nb(0,0), 1e-12);
            check("omega body-space dcm12   ", Rdot1_nb(0,1), Rdot2_nb(0,1), 1e-12);
            check("omega body-space dcm13   ", Rdot1_nb(0,2), Rdot2_nb(0,2), 1e-12);
            check("omega body-space dcm21   ", Rdot1_nb(1,0), Rdot2_nb(1,0), 1e-12);
            check("omega body-space dcm22   ", Rdot1_nb(1,1), Rdot2_nb(1,1), 1e-12);
            check("omega body-space dcm23   ", Rdot1_nb(1,2), Rdot2_nb(1,2), 1e-12);
            check("omega body-space dcm31   ", Rdot1_nb(2,0), Rdot2_nb(2,0), 1e-12);
            check("omega body-space dcm32   ", Rdot1_nb(2,1), Rdot2_nb(2,1), 1e-12);
            check("omega body-space dcm33   ", Rdot1_nb(2,2), Rdot2_nb(2,2), 1e-12);
            check("omega body-space rodr0   ", qdot1_nb(0), qdot2_nb(0), 1e-12);
            check("omega body-space rodr1   ", qdot1_nb(1), qdot2_nb(1), 1e-12);
            check("omega body-space rodr2   ", qdot1_nb(2), qdot2_nb(2), 1e-12);
            check("omega body-space rodr3   ", qdot1_nb(3), qdot2_nb(3), 1e-12);
            check("omega body-space euler1  ", eulerdot1_nb(0), eulerdot2_nb(0), 1e-12);
            check("omega body-space euler2  ", eulerdot1_nb(1), eulerdot2_nb(1), 1e-12);
            check("omega body-space euler3  ", eulerdot1_nb(2), eulerdot2_nb(2), 1e-12);

            check("omega space - 1a         ", w_nbn_rps()(0), w_nbn_rps_rotv()(0), 1e-10);
            check("omega space - 2a         ", w_nbn_rps()(1), w_nbn_rps_rotv()(1), 1e-10);
            check("omega space - 3a         ", w_nbn_rps()(2), w_nbn_rps_rotv()(2), 1e-10);
            check("omega space - 1b         ", w_nbn_rps()(0), w_nbn_rps_dcm()(0), 1e-12);
            check("omega space - 2b         ", w_nbn_rps()(1), w_nbn_rps_dcm()(1), 1e-12);
            check("omega space - 3b         ", w_nbn_rps()(2), w_nbn_rps_dcm()(2), 1e-12);
            check("omega space - 1c         ", w_nbn_rps()(0), w_nbn_rps_rodrigues()(0), 1e-12);
            check("omega space - 2c         ", w_nbn_rps()(1), w_nbn_rps_rodrigues()(1), 1e-12);
            check("omega space - 3c         ", w_nbn_rps()(2), w_nbn_rps_rodrigues()(2), 1e-12);
            check("omega space - 1d         ", w_nbn_rps()(0), w_nbn_rps_euler()(0), 1e-12);
            check("omega space - 2d         ", w_nbn_rps()(1), w_nbn_rps_euler()(1), 1e-12);
            check("omega space - 3d         ", w_nbn_rps()(2), w_nbn_rps_euler()(2), 1e-12);

            ang::rotv rv_bc(-0.1, 0.15, 0.7);
            ang::dcm R_bc(rv_bc);
            ang::rodrigues q_bc(rv_bc);
            ang::euler euler_bc(rv_bc);

            ang::rotv rv_nc, rv_nc_dcm, rv_nc_rodrigues;
            ang::dcm R_nc;
            ang::rodrigues q_nc;

            rv_nc  = rv1_nb * rv_bc;
            R_nc   = R1_nb * R_bc;
            q_nc   = q1_nb * q_bc;

            rv_nc_dcm       = ang::rotv(R_nc);
            rv_nc_rodrigues = ang::rotv(q_nc);

            check("combination rotv-dcm 1      ", rv_nc()(0), rv_nc_dcm()(0), 1e-12);
            check("combination rotv-dcm 2      ", rv_nc()(1), rv_nc_dcm()(1), 1e-12);
            check("combination rotv-dcm 3      ", rv_nc()(2), rv_nc_dcm()(2), 1e-12);
            check("combination rotv-rodrigues 1", rv_nc()(0), rv_nc_rodrigues()(0), 1e-12);
            check("combination rotv-rodrigues 2", rv_nc()(1), rv_nc_rodrigues()(1), 1e-12);
            check("combination rotv-rodrigues 3", rv_nc()(2), rv_nc_rodrigues()(2), 1e-12);

            ang::rotv rv_ncP, rv_ncP_dcm, rv_ncP_rodrigues;
            ang::dcm R_ncP;
            ang::rodrigues q_ncP;

            rv_ncP  = rv1_nb.inverse() / rv_bc;
            R_ncP   = R1_nb.inverse() / R_bc;
            q_ncP   = q1_nb.inverse() / q_bc;

            rv_ncP_dcm       = ang::rotv(R_ncP);
            rv_ncP_rodrigues = ang::rotv(q_ncP);

            check("combination / rotv-rotv 1     ", rv_ncP()(0), rv_nc()(0), 1e-12);
            check("combination / rotv-rotv 2     ", rv_ncP()(1), rv_nc()(1), 1e-12);
            check("combination / rotv-rotv 3     ", rv_ncP()(2), rv_nc()(2), 1e-12);
            check("combination / rotv-dcm 1      ", rv_ncP()(0), rv_ncP_dcm()(0), 1e-12);
            check("combination / rotv-dcm 2      ", rv_ncP()(1), rv_ncP_dcm()(1), 1e-12);
            check("combination / rotv-dcm 3      ", rv_ncP()(2), rv_ncP_dcm()(2), 1e-12);
            check("combination / rotv-rodrigues 1", rv_ncP()(0), rv_ncP_rodrigues()(0), 1e-12);
            check("combination / rotv-rodrigues 2", rv_ncP()(1), rv_ncP_rodrigues()(1), 1e-12);
            check("combination / rotv-rodrigues 3", rv_ncP()(2), rv_ncP_rodrigues()(2), 1e-12);

            ang::rotv Zrv2 = rv1_nb.exp_map_rodrigues().log_map();
            ang::rotv Zrv3 = rv1_nb.exp_map_dcm().log_map();

            check("exp is the opposite of log 1 ", rv1_nb()(0), Zrv2()(0), 1e-12);
            check("exp is the opposite of log 2 ", rv1_nb()(1), Zrv2()(1), 1e-12);
            check("exp is the opposite of log 3 ", rv1_nb()(2), Zrv2()(2), 1e-12);
            check("exp is the opposite of log 1 ", rv1_nb()(0), Zrv3()(0), 1e-12);
            check("exp is the opposite of log 2 ", rv1_nb()(1), Zrv3()(1), 1e-12);
            check("exp is the opposite of log 3 ", rv1_nb()(2), Zrv3()(2), 1e-12);
        }
    }
} // closes test_rotation

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_euler() {

    /* It checks the proper behavior of the "euler" class methods related with the obtainment
     * of different sets of Euler angles (bfsned, bfswfs, wfsned).
     */

    double d2r = math::constant::D2R();

    //ang::euler euler_wb(-10.0 * d2r, -5.0 * d2r, 0.0 * d2r);
    //ang::euler euler_nw(37.0 *d2r, 2.0*d2r, 11.0*d2r);

    //ang::euler euler_nb;
    //ang::euler::obtain_euler_nedbfs(euler_nb, euler_nw, euler_wb);
    //ang::euler euler_wb1;
    //ang::euler::obtain_euler_wfsbfs(euler_wb1, euler_nb, euler_nw);
    //ang::euler euler_nw1;
    //ang::euler::obtain_euler_nedwfs(euler_nw1, euler_nb, euler_wb);

    //check("bfswfs - 1       ", euler_wb.get_yaw_rad(),   euler_wb1.get_yaw_rad(), 1e-12);
    //check("bfswfs - 2       ", euler_wb.get_pitch_rad(), euler_wb1.get_pitch_rad(), 1e-12);
    //check("bfswfs - 3       ", euler_wb.get_bank_rad(),  euler_wb1.get_bank_rad(), 1e-12);

    //check("wfsned - 1       ", euler_nw.get_yaw_rad(),   euler_nw1.get_yaw_rad(), 1e-12);
    //check("wfsned - 2       ", euler_nw.get_pitch_rad(), euler_nw1.get_pitch_rad(), 1e-12);
    //check("wfsned - 3       ", euler_nw.get_bank_rad(),  euler_nw1.get_bank_rad(), 1e-12);

    Eigen::Vector3d v_n_mps(50, -20, 8);
    ang::euler Aeuler_nb(-6 * d2r, 8 * d2r, -3 * d2r);
    ang::euler euler_ng;
    ang::euler::obtain_euler_nedgrd(euler_ng, Aeuler_nb, v_n_mps);

    check("grdned - 1       ", euler_ng.get_yaw_rad(),   -0.380506377112365, 1e-12);
    check("grdned - 2       ", euler_ng.get_pitch_rad(), -0.147477689089775, 1e-12);
    check("grdned - 3       ", euler_ng.get_bank_rad(),  -0.010715602969010, 1e-12);

} // closes test_euler

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_exp_log_maps() {

    // the purpose of this test is to verify that the exponential and logarithmic functions relating
    // the rotation vector with both the rotation matrix and the quaternion comply with the common
    // arithmetic rules for exponential and logarithmic functions.

    ang::rotv Orotv_ba(0.1, 0.2, -0.15);
    ang::rotv Orotv_cb(-0.02, 0.12, 0.07);
    ang::rotv Orotv_ca = Orotv_cb * Orotv_ba;
    ang::rotv Orotv_3ba = Orotv_ba * Orotv_ba * Orotv_ba;

    // check that exponential equals the inverse of logarithm --> à = log(exp(a))
    ang::rodrigues Oq_ba      = Orotv_ba.exp_map_rodrigues();
    ang::dcm Odcm_ba          = Orotv_ba.exp_map_dcm();
    ang::rotv Orotv1_ba       = Oq_ba.log_map();
    ang::rotv Orotv2_ba       = Odcm_ba.log_map();
    check("exp-log quat 0       ", Orotv_ba()(0), Orotv1_ba()(0), 1e-12);
    check("exp-log quat 1       ", Orotv_ba()(1), Orotv1_ba()(1), 1e-12);
    check("exp-log quat 2       ", Orotv_ba()(2), Orotv1_ba()(2), 1e-12);
    check("exp-log dcm 0        ", Orotv_ba()(0), Orotv2_ba()(0), 1e-12);
    check("exp-log dcm 1        ", Orotv_ba()(1), Orotv2_ba()(1), 1e-12);
    check("exp-log dcm 2        ", Orotv_ba()(2), Orotv2_ba()(2), 1e-12);

    // check that exponential of a sum equals the product of exponentials --> exp(a+b) = exp(a) * exp(b)
    ang::rodrigues Oq_cb  = Orotv_cb.exp_map_rodrigues();
    ang::rodrigues Oq_ca  = Orotv_ca.exp_map_rodrigues();
    ang::rodrigues Oq1_ca = Oq_cb * Oq_ba;
    check("exp of sum - quat 0   ", Oq_ca()(0), Oq1_ca()(0), 1e-12);
    check("exp of sum - quat 1   ", Oq_ca()(1), Oq1_ca()(1), 1e-12);
    check("exp of sum - quat 2   ", Oq_ca()(2), Oq1_ca()(2), 1e-12);
    check("exp of sum - quat 3   ", Oq_ca()(3), Oq1_ca()(3), 1e-12);

    ang::dcm Odcm_cb  = Orotv_cb.exp_map_dcm();
    ang::dcm Odcm_ca  = Orotv_ca.exp_map_dcm();
    ang::dcm Odcm1_ca = Odcm_cb * Odcm_ba;
    check("exp of sum - dcm 00  ", Odcm_ca()(0,0), Odcm1_ca()(0,0), 1e-12);
    check("exp of sum - dcm 01  ", Odcm_ca()(0,1), Odcm1_ca()(0,1), 1e-12);
    check("exp of sum - dcm 02  ", Odcm_ca()(0,2), Odcm1_ca()(0,2), 1e-12);
    check("exp of sum - dcm 10  ", Odcm_ca()(1,0), Odcm1_ca()(1,0), 1e-12);
    check("exp of sum - dcm 11  ", Odcm_ca()(1,1), Odcm1_ca()(1,1), 1e-12);
    check("exp of sum - dcm 12  ", Odcm_ca()(1,2), Odcm1_ca()(1,2), 1e-12);
    check("exp of sum - dcm 20  ", Odcm_ca()(2,0), Odcm1_ca()(2,0), 1e-12);
    check("exp of sum - dcm 21  ", Odcm_ca()(2,1), Odcm1_ca()(2,1), 1e-12);
    check("exp of sum - dcm 22  ", Odcm_ca()(2,2), Odcm1_ca()(2,2), 1e-12);

    // check that exponential of multiple equals the power of exponentials --> exp(n*a) = (exp(a))^n
    ang::rodrigues Oq_3ba  = Orotv_3ba.exp_map_rodrigues();
    ang::rodrigues Oq1_3ba = Oq_ba * Oq_ba * Oq_ba;
    check("exp of multiple - quat 0 ", Oq_3ba()(0), Oq1_3ba()(0), 1e-12);
    check("exp of multiple - quat 1 ", Oq_3ba()(1), Oq1_3ba()(1), 1e-12);
    check("exp of multiple - quat 2 ", Oq_3ba()(2), Oq1_3ba()(2), 1e-12);
    check("exp of multiple - quat 3 ", Oq_3ba()(3), Oq1_3ba()(3), 1e-12);

    ang::dcm Odcm_3ba  = Orotv_3ba.exp_map_dcm();
    ang::dcm Odcm1_3ba = Odcm_ba * Odcm_ba * Odcm_ba;
    check("exp of multiple - dcm 00 ", Odcm_3ba()(0,0), Odcm1_3ba()(0,0), 1e-12);
    check("exp of multiple - dcm 01 ", Odcm_3ba()(0,1), Odcm1_3ba()(0,1), 1e-12);
    check("exp of multiple - dcm 02 ", Odcm_3ba()(0,2), Odcm1_3ba()(0,2), 1e-12);
    check("exp of multiple - dcm 10 ", Odcm_3ba()(1,0), Odcm1_3ba()(1,0), 1e-12);
    check("exp of multiple - dcm 11 ", Odcm_3ba()(1,1), Odcm1_3ba()(1,1), 1e-12);
    check("exp of multiple - dcm 12 ", Odcm_3ba()(1,2), Odcm1_3ba()(1,2), 1e-12);
    check("exp of multiple - dcm 20 ", Odcm_3ba()(2,0), Odcm1_3ba()(2,0), 1e-12);
    check("exp of multiple - dcm 21 ", Odcm_3ba()(2,1), Odcm1_3ba()(2,1), 1e-12);
    check("exp of multiple - dcm 22 ", Odcm_3ba()(2,2), Odcm1_3ba()(2,2), 1e-12);

    // check that logarithmic of product equals sum of logarithms --> log(a*b) = log(a) + log(b)
    ang::rotv Orotv5_ba = Oq_ba.log_map();
    ang::rotv Orotv5_cb = Oq_cb.log_map();
    ang::rotv Orotv5_ca = Orotv5_cb * Orotv5_ba;
    ang::rodrigues Oq5_ca = Oq_cb * Oq_ba;
    ang::rotv Orotv6_ca = Oq5_ca.log_map();
    check("log of product - quat 0 ", Orotv5_ca()(0), Orotv6_ca()(0), 1e-12);
    check("log of product - quat 1 ", Orotv5_ca()(1), Orotv6_ca()(1), 1e-12);
    check("log of product - quat 2 ", Orotv5_ca()(2), Orotv6_ca()(2), 1e-12);

    ang::rotv Orotv8_ba = Odcm_ba.log_map();
    ang::rotv Orotv8_cb = Odcm_cb.log_map();
    ang::rotv Orotv8_ca = Orotv8_cb * Orotv8_ba;
    ang::dcm Odcm8_ca = Odcm_cb * Odcm_ba;
    ang::rotv Orotv9_ca = Odcm8_ca.log_map();
    check("log of product - dcm 0  ", Orotv8_ca()(0), Orotv9_ca()(0), 1e-12);
    check("log of product - dcm 1  ", Orotv8_ca()(1), Orotv9_ca()(1), 1e-12);
    check("log of product - dcm 2  ", Orotv8_ca()(2), Orotv9_ca()(2), 1e-12);

    // check that logarithm of power equals multiple of logarithm
    ang::rodrigues Oq5_3ba = Oq_ba * Oq_ba * Oq_ba;
    ang::rotv Orotv5_3ba = Oq5_3ba.log_map();
    ang::rotv Orotv55_ba = Oq_ba.log_map();
    ang::rotv Orotv6_3ba = Orotv55_ba * Orotv55_ba * Orotv55_ba;
    check("log of power - quat 0  ", Orotv5_3ba()(0), Orotv6_3ba()(0), 1e-12);
    check("log of power - quat 1  ", Orotv5_3ba()(1), Orotv6_3ba()(1), 1e-12);
    check("log of power - quat 2  ", Orotv5_3ba()(2), Orotv6_3ba()(2), 1e-12);

    ang::dcm Odcm8_3ba = Odcm_ba * Odcm_ba * Odcm_ba;
    ang::rotv Orotv8_3ba = Odcm8_3ba.log_map();
    ang::rotv Orotv88_ba = Odcm_ba.log_map();
    ang::rotv Orotv9_3ba = Orotv88_ba * Orotv88_ba * Orotv88_ba;
    check("log of power - dcm 0   ", Orotv8_3ba()(0), Orotv9_3ba()(0), 1e-12);
    check("log of power - dcm 1   ", Orotv8_3ba()(1), Orotv9_3ba()(1), 1e-12);
    check("log of power - dcm 2   ", Orotv8_3ba()(2), Orotv9_3ba()(2), 1e-12);

} // closes test_exp_log_maps

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_exp_log_maps_small() {

    // The key to everything is the exponential and logarithmic maps between the rotation vector
    // and the quaternion. All others, when values are small, revert here.
    //
    // For the exponential map, there are two possible expressions:
    // E1. Based on sin and cos.
    // E2. Truncated for when rotation angle is small. The objective of this test is to determine
    //     where is the threshold.
    //
    // For the logarithmic map, I have three possible expressions:
    // L1. Based on arcos, that I used until I read Sola2017. It has two significant problems.
    //     1st is that it returns 0 very quick for small values and 2nd that it only returns
    //     between 0 and 180 deg. It should not be used under any circumstances.
    // L2. Based on atan, solves previous problems.
    // L3. Truncated for when rotation angle is small. The objective of this test is to determine
    //     where is the threshold.
    //
    // I believe a good criterion to evaluate the exponential map is the norm of the resulting
    // quaternion, which shall be one.
    //
    // I believe a good criterion to evaluate the logarithmic map is the norm of the resulting
    // rotation vector, which shall coincide with that of the initial rotation vector (before
    // the exponential map)

    // ===== ===== ===== RESULTS ===== ===== =====
    // ===========================================
    //              = Quat Norm Errors =    =====                      Rotation Vector Norm Errors          =====
    //      angle  q1vec norm     q1_norm     q2_norm   rv11_norm   rv21_norm   rv31_norm   rv12_norm   rv22_norm   rv23_norm  best exp()  best log()
    // +5.341e+00  +4.540e-01  +0.000e+00  -1.614e+00  +8.882e-16  +8.882e-16  +4.516e+00        +nan  -5.550e-01  +4.958e+00   E1 only     L1 & L2
    // +3.770e+00  +9.511e-01  +0.000e+00  -9.268e-02  +0.000e+00  +8.882e-16  -1.051e+01  -2.233e+00  -9.526e-01  +2.631e+00   E1 only     L1 & L2
    // +2.199e+00  +8.910e-01  +1.110e-16  +3.705e-02  +4.441e-16  +0.000e+00  +4.681e-01  -2.666e-02  -9.603e-02  -1.669e+00   E1 only     L1 & L2
    // +6.283e-01  +3.090e-01  +1.110e-16  +3.926e-04  +0.000e+00  +0.000e+00  +5.248e-02  -3.499e-05  -2.017e-04  +5.234e-02   E1 only     L1 & L2
    // +1.000e-03  +5.000e-04  +0.000e+00  +2.665e-15  +0.000e+00  +2.168e-19  +1.666e-07  -4.337e-19  -2.168e-18  +1.666e-07   E1 only     L1 & L2
    // +1.000e-05  +5.000e-06  +0.000e+00  +0.000e+00  +0.000e+00  +0.000e+00  +1.667e-11  +0.000e+00  +0.000e+00  +1.667e-11   E1 & E2     L1 & L2
    // +1.000e-07  +5.000e-08  +0.000e+00  +0.000e+00  +0.000e+00  +0.000e+00  +1.667e-15  +0.000e+00  +0.000e+00  +1.667e-15   E1 & E2     L1 & L2
    // +5.000e-08  +2.500e-08  +0.000e+00  +0.000e+00  +0.000e+00  +0.000e+00  +4.167e-16  +0.000e+00  +0.000e+00  +4.167e-16   E1 & E2     L1 & L2 & L3
    // +3.000e-08  +1.500e-08  +0.000e+00  +0.000e+00  +0.000e+00  +0.000e+00  +1.500e-16  +0.000e+00  +0.000e+00  +1.500e-16   E1 & E2     L1 & L2 & L3
    // +2.000e-08  +1.000e-08  +0.000e+00  +0.000e+00        -nan  -3.309e-24  +6.667e-17        -nan  -3.309e-24  +6.667e-17   E1 & E2     L2 & L3
    // +1.000e-08  +5.000e-09  +0.000e+00  +0.000e+00        -nan  -1.654e-24  +1.667e-17        -nan  -1.654e-24  +1.667e-17   E1 & E2     L2 & L3
    // +1.000e-09  +5.000e-10  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +1.667e-19        -nan  +0.000e+00  +1.667e-19   E1 & E2     L2 & L3
    // +1.000e-11  +5.000e-12  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +1.667e-23        -nan  +0.000e+00  +1.667e-23   E1 & E2     L2 & L3
    // +1.000e-13  +5.000e-14  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +1.666e-27        -nan  +0.000e+00  +1.666e-27   E1 & E2     L2 & L3
    // +1.000e-15  +5.000e-16  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +1.972e-31        -nan  +0.000e+00  +1.972e-31   E1 & E2     L2 & L3
    // +1.000e-16  +5.000e-17  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +0.000e+00        -nan  +0.000e+00  +0.000e+00   E1 & E2     L2 & L3
    // +0.000e+00  +0.000e+00        -nan  +0.000e+00        -nan        -nan        -nan        -nan        -nan  +0.000e+00   E2 only     L3 only

    // ===== ===== ===== CONCLUSIONS ===== ===== =====
    // ===============================================
    // Rotation                 angles (bigger /// smaller) than 3.0e-8 [rad] --> use (E1 /// E2)
    // Quaternion (vector part) angles (bigger /// smaller) than 1.5e-8 [rad] --> use (L2 /// L3)

    // ===== ===== ===== LATER CONCLUSIONS ===== ===== =====
    // =====================================================
    // Why do this if I can use the default expressions (E1 and L2) for every angle except when 0, where I have to use (E2 and L3).

    Eigen::Vector3d dir_full(1.,1.,1.);
    Eigen::Vector3d dir_unit = dir_full / dir_full.norm();

    // different rotation angles (less than 180 [deg] as otherwise shortest path is chosen)
    vector<double> Vangle_rad(17);
    Vangle_rad[0]  = 1.7 * math::constant::PI(); // more than  90 deg
    Vangle_rad[1]  = 1.2 * math::constant::PI(); // less than 90 deg
    Vangle_rad[2]  = 0.7 * math::constant::PI(); // more than  90 deg
    Vangle_rad[3]  = 0.2 * math::constant::PI(); // less than 90 deg
    Vangle_rad[4]  = 1e-3;
    Vangle_rad[5]  = 1e-5;
    Vangle_rad[6]  = 1e-7;
    Vangle_rad[7]  = 5e-8;
    Vangle_rad[8]  = 3e-8;
    Vangle_rad[9]  = 2e-8;
    Vangle_rad[10] = 1e-8;
    Vangle_rad[11] = 1e-9;
    Vangle_rad[12] = 1e-11;
    Vangle_rad[13] = 1e-13;
    Vangle_rad[14] = 1e-15;
    Vangle_rad[15] = 1e-16;
    Vangle_rad[16] = 0.;

    for (unsigned int j = 0; j != Vangle_rad.size(); ++j) {
        cout << endl << "Iteration j = " << j << endl;

        // small rotation vector
        double angle_rad = Vangle_rad[j];
        Eigen::Vector3d rv_nb = dir_unit * angle_rad;

        // compute quaternion based on exponential E1
        double cos_phihalf = cos(angle_rad/2.0);
        double sin_phihalf = sin(angle_rad/2.0);
        Eigen::Vector4d q1_nb; q1_nb << cos_phihalf, rv_nb * sin_phihalf / angle_rad;

        // compute quaternion based on exponential E2
        double a = 1.0 - pow(angle_rad,2.) / 8.;
        double b = 1.0 - pow(angle_rad,2.) / 24.;
        Eigen::Vector4d q2_nb; q2_nb << a, rv_nb * b * 0.5;

        // compute rotation vector based on log L1 (previously E1)
        double q1_nb_angle_rad = 2.0 * acos(q1_nb(0));
        Eigen::Vector3d rv11_nb = q1_nb.segment<3>(1) * 2.0 * q1_nb(0) * q1_nb_angle_rad / sin(q1_nb_angle_rad);

        // compute rotation vector based on log L2 (previously E1)
        double q1_nb_norm = q1_nb.segment<3>(1).norm();
        Eigen::Vector3d rv21_nb = q1_nb.segment<3>(1) * 2.0 * atan2(q1_nb_norm, q1_nb(0)) / q1_nb_norm;

        // compute rotation vector based on log L3 (previously E1)
        Eigen::Vector3d rv31_nb = q1_nb.segment<3>(1) * 2.0 / q1_nb(0) * (1.0 - q1_nb_norm / (3.0 * q1_nb(0) * q1_nb(0)));

        // compute rotation vector based on log L1 (previously L2)
        double q2_nb_angle_rad = 2.0 * acos(q2_nb(0));
        Eigen::Vector3d rv12_nb = q2_nb.segment<3>(1) * 2.0 * q2_nb(0) * q2_nb_angle_rad / sin(q2_nb_angle_rad);

        // compute rotation vector based on log L2 (previously E2)
        double q2_nb_norm = q2_nb.segment<3>(1).norm();
        Eigen::Vector3d rv22_nb = q2_nb.segment<3>(1) * 2.0 * atan2(q2_nb_norm, q2_nb(0)) / q2_nb_norm;

        // compute rotation vector based on log L3 (previously E1)
        Eigen::Vector3d rv32_nb = q2_nb.segment<3>(1) * 2.0 / q2_nb(0) * (1.0 - q2_nb_norm / (3.0 * q2_nb(0) * q2_nb(0)));

        cout << "ROTV --> QUAT" << endl;
        cout << "q1_nb    : " << showpos << scientific << setprecision(15) << setw(23) << q1_nb(0) << setprecision(15) << setw(23) << q1_nb(1) << setprecision(15) << setw(23) << q1_nb(2) << scientific << setprecision(15) << setw(23) << q1_nb(3)  << endl;
        cout << "q2_nb    : " << showpos << scientific << setprecision(15) << setw(23) << q2_nb(0) << setprecision(15) << setw(23) << q2_nb(1) << setprecision(15) << setw(23) << q2_nb(2) << scientific << setprecision(15) << setw(23) << q2_nb(3)  << endl;
        cout << "q norm   : " << showpos << scientific << setprecision(15) << setw(23) << 1.0 << endl;
        cout << "q1 norm  : " << showpos << scientific << setprecision(15) << setw(23) << q1_nb.norm() << endl;
        cout << "q2 norm  : " << showpos << scientific << setprecision(15) << setw(23) << q2_nb.norm() << endl; // shall be one
        cout << "BACK TO ROTV FROM q1_nb" << endl;
        cout << "rv_nb    : " << showpos << scientific << setprecision(15) << setw(23) << rv_nb(0)   << setprecision(15) << setw(23) << rv_nb(1)   << setprecision(15) << setw(23) << rv_nb(2)   << endl;
        cout << "rv11_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv11_nb(0) << setprecision(15) << setw(23) << rv11_nb(1) << setprecision(15) << setw(23) << rv11_nb(2) << endl;
        cout << "rv21_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv21_nb(0) << setprecision(15) << setw(23) << rv21_nb(1) << setprecision(15) << setw(23) << rv21_nb(2) << endl;
        cout << "rv31_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv31_nb(0) << setprecision(15) << setw(23) << rv31_nb(1) << setprecision(15) << setw(23) << rv31_nb(2) << endl;
        cout << "rv norm  : " << showpos << scientific << setprecision(15) << setw(20) << angle_rad << endl;
        cout << "rv11 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv11_nb.norm() << endl; // shall be equal to angle_rad
        cout << "rv21 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv21_nb.norm() << endl; // shall be equal to angle_rad
        cout << "rv31 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv31_nb.norm() << endl; // shall be equal to angle_rad
        cout << "BACK TO ROTV FROM q2_nb" << endl;
        cout << "rv_nb    : " << showpos << scientific << setprecision(15) << setw(23) << rv_nb(0)   << setprecision(15) << setw(23) << rv_nb(1)   << setprecision(15) << setw(23) << rv_nb(2)   << endl;
        cout << "rv12_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv12_nb(0) << setprecision(15) << setw(23) << rv12_nb(1) << setprecision(15) << setw(23) << rv12_nb(2)  << endl;
        cout << "rv22_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv22_nb(0) << setprecision(15) << setw(23) << rv22_nb(1) << setprecision(15) << setw(23) << rv22_nb(2)  << endl;
        cout << "rv32_nb  : " << showpos << scientific << setprecision(15) << setw(23) << rv32_nb(0) << setprecision(15) << setw(23) << rv32_nb(1) << setprecision(15) << setw(23) << rv32_nb(2)  << endl;
        cout << "rv norm  : " << showpos << scientific << setprecision(15) << setw(20) << angle_rad << endl;
        cout << "rv12 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv12_nb.norm() << endl; // shall be equal to angle_rad
        cout << "rv22 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv22_nb.norm() << endl; // shall be equal to angle_rad
        cout << "rv32 norm: " << showpos << scientific << setprecision(15) << setw(20) << rv32_nb.norm() << endl; // shall be equal to angle_rad
        cout << "IMPORTANT RESULTS" << endl;
        cout << "angle           : " << showpos << scientific << showpos << setprecision(3) << setw(11) << Vangle_rad[j] << endl;
        cout << "q1 vector norm  : " << showpos << scientific << showpos << setprecision(3) << setw(11) << q1_nb_norm    << endl;
        cout << "q2 vector norm  : " << showpos << scientific << showpos << setprecision(3) << setw(11) << q2_nb_norm    << endl;
        cout << "q1 norm error   : " << showpos << scientific << showpos << setprecision(3) << setw(11) << 1 - q1_nb.norm() << endl;
        cout << "q2 norm error   : " << showpos << scientific << showpos << setprecision(3) << setw(11) << 1 - q2_nb.norm() << endl;
        cout << "rv11 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv11_nb.norm() << endl;
        cout << "rv21 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv21_nb.norm() << endl;
        cout << "rv31 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv31_nb.norm() << endl;
        cout << "rv12 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv12_nb.norm() << endl;
        cout << "rv22 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv22_nb.norm() << endl;
        cout << "rv32 norm error : " << showpos << scientific << showpos << setprecision(3) << setw(11) << angle_rad - rv32_nb.norm() << endl;
        int A = 8;
    }
} // closes test_exp_log_maps_small

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_power() {

    // the purpose of this test is to verify that the power function of a rotation
    ang::rotv Orotv(0.2, -0.3, 0.12);
    ang::rodrigues Oq(Orotv);
    ang::dcm Odcm(Orotv);

    ang::rotv Orotv2 = Orotv.pow(0.3);
    ang::rodrigues Oq2 = Oq.pow(0.3);
    ang::dcm Odcm2 = Odcm.pow(0.3);

    ang::rotv Orotv22(Oq2);
    ang::rotv Orotv23(Odcm2);

    check("exp1 rotv - rodrigues - 1   ", Orotv2()(0), Orotv22()(0), 1e-12);
    check("exp1 rotv - rodrigues - 2   ", Orotv2()(1), Orotv22()(1), 1e-12);
    check("exp1 rotv - rodrigues - 3   ", Orotv2()(2), Orotv22()(2), 1e-12);
    check("exp1 rotv - dcm - 1         ", Orotv2()(0), Orotv23()(0), 1e-12);
    check("exp1 rotv - dcm - 2         ", Orotv2()(1), Orotv23()(1), 1e-12);
    check("exp1 rotv - dcm - 3         ", Orotv2()(2), Orotv23()(2), 1e-12);

    ang::rotv Orotv3 = Orotv.pow(3.3);
    ang::rodrigues Oq3 = Oq.pow(3.3);
    ang::dcm Odcm3 = Odcm.pow(3.3);
    ang::rotv Orotv32(Oq3);
    ang::rotv Orotv33(Odcm3);

    check("exp2 rotv - rodrigues - 1   ", Orotv3()(0), Orotv32()(0), 1e-12);
    check("exp2 rotv - rodrigues - 2   ", Orotv3()(1), Orotv32()(1), 1e-12);
    check("exp2 rotv - rodrigues - 3   ", Orotv3()(2), Orotv32()(2), 1e-12);
    check("exp2 rotv - dcm - 1         ", Orotv3()(0), Orotv33()(0), 1e-12);
    check("exp2 rotv - dcm - 2         ", Orotv3()(1), Orotv33()(1), 1e-12);
    check("exp2 rotv - dcm - 3         ", Orotv3()(2), Orotv33()(2), 1e-12);

    ang::rotv Orotv4 = Orotv.pow(-0.3);
    ang::rodrigues Oq4 = Oq.pow(-0.3);
    ang::dcm Odcm4 = Odcm.pow(-0.3);

    ang::rotv Orotv42(Oq4);
    ang::rotv Orotv43(Odcm4);

    check("exp3 rotv - rodrigues - 1   ", Orotv4()(0), Orotv42()(0), 1e-12);
    check("exp3 rotv - rodrigues - 2   ", Orotv4()(1), Orotv42()(1), 1e-12);
    check("exp3 rotv - rodrigues - 3   ", Orotv4()(2), Orotv42()(2), 1e-12);
    check("exp3 rotv - dcm - 1         ", Orotv4()(0), Orotv43()(0), 1e-12);
    check("exp3 rotv - dcm - 2         ", Orotv4()(1), Orotv43()(1), 1e-12);
    check("exp3 rotv - dcm - 3         ", Orotv4()(2), Orotv43()(2), 1e-12);

    ang::rotv Orotv5 = Orotv.pow(-3.3);
    ang::rodrigues Oq5 = Oq.pow(-3.3);
    ang::dcm Odcm5 = Odcm.pow(-3.3);

    ang::rotv Orotv52(Oq5);
    ang::rotv Orotv53(Odcm5);

    check("exp4 rotv - rodrigues - 1   ", Orotv5()(0), Orotv52()(0), 1e-12);
    check("exp4 rotv - rodrigues - 2   ", Orotv5()(1), Orotv52()(1), 1e-12);
    check("exp4 rotv - rodrigues - 3   ", Orotv5()(2), Orotv52()(2), 1e-12);
    check("exp4 rotv - dcm - 1         ", Orotv5()(0), Orotv53()(0), 1e-12);
    check("exp4 rotv - dcm - 2         ", Orotv5()(1), Orotv53()(1), 1e-12);
    check("exp4 rotv - dcm - 3         ", Orotv5()(2), Orotv53()(2), 1e-12);

    ang::rotv Orotv6 = Orotv.pow(0.);
    ang::rodrigues Oq6 = Oq.pow(0.);
    ang::dcm Odcm6 = Odcm.pow(0.);

    ang::rotv Orotv62(Oq6);
    ang::rotv Orotv63(Odcm6);

    check("exp5 rotv - rodrigues - 1   ", Orotv6()(0), Orotv62()(0), 1e-12);
    check("exp5 rotv - rodrigues - 2   ", Orotv6()(1), Orotv62()(1), 1e-12);
    check("exp5 rotv - rodrigues - 3   ", Orotv6()(2), Orotv62()(2), 1e-12);
    check("exp5 rotv - dcm - 1         ", Orotv6()(0), Orotv63()(0), 1e-12);
    check("exp5 rotv - dcm - 2         ", Orotv6()(1), Orotv63()(1), 1e-12);
    check("exp5 rotv - dcm - 3         ", Orotv6()(2), Orotv63()(2), 1e-12);

    check("exp5 rotv - rotv - 1        ", Orotv6()(0), 0., 1e-12);
    check("exp5 rotv - rotv - 2        ", Orotv6()(1), 0., 1e-12);
    check("exp5 rotv - rotv - 3        ", Orotv6()(2), 0., 1e-12);
    check("exp5 rotv - rodrigues - 1   ", Oq6()(0), 1., 1e-12);
    check("exp5 rotv - rodrigues - 2   ", Oq6()(1), 0., 1e-12);
    check("exp5 rotv - rodrigues - 3   ", Oq6()(2), 0., 1e-12);
    check("exp5 rotv - rodrigues - 4   ", Oq6()(3), 0., 1e-12);
    check("exp5 rotv - dcm - 11        ", Odcm6()(0,0), 1., 1e-12);
    check("exp5 rotv - dcm - 12        ", Odcm6()(0,1), 0., 1e-12);
    check("exp5 rotv - dcm - 13        ", Odcm6()(0,2), 0., 1e-12);
    check("exp5 rotv - dcm - 21        ", Odcm6()(1,0), 0., 1e-12);
    check("exp5 rotv - dcm - 22        ", Odcm6()(1,1), 1., 1e-12);
    check("exp5 rotv - dcm - 23        ", Odcm6()(1,2), 0., 1e-12);
    check("exp5 rotv - dcm - 31        ", Odcm6()(2,0), 0., 1e-12);
    check("exp5 rotv - dcm - 32        ", Odcm6()(2,1), 0., 1e-12);
    check("exp5 rotv - dcm - 33        ", Odcm6()(2,2), 1., 1e-12);

    ang::rotv Orotv7= Orotv.pow(1.);
    ang::rodrigues Oq7 = Oq.pow(1.);
    ang::dcm Odcm7 = Odcm.pow(1.);

    ang::rotv Orotv72(Oq7);
    ang::rotv Orotv73(Odcm7);

    check("exp6 rotv - rodrigues - 1   ", Orotv7()(0), Orotv72()(0), 1e-12);
    check("exp6 rotv - rodrigues - 2   ", Orotv7()(1), Orotv72()(1), 1e-12);
    check("exp6 rotv - rodrigues - 3   ", Orotv7()(2), Orotv72()(2), 1e-12);
    check("exp6 rotv - dcm - 1         ", Orotv7()(0), Orotv73()(0), 1e-12);
    check("exp6 rotv - dcm - 2         ", Orotv7()(1), Orotv73()(1), 1e-12);
    check("exp6 rotv - dcm - 3         ", Orotv7()(2), Orotv73()(2), 1e-12);

    check("exp6 rotv - rotv - 1        ", Orotv7()(0), Orotv()(0), 1e-12);
    check("exp6 rotv - rotv - 2        ", Orotv7()(1), Orotv()(1), 1e-12);
    check("exp6 rotv - rotv - 3        ", Orotv7()(2), Orotv()(2), 1e-12);
    check("exp6 rotv - rodrigues - 1   ", Oq7()(0), Oq()(0), 1e-12);
    check("exp6 rotv - rodrigues - 2   ", Oq7()(1), Oq()(1), 1e-12);
    check("exp6 rotv - rodrigues - 3   ", Oq7()(2), Oq()(2), 1e-12);
    check("exp6 rotv - rodrigues - 4   ", Oq7()(3), Oq()(3), 1e-12);
    check("exp6 rotv - dcm - 11        ", Odcm7()(0,0), Odcm()(0,0), 1e-12);
    check("exp6 rotv - dcm - 12        ", Odcm7()(0,1), Odcm()(0,1), 1e-12);
    check("exp6 rotv - dcm - 13        ", Odcm7()(0,2), Odcm()(0,2), 1e-12);
    check("exp6 rotv - dcm - 21        ", Odcm7()(1,0), Odcm()(1,0), 1e-12);
    check("exp6 rotv - dcm - 22        ", Odcm7()(1,1), Odcm()(1,1), 1e-12);
    check("exp6 rotv - dcm - 23        ", Odcm7()(1,2), Odcm()(1,2), 1e-12);
    check("exp6 rotv - dcm - 31        ", Odcm7()(2,0), Odcm()(2,0), 1e-12);
    check("exp6 rotv - dcm - 32        ", Odcm7()(2,1), Odcm()(2,1), 1e-12);
    check("exp6 rotv - dcm - 33        ", Odcm7()(2,2), Odcm()(2,2), 1e-12);

} // closes test_power

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_slerp() {

    // I make them multiples so I can manually compute the interpolated values
    ang::rotv Orotv00(0.05, 0.10, 0.15);
    ang::rotv Orotv22(0.10, 0.20, 0.30);
    ang::rotv Orotv11(0.06, 0.12, 0.18);

    ang::rotv Orotv0 = ang::rotv::slerp(Orotv00, Orotv22, 0.);
    ang::rotv Orotv1 = ang::rotv::slerp(Orotv00, Orotv22, 0.2);
    ang::rotv Orotv2 = ang::rotv::slerp(Orotv00, Orotv22, 1.);

    check("slerp rotv 00   ", Orotv00()(0), Orotv0()(0), 1e-12);
    check("slerp rotv 01   ", Orotv00()(1), Orotv0()(1), 1e-12);
    check("slerp rotv 02   ", Orotv00()(2), Orotv0()(2), 1e-12);
    check("slerp rotv 10   ", Orotv11()(0), Orotv1()(0), 1e-12);
    check("slerp rotv 11   ", Orotv11()(1), Orotv1()(1), 1e-12);
    check("slerp rotv 12   ", Orotv11()(2), Orotv1()(2), 1e-12);
    check("slerp rotv 20   ", Orotv22()(0), Orotv2()(0), 1e-12);
    check("slerp rotv 21   ", Orotv22()(1), Orotv2()(1), 1e-12);
    check("slerp rotv 22   ", Orotv22()(2), Orotv2()(2), 1e-12);

    ang::rodrigues Oq00(Orotv00);
    ang::rodrigues Oq22(Orotv22);
    ang::rodrigues Oq11(Orotv11);

    ang::rodrigues Oq0 = ang::rodrigues::slerp(Oq00, Oq22, 0.);
    ang::rodrigues Oq1 = ang::rodrigues::slerp(Oq00, Oq22, 0.2);
    ang::rodrigues Oq2 = ang::rodrigues::slerp(Oq00, Oq22, 1.);

    check("slerp rodrig 00 ", Oq00()(0), Oq0()(0), 1e-12);
    check("slerp rodrig 01 ", Oq00()(1), Oq0()(1), 1e-12);
    check("slerp rodrig 02 ", Oq00()(2), Oq0()(2), 1e-12);
    check("slerp rodrig 03 ", Oq00()(3), Oq0()(3), 1e-12);
    check("slerp rodrig 10 ", Oq11()(0), Oq1()(0), 1e-12);
    check("slerp rodrig 11 ", Oq11()(1), Oq1()(1), 1e-12);
    check("slerp rodrig 12 ", Oq11()(2), Oq1()(2), 1e-12);
    check("slerp rodrig 13 ", Oq11()(3), Oq1()(3), 1e-12);
    check("slerp rodrig 20 ", Oq22()(0), Oq2()(0), 1e-12);
    check("slerp rodrig 21 ", Oq22()(1), Oq2()(1), 1e-12);
    check("slerp rodrig 22 ", Oq22()(2), Oq2()(2), 1e-12);
    check("slerp rodrig 23 ", Oq22()(3), Oq2()(3), 1e-12);

    ang::dcm R00(Orotv00);
    ang::dcm R22(Orotv22);
    ang::dcm R11(Orotv11);

    ang::dcm R0 = ang::dcm::slerp(R00, R22, 0.);
    ang::dcm R1 = ang::dcm::slerp(R00, R22, 0.2);
    ang::dcm R2 = ang::dcm::slerp(R00, R22, 1.);

    check("slerp dcm 0 00 ", R00()(0,0), R0()(0,0), 1e-12);
    check("slerp dcm 0 01 ", R00()(0,1), R0()(0,1), 1e-12);
    check("slerp dcm 0 02 ", R00()(0,2), R0()(0,2), 1e-12);
    check("slerp dcm 0 10 ", R00()(1,0), R0()(1,0), 1e-12);
    check("slerp dcm 0 11 ", R00()(1,1), R0()(1,1), 1e-12);
    check("slerp dcm 0 12 ", R00()(1,2), R0()(1,2), 1e-12);
    check("slerp dcm 0 20 ", R00()(2,0), R0()(2,0), 1e-12);
    check("slerp dcm 0 21 ", R00()(2,1), R0()(2,1), 1e-12);
    check("slerp dcm 0 22 ", R00()(2,2), R0()(2,2), 1e-12);

    check("slerp dcm 1 00 ", R11()(0,0), R1()(0,0), 1e-12);
    check("slerp dcm 1 01 ", R11()(0,1), R1()(0,1), 1e-12);
    check("slerp dcm 1 02 ", R11()(0,2), R1()(0,2), 1e-12);
    check("slerp dcm 1 10 ", R11()(1,0), R1()(1,0), 1e-12);
    check("slerp dcm 1 11 ", R11()(1,1), R1()(1,1), 1e-12);
    check("slerp dcm 1 12 ", R11()(1,2), R1()(1,2), 1e-12);
    check("slerp dcm 1 20 ", R11()(2,0), R1()(2,0), 1e-12);
    check("slerp dcm 1 21 ", R11()(2,1), R1()(2,1), 1e-12);
    check("slerp dcm 1 22 ", R11()(2,2), R1()(2,2), 1e-12);

    check("slerp dcm 2 00 ", R22()(0,0), R2()(0,0), 1e-12);
    check("slerp dcm 2 01 ", R22()(0,1), R2()(0,1), 1e-12);
    check("slerp dcm 2 02 ", R22()(0,2), R2()(0,2), 1e-12);
    check("slerp dcm 2 10 ", R22()(1,0), R2()(1,0), 1e-12);
    check("slerp dcm 2 11 ", R22()(1,1), R2()(1,1), 1e-12);
    check("slerp dcm 2 12 ", R22()(1,2), R2()(1,2), 1e-12);
    check("slerp dcm 2 20 ", R22()(2,0), R2()(2,0), 1e-12);
    check("slerp dcm 2 21 ", R22()(2,1), R2()(2,1), 1e-12);
    check("slerp dcm 2 22 ", R22()(2,2), R2()(2,2), 1e-12);

} // closes test_slerp

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_plus_minus() {

    // the purpose of this test is to verify that the plus and minus functions of a rotation
    ang::rotv r_a(0.2, -0.3, 0.12);
    ang::rodrigues q_a(r_a);
    ang::dcm R_a(r_a);

    ang::rotv Delta_r(-0.01, 0.03, 0.05);

    ang::rodrigues q_b1 = q_a.plus_right(Delta_r);
    ang::dcm R_b1 = R_a.plus_right(Delta_r);
    ang::rotv r_b1 = r_a.plus_right(Delta_r);

    ang::rotv r_bq1(q_b1);
    ang::rotv r_bR1(R_b1);

    check("plus right rodrigues - dcm 1      ", r_bq1()(0), r_bR1()(0), 1e-12);
    check("plus right rodrigues - dcm 2      ", r_bq1()(1), r_bR1()(1), 1e-12);
    check("plus right rodrigues - dcm 3      ", r_bq1()(2), r_bR1()(2), 1e-12);
    check("plus right rodrigues - rv 1       ", r_bq1()(0), r_b1()(0), 1e-12);
    check("plus right rodrigues - rv 2       ", r_bq1()(1), r_b1()(1), 1e-12);
    check("plus right rodrigues - rv 3       ", r_bq1()(2), r_b1()(2), 1e-12);

    ang::rotv r_cq1 = q_b1.minus_right(q_a);
    ang::rotv r_cR1 = R_b1.minus_right(R_a);
    ang::rotv r_c1  = r_b1.minus_right(r_a);

    check("minus right rodrigues 1           ", r_cq1()(0), Delta_r()(0), 1e-12);
    check("minus right rodrigues 2           ", r_cq1()(1), Delta_r()(1), 1e-12);
    check("minus right rodrigues 3           ", r_cq1()(2), Delta_r()(2), 1e-12);
    check("minus right dcm 1                 ", r_cR1()(0), Delta_r()(0), 1e-12);
    check("minus right dcm 2                 ", r_cR1()(1), Delta_r()(1), 1e-12);
    check("minus right dcm 3                 ", r_cR1()(2), Delta_r()(2), 1e-12);
    check("minus right rv 1                  ", r_c1()(0),  Delta_r()(0), 1e-12);
    check("minus right rv 2                  ", r_c1()(1),  Delta_r()(1), 1e-12);
    check("minus right rv 3                  ", r_c1()(2),  Delta_r()(2), 1e-12);

    ang::rodrigues q_b2 = q_a.plus_left(Delta_r);
    ang::dcm R_b2 = R_a.plus_left(Delta_r);
    ang::rotv r_b2 = r_a.plus_left(Delta_r);

    ang::rotv r_bq2(q_b2);
    ang::rotv r_bR2(R_b2);

    check("plus left rodrigues - dcm 1       ", r_bq2()(0), r_bR2()(0), 1e-12);
    check("plus left rodrigues - dcm 2       ", r_bq2()(1), r_bR2()(1), 1e-12);
    check("plus left rodrigues - dcm 3       ", r_bq2()(2), r_bR2()(2), 1e-12);
    check("plus left rodrigues - rv 1        ", r_bq2()(0), r_b2()(0), 1e-12);
    check("plus left rodrigues - rv 2        ", r_bq2()(1), r_b2()(1), 1e-12);
    check("plus left rodrigues - rv 3        ", r_bq2()(2), r_b2()(2), 1e-12);

    ang::rotv r_cq2 = q_b2.minus_left(q_a);
    ang::rotv r_cR2 = R_b2.minus_left(R_a);
    ang::rotv r_c2  = r_b2.minus_left(r_a);

    check("minus left rodrigues 1           ", r_cq2()(0), Delta_r()(0), 1e-12);
    check("minus left rodrigues 2           ", r_cq2()(1), Delta_r()(1), 1e-12);
    check("minus left rodrigues 3           ", r_cq2()(2), Delta_r()(2), 1e-12);
    check("minus left dcm 1                 ", r_cR2()(0), Delta_r()(0), 1e-12);
    check("minus left dcm 2                 ", r_cR2()(1), Delta_r()(1), 1e-12);
    check("minus left dcm 3                 ", r_cR2()(2), Delta_r()(2), 1e-12);
    check("minus left rv 1                  ", r_c2()(0), Delta_r()(0), 1e-12);
    check("minus left rv 2                  ", r_c2()(1), Delta_r()(1), 1e-12);
    check("minus left rv 3                  ", r_c2()(2), Delta_r()(2), 1e-12);

} // closes test_plus_minus

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_quat_4_solutions() {
    double d2r = math::constant::D2R();

    // Although there is no way to check it other than manually, the four
    // next rotations fall in the four different dcm to quat or euler to quat
    // cases.

    // The other correspond to singular cases

    ang::euler Oeuler01(30 * d2r, 40 * d2r, 50 * d2r);
    ang::dcm   Odcm01(Oeuler01);
    ang::rodrigues q01a(Oeuler01);
    ang::rodrigues q01b(Odcm01);

    ang::euler Oeuler02(300 * d2r, 300 * d2r, 150 * d2r);
    ang::dcm   Odcm02(Oeuler02);
    ang::rodrigues q02a(Oeuler02);
    ang::rodrigues q02b(Odcm02);

    ang::euler Oeuler03(-300 * d2r, -140 * d2r, +300 * d2r);
    ang::dcm   Odcm03(Oeuler03);
    ang::rodrigues q03a(Oeuler03);
    ang::rodrigues q03b(Odcm03);

    ang::euler Oeuler04(300 * d2r, 140 * d2r, 150 * d2r);
    ang::dcm   Odcm04(Oeuler04);
    ang::rodrigues q04a(Oeuler04);
    ang::rodrigues q04b(Odcm04);

    check("1st quat case - 0  ", q01a()(0), q01b()(0), 1e-12);
    check("1st quat case - 1  ", q01a()(1), q01b()(1), 1e-12);
    check("1st quat case - 2  ", q01a()(2), q01b()(2), 1e-12);
    check("1st quat case - 3  ", q01a()(3), q01b()(3), 1e-12);

    check("2nd quat case - 0  ", q02a()(0), q02b()(0), 1e-12);
    check("2nd quat case - 1  ", q02a()(1), q02b()(1), 1e-12);
    check("2nd quat case - 2  ", q02a()(2), q02b()(2), 1e-12);
    check("2nd quat case - 3  ", q02a()(3), q02b()(3), 1e-12);

    check("3rd quat case - 0  ", q03a()(0), q03b()(0), 1e-12);
    check("3rd quat case - 1  ", q03a()(1), q03b()(1), 1e-12);
    check("3rd quat case - 2  ", q03a()(2), q03b()(2), 1e-12);
    check("3rd quat case - 3  ", q03a()(3), q03b()(3), 1e-12);

    check("4th quat case - 0  ", q04a()(0), q04b()(0), 1e-12);
    check("4th quat case - 1  ", q04a()(1), q04b()(1), 1e-12);
    check("4th quat case - 2  ", q04a()(2), q04b()(2), 1e-12);
    check("4th quat case - 3  ", q04a()(3), q04b()(3), 1e-12);

    ang::euler Oeuler05(180 * d2r, 0 * d2r, 0 * d2r);
    ang::dcm   Odcm05(Oeuler05);
    ang::rodrigues q05a(Oeuler05);
    ang::rodrigues q05b(Odcm05);

    ang::euler Oeuler06(0 * d2r, 180 * d2r, 0 * d2r);
    ang::dcm   Odcm06(Oeuler06);
    ang::rodrigues q06a(Oeuler06);
    ang::rodrigues q06b(Odcm06);

    ang::euler Oeuler07(0 * d2r, 0 * d2r, 180 * d2r);
    ang::dcm   Odcm07(Oeuler07);
    ang::rodrigues q07a(Oeuler07);
    ang::rodrigues q07b(Odcm07);

    ang::euler Oeuler08(180 * d2r, 180 * d2r, 180 * d2r);
    ang::dcm   Odcm08(Oeuler08);
    ang::rodrigues q08a(Oeuler08);
    ang::rodrigues q08b(Odcm08);

    check("case - 0  ", q05a()(0), q05b()(0), 1e-12);
    check("case - 1  ", q05a()(1), q05b()(1), 1e-12);
    check("case - 2  ", q05a()(2), q05b()(2), 1e-12);
    check("case - 3  ", q05a()(3), q05b()(3), 1e-12);

    check("case - 0  ", q06a()(0), q06b()(0), 1e-12);
    check("case - 1  ", q06a()(1), q06b()(1), 1e-12);
    check("case - 2  ", q06a()(2), q06b()(2), 1e-12);
    check("case - 3  ", q06a()(3), q06b()(3), 1e-12);

    check("case - 0  ", q07a()(0), q07b()(0), 1e-12);
    check("case - 1  ", q07a()(1), q07b()(1), 1e-12);
    check("case - 2  ", q07a()(2), q07b()(2), 1e-12);
    check("case - 3  ", q07a()(3), q07b()(3), 1e-12);

    check("case - 0  ", q08a()(0), q08b()(0), 1e-12);
    check("case - 1  ", q08a()(1), q08b()(1), 1e-12);
    check("case - 2  ", q08a()(2), q08b()(2), 1e-12);
    check("case - 3  ", q08a()(3), q08b()(3), 1e-12);

} // closes test_quat_4_solutions

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_rotv_euclidean_diff() {

    // this test is intended to verify if there is some relation between the Euclidean modulus of the
    // Euclidean difference between two rotation vectors and the amount of the rotation between the same
    // two rotation vectors. It is true, this test verifies it, and hence the Euclidean norm can be
    // considered a good minimization objective.

    // Important note: When obtaining the Euclidean differences between two rotation vectors, ensure
    // angle is less than PI, which is the maximum. This is a problem as rotation vectors are not continuous

    double d2r = math::constant::D2R();
    double yaw0_deg   = 15.0;
    double pitch0_deg = 22.0;
    double roll0_deg  = 41.0;
    double Delta_phi_deg = 0.0001; // play with this, using 0.001,0.01, 0.1, 1, 10

    ang::euler euler0(yaw0_deg * d2r, pitch0_deg * d2r, roll0_deg * d2r);
    ang::rotv r0(euler0);

    ang::euler euler1((yaw0_deg + Delta_phi_deg) * d2r, pitch0_deg * d2r, roll0_deg * d2r);
    ang::rotv r1(euler1);
    ang::rotv r01 = r0 / r1;
    Eigen::Vector3d x01 = r0() - r1();
    double phi01_deg = r01().norm() / d2r;
    double mod01_deg = x01.norm() / d2r;
    cout << "phi01_deg:   " << fixed << setw(13) << setprecision(8) << showpos << phi01_deg << endl;
    cout << "mod01_deg:   " << fixed << setw(13) << setprecision(8) << showpos << mod01_deg << endl;

    ang::euler euler2(yaw0_deg * d2r, (pitch0_deg + Delta_phi_deg) * d2r, roll0_deg * d2r);
    ang::rotv r2(euler2);
    ang::rotv r02 = r0 / r2;
    Eigen::Vector3d x02 = r0() - r2();
    double phi02_deg = r02().norm() / d2r;
    double mod02_deg = x02.norm() / d2r;
    cout << "phi02_deg:   " << fixed << setw(13) << setprecision(8) << showpos << phi02_deg << endl;
    cout << "mod02_deg:   " << fixed << setw(13) << setprecision(8) << showpos << mod02_deg << endl;

    ang::euler euler3(yaw0_deg * d2r, pitch0_deg * d2r, (roll0_deg + Delta_phi_deg) * d2r);
    ang::rotv r3(euler3);
    ang::rotv r03 = r0 / r3;
    Eigen::Vector3d x03 = r0() - r3();
    double phi03_deg = r03().norm() / d2r;
    double mod03_deg = x03.norm() / d2r;
    cout << "phi03_deg:   " << fixed << setw(13) << setprecision(8) << showpos << phi03_deg << endl;
    cout << "mod03_deg:   " << fixed << setw(13) << setprecision(8) << showpos << mod03_deg << endl;

    ang::euler euler4((yaw0_deg + Delta_phi_deg) * d2r, (pitch0_deg + Delta_phi_deg) * d2r, (roll0_deg + Delta_phi_deg) * d2r);
    ang::rotv r4(euler4);
    ang::rotv r04 = r0 / r4;
    Eigen::Vector3d x04 = r0() - r4();
    double phi04_deg = r04().norm() / d2r;
    double mod04_deg = x04.norm() / d2r;
    cout << "phi04_deg:   " << fixed << setw(13) << setprecision(8) << showpos << phi04_deg << endl;
    cout << "mod04_deg:   " << fixed << setw(13) << setprecision(8) << showpos << mod04_deg << endl;

} // closes test_rotv

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_adjoint() {
    ang::so3_tangent w_nbb_0(0.51, -0.38, -0.14);
    ang::so3_tangent_skew w_nbb_0_skew = ang::so3_tangent::hat_skew(w_nbb_0);
    ang::so3_tangent_quat w_nbb_0_quat = ang::so3_tangent::hat_quat(w_nbb_0);

    ang::euler euler_nb(0.12, 0.17, -0.03);
    ang::dcm R_nb(euler_nb);
    ang::rodrigues q_nb(euler_nb);
    ang::rotv r_nb(euler_nb);

    ang::so3_tangent w_nbn_1(R_nb * w_nbb_0());
    ang::so3_tangent w_nbn_2(q_nb * w_nbb_0());
    ang::so3_tangent w_nbn_3(r_nb * w_nbb_0());
    ang::so3_tangent w_nbn_4 = R_nb | w_nbb_0;;
    ang::so3_tangent w_nbn_5 = q_nb | w_nbb_0;
    ang::so3_tangent w_nbn_6 = r_nb | w_nbb_0;
    ang::so3_tangent w_nbn_7 = ang::so3_tangent::wedge(R_nb | w_nbb_0_skew);
    ang::so3_tangent w_nbn_8 = ang::so3_tangent::wedge(q_nb | w_nbb_0_quat);

    check("adjoint 21        ", w_nbn_1()(0), w_nbn_2()(0), 1e-12);
    check("adjoint 22        ", w_nbn_1()(1), w_nbn_2()(1), 1e-12);
    check("adjoint 23        ", w_nbn_1()(2), w_nbn_2()(2), 1e-12);

    check("adjoint 31        ", w_nbn_1()(0), w_nbn_3()(0), 1e-12);
    check("adjoint 32        ", w_nbn_1()(1), w_nbn_3()(1), 1e-12);
    check("adjoint 33        ", w_nbn_1()(2), w_nbn_3()(2), 1e-12);

    check("adjoint 41        ", w_nbn_1()(0), w_nbn_4()(0), 1e-12);
    check("adjoint 42        ", w_nbn_1()(1), w_nbn_4()(1), 1e-12);
    check("adjoint 43        ", w_nbn_1()(2), w_nbn_4()(2), 1e-12);

    check("adjoint 51        ", w_nbn_1()(0), w_nbn_5()(0), 1e-12);
    check("adjoint 52        ", w_nbn_1()(1), w_nbn_5()(1), 1e-12);
    check("adjoint 53        ", w_nbn_1()(2), w_nbn_5()(2), 1e-12);

    check("adjoint 61        ", w_nbn_1()(0), w_nbn_6()(0), 1e-12);
    check("adjoint 62        ", w_nbn_1()(1), w_nbn_6()(1), 1e-12);
    check("adjoint 63        ", w_nbn_1()(2), w_nbn_6()(2), 1e-12);

    check("adjoint 71        ", w_nbn_1()(0), w_nbn_7()(0), 1e-12);
    check("adjoint 72        ", w_nbn_1()(1), w_nbn_7()(1), 1e-12);
    check("adjoint 73        ", w_nbn_1()(2), w_nbn_7()(2), 1e-12);

    check("adjoint 81        ", w_nbn_1()(0), w_nbn_8()(0), 1e-12);
    check("adjoint 82        ", w_nbn_1()(1), w_nbn_8()(1), 1e-12);
    check("adjoint 83        ", w_nbn_1()(2), w_nbn_8()(2), 1e-12);

    ang::so3_tangent w_nbb_1(R_nb / w_nbn_1());
    ang::so3_tangent w_nbb_2(q_nb / w_nbn_2());
    ang::so3_tangent w_nbb_3(r_nb / w_nbn_2());
    ang::so3_tangent w_nbb_4 = R_nb % w_nbn_4;
    ang::so3_tangent w_nbb_5 = q_nb % w_nbn_5;
    ang::so3_tangent w_nbb_6 = r_nb % w_nbn_6;
    ang::so3_tangent w_nbb_7 = ang::so3_tangent::wedge(R_nb % ang::so3_tangent::hat_skew(w_nbn_7));
    ang::so3_tangent w_nbb_8 = ang::so3_tangent::wedge(q_nb % ang::so3_tangent::hat_quat(w_nbn_8));

    check("adjoint inv 11    ", w_nbb_0()(0), w_nbb_1()(0), 1e-12);
    check("adjoint inv 12    ", w_nbb_0()(1), w_nbb_1()(1), 1e-12);
    check("adjoint inv 13    ", w_nbb_0()(2), w_nbb_1()(2), 1e-12);

    check("adjoint inv 21    ", w_nbb_0()(0), w_nbb_2()(0), 1e-12);
    check("adjoint inv 22    ", w_nbb_0()(1), w_nbb_2()(1), 1e-12);
    check("adjoint inv 23    ", w_nbb_0()(2), w_nbb_2()(2), 1e-12);

    check("adjoint inv 31    ", w_nbb_0()(0), w_nbb_3()(0), 1e-12);
    check("adjoint inv 32    ", w_nbb_0()(1), w_nbb_3()(1), 1e-12);
    check("adjoint inv 33    ", w_nbb_0()(2), w_nbb_3()(2), 1e-12);

    check("adjoint inv 41    ", w_nbb_0()(0), w_nbb_4()(0), 1e-12);
    check("adjoint inv 42    ", w_nbb_0()(1), w_nbb_4()(1), 1e-12);
    check("adjoint inv 43    ", w_nbb_0()(2), w_nbb_4()(2), 1e-12);

    check("adjoint inv 51    ", w_nbb_0()(0), w_nbb_5()(0), 1e-12);
    check("adjoint inv 52    ", w_nbb_0()(1), w_nbb_5()(1), 1e-12);
    check("adjoint inv 53    ", w_nbb_0()(2), w_nbb_5()(2), 1e-12);

    check("adjoint inv 61    ", w_nbb_0()(0), w_nbb_6()(0), 1e-12);
    check("adjoint inv 62    ", w_nbb_0()(1), w_nbb_6()(1), 1e-12);
    check("adjoint inv 63    ", w_nbb_0()(2), w_nbb_6()(2), 1e-12);

    check("adjoint inv 71    ", w_nbb_0()(0), w_nbb_7()(0), 1e-12);
    check("adjoint inv 72    ", w_nbb_0()(1), w_nbb_7()(1), 1e-12);
    check("adjoint inv 73    ", w_nbb_0()(2), w_nbb_7()(2), 1e-12);

    check("adjoint inv 81    ", w_nbb_0()(0), w_nbb_8()(0), 1e-12);
    check("adjoint inv 82    ", w_nbb_0()(1), w_nbb_8()(1), 1e-12);
    check("adjoint inv 83    ", w_nbb_0()(2), w_nbb_8()(2), 1e-12);

    Eigen::Matrix3d Adf1 = R_nb.adjoint_matrix_forward();
    Eigen::Matrix3d Adf2 = q_nb.adjoint_matrix_forward();
    Eigen::Matrix3d Adf3 = r_nb.adjoint_matrix_forward();
    Eigen::Vector3d w_nbn_9 = Adf1 * w_nbb_0();

    check("adjoint for 11    ", Adf1(0,0), Adf2(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf2(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf2(0,2), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf2(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf2(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf2(1,2), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf2(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf2(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf2(2,2), 1e-12);

    check("adjoint for 11    ", Adf1(0,0), Adf3(0,0), 1e-12);
    check("adjoint for 12    ", Adf1(0,1), Adf3(0,1), 1e-12);
    check("adjoint for 13    ", Adf1(0,2), Adf3(0,2), 1e-12);
    check("adjoint for 21    ", Adf1(1,0), Adf3(1,0), 1e-12);
    check("adjoint for 22    ", Adf1(1,1), Adf3(1,1), 1e-12);
    check("adjoint for 23    ", Adf1(1,2), Adf3(1,2), 1e-12);
    check("adjoint for 31    ", Adf1(2,0), Adf3(2,0), 1e-12);
    check("adjoint for 32    ", Adf1(2,1), Adf3(2,1), 1e-12);
    check("adjoint for 33    ", Adf1(2,2), Adf3(2,2), 1e-12);

    check("adjoint for 1     ", w_nbn_9(0), w_nbn_1()(0), 1e-12);
    check("adjoint for 2     ", w_nbn_9(1), w_nbn_1()(1), 1e-12);
    check("adjoint for 3     ", w_nbn_9(2), w_nbn_1()(2), 1e-12);

    Eigen::Matrix3d Adb1 = R_nb.adjoint_matrix_backward();
    Eigen::Matrix3d Adb2 = q_nb.adjoint_matrix_backward();
    Eigen::Matrix3d Adb3 = r_nb.adjoint_matrix_backward();
    Eigen::Vector3d w_nbb_9 = Adb1 * w_nbn_9;

    check("adjoint bac 11    ", Adb1(0,0), Adb2(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb2(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb2(0,2), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb2(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb2(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb2(1,2), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb2(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb2(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb2(2,2), 1e-12);

    check("adjoint bac 11    ", Adb1(0,0), Adb3(0,0), 1e-12);
    check("adjoint bac 12    ", Adb1(0,1), Adb3(0,1), 1e-12);
    check("adjoint bac 13    ", Adb1(0,2), Adb3(0,2), 1e-12);
    check("adjoint bac 21    ", Adb1(1,0), Adb3(1,0), 1e-12);
    check("adjoint bac 22    ", Adb1(1,1), Adb3(1,1), 1e-12);
    check("adjoint bac 23    ", Adb1(1,2), Adb3(1,2), 1e-12);
    check("adjoint bac 31    ", Adb1(2,0), Adb3(2,0), 1e-12);
    check("adjoint bac 32    ", Adb1(2,1), Adb3(2,1), 1e-12);
    check("adjoint bac 33    ", Adb1(2,2), Adb3(2,2), 1e-12);

    check("adjoint bac 1     ", w_nbb_9(0), w_nbb_0()(0), 1e-12);
    check("adjoint bac 2     ", w_nbb_9(1), w_nbb_0()(1), 1e-12);
    check("adjoint bac 3     ", w_nbb_9(2), w_nbb_0()(2), 1e-12);

} // closes test_adjoint

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Tso3::test_point_velocity() {
    ang::so3_tangent w_nbb(0.51, -0.38, -0.14);

    ang::euler euler_nb(0.12, 0.17, -0.03);
    ang::dcm R_nb(euler_nb);
    ang::so3_tangent w_nbn = R_nb | w_nbb;

    Eigen::Vector3d p_b(7.1, -11.8, -9.4);
    Eigen::Vector3d p_n = R_nb * p_b;

    Eigen::Vector3d v_bpb_0 = w_nbb().cross(p_b);
    Eigen::Vector3d v_bpb_1 = w_nbb.point_velocity(p_b);

    check("velocity 11    ", v_bpb_0(0), v_bpb_1(0), 1e-12);
    check("velocity 12    ", v_bpb_0(1), v_bpb_1(1), 1e-12);
    check("velocity 13    ", v_bpb_0(2), v_bpb_1(2), 1e-12);

    Eigen::Vector3d v_npn_0 = R_nb * v_bpb_0; // origin of n coincides with origin of b
    Eigen::Vector3d v_npn_1 = w_nbn.point_velocity(p_n);

    check("velocity 21    ", v_npn_0(0), v_npn_1(0), 1e-12);
    check("velocity 22    ", v_npn_0(1), v_npn_1(1), 1e-12);
    check("velocity 23    ", v_npn_0(2), v_npn_1(2), 1e-12);

} // closes test_velocity

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



















