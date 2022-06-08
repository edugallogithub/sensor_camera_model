#ifndef ATT_TEST_INTEGRATION
#define ATT_TEST_INTEGRATION

#include "ang_test.h"
#include "ang/quat.h"
#include "ang/rotate/so3_tangent.h"
#include "ang/rotate/rodrigues.h"
#include "ang/transform/speu_rodrigues.h"
#include <jail/unit_test.h>
#include <Eigen/Core>
#include <vector>

/*
 * This text executes the same integration based on body linear and angular accelerations
 * employing all possible combinations of SO(3) and SE(3) representations, and either
 * employing space magnitudes or body magnitudes.
 * 01 --> space and rodrigues
 * 02 --> space and dcm
 * 03 --> space and rotv
 * 11 --> body  and rodrigues
 * 12 --> body  and dcm
 * 13 --> body  and rotv
 * 21 --> space and speu_rodrigues
 * 22 --> space and speu_dcm
 * 23 --> space and homogeneous
 * 24 --> space and trfv
 * 25 --> space and dual
 * 26 --> space and screw --- NOT YET NOT YET
 * 31 --> body  and speu_rodrigues
 * 32 --> body  and speu_dcm
 * 33 --> body  and homogeneous
 * 34 --> body  and trfv
 * 35 --> body  and dual
 * 36 --> body  and screw --- NOT YET NOT YET
*/

namespace ang {
    namespace test {

        class Tintegration: public ::jail::unit_test {
        public:
            /**< constructor based on counter */
            explicit Tintegration(jail::counter&);
            /**< execute tests and write results on console */
            void run() override;

            /**< specific tests */
            static void test_space_rodrigues     (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::rodrigues& q0_wb, const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_dcm           (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::dcm& R0_wb,       const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_rotv          (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::rotv& rv0_wb,     const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_rodrigues      (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::rodrigues& q0_wb, const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_dcm            (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::dcm& R0_wb,       const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_rotv           (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const Eigen::Vector3d& x0_wbw_m, const ang::rotv& rv0_wb,     const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_speu_rodrigues(std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::speu_rodrigues& Gq0_wb,                            const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_speu_dcm      (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::speu_dcm& GR0_wb,                                  const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_homogeneous   (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::homogeneous& M0_wb,                                const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_trfv          (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::trfv& tau0_wb,                                     const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpw_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_space_dual          (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::dual& Z0_wb,                                       const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_speu_rodrigues (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::speu_rodrigues& Gq0_wb,                            const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_speu_dcm       (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::speu_dcm& GR0_wb,                                  const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_homogeneous    (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::homogeneous& M0_wb,                                const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_trfv           (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::trfv& tau0_wb,                                     const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_body_dual           (std::vector<double>& res, const double& Deltat_sec, const double& t_sec_end, const ang::dual& Z0_wb,                                       const Eigen::Vector3d& v0_b_mps, const ang::so3_tangent& w0_bwb_rps, const Eigen::Vector3d& x_bpb_m, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);

            static void test_differences         (const ang::speu_dcm& GR0_wb, const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);
            static void test_analysis            (const ang::speu_dcm& GR0_wb, const Eigen::Vector3d& v0_w_mps, const ang::so3_tangent& w0_bww_rps, const Eigen::Vector3d& a_b_mps2, const Eigen::Vector3d& alpha_b_rps2);

        };

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


    } // closes namespace test
} // closes namespace ang

#endif

