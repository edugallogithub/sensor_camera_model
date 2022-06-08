#include "euler.h"
#include "rodrigues.h"
#include "dcm.h"
#include "so3_tangent.h"
#include "../tools.h"

using namespace ang;

// CLASS EULER
// ===========
// ===========

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

euler::euler(const rodrigues& q) {
    *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
    *(_ypr_rad + 1) = asin (-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
    *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
}
/* constructor based on Rodrigues parameters */

 euler::euler(const dcm& R) {
    *_ypr_rad       = atan2(R()(1,0), R()(0,0));
    *(_ypr_rad + 1) = asin(-R()(2,0));
    *(_ypr_rad + 2) = atan2(R()(2,1), R()(2,2));
}
/* constructor based on direction cosine matrix */

 euler::euler(const rotv& rv) {
	// there is no direct method;
	rodrigues q(rv);
    *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
    *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
    *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
}
/* constructor based on rotation vector */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

euler::euler(rodrigues&& q) {
     *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
     *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
     *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
}
/* move constructor based on Rodrigues parameters */

euler::euler(dcm&& R) {
     *_ypr_rad       = atan2(R()(1,0), R()(0,0));
     *(_ypr_rad + 1) = asin(-R()(2,0));
     *(_ypr_rad + 2) = atan2(R()(2,1), R()(2,2));
}
/* move constructor based on direction cosine matrix */

 euler::euler(rotv&& rv) {
     // there is no direct method;
     rodrigues q(rv);
     *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
     *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
     *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
}
/* move constructor based on rotation vector */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

euler& euler::operator=(const rodrigues& q) {
    *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
    *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
    *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
    return *this;
}
/* assignment operator based on Rodrigues parameters */

euler& euler::operator=(const dcm& R) {
    *_ypr_rad       = atan2(R()(1,0), R()(0,0));
    *(_ypr_rad + 1) = asin(-R()(2,0));
    *(_ypr_rad + 2) = atan2(R()(2,1), R()(2,2));
    return *this;
}
/* assignment operator based on direction cosine matrix */

euler& euler::operator=(const rotv& rv) {
    // there is no direct method;
    rodrigues q(rv);
    *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
    *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
    *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
    return *this;
}
/* assignment operator based on rotation vector */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

 euler& euler::operator=(rodrigues&& q) {
     *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
     *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
     *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
     return *this;
 }
/* move assignment operator based on Rodrigues parameters */

 euler& euler::operator=(dcm&& R) {
     *_ypr_rad       = atan2(R()(1,0), R()(0,0));
     *(_ypr_rad + 1) = asin(-R()(2,0));
     *(_ypr_rad + 2) = atan2(R()(2,1), R()(2,2));
     return *this;
 }
/* move assignment operator based on direction cosine matrix */

 euler& euler::operator=(rotv&& rv) {
     // there is no direct method;
     rodrigues q(rv);
     *_ypr_rad       = atan2(+2 * (+q()(1) * q()(2) + q()(0) * q()(3)), 1 - 2 * (pow(q()(2),2) + pow(q()(3),2)));
     *(_ypr_rad + 1) = asin(-2 * (-q()(0) * q()(2) + q()(1) * q()(3)));
     *(_ypr_rad + 2) = atan2(+2 * (+q()(2) * q()(3) + q()(0) * q()(1)), 1 - 2 * (pow(q()(1),2) + pow(q()(2),2)));
     return *this;
 }
/* move assignment operator based on rotation vector */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d euler::operator*(const Eigen::Vector3d& vecin) const {
    double sy = sin(*_ypr_rad);
    double cy = cos(*_ypr_rad);
    double sp = sin(*(_ypr_rad + 1));
    double cp = cos(*(_ypr_rad + 1));
    double sr = sin(*(_ypr_rad + 2));
    double cr = cos(*(_ypr_rad + 2));

    double vec11 = vecin(0);
    double vec12 = cr * vecin(1) - sr * vecin(2);
    double vec13 = cr * vecin(2) + sr * vecin(1);

    double vec21 = cp * vec11 + sp * vec13;
    double vec22 = vec12;
    double vec23 = cp * vec13 - sp * vec11;

    return {cy * vec21 - sy * vec22, cy * vec22 + sy * vec21, vec23};
}
/* overloaded operator * (forward rotation) */

Eigen::Vector3d euler::operator/(const Eigen::Vector3d& vecin) const {
    double sy = sin(*_ypr_rad);
    double cy = cos(*_ypr_rad);
    double sp = sin(*(_ypr_rad + 1));
    double cp = cos(*(_ypr_rad + 1));
    double sr = sin(*(_ypr_rad + 2));
    double cr = cos(*(_ypr_rad + 2));

    double vec11 = cy * vecin(0) + sy * vecin(1);
    double vec12 = cy * vecin(1) - sy * vecin(0);
    double vec13 = vecin(2);

    double vec21 = cp * vec11 - sp * vec13;
    double vec22 = vec12;
    double vec23 = cp * vec13 + sp * vec11;

    return {vec21, cr * vec22 + sr * vec23, cr * vec23 - sr * vec22};
}
/* overloaded operator / (backward rotation) */

/* ===== ===== Angular Velocity - Time Derivative ===== ===== */
/* ========================================================== */

 so3_tangent euler::dot2omegabody(const Eigen::Vector3d& eulerdot_rps) const {
    double sp = sin(*(_ypr_rad + 1));
    double cp = cos(*(_ypr_rad + 1));
    double sr = sin(*(_ypr_rad + 2));
    double cr = cos(*(_ypr_rad + 2));
	return {-sp * eulerdot_rps(0) + eulerdot_rps(2),
			+sr * cp * eulerdot_rps(0) + cr * eulerdot_rps(1),
			+cr * cp * eulerdot_rps(0) - sr * eulerdot_rps(1)};
}
/* obtains the body angular velocity from the Euler angles and their time derivatives. */

Eigen::Vector3d euler::omegabody2dot(const so3_tangent& w_body_rps) const {
	double tanp = tan(*(_ypr_rad + 1));
	double secp = 1 / cos(*(_ypr_rad + 1));
    double sinr = sin(*(_ypr_rad + 2));
    double cosr = cos(*(_ypr_rad + 2));
	return {sinr * secp * w_body_rps()(1) + cosr * secp * w_body_rps()(2),
		cosr * w_body_rps()(1) - sinr * w_body_rps()(2),
		w_body_rps()(0) + sinr * tanp * w_body_rps()(1) + cosr * tanp * w_body_rps()(2)};
}
/* obtains the Euler angles differentials with time based on the Euler angles and the body angular velocity. */

 so3_tangent euler::dot2omegaspace(const Eigen::Vector3d& eulerdot_rps) const {
    double sy = sin(*_ypr_rad);
    double cy = cos(*_ypr_rad);
    double sp = sin(*(_ypr_rad + 1));
    double cp = cos(*(_ypr_rad + 1));
    return {- sy * eulerdot_rps(1) + cy * cp * eulerdot_rps(2),
               + cy * eulerdot_rps(1) + sy * cp * eulerdot_rps(2),
               eulerdot_rps(0) - sp * eulerdot_rps(2)};
}
/* obtains the space angular velocity from the Euler angles and their time derivatives. */

Eigen::Vector3d euler::omegaspace2dot(const so3_tangent& w_space_rps) const {
    double siny = sin(*_ypr_rad);
    double cosy = cos(*_ypr_rad);
    double tanp = tan(*(_ypr_rad + 1));
    double secp = 1 / cos(*(_ypr_rad + 1));
    return {cosy * tanp * w_space_rps()(0) + siny * tanp * w_space_rps()(1) + w_space_rps()(2),
               - siny * w_space_rps()(0)      + cosy * w_space_rps()(1),
               cosy * secp * w_space_rps()(0) + siny * secp * w_space_rps()(1)};
}
/* obtains the Euler angles differentials with time based on the Euler angles and the space angular velocity. */

/* ===== ===== Other ===== ===== */
/* ============================= */

void euler::obtain_euler_nedgrd(euler& euler_ng, const euler& euler_nb, const Eigen::Vector3d& v_n_mps) {
	// absolute bearing angle
	*euler_ng._ypr_rad = atan2(v_n_mps(1), v_n_mps(0));
	// absolute path angle
	*(euler_ng._ypr_rad + 1) = atan2(-v_n_mps(2), sqrt(pow(v_n_mps(0),2) + pow(v_n_mps(1),2)));

	// yaw[-180, 180] and pitch[-90, 90]
	double sgamma = sin(*(euler_ng._ypr_rad + 1));
	double cgamma = cos(*(euler_ng._ypr_rad + 1));

	double stheta = sin(*(euler_nb._ypr_rad + 1));
	double ctheta = cos(*(euler_nb._ypr_rad + 1));

	double sxi = sin(*(euler_nb._ypr_rad + 2));
	double cxi = cos(*(euler_nb._ypr_rad + 2));
	double txi = sxi / cxi;

	double Deltachi = *euler_ng._ypr_rad - *euler_nb._ypr_rad;
	tools::correct_yaw_rad(Deltachi);
	double sDeltachi = sin(Deltachi);
	double cDeltachi = cos(Deltachi);

	*(euler_ng._ypr_rad + 2) = atan((+sgamma * sDeltachi + txi * stheta * sgamma * cDeltachi + txi * ctheta * cgamma) / (+cDeltachi - txi * stheta * sDeltachi));
}
/* obtains the Euler angles that define the rotation between NED and GRD (chi, gamma, mu) from the Euler angles between NED and BFS (psi,
 * theta, xi) and the absolute speed in NED. */

double euler::control_bearing_diff(const double& a_deg, const double& b_deg) {
    double diff = a_deg - b_deg;
    double ab = std::fabs(diff);
    if (ab > 180.0) {
        diff = (360.0 - ab) * (-1) * diff / ab;
    }
    if (std::fabs(diff) > 90.0) {
        diff = std::nan("");
    }
    return diff;
}
/* returns difference between two bearings in degrees, nan if higher than plus minus 90. */

void euler::obtain_cross_long_track_errors(double& error_cross_m, double& error_long_m, const double& error_north_m, const double& error_east_m, const double& chi_rad) {
    double error_hor_m = std::sqrt(std::pow(error_north_m,2) + std::pow(error_east_m,2));
    double alpha_rad   = std::atan2(error_east_m, error_north_m);
    double beta_rad    = alpha_rad - chi_rad;
    error_cross_m      = error_hor_m * sin(beta_rad);
    error_long_m       = error_hor_m * cos(beta_rad);
 }
/* fills up the cross track error (positive to the right, negative to the left) and the long
 * track errors (positive to the front, negative to the back) based on the North error (positive
 * to the North, negative to the South), the East error (positive to the East, negative to the
 * West), and the trajectory bearing angle */

/* ===== ===== Angles for Forward Looking Vectors (choose yaw and pitch) ===== ===== */
/* ================================================================================= */

 double euler::obtain_yaw_forward(const Eigen::Vector3d& v_ned) {
     return atan2(v_ned(1), v_ned(0));
 }
/* given a forward looking vector (generally speed) expressed in NED, it returns the yaw angle [rad] of that vector with respect to the NED frame */

 double euler::obtain_pitch_forward(const Eigen::Vector3d& v_ned) {
     return atan2(- v_ned(2), sqrt(pow(v_ned(0), 2) + pow(v_ned(1), 2)));
 }
/* given a forward looking vector (generally speed) expressed in NED, it returns the pitch angle [rad] of that vector with respect to the NED frame */

/* ===== ===== Angles for Downward Looking Vectors (choose pitch and roll) ===== ===== */
/* =================================================================================== */

double euler::obtain_pitch_downward(const Eigen::Vector3d& g_ned) {
     return atan2(- g_ned(2), g_ned(0));
}
/*< given a downward looking vector (generally gravitation) expressed in NED, it returns the pitch angle [rad] of that vector with respect to the NED frame */

double euler::obtain_bank_downward(const Eigen::Vector3d& g_ned) {
     return atan2(- g_ned(1), sqrt(pow(g_ned(0), 2) + pow(g_ned(2), 2)));
}
/* given a downward looking vector (generally gravitation) expressed in NED, it returns the bank angle [rad] of that vector with respect to the NED frame */

/* ===== ===== Individual Rotations ===== ===== */
/* ============================================ */

ang::dcm euler::R3_yaw() const {
    double sy = sin(*_ypr_rad);
    double cy = cos(*_ypr_rad);

    ang::dcm R3;
    R3()(0,0) = + cy;
    R3()(0,1) = - sy;
    R3()(0,2) = 0.;
    R3()(1,0) = + sy;
    R3()(1,1) = + cy;
    R3()(1,2) = 0.;
    R3()(2,0) = 0.;
    R3()(2,1) = 0.;
    R3()(2,2) = 1.;

    return R3;
}
/* return rotation matrix for 1st rotation (yaw around 3rd axis) */

ang::dcm euler::R2_pitch() const {
    double sp = sin(*(_ypr_rad + 1));
    double cp = cos(*(_ypr_rad + 1));

    ang::dcm R2;
    R2()(0,0) = + cp;
    R2()(0,1) = 0.;
    R2()(0,2) = + sp;
    R2()(1,0) = 0.;
    R2()(1,1) = 1.;
    R2()(1,2) = 0.;
    R2()(2,0) = - sp;
    R2()(2,1) = 0.;
    R2()(2,2) = + cp;

    return R2;
}
/* return rotation matrix for 2nd rotation (pitch around 2nd axis) */

ang::dcm euler::R1_roll() const {
    double sr = sin(*(_ypr_rad + 2));
    double cr = cos(*(_ypr_rad + 2));

    ang::dcm R1;
    R1()(0,0) = 1.;
    R1()(0,1) = 0.;
    R1()(0,2) = 0.;
    R1()(1,0) = 0.;
    R1()(1,1) = + cr;
    R1()(1,2) = - sr;
    R1()(2,0) = 0.;
    R1()(2,1) = + sr;
    R1()(2,2) = + cr;

    return R1;
}
/* return rotation matrix for 3rd rotation (roll around 1st axis) */



















