#include "rotv.h"
#include "euler.h"
#include "rodrigues.h"
#include "dcm.h"
#include "so3_tangent.h"
#include "../auxiliary.h"
#include "../tools.h"

#include <cmath>
#include <iostream>

using namespace ang;

// CLASS ROTV
// ==========
// ==========

void rotv::shortest_angle() {
    while (this->norm() > math::constant::PI()) {
        this->get() = this->get() - this->get() / this->norm() * math::constant::PI() * 2.0;
    }
}
/* modify parameters maintaining rotation so it follows shortest path (rotation angle less than 180 [deg]) */

void rotv::equivalent_rotv(Eigen::Vector3d& rv) {
    double angle_rad = rv.norm();
    Eigen::Vector3d n = rv / angle_rad;
    rv = (angle_rad - 2.0 * math::constant::PI()) * n;
}
/* replaces the input 3x1 vector (which should come from a rotation vector by employing the () operator or
 * the get() method) by an equivalent vector of angle [2PI - angle] instead of [angle] and direction [-n]
 * instead of [n]. This is completely equivalent, but the rotation angle is now between 180 and 360 [deg],
 * and that is the reason why this method employs Euclidean 3x1 vectors instead of rotation vector objects. */

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

rotv::rotv(double p1, double p2, double p3)
: Eigen::Vector3d(p1, p2, p3) {
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* normal constructor based on rotation vector components */

rotv::rotv(const Eigen::Vector3d& Ovec)
: Eigen::Vector3d(Ovec) {
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* constructor based on size 3 vector */

rotv::rotv(const Eigen::Vector3d& Ovec, const double& factor)
: Eigen::Vector3d(Ovec * factor) {
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* constructor based on size 3 vector multiplied by input factor */

rotv::rotv(const euler& Oeuler) {
	double sr = std::sin(Oeuler.get_bank_rad());
	double cr = std::cos(Oeuler.get_bank_rad());
	double sp = std::sin(Oeuler.get_pitch_rad());
	double cp = std::cos(Oeuler.get_pitch_rad());
	double sy = std::sin(Oeuler.get_yaw_rad());
	double cy = std::cos(Oeuler.get_yaw_rad());

    double angle_rad = acos(0.5 * (cp * cy + cr * cy + sr * sp * sy + cr * cp - 1)); // [0, PI]

    if (angle_rad > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (sr * cp + sr * cy - cr * sp * sy);
        (*this)()(1) = f * (sr * sy + cr * sp * cy + sp);
        (*this)()(2) = f * (cp * sy + cr * sy - sr * sp * cy);
    }
    else {
        // convert from euler angles to quaternion and then to rotation vector
        rodrigues temp(Oeuler);
        (*this) = temp;
    }
}
/* constructor based on Euler angles */

rotv::rotv(const rodrigues& q) {
    double q_vec_norm = q().segment<3>(1).norm();
    double angle_rad = 2.0 * std::atan2(q_vec_norm, q()(0)); // 2 * [-PI, PI] = [-2PI, 2PI]
    if (angle_rad > math::constant::PI()) {
        angle_rad = angle_rad - 2.0 * math::constant::PI();
    }
    else if (angle_rad < - math::constant::PI()) {
        angle_rad = angle_rad + 2.0 * math::constant::PI();
    }
    // now angle [-PI,PI]
    if (fabs(angle_rad) > math::constant::SMALL_ROT()) {
        (*this) = q().segment<3>(1) * angle_rad / q_vec_norm;
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        (*this) = q().segment<3>(1) * 2.0 / q()(0) * (1.0 - q_vec_norm / (3.0 * q()(0) * q()(0)));
    }
}
/* constructor based on Rodrigues parameters */

rotv::rotv(const dcm& R) {
	double angle_rad = acos(0.5 * (R().trace() - 1)); // [0, PI]
    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (R()(2, 1) - R()(1, 2));
        (*this)()(1) = f * (R()(0, 2) - R()(2, 0));
        (*this)()(2) = f * (R()(1, 0) - R()(0, 1));
    }
    else {
        // convert from rotation matrix to quaternion and then to rotation vector
        rodrigues temp(R);
        (*this) = temp;
    }
}
/* constructor based on direction cosine matrix */

rotv::rotv(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    *this = v1.cross(v2).normalized() * acos(v1.dot(v2) / v1.norm() / v2.norm());
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* VERY IMPORTANT constructor that obtains the rotation from v1 to v2, both normalized.
The direction is orthogonal to the plane formed by v1 and v2. The magnitude is the angle
between both input vectors in that plane. The modulus of the input vectors does not matter.
The result is such that v2.normalized() = this * v1.normalized(). */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

rotv::rotv(Eigen::Vector3d&& Ovec)
: Eigen::Vector3d(Ovec) {
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* move constructor based on size 3 vector */

rotv::rotv(const Eigen::Vector3d&& Ovec, const double& factor)
: Eigen::Vector3d(Ovec * factor) {
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* move constructor based on size 3 vector multiplied by input factor */

rotv::rotv(euler&& Oeuler) {
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());

    double angle_rad = acos(0.5 * (cp * cy + cr * cy + sr * sp * sy + cr * cp - 1)); // [0, PI]

    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (sr * cp + sr * cy - cr * sp * sy);
        (*this)()(1) = f * (sr * sy + cr * sp * cy + sp);
        (*this)()(2) = f * (cp * sy + cr * sy - sr * sp * cy);
    }
    else {
        // convert from euler angles to quaternion and then to rotation vector
        rodrigues temp(Oeuler);
        (*this) = temp;
    }
}
/* move constructor based on Euler angles */

rotv::rotv(rodrigues&& q) {
    double q_vec_norm = q().segment<3>(1).norm();
    double angle_rad = 2.0 * std::atan2(q_vec_norm, q()(0)); // 2 * [-PI, PI] = [-2PI, 2PI]
    if (angle_rad > math::constant::PI()) {
        angle_rad = angle_rad - 2.0 * math::constant::PI();
    }
    else if (angle_rad < - math::constant::PI()) {
        angle_rad = angle_rad + 2.0 * math::constant::PI();
    }
    // now angle [-PI,PI]
    if (angle_rad > math::constant::SMALL_ROT()) {
        (*this) = q().segment<3>(1) * angle_rad / q_vec_norm;
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        (*this) = q().segment<3>(1) * 2.0 / q()(0) * (1.0 - q_vec_norm / (3.0 * q()(0) * q()(0)));
    }
}
/* move constructor based on Rodrigues parameters */

rotv::rotv(dcm&& R) {
    double angle_rad = acos(0.5 * (R().trace() - 1)); // [0, PI]
    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (R()(2, 1) - R()(1, 2));
        (*this)()(1) = f * (R()(0, 2) - R()(2, 0));
        (*this)()(2) = f * (R()(1, 0) - R()(0, 1));
    }
    else {
        // convert from rotation matrix to quaternion and then to rotation vector
        rodrigues temp(R);
        (*this) = temp;
    }
}
/* move constructor based on direction cosine matrix */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

rotv& rotv::operator=(const Eigen::Vector3d& op2) {
    this->get() = op2;
    this->shortest_angle(); // ensure that rotation is less than PI
    return *this;
}
/* assignment operator based on size 3 vector */

rotv& rotv::operator=(const euler& Oeuler) {
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());

    double angle_rad = acos(0.5 * (cp * cy + cr * cy + sr * sp * sy + cr * cp - 1)); // [0, PI]

    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (sr * cp + sr * cy - cr * sp * sy);
        (*this)()(1) = f * (sr * sy + cr * sp * cy + sp);
        (*this)()(2) = f * (cp * sy + cr * sy - sr * sp * cy);
    }
    else {
        // convert from euler angles to quaternion and then to rotation vector
        rodrigues temp(Oeuler);
        (*this) = temp;
    }
    return *this;
}
/* assignment operator based on Euler angles */

rotv& rotv::operator=(const rodrigues& q) {
    double q_vec_norm = q().segment<3>(1).norm();
    double angle_rad = 2.0 * std::atan2(q_vec_norm, q()(0)); // 2 * [-PI, PI] = [-2PI, 2PI]
    if (angle_rad > math::constant::PI()) {
        angle_rad = angle_rad - 2.0 * math::constant::PI();
    }
    else if (angle_rad < - math::constant::PI()) {
        angle_rad = angle_rad + 2.0 * math::constant::PI();
    }
    // now angle [-PI,PI]
    if (angle_rad > math::constant::SMALL_ROT()) {
        (*this) = q().segment<3>(1) * angle_rad / q_vec_norm;
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        (*this) = q().segment<3>(1) * 2.0 / q()(0) * (1.0 - q_vec_norm / (3.0 * q()(0) * q()(0)));
    }
    return *this;
}
/* assignment operator based on Rodrigues parameters */

rotv& rotv::operator=(const dcm& R) {
    double angle_rad = acos(0.5 * (R().trace() - 1)); // [0, PI]
    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (R()(2, 1) - R()(1, 2));
        (*this)()(1) = f * (R()(0, 2) - R()(2, 0));
        (*this)()(2) = f * (R()(1, 0) - R()(0, 1));
    }
    else {
        // convert from rotation matrix to quaternion and then to rotation vector
        rodrigues temp(R);
        (*this) = temp;
    }
    return *this;
}
/* assignment operator based on direction cosine matrix */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

rotv& rotv::operator=(Eigen::Vector3d&& op2) {
    this->get() = op2;
    this->shortest_angle(); // ensure that rotation is less than PI
}
/* move assignment operator based on size 3 vector */

rotv& rotv::operator=(euler&& Oeuler) {
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());

    double angle_rad = acos(0.5 * (cp * cy + cr * cy + sr * sp * sy + cr * cp - 1)); // [0, PI]

    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (sr * cp + sr * cy - cr * sp * sy);
        (*this)()(1) = f * (sr * sy + cr * sp * cy + sp);
        (*this)()(2) = f * (cp * sy + cr * sy - sr * sp * cy);
    }
    else {
        // convert from euler angles to quaternion and then to rotation vector
        rodrigues temp(Oeuler);
        (*this) = temp;
    }
    return *this;
}
/* move assignment operator based on Euler angles */

rotv& rotv::operator=(rodrigues&& q) {
    double q_vec_norm = q().segment<3>(1).norm();
    double angle_rad = 2.0 * std::atan2(q_vec_norm, q()(0)); // 2 * [-PI, PI] = [-2PI, 2PI]
    if (angle_rad > math::constant::PI()) {
        angle_rad = angle_rad - 2.0 * math::constant::PI();
    }
    else if (angle_rad < - math::constant::PI()) {
        angle_rad = angle_rad + 2.0 * math::constant::PI();
    }
    // now angle [-PI,PI]
    if (angle_rad > math::constant::SMALL_ROT()) {
        (*this) = q().segment<3>(1) * angle_rad / q_vec_norm;
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        (*this) = q().segment<3>(1) * 2.0 / q()(0) * (1.0 - q_vec_norm / (3.0 * q()(0) * q()(0)));
    }
    return *this;
}
/* move assignment operator based on Rodrigues parameters */

rotv& rotv::operator=(dcm&& R) {
    double angle_rad = acos(0.5 * (R().trace() - 1)); // [0, PI]
    if (fabs(angle_rad) > math::constant::EPS()) {
        double f = 0.5 * angle_rad / std::sin(angle_rad);
        (*this)()(0) = f * (R()(2, 1) - R()(1, 2));
        (*this)()(1) = f * (R()(0, 2) - R()(2, 0));
        (*this)()(2) = f * (R()(1, 0) - R()(0, 1));
    }
    else {
        // convert from rotation matrix to quaternion and then to rotation vector
        rodrigues temp(R);
        (*this) = temp;
    }
    return *this;
}
/* move assignment operator based on direction cosine matrix */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

rotv rotv::operator*(const rotv& op2) const {
	double angle_rad = this->norm();
    double op2_angle_rad = op2.norm();

    if (angle_rad < math::constant::EPS()) {
        return op2;
    }
    else if (op2_angle_rad < math::constant::EPS()) {
        return *this;
    }
    else {
        double cos1 = std::cos(0.5 * angle_rad);
        double sin1 = std::sin(0.5 * angle_rad);
        double cos2 = std::cos(0.5 * op2_angle_rad);
        double sin2 = std::sin(0.5 * op2_angle_rad);

        double res_angle_rad = 2 * acos(cos1 * cos2 - sin1 * sin2 / angle_rad / op2_angle_rad * (op2.dot(*this)));
        double sin3 = std::sin(0.5 * res_angle_rad);

        rotv rv((this->get() / angle_rad * sin1 * cos2 / sin3
                          + op2() / op2_angle_rad * sin2 * cos1 / sin3
                          - op2.cross(*this) * sin1 * sin2 / sin3 / angle_rad / op2_angle_rad) * res_angle_rad);

        double fac = sin1 * cos2 / sin3 / angle_rad;
        Eigen::Vector3d PPP1 = this->get() * fac;
        Eigen::Vector3d PPP2 = op2() / op2_angle_rad * sin2 * cos1 / sin3;
        Eigen::Vector3d PPP3 = - op2.cross(*this) * sin1 * sin2 / sin3 / angle_rad / op2_angle_rad;

        rv.shortest_angle(); // ensure that rotation is less than PI
        return rv;
    }
}
/* overloaded operator * (combination of rotations) */

rotv rotv::operator/(const rotv& op2) const {
    double angle_rad = this->norm();
    double op2_angle_rad = op2.norm();

    if (angle_rad < math::constant::EPS()) {
        return op2;
    }
    else if (op2_angle_rad < math::constant::EPS()) {
        return this->inverse();
    }
    else {
        double cos1 = std::cos(0.5 * angle_rad);
        double sin1 = std::sin(0.5 * angle_rad);
        double cos2 = std::cos(0.5 * op2_angle_rad);
        double sin2 = std::sin(0.5 * op2_angle_rad);

        double res_angle_rad = 2 * acos(cos1 * cos2 + sin1 * sin2 / angle_rad / op2_angle_rad * (op2.dot(*this)));
        double sin3 = std::sin(0.5 * res_angle_rad);

        rotv rv((- this->get() / angle_rad * sin1 * cos2 / sin3
                          + op2() / op2_angle_rad * sin2 * cos1 / sin3
                          + op2.cross(*this) * sin1 * sin2 / sin3 / angle_rad / op2_angle_rad) * res_angle_rad);

        rv.shortest_angle(); // ensure that rotation is less than PI
        return rv;
    }
}
/* overloaded operator / (backward combination of rotations) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d rotv::operator*(const Eigen::Vector3d& vecin) const {
    double angle_rad    = this->norm();
    if (angle_rad < math::constant::EPS()) {return vecin;}
	double angle_rad_sq = angle_rad * angle_rad;
    return tools::itself_by_transpose(*this) * vecin / angle_rad_sq - this->cross(this->cross(vecin)) * std::cos(angle_rad) / angle_rad_sq + this->cross(vecin) * std::sin(angle_rad) / angle_rad;
}
/* overloaded operator * (forward rotation) */

Eigen::Vector3d rotv::operator/(const Eigen::Vector3d& vecin) const {
    double angle_rad    = this->norm();
    if (angle_rad < math::constant::EPS()) {return vecin;}
    double angle_rad_sq = angle_rad * angle_rad;
    return tools::itself_by_transpose(*this) * vecin / angle_rad_sq - this->cross(this->cross(vecin)) * std::cos(angle_rad) / angle_rad_sq - this->cross(vecin) * std::sin(angle_rad) / angle_rad;
}
/* overloaded operator / (backward rotation) */

rodrigues rotv::exp_map_rodrigues() const {
    return rodrigues(*this);
}
/* exponential map that returns rodrigues parameters */

dcm rotv::exp_map_dcm() const {
    return dcm(*this);
}
/* exponential map that returns direction cosine matrix */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

so3_tangent rotv::operator|(const so3_tangent& w) const {
    return so3_tangent((*this) * w());
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix3d rotv::adjoint_matrix_forward() const {
    return dcm(*this)();
}
/* returns forward adjoint matrix */

so3_tangent rotv::operator%(const so3_tangent& w) const {
    return so3_tangent((*this) / w());
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix3d rotv::adjoint_matrix_backward() const {
    return dcm(this->inverse())();
}
/* returns backward adjoint matrix */

/* ===== ===== Angular Velocity - Time Derivative ===== ===== */
/* ========================================================== */

so3_tangent rotv::dot2omegabody(const Eigen::Vector3d& rvdot) const {
    double angle_rad = this->norm();
    if (angle_rad < math::constant::EPS()) {return so3_tangent(rvdot);}
    return so3_tangent(rvdot - this->cross(rvdot) * (1 - std::cos(angle_rad)) / std::pow(angle_rad,2) + this->cross(this->cross(rvdot)) * (angle_rad - std::sin(angle_rad)) / std::pow(angle_rad,3));
}
/* obtains the body angular velocity from the rotation vector and its time derivative. */

Eigen::Vector3d rotv::omegabody2dot(const so3_tangent& w_body_rps) const {
    double angle_rad = this->norm();
    if (angle_rad < math::constant::EPS()) {return w_body_rps();}
    return w_body_rps() + this->cross(w_body_rps()) * 0.5 + this->cross(this->cross(w_body_rps())) * (1 / std::pow(angle_rad,2) - 0.5 / angle_rad / tan(0.5 * angle_rad));
}
/* obtains the rotation vector differential with time based on the rotation vector and the body angular velocity. */

so3_tangent rotv::dot2omegaspace(const Eigen::Vector3d& rvdot) const {
    double angle_rad = this->norm();
    if (angle_rad < math::constant::EPS()) {return so3_tangent(rvdot);}
    return so3_tangent(rvdot + this->cross(rvdot) * (1 - std::cos(angle_rad)) / std::pow(angle_rad,2) + this->cross(this->cross(rvdot)) * (angle_rad - std::sin(angle_rad)) / std::pow(angle_rad,3));
}
/* obtains the space angular velocity from the rotation vector and its time derivative. */

Eigen::Vector3d rotv::omegaspace2dot(const so3_tangent& w_space_rps) const {
    double angle_rad = this->norm();
    if (angle_rad < math::constant::EPS()) {return w_space_rps();}
    return w_space_rps() - this->cross(w_space_rps()) * 0.5 + this->cross(this->cross(w_space_rps())) * (1 / std::pow(angle_rad,2) - 0.5 / angle_rad / tan(0.5 * angle_rad));
}
/* obtains the rotation vector differential with time based on the rotation vector and the space angular velocity. */

/* ===== ===== Linear Algebra ===== ===== */
/* ====================================== */

Eigen::Vector3d rotv::euclidean_diff(const rotv& rv1, Eigen::Vector3d& r2) {
    Eigen::Vector3d res = rv1() - r2;
    double res_norm = res.norm();
    if (res_norm > math::constant::PI()) {
        rotv::equivalent_rotv(r2);
        res = rv1() - r2;
    }
    return res;
}
/* This function returns a very good approximation (specially when the two input rotations are close) to (rv1 / rv2),
 * computed as (rv1() - rv2()). It also modifies one of the inputs, so it can be used repeatedly with the same
 * rv2 without going into the conditional. *
 * This function shall be used exclusively when the intended difference between both inputs is relatively small.
 * The reason is NOT that this approximation is better the smaller their difference is (check the test_rotv_eclidean_diff
 * method), but that to compute the initial difference on which to apply the "if" condition (not the final result,
 * which in this case is the same), both inputs should be rotv (here one is Vector3d), and it would be (rv1 / rv2).norm().
 * But then we can not change the second one as the angle has to be less than 180 [deg], and the equivalent_rotv function
 * changes it to between 180 and 360 [deg]. */

/* ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================= */
Eigen::Matrix3d rotv::jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - dcm(*this)() * tools::skew3(vec);
}
/* returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(rv * v) / dDeltarB,
 * (rv plus DeltarB) * v = rv * v + J * DeltarB */

Eigen::Matrix3d rotv::jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    // return dcm(*this).inverse()() * tools::skew3(vec) * dcm(*this)(); // this result is identical but more expensive to compute
    return tools::skew3(dcm(*this).inverse()() * vec);
}
/* returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(rv / v) / dDeltarB,
 * (rv plus DeltarB) / v = rv / v + J * DeltarB */

Eigen::Matrix3d rotv::jac_right_log() const {
    return this->jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(rv)) / dDeltarv, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(rv plus Deltarv) ~= Log(rv) + (J * Deltarv) */

Eigen::Matrix3d rotv::jac_right_plus_wrt_first(const rotv& rv) const {
    return rv.adjoint_matrix_backward();
}
/* returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
 * equal to d(rv1 plus rv2) / dDeltar1B.
 * (rv1 plus Deltar1B) plus rv2 = (rv1 plus rv2) plus J * Deltar1B */

Eigen::Matrix3d rotv::jac_right_plus_wrt_second(const rotv& rv) const {
    return rv.jac_right();
}
/* returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
 * equal to d(rv1 plus rv2) / dDeltar2B.
 * rv1 plus (rv2 + Deltar2B) = (rv1 plus rv2) plus J * Deltar2B */

Eigen::Matrix3d rotv::jac_right_minus_wrt_first(const rotv& r) const {
    return this->minus_right(r).jac_right_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the first element,
 * equal to d(rv2 minus rv1) / dDeltar2B.
 * (rv2 plus Deltar2B) minus rv1 = (rv2 minus rv1) + J * Deltar2B */

Eigen::Matrix3d rotv::jac_right_minus_wrt_second(const rotv& r) const {
    return - this->minus_right(r).jac_left_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the second element,
 * equal to d(rv2 minus rv1) / dDeltar1B.
 * rv2 minus(rv1 plus Deltar1B) = (rv2 minus rv1) + J * Deltar1B */

Eigen::Matrix3d rotv::jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - dcm(*this)() * tools::skew3(w());
}
/* returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(Adr | w) / dDeltarB,
 * Ad(r plus DeltarB) | w = Adr | w + J * DeltarB */

Eigen::Matrix3d rotv::jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return tools::skew3(dcm(*this).inverse()() * w());
}
/* returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(Adr % w) / dDeltarB,
 * Ad(r plus DeltarB) % w = Adr % w + J * DeltarB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix3d rotv::jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - tools::skew3(dcm(*this)() * vec);
}
/* returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(rv * v) / dDeltarN,
 * (DeltarN plus rv) * v = rv * v + J * DeltarN */

Eigen::Matrix3d rotv::jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return dcm(*this)().transpose() * tools::skew3(vec);
}
/* returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(rv / v) / dDeltarN,
 * (DeltarN plus rv) / v = rv / v + J * DeltarN */

Eigen::Matrix3d rotv::jac_left_log() const {
    return this->jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(rv)) / dDeltarv, which coincides with the inverse
 * of the left jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(Deltarv plus rv) ~= Log(rv) + (J * Deltarv) */

Eigen::Matrix3d rotv::jac_left_plus_wrt_first(const rotv& rv) const {
    return rv.jac_left();
}
/* returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
 * equal to d(rv1 plus rv2) / dDeltar1N.
 * (rv1 + Deltar1N) plus rv2 = J * Deltar1N plus (rv1 plus rv2)  */

Eigen::Matrix3d rotv::jac_left_plus_wrt_second(const rotv& rv) const {
    return rv.adjoint_matrix_forward();
}
/* returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
 * equal to d(rv1 plus rv2) / dDeltar2N.
 * rv1 plus (Deltar2N plus rv2) = J * Deltar2N plus (rv1 plus rv2) */

Eigen::Matrix3d rotv::jac_left_minus_wrt_first(const rotv& rv) const {
    return this->minus_left(rv).jac_left_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the first element,
 * equal to d(rv2 minus rv1) / dDeltar2N.
 * (Deltar2N plus rv2) minus rv1 = (rv2 minus rv1) + J * Deltar2N */

Eigen::Matrix3d rotv::jac_left_minus_wrt_second(const rotv& rv) const {
    return - this->minus_left(rv).jac_right_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the second element,
 * equal to d(rv2 minus rv1) / dDeltar1N.
 * rv2 minus (Deltar1N plus rv1) = (rv2 minus rv1) + J * Deltar1N */

Eigen::Matrix3d rotv::jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - tools::skew3(dcm(*this)() * w());
}
/* returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(Adrv | w) / dDeltarN,
 * Ad(DeltarN plus rv) | w = Adrv | w + J * DeltarN */

Eigen::Matrix3d rotv::jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return dcm(*this)().transpose() * tools::skew3(w());
}
/* returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(Adrv % w) / dDeltarN,
 * Ad(DeltarN plus rv) % w = Adrv % w + J * DeltarN */

/* ===== ===== Exponential (Right and Left) Jacobians ===== ===== */
/* ============================================================== */

Eigen::Matrix3d rotv::jac_right() const {
    double rv_norm = this->norm();
    Eigen::Matrix3d rv_skew = tools::skew3(*this);
    if (rv_norm > math::constant::SMALL_ROT()) {
        return Eigen::Matrix3d::Identity() - rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2) + rv_skew * rv_skew * (rv_norm - std::sin(rv_norm)) / std::pow(rv_norm,3);
    }
    else { // easy Taylor expansion of above expression
        return Eigen::Matrix3d::Identity() - rv_skew * 0.5 + rv_skew * rv_skew / 6.0;
    }
}
/* returns the jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector itself,
 * also known as the right Jacobian of SO(3).
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(rv + Deltarv) ~= exp(rv) plus (JR * Deltarv) */

Eigen::Matrix3d rotv::jac_right_inv() const {
    double rv_norm = this->norm();
    Eigen::Matrix3d rv_skew = tools::skew3(*this);
    if (rv_norm > math::constant::SMALL_ROT()) {
        return Eigen::Matrix3d::Identity() + rv_skew * 0.5 + rv_skew * rv_skew * (1 / rv_norm / rv_norm - (1 + std::cos(rv_norm)) / (2 * rv_norm * std::sin(rv_norm)));
    }
    else { // perform a Taylor expansion of above // TODO Not done yet
        return Eigen::Matrix3d::Identity() + rv_skew * 0.5 + rv_skew * rv_skew * (1 / rv_norm / rv_norm - (1 + std::cos(rv_norm)) / (2 * rv_norm * std::sin(rv_norm)));
    }
}
/* returns the inverse jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector
 * itself, also known as the inverse of the right Jacobian of SO(3).
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * rv + JRinv * Deltarv ~= log[exp(rv) plus Deltarv] */

Eigen::Matrix3d rotv::jac_left() const {
    double rv_norm = this->norm();
    Eigen::Matrix3d rv_skew = tools::skew3(*this);
    if (rv_norm > math::constant::SMALL_ROT()) {
        return Eigen::Matrix3d::Identity() + rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2) + rv_skew * rv_skew * (rv_norm - std::sin(rv_norm)) / std::pow(rv_norm,3);
    }
    else { // easy Taylor expansion of above expression
        return Eigen::Matrix3d::Identity() + rv_skew * 0.5 + rv_skew * rv_skew / 6.0;
    }
}
/* returns the jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector itself,
 * also known as the left Jacobian of SO(3).
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(rv + Deltarv) ~= (JL * Deltarv) plus_left exp(rv) */

Eigen::Matrix3d rotv::jac_left_inv() const {
    double rv_norm = this->norm();
    Eigen::Matrix3d rv_skew = tools::skew3(*this);
    if (rv_norm > math::constant::SMALL_ROT()) {
        return Eigen::Matrix3d::Identity() - rv_skew * 0.5 + rv_skew * rv_skew * (1 / rv_norm / rv_norm - (1 + std::cos(rv_norm)) / (2 * rv_norm * std::sin(rv_norm)));
    }
    else { // perform a Taylor expansion of above // TODO Not done yet
        return Eigen::Matrix3d::Identity() - rv_skew * 0.5 + rv_skew * rv_skew * (1 / rv_norm / rv_norm - (1 + std::cos(rv_norm)) / (2 * rv_norm * std::sin(rv_norm)));
    }
}
/* returns the inverse jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector
 * itself, also known as the inverse of the left Jacobian of SO(3).
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * rv + jacLinv(rv) * Deltarv ~= log[exp(Deltarv) * exp(rv)] */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d rotv::jac_euclidean_forward_rotation_wrt_vector() const {
    return dcm(*this)();
}
/* returns the jacobian of the forward rotation action with respect to the vector, equal to d(rv * v) / dDeltav,
 * rv * v = J * v
 * rv * (v + Deltav) = rv * v + J * Deltav */

Eigen::Matrix3d rotv::jac_euclidean_backward_rotation_wrt_vector() const {
    return dcm(this->inverse())();
}
/* returns the jacobian of the backward rotation action with respect to the vector, equal to d(rv / v) / dDeltav,
 * rv / v = J * v
 * rv / (v + Deltav) = rv / v + J * Deltav */

Eigen::Matrix3d rotv::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return dcm(*this)();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adr | w) / dDeltaw,
 * Adr | w = J * w
 * Adr | (w + Deltaw) = Adr | w + J * Deltaw */

Eigen::Matrix3d rotv::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adr % w) / dDeltaw,
 * Adr % w = J * w
 * Adr % (w + Deltaw) = Adr % w + J * Deltaw */

Eigen::Matrix3d rotv::jac_euclidean_forward_rotation_wrt_rotv(const Eigen::Vector3d& v) const {
    return - this->exp_map_dcm()() * tools::skew3(v) * this->jac_right();
}
/* returns the jacobian [3x3] of a forward rotation with respect to the rotation vector, equal to d(r * v)/dr.
 * The forward rotation is NOT linear on the rotation vector, so
 * r * v != d(r * v)/dr |r*v * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) * v ~= r * v + d(r * v)/dr |r*v * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the right jacobian, but the result
 * with the left jacobian is the same. */

Eigen::Matrix3d rotv::jac_euclidean_forward_rotation_wrt_rotv_bis(const Eigen::Vector3d& v) const {
    return - tools::skew3(this->exp_map_dcm() * v) * this->jac_left();
}
/* returns the jacobian [3x3] of a forward rotation with respect to the rotation vector, equal to d(r * v)/dr.
 * The forward rotation is NOT linear on the rotation vector, so
 * r * v != d(r * v)/dr |r*v * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) * v ~= r * v + d(r * v)/dr |r*v * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the left jacobian, but the result
 * with the right jacobian is the same. */

Eigen::Matrix3d rotv::jac_euclidean_backward_rotation_wrt_rotv(const Eigen::Vector3d& v) const {
    //return tools::skew3(this->exp_map_dcm() / v) * this->jac_right(); // Jan 2022 -> this is better but I do not use it to avoid introducing changes (albeit ridiculously small, rounding error)
    return - this->inverse().jac_euclidean_forward_rotation_wrt_rotv(v);
}
/* returns the jacobian [3x3] of a backward rotation with respect to the rotation vector, equal to d[r / v]/dr.
 * The backward rotation is NOT linear on the rotation vector, so
 * r / v != d(r / v)/dr |r/v * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) / v ~= r / v + d(r / v)/dr |r/v * Delta r
 * Note that the increment is Delta r. */

Eigen::Matrix3d rotv::jac_euclidean_backward_rotation_wrt_rotv_bis(const Eigen::Vector3d& v) const {
    //return - this->exp_map_dcm()().transpose() * tools::skew3(v) * this->jac_left();// Jan 2022 -> this is better but I do not use it to avoid introducing changes (albeit ridiculously small, rounding error)
    return - this->inverse().jac_euclidean_forward_rotation_wrt_rotv_bis(v);
}
/* returns the jacobian [3x3] of a backward rotation with respect to the rotation vector, equal to d(r / v)/dr.
 * The backward rotation is NOT linear on the rotation vector, so
 * r / v != d(r / v)/dr |r/v * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) / v ~= r / v + d(r / v)/dr |r/v * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the left jacobian, but the result
 * with the right jacobian is the same. */

Eigen::Matrix3d rotv::jac_euclidean_forward_adjoint_wrt_rotv(const so3_tangent& w) const {
    return this->jac_euclidean_forward_rotation_wrt_rotv(w());
}
/* returns the jacobian [3x3] of a forward adjoint with respect to the rotation vector, equal to d(r | w)/dr.
 * The forward adjoint is NOT linear on the rotation vector, so
 * r | w != d(r | w)/dr |r|w * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) | w ~= r | w + d(r | w)/dr |r|w * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the right jacobian, but the result
 * with the left jacobian is the same. */

Eigen::Matrix3d rotv::jac_euclidean_forward_adjoint_wrt_rotv_bis(const so3_tangent& w) const {
    return this->jac_euclidean_forward_rotation_wrt_rotv_bis(w());
}
/* returns the jacobian [3x3] of a forward adjoint with respect to the rotation vector, equal to d(r | w)/dr.
 * The forward adjoint is NOT linear on the rotation vector, so
 * r | w != d(r | w)/dr |r|w * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) | w ~= r | w + d(r | w)/dr |r|w * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the left jacobian, but the result
 * with the right jacobian is the same. */

Eigen::Matrix3d rotv::jac_euclidean_backward_adjoint_wrt_rotv(const so3_tangent& w) const {
    return this->jac_euclidean_backward_rotation_wrt_rotv(w());
}
/* returns the jacobian [3x3] of a backward adjoint with respect to the rotation vector, equal to d[r % w]/dr.
 * The backward adjoint is NOT linear on the rotation vector, so
 * r % w != d(r % w)/dr |r%w * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) % w ~= r % w + d(r % w)/dr |r%w * Delta r
 * Note that the increment is Delta r. */

Eigen::Matrix3d rotv::jac_euclidean_backward_adjoint_wrt_rotv_bis(const so3_tangent& w) const {
    return this->jac_euclidean_backward_rotation_wrt_rotv_bis(w());
}
/* returns the jacobian [3x3] of a backward adjoint with respect to the rotation vector, equal to d(r % w)/dr.
 * The backward adjoint is NOT linear on the rotation vector, so
 * r % w != d(r % w)/dr |r%w * r
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * (r + Delta r) % w ~= r % w + d(r % w)/dr |r%w * Delta r
 * Note that the increment is Delta r.
 * Note that internally the expression is based on the left jacobian, but the result
 * with the right jacobian is the same. */
































