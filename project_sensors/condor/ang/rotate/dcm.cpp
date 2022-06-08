#include "dcm.h"
#include "euler.h"
#include "rodrigues.h"
#include "rotv.h"
#include "so3_tangent.h"
#include "../auxiliary.h"
#include "../tools.h"

#include <cmath>

using namespace ang;

// CLASS DCM
// =========
// =========

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

dcm::dcm(const euler& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    (*this) << + cp * cy, - cr*sy + sr*sp*cy, + sr*sy + cr*sp*cy,
               + cp * sy, + cr*cy + sr*sp*sy, - sr*cy + cr*sp*sy,
               - sp,      + sr * cp,          + cr * cp;
}
/* constructor based on Euler angles */

dcm::dcm(const rodrigues& q) {
    double q02 = std::pow(q()(0),2);
    double q12 = std::pow(q()(1),2);
    double q22 = std::pow(q()(2),2);
    double q32 = std::pow(q()(3),2);

    (*this)()(0,0) = q02 + q12 - q22 - q32;
    (*this)()(0,1) = 2 * (+q()(1) * q()(2) - q()(0) * q()(3));
    (*this)()(0,2) = 2 * (+q()(0) * q()(2) + q()(1) * q()(3));
    (*this)()(1,0) = 2 * (+q()(0) * q()(3) + q()(1) * q()(2));
    (*this)()(1,1) = q02 - q12 + q22 - q32;
    (*this)()(1,2) = 2 * (-q()(0) * q()(1) + q()(2) * q()(3));
    (*this)()(2,0) = 2 * (-q()(0) * q()(2) + q()(1) * q()(3));
    (*this)()(2,1) = 2 * (+q()(0) * q()(1) + q()(2) * q()(3));
    (*this)()(2,2) = q02 - q12 - q22 + q32;
}
/* constructor based on Rodrigues parameters */

dcm::dcm(const rotv& rv) {
    double rv_norm = rv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        Eigen::Matrix3d rv_skew = tools::skew3(rv());
        this->get() = Eigen::Matrix3d::Identity() + rv_skew * (std::sin(rv_norm) / rv_norm) + rv_skew * rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2);
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        // convert from rotation vector to quaternion and then to dcm
        rodrigues temp(rv);
        (*this) = temp;
    }
}
/* constructor based on rotation vector */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

dcm::dcm(euler&& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    (*this) << + cp * cy, - cr*sy + sr*sp*cy, + sr*sy + cr*sp*cy,
            + cp * sy, + cr*cy + sr*sp*sy, - sr*cy + cr*sp*sy,
            - sp,      + sr * cp,          + cr * cp;
}
/* move constructor based on Euler angles */

dcm::dcm(rodrigues&& q) {
    double q02 = std::pow(q()(0),2);
    double q12 = std::pow(q()(1),2);
    double q22 = std::pow(q()(2),2);
    double q32 = std::pow(q()(3),2);

    (*this)()(0,0) = q02 + q12 - q22 - q32;
    (*this)()(0,1) = -2 * (-q()(1) * q()(2) + q()(0) * q()(3));
    (*this)()(0,2) = -2 * (-q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(1,0) = -2 * (-q()(0) * q()(3) - q()(1) * q()(2));
    (*this)()(1,1) = q02 - q12 + q22 - q32;
    (*this)()(1,2) = -2 * (+q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,0) = -2 * (+q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(2,1) = -2 * (-q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,2) = q02 - q12 - q22 + q32;
}
/* move constructor based on Rodrigues parameters */

dcm::dcm(rotv&& rv) {
    double rv_norm = rv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        Eigen::Matrix3d rv_skew = tools::skew3(rv());
        this->get() = Eigen::Matrix3d::Identity() + rv_skew * (std::sin(rv_norm) / rv_norm) + rv_skew * rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2);
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        // convert from rotation vector to quaternion and then to dcm
        rodrigues temp(rv);
        (*this) = temp;
    }
}
/* move constructor based on rotation vector */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

dcm& dcm::operator=(const Eigen::Matrix3d& op2) {
    this->get() = op2;
    this->normalize();
    return *this;
}
/* assignemnt operator based on 3x3 matrix */

dcm& dcm::operator=(const euler& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    (*this) << + cp * cy, - cr*sy + sr*sp*cy, + sr*sy + cr*sp*cy,
            + cp * sy, + cr*cy + sr*sp*sy, - sr*cy + cr*sp*sy,
            - sp,      + sr * cp,          + cr * cp;
    return *this;
}
/* assignment operator based on Euler angles */

dcm& dcm::operator=(const rodrigues& q) {
    double q02 = std::pow(q()(0),2);
    double q12 = std::pow(q()(1),2);
    double q22 = std::pow(q()(2),2);
    double q32 = std::pow(q()(3),2);

    (*this)()(0,0) = q02 + q12 - q22 - q32;
    (*this)()(0,1) = -2 * (-q()(1) * q()(2) + q()(0) * q()(3));
    (*this)()(0,2) = -2 * (-q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(1,0) = -2 * (-q()(0) * q()(3) - q()(1) * q()(2));
    (*this)()(1,1) = q02 - q12 + q22 - q32;
    (*this)()(1,2) = -2 * (+q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,0) = -2 * (+q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(2,1) = -2 * (-q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,2) = q02 - q12 - q22 + q32;
    return *this;
}
/* assignment operator based on Rodrigues parameters */

dcm& dcm::operator=(const rotv& rv) {
    double rv_norm = rv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        Eigen::Matrix3d rv_skew = tools::skew3(rv());
        this->get() = Eigen::Matrix3d::Identity() + rv_skew * (std::sin(rv_norm) / rv_norm) + rv_skew * rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2);
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        // convert from rotation vector to quaternion and then to dcm
        rodrigues temp(rv);
        (*this) = temp;
    }
    return *this;
}
/* assignment operator based on rotation vector */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

dcm& dcm::operator=(Eigen::Matrix3d&& op2) {
    this->get() = op2;
    this->normalize();
    return *this;
}
/* move assignment operator based on 3x3 matrix */

dcm& dcm::operator=(euler&& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    (*this) << + cp * cy, - cr*sy + sr*sp*cy, + sr*sy + cr*sp*cy,
            + cp * sy, + cr*cy + sr*sp*sy, - sr*cy + cr*sp*sy,
            - sp,      + sr * cp,          + cr * cp;
    return *this;
}
/* move assignment operator based on Euler angles */

dcm& dcm::operator=(rodrigues&& q) {
    double q02 = std::pow(q()(0),2);
    double q12 = std::pow(q()(1),2);
    double q22 = std::pow(q()(2),2);
    double q32 = std::pow(q()(3),2);

    (*this)()(0,0) = q02 + q12 - q22 - q32;
    (*this)()(0,1) = -2 * (-q()(1) * q()(2) + q()(0) * q()(3));
    (*this)()(0,2) = -2 * (-q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(1,0) = -2 * (-q()(0) * q()(3) - q()(1) * q()(2));
    (*this)()(1,1) = q02 - q12 + q22 - q32;
    (*this)()(1,2) = -2 * (+q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,0) = -2 * (+q()(0) * q()(2) - q()(1) * q()(3));
    (*this)()(2,1) = -2 * (-q()(0) * q()(1) - q()(2) * q()(3));
    (*this)()(2,2) = q02 - q12 - q22 + q32;
    return *this;
}
/* move assignment operator based on Rodrigues parameters */

dcm& dcm::operator=(rotv&& rv) {
    double rv_norm = rv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        Eigen::Matrix3d rv_skew = tools::skew3(rv());
        this->get() = Eigen::Matrix3d::Identity() + rv_skew * (std::sin(rv_norm) / rv_norm) + rv_skew * rv_skew * (1 - std::cos(rv_norm)) / std::pow(rv_norm, 2);
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        // convert from rotation vector to quaternion and then to dcm
        rodrigues temp(rv);
        (*this) = temp;
    }
    return *this;
}
/* move assignment operator based on rotation vector */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

dcm dcm::pow(const double& t) const {
    return this->log_map().pow(t).exp_map_dcm();
}
/* executes object rotation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object logarithmic map. */

dcm dcm::slerp(const dcm& R0, const dcm& R1, const double& t) {
    rotv delta_rotv((rotv(R0.inverse() * R1)).pow(t));
    return R0 * dcm(delta_rotv);
}
/* spherical linear interpolation, returns R0 for t=0 and R1 for t=1 */

dcm dcm::plus_right(const rotv& rv) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (rv().norm() < math::constant::PI()) {
        return (*this) * rv.exp_map_dcm();
    }
    else {
        double new_norm = math::constant::PI() * 2. - rv().norm();
        rotv new_rv(rv / rv().norm() * new_norm * (-1));
        return (*this) * new_rv.exp_map_dcm();
    }
}
/* right plus operator (input rotation located in local tangent space). */

dcm dcm::plus_left(const rotv& rv) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (rv().norm() < math::constant::PI()) {
        return rv.exp_map_dcm() * (*this);
    }
    else {
        double new_norm = math::constant::PI() * 2. - rv().norm();
        rotv new_rv(rv / rv().norm() * new_norm * (-1));
        return new_rv.exp_map_dcm() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space). */

/* ===== ===== Operations ===== ===== */
/* ================================== */

rotv dcm::minus_right(const dcm& R) const {
    return (R.inverse() * (*this)).log_map();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output rotation located in local tangent space) */

ang::rotv dcm::minus_left(const dcm& R) const {
    return ((*this) * R.inverse()).log_map();
}
/* left minus operator (output rotation located in global tangent space) */

rotv dcm::log_map() const {
    return rotv(*this);
}
/* logarithmic map that returns the rotation vector */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

so3_tangent_skew dcm::operator|(const so3_tangent_skew& w_skew) const {
    return so3_tangent_skew(this->get() * w_skew() * this->inverse()());
}
/* overloaded operator | (forward adjoint) */

so3_tangent dcm::operator|(const so3_tangent& w) const {
    return so3_tangent((*this) * w());
}
/* overloaded operator | (forward adjoint) */

so3_tangent_skew dcm::operator%(const so3_tangent_skew& w_skew) const {
    return so3_tangent_skew(this->inverse()() * w_skew() * this->get());
}
/* overloaded operator % (backward adjoint) */

so3_tangent dcm::operator%(const so3_tangent& w) const {
    return so3_tangent((*this) / w());
}
/* overloaded operator % (backward adjoint) */

/* ===== ===== Angular Velocity - Time Derivative ===== ===== */
/* ========================================================== */
so3_tangent dcm::dot2omegabody(const Eigen::Matrix3d& dcmdot) const {
    return so3_tangent::wedge(so3_tangent_skew((*this).transpose() * dcmdot));
}
/* obtains the body angular velocity from the direction cosine matrix and its time derivative. */

Eigen::Matrix3d dcm::omegabody2dot(const so3_tangent& w_body_rps) const {
    return this->get() * so3_tangent::hat_skew(w_body_rps)();
}
/* obtains the direction cosine matrix differentials with time based on the direction cosine matrix and the body angular velocity. */

so3_tangent dcm::dot2omegaspace(const Eigen::Matrix3d& dcmdot) const {
    return so3_tangent::wedge(so3_tangent_skew(dcmdot * (*this).transpose()));
}
/* obtains the space angular velocity from the direction cosine matrix and its time derivative. */

Eigen::Matrix3d dcm::omegaspace2dot(const so3_tangent& w_space_rps) const {
    return so3_tangent::hat_skew(w_space_rps)() * (*this);
}
/* obtains the direction cosine matrix differentials with time based on the direction cosine matrix and the space angular velocity. */

void dcm::normalize() {
    // algorithm should be iterative, but this is more than enough for our purposes
    Eigen::Vector3d a0 = (this->row(0) + this->row(1).cross(this->row(2))) * 0.5;
    Eigen::Vector3d a1 = (this->row(1) + this->row(2).cross(this->row(0))) * 0.5;
    Eigen::Vector3d a2 = (this->row(2) + this->row(0).cross(this->row(1))) * 0.5;
    this->row(0) = a0 / a0.norm();
    this->row(1) = a1 / a1.norm();
    this->row(2) = a2 / a2.norm();
}
/* normalize the rotation matrix ensuring it is orthonormal */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */
Eigen::Vector9d dcm::wedge(const dcm& R) {
    Eigen::Vector9d res;
    res.segment<3>(0) = R().row(0);
    res.segment<3>(3) = R().row(1);
    res.segment<3>(6) = R().row(2);
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the nine components of the
 * rotation matrix in matrix form and returns them in vector form, ordered by row. */

Eigen::Vector9d dcm::wedge(const Eigen::Matrix3d& R) {
    Eigen::Vector9d res;
    res.segment<3>(0) = R.row(0);
    res.segment<3>(3) = R.row(1);
    res.segment<3>(6) = R.row(2);
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the nine components of the
 * rotation matrix in matrix form and returns them in vector form, ordered by row. */

dcm dcm::hat(const Eigen::Vector9d& v) {
    Eigen::Matrix3d R;
    R.row(0) = v.segment<3>(0);
    R.row(1) = v.segment<3>(3);
    R.row(2) = v.segment<3>(6);
    return dcm(R);
}
/* although the hat operator usually applies to the tangent space, here it takes the nine components of the
 * rotation matrix in vector form, ordered by row, and returns them in matrix form (object). */

/* ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================= */
Eigen::Matrix3d dcm::jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - (*this)() * tools::skew3(vec);
}
/* returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(R * v) / dDeltarB,
 * (R plus DeltarB) * v = R * v + J * DeltarB */

Eigen::Matrix3d dcm::jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    //return (*this).inverse()() * tools::skew3(vec) * (*this)(); // this result is identical but more expensive to compute
    return tools::skew3((*this).inverse()() * vec);
}
/* returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(R / v) / dDeltarB,
 * (R plus DeltarB) / v = R / v + J * DeltarB */

Eigen::Matrix3d dcm::jac_right_log() const {
    return this->log_map().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(R)) / dDeltarv, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(R plus Deltarv) ~= Log(R) + (J * Deltarv) */

Eigen::Matrix3d dcm::jac_right_plus_wrt_first(const rotv& rv) const {
    return rv.exp_map_dcm().adjoint_matrix_backward();
}
/* returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
 * equal to d(R1 plus rv2) / dDeltar1B.
 * (R1 plus Deltar1B) plus rv2 = (R1 plus rv2) plus J * Deltar1B */

Eigen::Matrix3d dcm::jac_right_plus_wrt_second(const rotv& rv) const {
    return rv.jac_right();
}
/* returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
 * equal to d(R1 plus rv2) / dDeltar2B.
 * R1 plus (rv2 + Deltar2B) = (R1 plus rv2) plus J * Deltar2B */

Eigen::Matrix3d dcm::jac_right_minus_wrt_first(const dcm& R) const {
    return this->minus_right(R).jac_right_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the first element,
 * equal to d(R2 minus R1) / dDeltar2B.
 * (R2 plus Deltar2B) minus R1 = (R2 minus R1) + J * Deltar2B */

Eigen::Matrix3d dcm::jac_right_minus_wrt_second(const dcm& R) const {
    return - this->minus_right(R).jac_left_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the second element,
 * equal to d(R2 minus R1) / dDeltar1B.
 * R2 minus(R1 plus Deltar1B) = (R2 minus R1) + J * Deltar1B */

Eigen::Matrix3d dcm::jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - (*this)() * tools::skew3(w());
}
/* returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(AdR | w) / dDeltarB,
 * Ad(R plus DeltarB) | w = AdR | w + J * DeltarB */

Eigen::Matrix3d dcm::jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return tools::skew3((*this).inverse()() * w());
}
/* returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(AdR % w) / dDeltarB,
 * Ad(R plus DeltarB) % w = AdR % w + J * DeltarB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix3d dcm::jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - tools::skew3((*this)() * vec);
}
/* returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(R * v) / dDeltarN,
 * (DeltarN plus R) * v = R * v + J * DeltarN */

Eigen::Matrix3d dcm::jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return (*this)().transpose() * tools::skew3(vec);
}
/* returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(R / v) / dDeltarN,
 * (DeltarN plus R) / v = R / v + J * DeltarN */

Eigen::Matrix3d dcm::jac_left_log() const {
    return this->log_map().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(R)) / dDeltarv, which coincides with the inverse
 * of the left jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(Deltarv plus R) ~= Log(R) + (J * Deltarv) */

Eigen::Matrix3d dcm::jac_left_plus_wrt_first(const rotv& rv) const {
    return rv.jac_left();
}
/* returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
 * equal to d(rv1 plus R2) / dDeltar1N.
 * (rv1 + Deltar1N) plus R2 = J * Deltar1N plus (rv1 plus R2)  */

Eigen::Matrix3d dcm::jac_left_plus_wrt_second(const rotv& rv) const {
    return rv.exp_map_dcm().adjoint_matrix_forward();
}
/* returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
 * equal to d(rv1 plus R2) / dDeltar2N.
 * rv1 plus (Deltar2N plus R2) = J * Deltar2N plus (rv1 plus R2) */

Eigen::Matrix3d dcm::jac_left_minus_wrt_first(const dcm& R) const {
    return this->minus_left(R).jac_left_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the first element,
 * equal to d(R2 minus R1) / dDeltar2N.
 * (Deltar2N plus R2) minus R1 = (R2 minus R1) + J * Deltar2N */

Eigen::Matrix3d dcm::jac_left_minus_wrt_second(const dcm& R) const {
    return - this->minus_left(R).jac_right_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the second element,
 * equal to d(R2 minus R1) / dDeltar1N.
 * R2 minus (Deltar1N plus R1) = (R2 minus R1) + J * Deltar1N */

Eigen::Matrix3d dcm::jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - tools::skew3((*this)() * w());
}
/* returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(AdR | w) / dDeltarN,
 * Ad(DeltarN plus R) | w = AdR | w + J * DeltarN */

Eigen::Matrix3d dcm::jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return (*this)().transpose() * tools::skew3(w());
}
/* returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(AdR % w) / dDeltarN,
 * Ad(DeltarN plus R) % w = AdR % w + J * DeltarN */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d dcm::jac_euclidean_forward_rotation_wrt_vector() const {
    return (*this)();
}
/* returns the jacobian of the forward rotation action with respect to the vector, equal to d(R * v) / dDeltav,
 * R * v = J * v
 * R * (v + Deltav) = R * v + J * Deltav */

Eigen::Matrix3d dcm::jac_euclidean_backward_rotation_wrt_vector() const {
    return this->inverse()();
}
/* returns the jacobian of the backward rotation action with respect to the vector, equal to d(R / v) / dDeltav,
 * R / v = J * v
 * R / (v + Deltav) = R / v + J * Deltav */

Eigen::Matrix3d dcm::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return (*this)();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdR | w) / dDeltaw,
 * AdR | w = J * w
 * AdR | (w + Deltaw) = AdR | w + J * Deltaw */

Eigen::Matrix3d dcm::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdR % w) / dDeltaw,
 * AdR % w = J * w
 * AdR % (w + Deltaw) = AdR % w + J * Deltaw */

Eigen::Matrix39d dcm::jac_euclidean_forward_rotation_wrt_dcm(const Eigen::Vector3d& v) const {
    Eigen::Matrix39d res = Eigen::Matrix39d::Zero();
    res.block<1,3>(0,0) = v.transpose();
    res.block<1,3>(1,3) = v.transpose();
    res.block<1,3>(2,6) = v.transpose();
    return res;
}
/* returns the jacobian [3x9] of a forward rotation with respect to the rotation matrix, equal to d(R * v)/dR.
 * The forward rotation is linear on the rotation matrix, so
 * R * v = d(R * v)/dR |R*v * R
 * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [R * exp(Delta r)] * v = R * v + d(R * v)/dR |R*v * {[R * exp(Delta r)] - R}
 * [exp(Delta r) * R] * v = R * v + d(R * v)/dR |R*v * {[exp(Delta r) * R] - R}
 * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
 * to each other, and never exp(Delta r) or Delta R. */

Eigen::Matrix39d dcm::jac_euclidean_backward_rotation_wrt_dcm(const Eigen::Vector3d& v) const {
    Eigen::Matrix39d res = Eigen::Matrix39d::Zero();
    res.block<3,3>(0,0).diagonal() = v(0) * Eigen::Vector3d::Ones();
    res.block<3,3>(0,3).diagonal() = v(1) * Eigen::Vector3d::Ones();
    res.block<3,3>(0,6).diagonal() = v(2) * Eigen::Vector3d::Ones();
    return res;
}
/* returns the jacobian [3x9] of a backward rotation with respect to the rotation matrix, equal to d(R / v)/dR.
 * The backward rotation is linear on the rotation matrix, so
 * R / v = d(R / v)/dR |R/v * R
 * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [R * exp(Delta r)] / v = R / v + d(R / v)/dR |R/v * {[R * exp(Delta r)] - R}
 * [exp(Delta r) * R] / v = R / v + d(R / v)/dR |R/v * {[exp(Delta r) * R] - R}
 * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
 * to each other, and never exp(Delta r) or Delta R. */

Eigen::Matrix39d dcm::jac_euclidean_forward_adjoint_wrt_dcm(const so3_tangent& w) const {
    return this->jac_euclidean_forward_rotation_wrt_dcm(w());
}
/* returns the jacobian [3x9] of a forward adjoint with respect to the rotation matrix, equal to d(R | w)/dR.
 * The forward adjoint is linear on the rotation matrix, so
 * R | w = d(R | w)/dR |R|w * R
 * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [R * exp(Delta r)] | w = R | w + d(R | w)/dR |R|w * {[R * exp(Delta r)] - R}
 * [exp(Delta r) * R] | w = R | w + d(R | w)/dR |R|w * {[exp(Delta r) * R] - R}
 * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
 * to each other, and never exp(Delta r) or Delta R. */

Eigen::Matrix39d dcm::jac_euclidean_backward_adjoint_wrt_dcm(const so3_tangent& w) const {
    return this->jac_euclidean_backward_rotation_wrt_dcm(w());
}
/* returns the jacobian [3x9] of a backward adjoint with respect to the rotation matrix, equal to d(R % w)/dR.
 * The backward adjoint is linear on the rotation matrix, so
 * R % w = d(R % w)/dR |R%w * R
 * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [R * exp(Delta r)] % w = R % w + d(R % w)/dR |R%w * {[R * exp(Delta r)] - R}
 * [exp(Delta r) * R] % w = R % w + d(R % w)/dR |R%w * {[exp(Delta r) * R] - R}
 * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
 * to each other, and never exp(Delta r) or Delta R. */



















