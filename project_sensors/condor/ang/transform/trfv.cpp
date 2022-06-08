#include "trfv.h"
#include "speu_rodrigues.h"
#include "speu_dcm.h"
#include "homogeneous.h"
#include "screw.h"
#include "dual.h"
#include "se3_tangent.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS TRANSFORM VECTOR
// ======================
// ======================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

trfv::trfv(double s1, double s2, double s3, double r1, double r2, double r3) {
    *this << s1, s2, s3, r1, r2, r3;
}
/* constructor based on transform vector components */

trfv::trfv(const rotv& rotv, const Eigen::Vector3d& T) {
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * T;
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(T);
    }
}
/* constructor based on rotation vector and translation vector */

trfv::trfv(const speu_rodrigues& gq) {
    *this = trfv(gq.get_rotv(), gq.get_T());
}
/* constructor based on special Euclidean (rodrigues) */

trfv::trfv(const speu_dcm& gR) {
    *this = trfv(gR.get_rotv(), gR.get_T());
}
/* constructor based on special Euclidean (dcm) */

trfv::trfv(const ang::homogeneous& M) {
    *this = trfv(M.get_rotv(), M.get_T());
}
/* constructor based on homogeneous */

trfv::trfv(const screw& S) {
    *this = trfv(S.get_rotv(), S.get_T());
}
/* constructor based on screw */

trfv::trfv(const dual& Z) {
    *this = trfv(Z.get_rotv(), Z.get_T());
}
/* constructor based on unit dual quaternion */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

trfv::trfv(rotv&& rotv, Eigen::Vector3d&& T) {
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * T;
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(T);
    }
}
/* move constructor based on rotation vector and translation vector */

trfv::trfv(speu_rodrigues&& gq) {
    *this = trfv(gq.get_rotv(), gq.get_T());
}
/* move constructor based on special Euclidean (rodrigues) */

trfv::trfv(speu_dcm&& gR) {
    *this = trfv(gR.get_rotv(), gR.get_T());
}
/* move constructor based on special Euclidean (dcm) */

trfv::trfv(ang::homogeneous&& M) {
    *this = trfv(M.get_rotv(), M.get_T());
}
/* move constructor based on homogeneous */

trfv::trfv(screw&& S) {
    *this = trfv(S.get_rotv(), S.get_T());
}
/* move constructor based on screw */

trfv::trfv(dual&& Z) {
    *this = trfv(Z.get_rotv(), Z.get_T());
}
/* move constructor based on unit dual quaternion */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

trfv& trfv::operator=(const Eigen::Vector6d& op2) {
    static_cast<Eigen::Vector6d&>(*this) = op2;
    return *this;
}
/* assignment operator based on size 6 vector */

trfv& trfv::operator=(const speu_rodrigues& gq) {
    rotv rotv(gq.get_rodrigues());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * gq.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(gq.get_T());
    }
    return *this;
}
/* assignment operator based on special Euclidean (rodrigues) */

trfv& trfv::operator=(const speu_dcm& gR) {
    rotv rotv(gR.get_dcm());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * gR.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(gR.get_T());
    }
    return *this;
}
/* assignment operatoro based on special Euclidean (dcm) */

trfv& trfv::operator=(const ang::homogeneous& M) {
    rotv rotv(M.get_dcm());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * M.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(M.get_T());
    }
    return *this;
}
/* assignment operator based on homogeneous */

trfv& trfv::operator=(const screw& S) {
    *this = trfv(S.get_rotv(), S.get_T());
    return *this;
}
/* assignment operator = based on screw */

trfv& trfv::operator=(const dual& Z) {
    *this = trfv(Z.get_rotv(), Z.get_T());
    return *this;
}
/* assignment operator based on unit dual quaternion */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */
trfv& trfv::operator=(Eigen::Vector6d&& op2) {
    static_cast<Eigen::Vector6d&>(*this) = op2;
    return *this;
}
/* move assignment operator based on size 6 vector */

trfv& trfv::operator=(speu_rodrigues&& gq) {
    rotv rotv(gq.get_rodrigues());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * gq.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(gq.get_T());
    }
    return *this;
}
/* move assignment operator based on special Euclidean (rodrigues) */

trfv& trfv::operator=(speu_dcm&& gR) {
    rotv rotv(gR.get_dcm());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * gR.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(gR.get_T());
    }
    return *this;
}
/* move assignment operator based on special Euclidean (dcm) */

trfv& trfv::operator=(ang::homogeneous&& M) {
    rotv rotv(M.get_dcm());
    this->tail<3>() = rotv();
    double theta = rotv().norm();
    Eigen::Matrix3d rm = tools::skew3(rotv());

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        this->head<3>() = (Eigen::Matrix3d::Identity()- 0.5 * rm + (1./12.) * (rm * rm)) * M.get_T();
    }
    else {
        double theta_sq = theta * theta;
        Eigen::Matrix3d V = ((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(rotv())) / theta_sq;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> Odecom(V); // this decomposition works, LDLT does not
        this->head<3>() = Odecom.solve(M.get_T());
    }
    return *this;
}
/* move assignment operator based on homogeneous */

trfv& trfv::operator=(screw&& S) {
    *this = trfv(S.get_rotv(), S.get_T());
    return *this;
}
/* move assignment operator = based on screw */

trfv& trfv::operator=(dual&& Z) {
    *this = trfv(Z.get_rotv(), Z.get_T());
    return *this;
}
/* move assignment operator based on unit dual quaternion */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

trfv trfv::pow(const double& t) const {
    // return trfv(this->get() * t); // THIS DOES NOT WORK, IT IS NOT THE SAME
    return this->explog_map_screw().pow(t).explog_map_trfv();
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns map of the power function applied to the object screw exponential-log map. */

trfv trfv::sclerp(const trfv& tau0, const trfv& tau1, const double& t) {
    screw delta_screw((screw(tau0.inverse() * tau1)).pow(t));
    return tau0 * trfv(delta_screw);
}
/* screw linear interpolation, returns tau0 for t=0 and tau1 for t=1 */

trfv trfv::plus_right(const screw& S) const {
    return (*this) * S.explog_map_trfv();
}
/* right plus operator (input rotation located in local tangent space) */

trfv trfv::plus_left(const screw& S) const {
    return S.explog_map_trfv() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d trfv::operator *(const Eigen::Vector3d& vecin) const {
    return speu_rodrigues(*this) * vecin;
}
/* overloaded operator * (forward transformation of point, not vector) */

Eigen::Vector3d trfv::operator /(const Eigen::Vector3d& vecin) const {
    return speu_rodrigues(*this) / vecin;
}
/* overloaded operator / (backward transformation of point, not vector) */

Eigen::Vector3d trfv::operator^(const Eigen::Vector3d& vecin) const {
    rotv rv(this->tail<3>());
    return rv * vecin;
}
/* overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */

Eigen::Vector3d trfv::operator&(const Eigen::Vector3d& vecin) const {
    rotv rv(this->tail<3>());
    return rv / vecin;
}
/* overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */

screw trfv::minus_right_screw(const trfv& tau) const {
    return (tau.inverse() * (*this)).explog_map_screw();
}
/* minus operator (output transformation located in local tangent space) */

screw trfv::minus_left_screw(const trfv& tau) const {
    return ((*this) * tau.inverse()).explog_map_screw();
}
/* left minus operator (output transformation located in global tangent space) */

speu_rodrigues trfv::exp_map_speu_rodrigues() const {
    return speu_rodrigues(*this);
}
/* exponential map that returns special euclidean rodrigues parameters */

speu_dcm trfv::exp_map_speu_dcm() const {
    return speu_dcm(*this);
}
/* exponential map that returns special euclidean direction cosine matrix */

homogeneous trfv::exp_map_homogeneous() const {
    return ang::homogeneous(*this);
}
/* exponential map that returns homogeneous matrix */

dual trfv::exp_map_dual() const {
    return dual(*this);
}
/* exponential map that returns the unit dual quaternion */

screw trfv::explog_map_screw() const {
    return screw(*this);
}
/* exponential and logarithmic map (they are the same in this case) that returns the screw */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

se3_tangent trfv::operator|(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = (*this) ^ xi().tail<3>();
    res().head<3>() = ((*this) ^ xi().head<3>()) + this->get_T().cross(res().tail<3>()); // do not remove parenthesis
    return res;
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix6d trfv::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;
    ang::dcm R(this->get_rotv());

    res.topLeftCorner<3,3>()     = R();
    res.topRightCorner<3,3>()    = tools::skew3(this->get_T()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R();

    return res;
}
/* returns forward adjoint matrix */

se3_tangent trfv::operator%(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = (*this) & xi().tail<3>();
    res().head<3>() = ((*this) & xi().head<3>()) - ((*this) & this->get_T().cross(xi().tail<3>())); // do not remove parenthesis
    return res;
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix6d trfv::adjoint_matrix_backward() const {
    Eigen::Matrix6d res;
    ang::dcm R(this->get_rotv());

    res.topLeftCorner<3,3>()     = R().transpose();
    res.topRightCorner<3,3>()    = - R().transpose() * tools::skew3(this->get_T());
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R().transpose();

    return res;
}
/* returns backward adjoint matrix */

/* ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
/* ================================================================== */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

Eigen::Vector6d trfv::wedge(const trfv& tau) {
    return tau();
}
/* although the wedge operator usually applies to the tangent space, here it takes the six components of the
 * transform vector and returns them in vector form. */

trfv trfv::hat(const Eigen::Vector6d& v) {
    return trfv(v);
}
/* although the hat operator usually applies to the tangent space, here it takes the six components of the
 * transform vector in vector form, and returns them in transform vector form (object). */

/* ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================================== */

Eigen::Matrix36d trfv::jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(this->get_rotv())();
    res.leftCols<3>()  = R;;
    res.rightCols<3>() = - R * tools::skew3(p);
    return res;
}
/* returns the right jacobian of the forward motion action with respect to the motion, equal to d(tau * p) / dDeltatauB,
 * (tau plus DeltatauB) * p = tau * p + J * DeltatauB */

Eigen::Matrix36d trfv::jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(this->get_rotv())();
    res.leftCols<3>()  = - Eigen::Matrix3d::Identity();
    //res.rightCols<3>() = R.transpose() * tools::skew3(p - this->get_T()) * R; // this result is identical but more expensive to compute
    res.rightCols<3>() = tools::skew3(dcm(this->get_inverse_rotv())() * (p - this->get_T()));
    return res;
}
/* returns the right jacobian of the backward motion action with respect to the motion, equal to d(tau / p) / dDeltatauB,
 * (tau plus DeltatauB) / p = tau / p + J * DeltatauB */

Eigen::Matrix6d trfv::jac_right_log() const {
    return this->jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(tau)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(tau plus Deltatau) ~= Log(tau) + (J * Deltatau) */

Eigen::Matrix6d trfv::jac_right_plus_wrt_first(const trfv& tau) const {
    return tau.adjoint_matrix_backward();
}
/* returns the right jacobian of the motion right plus operator with respect to the group object (first element),
 * equal to d(tau1 plus tau2) / dDeltatau1B.
 * (tau1 plus Deltatau1B) plus tau2 = (tau1 plus tau2) plus J * Deltatau1B */

Eigen::Matrix6d trfv::jac_right_plus_wrt_second(const trfv& tau) const {
    return tau.jac_right();
}
/* returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
 * equal to d(tau1 plus tau2) / dDeltatau2B.
 * tau1 plus (tau2 + Deltatau2B) = (tau1 plus tau2) plus J * Deltatau2B */

Eigen::Matrix6d trfv::jac_right_minus_wrt_first(const trfv& tau) const {
    return this->minus_right_trfv(tau).jac_right_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the first element,
 * equal to d(tau2 minus tau1) / dDeltatau2B.
 * (tau2 plus Deltatau2B) minus tau1 = (tau2 minus tau1) + J * Deltatau2B */

Eigen::Matrix6d trfv::jac_right_minus_wrt_second(const trfv& tau) const {
    return - this->minus_right_trfv(tau).jac_left_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the second element,
 * equal to d(tau2 minus tau1) / dDeltatau1B.
 * tau2 minus(tau1 plus Deltatau1B) = (tau2 minus tau1) + J * Deltatau1B */

Eigen::Matrix6d trfv::jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    dcm R(this->get_rotv());
    Eigen::Vector3d T = this->get_T();
    ang::so3_tangent w = xi.get_w();
    res.topLeftCorner<3,3>()     = - tools::skew3(R() * w()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - R() * tools::skew3(xi.get_vi()) - tools::skew3(T) * R() * tools::skew3(w());
    res.bottomRightCorner<3,3>() = - R() * tools::skew3(w());
    return res;
}
/* returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(Adtau | xi) / dDeltatauB,
 * Ad(tau plus DeltatauB) | xi = Adtau | xi + J * DeltatauB */

Eigen::Matrix6d trfv::jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    dcm R(this->get_rotv());
    Eigen::Matrix3d Rtr = R().transpose();
    Eigen::Vector3d T = this->get_T();

    res.topLeftCorner<3,3>()     = Rtr * tools::skew3(w()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = tools::skew3(Rtr * (xi.get_vi() - tools::skew3(T) * w()));
    res.bottomRightCorner<3,3>() = tools::skew3(Rtr * w());

    return res;
}
/* returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(Adtau % xi) / dDeltatauB,
 * Ad(tau plus DeltatauB) % xi = Adtau % xi + J * DeltatauB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix36d trfv::jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(this->get_rotv())();
    res.leftCols<3>()  = Eigen::Matrix3d::Identity();
    res.rightCols<3>() = - tools::skew3(R * p + this->get_T());
    return res;
}
/* returns the left jacobian of the forward motion action with respect to the motion, equal to d(tau * p) / dDeltatauE,
 * (DeltatauE plus tau) * p = tau * p + J * DeltatauE */

Eigen::Matrix36d trfv::jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d Rtr = dcm(this->get_rotv())().transpose();
    res.leftCols<3>()  = - Rtr;
    res.rightCols<3>() = Rtr * tools::skew3(p);
    return res;
}
/* returns the left jacobian of the backward motion action with respect to the motion, equal to d(tau / p) / dDeltatauE,
 * (DeltatauE plus tau) / p = tau / p + J * DeltatauE */

Eigen::Matrix6d trfv::jac_left_log() const {
    return this->jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(tau)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(Deltatau plus tau) ~= Log(tau) + (J * Deltatau) */

Eigen::Matrix6d trfv::jac_left_plus_wrt_first(const trfv& tau) const {
    return tau.jac_left();
}
/* returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
 * equal to d(tau1 plus tau2) / dDeltatau1E.
 * (tau1 + Deltatau1E) plus tau2 = J * Deltatau1E plus (tau1 plus tau2)  */

Eigen::Matrix6d trfv::jac_left_plus_wrt_second(const trfv& tau) const {
    return tau.adjoint_matrix_forward();
}
/* returns the left jacobian of the motion left plus operator with respect to the group object (second element),
 * equal to d(tau1 plus tau2) / dDeltatau2E.
 * tau1 plus (Deltatau2E plus tau2) = J * Deltatau2E plus (tau1 plus tau2) */

Eigen::Matrix6d trfv::jac_left_minus_wrt_first(const trfv& tau) const {
    return this->minus_left_trfv(tau).jac_left_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the first element,
 * equal to d(tau2 minus tau1) / dDeltatau2E.
 * (Deltatau2E plus tau2) minus tau1 = (tau2 minus tau1) + J * Deltatau2E */

Eigen::Matrix6d trfv::jac_left_minus_wrt_second(const trfv& tau) const {
    return - this->minus_left_trfv(tau).jac_right_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the second element,
 * equal to d(tau2 minus tau1) / dDeltatau1E.
 * tau2 minus (Deltatau1E plus tau1) = (tau2 minus tau1) + J * Deltatau1E */

Eigen::Matrix6d trfv::jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d A = tools::skew3(dcm(this->get_rotv())() * xi.get_w()());
    res.topLeftCorner<3,3>()     = - A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - tools::skew3(dcm(this->get_rotv())() * xi.get_vi()) - tools::skew3(this->get_T()) * A + A * tools::skew3(this->get_T());
    res.bottomRightCorner<3,3>() = - A;
    return res;
}
/* returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adtau | xi) / dDeltatauE,
 * Ad(DeltatauE plus tau) | xi = Adtau | xi + J * DeltatauE */

Eigen::Matrix6d trfv::jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = dcm(this->get_rotv())().transpose();
    Eigen::Matrix3d A = Rtr * tools::skew3(w());

    res.topLeftCorner<3,3>()     = A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = Rtr * (tools::skew3(xi.get_vi()) - tools::skew3(this->get_T()) * tools::skew3(w()));
    res.bottomRightCorner<3,3>() = A;

    return res;
}
/* returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adtau % xi) / dDeltatauE,
 * Ad(DeltatauE plus tau) % xi = Adtau % xi + J * DeltatauE */

/* ===== ===== Exponential (Right and Left) Jacobians ===== ===== */
/* ============================================================== */

Eigen::Matrix3d trfv::jac_block() const {
    rotv rv(this->tail<3>());
    Eigen::Vector3d rho(this->head<3>());

    double theta     = rv().norm();
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    double theta_sq  = theta * theta;
    double theta_cub = theta_sq * theta;
    double factor    = (1.0 - 0.5 * theta_sq - cos_theta) / (theta_sq * theta_sq);

    Eigen::Matrix3d rv_skew        = tools::skew3(rv());
    Eigen::Matrix3d rho_skew       = tools::skew3(rho);
    Eigen::Matrix3d rv_sq_skew     = rv_skew * rv_skew;
    Eigen::Matrix3d rv_rho_skew    = rv_skew * rho_skew;
    Eigen::Matrix3d rho_rv_skew    = rho_skew * rv_skew;
    Eigen::Matrix3d rv_rho_rv_skew = rv_rho_skew * rv_skew;

    return rho_skew * 0.5
           + (theta - sin_theta) / (theta_cub) * (rv_rho_skew + rho_rv_skew + rv_rho_rv_skew)
           - factor * (rv_skew * rv_rho_skew + rho_rv_skew * rv_skew - 3.0 * rv_rho_rv_skew)
           - 0.5 * (factor - 3.0 * (theta - sin_theta - theta_cub / 6.0) / (theta_sq * theta_cub)) * (rv_rho_skew * rv_sq_skew + rv_sq_skew * rho_rv_skew);
}
/* returns a [3x3] matrix which is part of the right and left jacobians as well as their inverses. Should be private. */

Eigen::Matrix6d trfv::jac_right() const {
    trfv tau_neg(- this->get());
    return tau_neg.jac_left();
}
/* returns the jacobian [6x6] of the exponential of the transform vector with respect to the transform vector itself,
 * also known as the right Jacobian of SE(3).
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * exp(tau + Deltatau) ~= exp(tau) plus (JR * Deltatau) */

Eigen::Matrix6d trfv::jac_right_inv() const {
    trfv tau_neg(- this->get());
    return tau_neg.jac_left_inv();
}
/* returns the inverse jacobian [6x6] of the exponential of the transform vector with respect to the transform vector
 * itself, also known as the inverse of the right Jacobian of SE(3).
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * tau + jacRinv * Deltatau ~= log[exp(tau) plus Deltatau] */

Eigen::Matrix6d trfv::jac_left() const {
    Eigen::Matrix6d J;
    J.block<3,3>(0,0) = this->get_rotv().jac_left();
    J.block<3,3>(0,3) = this->jac_block();
    J.block<3,3>(3,0) = Eigen::Matrix3d::Zero();
    J.block<3,3>(3,3) = J.block<3,3>(0,0);
    return J;
}
/* returns the jacobian [6x6] of the exponential of the transform vector with respect to the transform vector itself,
 * also known as the left Jacobian of SE(3).
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * exp(tau + Deltatau) ~= (JL * Deltatau) plus_left exp(tau) */

Eigen::Matrix6d trfv::jac_left_inv() const {
    Eigen::Matrix6d J;
    J.block<3,3>(0,0) = this->get_rotv().jac_left_inv();
    J.block<3,3>(0,3) = - J.block<3,3>(0,0) * this->jac_block() * J.block<3,3>(0,0);
    J.block<3,3>(3,0) = Eigen::Matrix3d::Zero();
    J.block<3,3>(3,3) = J.block<3,3>(0,0);
    return J;
}
/* returns the inverse jacobian [6x6] of the exponential of the transform vector with respect to the transform vector
 * itself, also known as the inverse of the left Jacobian of SE(3).
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * tau + JLinv * Deltatau ~= log[Deltatau plus_left exp(tau)] */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d trfv::jac_euclidean_forward_motion_wrt_point() const {
    return dcm(this->get_rotv())();
}
/* returns the jacobian of the forward motion action with respect to the point, equal to d(tau * p) / dDeltap,
 * tau * (p + Deltap) = tau * p + J * Deltap */

Eigen::Matrix3d trfv::jac_euclidean_backward_motion_wrt_point() const {
    return dcm(this->get_inverse_rotv())();
}
/* returns the jacobian of the backward motion action with respect to the point, equal to d(tau / p) / dDeltap,
 * tau / (p + Deltap) = tau / p + J * Deltap */

Eigen::Matrix6d trfv::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_forward();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adtau | xi) / dDeltaxi,
 * Adtau | (xi + Deltaxi) = Adtau | xi + J * Deltaxi */

Eigen::Matrix6d trfv::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adtau % xi) / dDeltaxi,
 * Adtau % (xi + Deltaxi) = Adtau % xi + J * Deltaxi */

Eigen::Matrix36d trfv::jac_euclidean_forward_motion_wrt_trfv(const Eigen::Vector3d& p) const {
    dcm R(this->get_rotv());
    Eigen::Matrix3d Jleft = this->get_rotv().jac_left();
    Eigen::Matrix36d J;
    J.block<3,3>(0,0) = Jleft;
    J.block<3,3>(0,3) = - tools::skew3(dcm(this->get_rotv())() * p + this->get_T()) * Jleft + this->jac_block();
    return J;
}
/* returns the jacobian [3x6] of a forward rotation with respect to the transform vector, equal to d(tau * p)/dtau.
 * The forward transformation is NOT linear on the transform vector, so
 * tau * p != d(tau * p)/dtau |tau*p * tau
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * (tau + Delta tau) * p ~= tau * p + d(tau * p)/dtau |tau*p * Delta tau
 * Note that the increment is Delta tau. */

Eigen::Matrix36d trfv::jac_euclidean_forward_motion_wrt_trfv_bis(const Eigen::Vector3d& p) const {
    dcm R(this->get_rotv());
    trfv tau_neg(- this->get());
    Eigen::Matrix3d Jleftneg = tau_neg.get_rotv().jac_left();
    Eigen::Matrix36d J;
    J.block<3,3>(0,0) = R() * Jleftneg;
    J.block<3,3>(0,3) = - R() * tools::skew3(p) * Jleftneg + R() * tau_neg.jac_block();
    return J;
}
/* returns the jacobian [3x6] of a forward transformation with respect to the transform vector, equal to d(tau * p)/dtau.
 * The forward transformation is NOT linear on the transform vector, so
 * tau * p != d(tau * p)/dtau |tau*p * tau
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * (tau + Delta tau) * p ~= tau * p + d(tau * p)/dtau |tau*p * Delta tau
 * Note that the increment is Delta tau. */

Eigen::Matrix36d trfv::jac_euclidean_backward_motion_wrt_trfv_zero(const Eigen::Vector3d& p) const {
    dcm R(this->get_rotv());
    trfv tau_neg(- this->get());
    Eigen::Matrix3d Jleftneg = tau_neg.get_rotv().jac_left();
    Eigen::Matrix36d J;
    J.block<3,3>(0,0) = - Jleftneg;
    J.block<3,3>(0,3) = - tau_neg.jac_block() - tools::skew3(R().transpose() * (this->get_T() - p)) * Jleftneg;
    return J;
}
/* returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d[tau / p]/dtau.
 * The backward transformation is NOT linear on the transform vector, so
 * tau / p != d(tau / p)/dtau |tau/p * tau
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
 * Note that the increment is Delta tau. */

Eigen::Matrix36d trfv::jac_euclidean_backward_motion_wrt_trfv(const Eigen::Vector3d& p) const {
    return - this->inverse().jac_euclidean_forward_motion_wrt_trfv(p);
}
/* returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d[tau / p]/dtau.
 * The backward transformation is NOT linear on the transform vector, so
 * tau / p != d(tau / p)/dtau |tau/p * tau
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
 * Note that the increment is Delta tau. */

Eigen::Matrix36d trfv::jac_euclidean_backward_motion_wrt_trfv_bis(const Eigen::Vector3d& p) const {
    return - this->inverse().jac_euclidean_forward_motion_wrt_trfv_bis(p);
}
/* returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d(tau / p)/dtau.
 * The backward transformation is NOT linear on the transform vector, so
 * tau / p != d(tau / p)/dtau |tau/p * tau
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
 * Note that the increment is Delta tau. */

Eigen::Matrix3d trfv::jac_euclidean_forward_motion_T_wrt_rotv() const {
    Eigen::Vector3d s(this->head<3>());
    rotv rv(this->tail<3>());

    double theta       = rv().norm();
    double theta_sq    = theta * theta;
    double theta_4     = theta_sq * theta_sq;
    Eigen::Vector3d rs = rv().cross(s);

    Eigen::Vector3d part_one   = tools::itself_by_transpose(rv()) * s;
    double part_three          = 1. / theta_sq;
    Eigen::Vector3d part_five  = rv().cross(s);
    Eigen::Vector3d part_nine  = rv.exp_map_dcm()() * rv().cross(s);

    Eigen::Matrix3d J_one      = rv() * s.transpose() + Eigen::Matrix3d::Identity() * rv().dot(s);
    Eigen::RowVector3d J_three = - 2 * rv().transpose() / theta_4;
    Eigen::Matrix3d J_five     = - tools::skew3(s);
    Eigen::Matrix3d J_nine     = rv.jac_euclidean_forward_rotation_wrt_rotv(rs) - rv.exp_map_dcm()() * tools::skew3(s);

    Eigen::Matrix3d J_four = J_one * part_three + part_one * J_three;
    // The second part is symmetric (looks like coincidence, but likely it is not), so it is the same as
    //Eigen::Matrix3d J_four = J_one * part_three + J_three.transpose() * part_one.transpose();

    Eigen::Matrix3d J_six = J_five * part_three + part_five * J_three;
    // The second part is NOT symmetric, so it is NOT the same as
    //Eigen::Matrix3d J_six = J_five * part_three + J_three.transpose() * part_five.transpose();

    Eigen::Matrix3d J_ten = J_nine * part_three + part_nine * J_three;

    return J_six - J_ten + J_four;
};
/* returns the jacobian [3x3] of the translation vector of a forward motion with respect to the rotation vector,
 * equal to d(exp(tau).get_T)/drv.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau + Deltatau).get_T() ~= exp(tau).get_T + d(exp(tau).get_T())/drv |[exp(tau)] * Deltarv.
 * Note that the jacobian is evaluated at |[exp(tau)]. */

Eigen::Matrix3d trfv::jac_euclidean_backward_motion_T_wrt_rotv() const {
    trfv tau_inv(this->inverse());
    return - tau_inv.jac_euclidean_forward_motion_T_wrt_rotv();
}
/* returns the jacobian [3x3] of the translation vector of a backward motion with respect to the rotation vector,
 * equal to d(exp(tau.inv()).get_T)/drv.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau.inv() + Deltatau).get_T() ~= exp(tau.inv()).get_T + d(exp(tau.inv()).get_T())/drv |[exp(tau.inv())] * Deltarv.
 * Note that the jacobian is evaluated at |[exp(tau.inv())]. */

Eigen::Matrix3d trfv::jac_euclidean_forward_motion_wrt_rotv(const Eigen::Vector3d& p) const {
    return this->jac_euclidean_forward_motion_T_wrt_rotv() + this->get_rotv().jac_euclidean_forward_rotation_wrt_rotv(p);
}
/* returns the jacobian [3x3] of a forward motion with respect to the rotation vector, equal to d(exp(tau) * p)/drv.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/drv |[exp(tau) * p] * Deltarv.
 * Note that the jacobian is evaluated at |[exp(tau)*p].*/

Eigen::Matrix3d trfv::jac_euclidean_backward_motion_wrt_rotv(const Eigen::Vector3d& p) const {
    return this->jac_euclidean_backward_motion_T_wrt_rotv() + this->get_rotv().jac_euclidean_backward_rotation_wrt_rotv(p);
}
/* returns the jacobian [3x3] of a backward motion with respect to the rotation vector, equal to d(exp(tau) / p)/drv.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) / p)/drv |[exp(tau)/p] * Deltarv.
 * Note that the jacobian is evaluated at |[exp(tau)/p].*/

Eigen::Matrix3d trfv::jac_euclidean_forward_motion_wrt_s(const Eigen::Vector3d& p) const {
    rotv rv(this->tail<3>());
    double theta = rv().norm();
    return ((Eigen::Matrix3d::Identity() - dcm(rv)()) * tools::skew3(rv()) + tools::itself_by_transpose(rv())) / theta / theta;
}
/* returns the jacobian [3x3] of a forward motion with respect to the translation vector s, equal to d(exp(tau) * p)/ds.
 * As the forward transformation IS linear on the translation vector, the following expression is true:
 * exp(tau + Deltatau) * p ==  exp(tau) * p + d(exp(tau) * p)/ds |[exp(tau)*p] * Deltas.
 * Note that the jacobian is evaluated at |[exp(tau)*p]. */

Eigen::Matrix3d trfv::jac_euclidean_backward_motion_wrt_s(const Eigen::Vector3d& p) const {
    trfv tau_inv(this->inverse());
    return - tau_inv.jac_euclidean_forward_motion_wrt_s(p);
}
/* returns the jacobian [3x3] of a backward motion with respect to the translation vector s, equal to d(exp(tau) / p)/ds.
 * As the forward transformation IS linear on the translation vector, the following expression is true:
 * exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) * p)/ds |[exp(tau)/p] * Deltas.
 * Note that the jacobian is evaluated at |[exp(tau)/p]. */

Eigen::Matrix36d trfv::jac_euclidean_forward_motion_wrt_trfv_tri(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d J;
    J.block<3,3>(0,0) = this->jac_euclidean_forward_motion_wrt_s(p);
    J.block<3,3>(0,3) = this->jac_euclidean_forward_motion_wrt_rotv(p);
    return J;
}
/* returns the jacobian [3x6] of a forward motion with respect to the transform vector, equal to d(exp(tau) * p)/dtau.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/dtau |[exp(tau)*p] * Deltatau.
 * Note that the jacobian is evaluated at |[exp(tau)*p]. */

Eigen::Matrix36d trfv::jac_euclidean_backward_motion_wrt_trfv_tri(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d J;
    J.block<3,3>(0,0) = this->jac_euclidean_backward_motion_wrt_s(p);
    J.block<3,3>(0,3) = this->jac_euclidean_backward_motion_wrt_rotv(p);
    return J;
}
/* returns the jacobian [3x6] of a backward motion with respect to the transform vector, equal to d(exp(tau) / p)/dtau.
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * exp(tau + Deltatau) / p ~=  exp(tau) / p + d(exp(tau) / p)/dtau |[exp(tau)/p] * Deltatau.
 * Note that the jacobian is evaluated at |[exp(tau)/p]. */

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/* ===== ===== Getters ===== ===== */
/* =============================== */

Eigen::Vector3d trfv::get_T() const {
    Eigen::Vector3d r = this->tail<3>();
    double theta = r.norm();
    rotv rotv(r);

    if (theta < math::constant::SMALL_ROT()) { // if rotation is very small (do not change threshold)
        return dcm(rotv) * this->head<3>(); // Note: That is an accurate expansion!
    }
    else {
        Eigen::Matrix3d rm = tools::skew3(r);
        double theta_sq = theta * theta;
        return (((Eigen::Matrix3d::Identity() - dcm(rotv)()) * rm + tools::itself_by_transpose(r)) / theta_sq) * this->head<3>();
    }
}
/* return translation vector to read */

rotv trfv::get_rotv() const {
    return rotv(this->tail<3>());
}
/* return rotation vector object */

Eigen::Vector3d trfv::get_s() const {
    return this->head<3>();
}
/* return d vector */

Eigen::Vector3d trfv::get_inverse_T() const {
    rotv rv(this->tail<3>());
    return - (rv.inverse() * this->get_T());
}
/* return translation vector of inverse or opposite transformation */

rotv trfv::get_inverse_rotv() const {
    return rotv(- this->tail<3>());
}
/* return rotation vector of inverse or opposite transformation */

void trfv::set(const Eigen::Vector6d& op2) {
    static_cast<Eigen::Vector6d&>(*this) = op2;
}
/* set the vector */

void trfv::set(Eigen::Vector6d&& op2) {
    static_cast<Eigen::Vector6d&>(*this) = op2;
}
/* set the vector */

