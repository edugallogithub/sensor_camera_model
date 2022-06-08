#include "dual.h"
#include "speu_rodrigues.h"
#include "speu_dcm.h"
#include "homogeneous.h"
#include "trfv.h"
#include "screw.h"
#include "se3_tangent.h"
#include "../dual_quat.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS UNIT DUAL QUATERNION
// ==========================
// ==========================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

dual::dual(const rodrigues& qr, const Eigen::Vector3d& T)
: _qr(qr) {
    quat qT = quat::convert_3dto4d(T);
    _qd = (qT * 0.5) * _qr();
}
/* constructor based on rodrigues parameters and translation vector */

dual::dual(const rotv& rv, const Eigen::Vector3d& T)
: _qr(rv) {
    quat qT = quat::convert_3dto4d(T);
    _qd = (qT * 0.5) * _qr();
}
/* constructor based on rotation vector and translation vector */

dual::dual(const speu_rodrigues& gq)
: _qr(gq.get_rodrigues()) {
    quat qT = quat::convert_3dto4d(gq.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* constructor based on special Euclidean (rodrigues) */

dual::dual(const speu_dcm& gR)
: _qr(gR.get_dcm()) {
    quat qT = quat::convert_3dto4d(gR.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* constructor based on special Euclidean (dcm) */

dual::dual(const homogeneous& M)
: _qr(M.get_dcm()) {
    quat qT = quat::convert_3dto4d(M.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* constructor based on homogeneous */

dual::dual(const trfv& tau){
    double angle_rad = tau.get_rotv()().norm();
    if (angle_rad < math::constant::EPS()) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., tau.get_T() * 0.5);
    }
    else {
        _qr = rodrigues(tau.get_rotv());
        quat qT = quat::convert_3dto4d(tau.get_T());
        _qd = (qT * 0.5) * _qr();
    }
}
/* constructor based on transform vector */

dual::dual(const screw& S) {
    if (std::isnan(S.get_h())) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., S.get_phi() * S.get_n() * 0.5);
    }
    else {
        _qr = rotv(S.get_n() * S.get_phi());
        double qd0          = - 0.5 * S.get_d() * sin(S.get_phi() * 0.5);
        Eigen::Vector3d qdv = sin(S.get_phi() * 0.5) * S.get_m() + 0.5 * S.get_d() * cos(S.get_phi() * 0.5) * S.get_n();
        _qd = quat(qd0, qdv);
    }
}
/* constructor based on screw */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

dual::dual(rodrigues&& qr, Eigen::Vector3d&& T)
        : _qr(qr) {
    quat qT = quat::convert_3dto4d(T);
    _qd = (qT * 0.5) * _qr();
}
/* move constructor based on rodrigues parameters and translation vector */

dual::dual(speu_rodrigues&& gq)
        : _qr(gq.get_rodrigues()) {
    quat qT = quat::convert_3dto4d(gq.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* move constructor based on special Euclidean (rodrigues) */

dual::dual(speu_dcm&& gR)
        : _qr(gR.get_dcm()) {
    quat qT = quat::convert_3dto4d(gR.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* move constructor based on special Euclidean (dcm) */

dual::dual(homogeneous&& M)
        : _qr(M.get_dcm()) {
    quat qT = quat::convert_3dto4d(M.get_T());
    _qd = (qT * 0.5) * _qr();
}
/* move constructor based on homogeneous */

dual::dual(trfv&& tau){
    double angle_rad = tau.get_rotv()().norm();
    if (angle_rad < math::constant::EPS()) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., tau.get_T() * 0.5);
    }
    else {
        _qr = rodrigues(tau.get_rotv());
        quat qT = quat::convert_3dto4d(tau.get_T());
        _qd = (qT * 0.5) * _qr();
    }
}
/* move constructor based on transform vector */

dual::dual(screw&& S) {
    if (std::isnan(S.get_h())) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., S.get_phi() * S.get_n() * 0.5);
    }
    else {
        _qr = rotv(S.get_n() * S.get_phi());
        double qd0          = - 0.5 * S.get_d() * sin(S.get_phi() * 0.5);
        Eigen::Vector3d qdv = sin(S.get_phi() * 0.5) * S.get_m() + 0.5 * S.get_d() * cos(S.get_phi() * 0.5) * S.get_n();
        _qd = quat(qd0, qdv);
    }
}
/* move constructor based on screw */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

dual& dual::operator=(const speu_rodrigues& gq) {
    _qr = gq.get_rodrigues();
    quat qT = quat::convert_3dto4d(gq.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* assignment operator based on special Euclidean (rodrigues) */

dual& dual::operator=(const speu_dcm& gR) {
    _qr = gR.get_dcm();
    quat qT = quat::convert_3dto4d(gR.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* assignment operator based on special Euclidean (dcm) */

dual& dual::operator=(const homogeneous& M) {
    _qr = M.get_dcm();
    quat qT = quat::convert_3dto4d(M.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* assignment operator based on homogeneous */

dual& dual::operator=(const trfv& tau){
    double angle_rad = tau.get_rotv()().norm();
    if (angle_rad < math::constant::EPS()) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., tau.get_T() * 0.5);
    }
    else {
        _qr = rodrigues(tau.get_rotv());
        quat qT = quat::convert_3dto4d(tau.get_T());
        _qd = (qT * 0.5) * _qr();
    }
    return *this;
}
/* assignment operator based on transform vector */

dual& dual::operator=(const screw& S) {
    if (std::isnan(S.get_h())) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., S.get_phi() * S.get_n() * 0.5);
    }
    else {
        _qr = rotv(S.get_n() * S.get_phi());
        double qd0          = - 0.5 * S.get_d() * sin(S.get_phi() * 0.5);
        Eigen::Vector3d qdv = sin(S.get_phi() * 0.5) * S.get_m() + 0.5 * S.get_d() * cos(S.get_phi() * 0.5) * S.get_n();
        _qd = quat(qd0, qdv);
    }
    return *this;
}
/* assignment operator based on screw */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

dual& dual::operator=(speu_rodrigues&& gq) {
    _qr = gq.get_rodrigues();
    quat qT = quat::convert_3dto4d(gq.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* move assignment operator based on special Euclidean (rodrigues) */

dual& dual::operator=(speu_dcm&& gR) {
    _qr = gR.get_dcm();
    quat qT = quat::convert_3dto4d(gR.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* move assignment operator based on special Euclidean (dcm) */

dual& dual::operator=(homogeneous&& M) {
    _qr = M.get_dcm();
    quat qT = quat::convert_3dto4d(M.get_T());
    _qd = (qT * 0.5) * _qr();
    return *this;
}
/* move assignment operator based on homogeneous */

dual& dual::operator=(trfv&& tau){
    double angle_rad = tau.get_rotv()().norm();
    if (angle_rad < math::constant::EPS()) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., tau.get_T() * 0.5);
    }
    else {
        _qr = rodrigues(tau.get_rotv());
        quat qT = quat::convert_3dto4d(tau.get_T());
        _qd = (qT * 0.5) * _qr();
    }
    return *this;
}
/* move assignment operator based on transform vector */

dual& dual::operator=(screw&& S) {
    if (std::isnan(S.get_h())) {
        _qr = rodrigues(1, 0, 0, 0);
        _qd = quat(0., S.get_phi() * S.get_n() * 0.5);
    }
    else {
        _qr = rotv(S.get_n() * S.get_phi());
        double qd0          = - 0.5 * S.get_d() * sin(S.get_phi() * 0.5);
        Eigen::Vector3d qdv = sin(S.get_phi() * 0.5) * S.get_m() + 0.5 * S.get_d() * cos(S.get_phi() * 0.5) * S.get_n();
        _qd = quat(qd0, qdv);
    }
    return *this;
}
/* move assignment operator based on screw */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

dual dual::operator*(const dual& op2) const {
    rodrigues qr(this->_qr() * op2._qr());
    quat qd = this->_qr() * op2._qd + this->_qd * op2._qr();
    return dual(qr, qd);
}
/* overloaded operator * (combination of transformations) */

dual dual::operator/(const dual& op2) const {
    rodrigues qr(this->_qr().adjoint() * op2._qr());
    quat qd = this->_qr().adjoint() * op2._qd + this->_qd.adjoint() * op2._qr();
    return dual(qr, qd);
}
/* overloaded operator / (backward combination of transformations) */

dual dual::pow(const double& t) const {
    return this->log_map_screw().pow(t).exp_map_dual();
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object screw logarithmic map. */

dual dual::sclerp(const dual& tau0, const dual& tau1, const double& t) {
    screw delta_screw((screw(tau0.inverse() * tau1)).pow(t));
    return tau0 * dual(delta_screw);
}
/* screw linear interpolation, returns Z0 for t=0 and Z1 for t=1 */

dual dual::plus_right(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return (*this) * tau.exp_map_dual();
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return (*this) * new_tau.exp_map_dual();
    }
}
/* right plus operator (input rotation located in local tangent space) */

dual dual::plus_right(const screw& S) const {
    return (*this) * S.exp_map_dual();
}
/* right plus operator (input rotation located in local tangent space) */

dual dual::plus_left(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return tau.exp_map_dual() * (*this);
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return new_tau.exp_map_dual() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space) */

dual dual::plus_left(const screw& S) const {
    return S.exp_map_dual() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

trfv dual::minus_right_trfv(const dual& Z) const {
    return (Z.inverse() * (*this)).log_map_trfv();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output transformation located in local tangent space) */

screw dual::minus_right_screw(const dual& Z) const {
    return (Z.inverse() * (*this)).log_map_screw();
}
/* right minus operator (output transformation located in local tangent space) */

trfv dual::minus_left_trfv(const dual& Z) const {
    return ((*this) * Z.inverse()).log_map_trfv();
}
/* left minus operator (output transformation located in global tangent space) */

screw dual::minus_left_screw(const dual& Z) const {
    return ((*this) * Z.inverse()).log_map_screw();
}
/* left minus operator (output transformation located in global tangent space) */

screw dual::log_map_screw() const {
    return screw(*this);
}
/* logarithmic map that returns the screw */

trfv dual::log_map_trfv() const {
    return trfv(*this);
}
/* logarithmic map that returns the transform vector */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

se3_tangent_dual dual::operator|(const se3_tangent_dual& xi_dual) const {
    quat qtemp_r = _qr() * xi_dual().get_qr();
    quat qtemp_d = _qr() * xi_dual().get_qd() + _qd * xi_dual().get_qr();
    return se3_tangent_dual(dual_quat(qtemp_r * _qr().adjoint(), qtemp_r * _qd.adjoint() + qtemp_d * _qr().adjoint()));
}
/* overloaded operator | (forward adjoint) */

se3_tangent dual::operator|(const se3_tangent& xi) const {
    return se3_tangent::wedge(*this | se3_tangent::hat_dual(xi));
    // reason is done this way is that dual is not a dual quaternion, and hence multiplication would not work
    // recoding not worth it in this case, at least for the time being
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix6d dual::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;
    dcm R(_qr);

    res.topLeftCorner<3,3>()     = R();
    res.topRightCorner<3,3>()    = tools::skew3(this->get_T()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R();

    return res;
}
/* returns forward adjoint matrix */

se3_tangent_dual dual::operator%(const se3_tangent_dual& xi_dual) const {
    quat qtemp_r = _qr().adjoint() * xi_dual().get_qr();
    quat qtemp_d = _qr().adjoint() * xi_dual().get_qd() + _qd.adjoint() * xi_dual().get_qr();
    return se3_tangent_dual(dual_quat(qtemp_r * _qr(), qtemp_r * _qd + qtemp_d * _qr()));
}
/* overloaded operator % (backward adjoint) */

se3_tangent dual::operator%(const se3_tangent& xi) const {
    return se3_tangent::wedge(*this % se3_tangent::hat_dual(xi));
    // reason is done this way is that dual is not a dual quaternion, and hence multiplication would not work
    // recoding not worth it in this case, at least for the time being
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix6d dual::adjoint_matrix_backward() const {
    Eigen::Matrix6d res;
    dcm R(_qr);

    res.topLeftCorner<3,3>()     = R().transpose();
    res.topRightCorner<3,3>()    = - R().transpose() * tools::skew3(this->get_T());
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R().transpose();

    return res;
}
/* returns backward adjoint matrix */

/* ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
/* ================================================================== */

se3_tangent dual::dot2xibody(const dual_quat& dualdot) const {
    return se3_tangent::wedge(se3_tangent_dual((  dual_quat(_qr().adjoint(), _qd.adjoint()) * dualdot) * 2.0));
}
/* obtains the body twist or motion velocity from the unit dual quaternion transformation and its time derivative */

dual_quat dual::xibody2dot(const se3_tangent& xi_body_mrps) const {
    return (dual_quat(_qr(), _qd) * se3_tangent::hat_dual(xi_body_mrps)()) * 0.5;
}
/* obtains the unit dual quaternion transformation derivative with time based on the unit dual quaternion
transformation and the body twist or motion velocity. */

se3_tangent dual::dot2xispace(const dual_quat& dualdot) const {
    return se3_tangent::wedge(se3_tangent_dual((dualdot * dual_quat(_qr().adjoint(), _qd.adjoint())) * 2.0));
}
/* obtains the space twist or motion velocity from the unit dual quaternion transformation and its time derivative */

dual_quat dual::xispace2dot(const se3_tangent& xi_space_mrps) const {
    return (se3_tangent::hat_dual(xi_space_mrps)() * dual_quat(_qr(), _qd)) * 0.5;
}
/* obtains the unit dual quaternions transformation derivative with time based on the unit dual quaternion
transformation and the space twist or motion velocity. */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

Eigen::Vector8d dual::wedge(const dual& zeta) {
    Eigen::Vector8d res;
    res.segment<4>(0) = zeta.get_qr()();
    res.segment<4>(4) = zeta.get_qd();
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the eight components of the
 * unit dual quaternion in dual form and returns them in vector form. */

dual dual::hat(const Eigen::Vector8d& v) {
    dual zeta;
    zeta._qr() = v.segment<4>(0);
    zeta._qd = v.segment<4>(4);
    zeta.normalize();
    return zeta;
}
/* although the hat operator usually applies to the tangent space, here it takes the eight components of the
 * unit dual quaternion in vector form, and returns them in unit dual quaternion form (object). */

/* ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================================== */

Eigen::Matrix36d dual::jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(_qr)();
    res.leftCols<3>()  = R;
    res.rightCols<3>() = - R * tools::skew3(p);
    return res;
}
/* returns the right jacobian of the forward motion action with respect to the motion, equal to d(zeta * p) / dDeltatauB,
 * (zeta plus DeltatauB) * p = zeta * p + J * DeltatauB */

Eigen::Matrix36d dual::jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(_qr)();
    res.leftCols<3>()  = - Eigen::Matrix3d::Identity();
    //res.rightCols<3>() = R.transpose() * tools::skew3(p - this->get_T()) * R; // this result is identical but more expensive to compute
    res.rightCols<3>() = tools::skew3(R.transpose() * (p - this->get_T()));
    return res;
}
/* returns the right jacobian of the backward motion action with respect to the motion, equal to d(zeta / p) / dDeltatauB,
 * (zeta plus DeltatauB) / p = zeta / p + J * DeltatauB */

Eigen::Matrix6d dual::jac_right_log() const {
    return this->log_map_trfv().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(zeta)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(zeta plus Deltatau) ~= Log(zeta) + (J * Deltatau) */

Eigen::Matrix6d dual::jac_right_plus_wrt_first(const trfv& tau) const {
    return tau.exp_map_dual().adjoint_matrix_backward();
}
/* returns the right jacobian of the motion right plus operator with respect to the group object (first element),
 * equal to d(zeta1 plus tau2) / dDeltatau1B.
 * (zeta1 plus Deltatau1B) plus tau2 = (zeta1 plus tau2) plus J * Deltatau1B */

Eigen::Matrix6d dual::jac_right_plus_wrt_second(const trfv& tau) const {
    return tau.jac_right();
}
/* returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
 * equal to d(zeta1 plus tau2) / dDeltatau2B.
 * zeta1 plus (tau2 + Deltatau2B) = (zeta1 plus tau2) plus J * Deltatau2B */

Eigen::Matrix6d dual::jac_right_minus_wrt_first(const dual& Z) const {
    return this->minus_right_trfv(Z).jac_right_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the first element,
 * equal to d(zeta2 minus zeta1) / dDeltatau2B.
 * (zeta2 plus Deltatau2B) minus zeta1 = (zeta2 minus zeta1) + J * Deltatau2B */

Eigen::Matrix6d dual::jac_right_minus_wrt_second(const dual& Z) const {
    return - this->minus_right_trfv(Z).jac_left_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the second element,
 * equal to d(zeta2 minus zeta1) / dDeltatau1B.
 * zeta2 minus(zeta1 plus Deltatau1B) = (zeta2 minus zeta1) + J * Deltatau1B */

Eigen::Matrix6d dual::jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d R = dcm(_qr)();
    ang::so3_tangent w = xi.get_w();
    res.topLeftCorner<3,3>()     = - tools::skew3(R * w()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - R * tools::skew3(xi.get_vi()) - tools::skew3(this->get_T()) * R * tools::skew3(w());
    res.bottomRightCorner<3,3>() = - R * tools::skew3(w());
    return res;
}
/* returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdZ | xi) / dDeltatauB,
 * Ad(Z plus DeltatauB) | xi = AdZ | xi + J * DeltatauB */

Eigen::Matrix6d dual::jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d R = dcm(_qr)();
    Eigen::Matrix3d Rtr = R.transpose();
    res.topLeftCorner<3,3>()     = Rtr * tools::skew3(w()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = tools::skew3(Rtr * (xi.get_vi() - tools::skew3(this->get_T()) * w()));
    res.bottomRightCorner<3,3>() = tools::skew3(Rtr * w());
    return res;
}
/* returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdZ % xi) / dDeltatauB,
 * Ad(Z plus DeltatauB) % xi = AdZ % xi + J * DeltatauB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix36d dual::jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    res.leftCols<3>()  = Eigen::Matrix3d::Identity();
    res.rightCols<3>() = - tools::skew3(dcm(_qr)() * p + this->get_T());
    return res;
}
/* returns the left jacobian of the forward motion action with respect to the motion, equal to d(z * p) / dDeltatauE,
 * (DeltatauE plus z) * p = z * p + J * DeltatauE */

Eigen::Matrix36d dual::jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d Rtr = dcm(_qr)().transpose();
    res.leftCols<3>()  = - Rtr;
    res.rightCols<3>() = Rtr * tools::skew3(p);
    return res;
}
/* returns the left jacobian of the backward motion action with respect to the motion, equal to d(z / p) / dDeltatauE,
 * (DeltatauE plus z) / p = z / p + J * DeltatauE */

Eigen::Matrix6d dual::jac_left_log() const {
    return this->log_map_trfv().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(z)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(Deltatau plus z) ~= Log(z) + (J * Deltatau) */

Eigen::Matrix6d dual::jac_left_plus_wrt_first(const trfv& tau) const {
    return tau.jac_left();
}
/* returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
 * equal to d(tau1 plus z2) / dDeltatau1E.
 * (tau1 + Deltatau1E) plus z2 = J * Deltatau1E plus (tau1 plus z2)  */

Eigen::Matrix6d dual::jac_left_plus_wrt_second(const trfv& tau) const {
    return tau.exp_map_dual().adjoint_matrix_forward();
}
/* returns the left jacobian of the motion left plus operator with respect to the group object (second element),
 * equal to d(tau1 plus z2) / dDeltatau2E.
 * tau1 plus (Deltatau2E plus z2) = J * Deltatau2E plus (tau1 plus z2) */

Eigen::Matrix6d dual::jac_left_minus_wrt_first(const dual& z) const {
    return this->minus_left_trfv(z).jac_left_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the first element,
 * equal to d(z2 minus z1) / dDeltatau2E.
 * (Deltatau2E plus z2) minus z1 = (z2 minus z1) + J * Deltatau2E */

Eigen::Matrix6d dual::jac_left_minus_wrt_second(const dual& z) const {
    return - this->minus_left_trfv(z).jac_right_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the second element,
 * equal to d(z2 minus z1) / dDeltatau1E.
 * z2 minus (Deltatau1E plus z1) = (z2 minus z1) + J * Deltatau1E */

Eigen::Matrix6d dual::jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d R = dcm(_qr)();
    Eigen::Matrix3d A = tools::skew3(R * xi.get_w()());
    res.topLeftCorner<3,3>()     = - A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - tools::skew3(R * xi.get_vi()) - tools::skew3(this->get_T()) * A + A * tools::skew3(this->get_T());
    res.bottomRightCorner<3,3>() = - A;
    return res;
}
/* returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adz | xi) / dDeltatauE,
 * Ad(DeltatauE plus z) | xi = Adz | xi + J * DeltatauE */

Eigen::Matrix6d dual::jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = dcm(_qr)().transpose();
    Eigen::Matrix3d A = Rtr * tools::skew3(w());
    res.topLeftCorner<3,3>()     = A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = Rtr * (tools::skew3(xi.get_vi()) - tools::skew3(this->get_T()) * tools::skew3(w()));
    res.bottomRightCorner<3,3>() = A;
    return res;
}
/* returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adz % xi) / dDeltatauE,
 * Ad(DeltatauE plus z) % xi = Adz % xi + J * DeltatauE */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d dual::jac_euclidean_forward_motion_wrt_point() const {
    return dcm(_qr)();
}
/* returns the jacobian of the forward motion action with respect to the point, equal to d(zeta * p) / dDeltap,
 * zeta * (p + Deltap) = zeta * p + J * Deltap */

Eigen::Matrix3d dual::jac_euclidean_backward_motion_wrt_point() const {
    return dcm(this->get_inverse_rotv())();
}
/* returns the jacobian of the backward motion action with respect to the point, equal to d(zeta / p) / dDeltap,
 * zeta / (p + Deltap) = zeta / p + J * Deltap */

Eigen::Matrix6d dual::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_forward();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdZ | xi) / dDeltaxi,
 * AdZ | (xi + Deltaxi) = AdZ | xi + J * Deltaxi */

Eigen::Matrix6d dual::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdZ % xi) / dDeltaxi,
 * AdZ % (xi + Deltaxi) = AdZ % xi + J * Deltaxi */

/* ===== ===== Setters and Getters ===== ===== */
/* =========================================== */

Eigen::Vector3d dual::get_T() const {
    return quat::convert_4dto3d(_qd * _qr().adjoint()) * 2.;
}
/* return translation vector to read */

rotv dual::get_rotv() const {
    return rotv(_qr);
}
/* return rotation vector object */

Eigen::Vector3d dual::get_inverse_T() const {
    return - (this->get_inverse_rotv() * this->get_T());
}
/* return translation vector of inverse or opposite transformation */

rotv dual::get_inverse_rotv() const {
    rodrigues qr_adj(_qr().adjoint());
    return rotv(qr_adj);
}
/* return rotation vector of inverse or opposite transformation */

void dual::set_T(const Eigen::Vector3d& T) {
    quat qT = quat::convert_3dto4d(T);
    _qd = (qT * 0.5) * _qr();
}
/* set translation component maintaining the rotation component based on the input translation vector.
 * Internally, this maintains the existing real quaternion (_qr) and updates the dual quaternion (_qd) */

void dual::set_rotv(const ang::rotv& r) {
    quat qThalf = (_qd * _qr().adjoint());
    _qr = r;
    _qd = qThalf * _qr();
}
/* set rotation component maintaining the translation component based on the input rotation vector.
 * Internally, this modifies both the real (_qr) and dual (_qd) quaternions. */

void dual::set_rodrigues(const ang::rodrigues& q) {
    quat qThalf = (_qd * _qr().adjoint());
    _qr = q;
    _qd = qThalf * _qr();
}
/* set rotation component maintaining the translation component based on the input unit quaternion.
 * Internally, this modifies both the real (_qr) and dual (_qd) quaternions. */

void dual::normalize() {
    _qr().normalize();
    _qd = _qd.get() - (_qr() * (_qd.dot(_qr()) / _qr().squaredNorm())).get();
}
/* normalize */





























