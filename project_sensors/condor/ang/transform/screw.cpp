#include "screw.h"
#include "speu_rodrigues.h"
#include "speu_dcm.h"
#include "homogeneous.h"
#include "trfv.h"
#include "dual.h"
#include "../rotate/rotv.h"
#include "../tools.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS SCREW
// ===========
// ===========

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

screw::screw(const Eigen::Vector3d& n, const Eigen::Vector3d& m, const double& phi, const double& d,const double& h)
: _n(n / n.norm()), _m(m), _phi(phi * n.norm()), _d(d), _h(h * n.norm()) {
}
/* constructor based on axis line, moment, magnitude, displacement, and pitch */

screw::screw(const rotv& rv, const Eigen::Vector3d& T) {
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = T.norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = T / _phi;
        _d     = T.dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (T.cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(T.cross(_n)));
        _d     = T.dot(_n);
        _h     = _d / _phi;
    }
}
/* constructor based on rotation vector and translation vector */

screw::screw(const speu_rodrigues& gq) {
    rotv rv(gq.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gq.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gq.get_T() / _phi;
        _d     = gq.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gq.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gq.get_T().cross(_n)));
        _d     = gq.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* constructor based on special Euclidean (rodrigues) */

screw::screw(const speu_dcm& gR) {
    rotv rv(gR.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gR.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gR.get_T() / _phi;
        _d     = gR.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gR.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gR.get_T().cross(_n)));
        _d     = gR.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* constructor based on special Euclidean (dcm) */

screw::screw(const homogeneous& M) {
    rotv rv(M.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = M.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = M.get_T() / _phi;
        _d     = M.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (M.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(M.get_T().cross(_n)));
        _d     = M.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* constructor based on homogeneous */

screw::screw(const trfv& tau) {
    double angle_rad = tau.get_rotv()().norm();
    Eigen::Vector3d T = tau.get_T();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = T.norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = T / _phi;
        _d     = T.dot(_n);
    }
    else {
        _n     = tau.get_rotv()() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (T.cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(T.cross(_n)));
        _d     = T.dot(_n);
        _h     = _d / _phi;
    }
}
/* constructor based on transform vector */

screw::screw(const dual& Z) {
    rotv rv(Z.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = Z.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = Z.get_T() / _phi;
        _d     = Z.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _d     = - 2. * Z.get_qd()(0) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _m     = (quat::convert_4dto3d(Z.get_qd()) - 0.5 * _d * Z.get_qr()()(0) * _n) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _h     = _d / _phi;
    }
}
/* constructor based on unit dual quaternion */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

screw::screw(speu_rodrigues&& gq) {
    rotv rv(gq.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gq.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gq.get_T() / _phi;
        _d     = gq.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gq.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gq.get_T().cross(_n)));
        _d     = gq.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* move constructor based on special Euclidean (rodrigues) */

screw::screw(speu_dcm&& gR) {
    rotv rv(gR.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gR.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gR.get_T() / _phi;
        _d     = gR.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gR.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gR.get_T().cross(_n)));
        _d     = gR.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* move constructor based on special Euclidean (dcm) */

screw::screw(homogeneous&& M) {
    rotv rv(M.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = M.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = M.get_T() / _phi;
        _d     = M.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (M.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(M.get_T().cross(_n)));
        _d     = M.get_T().dot(_n);
        _h     = _d / _phi;
    }
}
/* move constructor based on homogeneous */

screw::screw(trfv&& tau) {
    double angle_rad = tau.get_rotv()().norm();
    Eigen::Vector3d T = tau.get_T();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = T.norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = T / _phi;
        _d     = T.dot(_n);
    }
    else {
        _n     = tau.get_rotv()() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (T.cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(T.cross(_n)));
        _d     = T.dot(_n);
        _h     = _d / _phi;
    }
}
/* move constructor based on transform vector */

screw::screw(dual&& Z) {
    rotv rv(Z.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = Z.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = Z.get_T() / _phi;
        _d     = Z.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _d     = - 2. * Z.get_qd()(0) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _m     = (quat::convert_4dto3d(Z.get_qd()) - 0.5 * _d * Z.get_qr()()(0) * _n) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _h     = _d / _phi;
    }
}
/* move constructor based on unit dual quaternion */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

screw& screw::operator=(const speu_rodrigues& gq) {
    rotv rv(gq.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gq.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gq.get_T() / _phi;
        _d     = gq.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gq.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gq.get_T().cross(_n)));
        _d     = gq.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* assignment operator based on special Euclidean (rodrigues) */

screw& screw::operator=(const speu_dcm& gR) {
    rotv rv(gR.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gR.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gR.get_T() / _phi;
        _d     = gR.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gR.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gR.get_T().cross(_n)));
        _d     = gR.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* assignment operator based on special Euclidean (dcm) */

screw& screw::operator=(const homogeneous& M) {
    rotv rv(M.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = M.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = M.get_T() / _phi;
        _d     = M.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (M.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(M.get_T().cross(_n)));
        _d     = M.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* assignment operator based on homogeneous */

screw& screw::operator=(const trfv& tau) {
    double angle_rad = tau.get_rotv()().norm();
    Eigen::Vector3d T = tau.get_T();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = T.norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = T / _phi;
        _d     = T.dot(_n);
    }
    else {
        _n     = tau.get_rotv()() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (T.cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(T.cross(_n)));
        _d     = T.dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* assignment operator based on transform vector */

screw& screw::operator=(const dual& Z) {
    rotv rv(Z.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = Z.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = Z.get_T() / _phi;
        _d     = Z.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _d     = - 2. * Z.get_qd()(0) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _m     = (quat::convert_4dto3d(Z.get_qd()) - 0.5 * _d * Z.get_qr()()(0) * _n) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _h     = _d / _phi;
    }
    return *this;
}
/* assignment operators = based on unit dual quaternion */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

screw& screw::operator=(speu_rodrigues&& gq) {
    rotv rv(gq.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gq.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gq.get_T() / _phi;
        _d     = gq.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gq.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gq.get_T().cross(_n)));
        _d     = gq.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* move assignment operator based on special Euclidean (rodrigues) */

screw& screw::operator=(speu_dcm&& gR) {
    rotv rv(gR.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = gR.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = gR.get_T() / _phi;
        _d     = gR.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (gR.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(gR.get_T().cross(_n)));
        _d     = gR.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* move assignment operatoro based on special Euclidean (dcm) */

screw& screw::operator=(homogeneous&& M) {
    rotv rv(M.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = M.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = M.get_T() / _phi;
        _d     = M.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (M.get_T().cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(M.get_T().cross(_n)));
        _d     = M.get_T().dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* move assignment operator based on homogeneous */

screw& screw::operator=(trfv&& tau) {
    double angle_rad = tau.get_rotv()().norm();
    Eigen::Vector3d T = tau.get_T();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = T.norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = T / _phi;
        _d     = T.dot(_n);
    }
    else {
        _n     = tau.get_rotv()() / angle_rad;
        _phi   = angle_rad;
        _m     = 0.5 * (T.cross(_n) + 1 / std::tan(0.5 * angle_rad) * _n.cross(T.cross(_n)));
        _d     = T.dot(_n);
        _h     = _d / _phi;
    }
    return *this;
}
/* move assignment operator based on transform vector */

screw& screw::operator=(dual&& Z) {
    rotv rv(Z.get_rotv());
    double angle_rad = rv().norm();
    if (angle_rad < math::constant::EPS()) {
        _h     = std::nan("");
        _phi   = Z.get_T().norm();
        _m     = Eigen::Vector3d::Zero();
        _n     = Z.get_T() / _phi;
        _d     = Z.get_T().dot(_n);
    }
    else {
        _n     = rv() / angle_rad;
        _phi   = angle_rad;
        _d     = - 2. * Z.get_qd()(0) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _m     = (quat::convert_4dto3d(Z.get_qd()) - 0.5 * _d * Z.get_qr()()(0) * _n) /  quat::convert_4dto3d(Z.get_qr()()).norm();
        _h     = _d / _phi;
    }
    return *this;
}
/* move assignment operators = based on unit dual quaternion */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

screw screw::pow(const double& t) const {
    if (t < math::constant::EPS()) {
        screw res(_n, _m, _phi * t, _d * t, _h);
        res._h = std::nan("");
        return res;
    }
    else {
        return screw(_n, _m, _phi * t, _d * t, _h);
    }
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns input fraction (interpolation or extrapolation) of the object. */

screw screw::sclerp(const screw& S0, const screw& S1, const double& t) {
    // I changed this from the commented line to the new one in Jan 23 2020 for commonality with other sclerps
    // I can not understand how this was here before
    ////////////////return screw(S1 * S0.inverse()).pow(t) * S0;
    return S0 * screw(S0.inverse() * S1).pow(t); // WHAT SHOULD BE
}
/* screw linear interpolation, returns S0 for t=0 and S1 for t=1 */

screw screw::plus_right(const trfv& tau) const {
    return (*this) * tau.explog_map_screw();
}
/* right plus operator (input rotation located in local tangent space) */

screw screw::plus_left(const trfv& tau) const {
    return tau.explog_map_screw() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d screw::operator *(const Eigen::Vector3d& vecin) const {
    return dual(*this) * vecin;
}
/* overloaded operator * (forward transformation of point, not vector) */

Eigen::Vector3d screw::operator /(const Eigen::Vector3d& vecin) const {
    return dual(*this) / vecin;
}
/* overloaded operator / (backward transformation of point, not vector) */

Eigen::Vector3d screw::operator^(const Eigen::Vector3d& vecin) const {
    return dual(*this) ^ vecin;
}
/* overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */

Eigen::Vector3d screw::operator&(const Eigen::Vector3d& vecin) const {
    return dual(*this) & vecin;
}
/* overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */

trfv screw::minus_right_trfv(const screw& S) const {
    return (S.inverse() * (*this)).explog_map_trfv();
}
/* minus operator (output transformation located in local tangent space) */

trfv screw::minus_left_trfv(const screw& S) const {
    return ((*this) * S.inverse()).explog_map_trfv();
}
/* left minus operator (output transformation located in global tangent space) */

speu_dcm screw::exp_map_speu_dcm() const {
    return speu_dcm(*this);
}
/* exponential map that returns special euclidean direction cosine matrix */

speu_rodrigues screw::exp_map_speu_rodrigues() const {
    return speu_rodrigues(*this);
}
/* exponential map that returns special euclidean rodrigues parameters */

homogeneous screw::exp_map_homogeneous() const {
    return homogeneous(*this);
}
/* exponential map that returns homogeneous matrix */

trfv screw::explog_map_trfv() const {
    return trfv(*this);
}
/* exponential and logarithmic map (they are the same in this case) that returns the transform vector */

dual screw::exp_map_dual() const {
    return dual(*this);
}
/* exponential map that returns unit dual quaternion */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

Eigen::Matrix6d screw::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;
    ang::dcm R(this->get_rotv());

    res.topLeftCorner<3,3>()     = R();
    res.topRightCorner<3,3>()    = tools::skew3(this->get_T()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R();

    return res;
}
/* returns forward adjoint matrix */

Eigen::Matrix6d screw::adjoint_matrix_backward() const {
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

/* ===== ===== Jacobians ===== ===== */
/* ================================= */

/* ===== ===== Getters ===== ===== */
/* =============================== */

Eigen::Vector3d screw::get_T() const {
    dual Z(*this);
    return Z.get_T();

    // DO NOT DELETE ===== DO NOT DELETE ===== DO NOT DELETE
    // IMPORTANT NOTE: Never obtain directly the translation vector with either of the following formulas.
    // The first only works for significant angles. When they start becoming small, there are differences,
    // and I do not know how to expand it accordingly.
    // The second seems to be for the twist screw, which looks like a completely different thing
    // Eigen::Vector3d T = this->get_p() - sin(_phi) * _n.cross(this->get_p()) - cos(_phi) * this->get_p() + _d * _n;
    // Eigen::Vector3d T = _phi * _m + _d * _n;
    // DO NOT DELETE ===== DO NOT DELETE ===== DO NOT DELETE
}
/* return translation vector to read */

rotv screw::get_rotv() const {
    if (std::isnan(_h)) {
        return rotv(Eigen::Vector3d::Zero());
    }
    else {
        return rotv(_phi * _n);
    }
}
/* return rotation vector object */

Eigen::Vector3d screw::get_inverse_T() const {
    dual Z(*this);
    return Z.get_inverse_T();
}
/* return translation vector of inverse or opposite transformation */

rotv screw::get_inverse_rotv() const {
    if (std::isnan(_h)) {
        return rotv(Eigen::Vector3d::Zero());
    }
    else {
        return rotv(- _phi * _n);
    }
}
/* return rotation vector of inverse or opposite transformation */

Eigen::Vector3d screw::get_p() const {
    return _n.cross(_m);
}
/**< get point belonging to axis closest to origin */




