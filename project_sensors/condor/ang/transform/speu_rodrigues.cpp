#include "speu_rodrigues.h"
#include "speu_dcm.h"
#include "homogeneous.h"
#include "trfv.h"
#include "screw.h"
#include "dual.h"
#include "se3_tangent.h"
#include "../rotate/so3_tangent.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS SPECIAL EUCLIDEAN (based on RODRIGUES)
// ============================================
// ============================================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

speu_rodrigues::speu_rodrigues(const speu_dcm& gR)
: _q(gR.get_dcm()), _T(gR.get_T()) {
}
/* constructor based on special Euclidean (dcm) */

speu_rodrigues::speu_rodrigues(const homogeneous& M)
: _q(M.get_dcm()), _T(M.get_T()) {
}
/* constructor based on homogeneous */

speu_rodrigues::speu_rodrigues(const trfv& tau)
: _q(tau.get_rotv()), _T(tau.get_T()) {
}
/* constructor based on transform vector */

speu_rodrigues::speu_rodrigues(const screw& S)
: _q(S.get_rotv()), _T(S.get_T()) {
}
/* constructor based on transform vector */

speu_rodrigues::speu_rodrigues(const dual& Z)
: _q(Z.get_qr()), _T(Z.get_T()) {
}
/* constructor based on unit dual quaternion */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

speu_rodrigues::speu_rodrigues(speu_dcm&& gR)
        : _q(gR.get_dcm()), _T(gR.get_T()) {
}
/* move constructor based on special Euclidean (dcm) */

speu_rodrigues::speu_rodrigues(homogeneous&& M)
        : _q(M.get_dcm()), _T(M.get_T()) {
}
/* move constructor based on homogeneous */

speu_rodrigues::speu_rodrigues(trfv&& tau)
        : _q(tau.get_rotv()), _T(tau.get_T()) {
}
/* move constructor based on transform vector */

speu_rodrigues::speu_rodrigues(screw&& S)
        : _q(S.get_rotv()), _T(S.get_T()) {
}
/* constructor based on transform vector */

speu_rodrigues::speu_rodrigues(dual&& Z)
        : _q(Z.get_qr()), _T(Z.get_T()) {
}
/* move constructor based on unit dual quaternion */

/* ===== ===== Assignments ===== ===== */
/* ==================================== */

speu_rodrigues& speu_rodrigues::operator=(const speu_dcm& gR) {
    _q = gR.get_dcm();
    _T = gR.get_T();
    return *this;
}
/* assignment operator = based on special Euclidean (dcm) */

speu_rodrigues& speu_rodrigues::operator=(const homogeneous& M) {
    _q = M.get_dcm();
    _T = M.get_T();
    return *this;
}
/* assignment operator = based on homogeneous */

speu_rodrigues& speu_rodrigues::operator=(const trfv& tau) {
    _q = tau.get_rotv();
    _T = tau.get_T();
    return *this;
}
/* assignment operator = based on transform vector */

speu_rodrigues& speu_rodrigues::operator=(const screw& S) {
    _q = S.get_rotv();
    _T = S.get_T();
    return *this;
}
/* assignment operator = based on screw */

speu_rodrigues& speu_rodrigues::operator=(const dual& Z) {
    _q = Z.get_qr();
    _T = Z.get_T();
    return *this;
}
/* assignment operator = based on unit dual quaternion */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

speu_rodrigues& speu_rodrigues::operator=(speu_dcm&& gR) {
    _q = gR.get_dcm();
    _T = gR.get_T();
    return *this;
}
/* move assignment operator = based on special Euclidean (dcm) */

speu_rodrigues& speu_rodrigues::operator=(homogeneous&& M) {
    _q = M.get_dcm();
    _T = M.get_T();
    return *this;
}
/* move assignment operator = based on homogeneous */

speu_rodrigues& speu_rodrigues::operator=(trfv&& tau) {
    _q = tau.get_rotv();
    _T = tau.get_T();
    return *this;
}
/* move assignment operator = based on transform vector */

speu_rodrigues& speu_rodrigues::operator=(screw&& S) {
    _q = S.get_rotv();
    _T = S.get_T();
    return *this;
}
/* move assignment operator = based on screw */

speu_rodrigues& speu_rodrigues::operator=(dual&& Z) {
    _q = Z.get_qr();
    _T = Z.get_T();
    return *this;
}
/* move assignment operator = based on unit dual quaternion */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

speu_rodrigues speu_rodrigues::operator/(const speu_rodrigues& op2) const {
    rodrigues qinv = _q.inverse();
    return speu_rodrigues(qinv * op2._q, qinv * (op2._T - _T));
}
/* overloaded operator / (backward combination of transformations) */

speu_rodrigues speu_rodrigues::pow(const double& t) const {
    return this->log_map_screw().pow(t).exp_map_speu_rodrigues();
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object screw logarithmic map. */

speu_rodrigues speu_rodrigues::sclerp(const speu_rodrigues& gq0, const speu_rodrigues& gq1, const double& t) {
    /////////////////////////////////////// TODO TODO MOST LIKELY I NEED PROTECTION AGAINST NEGATIVE QUATERNIONS LIKE THE QUATERNION SLERP CASE
    screw delta_screw((screw(gq0.inverse() * gq1)).pow(t));
    return gq0 * speu_rodrigues(delta_screw);
}
/* screw linear interpolation, returns gq0 for t=0 and gq1 for t=1 */

speu_rodrigues speu_rodrigues::plus_right(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return (*this) * tau.exp_map_speu_rodrigues();
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return (*this) * new_tau.exp_map_speu_rodrigues();
    }
}
/* right plus operator (input rotation located in local tangent space) */

speu_rodrigues speu_rodrigues::plus_right(const screw& S) const {
    return (*this) * S.exp_map_speu_rodrigues();
}
/* right plus operator (input rotation located in local tangent space) */

speu_rodrigues speu_rodrigues::plus_left(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return tau.exp_map_speu_rodrigues() * (*this);
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return  new_tau.exp_map_speu_rodrigues() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space) */

speu_rodrigues speu_rodrigues::plus_left(const screw& S) const {
    return S.exp_map_speu_rodrigues() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

trfv speu_rodrigues::minus_right_trfv(const speu_rodrigues& gq) const {
    return (gq.inverse() * (*this)).log_map_trfv();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output transformation located in local tangent space) */

screw speu_rodrigues::minus_right_screw(const speu_rodrigues& M) const {
    return (M.inverse() * (*this)).log_map_screw();
}
/* right minus operator (output transformation located in local tangent space) */

trfv speu_rodrigues::minus_left_trfv(const speu_rodrigues& gq) const {
    return ((*this) * gq.inverse()).log_map_trfv();
}
/* left minus operator (output transformation located in global tangent space) */

screw speu_rodrigues::minus_left_screw(const speu_rodrigues& gq) const {
    return ((*this) * gq.inverse()).log_map_screw();
}
/* left minus operator (output transformation located in global tangent space) */

trfv speu_rodrigues::log_map_trfv() const {
    return trfv(*this);
}
/* logarithmic map that returns the transform vector */

screw speu_rodrigues::log_map_screw() const {
    return screw(*this);
}
/* logarithmic map that returns the screw */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

se3_tangent_homo speu_rodrigues::operator|(const se3_tangent_homo& xi_homo) const {
    homogeneous M(*this);
    return se3_tangent_homo(M() * xi_homo() * M.inverse()());
}
/* overloaded operator | (forward adjoint) */

se3_tangent speu_rodrigues::operator|(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = _q * xi().tail<3>();
    res().head<3>() = _q * xi().head<3>() + _T.cross(res().tail<3>());
    return res;
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix6d speu_rodrigues::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;
    dcm R(_q);

    res.topLeftCorner<3,3>()     = R();
    res.topRightCorner<3,3>()    = tools::skew3(_T) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R();

    return res;
}
/* returns forward adjoint matrix */

se3_tangent_homo speu_rodrigues::operator%(const se3_tangent_homo& xi_homo) const {
    homogeneous M(*this);
    return se3_tangent_homo(M.inverse()() * xi_homo() * M());
}
/* overloaded operator % (backward adjoint) */

se3_tangent speu_rodrigues::operator%(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = _q / xi().tail<3>();
    res().head<3>() = _q / xi().head<3>() - _q / (_T.cross(xi().tail<3>()));
    return res;
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix6d speu_rodrigues::adjoint_matrix_backward() const {
    Eigen::Matrix6d res;
    dcm R(_q);

    res.topLeftCorner<3,3>()     = R().transpose();
    res.topRightCorner<3,3>()    = - R().transpose() * tools::skew3(_T);
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R().transpose();

    return res;
}
/* returns backward adjoint matrix */

/* ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
/* ================================================================== */

se3_tangent speu_rodrigues::dot2xibody(const Eigen::Vector7d& speu_rodriguesdot) const {
    se3_tangent res;
    res.set_vi(_q / speu_rodriguesdot.tail<3>());
    res.set_w(_q.dot2omegabody(speu_rodriguesdot.head<4>()));
    return res;
}
/* obtains the body twist or motion velocity from the special euclidean Rodrigues parameters and its time derivative */

Eigen::Vector7d speu_rodrigues::xibody2dot(const se3_tangent& xi_body_mrps) const {
    Eigen::Vector7d res;
    res.head<4>() = _q.omegabody2dot(xi_body_mrps.get_w());
    res.tail<3>() = _q * xi_body_mrps.get_vi();
    return res;
}
/* obtains the special euclidean Rodrigues parameters derivative with time based on the homogeneous
transformation and the body twist or motion velocity. */

se3_tangent speu_rodrigues::dot2xispace(const Eigen::Vector7d& speu_rodriguesdot) const {
    se3_tangent res;
    res.set_w(_q.dot2omegaspace(speu_rodriguesdot.head<4>()));
    res.set_vi(speu_rodriguesdot.tail<3>() - res.get_w()().cross(_T));
    return res;
}
/* obtains the space twist or motion velocity from the special euclidean Rodrigues parameters and its time derivative */

Eigen::Vector7d speu_rodrigues::xispace2dot(const se3_tangent& xi_space_mrps) const {
    Eigen::Matrix<double,7,1> res;
    res.head<4>() = _q.omegaspace2dot(xi_space_mrps.get_w());
    res.tail<3>() = xi_space_mrps.get_vi() + xi_space_mrps.get_w()().cross(_T);
    return res;
}
/* obtains the special euclidean Rodrigues parameters derivative with time based on the homogeneous
transformation and the space twist or motion velocity. */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

Eigen::Vector7d speu_rodrigues::wedge(const speu_rodrigues& Gq) {
    Eigen::Vector7d res;
    res.head<4>() = Gq._q();
    res.tail<3>() = Gq._T;
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the seven components of the
 * special euclidean rodrigues form and returns them in vector form (quaternion first). */

speu_rodrigues speu_rodrigues::hat(const Eigen::Vector7d& v) {
    return {rodrigues(v.head<4>()), v.tail<3>()};
}
/* although the hat operator usually applies to the tangent space, here it takes the seven components of the
 * special euclidean rodrigues in vector from and returns them in object form. */

/* ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================================== */

Eigen::Matrix36d speu_rodrigues::jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(_q)();
    res.leftCols<3>()  = R;;
    res.rightCols<3>() = - R * tools::skew3(p);
    return res;
}
/* returns the right jacobian of the forward motion action with respect to the motion, equal to d(gq * p) / dDeltatauB,
 * (gq plus DeltatauB) * p = gq * p + J * DeltatauB */

Eigen::Matrix36d speu_rodrigues::jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(_q)();
    res.leftCols<3>()  = - Eigen::Matrix3d::Identity();
    //res.rightCols<3>() = R.transpose() * tools::skew3(p - _T) * R; // this result is identical but more expensive to compute
    res.rightCols<3>() = tools::skew3(R.transpose() * (p - _T));
    return res;
}
/* returns the right jacobian of the backward motion action with respect to the motion, equal to d(gq / p) / dDeltatauB,
 * (gq plus DeltatauB) / p = gq / p + J * DeltatauB */

Eigen::Matrix6d speu_rodrigues::jac_right_log() const {
    return this->log_map_trfv().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(gq)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(gq plus Deltatau) ~= Log(gq) + (J * Deltatau) */

Eigen::Matrix6d speu_rodrigues::jac_right_plus_wrt_first(const trfv& tau) const {
    return tau.exp_map_speu_rodrigues().adjoint_matrix_backward();
}
/* returns the right jacobian of the motion right plus operator with respect to the group object (first element),
 * equal to d(gq1 plus tau2) / dDeltatau1B.
 * (gq1 plus Deltatau1B) plus tau2 = (gq1 plus tau2) plus J * Deltatau1B */

Eigen::Matrix6d speu_rodrigues::jac_right_plus_wrt_second(const trfv& tau) const {
    return tau.jac_right();
}
/* returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
 * equal to d(gq1 plus tau2) / dDeltatau2B.
 * gq1 plus (tau2 + Deltatau2B) = (gq1 plus tau2) plus J * Deltatau2B */

Eigen::Matrix6d speu_rodrigues::jac_right_minus_wrt_first(const speu_rodrigues& gq) const {
    return this->minus_right_trfv(gq).jac_right_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the first element,
 * equal to d(gq2 minus gq1) / dDeltatau2B.
 * (gq2 plus Deltatau2B) minus gq1 = (gq2 minus gq1) + J * Deltatau2B */

Eigen::Matrix6d speu_rodrigues::jac_right_minus_wrt_second(const speu_rodrigues& gq) const {
    return - this->minus_right_trfv(gq).jac_left_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the second element,
 * equal to d(gq2 minus gq1) / dDeltatau1B.
 * gq2 minus(gq1 plus Deltatau1B) = (gq2 minus gq1) + J * Deltatau1B */

Eigen::Matrix6d speu_rodrigues::jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    dcm R(_q);
    ang::so3_tangent w = xi.get_w();
    res.topLeftCorner<3,3>()     = - tools::skew3(R() * w()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - R() * tools::skew3(xi.get_vi()) - tools::skew3(_T) * R() * tools::skew3(w());
    res.bottomRightCorner<3,3>() = - R() * tools::skew3(w());
    return res;
}
/* returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(Adgq | xi) / dDeltatauB,
 * Ad(gq plus DeltatauB) | xi = Adgq | xi + J * DeltatauB */

Eigen::Matrix6d speu_rodrigues::jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    dcm R(_q);
    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = R().transpose();

    res.topLeftCorner<3,3>()     = Rtr * tools::skew3(w()) * R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = tools::skew3(Rtr * (xi.get_vi() - tools::skew3(_T) * w()));
    res.bottomRightCorner<3,3>() = tools::skew3(Rtr * w());

    return res;
}
/* returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(Adgq % xi) / dDeltatauB,
 * Ad(gq plus DeltatauB) % xi = Adgq % xi + J * DeltatauB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix36d speu_rodrigues::jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R = dcm(_q)();
    res.leftCols<3>()  = Eigen::Matrix3d::Identity();
    res.rightCols<3>() = - tools::skew3(R * p + _T);
    return res;
}
/* returns the left jacobian of the forward motion action with respect to the motion, equal to d(gq * p) / dDeltatauE,
 * (DeltatauE plus gq) * p = gq * p + J * DeltatauE */

Eigen::Matrix36d speu_rodrigues::jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d Rtr = dcm(_q)().transpose();
    res.leftCols<3>()  = - Rtr;
    res.rightCols<3>() = Rtr * tools::skew3(p);
    return res;
}
/* returns the left jacobian of the backward motion action with respect to the motion, equal to d(gq / p) / dDeltatauE,
 * (DeltatauE plus gq) / p = gq / p + J * DeltatauE */

Eigen::Matrix6d speu_rodrigues::jac_left_log() const {
    return this->log_map_trfv().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(gq)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(Deltatau plus gq) ~= Log(gq) + (J * Deltatau) */

Eigen::Matrix6d speu_rodrigues::jac_left_plus_wrt_first(const trfv& tau) const {
    return tau.jac_left();
}
/* returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
 * equal to d(tau1 plus gq2) / dDeltatau1E.
 * (tau1 + Deltatau1E) plus gq2 = J * Deltatau1E plus (tau1 plus gq2)  */

Eigen::Matrix6d speu_rodrigues::jac_left_plus_wrt_second(const trfv& tau) const {
    return tau.exp_map_speu_rodrigues().adjoint_matrix_forward();
}
/* returns the left jacobian of the motion left plus operator with respect to the group object (second element),
 * equal to d(tau1 plus gq2) / dDeltatau2E.
 * tau1 plus (Deltatau2E plus gq2) = J * Deltatau2E plus (tau1 plus gq2) */

Eigen::Matrix6d speu_rodrigues::jac_left_minus_wrt_first(const speu_rodrigues& gq) const {
    return this->minus_left_trfv(gq).jac_left_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the first element,
 * equal to d(gq2 minus gq1) / dDeltatau2E.
 * (Deltatau2E plus gq2) minus gq1 = (gq2 minus gq1) + J * Deltatau2E */

Eigen::Matrix6d speu_rodrigues::jac_left_minus_wrt_second(const speu_rodrigues& gq) const {
    return - this->minus_left_trfv(gq).jac_right_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the second element,
 * equal to d(gq2 minus gq1) / dDeltatau1E.
 * gq2 minus (Deltatau1E plus gq1) = (gq2 minus gq1) + J * Deltatau1E */


Eigen::Matrix6d speu_rodrigues::jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d R = dcm(_q)();
    Eigen::Matrix3d A = tools::skew3(R * xi.get_w()());
    res.topLeftCorner<3,3>()     = - A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - tools::skew3(R * xi.get_vi()) - tools::skew3(_T) * A + A * tools::skew3(_T);
    res.bottomRightCorner<3,3>() = - A;
    return res;
}
/* returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adgq | xi) / dDeltatauE,
 * Ad(DeltatauE plus gq) | xi = Adgq | xi + J * DeltatauE */

Eigen::Matrix6d speu_rodrigues::jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = dcm(_q)().transpose();
    Eigen::Matrix3d A = Rtr * tools::skew3(w());

    res.topLeftCorner<3,3>()     = A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = Rtr * (tools::skew3(xi.get_vi()) - tools::skew3(_T) * tools::skew3(w()));
    res.bottomRightCorner<3,3>() = A;

    return res;
}
/* returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adgq % xi) / dDeltatauE,
 * Ad(DeltatauE plus gq) % xi = Adgq % xi + J * DeltatauE */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d speu_rodrigues::jac_euclidean_forward_motion_wrt_point() const {
    return dcm(_q)();
}
/* returns the jacobian of the forward motion action with respect to the point, equal to d(gq * p) / dDeltap,
 * gq * (p + Deltap) = gq * p + J * Deltap */

Eigen::Matrix3d speu_rodrigues::jac_euclidean_backward_motion_wrt_point() const {
    return dcm(this->get_inverse_rodrigues())();
}
/* returns the jacobian of the backward motion action with respect to the point, equal to d(gq / p) / dDeltap,
 * gq / (p + Deltap) = gq / p + J * Deltap */

Eigen::Matrix6d speu_rodrigues::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_forward();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adgq | xi) / dDeltaxi,
 * Adgq | (xi + Deltaxi) = Adgq | xi + J * Deltaxi */

Eigen::Matrix6d speu_rodrigues::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adgq % xi) / dDeltaxi,
 * Adgq % (xi + Deltaxi) = Adgq % xi + J * Deltaxi */

Eigen::Matrix37d speu_rodrigues::jac_euclidean_forward_motion_wrt_speu_rodrigues(const Eigen::Vector3d& p) const {
    Eigen::Matrix37d res;
    res.block<3,4>(0,0) = _q.jac_euclidean_forward_rotation_wrt_rodrigues(p);
    res.block<3,3>(0,4) = Eigen::Matrix3d::Identity();
    return res;
}
/* returns the jacobian [3x7] of a forward transformation with respect to the affine unit quaternion, equal to d(Gq * p)/dGq.
 * The 1st order Taylor approximation is valid for very small affine unit quaternion changes (the two following expressions are the same):
 * [Gq * exp(Delta tau)] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[Gq * exp(Delta tau)] - Gq}
 * [exp(Delta tau) * Gq] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[exp(Delta tau) * Gq] - Gq}
 * Note that the jacobian is evaluated at |Gq*p.
 * Note that the increment is {[Gq * exp(Delta tau)] - Gq} or [(exp(Delta tau) * Gq) - Gq], and not exp(Delta tau) or Delta Gq. */

Eigen::Matrix37d speu_rodrigues::jac_euclidean_backward_motion_wrt_speu_rodrigues(const Eigen::Vector3d& p) const {
    Eigen::Matrix37d res;
    res.block<3,4>(0,0) = _q.jac_euclidean_backward_rotation_wrt_rodrigues(p - _T);
    res.block<3,3>(0,4) = - _q.jac_euclidean_backward_rotation_wrt_vector();
    return res;
}
/* returns the jacobian [3x7] of a backward transformation with respect to the affine unit quaternion, equal to d(Gq / p)/dGq.
 * The 1st order Taylor approximation is valid for very small affine unit quaternion changes (the two following expressions are the same):
 * [Gq * exp(Delta tau)] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[Gq * exp(Delta tau)] - Gq}
 * [exp(Delta tau) * Gq] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[exp(Delta tau) * Gq] - Gq}
 * Note that the jacobian is evaluated at |Gq/p.
 * Note that the increment is {[Gq * exp(Delta tau)] - Gq} or [(exp(Delta tau) * Gq) - Gq], and not exp(Delta tau) or Delta Gq. */

/* ===== ===== Getters and Setters ===== ===== */
rotv speu_rodrigues::get_rotv() const {
    return rotv(_q);
}
/* return rotation vector object */

void speu_rodrigues::set(const rodrigues& q, const Eigen::Vector3d& T) {
    _q = q;
    _T = T;
}
/* modify the object based on the input rotation and translation */

void speu_rodrigues::set(const rodrigues& q_inv, const Eigen::Vector3d& T_inv, bool) {
    _q = q_inv.inverse();
    _T = _q * (-T_inv);
}
/* modify the object based on the input inverse rotation and inverse translation.
 * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */








