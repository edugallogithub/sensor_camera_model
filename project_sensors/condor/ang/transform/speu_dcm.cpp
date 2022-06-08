#include "speu_dcm.h"
#include "speu_rodrigues.h"
#include "homogeneous.h"
#include "trfv.h"
#include "screw.h"
#include "dual.h"
#include "se3_tangent.h"
#include "../rotate/so3_tangent.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS SPECIAL EUCLIDEAN (based on DCM)
// ======================================
// ======================================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

speu_dcm::speu_dcm(const speu_rodrigues& gq)
: _R(gq.get_rodrigues()), _T(gq.get_T()) {
}
/* constructor based on special Euclidean (rodrigues) */

speu_dcm::speu_dcm(const homogeneous& M)
: _R(M.get_dcm()), _T(M.get_T()) {
}
/* constructor based on homogeneous */

speu_dcm::speu_dcm(const trfv& tau)
: _R(tau.get_rotv()), _T(tau.get_T()) {
}
/* constructor based on transform vector */

speu_dcm::speu_dcm(const screw& S)
: _R(S.get_rotv()), _T(S.get_T()) {
}
/* constructor based on transform vector */

speu_dcm::speu_dcm(const dual& Z)
: _R(Z.get_qr()), _T(Z.get_T()) {
}
/* constructor based on unit dual quaternion */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

speu_dcm::speu_dcm(speu_rodrigues&& gR)
        : _R(gR.get_rodrigues()), _T(gR.get_T()) {
}
/* move constructor based on special Euclidean (dcm) */

speu_dcm::speu_dcm(homogeneous&& M)
        : _R(M.get_dcm()), _T(M.get_T()) {
}
/* move constructor based on homogeneous */

speu_dcm::speu_dcm(trfv&& tau)
        : _R(tau.get_rotv()), _T(tau.get_T()) {
}
/* move constructor based on transform vector */

speu_dcm::speu_dcm(screw&& S)
        : _R(S.get_rotv()), _T(S.get_T()) {
}
/* constructor based on transform vector */

speu_dcm::speu_dcm(dual&& Z)
: _R(Z.get_qr()), _T(Z.get_T()) {
}
/* constructor based on unit dual quaternion */

/* ===== ===== Assignments ===== ===== */
/* ==================================== */

speu_dcm& speu_dcm::operator=(const speu_rodrigues& gR) {
    _R = gR.get_rodrigues();
    _T = gR.get_T();
    return *this;
}
/* assignment operator = based on special Euclidean (rodrigues) */

speu_dcm& speu_dcm::operator=(const homogeneous& M) {
    _R = M.get_dcm();
    _T = M.get_T();
    return *this;
}
/* assignment operator = based on homogeneous */

speu_dcm& speu_dcm::operator=(const trfv& tau) {
    _R = tau.get_rotv();
    _T = tau.get_T();
    return *this;
}
/* assignment operator = based on transform vector */

speu_dcm& speu_dcm::operator=(const screw& S) {
    _R = S.get_rotv();
    _T = S.get_T();
    return *this;
}
/* assignment operator = based on screw */

speu_dcm& speu_dcm::operator=(const dual& Z) {
    _R = Z.get_qr();
    _T = Z.get_T();
    return *this;
}
/* assignment operator = based on unit dual quaternion */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

speu_dcm& speu_dcm::operator=(speu_rodrigues&& gR) {
    _R = gR.get_rodrigues();
    _T = gR.get_T();
    return *this;
}
/* move assignment operator = based on special Euclidean (rodrigues) */

speu_dcm& speu_dcm::operator=(homogeneous&& M) {
    _R = M.get_dcm();
    _T = M.get_T();
    return *this;
}
/* move assignment operator = based on homogeneous */

speu_dcm& speu_dcm::operator=(trfv&& tau) {
    _R = tau.get_rotv();
    _T = tau.get_T();
    return *this;
}
/* move assignment operator = based on transform vector */

speu_dcm& speu_dcm::operator=(screw&& S) {
    _R = S.get_rotv();
    _T = S.get_T();
    return *this;
}
/* move assignment operator = based on screw */

speu_dcm& speu_dcm::operator=(dual&& Z) {
    _R = Z.get_qr();
    _T = Z.get_T();
    return *this;
}
/* assignment operator = based on unit dual quaternion */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

speu_dcm speu_dcm::operator/(const speu_dcm& op2) const {
    dcm Rinv = _R.inverse();
    return speu_dcm(Rinv * op2._R, Rinv * (op2._T - _T));
}
/* overloaded operator / (backward combination of transformations) */

speu_dcm speu_dcm::pow(const double& t) const {
    return this->log_map_screw().pow(t).exp_map_speu_dcm();
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object screw logarithmic map. */

speu_dcm speu_dcm::sclerp(const speu_dcm& gR0, const speu_dcm& gR1, const double& t) {
    screw delta_screw((screw(gR0.inverse() * gR1)).pow(t));
    return gR0 * speu_dcm(delta_screw);
}
/* screw linear interpolation, returns gR0 for t=0 and gR1 for t=1 */

speu_dcm speu_dcm::plus_right(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return (*this) * tau.exp_map_speu_dcm();
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return (*this) * new_tau.exp_map_speu_dcm();
    }
}
/* right plus operator (input rotation located in local tangent space) */

speu_dcm speu_dcm::plus_right(const screw& S) const {
    return (*this) * S.exp_map_speu_dcm();
}
/* right plus operator (input rotation located in local tangent space) */

speu_dcm speu_dcm::plus_left(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return tau.exp_map_speu_dcm() * (*this);
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return new_tau.exp_map_speu_dcm() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space) */

speu_dcm speu_dcm::plus_left(const screw& S) const {
    return S.exp_map_speu_dcm() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

trfv speu_dcm::minus_right_trfv(const speu_dcm& gR) const {
    return (gR.inverse() * (*this)).log_map_trfv();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output transformation located in local tangent space) */

screw speu_dcm::minus_right_screw(const speu_dcm& M) const {
    return (M.inverse() * (*this)).log_map_screw();
}
/* right minus operator (output transformation located in local tangent space) */

trfv speu_dcm::minus_left_trfv(const speu_dcm& gR) const {
    return ((*this) * gR.inverse()).log_map_trfv();
}
/* left minus operator (output transformation located in global tangent space) */

screw speu_dcm::minus_left_screw(const speu_dcm& gR) const {
    return ((*this) * gR.inverse()).log_map_screw();
}
/* left minus operator (output transformation located in global tangent space) */

trfv speu_dcm::log_map_trfv() const {
    return trfv(*this);
}
/* logarithmic map that returns the transform vector */

screw speu_dcm::log_map_screw() const {
    return screw(*this);
}
/* logarithmic map that returns the screw */

/* ===== ===== Adjoint ===== ===== */
/* =============================== */

se3_tangent_homo speu_dcm::operator|(const se3_tangent_homo& xi_homo) const {
    homogeneous M(*this);
    return se3_tangent_homo(M() * xi_homo() * M.inverse()());
}
/* overloaded operator | (forward adjoint) */

se3_tangent speu_dcm::operator|(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = _R * xi().tail<3>();
    res().head<3>() = _R * xi().head<3>() + _T.cross(res().tail<3>());
    return res;
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix6d speu_dcm::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;

    res.topLeftCorner<3,3>()     = _R();
    res.topRightCorner<3,3>()    = tools::skew3(_T) * _R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = _R();

    return res;
}
/* returns forward adjoint matrix */

se3_tangent_homo speu_dcm::operator%(const se3_tangent_homo& xi_homo) const {
    homogeneous M(*this);
    return se3_tangent_homo(M.inverse()() * xi_homo() * M());
}
/* overloaded operator % (backward adjoint) */

se3_tangent speu_dcm::operator%(const se3_tangent& xi) const {
    se3_tangent res;
    res().tail<3>() = _R / xi().tail<3>();
    res().head<3>() = _R / xi().head<3>() - _R / (_T.cross(xi().tail<3>()));
    return res;
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix6d speu_dcm::adjoint_matrix_backward() const {
    Eigen::Matrix6d res;

    res.topLeftCorner<3,3>()     = _R().transpose();
    res.topRightCorner<3,3>()    = - _R().transpose() * tools::skew3(_T);
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = _R().transpose();

    return res;
}
/* returns backward adjoint matrix */

/* ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
/* ================================================================== */

se3_tangent speu_dcm::dot2xibody(const Eigen::Matrix34d& speu_dcmdot) const {
    se3_tangent res;
    res.set_vi(_R / speu_dcmdot.rightCols<1>());
    res.set_w(_R.dot2omegabody(speu_dcmdot.leftCols<3>()));
    return res;
}
/* obtains the body twist or motion velocity from the special euclidean direct cosine matrix and its time derivative */

Eigen::Matrix34d speu_dcm::xibody2dot(const se3_tangent& xi_body_mrps) const {
    Eigen::Matrix34d res;
    res.leftCols<3>()  = _R.omegabody2dot(xi_body_mrps.get_w());
    res.rightCols<1>() = _R * xi_body_mrps.get_vi();
    return res;
}
/* obtains the special euclidean direct cosine matrix derivative with time based on the homogeneous
transformation and the body twist or motion velocity. */

se3_tangent speu_dcm::dot2xispace(const Eigen::Matrix34d& speu_dcmdot) const {
    se3_tangent res;
    res.set_w(_R.dot2omegaspace(speu_dcmdot.leftCols<3>()));
    res.set_vi(speu_dcmdot.rightCols<1>() - res.get_w()().cross(_T));
    return res;
}
/* obtains the space twist or motion velocity from the special euclidean direct cosine matrix and its time derivative */

Eigen::Matrix34d speu_dcm::xispace2dot(const se3_tangent& xi_space_mrps) const {
    Eigen::Matrix34d res;
    res.leftCols<3>()  = _R.omegaspace2dot(xi_space_mrps.get_w());
    res.rightCols<1>() = xi_space_mrps.get_vi() + xi_space_mrps.get_w()().cross(_T);
    return res;
}
/* obtains the special euclidean direct cosine matrix derivative with time based on the homogeneous
transformation and the space twist or motion velocity. */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

Eigen::Vector12d speu_dcm::wedge(const speu_dcm& Gr) {
    Eigen::Vector12d res;
    res.head<9>() = dcm::wedge(Gr._R);
    res.tail<3>() = Gr._T;
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
 * special euclidean dcm form form and returns them in vector form, first those of the rotation matrix (ordered
 * by row), and last the translation vector. */

Eigen::Vector12d speu_dcm::wedge(const Eigen::Matrix34d& Gr) {
    Eigen::Vector12d res;
    res.head<9>() = dcm::wedge(Gr.leftCols<3>());
    res.tail<3>() = Gr.rightCols<1>();
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
 * special euclidean dcm form form and returns them in vector form, first those of the rotation matrix (ordered
 * by row), and last the translation vector. */

speu_dcm speu_dcm::hat(const Eigen::Vector12d& v) {
    return {dcm::hat(v.head<9>()), v.tail<3>()};
}
/* although the hat operator usually applies to the tangent space, here it takes the twelve components of the
 * vector form, first those of the rotation matrix (ordered by row), and last the translation vector, and returns
 * them in special euclidean dcm form. */

/* ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================================== */

Eigen::Matrix36d speu_dcm::jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    res.leftCols<3>()  = _R();
    res.rightCols<3>() = - _R() * tools::skew3(p);
    return res;
}
/* returns the right jacobian of the forward motion action with respect to the motion, equal to d(gR * p) / dDeltatauB,
 * (gR plus DeltatauB) * p = gR * p + J * DeltatauB */

Eigen::Matrix36d speu_dcm::jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    res.leftCols<3>()  = - Eigen::Matrix3d::Identity();
    //res.rightCols<3>() = _R.inverse()() * tools::skew3(p - _T) * _R(); // this result is identical but more expensive to compute
    res.rightCols<3>() = tools::skew3(_R.inverse()() * (p - _T));
    return res;
}
/* returns the right jacobian of the backward motion action with respect to the motion, equal to d(gR / p) / dDeltatauB,
 * (gR plus DeltatauB) / p = gR / p + J * DeltatauB */

Eigen::Matrix6d speu_dcm::jac_right_log() const {
    return this->log_map_trfv().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(gR)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(gR plus Deltatau) ~= Log(gR) + (J * Deltatau) */

Eigen::Matrix6d speu_dcm::jac_right_plus_wrt_first(const trfv& tau) const {
    return tau.exp_map_speu_dcm().adjoint_matrix_backward();
}
/* returns the right jacobian of the motion right plus operator with respect to the group object (first element),
 * equal to d(gR1 plus tau2) / dDeltatau1B.
 * (gR1 plus Deltatau1B) plus tau2 = (gR1 plus tau2) plus J * Deltatau1B */

Eigen::Matrix6d speu_dcm::jac_right_plus_wrt_second(const trfv& tau) const {
    return tau.jac_right();
}
/* returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
 * equal to d(gR1 plus tau2) / dDeltatau2B.
 * gR1 plus (tau2 + Deltatau2B) = (gR1 plus tau2) plus J * Deltatau2B */

Eigen::Matrix6d speu_dcm::jac_right_minus_wrt_first(const speu_dcm& gR) const {
    return this->minus_right_trfv(gR).jac_right_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the first element,
 * equal to d(gR2 minus gR1) / dDeltatau2B.
 * (gR2 plus Deltatau2B) minus gR1 = (gR2 minus gR1) + J * Deltatau2B */

Eigen::Matrix6d speu_dcm::jac_right_minus_wrt_second(const speu_dcm& gR) const {
    return - this->minus_right_trfv(gR).jac_left_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the second element,
 * equal to d(gR2 minus gR1) / dDeltatau1B.
 * gR2 minus(gR1 plus Deltatau1B) = (gR2 minus gR1) + J * Deltatau1B */

Eigen::Matrix6d speu_dcm::jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    ang::so3_tangent w = xi.get_w();
    res.topLeftCorner<3,3>()     = - tools::skew3(_R() * w()) * _R();
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - _R() * tools::skew3(xi.get_vi()) - tools::skew3(_T) * _R() * tools::skew3(w());
    res.bottomRightCorner<3,3>() = - _R() * tools::skew3(w());
    return res;
}
/* returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdgR | xi) / dDeltatauB,
 * Ad(gR plus DeltatauB) | xi = AdgR | xi + J * DeltatauB */

Eigen::Matrix6d speu_dcm::jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d R = _R();
    Eigen::Matrix3d Rtr = R.transpose();

    res.topLeftCorner<3,3>()     = Rtr * tools::skew3(w()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = tools::skew3(Rtr * (xi.get_vi() - tools::skew3(_T) * w()));
    res.bottomRightCorner<3,3>() = tools::skew3(Rtr * w());

    return res;
}
/* returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdgR % xi) / dDeltatauB,
 * Ad(gR plus DeltatauB) % xi = AdgR % xi + J * DeltatauB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix36d speu_dcm::jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    res.leftCols<3>()  = Eigen::Matrix3d::Identity();
    res.rightCols<3>() = - tools::skew3(_R() * p + _T);
    return res;
}
/* returns the left jacobian of the forward motion action with respect to the motion, equal to d(gR * p) / dDeltatauE,
 * (DeltatauE plus gR) * p = gR * p + J * DeltatauE */

Eigen::Matrix36d speu_dcm::jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d Rtr = _R().transpose();
    res.leftCols<3>()  = - Rtr;
    res.rightCols<3>() = Rtr * tools::skew3(p);
    return res;
}
/* returns the left jacobian of the backward motion action with respect to the motion, equal to d(gR / p) / dDeltatauE,
 * (DeltatauE plus gR) / p = gR / p + J * DeltatauE */

Eigen::Matrix6d speu_dcm::jac_left_log() const {
    return this->log_map_trfv().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(gR)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(Deltatau plus gR) ~= Log(gR) + (J * Deltatau) */

Eigen::Matrix6d speu_dcm::jac_left_plus_wrt_first(const trfv& tau) const {
    return tau.jac_left();
}
/* returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
 * equal to d(tau1 plus gR2) / dDeltatau1E.
 * (tau1 + Deltatau1E) plus gR2 = J * Deltatau1E plus (tau1 plus gR2)  */

Eigen::Matrix6d speu_dcm::jac_left_plus_wrt_second(const trfv& tau) const {
    return tau.exp_map_speu_dcm().adjoint_matrix_forward();
}
/* returns the left jacobian of the motion left plus operator with respect to the group object (second element),
 * equal to d(tau1 plus gR2) / dDeltatau2E.
 * tau1 plus (Deltatau2E plus gR2) = J * Deltatau2E plus (tau1 plus gR2) */

Eigen::Matrix6d speu_dcm::jac_left_minus_wrt_first(const speu_dcm& gR) const {
    return this->minus_left_trfv(gR).jac_left_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the first element,
 * equal to d(gR2 minus gR1) / dDeltatau2E.
 * (Deltatau2E plus gR2) minus gR1 = (gR2 minus gR1) + J * Deltatau2E */

Eigen::Matrix6d speu_dcm::jac_left_minus_wrt_second(const speu_dcm& gR) const {
    return - this->minus_left_trfv(gR).jac_right_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the second element,
 * equal to d(gR2 minus gR1) / dDeltatau1E.
 * gR2 minus (Deltatau1E plus gR1) = (gR2 minus gR1) + J * Deltatau1E */

Eigen::Matrix6d speu_dcm::jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d A = tools::skew3(_R() * xi.get_w()());
    res.topLeftCorner<3,3>()     = - A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - tools::skew3(_R() * xi.get_vi()) - tools::skew3(_T) * A + A * tools::skew3(_T);
    res.bottomRightCorner<3,3>() = - A;
    return res;
}
/* returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(AdgR | xi) / dDeltatauE,
 * Ad(DeltatauE plus gR) | xi = AdgR | xi + J * DeltatauE */

Eigen::Matrix6d speu_dcm::jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = _R().transpose();
    Eigen::Matrix3d A = Rtr * tools::skew3(w());

    res.topLeftCorner<3,3>()     = A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = Rtr * (tools::skew3(xi.get_vi()) - tools::skew3(_T) * tools::skew3(w()));
    res.bottomRightCorner<3,3>() = A;

    return res;
}
/* returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(AdgR % xi) / dDeltatauE,
 * Ad(DeltatauE plus gR) % xi = AdgR % xi + J * DeltatauE */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d speu_dcm::jac_euclidean_forward_motion_wrt_point() const {
    return _R();
}
/* returns the jacobian of the forward motion action with respect to the point, equal to d(gR * p) / dDeltap,
 * gR * (p + Deltap) = gR * p + J * Deltap */

Eigen::Matrix3d speu_dcm::jac_euclidean_backward_motion_wrt_point() const {
    return this->get_inverse_dcm()();
}
/* returns the jacobian of the backward motion action with respect to the point, equal to d(gR / p) / dDeltap,
 * gR / (p + Deltap) = gR / p + J * Deltap */

Eigen::Matrix6d speu_dcm::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_forward();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdgR | xi) / dDeltaxi,
 * AdgR | (xi + Deltaxi) = AdgR | xi + J * Deltaxi */

Eigen::Matrix6d speu_dcm::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdgR % xi) / dDeltaxi,
 * AdgR % (xi + Deltaxi) = AdgR % xi + J * Deltaxi */

Eigen::Matrix312d speu_dcm::jac_euclidean_forward_motion_wrt_speu_dcm(const Eigen::Vector3d& p) const {
    Eigen::Matrix312d res;
    res.block<3,9>(0,0) = _R.jac_euclidean_forward_rotation_wrt_dcm(p);
    res.block<3,3>(0,9) = Eigen::Matrix3d::Identity();
    return res;
}
/* returns the jacobian [3x12] of a forward transformation with respect to the affine rotation matrix, equal to d(Gr * p)/dGr.
 * The forward transformation (unlike the backward one) is linear on the affine rotation matrix, so the following expressions are true:
 * Gr * p = d(Gr * p)/dGr |Gr*p * Gr
 * [Gr * exp(Delta tau)] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[Gr * exp(Delta tau)] - Gr}
 * [exp(Delta tau) * Gr] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[exp(Delta tau) * Gr] - Gr}
 * Note that the jacobian is evaluated at |Gr*p.
 * Note that the increment is {[Gr * exp(Delta tau)] - Gr} or [(exp(Delta tau) * Gr) - Gr], and not exp(Delta tau) or Delta Gr. */

Eigen::Matrix312d speu_dcm::jac_euclidean_backward_motion_wrt_speu_dcm(const Eigen::Vector3d& p) const {
    Eigen::Matrix312d res;
    res.block<3,9>(0,0) = _R.jac_euclidean_backward_rotation_wrt_dcm(p - _T);
    res.block<3,3>(0,9) = - _R.jac_euclidean_backward_rotation_wrt_vector();
    return res;
}
/* returns the jacobian [3x12] of a backward transformation with respect to the affine rotation matrix, equal to d(Gr / p)/dGr.
 * The backward transformation (unlike the forward one) is NOT linear on the affine rotation matrix. The 1st order Taylor
 * approximation is valid for very small affine rotation matrix changes (the two following expressions are the same):
 * [Gr * exp(Delta tau)] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[Gr * exp(Delta tau)] - Gr}
 * [exp(Delta tau) * Gr] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[exp(Delta tau) * Gr] - Gr}
 * Note that the jacobian is evaluated at |Gr/p.
 * Note that the increment is {[Gr * exp(Delta tau)] - Gr} or [(exp(Delta tau) * Gr) - Gr], and not exp(Delta tau) or Delta Gr. */

/* ===== ===== Getters ===== ===== */
rotv speu_dcm::get_rotv() const {
    return rotv(_R);
}
/* return rotation vector object */



















