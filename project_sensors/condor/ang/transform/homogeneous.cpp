#include "homogeneous.h"
#include "speu_rodrigues.h"
#include "speu_dcm.h"
#include "trfv.h"
#include "screw.h"
#include "dual.h"
#include "se3_tangent.h"
#include "math/logic/constant.h"

using namespace ang;

// CLASS HOMOGENEOUS
// =================
// =================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

homogeneous::homogeneous(const Eigen::Matrix4d& Omat)
: Eigen::Matrix4d(Omat) {
    this->normalize();
}
/* constructor based on 4x4 matrix */

homogeneous::homogeneous(const dcm& R, const Eigen::Vector3d& T) {
    this->topLeftCorner(3,3) = R();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = T;
    (*this)()(3,3) = 1;
}
/* constructor based on direction cosine matrix and translation vector */

homogeneous::homogeneous(const speu_rodrigues& gq) {
    this->topLeftCorner(3,3) = dcm(gq.get_rodrigues())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = gq.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on special Euclidean (rodrigues) */

homogeneous::homogeneous(const speu_dcm& gR) {
    this->topLeftCorner(3,3) = gR.get_dcm()();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = gR.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on special Euclidean (dcm) */

homogeneous::homogeneous(const trfv& tau) {
    this->topLeftCorner(3,3) = dcm(tau.get_rotv())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = tau.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on transform vector */

homogeneous::homogeneous(const screw& S) {
    this->topLeftCorner(3,3) = dcm(S.get_rotv())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = S.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on transform vector */

homogeneous::homogeneous(const dual& Z) {
    this->topLeftCorner(3,3) = dcm(Z.get_qr())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = Z.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on unit dual quaternion */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

homogeneous::homogeneous(Eigen::Matrix4d&& Omat)
: Eigen::Matrix4d(Omat) {
    this->normalize();
}
/* move constructor based on 4x4 matrix */

homogeneous::homogeneous(dcm&& R, Eigen::Vector3d&& tr) {
    this->topLeftCorner(3,3) = R();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = tr;
    (*this)()(3,3) = 1;
}
/* move constructor based on direction cosine matrix and translation vector */

homogeneous::homogeneous(speu_rodrigues&& gq) {
    this->topLeftCorner(3,3) = dcm(gq.get_rodrigues())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = gq.get_T();
    (*this)()(3,3) = 1;
}
/* move constructor based on special Euclidean (rodrigues) */

homogeneous::homogeneous(speu_dcm&& gR) {
    this->topLeftCorner(3,3) = gR.get_dcm()();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = gR.get_T();
    (*this)()(3,3) = 1;
}
/* move constructor based on special Euclidean (dcm) */

homogeneous::homogeneous(trfv&& tau) {
    this->topLeftCorner(3,3) = dcm(tau.get_rotv())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = tau.get_T();
    (*this)()(3,3) = 1;
}
/* move constructor based on transform vector */

homogeneous::homogeneous(screw&& S) {
    this->topLeftCorner(3,3) = dcm(S.get_rotv())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = S.get_T();
    (*this)()(3,3) = 1;
}
/* constructor based on transform vector */

homogeneous::homogeneous(dual&& Z) {
    this->topLeftCorner(3,3) = dcm(Z.get_qr())();
    this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    this->topRightCorner(3,1) = Z.get_T();
    (*this)()(3,3) = 1;
}
/* move constructor based on unit dual quaternion */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

homogeneous& homogeneous::operator=(const Eigen::Matrix4d& op2) {
    static_cast<Eigen::Matrix4d&>(*this) = op2;
    this->normalize();
    return *this;
}
/* assignment operator based on 4x4 matrix */

homogeneous& homogeneous::operator=(const speu_rodrigues& gq) {
    this->topLeftCorner(3, 3) = dcm(gq.get_rodrigues())();
    this->topRightCorner(3, 1) = gq.get_T();
    return *this;
}
/* assignment operator = based on special Euclidean (rodrigues) */

homogeneous& homogeneous::operator=(const speu_dcm& gR) {
    this->topLeftCorner(3, 3) = gR.get_dcm()();
    this->topRightCorner(3, 1) = gR.get_T();
    return *this;
}
/* assignment operator = based on special Euclidean (dcm) */

homogeneous& homogeneous::operator=(const trfv& tau) {
    this->topLeftCorner(3, 3) = dcm(tau.get_rotv())();
    this->topRightCorner(3, 1) = tau.get_T();
    return *this;
}
/* assignment operator = based on transform vector */

homogeneous& homogeneous::operator=(const screw& S) {
    this->topLeftCorner(3,3) = dcm(S.get_rotv())();
    this->topRightCorner(3,1) = S.get_T();
    return *this;
}
/* assignment operator = based on screw */

homogeneous& homogeneous::operator=(const dual& Z) {
    this->topLeftCorner(3,3) = dcm(Z.get_qr())();
    this->topRightCorner(3,1) = Z.get_T();
    return *this;
}
/* assignment operator = based on unit dual quaternion */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */
homogeneous& homogeneous::operator=(Eigen::Matrix4d&& op2) {
    static_cast<Eigen::Matrix4d&>(*this) = op2;
    this->normalize();
    return *this;
}
/* move assignment operator based on 4x4 matrix */

homogeneous& homogeneous::operator=(speu_rodrigues&& gq) {
    this->topLeftCorner(3, 3) = dcm(gq.get_rodrigues())();
    this->topRightCorner(3, 1) = gq.get_T();
    return *this;
}
/* move assignment operator = based on special Euclidean (rodrigues) */

homogeneous& homogeneous::operator=(speu_dcm&& gR) {
    this->topLeftCorner(3, 3) = gR.get_dcm()();
    this->topRightCorner(3, 1) = gR.get_T();
    return *this;
}
/* move assignment operator = based on special Euclidean (dcm) */

homogeneous& homogeneous::operator=(trfv&& tau) {
    this->topLeftCorner(3, 3) = dcm(tau.get_rotv())();
    this->topRightCorner(3, 1) = tau.get_T();
    return *this;
}
/* move assignment operator = based on transform vector */

homogeneous& homogeneous::operator=(screw&& S) {
    this->topLeftCorner(3,3) = dcm(S.get_rotv())();
    this->topRightCorner(3,1) = S.get_T();
    return *this;
}
/* move assignment operator = based on screw */

homogeneous& homogeneous::operator=(dual&& Z) {
    this->topLeftCorner(3,3) = dcm(Z.get_qr())();
    this->topRightCorner(3,1) = Z.get_T();
    return *this;
}
/* mov assignment operator = based on unit dual quaternion */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */

homogeneous homogeneous::operator*(const homogeneous& op2) const {
    homogeneous res;
    static_cast<Eigen::Matrix4d&>(res) = static_cast<const Eigen::Matrix4d&>(*this) * static_cast<const Eigen::Matrix4d&>(op2);
    // NOTE: THERE SHOULD BE AN ORTHONORMALIZATION METHOD HERE FOR THE ROTATION PART
    return res;
}
/* overloaded operator * (combination of transformations) */

homogeneous homogeneous::operator/(const homogeneous& op2) const {
    homogeneous res;
    static_cast<Eigen::Matrix4d&>(res) = static_cast<const Eigen::Matrix4d&>(this->inverse()) * static_cast<const Eigen::Matrix4d&>(op2);
    // NOTE: THERE SHOULD BE AN ORTHONORMALIZATION METHOD HERE FOR THE ROTATION PART
    return res;
}
/* overloaded operator / (backward combination of transformations) */

homogeneous homogeneous::inverse() const {
    homogeneous res;
    res.topLeftCorner(3,3)    = this->topLeftCorner(3,3).transpose();
    res.topRightCorner(3,1)   = - this->topLeftCorner(3,3).transpose() * this->topRightCorner(3,1);
    res.bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
    res()(3,3)                  = 1;
    return res;
}
/* returns inverse or opposite transformation */

homogeneous homogeneous::pow(const double& t) const {
    return this->log_map_screw().pow(t).exp_map_homogeneous();
}
/* executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object screw logarithmic map. */

homogeneous homogeneous::sclerp(const homogeneous& M0, const homogeneous& M1, const double& t) {
    screw delta_screw((screw(M0.inverse() * M1)).pow(t));
    return M0 * homogeneous(delta_screw);
}
/* screw linear interpolation, returns M0 for t=0 and M1 for t=1 */

homogeneous homogeneous::plus_right(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return (*this) * tau.exp_map_homogeneous();
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return (*this) * new_tau.exp_map_homogeneous();
    }
}
/* right plus operator (input rotation located in local tangent space) */

homogeneous homogeneous::plus_right(const screw& S) const {
    return (*this) * S.exp_map_homogeneous();
}
/* right plus operator (input rotation located in local tangent space) */

homogeneous homogeneous::plus_left(const trfv& tau) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (tau.get_rotv()().norm() < math::constant::PI()) {
        return tau.exp_map_homogeneous() * (*this);
    }
    else {
        double old_norm = tau.get_rotv()().norm();
        double new_norm = math::constant::PI() * 2. - old_norm;
        //rotv new_rv(tau.get_rotv() / old_norm * new_norm * (-1));
        trfv new_tau(rotv(tau.get_rotv() / old_norm * new_norm * (-1)), tau.get_T());
        return new_tau.exp_map_homogeneous() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space) */

homogeneous homogeneous::plus_left(const screw& S) const {
    return S.exp_map_homogeneous() * (*this);
}
/* left plus operator (input rotation located in global tangent space) */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d homogeneous::operator*(const Eigen::Vector3d& vecin) const {
    Eigen::Vector4d qin;
    qin << vecin , 1.0;
    return (static_cast<const Eigen::Matrix4d&>(*this) * qin).head<3>();
}
/* overloaded operator * (forward transformation of point, not vector) */

Eigen::Vector3d homogeneous::operator/(const Eigen::Vector3d& vecin) const {
    Eigen::Vector4d qin;
    qin << vecin , 1.0;
    return (static_cast<const Eigen::Matrix4d&>(this->inverse()) * qin).head<3>();
}
/* overloaded operator / (backward transformation of point, not vector) */

trfv homogeneous::minus_right_trfv(const homogeneous& M) const {
    return (M.inverse() * (*this)).log_map_trfv();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output transformation located in local tangent space) */

screw homogeneous::minus_right_screw(const homogeneous& M) const {
    return (M.inverse() * (*this)).log_map_screw();
}
/* right minus operator (output transformation located in local tangent space) */

trfv homogeneous::minus_left_trfv(const homogeneous& M) const {
    return ((*this) * M.inverse()).log_map_trfv();
}
/* left minus operator (output transformation located in global tangent space) */

screw homogeneous::minus_left_screw(const homogeneous& M) const {
    return ((*this) * M.inverse()).log_map_screw();
}
/* left minus operator (output transformation located in global tangent space) */

trfv homogeneous::log_map_trfv() const {
    return trfv(*this);
}
/* logarithmic map that returns the transform vector */

screw homogeneous::log_map_screw() const {
    return screw(*this);
}
/* logarithmic map that returns the screw */

/* ===== ===== Adjoint ===== ===== */
/* =================================== */

se3_tangent_homo homogeneous::operator|(const se3_tangent_homo& xi_homo) const {
    return se3_tangent_homo(this->get() * xi_homo() * this->inverse()());
}
/* overloaded operator | (forward adjoint) */

se3_tangent homogeneous::operator|(const se3_tangent& xi) const {
    return se3_tangent::wedge(se3_tangent_homo(this->get() * se3_tangent::hat_homo(xi)() * this->inverse()()));
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix6d homogeneous::adjoint_matrix_forward() const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d R = this->topLeftCorner<3,3>();

    res.topLeftCorner<3,3>()     = R;
    res.topRightCorner<3,3>()    = tools::skew3(this->topRightCorner<3,1>()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = R;

    return res;
}
/* returns forward adjoint matrix */

se3_tangent_homo homogeneous::operator%(const se3_tangent_homo& xi_homo) const {
    return se3_tangent_homo(this->inverse()() * xi_homo() * this->get());
}
/* overloaded operator % (backward adjoint) */

se3_tangent homogeneous::operator%(const se3_tangent& xi) const {
    return se3_tangent::wedge(se3_tangent_homo(this->inverse()() * se3_tangent::hat_homo(xi)() * this->get()));
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix6d homogeneous::adjoint_matrix_backward() const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d RT = this->topLeftCorner<3,3>().transpose();

    res.topLeftCorner<3,3>()     = RT;
    res.topRightCorner<3,3>()    = - RT * tools::skew3(this->topRightCorner<3,1>());
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.bottomRightCorner<3,3>() = RT;

    return res;
}
/* returns backward adjoint matrix */

/* ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
/* ================================================================== */

se3_tangent homogeneous::dot2xibody(const Eigen::Matrix4d& homogeneousdot) const {
    return se3_tangent::wedge(se3_tangent_homo(this->inverse()() * homogeneousdot));
}
/* obtains the body twist or motion velocity from the homogeneous transformation and its time derivative */

Eigen::Matrix4d homogeneous::xibody2dot(const se3_tangent& xi_body_mrps) const {
    return (this->get() * se3_tangent::hat_homo(xi_body_mrps)());
}
/* obtains the homogeneous transformation derivative with time based on the homogeneous
transformation and the body twist or motion velocity. */

se3_tangent homogeneous::dot2xispace(const Eigen::Matrix4d& homogeneousdot) const {
    return se3_tangent::wedge(se3_tangent_homo(homogeneousdot * this->inverse()()));
}
/* obtains the space twist or motion velocity from the homogeneous transformation and its time derivative */

Eigen::Matrix4d homogeneous::xispace2dot(const se3_tangent& xi_space_mrps) const {
    return (se3_tangent::hat_homo(xi_space_mrps)() * this->get());
}
/* obtains the homogeneous transformation derivative with time based on the homogeneous
transformation and the space twist or motion velocity. */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

Eigen::Vector12d homogeneous::wedge(const homogeneous& M) {
    Eigen::Vector12d res;
    res.segment<3>(0) = M().block<1,3>(0,0).transpose();
    res.segment<3>(3) = M().block<1,3>(1,0).transpose();
    res.segment<3>(6) = M().block<1,3>(2,0).transpose();
    res.segment<3>(9) = M().block<3,1>(0,3);
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
 * homogeneous matrix in matrix form and returns them in vector form, first those of the rotation matrix (ordered
 * by row), and last the translation vector. */

Eigen::Vector12d homogeneous::wedge(const Eigen::Matrix4d& M) {
    Eigen::Vector12d res;
    res.segment<3>(0) = M.block<1,3>(0,0).transpose();
    res.segment<3>(3) = M.block<1,3>(1,0).transpose();
    res.segment<3>(6) = M.block<1,3>(2,0).transpose();
    res.segment<3>(9) = M.block<3,1>(0,3);
    return res;
}
/* although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
 * homogeneous matrix in matrix form and returns them in vector form, first those of the rotation matrix (ordered
 * by row), and last the translation vector. */

homogeneous homogeneous::hat(const Eigen::Vector12d& v) {
    Eigen::Matrix4d M;
    M.block<1,3>(0,0) = v.segment<3>(0);
    M.block<1,3>(1,0) = v.segment<3>(3);
    M.block<1,3>(2,0) = v.segment<3>(6);
    M.block<3,1>(0,3) = v.segment<3>(9);
    M.block<1,4>(3,0) = Eigen::RowVector4d::Zero();
    return homogeneous(M);
}
/* although the hat operator usually applies to the tangent space, here it takes the twelve components of the
 * homogeneous matrix in vector form, first those of the rotation matrix (ordered by row), and last the translation
 * vector. */

/* ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================================== */

Eigen::Matrix36d homogeneous::jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    res.leftCols<3>()  = this->get_dcm()();;
    res.rightCols<3>() = - this->get_dcm()() * tools::skew3(p);
    return res;
}
/* returns the right jacobian of the forward motion action with respect to the motion, equal to d(M * p) / dDeltatauB,
 * (M plus DeltatauB) * p = M * p + J * DeltatauB */

Eigen::Matrix36d homogeneous::jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R   = this->get_dcm()();
    res.leftCols<3>()  = - Eigen::Matrix3d::Identity();
    //res.rightCols<3>() = R.transpose() * tools::skew3(p - this->get_T()) * R; // this result is identical but more expensive to compute
    res.rightCols<3>() = tools::skew3(R.transpose() * (p - this->get_T()));
    return res;
}
/* returns the right jacobian of the backward motion action with respect to the motion, equal to d(M / p) / dDeltatauB,
 * (M plus DeltatauB) / p = M / p + J * DeltatauB */

Eigen::Matrix6d homogeneous::jac_right_log() const {
    return this->log_map_trfv().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(M)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(M plus Deltatau) ~= Log(M) + (J * Deltatau) */

Eigen::Matrix6d homogeneous::jac_right_plus_wrt_first(const trfv& tau) const {
    return tau.exp_map_homogeneous().adjoint_matrix_backward();
}
/* returns the right jacobian of the motion right plus operator with respect to the group object (first element),
 * equal to d(M1 plus tau2) / dDeltatau1B.
 * (M1 plus Deltatau1B) plus tau2 = (M1 plus tau2) plus J * Deltatau1B */

Eigen::Matrix6d homogeneous::jac_right_plus_wrt_second(const trfv& tau) const {
    return tau.jac_right();
}
/* returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
 * equal to d(M1 plus tau2) / dDeltatau2B.
 * M1 plus (tau2 + Deltatau2B) = (M1 plus tau2) plus J * Deltatau2B */

Eigen::Matrix6d homogeneous::jac_right_minus_wrt_first(const homogeneous& M) const {
    return this->minus_right_trfv(M).jac_right_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the first element,
 * equal to d(M2 minus M1) / dDeltatau2B.
 * (M2 plus Deltatau2B) minus M1 = (M2 minus M1) + J * Deltatau2B */

Eigen::Matrix6d homogeneous::jac_right_minus_wrt_second(const homogeneous& M) const {
    return - this->minus_right_trfv(M).jac_left_inv();
}
/* returns the right jacobian of the motion right minus operator with respect to the second element,
 * equal to d(M2 minus M1) / dDeltatau1B.
 * M2 minus(M1 plus Deltatau1B) = (M2 minus M1) + J * Deltatau1B */

Eigen::Matrix6d homogeneous::jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d R = this->get_dcm()();
    Eigen::Vector3d T = this->get_T();
    ang::so3_tangent w = xi.get_w();
    res.topLeftCorner<3,3>()     = - tools::skew3(R * w()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - R * tools::skew3(xi.get_vi()) - tools::skew3(T) * R * tools::skew3(w());
    res.bottomRightCorner<3,3>() = - R * tools::skew3(w());
    return res;
}
/* returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdM | xi) / dDeltatauB,
 * Ad(M plus DeltatauB) | xi = AdM | xi + J * DeltatauB */

Eigen::Matrix6d homogeneous::jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d R = this->get_dcm()();
    Eigen::Matrix3d Rtr = R.transpose();
    Eigen::Vector3d T = this->get_T();

    res.topLeftCorner<3,3>()     = Rtr * tools::skew3(w()) * R;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = tools::skew3(Rtr * (xi.get_vi() - tools::skew3(T) * w()));
    res.bottomRightCorner<3,3>() = tools::skew3(Rtr * w());

    return res;
}
/* returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdM % xi) / dDeltatauB,
 * Ad(M plus DeltatauB) % xi = AdM % xi + J * DeltatauB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix36d homogeneous::jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d R   = this->get_dcm()();
    res.leftCols<3>()  = Eigen::Matrix3d::Identity();
    res.rightCols<3>() = - tools::skew3(R * p + this->get_T());
    return res;
}
/* returns the left jacobian of the forward motion action with respect to the motion, equal to d(M * p) / dDeltatauE,
 * (DeltatauE plus M) * p = M * p + J * DeltatauE */

Eigen::Matrix36d homogeneous::jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const {
    Eigen::Matrix36d res;
    Eigen::Matrix3d Rtr = this->get_dcm()().transpose();
    res.leftCols<3>()  = - Rtr;
    res.rightCols<3>() = Rtr * tools::skew3(p);
    return res;
}
/* returns the left jacobian of the backward motion action with respect to the motion, equal to d(M / p) / dDeltatauE,
 * (DeltatauE plus M) / p = M / p + J * DeltatauE */

Eigen::Matrix6d homogeneous::jac_left_log() const {
    return this->log_map_trfv().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(M)) / dDeltatau, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small transform vector changes:
 * Log(Deltatau plus M) ~= Log(M) + (J * Deltatau) */

Eigen::Matrix6d homogeneous::jac_left_plus_wrt_first(const trfv& tau) const {
    return tau.jac_left();
}
/* returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
 * equal to d(tau1 plus M2) / dDeltatau1E.
 * (tau1 + Deltatau1E) plus M2 = J * Deltatau1E plus (tau1 plus M2)  */

Eigen::Matrix6d homogeneous::jac_left_plus_wrt_second(const trfv& tau) const {
    return tau.exp_map_homogeneous().adjoint_matrix_forward();
}
/* returns the left jacobian of the motion left plus operator with respect to the group object (second element),
 * equal to d(tau1 plus M2) / dDeltatau2E.
 * tau1 plus (Deltatau2E plus M2) = J * Deltatau2E plus (tau1 plus M2) */

Eigen::Matrix6d homogeneous::jac_left_minus_wrt_first(const homogeneous& M) const {
    return this->minus_left_trfv(M).jac_left_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the first element,
 * equal to d(M2 minus M1) / dDeltatau2E.
 * (Deltatau2E plus M2) minus M1 = (M2 minus M1) + J * Deltatau2E */

Eigen::Matrix6d homogeneous::jac_left_minus_wrt_second(const homogeneous& M) const {
    return - this->minus_left_trfv(M).jac_right_inv();
}
/* returns the left jacobian of the motion left minus operator with respect to the second element,
 * equal to d(M2 minus M1) / dDeltatau1E.
 * M2 minus (Deltatau1E plus M1) = (M2 minus M1) + J * Deltatau1E */

Eigen::Matrix6d homogeneous::jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;
    Eigen::Matrix3d A = tools::skew3(this->get_dcm()() * xi.get_w()());
    res.topLeftCorner<3,3>()     = - A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = - tools::skew3(this->get_dcm()() * xi.get_vi()) - tools::skew3(this->get_T()) * A + A * tools::skew3(this->get_T());
    res.bottomRightCorner<3,3>() = - A;
    return res;
}
/* returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(AdM | xi) / dDeltatauE,
 * Ad(DeltatauE plus M) | xi = AdM | xi + J * DeltatauE */

Eigen::Matrix6d homogeneous::jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const {
    Eigen::Matrix6d res;

    ang::so3_tangent w = xi.get_w();
    Eigen::Matrix3d Rtr = this->get_dcm()().transpose();
    Eigen::Matrix3d A = Rtr * tools::skew3(w());

    res.topLeftCorner<3,3>()     = A;
    res.bottomLeftCorner<3,3>()  = Eigen::Matrix3d::Zero();
    res.topRightCorner<3,3>()    = Rtr * (tools::skew3(xi.get_vi()) - tools::skew3(this->get_T()) * tools::skew3(w()));
    res.bottomRightCorner<3,3>() = A;

    return res;
}
/* returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(AdM % xi) / dDeltatauE,
 * Ad(DeltatauE plus M) % xi = AdM % xi + J * DeltatauE */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d homogeneous::jac_euclidean_forward_motion_wrt_point() const {
    return this->get_dcm()();
}
/* returns the jacobian of the forward motion action with respect to the point, equal to d(M * p) / dDeltap,
 * M * (p + Deltap) = M * p + J * Deltap */

Eigen::Matrix3d homogeneous::jac_euclidean_backward_motion_wrt_point() const {
    return this->get_inverse_dcm()();
}
/* returns the jacobian of the backward motion action with respect to the point, equal to d(M / p) / dDeltap,
 * M / (p + Deltap) = M / p + J * Deltap */

Eigen::Matrix6d homogeneous::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_forward();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdM | xi) / dDeltaxi,
 * AdM | (xi + Deltaxi) = AdM | xi + J * Deltaxi */

Eigen::Matrix6d homogeneous::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdM % xi) / dDeltaxi,
 * AdM % (xi + Deltaxi) = AdM % xi + J * Deltaxi */

Eigen::Matrix312d homogeneous::jac_euclidean_forward_motion_wrt_homogeneous(const Eigen::Vector3d& p) const {
    Eigen::Matrix312d res = Eigen::Matrix312d::Zero();
    res.block<1,3>(0,0) = p.transpose();
    res.block<1,3>(1,3) = p.transpose();
    res.block<1,3>(2,6) = p.transpose();
    res.block<3,3>(0,9).diagonal() = Eigen::Vector3d::Ones();
    return res;
}
/* returns the jacobian [3x12] of a forward transformation with respect to the homogeneous matrix, equal to d(M * p)/dM.
 * The forward transformation (unlike the backward one) is linear on the homogeneous matrix, so the following expressions are true:
 * M * p = d(M * p)/dM |M*p * M
 * [M * exp(Delta tau)] * p ~= M * p + d(M * p)/dM |M*p * {[M * exp(Delta tau)] - M}
 * [exp(Delta tau) * M] * p ~= M * p + d(M * p)/dM |M*p * {[exp(Delta tau) * M] - M}
 * Note that the jacobian is evaluated at |M*p.
 * Note that the increment is {[M * exp(Delta tau)] - M} or [(exp(Delta tau) * M) - M], and not exp(Delta tau) or Delta M.
 * Note that I discard the last row of the homogeneous matrix (it does not add anything) to obtain size 3x12. */

Eigen::Matrix312d homogeneous::jac_euclidean_backward_motion_wrt_homogeneous(const Eigen::Vector3d& p) const {
    Eigen::Matrix312d res;
    res.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * (p(0) - (*this)()(0,3));
    res.block<3,3>(0,3) = Eigen::Matrix3d::Identity() * (p(1) - (*this)()(1,3));
    res.block<3,3>(0,6) = Eigen::Matrix3d::Identity() * (p(2) - (*this)()(2,3));
    res.block<3,3>(0,9) = - this->block<3,3>(0,0).transpose();
    return res;
}
/* returns the jacobian [3x12] of a backward transformation with respect to the homogeneous matrix, equal to d(M / p)/dM.
 * The backward transformation (unlike the forward one) is NOT linear on the homogeneous matrix. The 1st order Taylor approximation
 * is valid for very small affine rotation matrix changes (the two following expressions are the same):
 * [M * exp(Delta tau)] / p ~= M / p + d(M / p)/dM |M/p * {[M * exp(Delta tau)] - M}
 * [exp(Delta tau) * M] / p ~= M / p + d(M / p)/dM |M/p * {[exp(Delta tau) * M] - M}
 * Note that the jacobian is evaluated at |M/p.
 * Note that the increment is {[M plus exp(Delta tau)] - M} or [(M plus Delta M) - M], and not exp(Delta tau) or Delta M.
 * Note that I discard the last row of the homogeneous matrix (it does not add anything) to obtain size 3x12. */

/* ===== ===== Other ===== ===== */
/* ============================= */

void homogeneous::normalize() {
    // Important: As I divide by norm, a0, a1, and a2 are double what they should be
    // algorithm should be iterative, but this is more than enough for our purposes
    Eigen::Vector3d a0((*this)()(0,0) + (*this)()(1,1) * (*this)()(2,2) - (*this)()(1,2) * (*this)()(2,1),
                       (*this)()(0,1) + (*this)()(1,2) * (*this)()(2,0) - (*this)()(1,0) * (*this)()(2,2),
                       (*this)()(0,2) + (*this)()(1,0) * (*this)()(2,1) - (*this)()(1,1) * (*this)()(2,0));
    Eigen::Vector3d a1((*this)()(1,0) + (*this)()(2,1) * (*this)()(0,2) - (*this)()(2,2) * (*this)()(0,1),
                       (*this)()(1,1) + (*this)()(2,2) * (*this)()(0,0) - (*this)()(2,0) * (*this)()(0,2),
                       (*this)()(1,2) + (*this)()(2,0) * (*this)()(0,1) - (*this)()(2,1) * (*this)()(0,0));
    Eigen::Vector3d a2((*this)()(2,0) + (*this)()(0,1) * (*this)()(1,2) - (*this)()(0,2) * (*this)()(1,1),
                       (*this)()(2,1) + (*this)()(0,2) * (*this)()(1,0) - (*this)()(0,0) * (*this)()(1,2),
                       (*this)()(2,2) + (*this)()(0,0) * (*this)()(1,1) - (*this)()(0,1) * (*this)()(1,0));
    this->block<1,3>(0,0) = a0 / a0.norm();
    this->block<1,3>(1,0) = a1 / a1.norm();
    this->block<1,3>(2,0) = a2 / a2.norm();

    //Eigen::Vector3d a0 = (this->row(0) + this->row(1).cross(this->row(2))) * 0.5;
    //Eigen::Vector3d a1 = (this->row(1) + this->row(2).cross(this->row(0))) * 0.5;
    //Eigen::Vector3d a2 = (this->row(2) + this->row(0).cross(this->row(1))) * 0.5;
    //this->row(0) = a0 / a0.norm();
    //this->row(1) = a1 / a1.norm();
    //this->row(2) = a2 / a2.norm();
}
/* normalize the rotation matrix block ensuring it is orthonormal */

/* ===== ===== Getters ===== ===== */
rotv homogeneous::get_rotv() const {
    dcm R (this->topLeftCorner(3,3));
    return rotv(R);
}
/* return rotation vector object */






