#include "se3_tangent.h"
#include "../rotate/so3_tangent.h"

using namespace ang;

// CLASS SE3 TANGENT
// =================
// =================

/* ===== ===== Constructors and Assignments ===== ===== */
/* ==================================================== */

se3_tangent::se3_tangent(const Eigen::Vector3d& vi, const so3_tangent& w) {
    this->head<3>() = vi;
    this->tail<3>() = w();
}
/* constructor based on linear velocity vector vi plus angular velocity */

se3_tangent::se3_tangent(const se3_tangent_homo& xihomo) {
    this->tail<3>() = tools::skew3_inverse(xihomo().block<3,3>(0,0));
    this->head<3>() = xihomo().block<3,1>(0,3);
}
/* constructor based on homogeneous se3 tangent object */

se3_tangent::se3_tangent(const se3_tangent_dual& xidual) {
    this->tail<3>() = quat::convert_4dto3d(xidual().get_qr());
    this->head<3>() = quat::convert_4dto3d(xidual().get_qd());
}
/* constructor based on pure dual quaternion se3 tangent object */

se3_tangent::se3_tangent(Eigen::Vector3d&& vi, so3_tangent&& w) {
    this->head<3>() = vi;
    this->tail<3>() = w();
}
/* move constructor based on linear velocity vector vi plus angular velocity */

se3_tangent::se3_tangent(se3_tangent_homo&& xihomo) {
    this->tail<3>() = tools::skew3_inverse(xihomo().block<3,3>(0,0));
    this->head<3>() = xihomo().block<3,1>(0,3);
}
/* move constructor based on homogeneous se3 tangent object */

se3_tangent::se3_tangent(se3_tangent_dual&& xidual) {
    this->tail<3>() = quat::convert_4dto3d(xidual().get_qr());
    this->head<3>() = quat::convert_4dto3d(xidual().get_qd());
}
/* move constructor based on pure dual quaternion se3 tangent object */

se3_tangent& se3_tangent::operator=(const Eigen::Vector6d& xi) {
    this->get() = xi;
    return *this;
}
/* assignment operator based size 6 vector */

se3_tangent& se3_tangent::operator=(Eigen::Vector6d&& xi) {
    this->get() = xi;
    return *this;
}
/* move assignment operator based size 6 vector */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

se3_tangent_homo se3_tangent::hat_homo(const se3_tangent& xi) {
    return se3_tangent_homo(xi);
}
/* takes the twist or motion velocity in vector form and returns it in matrix form (homogeneous) */

se3_tangent_dual se3_tangent::hat_dual(const se3_tangent& xi) {
    return se3_tangent_dual(xi);
}
/* takes the twist or motion velocity in vector form and returns it in pure dual quaternion form. */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SE3_TANGENT_HOMO
// ======================
// ======================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

se3_tangent_homo::se3_tangent_homo(const Eigen::Matrix4d& xihomo)
 : Eigen::Matrix4d(xihomo) {
}
/* constructor based on 4x4 homogeneous matrix (does not verify it) */

se3_tangent_homo::se3_tangent_homo(const se3_tangent& xi) {
    this->block<3,3>(0,0) = tools::skew3(xi().tail<3>());
    this->block<3,1>(0,3) = xi().head<3>();
    this->block<1,4>(3,0) = Eigen::RowVector4d::Zero();
}
/* constructor based on se3 tangent object */

se3_tangent_homo::se3_tangent_homo(se3_tangent&& xi) {
    this->block<3,3>(0,0) = tools::skew3(xi().tail<3>());
    this->block<3,1>(0,3) = xi().head<3>();
    this->block<1,4>(3,0) = Eigen::RowVector4d::Zero();
}
/* move constructor based on se3 tangent object */

se3_tangent_homo& se3_tangent_homo::operator=(const Eigen::Matrix4d& xihomo) {
    this->get() = xihomo;
    return *this;
}
/* assignment operator based on 4x4 homogeneous matrix (does not verify it) */

se3_tangent_homo& se3_tangent_homo::operator=(Eigen::Matrix4d&& xihomo) {
    this->get() = xihomo;
    return *this;
}
/* move assignment operator based on 4x4 homogeneous matrix (does not verify it) */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SE3_TANGENT_DUAL
// ======================
// ======================

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

se3_tangent_dual::se3_tangent_dual(const dual_quat& xidual)
: dual_quat(xidual) {
}
/* constructor based on pure dual quaternion (does not verify it) */

se3_tangent_dual::se3_tangent_dual(const se3_tangent& xi) {
    this->_qr = quat::convert_3dto4d(xi().tail<3>());
    this->_qd = quat::convert_3dto4d(xi().head<3>());
}
/* constructor based on se3 tangent object */

se3_tangent_dual::se3_tangent_dual(se3_tangent&& xi) {
    this->_qr = quat::convert_3dto4d(xi().tail<3>());
    this->_qd = quat::convert_3dto4d(xi().head<3>());
}
/* move constructor based on se3 tangent object */

se3_tangent_dual& se3_tangent_dual::operator=(const dual_quat& xidual) {
    this->get() = xidual;
    return *this;
}
/* assignment operator based on pure dual quaternion (does not verify it) */

se3_tangent_dual& se3_tangent_dual::operator=(dual_quat&& xidual) {
    this->get() = xidual;
    return *this;
}
/* move assignment operator based on pure dual quaternion (does not verify it) */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////











