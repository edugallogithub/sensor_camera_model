#include "so3_tangent.h"
#include "../tools.h"
#include "../quat.h"

using namespace ang;

// CLASS SO3 TANGENT
// =================
// =================

/* ===== ===== Constructors and Assignments ===== ===== */
/* ==================================================== */

so3_tangent::so3_tangent(const so3_tangent_skew& wskew)
: Eigen::Vector3d(tools::skew3_inverse(wskew())) {
}
/* constructor based on skew so3 tangent object */

so3_tangent::so3_tangent(const so3_tangent_quat& wquat)
       : Eigen::Vector3d(quat::convert_4dto3d(wquat())) {
}
/* constructor based on skew so3 tangent object */

so3_tangent::so3_tangent(so3_tangent_skew&& wskew)
: Eigen::Vector3d(tools::skew3_inverse(wskew())) {
}
/* move constructor based on skew so3 tangent object */

so3_tangent::so3_tangent(so3_tangent_quat&& wquat)
        : Eigen::Vector3d(quat::convert_4dto3d(wquat())) {
}
/* move constructor based on quat so3 tangent object */

so3_tangent& so3_tangent::operator=(const Eigen::Vector3d& w) {
    this->get() = w;
    return *this;
}
/* assignment operator based size 3 vector */

so3_tangent& so3_tangent::operator=(Eigen::Vector3d&& w) {
    this->get() = w;
    return *this;
}
/* move assignment operator based size 3 vector */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */

so3_tangent_skew so3_tangent::hat_skew(const so3_tangent& w) {
    return so3_tangent_skew(w);
}
/* takes the angular velocity in vector form and returns it in matrix form (skew symmetric) */

so3_tangent_quat so3_tangent::hat_quat(const so3_tangent& w) {
    return so3_tangent_quat(w);
}
/* takes the angular velocity in vector form and returns it in quaternion form (pure) */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SO3_TANGENT_SKEW
// ======================
// ======================

so3_tangent_skew& so3_tangent_skew::operator=(const Eigen::Matrix3d& wskew) {
    this->get() = wskew;
    return *this;
}
/* assignment operator based on 3x3 skew symmetric matrix (does not verify skew symmetric) */

so3_tangent_skew& so3_tangent_skew::operator=(Eigen::Matrix3d&& wskew) {
    this->get() = wskew;
    return *this;
}
/* move assignment operator based on 3x3 skew symmetric matrix (does not verify skew symmetric) */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SO3_TANGENT_QUAT
// ======================
// ======================

so3_tangent_quat& so3_tangent_quat::operator=(const quat& wquat) {
    this->get() = wquat;
    return *this;
}
/* assignment operator based on pure quaternion (does not verify) */

so3_tangent_quat& so3_tangent_quat::operator=(quat&& wquat) {
    this->get() = wquat;
    return *this;
}
/* move assignment operator based on pure quaternion (does not verify) */

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////












