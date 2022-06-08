#ifndef ANG_SO3_TANGENT
#define ANG_SO3_TANGENT

#include "../ang.h"
#include "../quat.h"
#include "../tools.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * This file contains the "so3_tangent" or angular velocity class, which contains the
 * 3-vector form, the skew-symmetric matrix form, and the pure quaternion form.
*/

namespace ang {

    class so3_tangent_skew;
    class so3_tangent_quat;

// CLASS SO3_TANGENT
// =================
// =================

class ANG_API so3_tangent : private Eigen::Vector3d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not initialized) */
    so3_tangent() = default;
    /**< normal constructor based on angular velocity components */
    so3_tangent(double w1, double w2, double w3) : Eigen::Vector3d(w1, w2, w3) {}
    /**< constructor based on size 3 vector */
    explicit so3_tangent(const Eigen::Vector3d& w) : Eigen::Vector3d(w) {}
    /**< constructor based on skew so3 tangent object */
    explicit so3_tangent(const ang::so3_tangent_skew& wskew);
    /**< constructor based on quat so3 tangent object */
    explicit so3_tangent(const ang::so3_tangent_quat& wquat);
    /**< copy constructor */
    so3_tangent(const so3_tangent&) = default;
    /**< destructor */
    ~so3_tangent() = default;

    /**< reverse angular velocity */
    so3_tangent reverse() const {return ang::so3_tangent(- this->get());}

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    so3_tangent(so3_tangent&&) = default;
    /**< move constructor based on size 3 vector */
    explicit so3_tangent(Eigen::Vector3d&& w) : Eigen::Vector3d(w) {}
    /**< move constructor based on skew so3 tangent object */
    explicit so3_tangent(ang::so3_tangent_skew&& wskew);
    /**< move constructor based on quat so3 tangent object */
    explicit so3_tangent(ang::so3_tangent_quat&& wquat);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    so3_tangent& operator=(const so3_tangent&) = default;
    /**< assignment operator based on size 3 vector */
    so3_tangent& operator=(const Eigen::Vector3d& w);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    so3_tangent& operator=(so3_tangent&&) = default;
    /**< move assignment operator based on size 3 vector */
    so3_tangent& operator=(Eigen::Vector3d&& w);

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator + (addition) */
    so3_tangent operator+(const ang::so3_tangent& o) const {return ang::so3_tangent{this->get() + o()};}
    /**< overloaded operator + (addition) */
    so3_tangent operator+(const Eigen::Vector3d& o) const {return ang::so3_tangent{this->get() + o};}
    /**< overloaded operator - (subtraction) */
    so3_tangent operator-(const ang::so3_tangent& o) const {return ang::so3_tangent{this->get() - o()};}
    /**< overloaded operator - (subtraction) */
    so3_tangent operator-(const Eigen::Vector3d& o) const {return ang::so3_tangent{this->get() - o};}

    /**< ===== ===== Point Velocity ===== ===== */
    /**< returns the velocity [mps] of a point attached to the body frame based on the point coordinates. All three
     * variables (object, input and one output) are viewed in the same frame, either spatial or local, so two possible uses:
     * - returns the space velocity of a point attached to the body frame based on the space angular velocity and the point space coordinates.
     * - returns the body velocity of a point attached to the body frame based on the body angular velocity and the point body coordinates (constant). */
    Eigen::Vector3d point_velocity(const Eigen::Vector3d& p_m) {return this->get().cross(p_m);}

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get size 3 vector */
    Eigen::Vector3d& operator()() {return *this;}
    const Eigen::Vector3d& operator()() const {return *this;}
    /**< get size 3 vector */
    Eigen::Vector3d& get() {return *this;}
    const Eigen::Vector3d& get() const {return *this;}

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< takes the angular velocity in matrix form (skew symmetric) and returns it in vector form. */
    static ang::so3_tangent wedge(const ang::so3_tangent_skew& wskew) {return ang::so3_tangent(wskew);}
    /**< takes the angular velocity in quaternion form (pure) and returns it in vector form. */
    static ang::so3_tangent wedge(const ang::so3_tangent_quat& wquat) {return ang::so3_tangent(wquat);}
    /**< takes the angular velocity in vector form and returns it in matrix form (skew symmetric) */
    static ang::so3_tangent_skew hat_skew(const ang::so3_tangent& w);
    /**< takes the angular velocity in vector form and returns it in quaternion form (pure) */
    static ang::so3_tangent_quat hat_quat(const ang::so3_tangent& w);
}; // closes class so3_tangent

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SO3_TANGENT_SKEW
// ======================
// ======================

class ANG_API so3_tangent_skew : private Eigen::Matrix3d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not implemented) */
    so3_tangent_skew() = delete;
    /**< constructor based on 3x3 skew symmetric matrix (does not verify skew symmetric) */
    explicit so3_tangent_skew(const Eigen::Matrix3d& wskew) : Eigen::Matrix3d(wskew) {}
    /**< constructor based on so3 tangent object */
    explicit so3_tangent_skew(const ang::so3_tangent& w) : Eigen::Matrix3d(ang::tools::skew3(w())) {}
    /**< copy constructor */
    so3_tangent_skew(const so3_tangent_skew&) = default;
    /**< destructor */
    ~so3_tangent_skew() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    so3_tangent_skew(so3_tangent_skew&&) = default;
    /**< move constructor based on 3x3 skew symmetric matrix (does not verify skew symmetric) */
    explicit so3_tangent_skew(Eigen::Matrix3d&& wskew) : Eigen::Matrix3d(wskew) {}
    /**< move constructor based on so3 tangent object */
    explicit so3_tangent_skew(ang::so3_tangent&& w) : Eigen::Matrix3d(ang::tools::skew3(w())) {}

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    so3_tangent_skew& operator=(const so3_tangent_skew&) = default;
    /**< assignment operator based on 3x3 skew symmetric matrix (does not verify skew symmetric) */
    so3_tangent_skew& operator=(const Eigen::Matrix3d& w_skew);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    so3_tangent_skew& operator=(so3_tangent_skew&&) = default;
    /**< move assignment operator based on 3x3 skew symmetric matrix (does not verify skew symmetric) */
    so3_tangent_skew& operator=(Eigen::Matrix3d&& wskew);

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 3x3 matrix */
    Eigen::Matrix3d& operator()() {return *this;}
    const Eigen::Matrix3d& operator()() const {return *this;}
    /**< get 3x3 matrix */
    Eigen::Matrix3d& get() {return *this;}
    const Eigen::Matrix3d& get() const {return *this;}
}; // closes class so3_tangent_skew

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SO3_TANGENT_QUAT
// ======================
// ======================

class ANG_API so3_tangent_quat : private ang::quat {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not implemented) */
    so3_tangent_quat() = delete;
    /**< constructor based on pure quaternion (does not verify) */
    explicit so3_tangent_quat(const ang::quat& wquat) : ang::quat(wquat) {}
    /**< constructor based on so3 tangent object */
    explicit so3_tangent_quat(const ang::so3_tangent& w) : ang::quat(ang::quat::convert_3dto4d(w())) {}
    /**< copy constructor */
    so3_tangent_quat(const so3_tangent_quat&) = default;
    /**< destructor */
    ~so3_tangent_quat() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    so3_tangent_quat(so3_tangent_quat&&) = default;
    /**< move constructor based on pure quaternion (does not verify) */
    explicit so3_tangent_quat(ang::quat&& wquat) : ang::quat(wquat) {}
    /**< move constructor based on so3 tangent object */
    explicit so3_tangent_quat(ang::so3_tangent&& w) : ang::quat(ang::quat::convert_3dto4d(w())) {}

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    so3_tangent_quat& operator=(const so3_tangent_quat&) = default;
    /**< assignment operator based on pure quaternion (does not verify) */
    so3_tangent_quat& operator=(const ang::quat& w_quat);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    so3_tangent_quat& operator=(so3_tangent_quat&&) = default;
    /**< move assignment operator based on pure quaternion (does not verify) */
    so3_tangent_quat& operator=(ang::quat&& wquat);

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get quaternion */
    ang::quat& operator()() {return *this;}
    const ang::quat& operator()() const {return *this;}
    /**< get quaternion */
    ang::quat& get() {return *this;}
    const ang::quat& get() const {return *this;}
}; // closes class so3_tangent_quat

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace ang

#endif














