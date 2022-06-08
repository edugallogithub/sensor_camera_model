#ifndef ANG_SE3_TANGENT
#define ANG_SE3_TANGENT

#include "../ang.h"
#include "../dual_quat.h"
#include "../auxiliary.h"
#include "../rotate/so3_tangent.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * This file contains the "se3_tangent", twist, or motion velocity class, which contains the
 * 6-vector form, the 4x4 homogeneous matrix form, and the pure dual quaternion form.
*/

namespace ang {

class se3_tangent_homo;
class se3_tangent_dual;

// CLASS SE3_TANGENT
// =================
// =================

class ANG_API se3_tangent : private Eigen::Vector6d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not initialized) */
    se3_tangent() = default;
    /**< normal constructor based on linear velocity components plus angular velocity components */
    se3_tangent(double vi1, double vi2, double vi3, double w1, double w2, double w3) {(*this) << vi1, vi2, vi3, w1, w2, w3;}
    /**< constructor based on size 6 vector */
    explicit se3_tangent(const Eigen::Vector6d& xi) : Eigen::Vector6d(xi) {}
    /**< constructor based on linear velocity vector vi plus angular velocity */
    se3_tangent(const Eigen::Vector3d& vi, const ang::so3_tangent& w);
    /**< constructor based on homogeneous se3 tangent object */
    explicit se3_tangent(const ang::se3_tangent_homo& xihomo);
    /**< constructor based on pure dual quaternion se3 tangent object */
    explicit se3_tangent(const ang::se3_tangent_dual& xidual);
    /**< copy constructor */
    se3_tangent(const se3_tangent&) = default;
    /**< destructor */
    ~se3_tangent() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    se3_tangent(se3_tangent&&) = default;
    /**< move constructor based on size 6 vector */
    explicit se3_tangent(Eigen::Vector6d&& xi) : Eigen::Vector6d(xi) {}
    /**< move constructor based on linear velocity vector vi plus angular velocity */
    se3_tangent(Eigen::Vector3d&& vi, ang::so3_tangent&& w);
    /**< move constructor based on homogeneous se3 tangent object */
    explicit se3_tangent(ang::se3_tangent_homo&& xihomo);
    /**< move constructor based on pure dual quaternion se3 tangent object */
    explicit se3_tangent(ang::se3_tangent_dual&& xidual);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    se3_tangent& operator=(const se3_tangent&) = default;
    /**< assignment operator based on size 6 vector */
    se3_tangent& operator=(const Eigen::Vector6d& xi);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    se3_tangent& operator=(se3_tangent&&) = default;
    /**< move assignment operator based on size 6 vector */
    se3_tangent& operator=(Eigen::Vector6d&& xi);

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator + (addition) */
    se3_tangent operator+(const ang::se3_tangent& o) const {return ang::se3_tangent{this->get() + o()};}
    /**< overloaded operator + (addition) */
    se3_tangent operator+(const Eigen::Vector6d& o) const {return ang::se3_tangent{this->get() + o};}
    /**< overloaded operator - (subtraction) */
    se3_tangent operator-(const ang::se3_tangent& o) const {return ang::se3_tangent{this->get() - o()};}
    /**< overloaded operator - (subtraction) */
    se3_tangent operator-(const Eigen::Vector6d& o) const {return ang::se3_tangent{this->get() - o};}

    /**< ===== ===== Point Velocity ===== ===== */
    /**< returns the velocity [mps] of a point attached to the body frame based on the point coordinates. All three
     * variables (object, input and one output) are viewed in the same frame, either spatial or local, so two possible uses:
     * - returns the space velocity of a point attached to the body frame based on the space twist and the point space coordinates.
     * - returns the body velocity of a point attached to the body frame based on the body twist and the point body coordinates (constant). */
    Eigen::Vector3d point_velocity(const Eigen::Vector3d& p_m) {return this->head<3>() + this->tail<3>().cross(p_m);}

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get size 6 vector */
    Eigen::Vector6d& operator()() {return *this;}
    const Eigen::Vector6d& operator()() const {return *this;}
    /**< get size 6 vector */
    Eigen::Vector6d& get() {return *this;}
    const Eigen::Vector6d& get() const {return *this;}

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< takes the twist of motion velocity in homogeneous form and returns it in vector form. */
    static ang::se3_tangent wedge(const ang::se3_tangent_homo& xihomo) {return ang::se3_tangent(xihomo);}
    /**< takes the twist of motion velocity in pure dual quaternion form and returns it in vector form. */
    static ang::se3_tangent wedge(const ang::se3_tangent_dual& xidual) {return ang::se3_tangent(xidual);}
    /**< takes the twist or motion velocity in vector form and returns it in homogeneous form. */
    static ang::se3_tangent_homo hat_homo(const ang::se3_tangent& xi);
    /**< takes the twist or motion velocity in vector form and returns it in pure dual quaternion form. */
    static ang::se3_tangent_dual hat_dual(const ang::se3_tangent& xi);

    /**< ===== ===== Getters and Setters ===== ===== */
    /**< get linear velocity */
    Eigen::Vector3d get_vi() const {return this->head<3>();}
    /**< set linear velocity */
    void set_vi(const Eigen::Vector3d& vi) {this->head<3>() = vi;}

    /**< get angular velocity */
    ang::so3_tangent get_w() const {return ang::so3_tangent(this->tail<3>());}
    /**< set angular velocity */
    void set_w(const ang::so3_tangent& w) {this->tail<3>() = w();}
}; // closes class se3_tangent

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SE3_TANGENT_HOMO
// ======================
// ======================

class ANG_API se3_tangent_homo : private Eigen::Matrix4d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not implemented) */
    se3_tangent_homo() = delete;
    /**< constructor based on 4x4 homogeneous matrix (does not verify it) */
    explicit se3_tangent_homo(const Eigen::Matrix4d& xihomo);
    /**< constructor based on se3 tangent object */
    explicit se3_tangent_homo(const ang::se3_tangent& xi);
    /**< copy constructor */
    se3_tangent_homo(const se3_tangent_homo&) = default;
    /**< destructor */
    ~se3_tangent_homo() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    se3_tangent_homo(se3_tangent_homo&&) = default;
    /**< move constructor based on 4x4 homogeneous matrix (does not verify it) */
    explicit se3_tangent_homo(Eigen::Matrix4d&& xihomo) : Eigen::Matrix4d(xihomo) {}
    /**< move constructor based on se3 tangent object */
    explicit se3_tangent_homo(ang::se3_tangent&& xi);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    se3_tangent_homo& operator=(const se3_tangent_homo&) = default;
    /**< assignment operator based on 4x4 homogeneous matrix (does not it) */
    se3_tangent_homo& operator=(const Eigen::Matrix4d& xihomo);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    se3_tangent_homo& operator=(se3_tangent_homo&&) = default;
    /**< move assignment operator based on 4x4 homogeneous matrix (does not verify it) */
    se3_tangent_homo& operator=(Eigen::Matrix4d&& xihomo);

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 4x4 matrix */
    Eigen::Matrix4d& operator()() {return *this;}
    const Eigen::Matrix4d& operator()() const {return *this;}
    /**< get 4x4 matrix */
    Eigen::Matrix4d& get() {return *this;}
    const Eigen::Matrix4d& get() const {return *this;}
}; // closes class se4_tangent_homo

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS SE3_TANGENT_DUAL
// ======================
// ======================

class ANG_API se3_tangent_dual : private ang::dual_quat {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not implemented) */
    se3_tangent_dual() = delete;
    /**< constructor based on pure dual quaternion (does not verify it) */
    explicit se3_tangent_dual(const ang::dual_quat& xidual);
    /**< constructor based on se3 tangent object */
    explicit se3_tangent_dual(const ang::se3_tangent& xi);
    /**< copy constructor */
    se3_tangent_dual(const se3_tangent_dual&) = default;
    /**< destructor */
    ~se3_tangent_dual() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    se3_tangent_dual(se3_tangent_dual&&) = default;
    /**< move constructor based on pure dual quaternion (does not verify it) */
    explicit se3_tangent_dual(ang::dual_quat&& xidual) : ang::dual_quat(xidual) {}
    /**< move constructor based on se3 tangent object */
    explicit se3_tangent_dual(ang::se3_tangent&& xi);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    se3_tangent_dual& operator=(const se3_tangent_dual&) = default;
    /**< assignment operator based on pure dual quaternion (does not verify it) */
    se3_tangent_dual& operator=(const ang::dual_quat& xidual);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    se3_tangent_dual& operator=(se3_tangent_dual&&) = default;
    /**< move assignment operator based on pure dual quaternion (does not verify it) */
    se3_tangent_dual& operator=(ang::dual_quat&& xidual);

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get dual quat */
    ang::dual_quat& operator()() {return *this;}
    const ang::dual_quat& operator()() const {return *this;}
    /**< get dual quat */
    ang::dual_quat& get() {return *this;}
    const ang::dual_quat& get() const {return *this;}
}; // closes class se4_tangent_dual

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace ang

#endif














