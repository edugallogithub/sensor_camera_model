#ifndef ANG_REPRESENTATION_ROTV
#define ANG_REPRESENTATION_ROTV

#include "../ang.h"
#include "../auxiliary.h" // required for the proper behavior of the << operator
#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * This file contains the "rotv" class that models the rotation vector representation of a rotation.
 * By convention, the modulus of the rotation is always between [0,180] deg.
 * NOTE: On Nov 19 2020 I discovered that this class sometimes returned rotations bigger than 180 degrees, which
 * could create very serious problems (errors) when optimizing with the exponentials. I corrected it.
 */

namespace ang {

class euler;
class rodrigues;
class dcm;
class so3_tangent;
class so3_tangent_skew;

// CLASS ROTV
// ==========
// ==========

class ANG_API rotv : private Eigen::Vector3d {
private:
    /**< modify parameters maintaining rotation so it follows shortest path (rotation angle less than 180 [deg]) */
    void shortest_angle();
    /**< replaces the input 3x1 vector (which should come from a rotation vector by employing the () operator or
     * the get() method) by an equivalent vector of angle [2PI - angle] instead of [angle] and direction [-n]
     * instead of [n]. This is completely equivalent, but the rotation angle is now between 180 and 360 [deg],
     * and that is the reason why this method employs Euclidean 3x1 vectors instead of rotation vector objects. */
    static void equivalent_rotv(Eigen::Vector3d& rv);
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
	/**< empty constructor (not initialized) */
	rotv() = default;
    /**< normal constructor based on rotation vector components */
    rotv(double p1, double p2, double p3);
    /**< constructor based on size 3 vector */
    explicit rotv(const Eigen::Vector3d& Ovec);
    /**< constructor based on size 3 vector multiplied by input factor */
    rotv(const Eigen::Vector3d& Ovec, const double& factor);
	/**< constructor based on Euler angles */
	explicit rotv(const euler&);
	/**< constructor based on Rodrigues parameters */
	explicit rotv(const rodrigues&);
	/**< constructor based on direction cosine matrix */
	explicit rotv(const dcm&);
    /**< VERY IMPORTANT constructor that obtains the rotation from v1 to v2, both normalized.
    The direction is orthogonal to the plane formed by v1 and v2. The magnitude is the angle
    between both input vectors in that plane. The modulus of the input vectors does not matter.
    The result is such that v2.normalized() = this * v1.normalized(). */
    rotv(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);
    /**< copy constructor */
	rotv(const rotv&) = default;
	/**< destructor */
	~rotv() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    rotv(rotv&&) = default;
    /**< move constructor based on size 3 vector */
    explicit rotv(Eigen::Vector3d&& Ovec);
    /**< move constructor based on size 3 vector multiplied by input factor */
    rotv(const Eigen::Vector3d&& Ovec, const double& factor);
    /**< move constructor based on Euler angles */
    explicit rotv(euler&&);
    /**< move constructor based on Rodrigues parameters */
    explicit rotv(rodrigues&&);
    /**< move constructor based on direction cosine matrix */
    explicit rotv(dcm&&);

    /**< ===== ===== Assignments ===== ===== */
	/**< copy assignment */
	rotv& operator=(const rotv&) = default;
    /**< assignment operator based on size 3 vector */
    rotv& operator=(const Eigen::Vector3d&);
    /**< assignment operator based on Euler angles */
    rotv& operator=(const euler&);
    /**< assignment operator based on Rodrigues parameters */
    rotv& operator=(const rodrigues&);
    /**< assignment operator based on direction cosine matrix */
    rotv& operator=(const dcm&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    rotv& operator=(rotv&&) = default;
    /**< move assignment operator based on size 3 vector */
    rotv& operator=(Eigen::Vector3d&&);
    /**< move assignment operator based on Euler angles */
    rotv& operator=(euler&&);
    /**< move assignment operator based on Rodrigues parameters */
    rotv& operator=(rodrigues&&);
    /**< move assignment operator based on direction cosine matrix */
    rotv& operator=(dcm&&);

    /**< ===== ===== Transformations ===== ===== */
	/**< overloaded operator * (combination of rotations) */
	rotv operator*(const rotv& op2) const;
    /**< overloaded operator / (backward combination of rotations) */
    rotv operator/(const rotv& op2) const;
    /**< returns inverse or opposite rotation */
    rotv inverse() const {return rotv(- *this);}
    /**< overloaded operator * (multiplication of rotation angle with no change in direction) */
    rotv operator*(const double& factor) const {return ang::rotv(this->get() * factor);} // already ensures that rotation is less than PI
    /**< overloaded operator / (division of rotation angle with no change in direction) */
    rotv operator/(const double& factor) const {return ang::rotv(this->get() / factor);} // already ensures that rotation is less than PI
    /**< executes object rotation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns input fraction (interpolation or extrapolation) of the object. */
    rotv pow(const double& t) const  {return ang::rotv(this->get() * t);} // already ensures that rotation is less than PI
    /**< spherical linear interpolation, returns rotv0 for t=0 and rotv1 for t=1 */
    static rotv slerp(const rotv& rotv0, const rotv& rotv1, const double& t) {return rotv0 * (rotv0 / rotv1).pow(t);} // already ensures that rotation is less than PI
    /**< right plus operator (input rotation located in local tangent space) */
    rotv plus_right(const rotv& rv) const {return (*this) * rv;}
    /**< left plus operator (input rotation located in global tangent space) */
    rotv plus_left(const rotv& rv) const {return rv * (*this);}

    /**< ===== ===== Operations ===== ===== */
	/**< overloaded operator * (forward rotation) */
	Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward rotation) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;
    /**< right minus operator (output rotation located in local tangent space). */
    ang::rotv minus_right(const rotv& rv) const {return rv.inverse() * (*this);}
    /**< left minus operator (output rotation located in global tangent space). */
    ang::rotv minus_left(const rotv& rv) const {return (*this) * rv.inverse();}
    /**< exponential map that returns rodrigues parameters */
    ang::rodrigues exp_map_rodrigues() const;
    /**< exponential map that returns direction cosine matrix */
    ang::dcm exp_map_dcm() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::so3_tangent operator|(const ang::so3_tangent& w) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_forward() const;
    /**< overloaded operator % (backward adjoint) */
    ang::so3_tangent operator%(const ang::so3_tangent& w) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_backward() const;

    /**< ===== ===== Angular Velocity - Time Derivative ===== ===== */
	/**< obtains the body angular velocity from the rotation vector and	its time derivative. */
    ang::so3_tangent dot2omegabody(const Eigen::Vector3d& rvdot) const;
    /**< obtains the rotation vector differential with time based on the rotation vector and the body angular velocity. */
    Eigen::Vector3d omegabody2dot(const ang::so3_tangent& w_body_rps) const;
    /**< obtains the space angular velocity from the rotation vector and its time derivative. */
    ang::so3_tangent dot2omegaspace(const Eigen::Vector3d& rvdot) const;
    /**< obtains the rotation vector differential with time based on the rotation vector and the space angular velocity. */
    Eigen::Vector3d omegaspace2dot(const ang::so3_tangent& w_space_rps) const;

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 3x1 vector */
    Eigen::Vector3d& operator()() {return *this;}
    const Eigen::Vector3d& operator()() const {return *this;}
    /**< get 3x1 vector */
    Eigen::Vector3d& get() {return *this;}
    const Eigen::Vector3d& get() const {return *this;}
    /**< This function returns a very good approximation (specially when the two input rotations are close) to (rv1 / rv2),
     * computed as (rv1() - rv2()). It also modifies one of the inputs, so it can be used repeatedly with the same
     * rv2 without going into the conditional. *
     * This function shall be used exclusively when the intended difference between both inputs is relatively small.
     * The reason is NOT that this approximation is better the smaller their difference is (check the test_rotv_eclidean_diff
     * method), but that to compute the initial difference on which to apply the "if" condition (not the final result,
     * which in this case is the same), both inputs should be rotv (here one is Vector3d), and it would be (rv1 / rv2).norm().
     * But then we can not change the second one as the angle has to be less than 180 [deg], and the equivalent_rotv function
     * changes it to between 180 and 360 [deg]. */
    static Eigen::Vector3d euclidean_diff(const ang::rotv& rv1, Eigen::Vector3d& r2);

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the three components of the
     * rotation vector and return them in vector form. */
    static Eigen::Vector3d wedge(const ang::rotv& r) {return r();}
    /**< although the hat operator usually applies to the tangent space, here it takes the three components of the
     * rotation vector in vector form, and returns them in rotation vector form (object). */
    static ang::rotv hat(const Eigen::Vector3d& v) {return ang::rotv(v);}

    /**< ===== ===== Obtain Individual Euler angles ===== ===== */

    /**< ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse rotation, equal to d(rv^-1)/dDeltarB, which coincides with the negative of the adjoint.
     * (rv plus DeltarB)^-1 = rv^-1 + J * DeltarB */
    Eigen::Matrix3d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(rv1 * rv2) / dDeltarB1, which coincides with the inverse of the 2nd rotation adjoint.
     * (rv1 plus DeltarB1) * rv2 = rv1 * rv2 plus J * DeltarB1 */
    Eigen::Matrix3d jac_right_composition_wrt_first(const rotv& r) const {return r.adjoint_matrix_backward();}
    /**< returns the right jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(rv1 * rv2) / dDeltarB2, which as a matter of fact is the identity.
     * rv1 * (rv2 plus DeltarB2) = rv1 * rv2 plus J * DeltarB2 */
    Eigen::Matrix3d jac_right_composition_wrt_second(const rotv& r) const {return Eigen::Matrix3d::Identity();}

    /**< returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(rv * v) / dDeltarB,
     * (rv plus DeltarB) * v = rv * v + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(rv / v) / dDeltarB,
     * (rv plus DeltarB) / v = rv / v + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(rv)) / dDeltarv, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(rv plus Deltarv) ~= Log(rv) + (J * Deltarv) */
    Eigen::Matrix3d jac_right_log() const;

    /**< returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
     * equal to d(rv1 plus rv2) / dDeltar1B.
     * (rv1 plus Deltar1B) plus rv2 = (rv1 plus rv2) plus J * Deltar1B */
    Eigen::Matrix3d jac_right_plus_wrt_first(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
     * equal to d(rv1 plus rv2) / dDeltar2B.
     * rv1 plus (rv2 + Deltar2B) = (rv1 plus rv2) plus J * Deltar2B */
    Eigen::Matrix3d jac_right_plus_wrt_second(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the first element,
     * equal to d(rv2 minus rv1) / dDeltar2B.
     * (rv2 plus Deltar2B) minus rv1 = (rv2 minus rv1) + J * Deltar2B */
    Eigen::Matrix3d jac_right_minus_wrt_first(const rotv& r) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the second element,
     * equal to d(rv2 minus rv1) / dDeltar1B.
     * rv2 minus(rv1 plus Deltar1B) = (rv2 minus rv1) + J * Deltar1B */
    Eigen::Matrix3d jac_right_minus_wrt_second(const rotv& r) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(Adr | w) / dDeltarB,
     * Ad(r plus DeltarB) | w = Adr | w + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(Adr % w) / dDeltarB,
     * Ad(r plus DeltarB) % w = Adr % w + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse rotation, equal to d(rv^-1)/dDeltarN, which coincides with the negative of the inverse adjoint.
     * (DeltarN plus rv)^-1 = J * DeltarN plus rv^-1 */
    Eigen::Matrix3d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(rv1 * rv2) / dDeltarN1, which as a matter of fact is the identity.
     * (DeltarN1 plus rv1) * rv2 = J * DeltarN1 plus rv1 * rv2 */
    Eigen::Matrix3d jac_left_composition_wrt_first(const rotv& rv) const {return Eigen::Matrix3d::Identity();}
    /**< returns the left jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(rv1 * rv2) / dDeltarN2, which coincides with the adjoint of the 1st rotation.
     * rv1 * (DeltarN2 plus rv2) = J * DeltarN2 plus rv1 * rv2 */
    Eigen::Matrix3d jac_left_composition_wrt_second(const rotv& rv) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(rv * v) / dDeltarN,
     * (DeltarN plus rv) * v = rv * v + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(rv / v) / dDeltarN,
     * (DeltarN plus rv) / v = rv / v + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(rv)) / dDeltarv, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(Deltarv plus rv) ~= Log(rv) + (J * Deltarv) */
    Eigen::Matrix3d jac_left_log() const;

    /**< returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
     * equal to d(rv1 plus rv2) / dDeltar1N.
     * (rv1 + Deltar1N) plus rv2 = J * Deltar1N plus (rv1 plus rv2) */
    Eigen::Matrix3d jac_left_plus_wrt_first(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
     * equal to d(rv1 plus rv2) / dDeltar2N.
     * rv1 plus (Deltar2N plus rv2) = J * Deltar2N plus (rv1 plus rv2) */
    Eigen::Matrix3d jac_left_plus_wrt_second(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the first element,
     * equal to d(rv2 minus rv1) / dDeltar2N.
     * (Deltar2N plus rv2) minus rv1 = (rv2 minus rv1) + J * Deltar2N */
    Eigen::Matrix3d jac_left_minus_wrt_first(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the second element,
     * equal to d(rv2 minus rv1) / dDeltar1N.
     * rv2 minus (Deltar1N plus rv1) = (rv2 minus rv1) + J * Deltar1N */
    Eigen::Matrix3d jac_left_minus_wrt_second(const rotv& rv) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(Adrv | w) / dDeltarN,
     * Ad(DeltarN plus rv) | w = Adrv | w + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(Adrv % w) / dDeltarN,
     * Ad(DeltarN plus rv) % w = Adrv % w + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Exponential (Right and Left) Jacobians ===== ===== */
    /**< returns the jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector itself,
     * also known as the right Jacobian of SO(3), equal to d(Exp(rv)) / dDeltarv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(rv + Deltarv) ~= exp(rv) plus (JR * Deltarv) */
    Eigen::Matrix3d jac_right() const;
    /**< returns the inverse jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector
     * itself, also known as the inverse of the right Jacobian of SO(3).
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * rv + JRinv * Deltarv ~= log[exp(rv) plus Deltarv] */
    Eigen::Matrix3d jac_right_inv() const;
    /**< returns the jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector itself,
     * also known as the left Jacobian of SO(3), equal to d(Exp(rv)) / dDeltarv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(rv + Deltarv) ~= (JL * Deltarv) plus_left exp(rv) */
    Eigen::Matrix3d jac_left() const;
    /**< returns the inverse jacobian [3x3] of the exponential of the rotation vector with respect to the rotation vector
     * itself, also known as the inverse of the left Jacobian of SO(3).
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * rv + JLinv * Deltarv ~= log[Deltarv plus_left exp(rv)] */
    Eigen::Matrix3d jac_left_inv() const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward rotation action with respect to the vector, equal to d(rv * v) / dDeltav,
     * rv * v = J * v
     * rv * (v + Deltav) = rv * v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_forward_rotation_wrt_vector() const;
    /**< returns the jacobian of the backward rotation action with respect to the vector, equal to d(rv / v) / dDeltav,
     * rv / v = J * v
     * rv / (v + Deltav) = rv / v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_backward_rotation_wrt_vector() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adr | w) / dDeltaw,
     * Adr | w = J * w
     * Adr | (w + Deltaw) = Adr | w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adr % w) / dDeltaw,
     * Adr % w = J * w
     * Adr % (w + Deltaw) = Adr % w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x3] of a forward rotation with respect to the rotation vector, equal to d(r * v)/dr.
     * The forward rotation is NOT linear on the rotation vector, so
     * r * v != d(r * v)/dr |r*v * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) * v ~= r * v + d(r * v)/dr |r*v * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the right jacobian, but the result
     * with the left jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_forward_rotation_wrt_rotv(const Eigen::Vector3d& v) const;
    /**< returns the jacobian [3x3] of a forward rotation with respect to the rotation vector, equal to d(r * v)/dr.
     * The forward rotation is NOT linear on the rotation vector, so
     * r * v != d(r * v)/dr |r*v * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) * v ~= r * v + d(r * v)/dr |r*v * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the left jacobian, but the result
     * with the right jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_forward_rotation_wrt_rotv_bis(const Eigen::Vector3d& v) const;
    /**< returns the jacobian [3x3] of a backward rotation with respect to the rotation vector, equal to d[r / v]/dr.
     * The backward rotation is NOT linear on the rotation vector, so
     * r / v != d(r / v)/dr |r/v * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) / v ~= r / v + d(r / v)/dr |r/v * Delta r
     * Note that the increment is Delta r. */
    Eigen::Matrix3d jac_euclidean_backward_rotation_wrt_rotv(const Eigen::Vector3d& v) const;
    /**< returns the jacobian [3x3] of a backward rotation with respect to the rotation vector, equal to d(r / v)/dr.
     * The backward rotation is NOT linear on the rotation vector, so
     * r / v != d(r / v)/dr |r/v * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) / v ~= r / v + d(r / v)/dr |r/v * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the left jacobian, but the result
     * with the right jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_backward_rotation_wrt_rotv_bis(const Eigen::Vector3d& v) const;

    /**< returns the jacobian [3x3] of a forward adjoint with respect to the rotation vector, equal to d(r | w)/dr.
     * The forward adjoint is NOT linear on the rotation vector, so
     * r | w != d(r | w)/dr |r|w * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) | w ~= r | w + d(r | w)/dr |r|w * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the right jacobian, but the result
     * with the left jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_forward_adjoint_wrt_rotv(const ang::so3_tangent& w) const;
    /**< returns the jacobian [3x3] of a forward adjoint with respect to the rotation vector, equal to d(r | w)/dr.
     * The forward adjoint is NOT linear on the rotation vector, so
     * r | w != d(r | w)/dr |r|w * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) | w ~= r | w + d(r | w)/dr |r|w * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the left jacobian, but the result
     * with the right jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_forward_adjoint_wrt_rotv_bis(const ang::so3_tangent& w) const;
    /**< returns the jacobian [3x3] of a backward adjoint with respect to the rotation vector, equal to d[r % w]/dr.
     * The backward adjoint is NOT linear on the rotation vector, so
     * r % w != d(r % w)/dr |r%w * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) % w ~= r % w + d(r % w)/dr |r%w * Delta r
     * Note that the increment is Delta r. */
    Eigen::Matrix3d jac_euclidean_backward_adjoint_wrt_rotv(const ang::so3_tangent& w) const;
    /**< returns the jacobian [3x3] of a backward adjoint with respect to the rotation vector, equal to d(r % w)/dr.
     * The backward adjoint is NOT linear on the rotation vector, so
     * r % w != d(r % w)/dr |r%w * r
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * (r + Delta r) % w ~= r % w + d(r % w)/dr |r%w * Delta r
     * Note that the increment is Delta r.
     * Note that internally the expression is based on the left jacobian, but the result
     * with the right jacobian is the same. */
    Eigen::Matrix3d jac_euclidean_backward_adjoint_wrt_rotv_bis(const ang::so3_tangent& w) const;



}; // closes class rotv
    
} // closes namespace ang

#endif


























