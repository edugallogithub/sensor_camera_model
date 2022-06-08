#ifndef ANG_REPRESENTATION_RODRIGUES
#define ANG_REPRESENTATION_RODRIGUES

#include "../ang.h"
#include "../quat.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * This file contains the "rodrigues" class that models the Rodrigues parameters (unit quaternion) representation of a rotation.
 */

namespace ang {

class euler;
class dcm;
class rotv;
class so3_tangent;
class so3_tangent_quat;

// CLASS RODRIGUES
// ===============
// ===============

class ANG_API rodrigues : private quat {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
	/**< default constructor (not initialized) */
	rodrigues() = default;
	/**< normal constructor */
	rodrigues(double p0, double p1, double p2, double p3) : ang::quat(p0, p1, p2, p3) {this->normalize();}
    /**< constructor based on quaternion */
    explicit rodrigues(const quat& q) : ang::quat(q) {this->normalize();}
    /**< constructor based on Euler angles */
	explicit rodrigues(const euler&);
	/**< constructor based on direction cosine matrix */
	explicit rodrigues(const dcm&);
	/**< constructor based on rotation vector */
	explicit rodrigues(const rotv&);
    /**< copy constructor */
    rodrigues(const rodrigues&) = default;
    /**< destructor */
	~rodrigues() override = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    rodrigues(rodrigues&&) = default;
    /**< move constructor based on quaternion */
    explicit rodrigues(quat&& q) : ang::quat(q) {this->normalize();}
    /**< move constructor based on Euler angles */
    explicit rodrigues(euler&&);
    /**< move constructor based on direction cosine matrix */
    explicit rodrigues(dcm&&);
    /**< move constructor based on rotation vector */
    explicit rodrigues(rotv&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    rodrigues& operator=(const rodrigues&) = default;
    /**< assignment operator = based on quaternion */
    rodrigues& operator=(const quat&);
    /**< assignment operator based on Euler angles */
    rodrigues& operator=(const euler&);
    /**< assignment operator based on direction cosine matrix */
    rodrigues& operator=(const dcm&);
    /**< assignment operator based on rotation vector */
    rodrigues& operator=(const rotv&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    rodrigues& operator=(rodrigues&&) = default;
    /**< move assignment operator = based on quaternion */
    rodrigues& operator=(quat&&);
    /**< move assignment operator based on Euler angles */
    rodrigues& operator=(euler&&);
    /**< move assignment operator based on direction cosine matrix */
    rodrigues& operator=(dcm&&);
    /**< move assignment operator based on rotation vector */
    rodrigues& operator=(rotv&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of rotations) */
	rodrigues operator*(const rodrigues& op2) const;
    /**< overloaded operator / (backward combination of rotations) */
    rodrigues operator/(const rodrigues& op2) const;
    /**< returns inverse or opposite rotation */
    rodrigues inverse() const {return ang::rodrigues(this->adjoint());}
    /**< returns the same rotation but with negative quaternion (all indexes reversed) */
    rodrigues negative() const {return ang::rodrigues(-(*this)()(0), -(*this)()(1), -(*this)()(2), -(*this)()(3));}
    /**< executes object rotation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns exponential map of the power function applied to the object logarithmic map. */
    rodrigues pow(const double& t) const;
    /**< spherical linear interpolation, returns q0 for t=0 and q1 for t=1 */
    static rodrigues slerp(const rodrigues& q0, const rodrigues& q1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    rodrigues plus_right(const rotv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    rodrigues plus_left(const rotv&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward rotation) */
    Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward rotation) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;
    /**< right minus operator (output rotation located in local tangent space). */
    ang::rotv minus_right(const rodrigues&) const;
    /**< left minus operator (output rotation located in global tangent space). */
    ang::rotv minus_left(const rodrigues&) const;
    /**< logarithmic map that returns the rotation vector */
    ang::rotv log_map() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::so3_tangent_quat operator|(const ang::so3_tangent_quat& w_quat) const;
    /**< overloaded operator | (forward adjoint) */
    ang::so3_tangent operator|(const ang::so3_tangent& w) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_forward() const;
    /**< overloaded operator % (backward adjoint) */
    ang::so3_tangent_quat operator%(const ang::so3_tangent_quat& w_quat) const;
    /**< overloaded operator % (backward adjoint) */
    ang::so3_tangent operator%(const ang::so3_tangent& w) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_backward() const;

    /**< ===== ===== Angular Velocity - Time Derivative ===== ===== */
	/**< obtains the body angular velocity from the Rodrigues parameters and their time derivative. */
    ang::so3_tangent dot2omegabody(const quat& rodriguesdot) const;
	/**< obtains the Rodrigues parameters differentials with time based on the Rodrigues parameters and the body angular velocity. */
	quat omegabody2dot(const ang::so3_tangent& w_body_rps) const;
    /**< obtains the space angular velocity from the Rodrigues parameters and their time derivative. */
    ang::so3_tangent dot2omegaspace(const quat& rodriguesdot) const;
    /**< obtains the Rodrigues parameters differentials with time based on the Rodrigues parameters and the space angular velocity. */
    quat omegaspace2dot(const ang::so3_tangent& w_space_rps) const;

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get quaternion */
    ang::quat& operator()() {return *this;}
    const ang::quat& operator()() const {return *this;}
    /**< get quaternion */
    ang::quat& get() {return *this;}
    const ang::quat& get() const {return *this;}

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the four components of the
     * unit quaternion in rodrigues form and returns them in vector form. */
    static Eigen::Vector4d wedge(const ang::rodrigues& q);
    /**< although the hat operator usually applies to the tangent space, here it takes the four components of the
     * unit quaternion in vector form, and returns them in unit quaternion form (object). */
    static ang::rodrigues hat(const Eigen::Vector4d& v);

    /**< ===== ===== Obtain Individual Euler angles ===== ===== */
    /**< returns yaw angle */
    double get_yaw_rad() const {return std::atan2(+2 * (+ (*this)()(1) * (*this)()(2) + (*this)()(0) * (*this)()(3)), 1 - 2 * (std::pow((*this)()(2),2) + std::pow((*this)()(3),2)));};
    /**< returns pitch angle */
    double get_pitch_rad() const {return std::asin(-2 * (- (*this)()(0) * (*this)()(2) + (*this)()(1) * (*this)()(3)));}
    /**< returns bank angle */
    double get_bank_rad() const {return std::atan2(+2 * (+ (*this)()(2) * (*this)()(3) + (*this)()(0) * (*this)()(1)), 1 - 2 * (std::pow((*this)()(1),2) + std::pow((*this)()(2),2)));}

    /**< ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse rotation, equal to d(q^-1)/dDeltarB, which coincides with the negative of the adjoint.
     * (q plus DeltarB)^-1 = q^-1 + J * DeltarB */
    Eigen::Matrix3d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(q1 * q2) / dDeltarB1, which coincides with the inverse of the 2nd rotation adjoint.
     * (q1 plus DeltarB1) * q2 = q1 * qq2 plus J * DeltarB1 */
    Eigen::Matrix3d jac_right_composition_wrt_first(const rodrigues& q) const {return q.adjoint_matrix_backward();}
    /**< returns the right jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(q1 * q2) / dDeltarB2, which as a matter of fact is the identity.
     * q1 * (q2 plus DeltarB2) = q1 * q2 plus J * DeltarB2 */
    Eigen::Matrix3d jac_right_composition_wrt_second(const rodrigues& q) const {return Eigen::Matrix3d::Identity();}

    /**< returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(q * v) / dDeltarB,
     * (q plus DeltarB) * v = q * v + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(q / v) / dDeltarB,
     * (q plus DeltarB) / v = q / v + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(q)) / dDeltarv, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(q plus Deltarv) ~= Log(q) + (J * Deltarv) */
    Eigen::Matrix3d jac_right_log() const;

    /**< returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
     * equal to d(q1 plus rv2) / dDeltar1B.
     * (q1 plus Deltar1B) plus rv2 = (q1 plus rv2) plus J * Deltar1B */
    Eigen::Matrix3d jac_right_plus_wrt_first(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
     * equal to d(q1 plus rv2) / dDeltar2B.
     * q1 plus (rv2 + Deltar2B) = (q1 plus rv2) plus J * Deltar2B */
    Eigen::Matrix3d jac_right_plus_wrt_second(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the first element,
     * equal to d(q2 minus q1) / dDeltar2B.
     * (q2 plus Deltar2B) minus q1 = (q2 minus q1) + J * Deltar2B */
    Eigen::Matrix3d jac_right_minus_wrt_first(const rodrigues& R) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the second element,
     * equal to d(q2 minus q1) / dDeltar1B.
     * q2 minus(q1 plus Deltar1B) = (q2 minus q1) + J * Deltar1B */
    Eigen::Matrix3d jac_right_minus_wrt_second(const rodrigues& q) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(Adq | w) / dDeltarB,
     * Ad(q plus DeltarB) | w = Adq | w + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(Adq % w) / dDeltarB,
     * Ad(q plus DeltarB) % w = Adq % w + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse rotation, equal to d(q^-1)/dDeltarN, which coincides with the negative of the inverse adjoint.
     * (DeltarN plus q)^-1 = J * DeltarN plus q^-1 */
    Eigen::Matrix3d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(q1 * q2) / dDeltarN1, which as a matter of fact is the identity.
     * (DeltarN1 plus q1) * q2 = J * DeltarN1 plus q1 * q2 */
    Eigen::Matrix3d jac_left_composition_wrt_first(const rodrigues& q) const {return Eigen::Matrix3d::Identity();}
    /**< returns the left jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(q1 * q2) / dDeltarN2, which coincides with the adjoint of the 1st rotation.
     * q1 * (DeltarN2 plus q2) = J * DeltarN2 plus q1 * q2 */
    Eigen::Matrix3d jac_left_composition_wrt_second(const rodrigues& q) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(q * v) / dDeltarN,
     * (DeltarN plus q) * v = q * v + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(q / v) / dDeltarN,
     * (DeltarN plus q) / v = q / v + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(q)) / dDeltarv, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(Deltarv plus q) ~= Log(q) + (J * Deltarv) */
    Eigen::Matrix3d jac_left_log() const;

    /**< returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
     * equal to d(rv1 plus q2) / dDeltar1N.
     * (rv1 + Deltar1N) plus q2 = J * Deltar1N plus (rv1 plus q2)  */
    Eigen::Matrix3d jac_left_plus_wrt_first(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
     * equal to d(rv1 plus q2) / dDeltar2N.
     * rv1 plus (Deltar2N plus q2) = J * Deltar2N plus (rv1 plus q2) */
    Eigen::Matrix3d jac_left_plus_wrt_second(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the first element,
     * equal to d(q2 minus q1) / dDeltar2N.
     * (Deltar2N plus q2) minus q1 = (q2 minus q1) + J * Deltar2N */
    Eigen::Matrix3d jac_left_minus_wrt_first(const rodrigues& q) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the second element,
     * equal to d(q2 minus q1) / dDeltar1N.
     * q2 minus (Deltar1N plus q1) = (q2 minus q1) + J * Deltar1N */
    Eigen::Matrix3d jac_left_minus_wrt_second(const rodrigues& q) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(Adq | w) / dDeltarN,
     * Ad(DeltarN plus q) | w = Adq | w + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(Adq % w) / dDeltarN,
     * Ad(DeltarN plus q) % w = Adq % w + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward rotation action with respect to the vector, equal to d(q * v) / dDeltav,
     * q * v = J * v
     * q * (v + Deltav) = q * v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_forward_rotation_wrt_vector() const;
    /**< returns the jacobian of the backward rotation action with respect to the vector, equal to d(q / v) / dDeltav,
     * q / v = J * v
     * q / (v + Deltav) = q / v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_backward_rotation_wrt_vector() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adq | w) / dDeltaw,
     * Adq | w = J * w
     * Adq | (w + Deltaw) = Adq | w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adq % w) / dDeltaw,
     * Adq % w = J * w
     * Adq % (w + Deltaw) = Adq % w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x4] of a forward rotation with respect to the quaternion, equal to d(q * v)/dq.
     * The forward rotation is BI linear on the quaternion, so
     * q * v = 0.5 * d(q * v)/dq |q*v * q
     * The 1st order Taylor approximation is valid for very small quaternion changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [q * exp(Delta r)] * v ~= q * v + d(q * v)/dq |q*v * {[q * exp(Delta r)] - q}
     * [exp(Delta r) * q] * v ~= q * v + d(q * v)/dq |q*v * {[exp(Delta r) * q] - q}
     * Note that the jacobian is evaluated at |q*v.
     * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
     * to each other, and never exp(Delta r) or Delta q. */
    Eigen::Matrix34d jac_euclidean_forward_rotation_wrt_rodrigues(const Eigen::Vector3d& v) const;
    /**< returns the jacobian [3x4] of a backward rotation with respect to the quaternion, equal to d(q / v)/dq.
     * The backward rotation is BI linear on the quaternion, so
     * q / v = 0.5 * d(q / v)/dq |q/v * q
     * The 1st order Taylor approximation is valid for very small quaternion changes, vut it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [q * exp(Delta r)] / v ~= q / v + d(q / v)/dq |q/v * {[q * exp(Delta r)] - q}
     * [exp(Delta r) * q] / v ~= q / v + d(q / v)/dq |q/v * {[exp(Delta r) * q] - q}
     * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
     * to each other, and never exp(Delta r) or Delta q. */
    Eigen::Matrix34d jac_euclidean_backward_rotation_wrt_rodrigues(const Eigen::Vector3d& v) const;

    /**< returns the jacobian [3x4] of a forward adjoint with respect to the quaternion, equal to d(q | w)/dq.
     * The forward adjoint is BI linear on the quaternion, so
     * q | w = 0.5 * d(q | w)/dq |q|w * q
     * The 1st order Taylor approximation is valid for very small quaternion changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [q * exp(Delta r)] | w ~= q | w + d(q | w)/dq |q|w * {[q * exp(Delta r)] - q}
     * [exp(Delta r) * q] | w ~= q | w + d(q | w)/dq |q|w * {[exp(Delta r) * q] - q}
     * Note that the jacobian is evaluated at |q|w.
     * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
     * to each other, and never exp(Delta r) or Delta q. */
    Eigen::Matrix34d jac_euclidean_forward_adjoint_wrt_rodrigues(const ang::so3_tangent& w) const;
    /**< returns the jacobian [3x4] of a backward adjoint with respect to the quaternion, equal to d(q % w)/dq.
     * The backward adjoint is BI linear on the quaternion, so
     * q % w = 0.5 * d(q % w)/dq |q%w * q
     * The 1st order Taylor approximation is valid for very small quaternion changes, vut it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [q * exp(Delta r)] % w ~= q % w + d(q % w)/dq |q%w * {[q * exp(Delta r)] - q}
     * [exp(Delta r) * q] % w ~= q % w + d(q % w)/dq |q%w * {[exp(Delta r) * q] - q}
     * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
     * to each other, and never exp(Delta r) or Delta q. */
    Eigen::Matrix34d jac_euclidean_backward_adjoint_wrt_rodrigues(const ang::so3_tangent& w) const;
}; // closes class rodrigues

} // closes namespace ang

#endif






















