#ifndef ANG_REPRESENTATION_DCM
#define ANG_REPRESENTATION_DCM

#include "../ang.h"
#include "../auxiliary.h" // required for the proper behavior of the << operator
#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * This file contains the "dcm" class that models the direction cosine matrix representation of a rotation.
 */

namespace ang {

class rodrigues;
class euler;
class rotv;
class so3_tangent;
class so3_tangent_skew;

// CLASS DCM
// =========
// =========

class ANG_API dcm : private Eigen::Matrix3d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not initialized) */
    dcm() = default;
    /**< constructor based on 3x3 matrix */
    explicit dcm(const Eigen::Matrix3d& Omat) : Eigen::Matrix3d(Omat) {this->normalize();}
    /**< constructor based on Euler angles */
    explicit dcm(const euler&);
    /**< constructor based on Rodrigues parameters */
    explicit dcm(const rodrigues&);
    /**< constructor based on rotation vector */
    explicit dcm(const rotv&);
    /**< copy constructor */
    dcm(const dcm&) = default;
    /**< destructor */
    ~dcm() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    dcm(dcm&&) = default;
    /**< move constructor based on 3x3 matrix */
    explicit dcm(Eigen::Matrix3d&& Omat) : Eigen::Matrix3d(Omat) {this->normalize();}
    /**< move constructor based on Euler angles */
    explicit dcm(euler&&);
    /**< move constructor based on Rodrigues parameters */
    explicit dcm(rodrigues&&);
    /**< move constructor based on rotation vector */
    explicit dcm(rotv&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    dcm& operator=(const dcm&) = default;
    /**< assignment operator based 3x3 matrix */
    dcm& operator=(const Eigen::Matrix3d&);
    /**< assignment operator based on Euler angles */
    dcm& operator=(const euler&);
    /**< assignment operator based on Rodrigues parameters */
    dcm& operator=(const rodrigues&);
    /**< assignment operator based on rotation vector */
    dcm& operator=(const rotv&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    dcm& operator=(dcm&&) = default;
    /**< move assignment operator based 3x3 matrix */
    dcm& operator=(Eigen::Matrix3d&&);
    /**< move assignment operator based on Euler angles */
    dcm& operator=(euler&&);
    /**< move assignment operator based on Rodrigues parameters */
    dcm& operator=(rodrigues&&);
    /**< move assignment operator based on rotation vector */
    dcm& operator=(rotv&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of rotations) */
    dcm operator*(const dcm& op2) const {return dcm{this->get() * op2()};}
    /**< overloaded operator / (backward combination of rotations) */
    dcm operator/(const dcm& op2) const {return dcm{this->inverse()()* op2()};}
    /**< returns inverse or opposite rotation */
    dcm inverse() const {return dcm{this->transpose()};}
    /**< executes object rotation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns exponential map of the power function applied to the object logarithmic map. */
    dcm pow(const double& t) const;
    /**< spherical linear interpolation, returns R0 for t=0 and R1 for t=1 */
    static dcm slerp(const dcm& R0, const dcm& R1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    dcm plus_right(const rotv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    dcm plus_left(const rotv&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward rotation) */
    Eigen::Vector3d operator*(const Eigen::Vector3d& vecin) const {return this->get() * vecin;}
    /**< overloaded operator / (backward rotation) */
    Eigen::Vector3d operator/(const Eigen::Vector3d& vecin) const {return this->inverse() * vecin;};
    /**< right minus operator (output rotation located in local tangent space). */
    ang::rotv minus_right(const dcm& R) const;
    /**< left minus operator (output rotation located in global tangent space). */
    ang::rotv minus_left(const dcm& R) const;
    /**< logarithmic map that returns the rotation vector. */
    ang::rotv log_map() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::so3_tangent_skew operator|(const ang::so3_tangent_skew& w_skew) const;
    /**< overloaded operator | (forward adjoint) */
    ang::so3_tangent operator|(const ang::so3_tangent& w) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_forward() const {return this->get();}
    /**< overloaded operator % (backward adjoint) */
    ang::so3_tangent_skew operator%(const ang::so3_tangent_skew& w_skew) const;
    /**< overloaded operator % (backward adjoint) */
    ang::so3_tangent operator%(const ang::so3_tangent& w) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix3d adjoint_matrix_backward() const {return this->inverse()();}

    /**< ===== ===== Angular Velocity - Time Derivative ===== ===== */
    /**< obtains the body angular velocity from the direction cosine matrix and its time derivative. */
    ang::so3_tangent dot2omegabody(const Eigen::Matrix3d& dcmdot) const;
    /**< obtains the direction cosine matrix differentials with time based on the direction cosine matrix and the body angular velocity. */
    Eigen::Matrix3d omegabody2dot(const ang::so3_tangent& w_body_rps) const;
    /**< obtains the space angular velocity from the direction cosine matrix and its time derivative. */
    ang::so3_tangent dot2omegaspace(const Eigen::Matrix3d& dcmdot) const;
    /**< obtains the direction cosine matrix differentials with time based on the direction cosine matrix and the space angular velocity. */
    Eigen::Matrix3d omegaspace2dot(const ang::so3_tangent& w_space_rps) const;

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 3x3 matrix */
    Eigen::Matrix3d& operator()() {return *this;}
    const Eigen::Matrix3d& operator()() const {return *this;}
    /**< get 3x3 matrix */
    Eigen::Matrix3d& get() {return *this;}
    const Eigen::Matrix3d& get() const {return *this;}

    /**< ===== ===== Other ===== ===== */
    /**< normalize the rotation matrix ensuring it is orthonormal */
    void normalize();

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the nine components of the
     * rotation matrix in matrix form and returns them in vector form, ordered by row. */
    static Eigen::Vector9d wedge(const ang::dcm& R);
    /**< although the wedge operator usually applies to the tangent space, here it takes the nine components of the
     * rotation matrix in matrix form and returns them in vector form, ordered by row. */
    static Eigen::Vector9d wedge(const Eigen::Matrix3d& R);
    /**< although the hat operator usually applies to the tangent space, here it takes the nine components of the
     * rotation matrix in vector form, ordered by row, and returns them in matrix form (object). */
    static ang::dcm hat(const Eigen::Vector9d& v);

    /**< ===== ===== Obtain Individual Euler angles ===== ===== */
    /**< returns yaw angle */
    double get_yaw_rad() const {return std::atan2((*this)()(1,0), (*this)()(0,0));}
    /**< returns pitch angle */
    double get_pitch_rad() const {return std::asin(- (*this)()(2,0));}
    /**< returns bank angle */
    double get_bank_rad() const {return std::atan2((*this)()(2,1), (*this)()(2,2));}

    /**< ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse rotation, equal to d(R^-1)/dDeltarB, which coincides with the negative of the adjoint.
     * (R plus DeltarB)^-1 = R^-1 plus J * DeltarB */
    Eigen::Matrix3d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(R1 * R2) / dDeltarB1, which coincides with the inverse of the 2nd rotation adjoint.
     * (R1 plus DeltarB1) * R2 = R1 * R2 plus J * DeltarB1 */
    Eigen::Matrix3d jac_right_composition_wrt_first(const dcm& R) const {return R.adjoint_matrix_backward();}
    /**< returns the right jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(R1 * R2) / dDeltarB2, which as a matter of fact is the identity.
     * R1 * (R2 plus DeltarB2) = R1 * R2 plus J * DeltarB2 */
    Eigen::Matrix3d jac_right_composition_wrt_second(const dcm& R) const {return Eigen::Matrix3d::Identity();}

    /**< returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(R * v) / dDeltarB,
     * (R plus DeltarB) * v = R * v + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(R / v) / dDeltarB,
     * (R plus DeltarB) / v = R / v + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(R)) / dDeltarv, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(R plus Deltarv) ~= Log(R) + (J * Deltarv) */
    Eigen::Matrix3d jac_right_log() const;

    /**< returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
     * equal to d(R1 plus rv2) / dDeltar1B.
     * (R1 plus Deltar1B) plus rv2 = (R1 plus rv2) plus J * Deltar1B */
    Eigen::Matrix3d jac_right_plus_wrt_first(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
     * equal to d(R1 plus rv2) / dDeltar2B.
     * R1 plus (rv2 + Deltar2B) = (R1 plus rv2) plus J * Deltar2B */
    Eigen::Matrix3d jac_right_plus_wrt_second(const rotv& rv) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the first element,
     * equal to d(R2 minus R1) / dDeltar2B.
     * (R2 plus Deltar2B) minus R1 = (R2 minus R1) + J * Deltar2B */
    Eigen::Matrix3d jac_right_minus_wrt_first(const dcm& R) const;
    /**< returns the right jacobian of the rotation right minus operator with respect to the second element,
     * equal to d(R2 minus R1) / dDeltar1B.
     * R2 minus(R1 plus Deltar1B) = (R2 minus R1) + J * Deltar1B */
    Eigen::Matrix3d jac_right_minus_wrt_second(const dcm& R) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(AdR | w) / dDeltarB,
     * Ad(R plus DeltarB) | w = AdR | w + J * DeltarB */
    Eigen::Matrix3d jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(AdR % w) / dDeltarB,
     * Ad(R plus DeltarB) % w = AdR % w + J * DeltarB */
    Eigen::Matrix3d jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse rotation, equal to d(R^-1)/dDeltarN, which coincides with the negative of the inverse adjoint.
     * (DeltarN plus R)^-1 = J * DeltarN plus R^-1 */
    Eigen::Matrix3d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the rotation composition with respect to the first (left) rotation, equal to
     * d(R1 * R2) / dDeltarN1, which as a matter of fact is the identity.
     * (DeltarN1 plus R1) * R2 = J * DeltarN1 plus R1 * R2 */
    Eigen::Matrix3d jac_left_composition_wrt_first(const dcm& R) const {return Eigen::Matrix3d::Identity();}
    /**< returns the left jacobian of the rotation composition with respect to the second (right) rotation, equal to
     * d(R1 * R2) / dDeltarN2, which coincides with the adjoint of the 1st rotation.
     * R1 * (DeltarN2 plus R2) = J * DeltarN2 plus R1 * R2 */
    Eigen::Matrix3d jac_left_composition_wrt_second(const dcm& R) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(R * v) / dDeltarN,
     * (DeltarN plus R) * v = R * v + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;
    /**< returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(R / v) / dDeltarN,
     * (DeltarN plus R) / v = R / v + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(R)) / dDeltarv, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * Log(Deltarv plus R) ~= Log(R) + (J * Deltarv) */
    Eigen::Matrix3d jac_left_log() const;

    /**< returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
     * equal to d(rv1 plus R2) / dDeltar1N.
     * (rv1 + Deltar1N) plus R2 = J * Deltar1N plus (rv1 plus R2)  */
    Eigen::Matrix3d jac_left_plus_wrt_first(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
     * equal to d(rv1 plus R2) / dDeltar2N.
     * rv1 plus (Deltar2N plus R2) = J * Deltar2N plus (rv1 plus R2) */
    Eigen::Matrix3d jac_left_plus_wrt_second(const rotv& rv) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the first element,
     * equal to d(R2 minus R1) / dDeltar2N.
     * (Deltar2N plus R2) minus R1 = (R2 minus R1) + J * Deltar2N */
    Eigen::Matrix3d jac_left_minus_wrt_first(const dcm& R) const;
    /**< returns the left jacobian of the rotation left minus operator with respect to the second element,
     * equal to d(R2 minus R1) / dDeltar1N.
     * R2 minus (Deltar1N plus R1) = (R2 minus R1) + J * Deltar1N */
    Eigen::Matrix3d jac_left_minus_wrt_second(const dcm& R) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(AdR | w) / dDeltarN,
     * Ad(DeltarN plus R) | w = AdR | w + J * DeltarN */
    Eigen::Matrix3d jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(AdR % w) / dDeltarN,
     * Ad(DeltarN plus R) % w = AdR % w + J * DeltarN */
    Eigen::Matrix3d jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward rotation action with respect to the vector, equal to d(R * v) / dDeltav,
     * R * v = J * v
     * R * (v + Deltav) = R * v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_forward_rotation_wrt_vector() const;
    /**< returns the jacobian of the backward rotation action with respect to the vector, equal to d(R / v) / dDeltav,
     * R / v = J * v
     * R / (v + Deltav) = R / v + J * Deltav */
    Eigen::Matrix3d jac_euclidean_backward_rotation_wrt_vector() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdR | w) / dDeltaw,
     * AdR | w = J * w
     * AdR | (w + Deltaw) = AdR | w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdR % w) / dDeltaw,
     * AdR % w = J * w
     * AdR % (w + Deltaw) = AdR % w + J * Deltaw */
    Eigen::Matrix3d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x9] of a forward rotation with respect to the rotation matrix, equal to d(R * v)/dR.
     * The forward rotation is linear on the rotation matrix, so
     * R * v = d(R * v)/dR |R*v * R
     * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [R * exp(Delta r)] * v = R * v + d(R * v)/dR |R*v * {[R * exp(Delta r)] - R}
     * [exp(Delta r) * R] * v = R * v + d(R * v)/dR |R*v * {[exp(Delta r) * R] - R}
     * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
     * to each other, and never exp(Delta r) or Delta R. */
    Eigen::Matrix39d jac_euclidean_forward_rotation_wrt_dcm(const Eigen::Vector3d& v) const;
    /**< returns the jacobian [3x9] of a backward rotation with respect to the rotation matrix, equal to d(R / v)/dR.
     * The backward rotation is linear on the rotation matrix, so
     * R / v = d(R / v)/dR |R/v * R
     * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [R * exp(Delta r)] / v = R / v + d(R / v)/dR |R/v * {[R * exp(Delta r)] - R}
     * [exp(Delta r) * R] / v = R / v + d(R / v)/dR |R/v * {[exp(Delta r) * R] - R}
     * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
     * to each other, and never exp(Delta r) or Delta R. */
    Eigen::Matrix39d jac_euclidean_backward_rotation_wrt_dcm(const Eigen::Vector3d& v) const;

    /**< returns the jacobian [3x9] of a forward adjoint with respect to the rotation matrix, equal to d(R | w)/dR.
     * The forward adjoint is linear on the rotation matrix, so
     * R | w = d(R | w)/dR |R|w * R
     * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [R * exp(Delta r)] | w = R | w + d(R | w)/dR |R|w * {[R * exp(Delta r)] - R}
     * [exp(Delta r) * R] | w = R | w + d(R | w)/dR |R|w * {[exp(Delta r) * R] - R}
     * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
     * to each other, and never exp(Delta r) or Delta R. */
    Eigen::Matrix39d jac_euclidean_forward_adjoint_wrt_dcm(const ang::so3_tangent& w) const;
    /**< returns the jacobian [3x9] of a backward adjoint with respect to the rotation matrix, equal to d(R % w)/dR.
     * The backward adjoint is linear on the rotation matrix, so
     * R % w = d(R % w)/dR |R%w * R
     * The 1st order Taylor approximation is valid for very small rotation matrix changes, but it does not matter if
     * the perturbation is local or global, as the two following expressions are both valid:
     * [R * exp(Delta r)] % w = R % w + d(R % w)/dR |R%w * {[R * exp(Delta r)] - R}
     * [exp(Delta r) * R] % w = R % w + d(R % w)/dR |R%w * {[exp(Delta r) * R] - R}
     * Note that the increment is either {[R * exp(Delta r)] - R} or {[exp(Delta r) * R] - R}, which are not equal
     * to each other, and never exp(Delta r) or Delta R. */
    Eigen::Matrix39d jac_euclidean_backward_adjoint_wrt_dcm(const ang::so3_tangent& w) const;
}; // closes class dcm

} // closes namespace ang

#endif






























