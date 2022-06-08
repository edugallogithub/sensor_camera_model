#ifndef ANG_REPRESENTATION_SPECIAL_EUCLIDEAN_DCM
#define ANG_REPRESENTATION_SPECIAL_EUCLIDEAN_DCM

#include "../ang.h"
#include "../rotate/dcm.h"

/*
 * This file contains the "speu_dcm" class that models the special Euclidean SE(3)
 * transformation composed by a rotation (based on dcm) and a translation.
 */

namespace ang {
    class speu_rodrigues;
    class homogeneous;
    class trfv;
    class screw;
    class dual;
    class se3_tangent;
    class se3_tangent_homo;

// CLASS SPECIAL EUCLIDEAN (based on DCM)
// ======================================
// ======================================

class ANG_API speu_dcm {
private:
    /**< direction cosine matrix representing a rotation */
    ang::dcm _R;
    /**< vector representing translation */
    Eigen::Vector3d _T;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< default constructor (not initialized) */
    speu_dcm() = default;
    /**< constructor based on direction cosine matrix and translation */
    speu_dcm(const ang::dcm& R, const Eigen::Vector3d& T) : _R(R), _T(T) {}
    /**< constructor based on inverse direction cosine matrix and inverse translation.
     * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */
    speu_dcm(const ang::dcm& R_inv, const Eigen::Vector3d& T_inv, bool not_used) : _R(R_inv.inverse()), _T(R_inv.inverse() * (-T_inv)) {}
    /**< constructor based on special Euclidean (rodrigues) */
    explicit speu_dcm(const speu_rodrigues&);
    /**< constructor based on homogeneous */
    explicit speu_dcm(const homogeneous&);
    /**< constructor based on transform vector */
    explicit speu_dcm(const trfv&);
    /**< constructor based on screw */
    explicit speu_dcm(const screw&);
    /**< constructor based on unit dual quaternion */
    explicit speu_dcm(const dual&);
    /**< copy constructor */
    speu_dcm(const speu_dcm&) = default;
    /**< destructor */
    virtual ~speu_dcm() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    speu_dcm(speu_dcm&&) = default;
    /**< move constructor based on direction cosine matrix and translation */
    speu_dcm(ang::dcm&& R, Eigen::Vector3d&& T) : _R(R), _T(T) {}
    /**< move constructor based on inverse direction cosine matrix and inverse translation.
     * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */
    speu_dcm(ang::dcm&& R_inv, Eigen::Vector3d&& T_inv, bool not_used)  : _R(R_inv.inverse()), _T(R_inv.inverse() * (-T_inv)) {}
    /**< move constructor based on special Euclidean (rodrigues) */
    explicit speu_dcm(speu_rodrigues&&);
    /**< move constructor based on homogeneous */
    explicit speu_dcm(homogeneous&&);
    /**< move constructor based on transform vector */
    explicit speu_dcm(trfv&&);
    /**< move constructor based on screw */
    explicit speu_dcm(screw&&);
    /**< move constructor based on unit dual quaternion */
    explicit speu_dcm(dual&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    speu_dcm& operator=(const speu_dcm&) = default;
    /**< assignment operator = based on special Euclidean (rodrigues) */
    speu_dcm& operator=(const speu_rodrigues&);
    /**< assignment operator = based on homogeneous */
    speu_dcm& operator=(const homogeneous&);
    /**< assignment operator = based on transform vector */
    speu_dcm& operator=(const trfv&);
    /**< assignment operator = based on screw */
    speu_dcm& operator=(const screw&);
    /**< assignment operator = based on unit dual quaternion */
    speu_dcm& operator=(const dual&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    speu_dcm& operator=(speu_dcm&&) = default;
    /**< move assignment operator = based on special Euclidean (rodrigues) */
    speu_dcm& operator=(speu_rodrigues&&);
    /**< move assignment operator = based on homogeneous */
    speu_dcm& operator=(homogeneous&&);
    /**< move assignment operator = based on transform vector */
    speu_dcm& operator=(trfv&&);
    /**< move assignment operator = based on screw */
    speu_dcm& operator=(screw&&);
    /**< move assignment operator = based on unit dual quaternion */
    speu_dcm& operator=(dual&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    speu_dcm operator*(const speu_dcm& op2) const {return ang::speu_dcm(_R * op2._R, _R * op2._T + _T);}
    /**< overloaded operator / (backward combination of transformations) */
    speu_dcm operator/(const speu_dcm& op2) const;
    /**< returns inverse or opposite transformation */
    speu_dcm inverse() const {return ang::speu_dcm(_R.inverse(), _R.inverse() * (- _T));}
    /**< executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns exponential map of the power function applied to the object screw logarithmic map. */
    speu_dcm pow(const double& t) const;
    /**< screw linear interpolation, returns gR0 for t=0 and gR1 for t=1 */
    static speu_dcm sclerp(const speu_dcm& gR0, const speu_dcm& gR1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    speu_dcm plus_right(const trfv&) const;
    /**< right plus operator (input rotation located in local tangent space) */
    speu_dcm plus_right(const screw&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    speu_dcm plus_left(const trfv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    speu_dcm plus_left(const screw&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d& vecin) const {return _R * vecin + _T;}
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d& vecin) const {return _R / vecin - _R / _T;}
    /**< overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -.*/
    Eigen::Vector3d operator^(const Eigen::Vector3d& vecin) const {return _R * vecin;}
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d& vecin) const {return _R / vecin;}
    /**< right minus operator (output transformation located in local tangent space) */
    trfv minus_right_trfv(const speu_dcm&) const;
    /**< right minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const speu_dcm&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const speu_dcm&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const speu_dcm&) const;
    /**< logarithmic map that returns the transform vector */
    ang::trfv log_map_trfv() const;
    /**< logarithmic map that returns the screw */
    ang::screw log_map_screw() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::se3_tangent_homo operator|(const ang::se3_tangent_homo& xi_homo) const;
    /**< overloaded operator | (forward adjoint) */
    ang::se3_tangent operator|(const ang::se3_tangent& xi) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_forward() const;
    /**< overloaded operator % (backward adjoint) */
    ang::se3_tangent_homo operator%(const ang::se3_tangent_homo& xi_homo) const;
    /**< overloaded operator % (backward adjoint) */
    ang::se3_tangent operator%(const ang::se3_tangent& xi) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_backward() const;

    /**< ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
    /**< obtains the body twist or motion velocity from the special euclidean direct cosine matrix and its time derivative */
    ang::se3_tangent dot2xibody(const Eigen::Matrix34d& speu_dcmdot) const;
    /**< obtains the special euclidean direct cosine matrix derivative with time based on the homogeneous
    transformation and the body twist or motion velocity. */
    Eigen::Matrix34d xibody2dot(const ang::se3_tangent& xi_body_mrps) const;
    /**< obtains the space twist or motion velocity from the special euclidean direct cosine matrix and its time derivative */
    ang::se3_tangent dot2xispace(const Eigen::Matrix34d& speu_dcmdot) const;
    /**< obtains the special euclidean direct cosine matrix derivative with time based on the homogeneous
    transformation and the space twist or motion velocity. */
    Eigen::Matrix34d xispace2dot(const ang::se3_tangent& xi_space_mrps) const;

    /**< ===== ===== Linear Algebra ===== ===== */

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
     * special euclidean dcm form form and returns them in vector form, first those of the rotation matrix (ordered
     * by row), and last the translation vector. */
    static Eigen::Vector12d wedge(const ang::speu_dcm& Gr);
    /**< although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
     * special euclidean dcm form form and returns them in vector form, first those of the rotation matrix (ordered
     * by row), and last the translation vector. */
    static Eigen::Vector12d wedge(const Eigen::Matrix34d& Gr);
    /**< although the hat operator usually applies to the tangent space, here it takes the twelve components of the
     * vector form, first those of the rotation matrix (ordered by row), and last the translation vector, and returns
     * them in special euclidean dcm form. */
    static ang::speu_dcm hat(const Eigen::Vector12d& v);

    /**< ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse motion, equal to d(gR^-1)/dDeltatauB, which coincides with the negative of the adjoint.
     * (gR plus DeltatauB)^-1 = gR^-1 + J * DeltatauB */
    Eigen::Matrix6d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(gR1 * gR2) / dDeltatauB1, which coincides with the inverse of the 2nd motion adjoint.
     * (gR1 plus DeltatauB1) * gR2 = gR1 * gR2 plus J * DeltatauB1 */
    Eigen::Matrix6d jac_right_composition_wrt_first(const speu_dcm& gR) const {return gR.adjoint_matrix_backward();}
    /**< returns the right jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(gR1 * gR2) / dDeltatauB2, which as a matter of fact is the identity.
     * gR1 * (gR2 plus DeltatauB2) = gR1 * gR2 plus J * DeltatauB2 */
    Eigen::Matrix6d jac_right_composition_wrt_second(const speu_dcm& gR) const {return Eigen::Matrix6d::Identity();}

    /**< returns the right jacobian of the forward motion action with respect to the motion, equal to d(gR * p) / dDeltatauB,
     * (gR plus DeltatauB) * p = gR * p + J * DeltatauB */
    Eigen::Matrix36d jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the right jacobian of the backward motion action with respect to the motion, equal to d(gR / p) / dDeltatauB,
     * (gR plus DeltatauB) / p = gR / p + J * DeltatauB */
    Eigen::Matrix36d jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(gR)) / dDeltatau, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(gR plus Deltatau) ~= Log(gR) + (J * Deltatau) */
    Eigen::Matrix6d jac_right_log() const;

    /**< returns the right jacobian of the motion right plus operator with respect to the group object (first element),
     * equal to d(gR1 plus tau2) / dDeltatau1B.
     * (gR1 plus Deltatau1B) plus tau2 = (gR1 plus tau2) plus J * Deltatau1B */
    Eigen::Matrix6d jac_right_plus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
     * equal to d(gR1 plus tau2) / dDeltatau2B.
     * gR1 plus (tau2 + Deltatau2B) = (gR1 plus tau2) plus J * Deltatau2B */
    Eigen::Matrix6d jac_right_plus_wrt_second(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the first element,
     * equal to d(gR2 minus gR1) / dDeltatau2B.
     * (gR2 plus Deltatau2B) minus gR1 = (gR2 minus gR1) + J * Deltatau2B */
    Eigen::Matrix6d jac_right_minus_wrt_first(const speu_dcm& gR) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the second element,
     * equal to d(gR2 minus gR1) / dDeltatau1B.
     * gR2 minus(gR1 plus Deltatau1B) = (gR2 minus gR1) + J * Deltatau1B */
    Eigen::Matrix6d jac_right_minus_wrt_second(const speu_dcm& gR) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdgR | xi) / dDeltatauB,
     * Ad(gR plus DeltatauB) | xi = AdgR | xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdgR % xi) / dDeltatauB,
     * Ad(gR plus DeltatauB) % xi = AdgR % xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse motion, equal to d(gR^-1)/dDeltatauE, which coincides with the negative of the inverse adjoint.
     * (DeltatauE plus gR)^-1 = J * DeltatauE plus gR^-1 */
    Eigen::Matrix6d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(gR1 * gR2) / dDeltatauE1, which as a matter of fact is the identity.
     * (DeltatauE1 plus gR1) * gR2 = J * DeltatauE1 plus gR1 * gR2 */
    Eigen::Matrix6d jac_left_composition_wrt_first(const speu_dcm& gR) const {return Eigen::Matrix6d::Identity();}
    /**< returns the left jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(gR1 * gR2) / dDeltatauE2, which coincides with the adjoint of the 1st motion.
     * gR1 * (DeltatauE2 plus gR2) = J * DeltatauE2 plus gR1 * gR2 */
    Eigen::Matrix6d jac_left_composition_wrt_second(const speu_dcm& gR) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward motion action with respect to the motion, equal to d(gR * p) / dDeltatauE,
     * (DeltatauE plus gR) * p = gR * p + J * DeltatauE */
    Eigen::Matrix36d jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the left jacobian of the backward motion action with respect to the motion, equal to d(gR / p) / dDeltatauE,
     * (DeltatauE plus gR) / p = gR / p + J * DeltatauE */
    Eigen::Matrix36d jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(gR)) / dDeltatau, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(Deltatau plus gR) ~= Log(gR) + (J * Deltatau) */
    Eigen::Matrix6d jac_left_log() const;

    /**< returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
     * equal to d(tau1 plus gR2) / dDeltatau1E.
     * (tau1 + Deltatau1E) plus gR2 = J * Deltatau1E plus (tau1 plus gR2)  */
    Eigen::Matrix6d jac_left_plus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left plus operator with respect to the group object (second element),
     * equal to d(tau1 plus gR2) / dDeltatau2E.
     * tau1 plus (Deltatau2E plus gR2) = J * Deltatau2E plus (tau1 plus gR2) */
    Eigen::Matrix6d jac_left_plus_wrt_second(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the first element,
     * equal to d(gR2 minus gR1) / dDeltatau2E.
     * (Deltatau2E plus gR2) minus gR1 = (gR2 minus gR1) + J * Deltatau2E */
    Eigen::Matrix6d jac_left_minus_wrt_first(const speu_dcm& gR) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the second element,
     * equal to d(gR2 minus gR1) / dDeltatau1E.
     * gR2 minus (Deltatau1E plus gR1) = (gR2 minus gR1) + J * Deltatau1E */
    Eigen::Matrix6d jac_left_minus_wrt_second(const speu_dcm& gR) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(AdgR | xi) / dDeltatauE,
     * Ad(DeltatauE plus gR) | xi = AdgR | xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(AdgR % xi) / dDeltatauE,
     * Ad(DeltatauE plus gR) % xi = AdgR % xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward motion action with respect to the point, equal to d(gR * p) / dDeltap,
     * gR * (p + Deltap) = gR * p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_point() const;
    /**< returns the jacobian of the backward motion action with respect to the point, equal to d(gR / p) / dDeltap,
     * gR / (p + Deltap) = gR / p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_point() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdgR | xi) / dDeltaxi,
     * AdgR | (xi + Deltaxi) = AdgR | xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdgR % xi) / dDeltaxi,
     * AdgR % (xi + Deltaxi) = AdgR % xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x12] of a forward transformation with respect to the affine rotation matrix, equal to d(Gr * p)/dGr.
     * The forward transformation (unlike the backward one) is linear on the affine rotation matrix, so the following expressions are true:
     * Gr * p = d(Gr * p)/dGr |Gr*p * p
     * [Gr * exp(Delta tau)] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[Gr * exp(Delta tau)] - Gr}
     * [exp(Delta tau) * Gr] * p ~= Gr * p + d(Gr * p)/dGr |Gr*p * {[exp(Delta tau) * Gr] - Gr}
     * Note that the jacobian is evaluated at |Gr*p.
     * Note that the increment is {[Gr * exp(Delta tau)] - Gr} or [(exp(Delta tau) * Gr) - Gr], and not exp(Delta tau) or Delta Gr. */
    Eigen::Matrix312d jac_euclidean_forward_motion_wrt_speu_dcm(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x12] of a backward transformation with respect to the affine rotation matrix, equal to d(Gr / p)/dGr.
     * The backward transformation (unlike the forward one) is NOT linear on the affine rotation matrix. The 1st order Taylor
     * approximation is valid for very small affine rotation matrix changes:
     * [Gr * exp(Delta tau)] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[Gr * exp(Delta tau)] - Gr}
     * [exp(Delta tau) * Gr] / p ~= Gr / p + d(Gr / p)/dGr |Gr/p * {[exp(Delta tau) * Gr] - Gr}
     * Note that the jacobian is evaluated at |Gr/p.
     * Note that the increment is {[Gr * exp(Delta tau)] - Gr} or [(exp(Delta tau) * Gr) - Gr], and not exp(Delta tau) or Delta Gr. */
    Eigen::Matrix312d jac_euclidean_backward_motion_wrt_speu_dcm(const Eigen::Vector3d& p) const;

    /**< ===== ===== Getters ===== ===== */
    /**< return translation vector */
    const Eigen::Vector3d& get_T() const {return _T;}
    /**< return direction cosine matrix rotation */
    const ang::dcm& get_dcm() const {return _R;}
    /**< return rotation vector object */
    ang::rotv get_rotv() const;

    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const {return _R.inverse() * (- _T);}
    /**< return direction cosine matrix rotation of inverse or opposite transformation */
    ang::dcm get_inverse_dcm() const {return _R.inverse();}
}; // closes class speu_dcm

} // closes namespace ang

#endif

















