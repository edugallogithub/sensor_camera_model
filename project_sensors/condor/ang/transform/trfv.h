#ifndef ANG_REPRESENTATION_TRFV
#define ANG_REPRESENTATION_TRFV

#include "../ang.h"
#include "../auxiliary.h"
#include "homogeneous.h"
#include "../rotate/rotv.h"

/*
 * This file contains the "trfv" or "exponential" representation of a transformation.
 */

namespace ang {
    class speu_rodrigues;
    class speu_dcm;
    class screw;
    class dual;
    class se3_tangent;
    class se3_tangent_homo;

// CLASS TRANSFORM VECTOR
// ======================
// ======================

class ANG_API trfv : private Eigen::Vector6d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< default constructor (not initialized) */
    trfv() = default;
    /**< constructor based on size 6 vector */
    explicit trfv(const Eigen::Vector6d& t) : Eigen::Vector6d(t) {}
    /**< constructor based on transform vector components */
    trfv(double s1, double s2, double s3, double r1, double r2, double r3);
    /**< constructor based on rotation vector and translation vector */
    trfv(const ang::rotv&, const Eigen::Vector3d& T);
    /**< constructor based on special Euclidean (rodrigues) */
    explicit trfv(const ang::speu_rodrigues&);
    /**< constructor based on special Euclidean (dcm) */
    explicit trfv(const ang::speu_dcm&);
    /**< constructor based on homogeneous */
    explicit trfv(const ang::homogeneous&);
    /**< constructor based on screw */
    explicit trfv(const ang::screw&);
    /**< constructor based on unit dual quaternion */
    explicit trfv(const ang::dual&);
    /**< copy constructor */
    trfv(const trfv&) = default;
    /**< destructor */
    ~trfv() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    trfv(trfv&&) = default;
    /**< move constructor based on size 6 vector */
    explicit trfv(Eigen::Vector6d&& t) : Eigen::Vector6d(t) {}
    /**< move constructor based on rotation vector and translation vector */
    trfv(ang::rotv&&, Eigen::Vector3d&& T);
    /**< move constructor based on special Euclidean (rodrigues) */
    explicit trfv(ang::speu_rodrigues&&);
    /**< move constructor based on special Euclidean (dcm) */
    explicit trfv(ang::speu_dcm&&);
    /**< move constructor based on homogeneous */
    explicit trfv(ang::homogeneous&&);
    /**< move constructor based on screw */
    explicit trfv(ang::screw&&);
    /**< move constructor based on unit dual quaternion */
    explicit trfv(ang::dual&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    trfv& operator=(const trfv&) = default;
    /**< assignment operator = based on size 6 vector */
    trfv& operator=(const Eigen::Vector6d&);
    /**< assignment operator = based on special Euclidean (rodrigues) */
    trfv& operator=(const ang::speu_rodrigues&);
    /**< assignment operator = based on special Euclidean (dcm) */
    trfv& operator=(const ang::speu_dcm&);
    /**< assignment operator = based on homogeneous */
    trfv& operator=(const ang::homogeneous&);
    /**< assignment operator = based on screw */
    trfv& operator=(const ang::screw&);
    /**< assignment operator = based on unit dual quaternion */
    trfv& operator=(const ang::dual&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    trfv& operator=(trfv&&) = default;
    /**< move assignment operator based on size 6 vector */
    trfv& operator=(Eigen::Vector6d&&);
    /**< move assignment operator = based on special Euclidean (rodrigues) */
    trfv& operator=(ang::speu_rodrigues&&);
    /**< move assignment operator = based on special Euclidean (dcm) */
    trfv& operator=(ang::speu_dcm&&);
    /**< move assignment operator = based on homogeneous */
    trfv& operator=(ang::homogeneous&&);
    /**< move assignment operator = based on screw */
    trfv& operator=(ang::screw&&);
    /**< move assignment operator = based on unit dual quaternion */
    trfv& operator=(ang::dual&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    trfv operator*(const trfv& op2) const {return trfv(this->get_rotv() * op2.get_rotv(), this->get_rotv() * op2.get_T() + this->get_T());}
    /**< overloaded operator / (backward combination of transformations) */
    trfv operator/(const trfv& op2) const {return trfv(this->get_rotv() / op2.get_rotv(), this->get_rotv() / (op2.get_T() - this->get_T()));}
    /**< returns inverse or opposite transformation. It can be replaced by the negative (as in the rotation vector) ONLY when it is small,
     * this is, when treated as the tangent space. I however sometimes treat is as a normal transformation representation accepting
     * significant pose changes, and then the negative does not provide accurate results. It may be a good idea to implement a
     * a negative method, but then how do I distinguish them? */
    trfv inverse() const {return trfv(this->get_inverse_rotv(), this->get_inverse_T());}
    /**< executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns map of the power function applied to the object screw exponential-log map. */
    trfv pow(const double& t) const;
    /**< screw linear interpolation, returns tau0 for t=0 and tau1 for t=1 */
    static trfv sclerp(const trfv& tau0, const trfv& tau1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    trfv plus_right(const trfv& tau) const {return (*this) * tau;}
    /**< right plus operator (input rotation located in local tangent space) */
    trfv plus_right(const screw&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    trfv plus_left(const trfv& tau) const {return tau * (*this);}
    /**< left plus operator (input rotation located in global tangent space) */
    trfv plus_left(const screw&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;
    /**< overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator^(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d&) const;
    /**< minus operator (output transformation located in local tangent space)  */
    trfv minus_right_trfv(const trfv& tau) const {return tau.inverse() * (*this);}
    /**< minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const trfv&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const trfv& tau) const {return (*this) * tau.inverse();}
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const trfv&) const;
    /**< exponential map that returns special euclidean rodrigues parameters */
    ang::speu_rodrigues exp_map_speu_rodrigues() const;
    /**< exponential map that returns special euclidean direction cosine matrix */
    ang::speu_dcm exp_map_speu_dcm() const;
    /**< exponential map that returns homogeneous matrix */
    ang::homogeneous exp_map_homogeneous() const;
    /**< exponential map that returns the unit dual quaternion */
    ang::dual exp_map_dual() const;
    /**< exponential and logarithmic map (they are the same in this case) that returns the screw */
    ang::screw explog_map_screw() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::se3_tangent operator|(const ang::se3_tangent& xi) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_forward() const;
    /**< overloaded operator % (backward adjoint) */
    ang::se3_tangent operator%(const ang::se3_tangent& xi) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_backward() const;

    /**< ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 6x1 vector */
    Eigen::Vector6d& operator()() {return *this;}
    const Eigen::Vector6d& operator()() const {return *this;}
    /**< get 6x1 vector */
    Eigen::Vector6d& get() {return *this;}
    const Eigen::Vector6d& get() const {return *this;}

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the six components of the
     * transform vector and returns them in vector form. */
    static Eigen::Vector6d wedge(const ang::trfv& tau);
    /**< although the hat operator usually applies to the tangent space, here it takes the six components of the
     * transform vector in vector form, and returns them in transform vector form (object). */
    static ang::trfv hat(const Eigen::Vector6d& v);

    /**< ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse motion, equal to d(tau^-1)/dDeltatauB, which coincides with the negative of the adjoint.
     * (tau plus DeltatauB)^-1 = tau^-1 + J * DeltatauB */
    Eigen::Matrix6d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(tau1 * tau2) / dDeltatauB1, which coincides with the inverse of the 2nd motion adjoint.
     * (tau1 plus DeltatauB1) * tau2 = tau1 * tau2 plus J * DeltatauB1 */
    Eigen::Matrix6d jac_right_composition_wrt_first(const trfv& tau) const {return tau.adjoint_matrix_backward();}
    /**< returns the right jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(tau1 * tau2) / dDeltatauB2, which as a matter of fact is the identity.
     * tau1 * (tau2 plus DeltatauB2) = tau1 * tau2 plus J * DeltatauB2 */
    Eigen::Matrix6d jac_right_composition_wrt_second(const trfv& tau) const {return Eigen::Matrix6d::Identity();}

    /**< returns the right jacobian of the forward motion action with respect to the motion, equal to d(tau * p) / dDeltatauB,
     * (tau plus DeltatauB) * p = tau * p + J * DeltatauB */
    Eigen::Matrix36d jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the right jacobian of the backward motion action with respect to the motion, equal to d(tau / p) / dDeltatauB,
     * (tau plus DeltatauB) / p = tau / p + J * DeltatauB */
    Eigen::Matrix36d jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(tau)) / dDeltatau, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(tau plus Deltatau) ~= Log(tau) + (J * Deltatau) */
    Eigen::Matrix6d jac_right_log() const;

    /**< returns the right jacobian of the motion right plus operator with respect to the group object (first element),
     * equal to d(tau1 plus tau2) / dDeltatau1B.
     * (tau1 plus Deltatau1B) plus tau2 = (tau1 plus tau2) plus J * Deltatau1B */
    Eigen::Matrix6d jac_right_plus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
     * equal to d(tau1 plus tau2) / dDeltatau2B.
     * tau1 plus (tau2 + Deltatau2B) = (tau1 plus tau2) plus J * Deltatau2B */
    Eigen::Matrix6d jac_right_plus_wrt_second(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the first element,
     * equal to d(tau2 minus tau1) / dDeltatau2B.
     * (tau2 plus Deltatau2B) minus tau1 = (tau2 minus tau1) + J * Deltatau2B */
    Eigen::Matrix6d jac_right_minus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the second element,
     * equal to d(tau2 minus tau1) / dDeltatau1B.
     * tau2 minus(tau1 plus Deltatau1B) = (tau2 minus tau1) + J * Deltatau1B */
    Eigen::Matrix6d jac_right_minus_wrt_second(const trfv& tau) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(Adtau | xi) / dDeltatauB,
     * Ad(tau plus DeltatauB) | xi = Adtau | xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(Adtau % xi) / dDeltatauB,
     * Ad(tau plus DeltatauB) % xi = Adtau % xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse motion, equal to d(tau^-1)/dDeltatauE, which coincides with the negative of the inverse adjoint.
     * (DeltatauE plus tau)^-1 = J * DeltatauE plus tau^-1 */
    Eigen::Matrix6d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(tau1 * tau2) / dDeltatauE1, which as a matter of fact is the identity.
     * (DeltatauE1 plus tau1) * tau2 = J * DeltatauE1 plus tau1 * tau2 */
    Eigen::Matrix6d jac_left_composition_wrt_first(const trfv& tau) const {return Eigen::Matrix6d::Identity();}
    /**< returns the left jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(tau1 * tau2) / dDeltatauE2, which coincides with the adjoint of the 1st motion.
     * tau1 * (DeltatauE2 plus tau2) = J * DeltatauE2 plus tau1 * tau2 */
    Eigen::Matrix6d jac_left_composition_wrt_second(const trfv& tau) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward motion action with respect to the motion, equal to d(tau * p) / dDeltatauE,
     * (DeltatauE plus tau) * p = tau * p + J * DeltatauE */
    Eigen::Matrix36d jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the left jacobian of the backward motion action with respect to the motion, equal to d(tau / p) / dDeltatauE,
     * (DeltatauE plus tau) / p = tau / p + J * DeltatauE */
    Eigen::Matrix36d jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(tau)) / dDeltatau, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(Deltatau plus tau) ~= Log(tau) + (J * Deltatau) */
    Eigen::Matrix6d jac_left_log() const;

    /**< returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
     * equal to d(tau1 plus tau2) / dDeltatau1E.
     * (tau1 + Deltatau1E) plus tau2 = J * Deltatau1E plus (tau1 plus tau2)  */
    Eigen::Matrix6d jac_left_plus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left plus operator with respect to the group object (second element),
     * equal to d(tau1 plus tau2) / dDeltatau2E.
     * tau1 plus (Deltatau2E plus tau2) = J * Deltatau2E plus (tau1 plus tau2) */
    Eigen::Matrix6d jac_left_plus_wrt_second(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the first element,
     * equal to d(tau2 minus tau1) / dDeltatau2E.
     * (Deltatau2E plus tau2) minus tau1 = (tau2 minus tau1) + J * Deltatau2E */
    Eigen::Matrix6d jac_left_minus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the second element,
     * equal to d(tau2 minus tau1) / dDeltatau1E.
     * tau2 minus (Deltatau1E plus tau1) = (tau2 minus tau1) + J * Deltatau1E */
    Eigen::Matrix6d jac_left_minus_wrt_second(const trfv& tau) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adtau | xi) / dDeltatauE,
     * Ad(DeltatauE plus tau) | xi = Adtau | xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adtau % xi) / dDeltatauE,
     * Ad(DeltatauE plus tau) % xi = Adtau % xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Exponential (Right and Left) Jacobians ===== ===== */
    /**< returns a [3x3] matrix which is part of the right and left jacobians as well as their inverses. Should be private. */
    Eigen::Matrix3d jac_block() const;
    /**< returns the jacobian [6x6] of the exponential of the transform vector with respect to the transform vector itself,
     * also known as the right Jacobian of SE(3), equal to d(Exp(tau)) / dDeltatau.
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * exp(tau + Deltatau) ~= exp(tau) plus (JR * Deltatau) */
    Eigen::Matrix6d jac_right() const;
    /**< returns the inverse jacobian [6x6] of the exponential of the transform vector with respect to the transform vector
     * itself, also known as the inverse of the right Jacobian of SE(3).
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * tau + jacRinv * Deltatau ~= log[exp(tau) plus Deltatau] */
    Eigen::Matrix6d jac_right_inv() const;
    /**< returns the jacobian [6x6] of the exponential of the transform vector with respect to the transform vector itself,
     * also known as the left Jacobian of SE(3), equal to d(Exp(tau)) / dDeltatau.
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * exp(tau + Deltatau) ~= (JL * Deltatau) plus_left exp(tau) */
    Eigen::Matrix6d jac_left() const;
    /**< returns the inverse jacobian [6x6] of the exponential of the transform vector with respect to the transform vector
     * itself, also known as the inverse of the left Jacobian of SE(3).
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * tau + JLinv * Deltatau ~= log[Deltatau plus_left exp(tau)] */
    Eigen::Matrix6d jac_left_inv() const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward motion action with respect to the point, equal to d(tau * p) / dDeltap,
     * tau * (p + Deltap) = tau * p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_point() const;
    /**< returns the jacobian of the backward motion action with respect to the point, equal to d(tau / p) / dDeltap,
     * tau / (p + Deltap) = tau / p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_point() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adtau | xi) / dDeltaxi,
     * Adtau | (xi + Deltaxi) = Adtau | xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adtau % xi) / dDeltaxi,
     * Adtau % (xi + Deltaxi) = Adtau % xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x6] of a forward transformation with respect to the transform vector, equal to d(tau * p)/dtau.
     * The forward transformation is NOT linear on the transform vector, so
     * tau * p != d(tau * p)/dtau |tau*p * tau
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * (tau + Delta tau) * p ~= tau * p + d(tau * p)/dtau |tau*p * Delta tau
     * Note that the increment is Delta tau. */
    Eigen::Matrix36d jac_euclidean_forward_motion_wrt_trfv(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x6] of a forward transformation with respect to the transform vector, equal to d(tau * p)/dtau.
     * The forward transformation is NOT linear on the transform vector, so
     * tau * p != d(tau * p)/dtau |tau*p * tau
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * (tau + Delta tau) * p ~= tau * p + d(tau * p)/dtau |tau*p * Delta tau
     * Note that the increment is Delta tau. */
    Eigen::Matrix36d jac_euclidean_forward_motion_wrt_trfv_bis(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d[tau / p]/dtau.
     * The backward transformation is NOT linear on the transform vector, so
     * tau / p != d(tau / p)/dtau |tau/p * tau
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
     * Note that the increment is Delta tau. */
    Eigen::Matrix36d jac_euclidean_backward_motion_wrt_trfv_zero(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d[tau / p]/dtau.
     * The backward transformation is NOT linear on the transform vector, so
     * tau / p != d(tau / p)/dtau |tau/p * tau
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
     * Note that the increment is Delta tau. */
    Eigen::Matrix36d jac_euclidean_backward_motion_wrt_trfv(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x6] of a backward transformation with respect to the transform vector, equal to d(tau / p)/dtau.
     * The backward transformation is NOT linear on the transform vector, so
     * tau / p != d(tau / p)/dtau |tau/p * tau
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * (tau + Delta tau) / p ~= tau / p + d(tau / p)/dtau |tau/p * Delta tau
     * Note that the increment is Delta tau. */
    Eigen::Matrix36d jac_euclidean_backward_motion_wrt_trfv_bis(const Eigen::Vector3d& p) const;

    /**< returns the jacobian [3x3] of the translation vector of a forward motion with respect to the rotation vector,
     * equal to d(exp(tau).get_T)/drv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau + Deltatau).get_T() ~= exp(tau).get_T + d(exp(tau).get_T())/drv |[exp(tau)] * Deltarv.
     * Note that the jacobian is evaluated at |[exp(tau)]. */
    Eigen::Matrix3d jac_euclidean_forward_motion_T_wrt_rotv() const;
    /**< returns the jacobian [3x3] of the translation vector of a backward motion with respect to the rotation vector,
     * equal to d(exp(tau.inv()).get_T)/drv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau.inv() + Deltatau).get_T() ~= exp(tau.inv()).get_T + d(exp(tau.inv()).get_T())/drv |[exp(tau.inv())] * Deltarv.
     * Note that the jacobian is evaluated at |[exp(tau.inv())]. */
    Eigen::Matrix3d jac_euclidean_backward_motion_T_wrt_rotv() const;

    /**< returns the jacobian [3x3] of a forward motion with respect to the rotation vector, equal to d(exp(tau) * p)/drv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/drv |[exp(tau) * p] * Deltarv.
     * Note that the jacobian is evaluated at |[exp(tau)*p].*/
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_rotv(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x3] of a backward motion with respect to the rotation vector, equal to d(exp(tau) / p)/drv.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) / p)/drv |[exp(tau)/p] * Deltarv.
     * Note that the jacobian is evaluated at |[exp(tau)/p].*/
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_rotv(const Eigen::Vector3d& p) const;

    /**< returns the jacobian [3x3] of a forward motion with respect to the translation vector s, equal to d(exp(tau) * p)/ds.
     * As the forward transformation IS linear on the translation vector, the following expression is true:
     * exp(tau + Deltatau) * p ==  exp(tau) * p + d(exp(tau) * p)/ds |[exp(tau)*p] * Deltas.
     * Note that the jacobian is evaluated at |[exp(tau)*p]. */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_s(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x3] of a backward motion with respect to the translation vector s, equal to d(exp(tau) / p)/ds.
     * As the forward transformation IS linear on the translation vector, the following expression is true:
     * exp(tau + Deltatau) / p ==  exp(tau) / p + d(exp(tau) * p)/ds |[exp(tau)/p] * Deltas.
     * Note that the jacobian is evaluated at |[exp(tau)/p]. */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_s(const Eigen::Vector3d& p) const;

    /**< returns the jacobian [3x6] of a forward motion with respect to the transform vector, equal to d(exp(tau) * p)/dtau.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau + Deltatau) * p ~=  exp(tau) * p + d(exp(tau) * p)/dtau |[exp(tau)*p] * Deltatau.
     * Note that the jacobian is evaluated at |[exp(tau)*p]. */
    Eigen::Matrix36d jac_euclidean_forward_motion_wrt_trfv_tri(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x6] of a backward motion with respect to the transform vector, equal to d(exp(tau) / p)/dtau.
     * The 1st order Taylor approximation is valid for very small rotation vector changes:
     * exp(tau + Deltatau) / p ~=  exp(tau) / p + d(exp(tau) / p)/dtau |[exp(tau)/p] * Deltatau.
     * Note that the jacobian is evaluated at |[exp(tau)/p]. */
    Eigen::Matrix36d jac_euclidean_backward_motion_wrt_trfv_tri(const Eigen::Vector3d& p) const;

    /**< ===== ===== Getters ===== ===== */
    /**< return translation vector */
    Eigen::Vector3d get_T() const;
    /**< return rotation vector object */
    ang::rotv get_rotv() const;
    /**< return d vector */
    Eigen::Vector3d get_s() const;
    
    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const;
    /**< return rotation vector of inverse or opposite transformation */
    ang::rotv get_inverse_rotv() const;

    /**< set the vector */
    void set(const Eigen::Vector6d& op2);
    /**< set the vector */
    void set(Eigen::Vector6d&& op2);
}; // closes class trfv

/**< adds the input transform vector object to the stream (can not rely on Vector6d << operator as auxiliary.h file not included) */
ANG_API inline std::ostream& operator <<(std::ostream & out_str, const trfv& t) {
    out_str << t;
    return out_str;
}
    
} // closes namespace ang

#endif

