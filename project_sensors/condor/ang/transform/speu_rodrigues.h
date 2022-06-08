#ifndef ANG_REPRESENTATION_SPECIAL_EUCLIDEAN_RODRIGUES
#define ANG_REPRESENTATION_SPECIAL_EUCLIDEAN_RODRIGUES

#include "../ang.h"
#include "../rotate/rodrigues.h"

/*
 * This file contains the "speu_rodrigues" class that models the special Euclidean SE(3)
 * transformation composed by a rotation (based on rodrigues) and a translation.
 */

namespace ang {
    class speu_dcm;
    class homogeneous;
    class trfv;
    class screw;
    class dual;
    class se3_tangent;
    class se3_tangent_homo;

// CLASS SPECIAL EUCLIDEAN (based on RODRIGUES)
// ============================================
// ============================================

class ANG_API speu_rodrigues {
private:
    /**< rodrigues representing a rotation */
    ang::rodrigues _q;
    /**< vector representing translation */
    Eigen::Vector3d _T;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< default constructor (not initialized) */
    speu_rodrigues() = default;
    /**< constructor based on rodrigues and translation */
    speu_rodrigues(const ang::rodrigues& q, const Eigen::Vector3d& T) : _q(q), _T(T) {}
    /**< constructor based on inverse rodrigues and inverse translation.
     * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */
    speu_rodrigues(const ang::rodrigues& q_inv, const Eigen::Vector3d& T_inv, bool not_used) : _q(q_inv.inverse()), _T(q_inv.inverse() * (-T_inv)) {}
    /**< constructor based on special Euclidean (dcm) */
    explicit speu_rodrigues(const speu_dcm&);
    /**< constructor based on homogeneous */
    explicit speu_rodrigues(const homogeneous&);
    /**< constructor based on transform vector */
    explicit speu_rodrigues(const trfv&);
    /**< constructor based on screw */
    explicit speu_rodrigues(const screw&);
    /**< constructor based on unit dual quaternion */
    explicit speu_rodrigues(const dual&);
    /**< copy constructor */
    speu_rodrigues(const speu_rodrigues&) = default;
    /**< destructor */
    virtual ~speu_rodrigues() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    speu_rodrigues(speu_rodrigues&&) = default;
    /**< move constructor based on rodrigues and translation */
    speu_rodrigues(ang::rodrigues&& q, Eigen::Vector3d&& T) : _q(q), _T(T) {}
    /**< move constructor based on inverse rodrigues and inverse translation.
     * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */
    speu_rodrigues(ang::rodrigues&& q_inv, Eigen::Vector3d&& T_inv, bool not_used) : _q(q_inv.inverse()), _T(q_inv.inverse() * (-T_inv)) {}
    /**< move constructor based on special Euclidean (dcm) */
    explicit speu_rodrigues(speu_dcm&&);
    /**< move constructor based on homogeneous */
    explicit speu_rodrigues(homogeneous&&);
    /**< move constructor based on transform vector */
    explicit speu_rodrigues(trfv&&);
    /**< move constructor based on screw */
    explicit speu_rodrigues(screw&&);
    /**< move constructor based on unit dual quaternion */
    explicit speu_rodrigues(dual&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    speu_rodrigues& operator=(const speu_rodrigues&) = default;
    /**< assignment operator = based on special Euclidean (dcm) */
    speu_rodrigues& operator=(const speu_dcm&);
    /**< assignment operator = based on homogeneous */
    speu_rodrigues& operator=(const homogeneous&);
    /**< assignment operator = based on transform vector */
    speu_rodrigues& operator=(const trfv&);
    /**< assignment operator = based on screw */
    speu_rodrigues& operator=(const screw&);
    /**< assignment operator = based on unit dual quaternion */
    speu_rodrigues& operator=(const dual&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    speu_rodrigues& operator=(speu_rodrigues&&) = default;
    /**< move assignment operator = based on special Euclidean (dcm) */
    speu_rodrigues& operator=(speu_dcm&&);
    /**< move assignment operator = based on homogeneous */
    speu_rodrigues& operator=(homogeneous&&);
    /**< move assignment operator = based on transform vector */
    speu_rodrigues& operator=(trfv&&);
    /**< move assignment operator = based on screw */
    speu_rodrigues& operator=(screw&&);
    /**< move assignment operator = based on unit dual quaternion */
    speu_rodrigues& operator=(dual&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    speu_rodrigues operator*(const speu_rodrigues& op2) const {return ang::speu_rodrigues(_q * op2._q, _q * op2._T + _T);}
    /**< overloaded operator / (backward combination of transformations) */
    speu_rodrigues operator/(const speu_rodrigues& op2) const;
    /**< returns inverse or opposite transformation */
    speu_rodrigues inverse() const {return ang::speu_rodrigues(_q.inverse(), _q.inverse() * (- _T));}
    /**< executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns exponential map of the power function applied to the object screw logarithmic map. */
    speu_rodrigues pow(const double& t) const;
    /**< screw linear interpolation, returns gq0 for t=0 and gq1 for t=1 */
    static speu_rodrigues sclerp(const speu_rodrigues& gq0, const speu_rodrigues& gq1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    speu_rodrigues plus_right(const trfv&) const;
    /**< right plus operator (input rotation located in local tangent space) */
    speu_rodrigues plus_right(const screw&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    speu_rodrigues plus_left(const trfv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    speu_rodrigues plus_left(const screw&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d& vecin) const {return _q * vecin + _T;}
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d& vecin) const {return _q / vecin - _q / _T;}
    /**< overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -.*/
    Eigen::Vector3d operator^(const Eigen::Vector3d& vecin) const {return _q * vecin;}
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d& vecin) const {return _q / vecin;}
    /**< right minus operator (output transformation located in local tangent space) */
    trfv minus_right_trfv(const speu_rodrigues&) const;
    /**< right minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const speu_rodrigues&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const speu_rodrigues&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const speu_rodrigues&) const;

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

    /**< ===== ===== Twist or Motion Velocity - Time Derivative  ===== ===== */
    /**< obtains the body twist or motion velocity from the special euclidean Rodrigues parameters and its time derivative */
    ang::se3_tangent dot2xibody(const Eigen::Vector7d& speu_rodriguesdot) const;
    /**< obtains the special euclidean Rodrigues parameters derivative with time based on the homogeneous
    transformation and the body twist or motion velocity. */
    Eigen::Vector7d xibody2dot(const ang::se3_tangent& xi_body_mrps) const;
    /**< obtains the space twist or motion velocity from the special euclidean Rodrigues parameters and its time derivative */
    ang::se3_tangent dot2xispace(const Eigen::Vector7d& speu_rodriguesdot) const;
    /**< obtains the special euclidean Rodrigues parameters derivative with time based on the homogeneous
    transformation and the space twist or motion velocity. */
    Eigen::Vector7d xispace2dot(const ang::se3_tangent& xi_space_mrps) const;

    /**< ===== ===== Linear Algebra ===== ===== */

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the seven components of the
     * special euclidean rodrigues form and returns them in vector form (quaternion first). */
    static Eigen::Vector7d wedge(const ang::speu_rodrigues& Gq);
    /**< although the hat operator usually applies to the tangent space, here it takes the seven components of the
     * special euclidean rodrigues in vector from and returns them in object form. */
    static ang::speu_rodrigues hat(const Eigen::Vector7d& v);

    /**< ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse motion, equal to d(gq^-1)/dDeltatauB, which coincides with the negative of the adjoint.
     * (gq plus DeltatauB)^-1 = gq^-1 + J * DeltatauB */
    Eigen::Matrix6d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(gq1 * gq2) / dDeltatauB1, which coincides with the inverse of the 2nd motion adjoint.
     * (gq1 plus DeltatauB1) * gq2 = gq1 * gq2 plus J * DeltatauB1 */
    Eigen::Matrix6d jac_right_composition_wrt_first(const speu_rodrigues& gq) const {return gq.adjoint_matrix_backward();}
    /**< returns the right jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(gq1 * gq2) / dDeltatauB2, which as a matter of fact is the identity.
     * gq1 * (gq2 plus DeltatauB2) = gq1 * gq2 plus J * DeltatauB2 */
    Eigen::Matrix6d jac_right_composition_wrt_second(const speu_rodrigues& gq) const {return Eigen::Matrix6d::Identity();}

    /**< returns the right jacobian of the forward motion action with respect to the motion, equal to d(gq * p) / dDeltatauB,
     * (gq plus DeltatauB) * p = gq * p + J * DeltatauB */
    Eigen::Matrix36d jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the right jacobian of the backward motion action with respect to the motion, equal to d(gq / p) / dDeltatauB,
     * (gq plus DeltatauB) / p = gq / p + J * DeltatauB */
    Eigen::Matrix36d jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(gq)) / dDeltatau, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(gq plus Deltatau) ~= Log(gq) + (J * Deltatau) */
    Eigen::Matrix6d jac_right_log() const;

    /**< returns the right jacobian of the motion right plus operator with respect to the group object (first element),
     * equal to d(gq1 plus tau2) / dDeltatau1B.
     * (gq1 plus Deltatau1B) plus tau2 = (gq1 plus tau2) plus J * Deltatau1B */
    Eigen::Matrix6d jac_right_plus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
     * equal to d(gq1 plus tau2) / dDeltatau2B.
     * gq1 plus (tau2 + Deltatau2B) = (gq1 plus tau2) plus J * Deltatau2B */
    Eigen::Matrix6d jac_right_plus_wrt_second(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the first element,
     * equal to d(gq2 minus gq1) / dDeltatau2B.
     * (gq2 plus Deltatau2B) minus gq1 = (gq2 minus gq1) + J * Deltatau2B */
    Eigen::Matrix6d jac_right_minus_wrt_first(const speu_rodrigues& gq) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the second element,
     * equal to d(gq2 minus gq1) / dDeltatau1B.
     * gq2 minus(gq1 plus Deltatau1B) = (gq2 minus gq1) + J * Deltatau1B */
    Eigen::Matrix6d jac_right_minus_wrt_second(const speu_rodrigues& gq) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(Adgq | xi) / dDeltatauB,
     * Ad(gq plus DeltatauB) | xi = Adgq | xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(Adgq % xi) / dDeltatauB,
     * Ad(gq plus DeltatauB) % xi = Adgq % xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse motion, equal to d(gq^-1)/dDeltatauE, which coincides with the negative of the inverse adjoint.
     * (DeltatauE plus gq)^-1 = J * DeltatauE plus gq^-1 */
    Eigen::Matrix6d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(gq1 * gq2) / dDeltatauE1, which as a matter of fact is the identity.
     * (DeltatauE1 plus gq1) * gq2 = J * DeltatauE1 plus gq1 * gq2 */
    Eigen::Matrix6d jac_left_composition_wrt_first(const speu_rodrigues& gq) const {return Eigen::Matrix6d::Identity();}
    /**< returns the left jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(gq1 * gq2) / dDeltatauE2, which coincides with the adjoint of the 1st motion.
     * gq1 * (DeltatauE2 plus gq2) = J * DeltatauE2 plus gq1 * gq2 */
    Eigen::Matrix6d jac_left_composition_wrt_second(const speu_rodrigues& gq) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward motion action with respect to the motion, equal to d(gq * p) / dDeltatauE,
     * (DeltatauE plus gq) * p = gq * p + J * DeltatauE */
    Eigen::Matrix36d jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the left jacobian of the backward motion action with respect to the motion, equal to d(gq / p) / dDeltatauE,
     * (DeltatauE plus gq) / p = gq / p + J * DeltatauE */
    Eigen::Matrix36d jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(gq)) / dDeltatau, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(Deltatau plus gq) ~= Log(gq) + (J * Deltatau) */
    Eigen::Matrix6d jac_left_log() const;

    /**< returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
     * equal to d(tau1 plus gq2) / dDeltatau1E.
     * (tau1 + Deltatau1E) plus gq2 = J * Deltatau1E plus (tau1 plus gq2)  */
    Eigen::Matrix6d jac_left_plus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left plus operator with respect to the group object (second element),
     * equal to d(tau1 plus gq2) / dDeltatau2E.
     * tau1 plus (Deltatau2E plus gq2) = J * Deltatau2E plus (tau1 plus gq2) */
    Eigen::Matrix6d jac_left_plus_wrt_second(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the first element,
     * equal to d(gq2 minus gq1) / dDeltatau2E.
     * (Deltatau2E plus gq2) minus gq1 = (gq2 minus gq1) + J * Deltatau2E */
    Eigen::Matrix6d jac_left_minus_wrt_first(const speu_rodrigues& gq) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the second element,
     * equal to d(gq2 minus gq1) / dDeltatau1E.
     * gq2 minus (Deltatau1E plus gq1) = (gq2 minus gq1) + J * Deltatau1E */
    Eigen::Matrix6d jac_left_minus_wrt_second(const speu_rodrigues& gq) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adgq | xi) / dDeltatauE,
     * Ad(DeltatauE plus gq) | xi = Adgq | xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adgq % xi) / dDeltatauE,
     * Ad(DeltatauE plus gq) % xi = Adgq % xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward motion action with respect to the point, equal to d(gq * p) / dDeltap,
     * gq * (p + Deltap) = gq * p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_point() const;
    /**< returns the jacobian of the backward motion action with respect to the point, equal to d(gq / p) / dDeltap,
     * gq / (p + Deltap) = gq / p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_point() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adgq | xi) / dDeltaxi,
     * Adgq | (xi + Deltaxi) = Adgq | xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adgq % xi) / dDeltaxi,
     * Adgq % (xi + Deltaxi) = Adgq % xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x7] of a forward transformation with respect to the affine unit quaternion, equal to d(Gq * p)/dGq.
     * The 1st order Taylor approximation is valid for very small affine unit quaternion changes:
     * [Gq * exp(Delta tau)] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[Gq * exp(Delta tau)] - Gq}
     * [exp(Delta tau) * Gq] * p ~= Gq * p + d(Gq * p)/dGq |Gq*p * {[exp(Delta tau) * Gq] - Gq}
     * Note that the jacobian is evaluated at |Gq*p.
     * Note that the increment is {[Gq * exp(Delta tau)] - Gq} or [(exp(Delta tau) * Gq) - Gq], and not exp(Delta tau) or Delta Gq. */
    Eigen::Matrix37d jac_euclidean_forward_motion_wrt_speu_rodrigues(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x7] of a backward transformation with respect to the affine unit quaternion, equal to d(Gq / p)/dGq.
     * The 1st order Taylor approximation is valid for very small affine unit quaternion changes):
     * [Gq * exp(Delta tau)] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[Gq * exp(Delta tau)] - Gq}
     * [exp(Delta tau) * Gq] / p ~= Gq / p + d(Gq / p)/dGq |Gq/p * {[exp(Delta tau) * Gq] - Gq}
     * Note that the jacobian is evaluated at |Gq/p.
     * Note that the increment is {[Gq * exp(Delta tau)] - Gq} or [(exp(Delta tau) * Gq) - Gq], and not exp(Delta tau) or Delta Gq. */
    Eigen::Matrix37d jac_euclidean_backward_motion_wrt_speu_rodrigues(const Eigen::Vector3d& p) const;

    /**< ===== ===== Getters and Setters ===== ===== */
    /**< return translation vector */
    const Eigen::Vector3d& get_T() const {return _T;}
    /**< return rodrigues rotation */
    const ang::rodrigues& get_rodrigues() const {return _q;}
    /**< return rotation vector object */
    ang::rotv get_rotv() const;
    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const {return _q.inverse() * (- _T);}
    /**< return rodrigues rotation of inverse or opposite transformation */
    ang::rodrigues get_inverse_rodrigues() const {return _q.inverse();}
    /**< modify the object based on the input rotation and translation */
    void set(const ang::rodrigues& q, const Eigen::Vector3d& T);
    /**< modify the object based on the input inverse rotation and inverse translation.
     * Boolean is meaningless, only purpose is to remind yourself that inputs are inverses. */
    void set(const ang::rodrigues& q_inv, const Eigen::Vector3d& T_inv, bool);
}; // closes class speu_rodrigues

/**< adds the input special euclidean object to the stream (can not rely on Vector3d << operator as auxiliary.h file not included) */
ANG_API inline std::ostream& operator <<(std::ostream & out_str, const speu_rodrigues& op2) {
    out_str << op2.get_rodrigues()() << " | " << op2.get_T();
    return out_str;
}
    
} // closes namespace ang

#endif

