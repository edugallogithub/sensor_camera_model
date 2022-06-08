#ifndef ANG_REPRESENTATION_DUAL
#define ANG_REPRESENTATION_DUAL

#include "../ang.h"
#include "../rotate/rodrigues.h"

/*
 * This file contains the "unit dual quaternion" representation of a transformation.
 */

namespace ang {
    class rodrigues;
    class speu_rodrigues;
    class speu_dcm;
    class homogeneous;
    class trfv;
    class screw;
    class se3_tangent;
    class se3_tangent_dual;
    class dual_quat;

// CLASS UNIT DUAL QUATERNION
// ==========================
// ==========================

class ANG_API dual {
private:
    /**< unit quaternion representing rotation */
    ang::rodrigues _qr;
    /**< dual quaternion representing translation */
    ang::quat _qd;

    /**< normalize */
    void normalize();
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< default constructor (not initialized) */
    dual() = default;
    /**< constructor based on rodrigues parameters and quaternion */
    dual(const ang::rodrigues& qr, const ang::quat& qd) : _qr(qr), _qd(qd) {}
    /**< constructor based on rodrigues parameters and translation vector */
    dual(const ang::rodrigues&, const Eigen::Vector3d&);
    /**< constructor based on rotation vector and translation vector */
    dual(const ang::rotv&, const Eigen::Vector3d&);
    /**< constructor based on special Euclidean (rodrigues) */
    explicit dual(const ang::speu_rodrigues&);
    /**< constructor based on special Euclidean (dcm) */
    explicit dual(const ang::speu_dcm&);
    /**< constructor based on homogeneous */
    explicit dual(const ang::homogeneous&);
    /**< constructor based on transform vector */
    explicit dual(const ang::trfv&);
    /**< constructor based on screw */
    explicit dual(const ang::screw&);
    /**< copy constructor */
    dual(const dual&) = default;
    /**< destructor */
    ~dual() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    dual(ang::dual&&) = default;
    /**< move constructor based on rodrigues parameters and quaternion */
    dual(ang::rodrigues&& qr, ang::quat&& qd) : _qr(qr), _qd(qd) {}
    /**< move constructor based on rodrigues parameters and translation vector */
    dual(ang::rodrigues&&, Eigen::Vector3d&&);
    /**< move constructor based on special Euclidean (rodrigues) */
    explicit dual(ang::speu_rodrigues&&);
    /**< move constructor based on special Euclidean (dcm) */
    explicit dual(ang::speu_dcm&&);
    /**< move constructor based on homogeneous */
    explicit dual(ang::homogeneous&&);
    /**< move constructor based on transform vector */
    explicit dual(ang::trfv&&);
    /**< move constructor based on screw */
    explicit dual(ang::screw&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    dual& operator=(const ang::dual&) = default;
    /**< assignment operator = based on special Euclidean (rodrigues) */
    dual& operator=(const ang::speu_rodrigues&);
    /**< assignment operator = based on special Euclidean (dcm) */
    dual& operator=(const ang::speu_dcm&);
    /**< assignment operator = based on homogeneous */
    dual& operator=(const ang::homogeneous&);
    /**< assignment operator = based on transform vector */
    dual& operator=(const ang::trfv&);
    /**< assignment operator = based on screw */
    dual& operator=(const ang::screw&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    dual& operator=(ang::dual&&) = default;
    /**< move assignment operator = based on special Euclidean (rodrigues) */
    dual& operator=(ang::speu_rodrigues&&);
    /**< move assignment operator = based on special Euclidean (dcm) */
    dual& operator=(ang::speu_dcm&&);
    /**< move assignment operator = based on homogeneous */
    dual& operator=(ang::homogeneous&&);
    /**< move assignment operator = based on transform vector */
    dual& operator=(ang::trfv&&);
    /**< move assignment operator = based on screw */
    dual& operator=(screw&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    dual operator*(const dual& op2) const;
    /**< overloaded operator / (backward combination of transformations) */
    dual operator/(const dual& op2) const;
    /**< returns inverse or opposite transformation */
    dual inverse() const {return dual(rodrigues(_qr().adjoint()), _qd.adjoint());}
    /**< executes object transformation fraction (t < 1) or a multiple (t > 1) of times.
    * Returns exponential map of the power function applied to the object screw logarithmic map. */
    dual pow(const double& t) const;
    /**< screw linear interpolation, returns Z0 for t=0 and Z1 for t=1 */
    static dual sclerp(const dual& Z0, const dual& Z1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    dual plus_right(const trfv&) const;
    /**< right plus operator (input rotation located in local tangent space) */
    dual plus_right(const screw&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    dual plus_left(const trfv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    dual plus_left(const screw&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d& vecin) const {return _qr * vecin + quat::convert_4dto3d(_qd * _qr().adjoint()) - quat::convert_4dto3d(_qr() * _qd.adjoint());}
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d& vecin) const {return _qr / vecin + quat::convert_4dto3d(_qd.adjoint() * _qr()) - quat::convert_4dto3d(_qr().adjoint() * _qd);}
    /**< overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator^(const Eigen::Vector3d& vecin) const {return _qr * vecin;}
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d& vecin) const {return _qr / vecin;}
    /**< right minus operator (output transformation located in local tangent space) */
    trfv minus_right_trfv(const dual&) const;
    /**< right minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const dual&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const dual&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const dual&) const;
    /**< logarithmic map that returns the screw */
    ang::screw log_map_screw() const;
    /**< logarithmic map that returns the transform vector */
    ang::trfv log_map_trfv() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< overloaded operator | (forward adjoint) */
    ang::se3_tangent_dual operator|(const ang::se3_tangent_dual& xi_dual) const;
    /**< overloaded operator | (forward adjoint) */
    ang::se3_tangent operator|(const ang::se3_tangent& xi) const;
    /**< returns forward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_forward() const;
    /**< overloaded operator % (backward adjoint) */
    ang::se3_tangent_dual operator%(const ang::se3_tangent_dual& xi_dual) const;
    /**< overloaded operator % (backward adjoint) */
    ang::se3_tangent operator%(const ang::se3_tangent& xi) const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_backward() const;

    /**< ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */
    /**< obtains the body twist or motion velocity from the unit dual quaternion transformation and its time derivative */
    ang::se3_tangent dot2xibody(const dual_quat& dualdot) const;
    /**< obtains the unit dual quaternion transformation derivative with time based on the unit dual quaternion
    transformation and the body twist or motion velocity. */
    ang::dual_quat xibody2dot(const ang::se3_tangent& xi_body_mrps) const;
    /**< obtains the space twist or motion velocity from the unit dual quaternion transformation and its time derivative */
    ang::se3_tangent dot2xispace(const dual_quat& dualdot) const;
    /**< obtains the unit dual quaternions transformation derivative with time based on the unit dual quaternion
    transformation and the space twist or motion velocity. */
    ang::dual_quat xispace2dot(const ang::se3_tangent& xi_space_mrps) const;

    /**< ===== ===== Linear Algebra ===== ===== */

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the eight components of the
     * unit dual quaternion in dual form and returns them in vector form. */
    static Eigen::Vector8d wedge(const ang::dual& zeta);
    /**< although the hat operator usually applies to the tangent space, here it takes the eight components of the
     * unit dual quaternion in vector form, and returns them in unit dual quaternion form (object). */
    static ang::dual hat(const Eigen::Vector8d& v);

    /**< ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse motion, equal to d(zeta^-1)/dDeltatauB, which coincides with the negative of the adjoint.
     * (zeta plus DeltatauB)^-1 = zeta^-1 + J * DeltatauB */
    Eigen::Matrix6d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(zeta1 * zeta2) / dDeltatauB1, which coincides with the inverse of the 2nd motion adjoint.
     * (zeta1 plus DeltatauB1) * zeta2 = zeta1 * zeta2 plus J * DeltatauB1 */
    Eigen::Matrix6d jac_right_composition_wrt_first(const dual& zeta) const {return zeta.adjoint_matrix_backward();}
    /**< returns the right jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(zeta1 * zeta2) / dDeltatauB2, which as a matter of fact is the identity.
     * zeta1 * (zeta2 plus DeltatauB2) = zeta1 * zeta2 plus J * DeltatauB2 */
    Eigen::Matrix6d jac_right_composition_wrt_second(const dual& zeta) const {return Eigen::Matrix6d::Identity();}

    /**< returns the right jacobian of the forward motion action with respect to the motion, equal to d(zeta * p) / dDeltatauB,
     * (zeta plus DeltatauB) * p = zeta * p + J * DeltatauB */
    Eigen::Matrix36d jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the right jacobian of the backward motion action with respect to the motion, equal to d(zeta / p) / dDeltatauB,
     * (zeta plus DeltatauB) / p = zeta / p + J * DeltatauB */
    Eigen::Matrix36d jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(zeta)) / dDeltatau, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(zeta plus Deltatau) ~= Log(zeta) + (J * Deltatau) */
    Eigen::Matrix6d jac_right_log() const;

    /**< returns the right jacobian of the motion right plus operator with respect to the group object (first element),
     * equal to d(zeta1 plus tau2) / dDeltatau1B.
     * (zeta1 plus Deltatau1B) plus tau2 = (zeta1 plus tau2) plus J * Deltatau1B */
    Eigen::Matrix6d jac_right_plus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
     * equal to d(zeta1 plus tau2) / dDeltatau2B.
     * zeta1 plus (tau2 + Deltatau2B) = (zeta1 plus tau2) plus J * Deltatau2B */
    Eigen::Matrix6d jac_right_plus_wrt_second(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the first element,
     * equal to d(zeta2 minus zeta1) / dDeltatau2B.
     * (zeta2 plus Deltatau2B) minus zeta1 = (zeta2 minus zeta1) + J * Deltatau2B */
    Eigen::Matrix6d jac_right_minus_wrt_first(const dual& Z) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the second element,
     * equal to d(zeta2 minus zeta1) / dDeltatau1B.
     * zeta2 minus(zeta1 plus Deltatau1B) = (zeta2 minus zeta1) + J * Deltatau1B */
    Eigen::Matrix6d jac_right_minus_wrt_second(const dual& Z) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdZ | xi) / dDeltatauB,
     * Ad(Z plus DeltatauB) | xi = AdZ | xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdZ % xi) / dDeltatauB,
     * Ad(Z plus DeltatauB) % xi = AdZ % xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse motion, equal to d(z^-1)/dDeltatauE, which coincides with the negative of the inverse adjoint.
     * (DeltatauE plus z)^-1 = J * DeltatauE plus z^-1 */
    Eigen::Matrix6d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(z1 * z2) / dDeltatauE1, which as a matter of fact is the identity.
     * (DeltatauE1 plus z1) * z2 = J * DeltatauE1 plus z1 * z2 */
    Eigen::Matrix6d jac_left_composition_wrt_first(const dual& z) const {return Eigen::Matrix6d::Identity();}
    /**< returns the left jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(z1 * z2) / dDeltatauE2, which coincides with the adjoint of the 1st motion.
     * z1 * (DeltatauE2 plus z2) = J * DeltatauE2 plus z1 * z2 */
    Eigen::Matrix6d jac_left_composition_wrt_second(const dual& z) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward motion action with respect to the motion, equal to d(z * p) / dDeltatauE,
     * (DeltatauE plus z) * p = z * p + J * DeltatauE */
    Eigen::Matrix36d jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the left jacobian of the backward motion action with respect to the motion, equal to d(z / p) / dDeltatauE,
     * (DeltatauE plus z) / p = z / p + J * DeltatauE */
    Eigen::Matrix36d jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(z)) / dDeltatau, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(Deltatau plus z) ~= Log(z) + (J * Deltatau) */
    Eigen::Matrix6d jac_left_log() const;

    /**< returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
     * equal to d(tau1 plus z2) / dDeltatau1E.
     * (tau1 + Deltatau1E) plus z2 = J * Deltatau1E plus (tau1 plus z2)  */
    Eigen::Matrix6d jac_left_plus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left plus operator with respect to the group object (second element),
     * equal to d(tau1 plus z2) / dDeltatau2E.
     * tau1 plus (Deltatau2E plus z2) = J * Deltatau2E plus (tau1 plus z2) */
    Eigen::Matrix6d jac_left_plus_wrt_second(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the first element,
     * equal to d(z2 minus z1) / dDeltatau2E.
     * (Deltatau2E plus z2) minus z1 = (z2 minus z1) + J * Deltatau2E */
    Eigen::Matrix6d jac_left_minus_wrt_first(const dual& z) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the second element,
     * equal to d(z2 minus z1) / dDeltatau1E.
     * z2 minus (Deltatau1E plus z1) = (z2 minus z1) + J * Deltatau1E */
    Eigen::Matrix6d jac_left_minus_wrt_second(const dual& z) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(Adz | xi) / dDeltatauE,
     * Ad(DeltatauE plus z) | xi = Adz | xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(Adz % xi) / dDeltatauE,
     * Ad(DeltatauE plus z) % xi = Adz % xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward motion action with respect to the point, equal to d(zeta * p) / dDeltap,
     * zeta * (p + Deltap) = zeta * p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_point() const;
    /**< returns the jacobian of the backward motion action with respect to the point, equal to d(zeta / p) / dDeltap,
       * zeta / (p + Deltap) = zeta / p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_point() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdZ | xi) / dDeltaxi,
     * AdZ | (xi + Deltaxi) = AdZ | xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdZ % xi) / dDeltaxi,
     * AdZ % (xi + Deltaxi) = AdZ % xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< ===== ===== Setters and  Getters ===== ===== */
    /**< return translation vector */
    Eigen::Vector3d get_T() const;
    /**< return rotation vector object */
    ang::rotv get_rotv() const;
    /**< return rodrigues rotation */
    const ang::rodrigues& get_rodrigues() const {return _qr;}

    /**< set translation component maintaining the rotation component based on the input translation vector.
     * Internally, this maintains the existing real quaternion (_qr) and updates the dual quaternion (_qd) */
    void set_T(const Eigen::Vector3d& T);
    /**< set rotation component maintaining the translation component based on the input rotation vector.
     * Internally, this modifies both the real (_qr) and dual (_qd) quaternions. */
    void set_rotv(const ang::rotv& r);
    /**< set rotation component maintaining the translation component based on the input unit quaternion.
     * Internally, this modifies both the real (_qr) and dual (_qd) quaternions. */
    void set_rodrigues(const ang::rodrigues& q);




    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const;
    /**< return rotation vector of inverse or opposite transformation */
    ang::rotv get_inverse_rotv() const;

    /**< get real quaternion representing rotation */
    const ang::rodrigues& get_qr() const {return _qr;}
    /**< get dual quaternion representing translation */
    const ang::quat& get_qd() const {return _qd;}



}; // closes class dual

} // closes namespace ang

#endif



























