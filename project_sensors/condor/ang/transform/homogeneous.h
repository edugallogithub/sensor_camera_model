#ifndef ANG_REPRESENTATION_HOMOGENEOUS
#define ANG_REPRESENTATION_HOMOGENEOUS

#include "../ang.h"
#include "../rotate/dcm.h"

/*
 * This file contains the "homogeneous" representation of a transformation.
 */

namespace ang {
    class speu_rodrigues;
    class speu_dcm;
    class trfv;
    class screw;
    class dual;
    class se3_tangent;
    class se3_tangent_homo;

// CLASS HOMOGENEOUS
// =================
// =================

class ANG_API homogeneous : private Eigen::Matrix4d {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not initialized) */
    homogeneous() = default;
    /**< constructor based on 4x4 matrix */
    explicit homogeneous(const Eigen::Matrix4d& Omat);
    /**< constructor based on direction cosine matrix and translation vector */
    homogeneous(const ang::dcm&, const Eigen::Vector3d&);
    /**< constructor based on special Euclidean (rodrigues) */
    explicit homogeneous(const ang::speu_rodrigues&);
    /**< constructor based on special Euclidean (dcm) */
    explicit homogeneous(const ang::speu_dcm&);
    /**< constructor based on transform vector */
    explicit homogeneous(const ang::trfv&);
    /**< constructor based on screw */
    explicit homogeneous(const ang::screw&);
    /**< constructor based on unit dual quaternion */
    explicit homogeneous(const ang::dual&);
    /**< copy constructor */
    homogeneous(const homogeneous&) = default;
    /**< destructor */
    ~homogeneous() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    homogeneous(homogeneous&&) = default;
    /**< move constructor based on 4x4 matrix */
    explicit homogeneous(Eigen::Matrix4d&& Omat);
    /**< move constructor based on direction cosine matrix and translation vector */
    homogeneous(ang::dcm&&, Eigen::Vector3d&&);
    /**< move constructor based on special Euclidean (rodrigues) */
    explicit homogeneous(ang::speu_rodrigues&&);
    /**< move constructor based on special Euclidean (dcm) */
    explicit homogeneous(ang::speu_dcm&&);
    /**< move constructor based on transform vector */
    explicit homogeneous(ang::trfv&&);
    /**< move constructor based on screw */
    explicit homogeneous(ang::screw&&);
    /**< move constructor based on unit dual quaternion */
    explicit homogeneous(ang::dual&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    homogeneous& operator=(const homogeneous&) = default;
    /**< assignment operator based on 4x4 matrix */
    homogeneous& operator=(const Eigen::Matrix4d&);
    /**< assignment operator = based on special Euclidean (rodrigues) */
    homogeneous& operator=(const ang::speu_rodrigues&);
    /**< assignment operator = based on special Euclidean (dcm) */
    homogeneous& operator=(const ang::speu_dcm&);
    /**< assignment operator = based on transform vector */
    homogeneous& operator=(const ang::trfv&);
    /**< assignment operator = based on screw */
    homogeneous& operator=(const ang::screw&);
    /**< assignment operator = based on unit dual quaternion */
    homogeneous& operator=(const ang::dual&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    homogeneous& operator=(homogeneous&&) = default;
    /**< move assignment operator based on 4x4 matrix */
    homogeneous& operator=(Eigen::Matrix4d&&);
    /**< move assignment operator = based on special Euclidean (rodrigues) */
    homogeneous& operator=(ang::speu_rodrigues&&);
    /**< move assignment operator = based on special Euclidean (dcm) */
    homogeneous& operator=(ang::speu_dcm&&);
    /**< move assignment operator = based on transform vector */
    homogeneous& operator=(ang::trfv&&);
    /**< move assignment operator = based on screw */
    homogeneous& operator=(ang::screw&&);
    /**< move assignment operator = based on unit dual quaternion */
    homogeneous& operator=(ang::dual&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    homogeneous operator*(const homogeneous& op2) const;
    /**< overloaded operator / (backward combination of transformations) */
    homogeneous operator/(const homogeneous& op2) const;
    /**< returns inverse or opposite transformation */
    homogeneous inverse() const;
    /**< executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns exponential map of the power function applied to the object screw logarithmic map. */
    homogeneous pow(const double& t) const;
    /**< screw linear interpolation, returns M0 for t=0 and M1 for t=1 */
    static homogeneous sclerp(const homogeneous& M0, const homogeneous& M1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    homogeneous plus_right(const trfv&) const;
    /**< right plus operator (input rotation located in local tangent space) */
    homogeneous plus_right(const screw&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    homogeneous plus_left(const trfv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    homogeneous plus_left(const screw&) const;

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;
    /**< overloaded operator ^ (forward transformation of vector, not point). WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator^(const Eigen::Vector3d& vecin) const {return this->topLeftCorner(3,3) * vecin;}
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d& vecin) const {return this->topLeftCorner(3,3).transpose() *  vecin;}
    /**< overloaded operator * (forward transformation of homogeneous point OR vector) */
    Eigen::Vector4d operator*(const Eigen::Vector4d& vecin) const {return static_cast<const Eigen::Matrix4d&>(*this) * vecin;}
    /**< overloaded operator / (backward transformation of homogeneous point OR vector) */
    Eigen::Vector4d operator/(const Eigen::Vector4d& vecin) const {return this->inverse() * vecin;}
    /**< right minus operator (output transformation located in local tangent space) */
    trfv minus_right_trfv(const homogeneous&) const;
    /**< right minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const homogeneous&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const homogeneous&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const homogeneous&) const;
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
    /**< obtains the body twist or motion velocity from the homogeneous transformation and its time derivative */
    ang::se3_tangent dot2xibody(const Eigen::Matrix4d& homogeneousdot) const;
    /**< obtains the homogeneous transformation derivative with time based on the homogeneous
    transformation and the body twist or motion velocity. */
    Eigen::Matrix4d xibody2dot(const ang::se3_tangent& xi_body_mrps) const;
    /**< obtains the space twist or motion velocity from the homogeneous transformation and its time derivative */
    ang::se3_tangent dot2xispace(const Eigen::Matrix4d& homogeneousdot) const;
    /**< obtains the homogeneous transformation derivative with time based on the homogeneous
    transformation and the space twist or motion velocity. */
    Eigen::Matrix4d xispace2dot(const ang::se3_tangent& xi_space_mrps) const;

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get 4x4 matrix */
    Eigen::Matrix4d& operator()() {return *this;}
    const Eigen::Matrix4d& operator()() const {return *this;}
    /**< get 4x4 matrix */
    Eigen::Matrix4d& get() {return *this;}
    const Eigen::Matrix4d& get() const {return *this;}

    /**< ===== ===== Hat and Wedge Operators ===== ===== */
    /**< although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
     * homogeneous matrix in matrix form and returns them in vector form, first those of the rotation matrix (ordered
     * by row), and last the translation vector. */
    static Eigen::Vector12d wedge(const ang::homogeneous& M);
    /**< although the wedge operator usually applies to the tangent space, here it takes the twelve components of the
     * homogeneous matrix in matrix form and returns them in vector form, first those of the rotation matrix (ordered
     * by row), and last the translation vector. */
    static Eigen::Vector12d wedge(const Eigen::Matrix4d& M);
    /**< although the hat operator usually applies to the tangent space, here it takes the twelve components of the
     * homogeneous matrix in vector form, first those of the rotation matrix (ordered by row), and last the translation
     * vector. */
    static ang::homogeneous hat(const Eigen::Vector12d& v);

    /**< ===== ===== Right Transform Vector (Local Tangent Space) Jacobians ===== ===== */
    /**< returns the right jacobian of the inverse motion, equal to d(M^-1)/dDeltatauB, which coincides with the negative of the adjoint.
     * (M plus DeltatauB)^-1 = M^-1 + J * DeltatauB */
    Eigen::Matrix6d jac_right_inverse() const {return - this->adjoint_matrix_forward();}

    /**< returns the right jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(M1 * M2) / dDeltatauB1, which coincides with the inverse of the 2nd motion adjoint.
     * (M1 plus DeltatauB1) * M2 = M1 * M2 plus J * DeltatauB1 */
    Eigen::Matrix6d jac_right_composition_wrt_first(const homogeneous& M) const {return M.adjoint_matrix_backward();}
    /**< returns the right jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(M1 * M2) / dDeltatauB2, which as a matter of fact is the identity.
     * M1 * (M2 plus DeltatauB2) = M1 * M2 plus J * DeltatauB2 */
    Eigen::Matrix6d jac_right_composition_wrt_second(const homogeneous& M) const {return Eigen::Matrix6d::Identity();}

    /**< returns the right jacobian of the forward motion action with respect to the motion, equal to d(M * p) / dDeltatauB,
     * (M plus DeltatauB) * p = M * p + J * DeltatauB */
    Eigen::Matrix36d jac_right_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the right jacobian of the backward motion action with respect to the motion, equal to d(M / p) / dDeltatauB,
     * (M plus DeltatauB) / p = M / p + J * DeltatauB */
    Eigen::Matrix36d jac_right_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the right jacobian of the Log function, equal to d(Log(M)) / dDeltatau, which coincides with the inverse
     * of the right jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(M plus Deltatau) ~= Log(M) + (J * Deltatau) */
    Eigen::Matrix6d jac_right_log() const;

    /**< returns the right jacobian of the motion right plus operator with respect to the group object (first element),
     * equal to d(M1 plus tau2) / dDeltatau1B.
     * (M1 plus Deltatau1B) plus tau2 = (M1 plus tau2) plus J * Deltatau1B */
    Eigen::Matrix6d jac_right_plus_wrt_first(const trfv& tau) const;
    /**< returns the right jacobian of the motion right plus operator with respect to the tangent space object (second element),
     * equal to d(M1 plus tau2) / dDeltatau2B.
     * M1 plus (tau2 + Deltatau2B) = (M1 plus tau2) plus J * Deltatau2B */
    Eigen::Matrix6d jac_right_plus_wrt_second(const trfv& tau) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the first element,
     * equal to d(M2 minus M1) / dDeltatau2B.
     * (M2 plus Deltatau2B) minus M1 = (M2 minus M1) + J * Deltatau2B */
    Eigen::Matrix6d jac_right_minus_wrt_first(const homogeneous& M) const;
    /**< returns the right jacobian of the motion right minus operator with respect to the second element,
     * equal to d(M2 minus M1) / dDeltatau1B.
     * M2 minus(M1 plus Deltatau1B) = (M2 minus M1) + J * Deltatau1B */
    Eigen::Matrix6d jac_right_minus_wrt_second(const homogeneous& M) const;

    /**< returns the right jacobian of the forward adjoint action with respect to the motion, equal to d(AdM | xi) / dDeltatauB,
     * Ad(M plus DeltatauB) | xi = AdM | xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the right jacobian of the backward adjoint action with respect to the motion, equal to d(AdM % xi) / dDeltatauB,
     * Ad(M plus DeltatauB) % xi = AdM % xi + J * DeltatauB */
    Eigen::Matrix6d jac_right_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
    /**< returns the left jacobian of the inverse motion, equal to d(M^-1)/dDeltatauE, which coincides with the negative of the inverse adjoint.
     * (DeltatauE plus M)^-1 = J * DeltatauE plus M^-1 */
    Eigen::Matrix6d jac_left_inverse() const {return - this->adjoint_matrix_backward();}

    /**< returns the left jacobian of the motion composition with respect to the first (left) motion, equal to
     * d(M1 * M2) / dDeltatauE1, which as a matter of fact is the identity.
     * (DeltatauE1 plus M1) * M2 = J * DeltatauE1 plus M1 * M2 */
    Eigen::Matrix6d jac_left_composition_wrt_first(const homogeneous& M) const {return Eigen::Matrix6d::Identity();}
    /**< returns the left jacobian of the motion composition with respect to the second (right) motion, equal to
     * d(M1 * M2) / dDeltatauE2, which coincides with the adjoint of the 1st motion.
     * M1 * (DeltatauE2 plus M2) = J * DeltatauE2 plus M1 * M2 */
    Eigen::Matrix6d jac_left_composition_wrt_second(const homogeneous& M) const {return this->adjoint_matrix_forward();}

    /**< returns the left jacobian of the forward motion action with respect to the motion, equal to d(M * p) / dDeltatauE,
     * (DeltatauE plus M) * p = M * p + J * DeltatauE */
    Eigen::Matrix36d jac_left_forward_motion_wrt_motion(const Eigen::Vector3d& p) const;
    /**< returns the left jacobian of the backward motion action with respect to the motion, equal to d(M / p) / dDeltatauE,
     * (DeltatauE plus M) / p = M / p + J * DeltatauE */
    Eigen::Matrix36d jac_left_backward_motion_wrt_motion(const Eigen::Vector3d& p) const;

    /**< returns the left jacobian of the Log function, equal to d(Log(M)) / dDeltatau, which coincides with the inverse
     * of the left jacobian..
     * The 1st order Taylor approximation is valid for very small transform vector changes:
     * Log(Deltatau plus M) ~= Log(M) + (J * Deltatau) */
    Eigen::Matrix6d jac_left_log() const;

    /**< returns the left jacobian of the motion left plus operator with respect to the tangent space object (first element),
     * equal to d(tau1 plus M2) / dDeltatau1E.
     * (tau1 + Deltatau1E) plus M2 = J * Deltatau1E plus (tau1 plus M2)  */
    Eigen::Matrix6d jac_left_plus_wrt_first(const trfv& tau) const;
    /**< returns the left jacobian of the motion left plus operator with respect to the group object (second element),
     * equal to d(tau1 plus M2) / dDeltatau2E.
     * tau1 plus (Deltatau2E plus M2) = J * Deltatau2E plus (tau1 plus M2) */
    Eigen::Matrix6d jac_left_plus_wrt_second(const trfv& tau) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the first element,
     * equal to d(M2 minus M1) / dDeltatau2E.
     * (Deltatau2E plus M2) minus M1 = (M2 minus M1) + J * Deltatau2E */
    Eigen::Matrix6d jac_left_minus_wrt_first(const homogeneous& M) const;
    /**< returns the left jacobian of the motion left minus operator with respect to the second element,
     * equal to d(M2 minus M1) / dDeltatau1E.
     * M2 minus (Deltatau1E plus M1) = (M2 minus M1) + J * Deltatau1E */
    Eigen::Matrix6d jac_left_minus_wrt_second(const homogeneous& M) const;

    /**< returns the left jacobian of the forward adjoint action with respect to the motion, equal to d(AdM | xi) / dDeltatauE,
     * Ad(DeltatauE plus M) | xi = AdM | xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_forward_adjoint_wrt_motion(const se3_tangent& xi) const;
    /**< returns the left jacobian of the backward adjoint action with respect to the motion, equal to d(AdM % xi) / dDeltatauE,
     * Ad(DeltatauE plus M) % xi = AdM % xi + J * DeltatauE */
    Eigen::Matrix6d jac_left_backward_adjoint_wrt_motion(const se3_tangent& xi) const;

    /**< ===== ===== Euclidean Jacobians ===== ===== */
    /**< returns the jacobian of the forward motion action with respect to the point, equal to d(M * p) / dDeltap,
     * M * (p + Deltap) = M * p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_forward_motion_wrt_point() const;
    /**< returns the jacobian of the backward motion action with respect to the point, equal to d(M / p) / dDeltap,
     * M / (p + Deltap) = M / p + J * Deltap */
    Eigen::Matrix3d jac_euclidean_backward_motion_wrt_point() const;

    /**< returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(AdM | xi) / dDeltaxi,
     * AdM | (xi + Deltaxi) = AdM | xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_forward_adjoint_wrt_tangent() const;
    /**< returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(AdM % xi) / dDeltaxi,
     * AdM % (xi + Deltaxi) = AdM % xi + J * Deltaxi */
    Eigen::Matrix6d jac_euclidean_backward_adjoint_wrt_tangent() const;

    /**< returns the jacobian [3x12] of a forward transformation with respect to the homogeneous matrix, equal to d(M * p)/dM.
     * The forward transformation (unlike the backward one) is linear on the homogeneous matrix, so the following expressions are true:
     * M * p = d(M * p)/dM |M*p * M
     * [M * exp(Delta tau)] * p ~= M * p + d(M * p)/dM |M*p * {[M * exp(Delta tau)] - M}
     * [exp(Delta tau) * M] * p ~= M * p + d(M * p)/dM |M*p * {[exp(Delta tau) * M] - M}
     * Note that the jacobian is evaluated at |M*p.
     * Note that the increment is {[M * exp(Delta tau)] - M} or [(exp(Delta tau) * M) - M], and not exp(Delta tau) or Delta M.
     * Note that I discard the last row of the homogeneous matrix (it does not add anything) to obtain size 3x12. */
    Eigen::Matrix312d jac_euclidean_forward_motion_wrt_homogeneous(const Eigen::Vector3d& p) const;
    /**< returns the jacobian [3x12] of a backward transformation with respect to the homogeneous matrix, equal to d(M / p)/dM.
     * The backward transformation (unlike the forward one) is NOT linear on the homogeneous matrix. The 1st order Taylor approximation
     * is valid for very small affine rotation matrix changes:
     * [M * exp(Delta tau)] / p ~= M / p + d(M / p)/dM |M/p * {[M * exp(Delta tau)] - M}
     * [exp(Delta tau) * M] / p ~= M / p + d(M / p)/dM |M/p * {[exp(Delta tau) * M] - M}
     * Note that the jacobian is evaluated at |M/p.
     * Note that the increment is {[M plus exp(Delta tau)] - M} or [(M plus Delta M) - M], and not exp(Delta tau) or Delta M.
     * Note that I discard the last row of the homogeneous matrix (it does not add anything) to obtain size 3x12. */
    Eigen::Matrix312d jac_euclidean_backward_motion_wrt_homogeneous(const Eigen::Vector3d& p) const;

    /**< ===== ===== Getters ===== ===== */
    /**< return translation vector */
    Eigen::Vector3d get_T() const {return this->topRightCorner(3,1);}
    /**< return dcm object */
    ang::dcm get_dcm() const {return ang::dcm(this->topLeftCorner(3,3));}
    /**< return rotation vector object */
    ang::rotv get_rotv() const;

    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const {return - this->topLeftCorner(3,3).transpose() * this->topRightCorner(3,1);}
    /**< return direction cosine matrix rotation of inverse or opposite transformation */
    ang::dcm get_inverse_dcm() const {return ang::dcm(this->topLeftCorner(3,3).transpose());}

    /**< ===== ===== Setters ===== ===== */
    /**< set the object equal to the input matrix (I SHOULD ADD AN ORTHONORMALIZATION METHOD HERE FOR THE ROTATION PART) */
    void set(const Eigen::Matrix4d& Omat) {static_cast<Eigen::Matrix4d&>(*this) = Omat;}
    /**< set the object equal to the input direction cosine matrix and translation vector */
    void set(const ang::dcm& R, const Eigen::Vector3d& tr) {
        this->topLeftCorner(3, 3) = R();
        this->topRightCorner(3, 1) = tr;
        this->bottomLeftCorner(1,3) = Eigen::RowVector3d::Zero(1,3);
        (*this)()(3,3) = 1;
    }

    /**< ===== ===== Other ===== ===== */
    /**< normalize the rotation matrix block ensuring it is orthonormal */
    void normalize();
}; // closes class homogeneous

} // closes namespace ang

#endif

