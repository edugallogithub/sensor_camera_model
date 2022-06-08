#ifndef ANG_REPRESENTATION_SCREW
#define ANG_REPRESENTATION_SCREW

#include "../ang.h"
#include "../auxiliary.h"
#include "../rotate/rotv.h"

/*
 * This file contains the "screw" representation of a transformation.
 */

namespace ang {
    class speu_rodrigues;
    class speu_dcm;
    class homogeneous;
    class trfv;
    class dual;

// CLASS SCREW
// ===========
// ===========

class ANG_API screw {
private:
    /**< screw axis lines */
    Eigen::Vector3d _n;
    /**< screw axis moment */
    Eigen::Vector3d _m;
    /**< screw magnitude */
    double _phi;
    /**< screw displacement */
    double _d;
    /**< screw pitch */
    double _h;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< default constructor (not initialized) */
    screw() = default;
    /**< constructor based on axis line, moment, magnitude, displacement, and pitch */
    screw(const Eigen::Vector3d& n, const Eigen::Vector3d& m, const double& phi, const double& d, const double& h);
    /**< constructor based on rotation vector and translation vector */
    screw(const ang::rotv&, const Eigen::Vector3d&);
    /**< constructor based on special Euclidean (rodrigues) */
    explicit screw(const speu_rodrigues&);
    /**< constructor based on special Euclidean (dcm) */
    explicit screw(const speu_dcm&);
    /**< constructor based on homogeneous */
    explicit screw(const ang::homogeneous&);
    /**< constructor based on transform vector */
    explicit screw(const ang::trfv&);
    /**< constructor based on unit dual quaternion */
    explicit screw(const ang::dual&);
    /**< copy constructor */
    screw(const screw&) = default;
    /**< destructor */
    ~screw() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    screw(screw&&) = default;
    /**< move constructor based on special Euclidean (rodrigues) */
    explicit screw(ang::speu_rodrigues&&);
    /**< move constructor based on special Euclidean (dcm) */
    explicit screw(ang::speu_dcm&&);
    /**< move constructor based on homogeneous */
    explicit screw(ang::homogeneous&&);
    /**< move constructor based on transform vector */
    explicit screw(ang::trfv&&);
    /**< move constructor based on unit dual quaternion */
    explicit screw(ang::dual&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
    screw& operator=(const screw&) = default;
    /**< assignment operator = based on special Euclidean (rodrigues) */
    screw& operator=(const ang::speu_rodrigues&);
    /**< assignment operator = based on special Euclidean (dcm) */
    screw& operator=(const ang::speu_dcm&);
    /**< assignment operator = based on homogeneous */
    screw& operator=(const ang::homogeneous&);
    /**< assignment operator = based on transform vector */
    screw& operator=(const ang::trfv&);
    /**< assignment operators = based on unit dual quaternion */
    screw& operator=(const ang::dual&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    screw& operator=(screw&&) = default;
    /**< move assignment operator = based on special Euclidean (rodrigues) */
    screw& operator=(ang::speu_rodrigues&&);
    /**< move assignment operator = based on special Euclidean (dcm) */
    screw& operator=(ang::speu_dcm&&);
    /**< move assignment operator = based on homogeneous */
    screw& operator=(ang::homogeneous&&);
    /**< move assignment operator = based on transform vector */
    screw& operator=(ang::trfv&&);
    /**< move assignment operators = based on unit dual quaternion */
    screw& operator=(ang::dual&&);

    /**< ===== ===== Transformations ===== ===== */
    /**< overloaded operator * (combination of transformations) */
    screw operator*(const screw& op2) const {return {this->get_rotv() * op2.get_rotv(), this->get_rotv() * op2.get_T() + this->get_T()};}
    /**< overloaded operator / (backward combination of transformations) */
    screw operator/(const screw& op2) const {return {this->get_rotv() / op2.get_rotv(), this->get_rotv() / (op2.get_T() - this->get_T())};}
    /**< returns inverse or opposite transformation */
    screw inverse() const {return screw(_n, _m, - _phi, - _d, _h);}
    /**< executes object transformation a fraction (t < 1) or a multiple (t > 1) of times.
     * Returns input fraction (interpolation or extrapolation) of the screw. */
    screw pow(const double& t) const;
    /**< screw linear interpolation, returns S0 for t=0 and S1 for t=1 */
    static screw sclerp(const screw& S0, const screw& S1, const double& t);
    /**< right plus operator (input rotation located in local tangent space) */
    screw plus_right(const trfv&) const;
    /**< right plus operator (input rotation located in local tangent space) */
    screw plus_right(const screw& S) const {return (*this) * S;}
    /**< left plus operator (input rotation located in global tangent space) */
    screw plus_left(const trfv&) const;
    /**< left plus operator (input rotation located in global tangent space) */
    screw plus_left(const screw& S) const {return S * (*this);}

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward transformation of point, not vector) */
    Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward transformation of point, not vector) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;
    /**< overloaded operator ^ (forward transformation of vector, not point) WATCH OUT -> ^ DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator^(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward transformation of vector, not point) WATCH OUT -> & DOES NOT HAVE PRECEDENCE OVER + or -. */
    Eigen::Vector3d operator&(const Eigen::Vector3d&) const;
    /**< minus operator (output transformation located in local tangent space) */
    trfv minus_right_trfv(const screw&) const;
    /**< minus operator (output transformation located in local tangent space) */
    screw minus_right_screw(const screw& S) const {return S.inverse() * (*this);}
    /**< left minus operator (output transformation located in global tangent space) */
    trfv minus_left_trfv(const screw&) const;
    /**< left minus operator (output transformation located in global tangent space) */
    screw minus_left_screw(const screw& S) const {return (*this) * S.inverse();}
    /**< exponential map that returns special euclidean direction cosine matrix */
    ang::speu_dcm exp_map_speu_dcm() const;
    /**< exponential map that returns special euclidean rodrigues parameters */
    ang::speu_rodrigues exp_map_speu_rodrigues() const;
    /**< exponential map that returns homogeneous matrix */
    ang::homogeneous exp_map_homogeneous() const;
    /**< exponential and logarithmic map (they are the same in this case) that returns the transform vector */
    ang::trfv explog_map_trfv() const;
    /**< exponential map that returns unit dual quaternion */
    ang::dual exp_map_dual() const;

    /**< ===== ===== Adjoint ===== ===== */
    /**< returns forward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_forward() const;
    /**< returns backward adjoint matrix */
    Eigen::Matrix6d adjoint_matrix_backward() const;

    /**< ===== ===== Twist or Motion Velocity - Time Derivative ===== ===== */

    /**< ===== ===== Linear Algebra ===== ===== */

    /**< ===== ===== Hat and Wedge Operators ===== ===== */

    /**< ===== ===== Jacobians ===== ===== */

    /**< ===== ===== Getters ===== ===== */
    /**< return translation vector */
    Eigen::Vector3d get_T() const;
    /**< return rotation vector object */
    ang::rotv get_rotv() const;
    /**< return translation vector of inverse or opposite transformation */
    Eigen::Vector3d get_inverse_T() const;
    /**< return rotation vector of inverse or opposite transformation */
    ang::rotv get_inverse_rotv() const;

    /**< get point belonging to axis closest to origin */
    Eigen::Vector3d get_p() const;
    /**< get screw axis lines */
    const Eigen::Vector3d& get_n() const {return _n;}
    /**< get screw axis moment */
    const Eigen::Vector3d& get_m() const {return _m;}
    /**< get screw pitch */
    const double& get_h() const {return _h;}
    /**< get screw magnitude */
    const double& get_phi() const {return _phi;}
    /**< get screw displacement */
    const double& get_d() const {return _d;}
}; // closes class screw

/**< adds the input transform vector object to the stream (can not rely on Vector6d << operator as auxiliary.h file not included) */

} // closes namespace ang

#endif

































