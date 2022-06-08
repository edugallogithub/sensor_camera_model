#ifndef ANG_REPRESENTATION_EULER
#define ANG_REPRESENTATION_EULER

#include "../ang.h"
#include "../auxiliary.h"
#include "math/logic/constant.h"
#include <Eigen/Core>
#include <ostream>
#include <iomanip>

/*
 * This file contains the "euler" class that models the Euler angles representation of a rotation.
 */

namespace ang {

class rodrigues;
class dcm;
class rotv;
class so3_tangent;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// CLASS EULER
// ===========
// ===========

class ANG_API euler {
private:
    /**< yaw, pitch, and roll */
    double _ypr_rad[3];
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< ===== ===== Constructors ===== ===== */
    /**< empty constructor (not initialized) */
	euler() = default;
	/**< constructor based on yaw, pitch, and bank angles */
	euler(double yaw_rad, double pitch_rad, double bank_rad) : _ypr_rad{yaw_rad, pitch_rad, bank_rad} {}
	/**< constructor based on Rodrigues parameters */
	explicit euler(const rodrigues&);
	/**< constructor based on direction cosine matrix */
	explicit euler(const dcm&);
	/**< constructor based on rotation vector */
	explicit euler(const rotv&);
    /**< copy constructor */
    euler(const euler&) = default;
    /**< destructor */
    virtual ~euler() = default;

    /**< ===== ===== Move Constructors ===== ===== */
    /**< move constructor */
    euler(euler&&) = default;
    /**< move constructor based on Rodrigues parameters */
    explicit euler(rodrigues&&);
    /**< move constructor based on direction cosine matrix */
    explicit euler(dcm&&);
    /**< move constructor based on rotation vector */
    explicit euler(rotv&&);

    /**< ===== ===== Assignments ===== ===== */
    /**< copy assignment */
	euler& operator=(const euler&) = default;
    /**< assignment operator based on Rodrigues parameters */
    euler& operator=(const rodrigues&);
    /**< assignment operator based on direction cosine matrix */
    euler& operator=(const dcm&);
    /**< assignment operator based on rotation vector */
    euler& operator=(const rotv&);

    /**< ===== ===== Move Assignments ===== ===== */
    /**< move assignment */
    euler& operator=(euler&&) = default;
    /**< move assignment operator based on Rodrigues parameters */
    euler& operator=(rodrigues&&);
    /**< move assignment operator based on direction cosine matrix */
    euler& operator=(dcm&&);
    /**< move assignment operator based on rotation vector */
    euler& operator=(rotv&&);

    /**< ===== ===== Transformations ===== ===== */

    /**< ===== ===== Operations ===== ===== */
    /**< overloaded operator * (forward rotation) */
    Eigen::Vector3d operator*(const Eigen::Vector3d&) const;
    /**< overloaded operator / (backward rotation) */
    Eigen::Vector3d operator/(const Eigen::Vector3d&) const;

    /**< ===== ===== Adjoint ===== ===== */

    /**< ===== ===== Angular Velocity - Time Derivative ===== ===== */
    /**< obtains the body angular velocity from the Euler angles and their time derivatives. */
    ang::so3_tangent dot2omegabody(const Eigen::Vector3d& eulerdot_rps) const;
    /**< obtains the Euler angles differentials with time based on the Euler angles and the body angular velocity. */
    Eigen::Vector3d omegabody2dot(const ang::so3_tangent& w_body_rps) const;

    /**< obtains the space angular velocity from the Euler angles and their time derivatives. */
    ang::so3_tangent dot2omegaspace(const Eigen::Vector3d& eulerdot_rps) const;
    /**< obtains the Euler angles differentials with time based on the Euler angles and the space angular velocity. */
    Eigen::Vector3d omegaspace2dot(const ang::so3_tangent& w_space_rps) const;

    /**< ===== ===== Linear Algebra ===== ===== */

    /**< ===== ===== Hat and Wedge Operators ===== ===== */

    /**< ===== ===== Obtain (and Set) Individual Euler angles ===== ===== */
    /**< access to yaw angle to read or write */
    inline const double& get_yaw_rad() const {return *_ypr_rad;}
    inline double& get_yaw_rad() {return *_ypr_rad;}
    inline void set_yaw_rad(double yaw_rad) {*_ypr_rad = yaw_rad;}
    /**< access to pitch angle to read or write */
    inline const double& get_pitch_rad() const {return *(_ypr_rad + 1);}
    inline double& get_pitch_rad() {return *(_ypr_rad + 1);}
    inline void set_pitch_rad(double pitch_rad) {*(_ypr_rad + 1) = pitch_rad;}
    /**< access to bank angle to read or write*/
    inline const double& get_bank_rad() const {return *(_ypr_rad + 2);}
    inline double& get_bank_rad() {return *(_ypr_rad + 2);}
    inline void set_bank_rad(double bank_rad) {*(_ypr_rad + 2) = bank_rad;}
    /**< set all three angles */
    inline void set_rad(double yaw_rad, double pitch_rad, double bank_rad)
    {*_ypr_rad = yaw_rad; *(_ypr_rad + 1) = pitch_rad; *(_ypr_rad + 2) = bank_rad;}

    /**< ===== ===== Rotation Jacobians ===== ===== */

    /**< ===== ===== Inverse Rotation Jacobians ===== ===== */





    /**< ===== ===== Adjoint Jacobians ===== ===== */

    /**< ===== ===== Other ===== ===== */
	/**< obtains the Euler angles that define the rotation between NED and GRD
	(chi, gamma, mu) from the Euler angles between NED and BFS (psi,
	# theta, xi) and the absolute speed in NED. */
	static void obtain_euler_nedgrd(euler& euler_ng, const euler& euler_nb, const Eigen::Vector3d& v_n_mps);
    /**< returns difference between two bearings in degrees, nan if higher than plus minus 90. */
    static double control_bearing_diff(const double& a_deg, const double& b_deg);
    /**< fills up the cross track error (positive to the right, negative to the left) and the long
     * track errors (positive to the front, negative to the back) based on the North error (positive
     * to the North, negative to the South), the East error (positive to the East, negative to the
     * West), and the trajectory bearing angle */
    static void obtain_cross_long_track_errors(double& error_cross_m, double& error_long_m, const double& error_north_m, const double& error_east_m, const double& chi_rad);

    /**< ===== ===== Angles for Forward Looking Vectors (choose yaw and pitch) ===== ===== */
    /**< given a forward looking vector (generally speed) expressed in NED, it returns the yaw angle [rad] of that vector with respect to the NED frame */
    static double obtain_yaw_forward(const Eigen::Vector3d& v_ned);
    /**< given a forward looking vector (generally speed) expressed in NED, it returns the pitch angle [rad] of that vector with respect to the NED frame */
    static double obtain_pitch_forward(const Eigen::Vector3d& v_ned);

    /**< ===== ===== Angles for Downward Looking Vectors (choose pitch and roll) ===== ===== */
    /**< given a downward looking vector (generally gravitation) expressed in NED, it returns the pitch angle [rad] of that vector with respect to the NED frame */
    static double obtain_pitch_downward(const Eigen::Vector3d& g_ned);
    /**< given a downward looking vector (generally gravitation) expressed in NED, it returns the bank angle [rad] of that vector with respect to the NED frame */
    static double obtain_bank_downward(const Eigen::Vector3d& g_ned);

    /**< ===== ===== Individual Rotations ===== ===== */
    /**< ============================================ */
    /**< return rotation matrix for 1st rotation (yaw around 3rd axis) */
    ang::dcm R3_yaw() const;
    /**< return rotation matrix for 2nd rotation (pitch around 2nd axis) */
    ang::dcm R2_pitch() const;
    /**< return rotation matrix for 3rd rotation (roll around 1st axis) */
    ang::dcm R1_roll() const;

}; // closes class euler

/**< adds the input euler object to the stream */
ANG_API inline std::ostream& operator <<(std::ostream & out_str, const euler& e) {
    out_str << std::fixed << std::showpos
            << std::setprecision(10) << std::setw(16) << e.get_yaw_rad()   * math::constant::R2D()
            << std::setprecision(10) << std::setw(16) << e.get_pitch_rad() * math::constant::R2D()
            << std::setprecision(10) << std::setw(16) << e.get_bank_rad()  * math::constant::R2D();
	return out_str;
}

/**< takes the stream data into the euler angles */
ANG_API inline std::istream& operator >>(std::istream & in_str, euler& O) {
in_str >> O.get_yaw_rad() >> O.get_pitch_rad() >> O.get_bank_rad();
return in_str;
}

} // closes namespace ang

#endif

