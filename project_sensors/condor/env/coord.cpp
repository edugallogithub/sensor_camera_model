#include "coord.h"
#include "ang/tools.h"
#include "ang/rotate/euler.h"
#include <iostream>

// CLASS GEOCENTRIC_COORD
// ======================
// ======================

void env::geocentric_coord::set(const double& theta_rad, const double& lambda_rad, const double& r_m) {
    _data(0) = theta_rad;
    _data(1) = lambda_rad;
    _data(2) = r_m;
}
/* modify all attributes simultaneously */

void env::geocentric_coord::set(const Eigen::Array3d& data) {
    _data = data;
}
/* modify all attributes simultaneously */

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CLASS GEODETIC_COORD
// ====================
// ====================

env::geodetic_coord& env::geodetic_coord::operator=(const Eigen::Array3d&& data) {
    _data = data;
    return *this;
}
/* copy assignment based on an array containing longitude, latitude, and geometric altitude (in this order) */

env::geodetic_coord& env::geodetic_coord::operator=(Eigen::Array3d&& data) {
    _data = data;
    return *this;
}
/* move assignment based on an array containing longitude, latitude, and geometric altitude (in this order) */

env::geodetic_coord env::geodetic_coord::operator+(const geodetic_coord& op2) const {
	return {_data(0) + op2._data(0), _data(1) + op2._data(1), _data(2) + op2._data(2)};
}
/* overloaded operator + (addition) */

env::geodetic_coord env::geodetic_coord::operator+(const Eigen::Array3d& op2) const {
    return {_data(0) + op2(0), _data(1) + op2(1), _data(2) + op2(2)};
}
/* overloaded operator + (addition) */

env::geodetic_coord env::geodetic_coord::operator-(const geodetic_coord& op2) const {
    return {ang::tools::angle_diff_rad(_data(0), op2._data(0)), ang::tools::angle_diff_rad(_data(1), op2._data(1)), _data(2) - op2._data(2)};
}
/* overloaded operator - (substraction) */

env::geodetic_coord env::geodetic_coord::operator-(const Eigen::Array3d& op2) const {
    return {ang::tools::angle_diff_rad(_data(0), op2(0)), ang::tools::angle_diff_rad(_data(1), op2(1)), _data(2) - op2(2)};
}
/* overloaded operator - (substraction) */

env::geodetic_coord env::geodetic_coord::operator*(const double& op2) const {
	return {op2 * _data(0), op2 * _data(1), op2 * _data(2)};
}
/* overloaded operator * (scalar product) */

ang::euler env::geodetic_coord::obtain_euler_en() const {
    return {_data(0), - _data(1) - math::constant::PIHALF(), 0.};
}
/* obtain Euler angles of rotation from ECEF to NED */

env::geodetic_coord env::geodetic_coord::zero_altitude() const {
    return {_data(0), _data(1), 0.};
}
/* return same geodetic coordinates but with zero altitude */

void env::geodetic_coord::set(const double& lambda_rad, const double& phi_rad, const double& h_m) {
    _data(0) = lambda_rad;
    _data(1) = phi_rad;
    _data(2) = h_m;
}
/* modify all attributes simultaneously */

void env::geodetic_coord::set(const Eigen::Array3d& data) {
    _data = data;
}
/* modify all attributes simultaneously */

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CLASS CARTESIAN_COORD
// =====================
// =====================

void env::cartesian_coord::set(const double& x1_m, const double& x2_m, const double& x3_m) {
    _data(0) = x1_m;
    _data(1) = x2_m;
    _data(2) = x3_m;
}
/* modify all attributes simultaneously */

void env::cartesian_coord::set(const Eigen::Array3d& data) {
    _data = data;
}
/* modify all attributes simultaneously */

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


