#include "rodrigues.h"
#include "euler.h"
#include "dcm.h"
#include "rotv.h"
#include "so3_tangent.h"
#include <cmath>
#include "../auxiliary.h"
#include "../tools.h"
#include "../quat.h"
#include <iostream>

using namespace ang;

// CLASS RODRIGUES
// ===============
// ===============

/* ===== ===== Constructors ===== ===== */
/* ==================================== */

rodrigues::rodrigues(const euler& Oeuler)  {
	double sy = std::sin(Oeuler.get_yaw_rad());
	double cy = std::cos(Oeuler.get_yaw_rad());
	double sp = std::sin(Oeuler.get_pitch_rad());
	double cp = std::cos(Oeuler.get_pitch_rad());
	double sr = std::sin(Oeuler.get_bank_rad());
	double cr = std::cos(Oeuler.get_bank_rad());
    double tr = + cp * cy + cr * cy + sr * sp * sy + cr * cp;
    if (tr > 0) {
        (*this)()(0) = 0.5  * std::sqrt(1 + tr);
        (*this)()(1) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(0);
        (*this)()(2) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(0);
        (*this)()(3) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(0);
    }
    else if ((+ cp * cy - cr * cy - sr * sp * sy - cr * cp) > 0) {
        (*this)()(1) = 0.5  * std::sqrt(1 + cp * cy - cr * cy - sr * sp * sy - cr * cp);
        (*this)()(0) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(1);
        (*this)()(2) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(1);
        (*this)()(3) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(1);
    }
    else if ((- cp * cy + cr * cy + sr * sp * sy - cr * cp) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - cp * cy + cr * cy + sr * sp * sy - cr * cp );
        (*this)()(0) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(2);
        (*this)()(1) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(2);
        (*this)()(3) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(2);
    }
    else if ((- cp * cy - cr * cy - sr * sp * sy + cr * cp) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - cp * cy - cr * cy - sr * sp * sy + cr * cp);
        (*this)()(0) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(3);
        (*this)()(1) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(3);
        (*this)()(2) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(3);
    }

    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
}
/* constructor based on Euler angles */

rodrigues::rodrigues(const dcm& R) {
	double tr = R().trace();
	if (tr > 0) {
        (*this)()(0) = 0.5  * std::sqrt(1 + tr);
        (*this)()(1) = 0.25 * (R()(2,1) - R()(1,2)) / (*this)()(0);
        (*this)()(2) = 0.25 * (R()(0,2) - R()(2,0)) / (*this)()(0);
        (*this)()(3) = 0.25 * (R()(1,0) - R()(0,1)) / (*this)()(0);
    }
    else if ((R()(0, 0) - R()(1, 1) - R()(2, 2)) > 0) {
        (*this)()(1) = 0.5  * std::sqrt(1 + R()(0,0) - R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(2,1) - R()(1,2)) / (*this)()(1);
        (*this)()(2) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(1);
        (*this)()(3) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(1);
    }
    else if ((- R()(0,0) + R()(1,1) - R()(2,2)) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - R()(0,0) + R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(0,2) - R()(2,0)) / (*this)()(2);
        (*this)()(1) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(2);
        (*this)()(3) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(2);
    }
    else if ((- R()(0,0) - R()(1,1) + R()(2,2)) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - R()(0,0) - R()(1,1) + R()(2,2));
        (*this)()(0) = 0.25 * (R()(1,0) - R()(0,1)) / (*this)()(3);
        (*this)()(1) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(3);
        (*this)()(2) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(3);
    }
    else {
        throw std::runtime_error("Incorrect rotation matrix.");
    }

	// ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
}
/* constructor based on direction cosine matrix */

rodrigues::rodrigues(const rotv& rotv) {
    double rv_norm = rotv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        double factor = std::sin(rv_norm/2.0) / rv_norm;
        (*this) << std::cos(rv_norm/2.0), factor * rotv();
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        double a = std::pow(rv_norm,2.);
        (*this) << (1.0 - a / 8.), (1.0 - a / 24.) * 0.5 * rotv();
    }
}
/* constructor based on rotation vector */

/* ===== ===== Move Constructors ===== ===== */
/* ========================================= */

rodrigues::rodrigues(euler&& Oeuler)  {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double tr = + cp * cy + cr * cy + sr * sp * sy + cr * cp;
    if (tr > 0) {
        (*this)()(0) = 0.5  * std::sqrt(1 + tr);
        (*this)()(1) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(0);
        (*this)()(2) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(0);
        (*this)()(3) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(0);
    }
    else if ((+ cp * cy - cr * cy - sr * sp * sy - cr * cp) > 0) {
        (*this)()(1) = 0.5  * std::sqrt(1 + cp * cy - cr * cy - sr * sp * sy - cr * cp);
        (*this)()(0) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(1);
        (*this)()(2) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(1);
        (*this)()(3) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(1);
    }
    else if ((- cp * cy + cr * cy + sr * sp * sy - cr * cp) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - cp * cy + cr * cy + sr * sp * sy - cr * cp );
        (*this)()(0) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(2);
        (*this)()(1) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(2);
        (*this)()(3) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(2);
    }
    else if ((- cp * cy - cr * cy - sr * sp * sy + cr * cp) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - cp * cy - cr * cy - sr * sp * sy + cr * cp);
        (*this)()(0) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(3);
        (*this)()(1) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(3);
        (*this)()(2) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(3);
    }

    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
}
/* move constructor based on Euler angles */

rodrigues::rodrigues(dcm&& R) {
    double tr = R().trace();
    if (tr > 0) {
        (*this)()(0) = 0.5 * std::sqrt(1 + tr);
        (*this)()(1) = (R()(2,1) - R()(1,2)) / (4 * (*this)()(0));
        (*this)()(2) = (R()(0,2) - R()(2,0)) / (4 * (*this)()(0));
        (*this)()(3) = (R()(1,0) - R()(0,1)) / (4 * (*this)()(0));
    }
    else if ((R()(0, 0) - R()(1, 1) - R()(2, 2)) > 0) {
        (*this)()(1) = 0.5 * std::sqrt(1 + R()(0,0) - R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(2,1) - R()(1,2)) / (*this)()(1);
        (*this)()(2) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(1);
        (*this)()(3) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(1);
    }
    else if ((- R()(0,0) + R()(1,1) - R()(2,2)) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - R()(0,0) + R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(0,2) - R()(2,0)) / (*this)()(2);
        (*this)()(1) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(2);
        (*this)()(3) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(2);
    }
    else if ((- R()(0,0) - R()(1,1) + R()(2,2)) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - R()(0,0) - R()(1,1) + R()(2,2));
        (*this)()(0) = 0.25 * (R()(1,0) - R()(0,1)) / (*this)()(3);
        (*this)()(1) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(3);
        (*this)()(2) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(3);
    }
    else {
        throw std::runtime_error("Incorrect rotation matrix.");
    }
    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
}
/* move constructor based on direction cosine matrix */

rodrigues::rodrigues(rotv&& rotv) {
    double rv_norm = rotv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        double factor = std::sin(rv_norm/2.0) / rv_norm;
        (*this) << std::cos(rv_norm/2.0), factor * rotv();
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        double a = std::pow(rv_norm,2.);
        (*this) << (1.0 - a / 8.), (1.0 - a / 24.) * 0.5 * rotv();
    }
}
/* move constructor based on rotation vector */

/* ===== ===== Assignments ===== ===== */
/* =================================== */

rodrigues& rodrigues::operator=(const quat& op2) {
    this->get() = op2;
    this->normalize();
    return *this;
}
/* assignment operator based on quaternion */

rodrigues& rodrigues::operator=(const euler& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double tr = + cp * cy + cr * cy + sr * sp * sy + cr * cp;
    if (tr > 0) {
        (*this)()(0) = 0.5  * std::sqrt(1 + tr);
        (*this)()(1) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(0);
        (*this)()(2) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(0);
        (*this)()(3) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(0);
    }
    else if ((+ cp * cy - cr * cy - sr * sp * sy - cr * cp) > 0) {
        (*this)()(1) = 0.5  * std::sqrt(1 + cp * cy - cr * cy - sr * sp * sy - cr * cp);
        (*this)()(0) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(1);
        (*this)()(2) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(1);
        (*this)()(3) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(1);
    }
    else if ((- cp * cy + cr * cy + sr * sp * sy - cr * cp) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - cp * cy + cr * cy + sr * sp * sy - cr * cp );
        (*this)()(0) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(2);
        (*this)()(1) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(2);
        (*this)()(3) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(2);
    }
    else if ((- cp * cy - cr * cy - sr * sp * sy + cr * cp) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - cp * cy - cr * cy - sr * sp * sy + cr * cp);
        (*this)()(0) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(3);
        (*this)()(1) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(3);
        (*this)()(2) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(3);
    }

    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
    return *this;
}
/* assignment operator based on Euler angles */

rodrigues& rodrigues::operator=(const dcm& R) {
    double tr = R().trace();
    if (tr > 0) {
        (*this)()(0) = 0.5 * std::sqrt(1 + tr);
        (*this)()(1) = (R()(2,1) - R()(1,2)) / (4 * (*this)()(0));
        (*this)()(2) = (R()(0,2) - R()(2,0)) / (4 * (*this)()(0));
        (*this)()(3) = (R()(1,0) - R()(0,1)) / (4 * (*this)()(0));
    }
    else if ((R()(0, 0) - R()(1, 1) - R()(2, 2)) > 0) {
        (*this)()(1) = 0.5 * std::sqrt(1 + R()(0,0) - R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(2,1) - R()(1,2)) / (*this)()(1);
        (*this)()(2) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(1);
        (*this)()(3) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(1);
    }
    else if ((- R()(0,0) + R()(1,1) - R()(2,2)) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - R()(0,0) + R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(0,2) - R()(2,0)) / (*this)()(2);
        (*this)()(1) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(2);
        (*this)()(3) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(2);
    }
    else if ((- R()(0,0) - R()(1,1) + R()(2,2)) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - R()(0,0) - R()(1,1) + R()(2,2));
        (*this)()(0) = 0.25 * (R()(1,0) - R()(0,1)) / (*this)()(3);
        (*this)()(1) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(3);
        (*this)()(2) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(3);
    }
    else {
        throw std::runtime_error("Incorrect rotation matrix.");
    }
    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
    return *this;
}
/* assignment operator based on direction cosine matrix */

rodrigues& rodrigues::operator=(const rotv& rotv) {
    double rv_norm = rotv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        double factor = std::sin(rv_norm/2.0) / rv_norm;
        (*this) << std::cos(rv_norm/2.0), factor * rotv();
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        double a = std::pow(rv_norm,2.);
        (*this) << (1.0 - a / 8.), (1.0 - a / 24.) * 0.5 * rotv();
    }
    return *this;
}
/* assignment operator based on rotation vector */

/* ===== ===== Move Assignments ===== ===== */
/* ======================================== */

rodrigues& rodrigues::operator=(quat&& op2) {
    this->get() = op2;
    this->normalize();
    return *this;
}
/* move assignment operator based on quaternion */

rodrigues& rodrigues::operator=(euler&& Oeuler) {
    double sy = std::sin(Oeuler.get_yaw_rad());
    double cy = std::cos(Oeuler.get_yaw_rad());
    double sp = std::sin(Oeuler.get_pitch_rad());
    double cp = std::cos(Oeuler.get_pitch_rad());
    double sr = std::sin(Oeuler.get_bank_rad());
    double cr = std::cos(Oeuler.get_bank_rad());
    double tr = + cp * cy + cr * cy + sr * sp * sy + cr * cp;
    if (tr > 0) {
        (*this)()(0) = 0.5  * std::sqrt(1 + tr);
        (*this)()(1) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(0);
        (*this)()(2) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(0);
        (*this)()(3) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(0);
    }
    else if ((+ cp * cy - cr * cy - sr * sp * sy - cr * cp) > 0) {
        (*this)()(1) = 0.5  * std::sqrt(1 + cp * cy - cr * cy - sr * sp * sy - cr * cp);
        (*this)()(0) = 0.25 * (+ sr * cp + sr * cy - cr * sp * sy) / (*this)()(1);
        (*this)()(2) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(1);
        (*this)()(3) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(1);
    }
    else if ((- cp * cy + cr * cy + sr * sp * sy - cr * cp) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - cp * cy + cr * cy + sr * sp * sy - cr * cp );
        (*this)()(0) = 0.25 * (+ sr * sy + cr * sp * cy + sp) / (*this)()(2);
        (*this)()(1) = 0.25 * (+ cp * sy - cr * sy + sr * sp * cy) / (*this)()(2);
        (*this)()(3) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(2);
    }
    else if ((- cp * cy - cr * cy - sr * sp * sy + cr * cp) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - cp * cy - cr * cy - sr * sp * sy + cr * cp);
        (*this)()(0) = 0.25 * (+ cp * sy + cr * sy - sr * sp * cy) / (*this)()(3);
        (*this)()(1) = 0.25 * (+ sr * sy + cr * sp * cy - sp) / (*this)()(3);
        (*this)()(2) = 0.25 * (+ sr * cp - sr * cy + cr * sp * sy) / (*this)()(3);
    }

    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
    return *this;
}
/* move assignment operator based on Euler angles */

rodrigues& rodrigues::operator=(dcm&& R) {
    double tr = R().trace();
    if (tr > 0) {
        (*this)()(0) = 0.5 * std::sqrt(1 + tr);
        (*this)()(1) = (R()(2,1) - R()(1,2)) / (4 * (*this)()(0));
        (*this)()(2) = (R()(0,2) - R()(2,0)) / (4 * (*this)()(0));
        (*this)()(3) = (R()(1,0) - R()(0,1)) / (4 * (*this)()(0));
    }
    else if ((R()(0, 0) - R()(1, 1) - R()(2, 2)) > 0) {
        (*this)()(1) = 0.5 * std::sqrt(1 + R()(0,0) - R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(2,1) - R()(1,2)) / (*this)()(1);
        (*this)()(2) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(1);
        (*this)()(3) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(1);
    }
    else if ((- R()(0,0) + R()(1,1) - R()(2,2)) > 0) {
        (*this)()(2) = 0.5 * std::sqrt(1 - R()(0,0) + R()(1,1) - R()(2,2));
        (*this)()(0) = 0.25 * (R()(0,2) - R()(2,0)) / (*this)()(2);
        (*this)()(1) = 0.25 * (R()(1,0) + R()(0,1)) / (*this)()(2);
        (*this)()(3) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(2);
    }
    else if ((- R()(0,0) - R()(1,1) + R()(2,2)) > 0) {
        (*this)()(3) = 0.5 * std::sqrt(1 - R()(0,0) - R()(1,1) + R()(2,2));
        (*this)()(0) = 0.25 * (R()(1,0) - R()(0,1)) / (*this)()(3);
        (*this)()(1) = 0.25 * (R()(2,0) + R()(0,2)) / (*this)()(3);
        (*this)()(2) = 0.25 * (R()(2,1) + R()(1,2)) / (*this)()(3);
    }
    else {
        throw std::runtime_error("Incorrect rotation matrix.");
    }
    // ensure first member is always positive
    if ((*this)()(0) < 0.0) {
        this->get() = - this->get();
    }
    return *this;
}
/* move assignment operator based on direction cosine matrix */

rodrigues& rodrigues::operator=(rotv&& rotv) {
    double rv_norm = rotv().norm();
    if (rv_norm > math::constant::SMALL_ROT()) {
        double factor = std::sin(rv_norm/2.0) / rv_norm;
        (*this) << std::cos(rv_norm/2.0), factor * rotv();
    }
    else { // truncated expression --> refer to Tso3:test_exp_log_small
        double a = std::pow(rv_norm,2.);
        (*this) << (1.0 - a / 8.), (1.0 - a / 24.) * 0.5 * rotv();
    }
    return *this;
}
/* move assignment operator based on rotation vector */

/* ===== ===== Transformations ===== ===== */
/* ======================================= */
rodrigues rodrigues::operator*(const rodrigues& op2) const {
	rodrigues res(this->get() * op2());
    res.normalize();
	return res;
}
/* overloaded operator * (combination of rotations) */

rodrigues rodrigues::operator/(const rodrigues& op2) const {
    rodrigues res(this->get().adjoint() * op2());
    res.normalize();
    return res;
}
/* overloaded operator / (backward combination of rotations) */

rodrigues rodrigues::pow(const double& t) const {
    return this->log_map().pow(t).exp_map_rodrigues();
}
/* executes object rotation a fraction (t < 1) or a multiple (t > 1) of times.
 * Returns exponential map of the power function applied to the object logarithmic map. */

rodrigues rodrigues::slerp(const rodrigues& q0, const rodrigues& q1, const double& t) {
    if (q0.dot(q1) < 0.) { // phi > PI, theta > PI/2
        return q0 * (q0.inverse() * q1.negative()).pow(t);
    }
    else {
        return q0 * (q0.inverse() * q1).pow(t);
    }
}
/* spherical linear interpolation, returns q0 for t=0 and q1 for t=1 */

rodrigues rodrigues::plus_right(const rotv& rv) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (rv().norm() < math::constant::PI()) {
        return (*this) * rv.exp_map_rodrigues();
    }
    else {
        double new_norm = math::constant::PI() * 2. - rv().norm();
        rotv new_rv(rv / rv().norm() * new_norm * (-1));
        return (*this) * new_rv.exp_map_rodrigues();
    }
}
/* right plus operator (input rotation located in local tangent space). */

rodrigues rodrigues::plus_left(const rotv& rv) const {
    // plus operator applies to small rotations, but valid for first manifold covering (< PI)
    if (rv().norm() < math::constant::PI()) {
        return rv.exp_map_rodrigues() * (*this);
    }
    else {
        double new_norm = math::constant::PI() * 2. - rv().norm();
        rotv new_rv(rv / rv().norm() * new_norm * (-1));
        return new_rv.exp_map_rodrigues() * (*this);
    }
}
/* left plus operator (input rotation located in global tangent space). */

/* ===== ===== Operations ===== ===== */
/* ================================== */

Eigen::Vector3d rodrigues::operator*(const Eigen::Vector3d& vecin) const {
	quat quatin = quat::convert_3dto4d(vecin);
    return quat::convert_4dto3d(this->get() * (quatin * this->adjoint()) );
}
/* overloaded operator * (forward rotation) */

Eigen::Vector3d rodrigues::operator/(const Eigen::Vector3d& vecin) const {
    quat quatin = quat::convert_3dto4d(vecin);
    return quat::convert_4dto3d( this->adjoint() * (quatin * (*this)) );
}
/* overloaded operator / (backward rotation) */

rotv rodrigues::minus_right(const rodrigues& q) const {
    return (q.inverse() * (*this)).log_map();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* right minus operator (output rotation located in local tangent space). */

rotv rodrigues::minus_left(const rodrigues& q) const {
    return ((*this) * q.inverse()).log_map();
    // the result applies to small rotations, but valid for first manifold covering (<PI). Not verified though.
}
/* left minus operator (output rotation located in global tangent space). */

rotv rodrigues::log_map() const {
    return rotv(*this);
}
/* logarithmic map that returns the rotation vector */

/* ===== ===== Adjoint ===== ===== */
/* =================================== */

so3_tangent_quat rodrigues::operator|(const so3_tangent_quat& w_quat) const {
    return so3_tangent_quat(quat::convert_3dto4d((*this) * quat::convert_4dto3d(w_quat())));
}
/* overloaded operator | (forward adjoint) */

so3_tangent rodrigues::operator|(const so3_tangent& w) const {
    return so3_tangent((*this) * w());
}
/* overloaded operator | (forward adjoint) */

Eigen::Matrix3d rodrigues::adjoint_matrix_forward() const {
    return dcm(*this)();
}
/* returns forward adjoint matrix */

so3_tangent_quat rodrigues::operator%(const so3_tangent_quat& w_quat) const {
    return so3_tangent_quat(quat::convert_3dto4d((*this) / quat::convert_4dto3d(w_quat())));
}
/* overloaded operator % (backward adjoint) */

so3_tangent rodrigues::operator%(const so3_tangent& w) const {
    return so3_tangent((*this) / w());
}
/* overloaded operator % (backward adjoint) */

Eigen::Matrix3d rodrigues::adjoint_matrix_backward() const {
    return dcm(this->inverse())();
}
/* returns backward adjoint matrix */

/* ===== ===== Angular Velocity - Time Derivative ===== ===== */
/* ========================================================== */
so3_tangent rodrigues::dot2omegabody(const quat& rodriguesdot) const {
    return so3_tangent(quat::convert_4dto3d( 2.0 * (this->adjoint() * rodriguesdot)));
}
/* obtains the body angular velocity from the Rodrigues parameters and their time derivative. */

quat rodrigues::omegabody2dot(const so3_tangent& w_body_rps) const {
	return 0.5 * quat(this->get() * quat::convert_3dto4d(w_body_rps()));
}
/* obtains the Rodrigues parameters differentials with time based on the Rodrigues parameters and the body angular velocity. */

so3_tangent rodrigues::dot2omegaspace(const quat& rodriguesdot) const {
    return so3_tangent(quat::convert_4dto3d(2.0 * (rodriguesdot * this->adjoint())));
}
/* obtains the space angular velocity from the Rodrigues parameters and their time derivative. */

quat rodrigues::omegaspace2dot(const so3_tangent& w_space_rps) const {
    return 0.5 * (quat(quat::convert_3dto4d(w_space_rps())) * this->get());
}
/* obtains the Rodrigues parameters differentials with time based on the Rodrigues parameters and the space angular velocity. */

/* ===== ===== Hat and Wedge Operators ===== ===== */
/* =============================================== */
Eigen::Vector4d rodrigues::wedge(const rodrigues& q) {
    return q();
}
/* although the wedge operator usually applies to the tangent space, here it takes the four components of the
 * unit quaternion in rodrigues form and returns them in vector form. */

rodrigues rodrigues::hat(const Eigen::Vector4d& v) {
    return rodrigues(v);
}
/* although the hat operator usually applies to the tangent space, here it takes the four components of the
 * unit quaternion in vector form, and returns them in unit quaternion form (object). */

/* ===== ===== Right (Local Tangent Space) Jacobians ===== ===== */
/* ============================================================= */
Eigen::Matrix3d rodrigues::jac_right_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - dcm(*this)() *tools::skew3(vec);
}
/* returns the right jacobian of the forward rotation action with respect to the rotation, equal to d(R * v) / dDeltarB,
 * (R plus DeltarB) * v = R * v + J * DeltarB */

Eigen::Matrix3d rodrigues::jac_right_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    //return dcm(*this).inverse()() * tools::skew3(vec) * dcm(*this)(); // this result is identical but more expensive to compute
    return tools::skew3(dcm(*this).inverse()() * vec);
}
/* returns the right jacobian of the backward rotation action with respect to the rotation, equal to d(q / v) / dDeltarB,
 * (q plus DeltarB) / v = q / v + J * DeltarB */

Eigen::Matrix3d rodrigues::jac_right_log() const {
    return this->log_map().jac_right_inv();
}
/* returns the right jacobian of the Log function, equal to d(Log(q)) / dDeltarv, which coincides with the inverse
 * of the right jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(q plus Deltarv) ~= Log(q) + (J * Deltarv) */

Eigen::Matrix3d rodrigues::jac_right_plus_wrt_first(const rotv& rv) const {
    return rv.exp_map_rodrigues().adjoint_matrix_backward();
}
/* returns the right jacobian of the rotation right plus operator with respect to the group object (first element),
 * equal to d(q1 plus rv2) / dDeltar1B.
 * (q1 plus Deltar1B) plus rv2 = (q1 plus rv2) plus J * Deltar1B */

Eigen::Matrix3d rodrigues::jac_right_plus_wrt_second(const rotv& rv) const {
    return rv.jac_right();
}
/* returns the right jacobian of the rotation right plus operator with respect to the tangent space object (second element),
 * equal to d(q1 plus rv2) / dDeltar2B.
 * q1 plus (rv2 + Deltar2B) = (q1 plus rv2) plus J * Deltar2B */

Eigen::Matrix3d rodrigues::jac_right_minus_wrt_first(const rodrigues& q) const {
    return this->minus_right(q).jac_right_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the first element,
 * equal to d(q2 minus q1) / dDeltar2B.
 * (q2 plus Deltar2B) minus q1 = (q2 minus q1) + J * Deltar2B */

Eigen::Matrix3d rodrigues::jac_right_minus_wrt_second(const rodrigues& q) const {
    return - this->minus_right(q).jac_left_inv();
}
/* returns the right jacobian of the rotation right minus operator with respect to the second element,
 * equal to d(q2 minus q1) / dDeltar1B.
 * q2 minus(q1 plus Deltar1B) = (q2 minus q1) + J * Deltar1B */

Eigen::Matrix3d rodrigues::jac_right_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - dcm(*this)() * tools::skew3(w());
}
/* returns the right jacobian of the forward adjoint action with respect to the rotation, equal to d(Adq | w) / dDeltarB,
 * Ad(q plus DeltarB) | w = Adq | w + J * DeltarB */

Eigen::Matrix3d rodrigues::jac_right_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return tools::skew3(dcm(*this).inverse()() * w());
}
/* returns the right jacobian of the backward adjoint action with respect to the rotation, equal to d(Adq % w) / dDeltarB,
 * Ad(q plus DeltarB) % w = Adq % w + J * DeltarB */

/* ===== ===== Left (Global Tangent Space) Jacobians ===== ===== */
/* ============================================================= */

Eigen::Matrix3d rodrigues::jac_left_forward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return - tools::skew3(dcm(*this)() * vec);
}
/* returns the left jacobian of the forward rotation action with respect to the rotation, equal to d(q * v) / dDeltarN,
 * (DeltarN plus q) * v = q * v + J * DeltarN */

Eigen::Matrix3d rodrigues::jac_left_backward_rotation_wrt_rotation(const Eigen::Vector3d& vec) const {
    return dcm(*this)().transpose() * tools::skew3(vec);
}
/* returns the left jacobian of the backward rotation action with respect to the rotation, equal to d(q / v) / dDeltarN,
 * (DeltarN plus q) / v = q / v + J * DeltarN */

Eigen::Matrix3d rodrigues::jac_left_log() const {
    return this->log_map().jac_left_inv();
}
/* returns the left jacobian of the Log function, equal to d(Log(q)) / dDeltarv, which coincides with the inverse
 * of the left jacobian..
 * The 1st order Taylor approximation is valid for very small rotation vector changes:
 * Log(Deltarv plus q) ~= Log(q) + (J * Deltarv) */

Eigen::Matrix3d rodrigues::jac_left_plus_wrt_first(const rotv& rv) const {
    return rv.jac_left();
}
/* returns the left jacobian of the rotation left plus operator with respect to the tangent space object (first element),
 * equal to d(rv1 plus q2) / dDeltar1N.
 * (rv1 + Deltar1N) plus q2 = J * Deltar1N plus (rv1 plus q2)  */

Eigen::Matrix3d rodrigues::jac_left_plus_wrt_second(const rotv& rv) const {
    return rv.exp_map_rodrigues().adjoint_matrix_forward();
}
/* returns the left jacobian of the rotation left plus operator with respect to the group object (second element),
 * equal to d(rv1 plus q2) / dDeltar2N.
 * rv1 plus (Deltar2N plus q2) = J * Deltar2N plus (rv1 plus q2) */

Eigen::Matrix3d rodrigues::jac_left_minus_wrt_first(const rodrigues& q) const {
    return this->minus_left(q).jac_left_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the first element,
 * equal to d(q2 minus q1) / dDeltar2N.
 * (Deltar2N plus q2) minus q1 = (q2 minus q1) + J * Deltar2N */

Eigen::Matrix3d rodrigues::jac_left_minus_wrt_second(const rodrigues& q) const {
    return - this->minus_left(q).jac_right_inv();
}
/* returns the left jacobian of the rotation left minus operator with respect to the second element,
 * equal to d(q2 minus q1) / dDeltar1N.
 * q2 minus (Deltar1N plus q1) = (q2 minus q1) + J * Deltar1N */

Eigen::Matrix3d rodrigues::jac_left_forward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return - tools::skew3(dcm(*this)() * w());
}
/* returns the left jacobian of the forward adjoint action with respect to the rotation, equal to d(Adq | w) / dDeltarN,
 * Ad(DeltarN plus q) | w = Adq | w + J * DeltarN */

Eigen::Matrix3d rodrigues::jac_left_backward_adjoint_wrt_rotation(const so3_tangent& w) const {
    return dcm(*this)().transpose() * tools::skew3(w());
}
/* returns the left jacobian of the backward adjoint action with respect to the rotation, equal to d(Adq % w) / dDeltarN,
 * Ad(DeltarN plus q) % w = Adq % w + J * DeltarN */

/* ===== ===== Euclidean Jacobians ===== ===== */
/* =========================================== */

Eigen::Matrix3d rodrigues::jac_euclidean_forward_rotation_wrt_vector() const {
    return dcm(*this)();
}
/* returns the jacobian of the forward rotation action with respect to the vector, equal to d(q * v) / dDeltav,
 * q * v = J * v
 * q * (v + Deltav) = q * v + J * Deltav */

Eigen::Matrix3d rodrigues::jac_euclidean_backward_rotation_wrt_vector() const {
    return dcm(this->inverse())();
}
/* returns the jacobian of the backward rotation action with respect to the vector, equal to d(q / v) / dDeltav,
 * q / v = J * v
 * q / (v + Deltav) = q / v + J * Deltav */

Eigen::Matrix3d rodrigues::jac_euclidean_forward_adjoint_wrt_tangent() const {
    return dcm(*this)();
}
/* returns the jacobian of the forward adjoint action with respect to the tangent space component, equal to d(Adq | w) / dDeltaw,
 * Adq | w = J * w
 * Adq | (w + Deltaw) = Adq | w + J * Deltaw */

Eigen::Matrix3d rodrigues::jac_euclidean_backward_adjoint_wrt_tangent() const {
    return this->adjoint_matrix_backward();
}
/* returns the jacobian of the backward adjoint action with respect to the tangent space component, equal to d(Adq % w) / dDeltaw,
 * Adq % w = J * w
 * Adq % (w + Deltaw) = Adq % w + J * Deltaw */

Eigen::Matrix34d rodrigues::jac_euclidean_forward_rotation_wrt_rodrigues(const Eigen::Vector3d& v) const {
    Eigen::Matrix34d res;
    res << + (*this)()(0) * v(0) - (*this)()(3) * v(1) + (*this)()(2) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2),
           - (*this)()(2) * v(0) + (*this)()(1) * v(1) + (*this)()(0) * v(2),
           - (*this)()(3) * v(0) - (*this)()(0) * v(1) + (*this)()(1) * v(2),

           + (*this)()(3) * v(0) + (*this)()(0) * v(1) - (*this)()(1) * v(2),
           + (*this)()(2) * v(0) - (*this)()(1) * v(1) - (*this)()(0) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2),
           + (*this)()(0) * v(0) - (*this)()(3) * v(1) + (*this)()(2) * v(2),

           - (*this)()(2) * v(0) + (*this)()(1) * v(1) + (*this)()(0) * v(2),
           + (*this)()(3) * v(0) + (*this)()(0) * v(1) - (*this)()(1) * v(2),
           - (*this)()(0) * v(0) + (*this)()(3) * v(1) - (*this)()(2) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2);
    return res * 2.0;
}
/* returns the jacobian [3x4] of a forward rotation with respect to the quaternion, equal to d(q * v)/dq.
 * The forward rotation is BI linear on the quaternion, so
 * q * v = 0.5 * d(q * v)/dq |q*v * q
 * The 1st order Taylor approximation is valid for very small quaternion changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [q * exp(Delta r)] * v ~= q * v + d(q * v)/dq |q*v * {[q * exp(Delta r)] - q}
 * [exp(Delta r) * q] * v ~= q * v + d(q * v)/dq |q*v * {[exp(Delta r) * q] - q}
 * Note that the jacobian is evaluated at |q*v.
 * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
 * to each other, and never exp(Delta r) or Delta q. */

Eigen::Matrix34d rodrigues::jac_euclidean_backward_rotation_wrt_rodrigues(const Eigen::Vector3d& v) const {
    Eigen::Matrix34d res;
    res << + (*this)()(0) * v(0) + (*this)()(3) * v(1) - (*this)()(2) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2),
           - (*this)()(2) * v(0) + (*this)()(1) * v(1) - (*this)()(0) * v(2),
           - (*this)()(3) * v(0) + (*this)()(0) * v(1) + (*this)()(1) * v(2),

           - (*this)()(3) * v(0) + (*this)()(0) * v(1) + (*this)()(1) * v(2),
           + (*this)()(2) * v(0) - (*this)()(1) * v(1) + (*this)()(0) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2),
           - (*this)()(0) * v(0) - (*this)()(3) * v(1) + (*this)()(2) * v(2),

           + (*this)()(2) * v(0) - (*this)()(1) * v(1) + (*this)()(0) * v(2),
           + (*this)()(3) * v(0) - (*this)()(0) * v(1) - (*this)()(1) * v(2),
           + (*this)()(0) * v(0) + (*this)()(3) * v(1) - (*this)()(2) * v(2),
           + (*this)()(1) * v(0) + (*this)()(2) * v(1) + (*this)()(3) * v(2);
    return res * 2.0;
}
/* returns the jacobian [3x4] of a backward rotation with respect to the quaternion, equal to d(q / v)/dq.
 * The backward rotation is BI linear on the quaternion, so
 * q / v = 0.5 * d(q / v)/dq |q/v * q
 * The 1st order Taylor approximation is valid for very small quaternion changes, vut it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [q * exp(Delta r)] / v ~= q / v + d(q / v)/dq |q/v * {[q * exp(Delta r)] - q}
 * [exp(Delta r) * q] / v ~= q / v + d(q / v)/dq |q/v * {[exp(Delta r) * q] - q}
 * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
 * to each other, and never exp(Delta r) or Delta q. */

Eigen::Matrix34d rodrigues::jac_euclidean_forward_adjoint_wrt_rodrigues(const so3_tangent& w) const {
    return this->jac_euclidean_forward_rotation_wrt_rodrigues(w());
}
/* returns the jacobian [3x4] of a forward adjoint with respect to the quaternion, equal to d(q | w)/dq.
 * The forward adjoint is BI linear on the quaternion, so
 * q | w = 0.5 * d(q | w)/dq |q|w * q
 * The 1st order Taylor approximation is valid for very small quaternion changes, but it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [q * exp(Delta r)] | w ~= q | w + d(q | w)/dq |q|w * {[q * exp(Delta r)] - q}
 * [exp(Delta r) * q] | w ~= q | w + d(q | w)/dq |q|w * {[exp(Delta r) * q] - q}
 * Note that the jacobian is evaluated at |q|w.
 * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
 * to each other, and never exp(Delta r) or Delta q. */

Eigen::Matrix34d rodrigues::jac_euclidean_backward_adjoint_wrt_rodrigues(const so3_tangent& w) const {
    return this->jac_euclidean_backward_rotation_wrt_rodrigues(w());
}
/* returns the jacobian [3x4] of a backward adjoint with respect to the quaternion, equal to d(q % w)/dq.
 * The backward adjoint is BI linear on the quaternion, so
 * q % w = 0.5 * d(q % w)/dq |q%w * q
 * The 1st order Taylor approximation is valid for very small quaternion changes, vut it does not matter if
 * the perturbation is local or global, as the two following expressions are both valid:
 * [q * exp(Delta r)] % w ~= q % w + d(q % w)/dq |q%w * {[q * exp(Delta r)] - q}
 * [exp(Delta r) * q] % w ~= q % w + d(q % w)/dq |q%w * {[exp(Delta r) * q] - q}
 * Note that the increment is either {[q * exp(Delta r)] - q} or {[exp(Delta r) * q] - q}, which are not equal
 * to each other, and never exp(Delta r) or Delta q. */








































