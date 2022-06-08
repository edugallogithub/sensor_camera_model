#ifndef ANG_ALGEBRA_DUAL_QUAT
#define ANG_ALGEBRA_DUAL_QUAT

#include "ang.h"
#include "quat.h"
#include "auxiliary.h" // required for the proper behavior of the << operator
#include <Eigen/Core>
#include <ostream>

/*
This file contains the "dual quat" class (models a dual quaternion). It
contains overloaded operators and the most typical algebraic operations. */

namespace ang {

// CLASS DUAL QUAT
// ===============
// ===============

class dual_quat {
protected:
    /**< real quaternion */
    ang::quat _qr;
    /**< dual quaternion */
    ang::quat _qd;

    friend class dual; // dual has access to the dual_quat private attributes
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor (does not initialize) */
    dual_quat() = default;
    /**< constructor based on two quaternions */
    dual_quat(const quat& qr, const quat& qd) : _qr(qr), _qd(qd) {}
    /**< move constructor based on two quaternions */
    dual_quat(quat&& qr, quat&& qd) : _qr(qr), _qd(qd) {}
    /**< copy constructor */
    dual_quat(const dual_quat& z) = default;
    /**< move constructor */
    dual_quat(const dual_quat&& z) noexcept: _qr(z._qr), _qd(z._qd) {}
    /**< copy assignment */
    dual_quat& operator=(const dual_quat& z) {
        this->_qr = z._qr;
        this->_qd = z._qd;
        return *this;
    }
    /**< move assignment */
    dual_quat& operator=(dual_quat&& z) {
        this->_qr = z._qr;
        this->_qd = z._qd;
        return *this;
    }

    /**< ===== ===== Overloaded Operators ===== ===== */
    /**< overloaded operator * (quaternion product) */
    inline dual_quat operator*(const dual_quat& o) const {
        dual_quat res;
        res._qr = this->_qr * o._qr;
        res._qd = this->_qr * o._qd + this->_qd * o._qr;
        return res;
    }
    /**< overloaded operator * (scalar product) */
    dual_quat operator*(const double& o) const {
        return {_qr * o, _qd * o};
    }
    /**< overloaded operator / (scalar division) */
    dual_quat operator/(const double& o) const {
        return {_qr / o, _qd / o};
    }

    /**< ===== ===== Other Methods ===== ===== */
    /**< returns circle adjoint dual quaternion */
    inline dual_quat circle_adjoint() const {
        return {this->_qr, -this->_qd};
    }
    /**< returns asterisk adjoint dual quaternion */
    inline dual_quat ast_adjoint() const {
        return {this->_qr.adjoint(), this->_qd.adjoint()};
    }
    /**< returns bullet adjoint dual quaternion */
    inline dual_quat bullet_adjoint() const {
        return {this->_qr.adjoint(), - this->_qd.adjoint()};
    }

    /**< ===== ===== Linear Algebra ===== ===== */
    /**< get real quaternion */
    ang::quat& get_qr() {return this->_qr;}
    const ang::quat& get_qr() const {return this->_qr;}
    /**< get dual quaternion */
    ang::quat& get_qd() {return this->_qd;}
    const ang::quat& get_qd() const {return this->_qd;}
    /**< get 8x1 vector (do not use operator () as then it is implicit I have to use everywhere) */
    Eigen::Vector8d get() const {Eigen::Vector8d res; res << _qr.get(), _qd.get(); return res;}

}; // closes class quat

} // closes namespace ang

#endif
