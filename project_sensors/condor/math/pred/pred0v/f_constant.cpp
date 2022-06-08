#include "f_constant.h"

// CLASS F_CONSTANT
// ================
// ================

const math::logic::PRED_NAME math::f_constant::_name = math::logic::f_constant;
/* predicate name */

const std::string math::f_constant::_st_name = "f_constant";
/* predicate name string */

math::f_constant::f_constant(const double& f0)
: _f0(f0) {}
/* constructor based on constant value for output in standard units */

math::f_constant::f_constant(const f_constant& other)
: _f0(other._f0) {}
/* copy constructor */

math::f_constant* math::f_constant::clone() const {
	return new f_constant(*this);
}
/* cloner */

bool math::f_constant::operator==(const pred0v& op2) const {
	return (op2.get_name() == math::logic::f_constant) ?
		(*this == static_cast<const f_constant&>(op2)) : false;	
}
bool math::f_constant::operator==(const f_constant& op2) const {
	return (fabs(_f0 - op2._f0) < math::constant::TOL());
}
/* overloaded operator == (equal) */





