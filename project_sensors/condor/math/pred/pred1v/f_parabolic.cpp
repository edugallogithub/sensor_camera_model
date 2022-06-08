#include "f_parabolic.h"

// CLASS F_PARABOLIC
// =================
// =================

const math::logic::PRED_NAME math::f_parabolic::_name = math::logic::f_parabolic;
/* predicate name */

const std::string math::f_parabolic::_st_name = "f_parabolic";
/* predicate name string */

math::f_parabolic::f_parabolic(const double f0,
							  const double f1,
							  const double f2)
: _f0(f0), _f1(f1), _f2(f2) {}
/* constructor based on parabolic function coefficients for input
and output in standard units */

math::f_parabolic::f_parabolic(const f_parabolic& other)
: _f0(other._f0), _f1(other._f1), _f2(other._f2) {}
/* copy constructor */

bool math::f_parabolic::operator==(const pred1v& op2) const {
	return (op2.get_name() == math::logic::f_parabolic) ?
		(*this == static_cast<const f_parabolic&>(op2)) : false;	
}
bool math::f_parabolic::operator==(const f_parabolic& op2) const {
	return ((fabs(_f0 - op2._f0) < constant::TOL()) &&
		    (fabs(_f1 - op2._f1) < constant::TOL()) &&
			(fabs(_f2 - op2._f2) < constant::TOL()));
}
/* overloaded operator == (equal) */

math::f_parabolic* math::f_parabolic::clone() const {
	return new f_parabolic(*this);
}
/* cloner */

double math::f_parabolic::value(const double& input) const {
	return _f0 + input * _f1 + pow(input,2) * _f2;
}
/* evaluates the function at the reference magnitude input, and writes the
result at the reference magnitude result. Only the magnitude value is
inserted into result, the units are assummed to be OK. */

double math::f_parabolic::d_dt(const double& input, const double& input_dt) const {
	return input_dt * _f1 + 2 * input * input_dt * _f2;
};
/* evaluates the function differential with time at the reference
magnitude input and its differential with time input_dt, and writes the
result at the reference magnitude result. */



