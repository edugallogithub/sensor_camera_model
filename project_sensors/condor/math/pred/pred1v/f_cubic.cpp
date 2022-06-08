#include "f_cubic.h"

// CLASS F_CUBIC
// =============
// =============

const math::logic::PRED_NAME math::f_cubic::_name = math::logic::f_cubic;
/* predicate name */

const std::string math::f_cubic::_st_name = "f_cubic";
/* predicate name string */

math::f_cubic::f_cubic(const double f0,
					  const double f1,
					  const double f2,
					  const double f3)
: _f0(f0), _f1(f1), _f2(f2), _f3(f3) {}
/* constructor based on cubic function coefficients for input
and output in standard units */

math::f_cubic::f_cubic(const f_cubic& other)
: _f0(other._f0), _f1(other._f1), _f2(other._f2), _f3(other._f3) {}
/* copy constructor */

bool math::f_cubic::operator==(const pred1v& op2) const {
	return (op2.get_name() == math::logic::f_cubic) ?
		(*this == static_cast<const f_cubic&>(op2)) : false;
}
bool math::f_cubic::operator==(const f_cubic& op2) const {
	return ((fabs(_f0 - op2._f0) < constant::TOL()) &&
		    (fabs(_f1 - op2._f1) < constant::TOL()) &&
		    (fabs(_f2 - op2._f2) < constant::TOL()) &&
			(fabs(_f3 - op2._f3) < constant::TOL()));
}
/* overloaded operator == (equal) */

math::f_cubic* math::f_cubic::clone() const {
	return new f_cubic(*this);
}
/* cloner */

double math::f_cubic::value(const double& input) const {
    return _f0 + input * _f1 + pow(input,2) * _f2 + pow(input,3) * _f3;
}
/* evaluates the function at the reference magnitude input, and writes the
result at the reference magnitude result. Only the magnitude value is
inserted into result, the units are assummed to be OK. */

double math::f_cubic::d_dt(const double& input, const double& input_dt) const {
	return input_dt * _f1 + 2 * input * input_dt * _f2 + 3 * pow(input,2) * input_dt * _f3;
};
/* evaluates the function differential with time at the reference
magnitude input and its differential with time input_dt, and writes the
result at the reference magnitude result. */

