#include "f_lineal.h"

// CLASS F_LINEAL
// ==============
// ==============

const math::logic::PRED_NAME math::f_lineal::_name = math::logic::f_lineal;
/* predicate name */

const std::string math::f_lineal::_st_name = "f_lineal";
/* predicate name string */

math::f_lineal::f_lineal(const double f0,
						const double f1)
: _f0(f0), _f1(f1) {
	fill_identity();
}
/* constructor based on lineal function coefficients for input and
output in standard units */

math::f_lineal::f_lineal(const double x1,
						const double x2,
						const double y1,
						const double y2)
: _f1((y2 - y1) / (x2 - x1)) {
	_f0 = y1 - _f1 * x1;
	fill_identity();
}
/* constructor based on input variable values x1 and x2 and output
variable values y1 and y2 in standard units */

math::f_lineal::f_lineal(const math::f_lineal& other)
: _f0(other._f0), _f1(other._f1), _identity(other._identity) {}
/* copy constructor */

bool math::f_lineal::operator==(const pred1v& op2) const {
	return (op2.get_name() == math::logic::f_lineal) ?
		(*this == static_cast<const f_lineal&>(op2)) : false;
}
bool math::f_lineal::operator==(const f_lineal& op2) const {
	return ((fabs(_f0 - op2._f0) < constant::TOL()) &&
		    (fabs(_f1 - op2._f1) < constant::TOL()));
}
/* overloaded operator == (equal) */

math::f_lineal* math::f_lineal::clone() const {
	return new f_lineal(*this);
}
/* cloner */

double math::f_lineal::value(const double& input) const {
	return _f0 + input * _f1;
}
/* evaluates the function at the reference magnitude input, and writes the
result at the reference magnitude result. Only the magnitude value is
inserted into result, the units are assummed to be OK. */

double math::f_lineal::d_dt(const double& input, const double& input_dt) const {
	return input_dt * _f1;
}
/* evaluates the function differential with time at the reference
magnitude input and its differential with time input_dt, and writes the
result at the reference magnitude result. */

void math::f_lineal::fill_identity() {
	_identity = (fabs(_f0 - 0.) < constant::TOL()) && (fabs(_f1 - 1.) < constant::TOL());
}
/* sets the _identity flag */



