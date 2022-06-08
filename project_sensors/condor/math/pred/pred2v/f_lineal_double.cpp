#include "f_lineal_double.h"

// CLASS F_LINEAL_DOUBLE
// =====================
// =====================

const math::logic::PRED_NAME math::f_lineal_double::_name = math::logic::f_lineal_double;
/* predicate name */

const std::string math::f_lineal_double::_st_name = "f_lineal_double";
/* predicate name string */

math::f_lineal_double::f_lineal_double(const double f0,
									  const double f1x,
									  const double f1y,
									  const double f1xy)
: _f0(f0), _f1x(f1x), _f1y(f1y), _f1xy(f1xy) {}
/* constructor based on lineal function coefficients for input and
output in standard units */

math::f_lineal_double::f_lineal_double(const f_lineal_double& other)
: _f0(other._f0), _f1x(other._f1x), _f1y(other._f1y), _f1xy(other._f1xy) {}
/* copy constructor */

double math::f_lineal_double::value(const double& input1,
								 const double& input2) const {
	double x = input1;
	double y = input2;
	return _f0 + _f1x * x + _f1y * y + _f1xy * x * y;
}
/* evaluates the function at the reference magnitudes input1 and input2,
and	writes the result at the reference magnitude result. Only the magnitude
value is inserted into result, the units are assummed to be OK. */

double math::f_lineal_double::d_dt(const double& input1,
                             const double& input2,
                             const double& input1_dt,
                             const double& input2_dt) const {
	return input1_dt * (_f1x + input2 * _f1xy) + input2_dt * (_f1y + input1 * _f1xy);
}
/* evaluates the function differential with time at the reference
magnitudes input1 and input2 and their differentials with time input1_dt
and input2_dt, and writes the result at the reference magnitude result. */

math::f_lineal_double* math::f_lineal_double::clone() const {
	return new f_lineal_double(*this);
}
/* cloner */

bool math::f_lineal_double::operator==(const pred2v& op2) const {
	return (op2.get_name() == math::logic::f_lineal_double) ?
		(*this == static_cast<const f_lineal_double&>(op2)) : false;
}
bool math::f_lineal_double::operator==(const f_lineal_double& op2) const {
	return ((fabs(_f0   - op2._f0)   < constant::TOL()) &&
		    (fabs(_f1x  - op2._f1x)  < constant::TOL()) &&
			(fabs(_f1y  - op2._f1y)  < constant::TOL()) &&
			(fabs(_f1xy - op2._f1xy) < constant::TOL()));
}
/* overloaded operator == (equal) */


