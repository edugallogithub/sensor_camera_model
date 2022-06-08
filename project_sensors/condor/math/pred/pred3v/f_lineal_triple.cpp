#include "f_lineal_triple.h"

// CLASS F_LINEAL_TRIPLE
// =====================
// =====================

const math::logic::PRED_NAME math::f_lineal_triple::_name = math::logic::f_lineal_triple;
/* predicate name */

const std::string math::f_lineal_triple::_st_name = "f_lineal_triple";
/* predicate name string */

math::f_lineal_triple::f_lineal_triple(const double f0,
									  const double f1x,
									  const double f1y, 
									  const double f1z,
									  const double f1xy, 
									  const double f1xz, 
									  const double f1yz, 
									  const double f1xyz)
 :  _f0(f0), _f1x(f1x), _f1y(f1y), _f1z(f1z),
 _f1xy(f1xy), _f1xz(f1xz), _f1yz(f1yz), _f1xyz(f1xyz) {}
/* constructor based on lineal function coefficients for input and
output in standard units */

math::f_lineal_triple::f_lineal_triple(const f_lineal_triple& other)
: _f0(other._f0), _f1x(other._f1x), _f1y(other._f1y), _f1z(other._f1z),
_f1xy(other._f1xy), _f1xz(other._f1xz), _f1yz(other._f1yz), _f1xyz(other._f1xyz) {}
/* copy constructor */

math::f_lineal_triple* math::f_lineal_triple::clone() const {
	return new f_lineal_triple(*this);
}
/* cloner */

double math::f_lineal_triple::value(const double& input1,
								 const double& input2,
								 const double& input3) const {
	double x = input1;
	double y = input2;
	double z = input3;
	return _f0 + _f1x * x + _f1y * y + _f1z * z +  _f1xy * x * y + _f1xz * x * z +
			_f1yz * y * z + _f1xyz * x * y * z;
}
/* evaluates the function at the reference magnitudes input1, input2, and input3,
and	writes the result at the reference magnitude result. Only the magnitude
value is inserted into result, the units are assummed to be OK. */

double math::f_lineal_triple::d_dt(const double& input1, const double& input2, const double& input3,
                             const double& input1_dt, const double& input2_dt, const double& input3_dt) const {
	return
			input1_dt * (_f1x + input2 * _f1xy + input3 * _f1xz + input2 * input3 * _f1xyz) +
			input2_dt * (_f1y + input1 * _f1xy + input3 * _f1yz + input1 * input3 * _f1xyz) +
			input3_dt * (_f1z + input1 * _f1xz + input2 * _f1yz + input1 * input2 * _f1xyz);
}
/* evaluates the function differential with time at the reference
magnitudes input1, input2, and input3 and their differentials with time
input1_dt, input2_dt, and input3_dt, and writes the result at the reference
magnitude result. */

bool math::f_lineal_triple::operator==(const pred3v& op2) const {
	return (op2.get_name() == math::logic::f_lineal_triple) ?
		(*this == static_cast<const f_lineal_triple&>(op2)) : false;	
}
bool math::f_lineal_triple::operator==(const f_lineal_triple& op2) const {
	return ((fabs(_f0    - op2._f0)    < constant::TOL()) &&
		    (fabs(_f1x   - op2._f1x)   < constant::TOL()) &&
			(fabs(_f1y   - op2._f1y)   < constant::TOL()) &&
			(fabs(_f1z   - op2._f1z)   < constant::TOL()) &&
		    (fabs(_f1xy  - op2._f1xy)  < constant::TOL()) &&
			(fabs(_f1xz  - op2._f1xz)  < constant::TOL()) &&
			(fabs(_f1yz  - op2._f1yz)  < constant::TOL()) &&
			(fabs(_f1xyz - op2._f1xyz) < constant::TOL()));
}
/* overloaded operator == (equal) */








