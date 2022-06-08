#include "f_table1V_spl.h"

// CLASS F_TABLE1V_SPL
// ===================
// ===================

const math::logic::PRED_NAME math::f_table1V_spl::_name = math::logic::f_table1V_spl;
/* predicate name */

const std::string math::f_table1V_spl::_st_name = "f_table1V_spl";
/* predicate name string */

math::f_table1V_spl::f_table1V_spl(std::vector<double>* points1,
								  std::vector<double>* values)
: _points1(points1), _values(values), _diffdiff(0), _dist(0) {
	unsigned short n_points = _points1->size();
	// check size compatibility between _points and _steps
	assert(n_points == _values->size());
	// check that points vector has at least two points
	assert(n_points >= 4);
	// check that points vector is ordered
	for (unsigned short i = 0; i != (n_points-1); ++i) {
		assert((*_points1)[i] < (*_points1)[i+1]);
	}
	// precompute spline
	this->complete_spline(*_points1, *_values);
	// compute differences
	_dist = new std::vector<double>(_diffdiff->size() - 1, 0.);
	for (int i = 0; i != _dist->size(); ++i) {
		(*_dist)[i] = (*_points1)[i+1] - (*_points1)[i];
	}
}
/* constructor based on pointers to the vectors representing the 
table inputs and outputs. Assumes standard units. */

math::f_table1V_spl::f_table1V_spl(const f_table1V_spl& other)
: _points1(new std::vector<double>(*(other._points1))),
  _values(new std::vector<double>(*(other._values))),
  _diffdiff(new std::vector<double>(*(other._diffdiff))),
  _dist(new std::vector<double>(*(other._dist))) {
}
/* copy constructor */

math::f_table1V_spl::~f_table1V_spl() {
	delete _points1;
	delete _values;
	delete _diffdiff;
	delete _dist;
}
/* destructor */

math::f_table1V_spl* math::f_table1V_spl::clone() const {
	return new f_table1V_spl(*this);
}
/* cloner */

bool math::f_table1V_spl::operator==(const pred1v& op2) const {
	return (op2.get_name() == math::logic::f_table1V_spl) ?
		(*this == static_cast<const f_table1V_spl&>(op2)) : false;	
}
bool math::f_table1V_spl::operator==(const f_table1V_spl& op2) const {
	if ((*_points1 != *(op2._points1)) || (*_values != *(op2._values))) {
		return false;
	}
	return true;	
}
/* overloaded operator == (equal) */

void math::f_table1V_spl::complete_spline(const std::vector<double>& points1,
										 const std::vector<double>& values) {
	_diffdiff = new std::vector<double>(points1.size(), 0.);
	int n = points1.size() - 1;
	double yp1 = values[0] / (points1[0] - points1[1]) + 
				 values[0] / (points1[0] - points1[2]) + 
				 values[0] / (points1[0] - points1[3]) + 
				 values[1] * (points1[0] - points1[2]) * (points1[0] - points1[3]) / (points1[1] - points1[0]) / (points1[1] - points1[2]) / (points1[1] - points1[3]) + 
				 values[2] * (points1[0] - points1[1]) * (points1[0] - points1[3]) / (points1[2] - points1[0]) / (points1[2] - points1[1]) / (points1[2] - points1[3]) + 
				 values[3] * (points1[0] - points1[1]) * (points1[0] - points1[2]) / (points1[3] - points1[0]) / (points1[3] - points1[1]) / (points1[3] - points1[2]);
	double ypn = values[n] / (points1[n] - points1[n-3]) + 
				 values[n] / (points1[n] - points1[n-2]) + 
				 values[n] / (points1[n] - points1[n-1]) + 
				 values[n-3] * (points1[n] - points1[n-2]) * (points1[n] - points1[n-1]) / (points1[n-3] - points1[n-2]) / (points1[n-3] - points1[n-1]) / (points1[n-3] - points1[n]) +
				 values[n-2] * (points1[n] - points1[n-3]) * (points1[n] - points1[n-1]) / (points1[n-2] - points1[n-3]) / (points1[n-2] - points1[n-1]) / (points1[n-2] - points1[n]) + 
				 values[n-1] * (points1[n] - points1[n-3]) * (points1[n] - points1[n-2]) / (points1[n-1] - points1[n-3]) / (points1[n-1] - points1[n-2]) / (points1[n-1] - points1[n]); 

	n = points1.size();
	std::vector<double> u(n - 1, 0.);		
	
	(*_diffdiff)[0] = -0.5;
	u[0]  = (3.0 / (points1[1] - points1[0])) * ((values[1] - values[0]) /	
			(points1[1] - points1[0]) - yp1);

	double sig = 0., p = 0.;
	for (int i = 1; i < (n-1); i++) {
		sig = (points1[i] - points1[i-1]) / (points1[i+1] - points1[i-1]);
		p = sig * (*_diffdiff)[i-1] + 2.0;
		(*_diffdiff)[i] = (sig - 1.0) / p;
		u[i] = (values[i+1] - values[i]) / (points1[i+1] - points1[i]) - 
			   (values[i] - values[i-1]) / (points1[i] - points1[i-1]);
		u[i] = (6.0 * u[i] / (points1[i+1] - points1[i-1]) - sig * u[i-1]) / p;
	}
	double qn = 0., un = 0.;

	qn = 0.5;
	un = (3.0 / (points1[n-1] - points1[n-2])) * (ypn - (values[n-1] - values[n-2]) /
		 (points1[n-1] - points1[n-2]));
	
	(*_diffdiff)[n-1] = (un - qn * u[n-2]) / (qn * (*_diffdiff)[n-2] + 1.0);
	for (int k = n-2; k >= 0; k--) {
		(*_diffdiff)[k] = (*_diffdiff)[k] * (*_diffdiff)[k+1] + u[k];
	}								 
}
/* complements the constructor based on two vectors of the same size (n), 
one containing the inputs and another the outputs of the function. */

int math::f_table1V_spl::compute_pos1(const double& input1) const {
	int j = upper_bound(_points1->begin(), _points1->end(), input1) - _points1->begin() - 1;
	int n = _points1->size();
	return std::min(std::max((j - (2 - 1) / 2), 0), n - 2);
}
/* Returns first position within _points1 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

math::ratio* math::f_table1V_spl::compute_ratio1(const double& input1,
													  const int& pos1) const {
	math::ratio_spline* result = math::ratio_mgr::from_pool_spline();
	result->set_ratio12((input1 - (*_points1)[pos1]) / ((*_points1)[pos1+1] - (*_points1)[pos1]));
	return result;
}
/* Returns ratio of input1 with respect to the two points1 identified by pos1 and
pos1+1. */

void math::f_table1V_spl::compute_value(double& result,
									   const int& pos1,
									   const math::ratio& rat1) const {
	double b = static_cast<const math::ratio_spline&>(rat1)._ratio12;
	double a = 1 - b;
	result = a * (*_values)[pos1] + b * (*_values)[pos1+1] +
		   ((a*a*a - a) * (*_diffdiff)[pos1] + (b*b*b - b) * (*_diffdiff)[pos1+1]) * ((*_dist)[pos1] * (*_dist)[pos1]) / 6.0;
}
/* Fills up result magnitude by interpolating based on the input position and 
ratio */

void math::f_table1V_spl::compute_diff(double& result,
									  const int& pos1, 
									  const double& input1,
									  const double& input1_dt) const {
	math::ratio* ratio1 = compute_ratio1(input1, pos1);
    double Omag1; // there is no way of knowing the type of magnitude
    this->compute_value(Omag1, pos1, *ratio1);
    math::ratio_spline* Ratio1 = static_cast<math::ratio_spline*>(ratio1);
    Ratio1->set_ratio12(Ratio1->_ratio12 + 1e-4);
    double difX = 1e-4 * ((*_points1)[pos1+1] - (*_points1)[pos1]);
    double OmagX;	// there is no way of knowing the type of magnitude
    this->compute_value(OmagX, pos1, *ratio1);
    result = (OmagX - Omag1) * input1_dt / difX;
	// release memory
	math::ratio_mgr::to_pool_spline(static_cast<math::ratio_spline*>(ratio1));
}
/* Fills up result differential basedon position and ratio */

double math::f_table1V_spl::value(const double& input1) const {
	int pos1 = compute_pos1(input1);
	math::ratio* ratio1 = compute_ratio1(input1, pos1);
    double res;
	this->compute_value(res, pos1, *ratio1);
	math::ratio_mgr::to_pool_spline(static_cast<math::ratio_spline*>(ratio1));
    return res;
}
/* evaluates the function at the reference magnitude input, and writes the
result at the reference magnitude result. Only the magnitude value is
inserted into result, the units are assummed to be OK. */

double math::f_table1V_spl::d_dt(const double& input1, const double& input1_dt) const {
	int pos1 = compute_pos1(input1);
    double res;
	compute_diff(res, pos1, input1, input1_dt);
    return res;
}
/* evaluates the function differential with time at the reference
magnitude input and its differential with time input_dt, and writes the
result at the reference magnitude result. */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////




