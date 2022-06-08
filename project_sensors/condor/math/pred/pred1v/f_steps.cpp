#include "f_steps.h"

// CLASS F_STEPS
// =============
// =============

const math::logic::PRED_NAME math::f_steps::_name = math::logic::f_steps;
/* predicate name */

const std::string math::f_steps::_st_name = "f_steps";
/* predicate name string */

math::f_steps::f_steps(vec1* points,
					  vec1* values)
: _points(points), _values(values), _del_flag(true) {
	//check size compatibility between _points and _steps
	if ((points->size1() - 1) != values->size1()) {
		delete _points;
		delete _values;
		throw std::runtime_error("Incorrect size.");
	}
	// points vector has at least two points	
	if (points->size1() < 2) {
		std::string st_name = "a";
		delete _points;
		delete _values;
        throw std::runtime_error("Incorrect size.");
	}
	// check that points vector is ordered
	for (unsigned short i = 0; i != points->size1() - 1; ++i) {
		if (_points->get(i) >= _points->get(i+1)) {
			std::string st_name = "a";
			delete _points;
			delete _values;
            throw std::runtime_error("Incorrect size.");
		}
	}
}
/* constructor based on pointers - deleted by destructor */

math::f_steps::f_steps(vec1& points,
					  vec1& values)
: _points(&points), _values(&values), _del_flag(false) {
	//check size compatibility between _points and _steps
	if ((points.size1() - 1) != values.size1()) {
        throw std::runtime_error("Incorrect size.");
	}
 	// points vector has at least two points	
	if (points.size1() < 2) {
        throw std::runtime_error("Incorrect size.");
	}
	// check that points vector is ordered
	for (unsigned short i = 0; i != points.size1() - 1; ++i) {
		if (_points->get(i) >= _points->get(i+1)) {
            throw std::runtime_error("Incorrect size.");
		}
	}
}
/* constructor based on references - not deleted by destructor */

math::f_steps::f_steps(const f_steps& other)
: _del_flag(other._del_flag) {
	if (_del_flag == true) {
		_points = other._points->clone();
		_values = other._values->clone();
	}
	else {
		_points = other._points;
		_values = other._values;
	}
}
/* copy constructor */

bool math::f_steps::operator==(const pred1v& op2) const {
	return (op2.get_name() == math::logic::f_steps) ?
		(*this == static_cast<const f_steps&>(op2)) : false;	
}
bool math::f_steps::operator==(const f_steps& op2) const {
	return ((*_points == *op2._points) && (*_values == *op2._values));	
}
/* overloaded operator == (equal) */

math::f_steps::~f_steps() {
	if (_del_flag == true) {
		delete _points;
		delete _values;
	}
}
/* destructor */

math::f_steps* math::f_steps::clone() const {
	return new f_steps(*this);
}
/* cloner */

double math::f_steps::value(const double& input) const {
	// Look for proper _point immediately superior to input
	std::vector<double>::const_iterator I = upper_bound(_points->begin(), _points->end(), input);
	int pos = I - _points->begin() - 1;
	if (pos == (-1)) { // out of range below
        return _values->front();
    }
    else if (pos != (_points->size1() - 1)) { // normal situation
        return _values->get(pos);
    }
    else { // out of range above and input coincides with last point
        return _values->back();
    }
}
/* evaluates the function at the reference magnitude input, and writes the
result at the reference magnitude result. Only the magnitude value is
inserted into result, the units are assummed to be OK. */

double math::f_steps::d_dt(const double& input, const double& input_dt) const {
	return 0.;
};
/* evaluates the function differential with time at the reference
magnitude input and its differential with time input_dt, and writes the
result at the reference magnitude result. */




