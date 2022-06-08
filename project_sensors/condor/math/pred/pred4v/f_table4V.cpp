#include "f_table4V.h"

// CLASS F_TABLE4V
// ===============
// ===============

const math::logic::PRED_NAME math::f_table4V::_name = math::logic::f_table4V;
/* predicate name */

const std::string math::f_table4V::_st_name = "f_table4V";
/* predicate name string */

math::f_table4V::f_table4V(vec1* points4,
							  vec1* points3,
							  vec1* points2,
							  vec1* points1,
							  vec4* VVValues,
							  math::logic::INTERP_MODE interp_mode)
: _points4(points4), _points3(points3), _points2(points2), _points1(points1),
_VVValues(VVValues), _del_flag(true), _interp_mode(interp_mode),
 _interp(math::interp::get_interp(interp_mode)),
_equi4(true), _equi3(true), _equi2(true), _equi1(true), _functor_diff(0),
_slopes_d4(0), _slopes_d3(0), _slopes_d2(0), _slopes_d1(0), _herm(0), 
_finder4(0), _finder3(0), _finder2(0), _finder1(0),
_checker4(new math::range_checker_inactive()),
_checker3(new math::range_checker_active()),
_checker2(new math::range_checker_active()),
_checker1(new math::range_checker_active()) {

	// check size compatibility between _points4, _points3, _points2, _points1, and _Values
	if ((points4->size1() != VVValues->size4()) ||
		(points3->size1() != VVValues->size3()) ||
		(points2->size1() != VVValues->size2()) ||
		(points1->size1() != VVValues->size1())) {
		delete _points4;	delete _points3;	delete _points2;	delete _points1;
		delete _VVValues;
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	// points vectors have at least minimum number of points
	if (points4->size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _points4;	delete _points3;	delete _points2;	delete _points1;
		delete _VVValues;
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points3->size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _points4;	delete _points3;	delete _points2;	delete _points1;
		delete _VVValues;
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points2->size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _points4;	delete _points3;	delete _points2;	delete _points1;
		delete _VVValues;
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points1->size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _points4;	delete _points3;	delete _points2;	delete _points1;
		delete _VVValues;
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	// difference between consecutive _points4 and _points3
	_points4_diff = _points4->get(1) - _points4->get(0);
	_points3_diff = _points3->get(1) - _points3->get(0);
	_points2_diff = _points2->get(1) - _points2->get(0);
	_points1_diff = _points1->get(1) - _points1->get(0);
	double temp = 0.;
	// difference between consecutive points4
	for (unsigned short i = 1; i != _points4->size1() - 1; ++i) {
		if (_equi4 == true) { // so far equispaced
			temp = _points4->get(i+1) - _points4->get(i);
			if (fabs((temp-_points4_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi4 = false;
				_points4_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _points4;	delete _points3;	delete _points2;	delete _points1;
					delete _VVValues;
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if (_points4->get(i+1) - _points4->get(i) <= 0.) {
				std::string st_name = "a";
				delete _points4;		delete _points3;	delete _points2;	delete _points1;
				delete _VVValues;
				delete _interp;	
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points3
	for (unsigned short i = 1; i != _points3->size1() - 1; ++i) {
		if (_equi3 == true) { // so far equispaced
			temp = _points3->get(i+1) - _points3->get(i);
			if (fabs((temp-_points3_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi3 = false;
				_points3_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _points4;	delete _points3;	delete _points2;	delete _points1;
					delete _VVValues;
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if (_points3->get(i+1) - _points3->get(i) <= 0.) {
				std::string st_name = "a";
				delete _points4;		delete _points3;	delete _points2;	delete _points1;
				delete _VVValues;
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points2
	for (unsigned short i = 1; i != _points2->size1() - 1; ++i) {
		if (_equi2 == true) { // so far equispaced
			temp = _points2->get(i+1) - _points2->get(i);
			if (fabs((temp-_points2_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi2 = false;
				_points2_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _points4;	delete _points3;	delete _points2;	delete _points1;
					delete _VVValues;
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if (_points2->get(i+1) - _points2->get(i) <= 0.) {
				std::string st_name = "a";
				delete _points4;		delete _points3;	delete _points2;	delete _points1;
				delete _VVValues;
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points1
	for (unsigned short i = 1; i != _points1->size1() - 1; ++i) {
		if (_equi1 == true) { // so far equispaced
			temp = _points1->get(i+1) - _points1->get(i);
			if (fabs((temp-_points1_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi1 = false;
				_points1_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _points4;	delete _points3;	delete _points2;	delete _points1;
					delete _VVValues;
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if (_points1->get(i+1) - _points1->get(i) <= 0.) {
				std::string st_name = "a";
				delete _points4;		delete _points3;	delete _points2;	delete _points1;
				delete _VVValues;
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}

	if (_equi4 == true) {_finder4 = new math::pos_finder_equispaced(*_points4, *_interp, _points4_diff);}
	else                {_finder4 = new math::pos_finder_binary(*_points4, *_interp);}
	if (_equi3 == true) {_finder3 = new math::pos_finder_equispaced(*_points3, *_interp, _points3_diff);}
	else                {_finder3 = new math::pos_finder_binary(*_points3, *_interp);}
	if (_equi2 == true) {_finder2 = new math::pos_finder_equispaced(*_points2, *_interp, _points2_diff);}
	else                {_finder2 = new math::pos_finder_binary(*_points2, *_interp);}
	if (_equi1 == true) {_finder1 = new math::pos_finder_equispaced(*_points1, *_interp, _points1_diff);}
	else                {_finder1 = new math::pos_finder_binary(*_points1, *_interp);}

	switch (_interp_mode) {
		case math::logic::lagrange_first_precompute:
			fill_up_slopes_lagrange_first_precompute();
			_herm = new hermite4v();
			_functor_diff = new math::table4V_diff_prec(*this);
			break;
		case math::logic::hermite_first:
		case math::logic::hermite_second:
			fill_up_slopes_hermite();
			_herm = new hermite4v(*_points4, *_points3, *_points2, *_points1, *_VVValues, *_slopes_d4, *_slopes_d3, *_slopes_d2, *_slopes_d1, _interp_mode);
			_functor_diff = new math::table4V_diff_real(*this);
			break;
		case math::logic::spline: // splines not allowed in 3 dimensions
			delete _points4;	delete _points3;	delete _points2;	delete _points1;
			delete _VVValues;
			delete _interp;
			delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
			delete _finder4;	delete _finder3;	delete _finder2;	delete _finder1;
            throw std::runtime_error("Incorrect interpolation method.");
			break;
		default:
			_slopes_d4 = 0;
			_slopes_d3 = 0;
			_slopes_d2 = 0;
			_slopes_d1 = 0;
			_herm = new hermite4v();
			_functor_diff = new math::table4V_diff_real(*this);
			break;
	}
}
/* constructor based on pointer to size l vector points4, pointer	to a
size m vector points3, pointer to a size n vector points2, pointer to a
size q vector point3, and pointer to a size l vector of size m vectors of size n
vectors of size q vectors VVValues - deleted by destructor. */

math::f_table4V::f_table4V(vec1& points4,
							  vec1& points3,
							  vec1& points2,
							  vec1& points1,
							  vec4& VVValues,
							  math::logic::INTERP_MODE interp_mode)
: _points4(&points4), _points3(&points3), _points2(&points2), _points1(&points1),
_VVValues(&VVValues), _del_flag(false), _interp_mode(interp_mode), 
_interp(math::interp::get_interp(interp_mode)),
_equi4(true), _equi3(true), _equi2(true), _equi1(true), _functor_diff(0),
_slopes_d4(0), _slopes_d3(0), _slopes_d2(0), _slopes_d1(0), _herm(0),
_finder4(0), _finder3(0), _finder2(0), _finder1(0),
_checker4(new math::range_checker_inactive()),
_checker3(new math::range_checker_active()),
_checker2(new math::range_checker_active()),
_checker1(new math::range_checker_active()) {
	// check size compatibility between _points4, _points3, _points2, _points1 and _Values
	if ((points4.size1() != VVValues.size4()) ||
		(points3.size1() != VVValues.size3()) ||
		(points2.size1() != VVValues.size2()) ||
		(points1.size1() != VVValues.size1())) {
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	// points vectors have at least minimum number of points
	if (points4.size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points3.size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points2.size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	if (points1.size1() < _interp->get_min_points()) {
		int n = _interp->get_min_points();
		delete _interp;
		delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
        throw std::runtime_error("Incorrect size.");
	}
	// difference between consecutive _points4 and _points3
	_points4_diff = _points4->get(1) - _points4->get(0);
	_points3_diff = _points3->get(1) - _points3->get(0);
	_points2_diff = _points2->get(1) - _points2->get(0);
	_points1_diff = _points1->get(1) - _points1->get(0);
	double temp = 0.;
	// difference between consecutive points4
	for (unsigned short i = 1; i != _points4->size1() - 1; ++i) {
		if (_equi4 == true) { // so far equispaced
			temp = _points4->get(i+1) - _points4->get(i);
			if (fabs((temp-_points4_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi4 = false;
				_points4_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if ((_points4->get(i+1) - _points4->get(i)) <= 0.) {
				std::string st_name = "a";
				delete _interp;	
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points3
	for (unsigned short i = 1; i != _points3->size1() - 1; ++i) {
		if (_equi3 == true) { // so far equispaced
			temp = _points3->get(i+1) - _points3->get(i);
			if (fabs((temp-_points3_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi3 = false;
				_points3_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if ((_points3->get(i+1) - _points3->get(i)) <= 0.) {
				std::string st_name = "a";
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points2
	for (unsigned short i = 1; i != _points2->size1() - 1; ++i) {
		if (_equi2 == true) { // so far equispaced
			temp = _points2->get(i+1) - _points2->get(i);
			if (fabs((temp-_points2_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi2 = false;
				_points2_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if ((_points2->get(i+1) - _points2->get(i)) <= 0.) {
				std::string st_name = "a";
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}
	// difference between consecutive points1
	for (unsigned short i = 1; i != _points1->size1() - 1; ++i) {
		if (_equi1 == true) { // so far equispaced
			temp = _points1->get(i+1) - _points1->get(i);
			if (fabs((temp-_points1_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi1 = false;
				_points1_diff = 0.; 
				if (temp <= 0.) {
					std::string st_name = "a";
					delete _interp;
					delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
                    throw std::runtime_error("Incorrect size.");
				}
			}	
		}
		else { // so far not equispaced and will remain that way
			if ((_points1->get(i+1) - _points1->get(i)) <= 0.) {
				std::string st_name = "a";
				delete _interp;
				delete _checker4;		delete _checker3;	delete _checker2;	delete _checker1;
                throw std::runtime_error("Incorrect size.");
			}
		}
	}

	if (_equi4 == true) {_finder4 = new math::pos_finder_equispaced(*_points4, *_interp, _points4_diff);}
	else                {_finder4 = new math::pos_finder_binary(*_points4, *_interp);}
	if (_equi3 == true) {_finder3 = new math::pos_finder_equispaced(*_points3, *_interp, _points3_diff);}
	else                {_finder3 = new math::pos_finder_binary(*_points3, *_interp);}
	if (_equi2 == true) {_finder2 = new math::pos_finder_equispaced(*_points2, *_interp, _points2_diff);}
	else                {_finder2 = new math::pos_finder_binary(*_points2, *_interp);}
	if (_equi1 == true) {_finder1 = new math::pos_finder_equispaced(*_points1, *_interp, _points1_diff);}
	else                {_finder1 = new math::pos_finder_binary(*_points1, *_interp);}

	switch (_interp_mode) {
		case math::logic::lagrange_first_precompute:
			fill_up_slopes_lagrange_first_precompute();
			_herm = new hermite4v();
			_functor_diff = new math::table4V_diff_prec(*this);
			break;
		case math::logic::hermite_first:
		case math::logic::hermite_second:
			fill_up_slopes_hermite();
			_herm = new hermite4v(*_points4, *_points3, *_points2, *_points1, *_VVValues, *_slopes_d4, *_slopes_d3, *_slopes_d2, *_slopes_d1, _interp_mode);
			_functor_diff = new math::table4V_diff_real(*this);
			break;
		case math::logic::spline: // splines not allowed in 3 dimensions
			delete _interp;
			delete _checker4;	delete _checker3;	delete _checker2;	delete _checker1;
			delete _finder4;	delete _finder3;	delete _finder2;	delete _finder1;
            throw std::runtime_error("Incorrect interpolation method.");
			break;
		default:
			_slopes_d4 = 0;
			_slopes_d3 = 0;
			_slopes_d2 = 0;
			_slopes_d1 = 0;
			_herm = new hermite4v();
			_functor_diff = new math::table4V_diff_real(*this);
			break;
	}
}
/* constructor based on reference to size l vector points4, reference	to a
size m vector points3, reference to a size n vector points2, reference to a
size q vector points1, and reference to	a size l vector of size m vectors of 
size n vectos of size q vectors VValues - not deleted by destructor. */

math::f_table4V::f_table4V(const f_table4V& other)
: _del_flag(other._del_flag),
_points4_diff(other._points4_diff), _equi4(other._equi4),
_points3_diff(other._points3_diff), _equi3(other._equi3),
_points2_diff(other._points2_diff), _equi2(other._equi2),
_points1_diff(other._points1_diff), _equi1(other._equi1),
_slopes_d4(0), _slopes_d3(0), _slopes_d2(0), _slopes_d1(0), _herm(0), 
_interp_mode(other._interp_mode),
_interp(math::interp::get_interp(_interp_mode)),
_checker4(other._checker4->clone()),
_checker3(other._checker3->clone()),
_checker2(other._checker2->clone()), 
_checker1(other._checker1->clone()),
_finder4(other._finder4->clone()),
_finder3(other._finder3->clone()),
_finder2(other._finder2->clone()),
_finder1(other._finder1->clone()) {
	if (other._herm != 0) {
		_herm = new hermite4v(*other._herm);
	}
	if (other._slopes_d4 != 0) {
        _slopes_d4 = new vec4(*other._slopes_d4);
	}
	if (other._slopes_d3 != 0) {
		_slopes_d3 = new vec4(*other._slopes_d3);
	}
	if (other._slopes_d2 != 0) {
		_slopes_d2 = new vec4(*other._slopes_d2);
	}
	if (other._slopes_d1 != 0) {
		_slopes_d1 = new vec4(*other._slopes_d1);
	}		
	if (_del_flag == true) {
		_points4=other._points4->clone();
		_points3=other._points3->clone();
		_points2=other._points2->clone();
		_points1=other._points1->clone();
		_VVValues=other._VVValues->clone();
	}
	else {
		_points4 = other._points4;
		_points3 = other._points3;
		_points2 = other._points2;
		_points1 = other._points1;
		_VVValues = other._VVValues;
	}

	switch (_interp_mode) {
		case math::logic::lagrange_first_precompute:
			_functor_diff = new math::table4V_diff_prec(*this);
			break;
		case math::logic::hermite_first:
		case math::logic::hermite_second:
			_functor_diff = new math::table4V_diff_real(*this);
			break;
		default:
			_functor_diff = new math::table4V_diff_real(*this);
			break;
	}
}
/* copy constructor */

math::f_table4V* math::f_table4V::clone() const {
	return new f_table4V(*this);
}
/* cloner */

bool math::f_table4V::operator==(const pred4v& op2) const {
	return (op2.get_name() == math::logic::f_table4V) ?
		(*this == static_cast<const f_table4V&>(op2)) : false;	
}
bool math::f_table4V::operator==(const f_table4V& op2) const {
	return ((*_points4 == *op2._points4) &&	(*_points3 == *op2._points3) &&
			(*_points2 == *op2._points2) &&	(*_points1 == *op2._points1) &&
			(*_VVValues == *op2._VVValues) && (_interp_mode == op2._interp_mode));	
}
/* overloaded operator == (equal) */

math::f_table4V::~f_table4V() {
	if (_del_flag == true) {
		delete _points4;
		delete _points3;
		delete _points2;
		delete _points1;
		delete _VVValues;
	}
	delete _interp;
	delete _functor_diff;
	delete _slopes_d4;
	delete _slopes_d3;
	delete _slopes_d2;
	delete _slopes_d1;
	delete _herm;
	delete _checker4;
	delete _checker3;
	delete _checker2;
	delete _checker1;
	delete _finder4;
	delete _finder3;
	delete _finder2;
	delete _finder1;
}
/* destructor */

int math::f_table4V::compute_pos4(const double& input4) const	{
	return _interp->find_index(_finder4->search(*_points4, input4, _points4_diff),
								_points4->size1());
}
/* Returns first position within _points4 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

int math::f_table4V::compute_pos3(const double& input3) const {
	return _interp->find_index(_finder3->search(*_points3, input3, _points3_diff),
								_points3->size1());
}
/* Returns first position within _points3 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

int math::f_table4V::compute_pos2(const double& input2) const {
	return _interp->find_index(_finder2->search(*_points2, input2, _points2_diff),
								_points2->size1());
}
/* Returns first position within _points2 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

int math::f_table4V::compute_pos1(const double& input1) const {
	return _interp->find_index(_finder1->search(*_points1, input1, _points1_diff),
								_points1->size1());
}
/* Returns first position within _points1 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

math::ratio* math::f_table4V::compute_ratio4(const double& input4,
													const int& pos4) const {
	return _interp->compute_ratio(input4, *_points4, pos4);
}
/* Returns ratio of input4 with respect to the two points4 identified by pos4 and pos4+1. */

math::ratio* math::f_table4V::compute_ratio3(const double& input3,
													const int& pos3) const {
	return _interp->compute_ratio(input3, *_points3, pos3);
}
/* Returns ratio of input3 with respect to the two points3 identified by pos3 and pos3+1. */

math::ratio* math::f_table4V::compute_ratio2(const double& input2,
													const int& pos2) const {
	return _interp->compute_ratio(input2, *_points2, pos2);
}
/* Returns ratio of input2 with respect to the two points2 identified by pos2 and pos2+1. */

math::ratio* math::f_table4V::compute_ratio1(const double& input1,
													const int& pos1) const {
	return _interp->compute_ratio(input1, *_points1, pos1);
}
/* Returns ratio of input1 with respect to the two points1 identified by pos1 and pos1+1. */

double math::f_table4V::compute_value(const int& pos4,
									 const int& pos3,
									 const int& pos2,
									 const int& pos1,
									 const math::ratio& ratio4,
									 const math::ratio& ratio3,
									 const math::ratio& ratio2,
									 const math::ratio& ratio1) const {
	double res;
	_interp->interp4(res, *_points4, *_points3, *_points2, *_points1, *_VVValues, pos4, pos3, pos2, pos1, ratio4, ratio3, ratio2, ratio1, *_herm);
    return res;
}
/* Fills up result magnitude by interpolating based on the input positions
and ratios */

double math::f_table4V::compute_diff(const int& pos4,
									const int& pos3,
									const int& pos2,
									const int& pos1,
									const double& input4,
									const double& input3,
									const double& input2,
									const double& input1,
									const double& input4_dt,
									const double& input3_dt,
									const double& input2_dt,
									const double& input1_dt) const {
	double res;
	_functor_diff->compute_diff(res, pos4, pos3, pos2, pos1,
						input4, input3, input2, input1, input4_dt, input3_dt, input2_dt, input1_dt);
    return res;
}
/* Fills up result differential based on positions and ratios */

double math::f_table4V::value(const double& input4,
							 const double& input3,
							 const double& input2,
							 const double& input1) const {
	// verify inputs are within range
	_checker4->check_range(*_points4, input4);
	_checker3->check_range(*_points3, input3);
	_checker2->check_range(*_points2, input2);
	_checker1->check_range(*_points1, input1);
							 
	// points4, points3, points2, and points1 are equispaced
	// pos4 provides the lower 3Dmatrix position
	// pos3 provides, for each 3Dmatrix, the lower 2D matrix position
	// pos2 provides, for each 2Dmatrix, the lower vector position
	// pos1 provides, for each vector, the lower position
	int pos4 = compute_pos4(input4);
	int pos3 = compute_pos3(input3);
	int pos2 = compute_pos2(input2);
	int pos1 = compute_pos1(input1);
	math::ratio* ratio4 = compute_ratio4(input4, pos4);
	math::ratio* ratio3 = compute_ratio3(input3, pos3);
	math::ratio* ratio2 = compute_ratio2(input2, pos2);
	math::ratio* ratio1 = compute_ratio1(input1, pos1);
    double res;
	this->_interp->interp4(res, *_points4, *_points3, *_points2, *_points1, *_VVValues, pos4, pos3, pos2, pos1, *ratio4, *ratio3, *ratio2, *ratio1, *_herm);
	_interp->to_pool(ratio4);
	_interp->to_pool(ratio3);
	_interp->to_pool(ratio2);
	_interp->to_pool(ratio1);
    return res;
}
/* evaluates the function at the reference magnitudes input4, input3, input2,
and input1, and writes the result at the reference magnitude result. Only the
magnitude value is inserted into result, the units are assummed to be OK. */

double math::f_table4V::d_dt(const double& input4,
							const double& input3,
							const double& input2,
							const double& input1,
							const double& input4_dt,
							const double& input3_dt,
							const double& input2_dt,
							const double& input1_dt) const {
	////////////////////////////////////////////////////////////////////
	// The computation of the four positions, which are quite expensive,
	// in theory have already been done before in the value method, and
	// should not be repeated here. However, there is no way to avoid it.
	////////////////////////////////////////////////////////////////////
	// verify inputs are within range
	_checker4->check_range(*_points4, input4);
	_checker3->check_range(*_points3, input3);
	_checker2->check_range(*_points2, input2);
	_checker1->check_range(*_points1, input1);

	// points4, points3, points2 and points1 are equispaced
	// pos4 provides the lower 3Dmatrix position
	// pos3 provides, for each 3Dmatrix, the lower 2D matrix position
	// pos2 provides, for each 2Dmatrix, the lower vector position
	// pos1 provides, for each vector, the lower position
	int pos4 = compute_pos4(input4);
	int pos3 = compute_pos3(input3);
	int pos2 = compute_pos2(input2);
	int pos1 = compute_pos1(input1);
    double res;
	this->_functor_diff->compute_diff(res, pos4, pos3, pos2, pos1,
						input4, input3, input2, input1, input4_dt, input3_dt, input2_dt, input1_dt);
    return res;
}
/* evaluates the function differential with time at the reference
magnitudes input4, input3, input2, and input1 and their differentials with time
input4_dt, input3_dt, input2_dt, and input1_dt and writes the result at the reference
magnitude result. */

void math::f_table4V::fill_up_slopes_lagrange_first_precompute() {
	// differentials with respect to first independent magnitude
	// has size (l) x (m) x (n) x (q)
	_slopes_d4 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	for (unsigned short j = 0; j != _slopes_d4->size4()-1; ++j) { // last one is meaningless and remains 0
		fill_up_slopes_aux4_lagrange_first_precompute(*_slopes_d4, j, _points3->size1(), _points2->size1(), _points1->size1());
	}
	for (unsigned short i = 0; i != _points3->size1(); ++i) {
		for (unsigned short j = 0; j != _points2->size1(); ++j) {
			for (unsigned short k = 0; k != _points1->size1(); ++k) {
				_slopes_d4->set(k, j, i, _slopes_d4->size4()-1, 0.);
			}
		}
	}

	// differentials with respect to second independent magnitude
	// has size (l) x (m) x (n) x (q)
	_slopes_d3 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	for (unsigned short j = 0; j != _slopes_d3->size4(); ++j) {
		fill_up_slopes_aux3_lagrange_first_precompute(*_slopes_d3, j, _points3->size1(), _points2->size1(), _points1->size1());
	}
	// differentials with respect to third independent magnitude
	// has size (l) x (m) x (n) x (q)
	_slopes_d2 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	for (unsigned short j = 0; j != _slopes_d2->size4(); ++j) {
		fill_up_slopes_aux2_lagrange_first_precompute(*_slopes_d2, j, _points3->size1(), _points2->size1(), _points1->size1());
	}
	// differentials with respect to fourth independent magnitude
	// has size (l) x (m) x (n) x (q)
	_slopes_d1 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	for (unsigned short j = 0; j != _slopes_d1->size4(); ++j) {
		fill_up_slopes_aux1_lagrange_first_precompute(*_slopes_d1, j, _points3->size1(), _points2->size1(), _points1->size1());
	}
}
/* fills up the _slopes_d4, _slopes_d3, _slopes_d2, and _slopes_d1 attributes based on
_points4, _points3, _points2, _points1, and _FVValues */

void math::f_table4V::fill_up_slopes_aux4_lagrange_first_precompute
				(math::vec4& slopes_d4,
				 unsigned short pos,
				 unsigned short siz2,
				 unsigned short siz3,
				 unsigned short siz4) {
    for (unsigned short i = 0; i != siz2; ++i) {
        for (unsigned short j = 0; j != siz3; ++j) {
            for (unsigned short k = 0; k != siz4; ++k) {
                slopes_d4.set(k, j, i, pos, (_VVValues->get(k, j, i, pos+1) - _VVValues->get(k, j, i, pos)) / _finder4->diff(pos+1, pos));
            }
        }
    }
}
/* fills up order "pos" of the _slopes_d4 4Dmatrix of differentials with respect to the
first independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux3_lagrange_first_precompute
						(math::vec4& slopes_d3,
						 unsigned short pos,
						 unsigned short siz2,
						 unsigned short siz3,
						 unsigned short siz4) {
    for (unsigned short i = 0; i != siz2-1; ++i) {
        for (unsigned short j = 0; j != siz3; ++j) {
            for (unsigned short k = 0; k != siz4; ++k) {
                slopes_d3.set(k, j, i, pos,
                              (_VVValues->get(k, j, i+1, pos) - _VVValues->get(k, j, i, pos)) / _finder3->diff(i+1,i));
            }
        }
    }
	for (unsigned short j = 0; j != siz3; ++j) {
		for (unsigned short k = 0; k != siz4; ++k) {
			slopes_d3.set(k, j, siz2-1, pos, 0.); // last value meaningless
		}
	}							 
}
/* fills up order "pos" of the _slopes_d3 4Dmatrix of differentials with respect
to the second independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux2_lagrange_first_precompute
						(math::vec4& slopes_d2,
						 unsigned short pos, 
						 unsigned short siz2, 
						 unsigned short siz3,
						 unsigned short siz4) {
    for (unsigned short i = 0; i != siz2; ++i) {
        for (unsigned short j = 0; j != siz3-1; ++j) {
            for (unsigned short k = 0; k != siz4; ++k) {
                slopes_d2.set(k, j, i, pos, (_VVValues->get(k, j+1, i, pos) - _VVValues->get(k, j, i, pos)) / _finder2->diff(j+1,j));
            }
        }
        for (unsigned short k = 0; k != siz4; ++k) {
            slopes_d2.set(k, siz3-1, i, pos, 0.); // last value meaningless
        }
    }
}
/* fills up order "pos" of the _slopes_d2 4Dmatrix of differentials with respect
to the third independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux1_lagrange_first_precompute
						(math::vec4& slopes_d1,
						 unsigned short pos, 
						 unsigned short siz2, 
						 unsigned short siz3,
						 unsigned short siz4) {
	for (unsigned short i = 0; i != siz2; ++i) {
        for (unsigned short j = 0; j != siz3; ++j) {
            for (unsigned short k = 0; k != siz4-1; ++k) {
                slopes_d1.set(k, j, i, pos, (_VVValues->get(k+1, j, i, pos) - _VVValues->get(k, j, i, pos)) / _finder1->diff(k+1,k));
            }
            slopes_d1.set(siz4-1, j, i, pos, 0.); // last value is meaningless
        }
    }
}
/* fills up order "pos" of the _slopes_d1 4Dmatrix of differentials with respect
to the fourth independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_hermite() {
	// differentials with respect to each of the four independent magnitudes
	// have size (l) x (m) x (n) x (q)
	_slopes_d4 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	_slopes_d3 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	_slopes_d2 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());
	_slopes_d1 = new vec4(_points1->size1(), _points2->size1(), _points3->size1(), _points4->size1());

	for (unsigned short i = 0; i != _points3->size1(); ++i) {
		for (unsigned short j = 0; j != _points2->size1(); ++j) {
			for (unsigned short k = 0; k != _points1->size1(); ++k) {
				fill_up_slopes_aux4_hermite(*_slopes_d4, i, j, k); 
			}
		}
	}

	for (unsigned short i = 0; i != _points4->size1(); ++i) {
		for (unsigned short j = 0; j != _points2->size1(); ++j) {
			for (unsigned short k = 0; k != _points1->size1(); ++k) {
				fill_up_slopes_aux3_hermite(*_slopes_d3, i, j, k); 
			}
		}
	}

	for (unsigned short i = 0; i != _points4->size1(); ++i) {
		for (unsigned short j = 0; j != _points3->size1(); ++j) {
			for (unsigned short k = 0; k != _points1->size1(); ++k) {
				fill_up_slopes_aux2_hermite(*_slopes_d2, i, j, k);
			}
		}
	}

	for (unsigned short i = 0; i != _points4->size1(); ++i) {
		for (unsigned short j = 0; j != _points3->size1(); ++j) {
			for (unsigned short k = 0; k != _points1->size1(); ++k) {
				fill_up_slopes_aux1_hermite(*_slopes_d1, i, j, k);
			}
		}
	}
}
/* fills up the _slopes_d4, _slopes_d3, _slopes_d2, and _slopes_d1 attributes based on
_points4, _points3, _points2, _points1, _and _VVValues */

void math::f_table4V::fill_up_slopes_aux4_hermite(math::vec4& slopes_d4,
												   unsigned short index3,
												   unsigned short index2,
												   unsigned short index1) {
	switch (_interp_mode) {
		case math::logic::hermite_first: {
			// compute temporary slopes as in Lagrange first order for each interval 
			vec1 temp(_points4->size1()-1);
			for (unsigned short k = 0; k != _points4->size1()-1; ++k) {
				temp[k] = (_VVValues->get(index1, index2, index3, k+1) - _VVValues->get(index1, index2, index3, k)) / _finder4->diff(k+1,k);
			}
			// compute slopes by average of previous slopes
			slopes_d4.set(index1, index2, index3, 0, temp[0]);
			for (unsigned short k = 1; k != _points4->size1()-1; ++k) {
				slopes_d4.set(index1, index2, index3, k, 0.5 * (temp[k-1] + temp[k]));
			}
			slopes_d4.set(index1, index2, index3, _points4->size1()-1, temp.back());
			break; }
		case math::logic::hermite_second: {
			// compute slopes by results of Lagrange 2nd order interpolation
			unsigned short n = _points4->size1();
			slopes_d4.set(index1, index2, index3, 0,
				 _VVValues->get(index1, index2, index3, 1) / _finder4->diff(2, 0, 1, 0, 2, 1) -
				 _VVValues->get(index1, index2, index3, 2) / _finder4->diff(1, 0, 2, 0, 2, 1) -
				 _VVValues->get(index1, index2, index3, 0) / _finder4->diff(2 ,0) - 
				 _VVValues->get(index1, index2, index3, 0) / _finder4->diff(1 ,0));
			for (unsigned short i = 1; i < n-1; ++i) {
				slopes_d4.set(index1, index2, index3, i,
					_VVValues->get(index1, index2, index3, i+1) / _finder4->diff(i,   i-1, i+1, i-1, i+1, i) - 
					_VVValues->get(index1, index2, index3, i-1) / _finder4->diff(i+1, i,   i,   i-1, i+1, i-1) + 			   
					_VVValues->get(index1, index2, index3, i)   / _finder4->diff(i,   i-1) -  
					_VVValues->get(index1, index2, index3, i)   / _finder4->diff(i+1, i));
			}
			slopes_d4.set(index1, index2, index3, n-1, 
				_VVValues->get(index1, index2, index3, n-3) / _finder4->diff(n-1, n-2, n-2, n-3, n-1, n-3) -
				_VVValues->get(index1, index2, index3, n-2) / _finder4->diff(n-1, n-3, n-2, n-3, n-1, n-2) +
				_VVValues->get(index1, index2, index3, n-1) / _finder4->diff(n-1, n-3) +
				_VVValues->get(index1, index2, index3, n-1) / _finder4->diff(n-1, n-2));	
			break; }
		default:
			throw std::runtime_error("Incorrect interpolation method.");
	}
}
/* fills up order "pos" of the _slopes_d4 4Dmatrix of differentials with respect to the
first independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux3_hermite(math::vec4& slopes_d3,
												   unsigned short index4,
												   unsigned short index2,
												   unsigned short index1) {
	switch (_interp_mode) {
		case math::logic::hermite_first: {
			// compute temporary slopes as in Lagrange first order for each interval 
			vec1 temp(_points3->size1()-1);
			for (unsigned short k = 0; k != _points3->size1()-1; ++k) {
				temp[k] = (_VVValues->get(index1, index2, k+1, index4) - _VVValues->get(index1, index2, k, index4)) / _finder3->diff(k+1,k);
			}
			// compute slopes by average of previous slopes
			slopes_d3.set(index1, index2, 0, index4, temp[0]);
			for (unsigned short k = 1; k != _points3->size1()-1; ++k) {
				slopes_d3.set(index1, index2, k, index4, 0.5 * (temp[k-1] + temp[k]));
			}
			slopes_d3.set(index1, index2, _points3->size1()-1, index4, temp.back());
			break; }
		case math::logic::hermite_second: {
			// compute slopes by results of Lagrange 2nd order interpolation
			unsigned short n = _points3->size1();
			slopes_d3.set(index1, index2, 0, index4,
				 _VVValues->get(index1, index2, 1, index4) / _finder3->diff(2, 0, 1, 0, 2, 1) -
				 _VVValues->get(index1, index2, 2, index4) / _finder3->diff(1, 0, 2, 0, 2, 1) -
				 _VVValues->get(index1, index2, 0, index4) / _finder3->diff(2, 0) - 
				 _VVValues->get(index1, index2, 0, index4) / _finder3->diff(1, 0));
			for (unsigned short i = 1; i < n-1; ++i) {
				slopes_d3.set(index1, index2, i, index4,
				_VVValues->get(index1, index2, i+1, index4) / _finder3->diff(i,   i-1, i+1, i-1, i+1, i) -  
				_VVValues->get(index1, index2, i-1, index4) / _finder3->diff(i+1, i,   i,   i-1, i+1, i-1) + 		 			   
				_VVValues->get(index1, index2, i, index4)   / _finder3->diff(i,   i-1) -  
				_VVValues->get(index1, index2, i, index4)   / _finder3->diff(i+1, i));
			}
			slopes_d3.set(index1, index2, n-1, index4,
				_VVValues->get(index1, index2, n-3, index4) / _finder3->diff(n-1, n-2, n-2, n-3, n-1, n-3) -
				_VVValues->get(index1, index2, n-2, index4) / _finder3->diff(n-1, n-3, n-2, n-3, n-1, n-2) +
				_VVValues->get(index1, index2, n-1, index4) / _finder3->diff(n-1, n-3) +
				_VVValues->get(index1, index2, n-1, index4) / _finder3->diff(n-1, n-2));
			break; }
		default:
			throw std::runtime_error("Incorrect interpolation method.");
	}
}
/* fills up order "pos" of the _slopes_d3 4Dmatrix of differentials with respect
to the second independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux2_hermite(math::vec4& slopes_d2,
												   unsigned short index4,
												   unsigned short index3,
												   unsigned short index1) {
	switch (_interp_mode) {
		case math::logic::hermite_first: {
			// compute temporary slopes as in Lagrange first order for each interval 
			vec1 temp(_points2->size1()-1);
			for (unsigned short k = 0; k != _points2->size1()-1; ++k) {
				temp[k] = (_VVValues->get(index1, k+1, index3, index4) - _VVValues->get(index1, k, index3, index4)) / _finder2->diff(k+1,k);
			}
			// compute slopes by average of previous slopes
			slopes_d2.set(index1, 0, index3, index4, temp[0]);
			for (unsigned short k = 1; k != _points2->size1()-1; ++k) {
				slopes_d2.set(index1, k, index3, index4, 0.5 * (temp[k-1] + temp[k]));
			}
			slopes_d2.set(index1, _points2->size1()-1, index3, index4, temp.back());
		break; }
		case math::logic::hermite_second: {
			// compute slopes by results of Lagrange 2nd order interpolation
			unsigned short n = _points2->size1();
			slopes_d2.set(index1, 0, index3, index4, 
				 _VVValues->get(index1, 1, index3, index4) / _finder2->diff(2, 0, 1, 0, 2, 1) -
				 _VVValues->get(index1, 2, index3, index4) / _finder2->diff(1, 0, 2, 0, 2, 1) -
				 _VVValues->get(index1, 0, index3, index4) / _finder2->diff(2, 0) -  
				 _VVValues->get(index1, 0, index3, index4) / _finder2->diff(1, 0));
			for (unsigned short i = 1; i < n-1; ++i) {
				slopes_d2.set(index1, i, index3, index4, 
				_VVValues->get(index1, i+1, index3, index4) / _finder2->diff(i,   i-1, i+1, i-1, i+1, i) -  
				_VVValues->get(index1, i-1, index3, index4) / _finder2->diff(i+1, i,   i,   i-1, i+1, i-1) + 		 			   
				_VVValues->get(index1, i, index3, index4)   / _finder2->diff(i,   i-1) -  
				_VVValues->get(index1, i, index3, index4)   / _finder2->diff(i+1, i));
			}
			slopes_d2.set(index1, n-1, index3, index4, 
				_VVValues->get(index1, n-3, index3, index4) / _finder2->diff(n-1, n-2, n-2, n-3, n-1, n-3) -
				_VVValues->get(index1, n-2, index3, index4) / _finder2->diff(n-1, n-3, n-2, n-3, n-1, n-2) +
				_VVValues->get(index1, n-1, index3, index4) / _finder2->diff(n-1, n-3) +
				_VVValues->get(index1, n-1, index3, index4) / _finder2->diff(n-1, n-2));
			break; }
		default:
			throw std::runtime_error("Incorrect interpolation method.");
	}		 		
}
/* fills up order "pos" of the _slopes_d2 4Dmatrix of differentials with respect
to the third independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

void math::f_table4V::fill_up_slopes_aux1_hermite(math::vec4& slopes_d1,
												   unsigned short index4,
												   unsigned short index3,
												   unsigned short index2) {
	switch (_interp_mode) {
		case math::logic::hermite_first: {
			// compute temporary slopes as in Lagrange first order for each interval 
			vec1 temp(_points1->size1()-1);
			for (unsigned short k = 0; k != _points1->size1()-1; ++k) {
				temp[k] = (_VVValues->get(k+1, index2, index3, index4) - _VVValues->get(k, index2, index3, index4)) / _finder1->diff(k+1,k);
			}
			// compute slopes by average of previous slopes
			slopes_d1.set(0, index2, index3, index4, temp[0]);
			for (unsigned short k = 1; k != _points1->size1()-1; ++k) {
				slopes_d1.set(k, index2, index3, index4, 0.5 * (temp[k-1] + temp[k]));
			}
			slopes_d1.set(_points1->size1()-1, index2, index3, index4, temp.back());
			break; }
		case math::logic::hermite_second: {
			// compute slopes by results of Lagrange 2nd order interpolation
			unsigned short n = _points1->size1();
			slopes_d1.set(0, index2, index3, index4,
				 _VVValues->get(1, index2, index3, index4) / _finder1->diff(2, 0, 1, 0, 2, 1) -
				 _VVValues->get(2, index2, index3, index4) / _finder1->diff(1, 0, 2, 0, 2, 1) -
				 _VVValues->get(0, index2, index3, index4) / _finder1->diff(2, 0) -   
				 _VVValues->get(0, index2, index3, index4) / _finder1->diff(1, 0));
			for (unsigned short i = 1; i < n-1; ++i) {
				slopes_d1.set(i, index2, index3, index4,
				_VVValues->get(i+1, index2, index3, index4) / _finder1->diff(i,   i-1, i+1, i-1, i+1, i) -   
				_VVValues->get(i-1, index2, index3, index4) / _finder1->diff(i+1, i,   i,   i-1, i+1, i-1) +  			   
				_VVValues->get(i, index2, index3, index4)   / _finder1->diff(i,   i-1) -   
				_VVValues->get(i, index2, index3, index4)   / _finder1->diff(i+1, i));
			}
			slopes_d1.set(n-1, index2, index3, index4,
				_VVValues->get(n-3, index2, index3, index4) / _finder1->diff(n-1, n-2, n-2, n-3, n-1, n-3) -
				_VVValues->get(n-2, index2, index3, index4) / _finder1->diff(n-1, n-3, n-2, n-3, n-1, n-2) +
				_VVValues->get(n-1, index2, index3, index4) / _finder1->diff(n-1, n-3) +
				_VVValues->get(n-1, index2, index3, index4) / _finder1->diff(n-1, n-2));
			break; }
		default:
			throw std::runtime_error("Incorrect interpolation method.");
		}	
}
/* fills up order "pos" of the _slopes_d1 4Dmatrix of differentials with respect
to the fourth independent magnitude with a matrix of size "siz2" by "siz3" by "siz4"*/

void math::f_table4V::activate_checker4() {
	delete _checker4;
	_checker4 = new math::range_checker_active();
}
void math::f_table4V::deactivate_checker4() {
	delete _checker4;
	_checker4 = new math::range_checker_inactive();
}
/* Activates or deactivates the out of range verification for _points4, which
is inactive by default */

void math::f_table4V::activate_checker3() {
	delete _checker3;
	_checker3 = new math::range_checker_active();
}
void math::f_table4V::deactivate_checker3() {
	delete _checker3;
	_checker3 = new math::range_checker_inactive();
}
/* Activates or deactivates the out of range verification for _points3, which
is active by default */

void math::f_table4V::activate_checker2() {
	delete _checker2;
	_checker2 = new math::range_checker_active();
}
void math::f_table4V::deactivate_checker2() {
	delete _checker2;
	_checker2 = new math::range_checker_inactive();
}
/* Activates or deactivates the out of range verification for _points2, which
is active by default */

void math::f_table4V::activate_checker1() {
	delete _checker1;
	_checker1 = new math::range_checker_active();
}
void math::f_table4V::deactivate_checker1() {
	delete _checker1;
	_checker1 = new math::range_checker_inactive();
}
/* Activates or deactivates the out of range verification for _points1, which
is active by default */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABLE4V_DIFF_PREC
// =======================
// =======================
	
math::table4V_diff_prec::table4V_diff_prec(math::f_table4V& pred)
: _pred(&pred) {}
/* constructor based on equispaced four dimensional table */

math::table4V_diff_prec::~table4V_diff_prec() {}
/* destructor */

void math::table4V_diff_prec::compute_diff(double& result,
												const int& pos4,
												const int& pos3,
												const int& pos2,
												const int& pos1,
												const double& input4,
												const double& input3,
												const double& input2,
												const double& input1,
												const double& input4_dt,
												const double& input3_dt,
												const double& input2_dt,
												const double& input1_dt) const {
	math::ratio_linear* ratio4 = static_cast<math::ratio_linear*>(_pred->compute_ratio4(input4, pos4));
	math::ratio_linear* ratio3 = static_cast<math::ratio_linear*>(_pred->compute_ratio3(input3, pos3));
	math::ratio_linear* ratio2 = static_cast<math::ratio_linear*>(_pred->compute_ratio2(input2, pos2));
	math::ratio_linear* ratio1 = static_cast<math::ratio_linear*>(_pred->compute_ratio1(input1, pos1));
    // Differential with respect to 1st independent magnitude
    double temp, temp11, temp12, temp21, temp22;
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11,
                _pred->_slopes_d4->get(pos1,   pos2, pos3, pos4),
                _pred->_slopes_d4->get(pos1+1, pos2, pos3, pos4),  ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp12,
                _pred->_slopes_d4->get(pos1,   pos2+1, pos3, pos4),
                _pred->_slopes_d4->get(pos1+1, pos2+1, pos3, pos4), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21,
                _pred->_slopes_d4->get(pos1,   pos2, pos3+1, pos4),
                _pred->_slopes_d4->get(pos1+1, pos2, pos3+1, pos4),  ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp22,
                _pred->_slopes_d4->get(pos1,   pos2+1, pos3+1, pos4),
                _pred->_slopes_d4->get(pos1+1, pos2+1, pos3+1, pos4), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11, temp11, temp12, ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21, temp21, temp22, ratio2->_ratio12);
    math::interp::basic_interp_lagrange_first(temp, temp11, temp21, ratio3->_ratio12);
    result = temp * input4_dt;
    // Differential with respect to 2nd independent magnitude
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11,
                _pred->_slopes_d3->get(pos1,   pos2, pos3, pos4),
                _pred->_slopes_d3->get(pos1+1, pos2, pos3, pos4),  ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp12,
                _pred->_slopes_d3->get(pos1,   pos2+1, pos3, pos4),
                _pred->_slopes_d3->get(pos1+1, pos2+1, pos3, pos4), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21,
                _pred->_slopes_d3->get(pos1,   pos2, pos3, pos4+1),
                _pred->_slopes_d3->get(pos1+1, pos2, pos3, pos4+1), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp22,
                _pred->_slopes_d3->get(pos1,   pos2+1, pos3, pos4+1),
                _pred->_slopes_d3->get(pos1+1, pos2+1, pos3, pos4+1), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11, temp11, temp12, ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21, temp21, temp22, ratio2->_ratio12);
    math::interp::basic_interp_lagrange_first(temp, temp11, temp21, ratio4->_ratio12);
    result = result + temp * input3_dt;
    // Differential with respect to 3rd independent magnitude
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11,
                _pred->_slopes_d2->get(pos1,   pos2, pos3, pos4),
                _pred->_slopes_d2->get(pos1+1, pos2, pos3, pos4),  ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp12,
                _pred->_slopes_d2->get(pos1,   pos2, pos3+1, pos4),
                _pred->_slopes_d2->get(pos1+1, pos2, pos3+1, pos4),	ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21,
                _pred->_slopes_d2->get(pos1,   pos2, pos3, pos4+1),
                _pred->_slopes_d2->get(pos1+1, pos2, pos3, pos4+1), ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp22,
                _pred->_slopes_d2->get(pos1,   pos2, pos3+1, pos4+1),
                _pred->_slopes_d2->get(pos1+1, pos2, pos3+1, pos4+1),  ratio1->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11, temp11, temp12, ratio3->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21, temp21, temp22, ratio3->_ratio12);
    math::interp::basic_interp_lagrange_first(temp, temp11, temp21, ratio4->_ratio12);
    result = result + temp * input2_dt;
    // Differential with respect to 4th independent magnitude
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11,
                _pred->_slopes_d1->get(pos1, pos2,   pos3, pos4),
                _pred->_slopes_d1->get(pos1, pos2+1, pos3, pos4), ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp12,
                _pred->_slopes_d1->get(pos1, pos2,   pos3+1, pos4),
                _pred->_slopes_d1->get(pos1, pos2+1, pos3+1, pos4), ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21,
                _pred->_slopes_d1->get(pos1, pos2,   pos3, pos4+1),
                _pred->_slopes_d1->get(pos1, pos2+1, pos3, pos4+1), ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp22,
                _pred->_slopes_d1->get(pos1, pos2,   pos3+1, pos4+1),
                _pred->_slopes_d1->get(pos1, pos2+1, pos3+1, pos4+1), ratio2->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp11, temp11, temp12, ratio3->_ratio12);
    math::interp_lagrange_first::basic_interp_lagrange_first(temp21, temp21, temp22, ratio3->_ratio12);
    math::interp::basic_interp_lagrange_first(temp, temp11, temp21, ratio4->_ratio12);
    result = result + temp * input1_dt;
	math::ratio_mgr::to_pool_linear(ratio4);
	math::ratio_mgr::to_pool_linear(ratio3);
	math::ratio_mgr::to_pool_linear(ratio2);
	math::ratio_mgr::to_pool_linear(ratio1);
}
/* fill up result differential based on positions, inputs, and its partial
differentials with time */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
												
// CLASS TABLE4V_DIFF_REAL
// =======================
// =======================

math::table4V_diff_real::table4V_diff_real(math::f_table4V& pred)
: _pred(&pred) {}
/* constructor based on equispaced four dimensional table */

math::table4V_diff_real::~table4V_diff_real() {}
/* destructor */

void math::table4V_diff_real::compute_diff(double& result,
												const int& pos4,
												const int& pos3,
												const int& pos2,
												const int& pos1,
												const double& input4,
												const double& input3,
												const double& input2,
												const double& input1,
												const double& input4_dt,
												const double& input3_dt,
												const double& input2_dt,
												const double& input1_dt) const {
	math::ratio* ratio4 = _pred->compute_ratio4(input4, pos4);
	math::ratio* ratio3 = _pred->compute_ratio3(input3, pos3);
	math::ratio* ratio2 = _pred->compute_ratio2(input2, pos2);
	math::ratio* ratio1 = _pred->compute_ratio1(input1, pos1);
	math::ratio* ratioX4 = _pred->_interp->copy_ratio(ratio4);
	math::ratio* ratioX3 = _pred->_interp->copy_ratio(ratio3);
	math::ratio* ratioX2 = _pred->_interp->copy_ratio(ratio2);
	math::ratio* ratioX1 = _pred->_interp->copy_ratio(ratio1);
	double Pmag1, PmagX;
    // compute result at point
    _pred->_interp->interp4(Pmag1, *_pred->_points4, *_pred->_points3, *_pred->_points2, *_pred->_points1, *_pred->_VVValues, pos4, pos3, pos2, pos1, *ratio4, *ratio3, *ratio2, *ratio1, *_pred->_herm);
    // compute result a short interval after point in direction 1
    double difX = _pred->_finder4->compute_Dratio(*ratioX4, input4, pos4);
    _pred->_interp->interp4(PmagX, *_pred->_points4, *_pred->_points3, *_pred->_points2, *_pred->_points1, *_pred->_VVValues, pos4, pos3, pos2, pos1, *ratioX4, *ratio3, *ratio2, *ratio1, *_pred->_herm);
    result = (PmagX - Pmag1) * _pred->_finder4->compute_final_diff(input4_dt, difX);
    // compute result a short interval after point in direction 2
    difX = _pred->_finder3->compute_Dratio(*ratioX3, input3, pos3);
    _pred->_interp->interp4(PmagX, *_pred->_points4, *_pred->_points3, *_pred->_points2, *_pred->_points1, *_pred->_VVValues, pos4, pos3, pos2, pos1, *ratio4, *ratioX3, *ratio2, *ratio1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * _pred->_finder3->compute_final_diff(input3_dt, difX);
    // compute result a short interval after point in direction 3
    difX = _pred->_finder2->compute_Dratio(*ratioX2, input2, pos2);
    _pred->_interp->interp4(PmagX, *_pred->_points4, *_pred->_points3, *_pred->_points2, *_pred->_points1, *_pred->_VVValues, pos4, pos3, pos2, pos1, *ratio4, *ratio3, *ratioX2, *ratio1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * _pred->_finder2->compute_final_diff(input2_dt, difX);
    // compute result a short interval after point in direction 4
    difX = _pred->_finder1->compute_Dratio(*ratioX1, input1, pos1);
    _pred->_interp->interp4(PmagX, *_pred->_points4, *_pred->_points3, *_pred->_points2, *_pred->_points1, *_pred->_VVValues, pos4, pos3, pos2, pos1, *ratio4, *ratio3, *ratio2, *ratioX1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * _pred->_finder1->compute_final_diff(input1_dt, difX);
	// release memory
	_pred->_interp->to_pool(ratio4);
	_pred->_interp->to_pool(ratio3);
	_pred->_interp->to_pool(ratio2);
	_pred->_interp->to_pool(ratio1);		
	_pred->_interp->to_pool(ratioX4);
	_pred->_interp->to_pool(ratioX3);
	_pred->_interp->to_pool(ratioX2);
	_pred->_interp->to_pool(ratioX1);
}
/* fill up result differential based on positions, inputs, and its partial
differentials with time */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////












