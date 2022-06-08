#include "f_tabular2V.h"
#include "../pred1v/f_table1V.h"

// CLASS F_TABULAR2V
// =================
// =================

const math::logic::PRED_NAME math::f_tabular2V::_name = math::logic::f_tabular2V;
/* predicate name */

const std::string math::f_tabular2V::_st_name = "f_tabular2V";
/* predicate name string */

math::f_tabular2V::f_tabular2V(vec1* points2,
							  std::vector<math::f_table1V*>* tables,
							  math::logic::INTERP_MODE interp_mode)
: _points2(points2), _tables(tables),
_del_flag(true), _interp_mode(interp_mode),
_interp(math::interp::get_interp(interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_equi2(true), _functor_diff(0), _herm(0), _finder2(0),
_checker2(new math::range_checker_inactive()) {
	initialize();
}
/* constructor based on pointer to size n vector points2 and pointer to 
size n vector of pointers to unidimensional tables - deleted by destructor. */
math::f_tabular2V::f_tabular2V(vec1& points2,
							  std::vector<math::f_table1V*>& tables,
							  math::logic::INTERP_MODE interp_mode)
: _points2(&points2), _tables(&tables),
_del_flag(false), _interp_mode(interp_mode),
_interp(math::interp::get_interp(interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_equi2(true), _functor_diff(0),
_herm(0), _finder2(0), 
_checker2(new math::range_checker_inactive()) {
	initialize();
}
/* constructor based on reference to size n vector points2 and reference to
size n vector of pointers to unidimensional tables - not deleted by destructor. */ 

void math::f_tabular2V::initialize() {
	// check size compatibility between _points2 and _tables
	int size2 = _points2->size1();
	if (size2 != _tables->size()) {
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// points vectors have at least minimum number of points
	int n = _interp->get_min_points();
	if (size2 < n) {
        std::string st_name = "a";
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// difference between first two points of _points2
	_points2_diff = _points2->get(1) - _points2->get(0);
	if (_points2_diff <= 0.) { // they need to be in increasing order
		std::string st_name = "a";
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// difference between consecutive _points2
	double temp = 0.;
	for (unsigned short i = 1; i != size2 - 1; ++i) {
		temp = _points2->get(i+1) - _points2->get(i);
		if (temp <= 0.) { // they need to be in increasing order
			std::string st_name = "a";
			destroy();
            throw std::runtime_error("Incorrect size.");
		}		
		if (_equi2 == true) { // so far equispaced
			if (fabs((temp-_points2_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi2 = false;
				_points2_diff = 0.; 
			}	
		}
		// else: so far not equispaced and will remain that way
	}
	if (_equi2 == true) {_finder2 = new math::pos_finder_equispaced(*_points2, *_interp, _points2_diff);}
	else                {_finder2 = new math::pos_finder_binary(*_points2, *_interp);}

	switch (_interp_mode) {
		case math::logic::lagrange_first:
		case math::logic::lagrange_second:
		case math::logic::lagrange_third:
		case math::logic::biparabolic:
		    _herm = new hermite2v(); // dummy
			_functor_diff = new math::tabular2V_diff_real(*this);
			break;
		default: // none other cases allowed
			// it is impossible to precompute anything as size of tables varies
			destroy();
			throw std::runtime_error("Incorrect interpolation method.");
			break;
	}
}
/**< initialization for constructors */

math::f_tabular2V::f_tabular2V(const f_tabular2V& other)
: _del_flag(other._del_flag),
_points2_diff(other._points2_diff), _equi2(other._equi2),
_interp_mode(other._interp_mode), _herm(0), 
_interp(math::interp::get_interp(_interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_checker2(other._checker2->clone()),
_finder2(other._finder2->clone()) {
	
	if (other._herm != 0) _herm = new hermite2v(*other._herm);
	
	if (_del_flag == true) {
		_points2 = other._points2->clone();
		_tables = new std::vector<math::f_table1V*>(other._tables->size(), 0);
		for (unsigned short i = 0; i != _tables->size(); ++i) 
			(*_tables)[i] = (*other._tables)[i]->clone();
	}
	else {
		_points2 = other._points2;
		_tables  = other._tables;
	}

	switch (_interp_mode) {
		case math::logic::lagrange_first:
		case math::logic::lagrange_second:
		case math::logic::lagrange_third:
		case math::logic::biparabolic:
			_functor_diff = new math::tabular2V_diff_real(*this);
			break;
		default:
            throw std::runtime_error("Incorrect interpolation method.");
			break;
	}
}
/* copy constructor */

void math::f_tabular2V::destroy() {
	if (_del_flag == true) {
		delete _points2;
		for (unsigned short i = 0; i != _tables->size(); ++i) delete (*_tables)[i];
		delete _tables;
	}
	delete _interp;
	delete _interp_diff;
	delete _herm;
	delete _checker2;
	delete _finder2;
	delete _functor_diff;
}
/* destructor */

math::f_tabular2V* math::f_tabular2V::clone() const {
	return new f_tabular2V(*this);
}
/* cloner */

bool math::f_tabular2V::operator==(const pred2v& op2) const {
	return (op2.get_name() == math::logic::f_tabular2V) ?
		(*this == static_cast<const f_tabular2V&>(op2)) : false;	
}
bool math::f_tabular2V::operator==(const f_tabular2V& op2) const {
	if ((*_points2 != *op2._points2) || (_interp_mode != op2._interp_mode)) {
		return false;
	}
	for (unsigned short i = 0; i != _tables->size(); ++i) {
		if (*(*_tables)[i] != *(*op2._tables)[i]) {
			return false;
		}
	}
	return true;
}
/* overloaded operator == (equal) */

int math::f_tabular2V::compute_pos2(const double& input2) const {
	return _interp->find_index(_finder2->search(*_points2, input2, _points2_diff),
								_points2->size1());
}
/* Returns first position within _points2 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

math::ratio* math::f_tabular2V::compute_ratio2(const double& input2,
												    const int& pos2) const {
	return _interp->compute_ratio(input2, *_points2, pos2);
}
/* Returns ratio of input2 with respect to the two points2 identified by pos2 and pos2+1. */

void math::f_tabular2V::compute_value(double& result,
									 const int& pos2, 
									 const math::ratio& ratio2,
									 const double& input1) const {
	_interp->interp2(result, *_points2, *_tables, pos2, ratio2, input1, *_herm);
}
/* Fills up result magnitude by interpolating based on the second dimension
position and ratio plus the other input. */

void math::f_tabular2V::compute_diff(double& result,
								    const int& pos2,
								    const double& input2,
								    const double& input1,
								    const double& input2_dt,
								    const double& input1_dt) const {
	_functor_diff->compute_diff(result, pos2,
						input2, input1, input2_dt, input1_dt);	
}
/* Fills up result differential based on positions and ratios */

double math::f_tabular2V::value(const double& input2,
                               const double& input1) const {
	// verify inputs are within range
	_checker2->check_range(*_points2, input2);
	// pos2 provides the upper vector position
	int pos2 = compute_pos2(input2);
	math::ratio* ratio2 = compute_ratio2(input2, pos2);
    double res;
	_interp->interp2(res, *_points2, *_tables, pos2, *ratio2, input1, *_herm);
	_interp->to_pool(ratio2);
    return res;
}
/* evaluates the function at the reference magnitudes input2, and
input1, and writes the result at the reference magnitude result. Only the
magnitude value is inserted into result, the units are assummed to be OK. */

double math::f_tabular2V::d_dt(const double& input2,
                         const double& input1,
                         const double& input2_dt,
                         const double& input1_dt) const {
	////////////////////////////////////////////////////////////////////
	// The computation of the two positions, which are quite expensive,
	// in theory have already been done before in the value method, and 
	// should not be repeated here. However, there is no way to avoid it.
	////////////////////////////////////////////////////////////////////
	// verify inputs are within range
	_checker2->check_range(*_points2, input2);
	// pos2 provides the upper vector position
	int pos2 = compute_pos2(input2);
    double res;
	_functor_diff->compute_diff(res, pos2, input2, input1, input2_dt, input1_dt);
    return res;
}
/* evaluates the function differential with time at the reference
magnitudes input2, and input1 and their differentials with time
input2_dt, and input1_dt, and writes the result at the reference
magnitude result. */

void math::f_tabular2V::activate_checker2() {
	delete _checker2;
	_checker2 = new math::range_checker_active();
}
void math::f_tabular2V::deactivate_checker2() {
	delete _checker2;
	_checker2 = new math::range_checker_inactive();
}
/* Activates or deactives the out of range verification for _points2, which
is inactive by default */


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR2V_DIFF_REAL
// =========================
// =========================

math::tabular2V_diff_real::tabular2V_diff_real(math::f_tabular2V& pred)
: _pred(&pred) {
}
/* constructor based on two dimensional table */

math::tabular2V_diff_real::~tabular2V_diff_real() {
}
/* destructor */

void math::tabular2V_diff_real::compute_diff(double& result,
											    const int& pos2,
											    const double& input2,
											    const double& input1,
											    const double& input2_dt,
											    const double& input1_dt) const {

	math::ratio* ratio2 = _pred->compute_ratio2(input2, pos2);
	math::ratio* ratioX2 = _pred->_interp->copy_ratio(ratio2);
	double Pmag1, PmagX;
	double Pinput1 = input1;
    // compute result at point
    _pred->_interp->interp2(Pmag1, *_pred->_points2, *_pred->_tables, pos2, *ratio2, input1, *_pred->_herm);

    // compute result a short interval after point in direction 1
    double diff2 = _pred->_finder2->compute_Dratio(*ratioX2, input2, pos2);
    _pred->_interp->interp2(PmagX, *_pred->_points2, *_pred->_tables, pos2, *ratioX2, input1, *_pred->_herm);
    result = (PmagX - Pmag1) * _pred->_finder2->compute_final_diff(input2_dt, diff2);

    // find unidimensional table inmediately below the input2 (lagrange first interpolation for that)
    int Xpos2 = _pred->_interp_diff->find_index(
                _pred->get_finder2().search(*_pred->_points2, input2, _pred->_points2_diff),
                _pred->_points2->size1());
    const math::f_table1V& Otable1V = *(*_pred->_tables)[Xpos2];

    // compute results a short interval after point in direction 2
    // to obtain that input value, use the unidimensional table for position Xpos2
    // find position inmediately below in that table (lagrange first interpolation for that only)
    const math::vec1& Opoints1 = Otable1V.get_points1();
    int Xpos1 = _pred->_interp_diff->find_index(
                Otable1V.get_finder1().search(Opoints1, input1, Otable1V.get_points1_diff()), Opoints1.size1());
    double diff1 = 1e-4 * (Opoints1[Xpos1+1] - Opoints1[Xpos1]);
    Pinput1 = input1 + diff1;
    _pred->_interp->interp2(PmagX, *_pred->_points2, *_pred->_tables, pos2, *ratio2, Pinput1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * input1_dt / diff1;
	// release memory
	_pred->_interp->to_pool(ratio2);
	_pred->_interp->to_pool(ratioX2);
}
/* fill up result differential based on positions, inputs, and its partial
differentials with time */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

