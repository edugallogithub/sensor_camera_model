#include "f_tabular3V.h"
#include "../pred2v/f_table2V.h"

// CLASS F_TABULAR3V
// =================
// =================

const math::logic::PRED_NAME math::f_tabular3V::_name = math::logic::f_tabular3V;
/* predicate name */

const std::string math::f_tabular3V::_st_name = "f_tabular3V";
/* predicate name string */

math::f_tabular3V::f_tabular3V(vec1* points3,
							  std::vector<math::f_table2V*>* tables,
							  math::logic::INTERP_MODE interp_mode)
: _points3(points3), _tables(tables),
_del_flag(true), _interp_mode(interp_mode),
_interp(math::interp::get_interp(interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_equi3(true), _functor_diff(0), _herm(0), _finder3(0),
_checker3(new math::range_checker_inactive()) {
	initialize();
}
/* constructor based on pointer to size l vector points3 and pointer to 
size l vector of pointers to bidimensional tables - deleted	by destructor. */

math::f_tabular3V::f_tabular3V(vec1& points3,
							  std::vector<math::f_table2V*>& tables,
							  math::logic::INTERP_MODE interp_mode)
: _points3(&points3), _tables(&tables),
_del_flag(false), _interp_mode(interp_mode),
_interp(math::interp::get_interp(interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_equi3(true), _functor_diff(0),
_herm(0), _finder3(0), 
_checker3(new math::range_checker_inactive()) {
	initialize();
}
/* constructor based on reference to size l vector points3 and reference to
size l vector of pointers to bidimensional tables - not deleted by destructor. */ 

void math::f_tabular3V::initialize() {
	// check size compatibility between _points3 and _tables
	int size3 = _points3->size1();
	if (size3 != _tables->size()) {
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// points vectors have at least minimum number of points
	int n = _interp->get_min_points();
	if (size3 < n) {
        std::string st_name = "a";
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// difference between first two points of _points3
	_points3_diff = _points3->get(1) - _points3->get(0);
	if (_points3_diff <= 0.) { // they need to be in increasing order
		std::string st_name = "a";
		destroy();
        throw std::runtime_error("Incorrect size.");
	}
	// difference between consecutive _points3
	double temp = 0.;
	for (unsigned short i = 1; i != size3 - 1; ++i) {
		temp = _points3->get(i+1) - _points3->get(i);
		if (temp <= 0.) { // they need to be in increasing order
			std::string st_name = "a";
			destroy();
            throw std::runtime_error("Incorrect size.");
		}
		if (_equi3 == true) { // so far equispaced
			if (fabs((temp-_points3_diff)/temp) > constant::DIFF()) {
				// not equiespaced
				_equi3 = false;
				_points3_diff = 0.; 
			}	
		}
		// else: so far not equispaced and will remain that way
	}
	if (_equi3 == true) {_finder3 = new math::pos_finder_equispaced(*_points3, *_interp, _points3_diff);}
	else                {_finder3 = new math::pos_finder_binary(*_points3, *_interp);}

	switch (_interp_mode) {
		case math::logic::lagrange_first:
		case math::logic::lagrange_second:
		case math::logic::lagrange_third:
		case math::logic::biparabolic:
			_herm = new hermite3v(); // dummy
			_functor_diff = new math::tabular3V_diff_real(*this);
			break;
		default: // none other cases allowed
			// it is impossible to precompute anything as size of tables varies
			destroy();
            throw std::runtime_error("Incorrect interpolation method.");
			break;
	}
}
/**< initialization for constructors */

math::f_tabular3V::f_tabular3V(const f_tabular3V& other)
: _del_flag(other._del_flag),
_points3_diff(other._points3_diff), _equi3(other._equi3),
_interp_mode(other._interp_mode), _herm(0), 
_interp(math::interp::get_interp(_interp_mode)),
_interp_diff(math::interp::get_interp(math::logic::lagrange_first)),
_checker3(other._checker3->clone()),
_finder3(other._finder3->clone()) {
	
	if (other._herm != 0) _herm = new hermite3v(*other._herm);

	if (_del_flag == true) {
		_points3 = other._points3->clone();
		_tables = new std::vector<math::f_table2V*>(other._tables->size(), 0);
		for (unsigned short i = 0; i != _tables->size(); ++i)
			(*_tables)[i] = (*other._tables)[i]->clone();
	}
	else {
		_points3 = other._points3;
		_tables  = other._tables;
	}

	switch (_interp_mode) {
		case math::logic::lagrange_first:
		case math::logic::lagrange_second:
		case math::logic::lagrange_third:
		case math::logic::biparabolic:
			_functor_diff = new math::tabular3V_diff_real(*this);
			break;
		default:
			throw std::runtime_error("Incorrect interpolation mode.");
    }
}
/* copy constructor */

void math::f_tabular3V::destroy() {
	if (_del_flag == true) {
		delete _points3;
		for (unsigned short i = 0; i != _tables->size(); ++i) delete (*_tables)[i];
		delete _tables;
	}
	delete _interp;
	delete _interp_diff;
	delete _herm;
	delete _checker3;
	delete _finder3;
	delete _functor_diff;
}
/* destructor */

math::f_tabular3V* math::f_tabular3V::clone() const {
	return new f_tabular3V(*this);
}
/* cloner */

bool math::f_tabular3V::operator==(const pred3v& op2) const {
	return (op2.get_name() == math::logic::f_tabular3V) ?
		(*this == static_cast<const f_tabular3V&>(op2)) : false;	
}
bool math::f_tabular3V::operator==(const f_tabular3V& op2) const {
	if ((*_points3 != *op2._points3) || (_interp_mode != op2._interp_mode)) {
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

int math::f_tabular3V::compute_pos3(const double& input3) const {
	return _interp->find_index(_finder3->search(*_points3, input3, _points3_diff),
								_points3->size1());
}
/* Returns first position within _points3 vector that shall be employed when
interpolating to obtain the result corresponding to the input magnitude. */

math::ratio* math::f_tabular3V::compute_ratio3(const double& input3,
												    const int& pos3) const {
	return _interp->compute_ratio(input3, *_points3, pos3);
}
/* Returns ratio of input3 with respect to the two points3 identified by pos3 and pos3+1. */

void math::f_tabular3V::compute_value(double& result,
                                     const int& pos3,
									 const math::ratio& ratio3,
									 const double& input2,
									 const double& input1) const {
	_interp->interp3(result, *_points3, *_tables, pos3, ratio3, input2, input1, *_herm);
}
/* Fills up result magnitude by interpolating based on the third dimension
position and ratio plus the other two inputs. */

void math::f_tabular3V::compute_diff(double& result,
									const int& pos3, 
									const double& input3,
									const double& input2,
									const double& input1,
									const double& input3_dt,
									const double& input2_dt,
									const double& input1_dt) const {
	_functor_diff->compute_diff(result, pos3, 
						input3, input2, input1, input3_dt, input2_dt, input1_dt);	
}
/* Fills up result differential based on positions and ratios */

double math::f_tabular3V::value(const double& input3,
						     const double& input2,
						     const double& input1) const {
	// verify inputs are within range
	_checker3->check_range(*_points3, input3);
	// pos3 provides the upper vector position
	int pos3 = compute_pos3(input3);
	math::ratio* ratio3 = compute_ratio3(input3, pos3);
    double res;
	this->_interp->interp3(res, *_points3, *_tables, pos3, *ratio3, input2, input1, *_herm);
	_interp->to_pool(ratio3);
    return res;
}
/* evaluates the function at the reference magnitudes input3, input2, and
input1, and writes the result at the reference magnitude result. Only the
magnitude value is inserted into result, the units are assummed to be OK. */

double math::f_tabular3V::d_dt(const double& input3, const double& input2, const double& input1,
                         const double& input3_dt, const double& input2_dt, const double& input1_dt) const {
	////////////////////////////////////////////////////////////////////
	// The computation of the three positions, which are quite expensive,
	// in theory have already been done before in the value method, and 
	// should not be repeated here. However, there is no way to avoid it.
	////////////////////////////////////////////////////////////////////
	// verify inputs are within range
	_checker3->check_range(*_points3, input3);
	// pos3 provides the upper vector position
	int pos3 = compute_pos3(input3);
    double res;
    _functor_diff->compute_diff(res, pos3,
						input3, input2, input1, input3_dt, input2_dt, input1_dt);
    return res;
}
/* evaluates the function differential with time at the reference
magnitudes input3, input2, and input1 and their differentials with time
input3_dt, input2_dt, and input1_dt, and writes the result at the reference
magnitude result. */

void math::f_tabular3V::activate_checker3() {
	delete _checker3;
	_checker3 = new math::range_checker_active();
}
void math::f_tabular3V::deactivate_checker3() {
	delete _checker3;
	_checker3 = new math::range_checker_inactive();
}
/* Activates or deactives the out of range verification for _points3, which
is inactive by default */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR3V_DIFF_REAL
// =========================
// =========================

math::tabular3V_diff_real::tabular3V_diff_real(math::f_tabular3V& pred)
: _pred(&pred) {
}
/* constructor based on three dimensional table */

math::tabular3V_diff_real::~tabular3V_diff_real() {
}
/* destructor */

void math::tabular3V_diff_real::compute_diff(double& result,
											    const int& pos3,
											    const double& input3,
											    const double& input2,
											    const double& input1,
											    const double& input3_dt,
											    const double& input2_dt,
											    const double& input1_dt) const {
	math::ratio* ratio3 = _pred->compute_ratio3(input3, pos3);
	math::ratio* ratioX3 = _pred->_interp->copy_ratio(ratio3);
	double Pmag1, PmagX;
	double Pinput2, Pinput1;
    // compute result at point
    _pred->_interp->interp3(Pmag1, *_pred->_points3, *_pred->_tables, pos3, *ratio3, input2, input1, *_pred->_herm);

    // compute result a short interval after point in direction 1
    double diff3 = _pred->_finder3->compute_Dratio(*ratioX3, input3, pos3);
    _pred->_interp->interp3(PmagX, *_pred->_points3, *_pred->_tables, pos3, *ratioX3, input2, input1, *_pred->_herm);
    result = (PmagX - Pmag1) * _pred->_finder3->compute_final_diff(input3_dt, diff3);

    // find bidimensional table inmediately below the input3 (lagrange first interpolation for that)
    int Xpos3 = _pred->_interp_diff->find_index(
                _pred->get_finder3().search(*_pred->_points3, input3, _pred->_points3_diff),
                _pred->_points3->size1());
    const math::f_table2V& Otable2V = *(*_pred->_tables)[Xpos3];

    // compute results a short interval after point in direction 2
    // to obtain that input value, use the bidimensional table for position Xpos3
    // find position inmediately below in that table (lagrange first interpolation for that only)
    const math::vec1& Opoints2 = Otable2V.get_points2();
    int Xpos2 = _pred->_interp_diff->find_index(
                Otable2V.get_finder2().search(Opoints2, input2, Otable2V.get_points2_diff()), Opoints2.size1());
    double diff2 = 1e-4 * (Opoints2[Xpos2+1] - Opoints2[Xpos2]);
    Pinput2 = input2 + diff2;
    _pred->_interp->interp3(PmagX, *_pred->_points3, *_pred->_tables, pos3, *ratio3, Pinput2, input1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * input2_dt / diff2;

    // compute results a short interval after point in direction 3
    // to obtain that input value, use the bidimensional table for position Xpos3
    // find position inmediately below in that table (lagrange first interpolation for that only)
    const math::vec1& Opoints1 = Otable2V.get_points1();
    int Xpos1 = _pred->_interp_diff->find_index(
                Otable2V.get_finder1().search(Opoints1, input1, Otable2V.get_points1_diff()), Opoints1.size1());
    double diff1 = 1e-4 * (Opoints1[Xpos1+1] - Opoints1[Xpos1]);
    Pinput1 = input1 + diff1;
    _pred->_interp->interp3(PmagX, *_pred->_points3, *_pred->_tables, pos3, *ratio3, input2, Pinput1, *_pred->_herm);
    result = result + (PmagX - Pmag1) * input1_dt / diff1;
	// release memory
	_pred->_interp->to_pool(ratio3);
	_pred->_interp->to_pool(ratioX3);
}
/* fill up result differential based on positions, inputs, and its partial
differentials with time */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////





