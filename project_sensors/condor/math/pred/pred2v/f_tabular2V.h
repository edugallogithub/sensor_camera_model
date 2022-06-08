
#ifndef MATH_F_TABULAR2V_H
#define TABLE_F_TABULAR2V_H

#include "../../math.h"
#include "pred2v.h"

/*
Derivate class of pred2v containing f_tabular2V objects. Refer to pred2v
for more info. Contains a vector of one dimensional tables of class f_table1V.

Extrapolation when outside limits. 

NOTE: Watch out when input magnitude is either a longitude or a bearing,
as a discontinuity in the results appears below the first input (inputs below
the first are considered as being after the last).

LEFT:
- Binary search could be accelerated if input from previous execution
  were added to the input. It would need to be stored however.
- Improve accuracy of differntial computation by doing before and after point,
  instead of only after.
- For Hermite interpolation, there exists an explicit expression for the
  differential computation.

SLOW:
- Repetition of binary searches.
- The linear_precompute interpolation method requires precomputing slopes at
  construction time (slower) but the computation of differentials is faster. The
  spline interpolation method requires computing the spline second differentials
  at construction time may times (slow), but also once at execution time (very slow). 
  Hermite also requires precomputation.

*/

namespace math {

	class f_table1V;
	class tabular2V_diff;

// CLASS F_TABULAR2V
// =================
// =================

class MATH_API f_tabular2V : public pred2v {
private:
	friend class tabular2V_diff_real;
	/**< classes that have access to private members */

	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */
	math::vec1* _points2;
	/**< pointer to object of class vec1 (ordered magnitude vector of size
	n defining the table 1st inputs). */
	std::vector<math::f_table1V*>* _tables;
	/**< pointer to vector of size n (of same size as _points2) of pointers
	to one dimensional tables f_table1V, each one corresponding to a _points2 member */
	bool _del_flag;
	/**< flag that indicates if the _points2 and _tables vectors
	shall be deleted by the destructor (true) or not (false) */

	math::logic::INTERP_MODE _interp_mode;
	/**< interpolation method name */
	math::interp* _interp;
	/**< pointer to interpolation method */
	math::interp* _interp_diff;
	/**< pointer to interpolation method only employed to compute differentials
	dealing with points2. Use only if you know what you are doing. */
	
	bool _equi2;
	/**< true if _points2 vector equispaced, false otherwise */
	double _points2_diff;
	/**< difference between consecutive _points2, 0. if not equispaced */

	math::tabular2V_diff* _functor_diff;
	/**< pointer to functor that solves differentials */
	math::hermite2v* _herm;
	/**< pointer to Hermite coefficients. It is a useless dummy here but required for
	compatibility with other tables */

	math::range_checker* _checker2;
	/**< object controlling out of range for _points2, inactive by default */
	math::pos_finder* _finder2;
	/**< object determining position of input with respect to _points2 */
	
	f_tabular2V();
	/**< empty constructor not implemented */
	f_tabular2V& operator=(const f_tabular2V& op2);
	/**< overloaded operator = (assignment) not implemented */

	void initialize();
	/**< initialization for constructors */
	void destroy();
	/**< destructor */
public:
	f_tabular2V(math::vec1* points2,
				std::vector<math::f_table1V*>* tables,
				math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on pointer to size n vector points2 and pointer to 
	size n vector of pointers to unidimensional tables - deleted by destructor. */ 
	f_tabular2V(math::vec1& points2,
				std::vector<math::f_table1V*>& tables,
			    math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on reference to size n vector points2 and reference to
	size n vector of pointers to unidimensional tables - not deleted by destructor. */ 

	f_tabular2V(const f_tabular2V&);
	/**< copy constructor */
	virtual ~f_tabular2V() {destroy();}
	/**< destructor */
	virtual f_tabular2V* clone() const;
	/**< cloner */
	bool operator==(const pred2v& op2) const;
	bool operator==(const f_tabular2V& op2) const;
	/**< overloaded operator == (equal) */
	bool operator!=(const f_tabular2V& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

	const math::vec1& get_points2() const {return *_points2;}
	/**< get the _points2 object as constant reference to read */
	const std::vector<math::f_table1V*>& get_tables() const {return *_tables;}
	/**< get the _tables vector as constant reference to read */
	const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred2v above */
	const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred2v above */
	const math::logic::INTERP_MODE& get_interp_mode() const {return _interp_mode;}
	/**< get the interpolation mode to read */
	const math::interp& get_interp() const {return *_interp;}
	/**< get constant reference to interpolation method */
	bool get_equi2() const {return _equi2;}
	/**< returns true if _points2 vector equispaced, false otherwise */
	const math::pos_finder& get_finder2() const {return *_finder2;}
	/**< returns constant reference to finder object */
	const double& get_points2_diff() const {return _points2_diff;}
	/**< returns difference between consecutive _points2, 0. if not equispaced */

	void activate_checker2();
	void deactivate_checker2();
	/**< Activates or deactives the out of range verification for _points2, which
	is inactive by default */

	int compute_pos2(const double& input2) const;
	/**< Returns first position within _points2 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	math::ratio* compute_ratio2(const double& input2,
									  const int& pos2) const;
	/**< Returns ratio of input2 with respect to the two points2 identified by pos2 and
	pos2+1. */
	void compute_value(double& result,
					   const int& pos2, 
					   const math::ratio& ratio2,
					   const double& input1) const;
	/**< Fills up result magnitude by interpolating based on the second dimension
	position and ratio plus the other input. */
	void compute_diff(double& result,
					  const int& pos2, 
					  const double& input2,
					  const double& input1,
					  const double& input2_dt,
					  const double& input1_dt) const;
	/**< Fills up result differential based on positions and ratios */

    double value(const double& input2,
                 const double& input1) const;
	/**< see virtual function of class predicate2v above. */
	double d_dt(const double& input2,
                const double& input1,
                const double& input2_dt,
                const double& input1_dt) const;
	/**< see virtual function of class pred2v above. */
}; // closes class f_tabular2V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR2V_DIFF
// ====================
// ====================

class MATH_API tabular2V_diff {
public:
	virtual void compute_diff(double& result,
							  const int& pos2,
							  const double& input2,
							  const double& input1,
							  const double& input2_dt,
							  const double& input1_dt) const = 0;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class tabular2V_diff

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR2V_DIFF_REAL
// =========================
// =========================

class MATH_API tabular2V_diff_real: public tabular2V_diff {
private:
	math::f_tabular2V* _pred;
	/**< weak pointer to f_tabular2V class */
	tabular2V_diff_real();
	/**< empty constructor not implemented */
	tabular2V_diff_real(const tabular2V_diff_real&);
	/**< copy constructor not implemented */
	tabular2V_diff_real& operator=(const tabular2V_diff_real& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	tabular2V_diff_real(math::f_tabular2V& pred);
	/**< constructor based on two dimensional table */
	~tabular2V_diff_real();
	/**< destructor */
	void compute_diff(double& result,
					  const int& pos2,
					  const double& input2,
					  const double& input1,
					  const double& input2_dt,
					  const double& input1_dt) const;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class tabular2V_diff_real

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

} // closes namespace math

#endif

