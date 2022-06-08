#ifndef MATH_F_TABULAR3V_H
#define MATH_F_TABULAR3V_H

#include "../../math.h"
#include "pred3v.h"

/*
Derivate class of pred3v containing f_tabular3V objects. Refer to pred3v
for more info. Contains a vector of two dimensional tables of class f_table2V.

Extrapolation when outside limits. However, tables are protected for 
all dimensions except first (generally altitude), so an exception is 
launched if extrapolation is required in the 2nd or 3rd dimensions.

NOTE: Watch out when input magnitude is either a longitude or a bearing,
as a discontinuity in the results appears below the first input (inputs below
the first are considered as being after the last).

LEFT:
- Binary search could be accelerated if input from previous execution
  were added to the input. It would need to be stored however.
- Improve accuracy of differential computation by doing before and after point,
  instead of only after.
- For Hermite interpolation, there exists an explicit expression for the
  differential computation.

SLOW:
- Repetition of binary searches.
- The linear_precompute interpolation method requires precomputing slopes at
  construction time (slower) but the computation of differentials is faster. The
  spline interpolation method is not allowed for three dimensions.
  Hermite also requires precomputation.
*/

namespace math {

	class f_table2V;
	class tabular3V_diff;

// CLASS F_TABULAR3V
// =================
// =================

class MATH_API f_tabular3V : public pred3v {
private:
	friend class tabular3V_diff_real;
	/**< classes that have access to private members */

	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */
	math::vec1* _points3;
	/**< pointer to object of class vec1 (ordered magnitude vector of size
	n defining the table 1st inputs). */
	std::vector<math::f_table2V*>* _tables;
	/**< pointer to vector of size n (of same size as _points3) of pointers
	to bidimensional tables f_table2V, each one corresponding to a _points3 member */
	bool _del_flag;
	/**< flag that indicates if the _points3 and _tables vectors
	shall be deleted by the destructor (true) or not (false) */

	math::logic::INTERP_MODE _interp_mode;
	/**< interpolation method name */
	math::interp* _interp;
	/**< pointer to interpolation method */
	math::interp* _interp_diff;
	/**< pointer to interpolation method only employed to compute differentials
	dealing with points1 and points2. Use only if you know what you are doing. */

	bool _equi3;
	/**< true if _points3 vector equispaced, false otherwise */
	double _points3_diff;
	/**< difference between consecutive _points3, 0. if not equispaced */
	math::tabular3V_diff* _functor_diff;
	/**< pointer to functor that solves differentials */
	math::hermite3v* _herm;
	/**< pointer to Hermite coefficients. It is a useless dummy here but required for
	compatibility with other tables */
	math::range_checker* _checker3;
	/**< object controlling out of range for _points3, inactive by default */
	math::pos_finder* _finder3;
	/**< object determining position of input with respect to _points3 */

	f_tabular3V();
	/**< empty constructor not implemented */
	f_tabular3V& operator=(const f_tabular3V& op2);
	/**< overloaded operator = (assignment) not implemented */

	void initialize();
	/**< initialization for constructors */
	void destroy();
	/**< destructor */
public:
	f_tabular3V(math::vec1* points3,
				std::vector<math::f_table2V*>* tables,
				math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on pointer to size l vector points3 and pointer to 
	size l vector of pointers to bidimensional tables - deleted	by destructor. */ 
	f_tabular3V(math::vec1& points3,
				std::vector<math::f_table2V*>& tables,
			    math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on reference to size l vector points3 and reference to
	size l vector of pointers to bidimensional tables - not deleted by destructor. */ 
	f_tabular3V(const f_tabular3V&);
	/**< copy constructor */
	virtual ~f_tabular3V() {destroy();}
	/**< destructor */
	virtual f_tabular3V* clone() const;
	/**< cloner */
	bool operator==(const pred3v& op2) const;
	bool operator==(const f_tabular3V& op2) const;
	/**< overloaded operator == (equal) */
	bool operator!=(const f_tabular3V& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

	const math::vec1& get_points3() const {return *_points3;}
	/**< get the _points3 object as constant reference to read */
	const std::vector<math::f_table2V*>& get_tables() const {return *_tables;}
	/**< get the _tables vector as constant reference to read */
	const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred3v above */
	const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred3v above */
	const math::logic::INTERP_MODE& get_interp_mode() const {return _interp_mode;}
	/**< get the interpolation mode to read */
	const math::interp& get_interp() const {return *_interp;}
	/**< get constant reference to interpolation method */
	bool get_equi3() const {return _equi3;}
	/**< returns true if _points3 vector equispaced, false otherwise */
	const math::pos_finder& get_finder3() const {return *_finder3;}
	/**< returns constant reference to finder object */
	const double& get_points3_diff() const {return _points3_diff;}
	/**< returns difference between consecutive _points3, 0. if not equispaced */

	void activate_checker3();
	void deactivate_checker3();
	/**< Activates or deactives the out of range verification for _points3, which
	is inactive by default */

	int compute_pos3(const double& input3) const;
	/**< Returns first position within _points3 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	math::ratio* compute_ratio3(const double& input3,
									  const int& pos3) const;
	/**< Returns ratio of input3 with respect to the two points3 identified by pos3
	and pos3+1. */
	void compute_value(double& result,
					   const int& pos3, 
					   const math::ratio& ratio3,
					   const double& input2,
					   const double& input1) const;
	/**< Fills up result magnitude by interpolating based on the third dimension
	position and ratio plus the other two inputs. */
	void compute_diff(double& result,
					  const int& pos3, 
					  const double& input3,
					  const double& input2,
					  const double& input1,
					  const double& input3_dt,
					  const double& input2_dt,
					  const double& input1_dt) const;
	/**< Fills up result differential based on positions and ratios */

    double value(const double& input3,
			   const double& input2,
			   const double& input1) const;
	/**< see virtual function of class predicate3v above. */
    double d_dt(const double& input3, const double& input2, const double& input1,
			  const double& input3_dt, const double& input2_dt, const double& input1_dt) const;
	/**< see virtual function of class pred3v above. */

}; // closes class f_tabular3V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR3V_DIFF
// ====================
// ====================

class MATH_API tabular3V_diff {
public:
	virtual void compute_diff(double& result,
							  const int& pos3,
							  const double& input3,
							  const double& input2,
							  const double& input1,
							  const double& input3_dt,
							  const double& input2_dt,
							  const double& input1_dt) const = 0;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class tabular3V_diff

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABULAR3V_DIFF_REAL
// =========================
// =========================

class MATH_API tabular3V_diff_real: public tabular3V_diff {
private:
	math::f_tabular3V* _pred;
	/**< weak pointer to f_tabular3V class */
	tabular3V_diff_real();
	/**< empty constructor not implemented */
	tabular3V_diff_real(const tabular3V_diff_real&);
	/**< copy constructor not implemented */
	tabular3V_diff_real& operator=(const tabular3V_diff_real& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	tabular3V_diff_real(math::f_tabular3V& pred);
	/**< constructor based on three dimensional table */
	~tabular3V_diff_real();
	/**< destructor */
	void compute_diff(double& result,
					  const int& pos3,
					  const double& input3,
					  const double& input2,
					  const double& input1,
					  const double& input3_dt,
					  const double& input2_dt,
					  const double& input1_dt) const;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class tabular3V_diff_real

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

} // closes namespace math

#endif





