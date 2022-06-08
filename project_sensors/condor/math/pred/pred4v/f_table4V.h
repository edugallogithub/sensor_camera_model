#ifndef MATH_F_TABLE4V_H
#define MATH_F_TABLE4V_H

#include "../../math.h"
#include "pred4v.h"

/*
Derivate class of pred4v containing f_table4Veq objects. Refer to pred4v
for more info. Contains a four dimensional table.

Extrapolation when outside limits. However, tables are protected for 
all dimensions except first (generally altitude), so an exception is 
launched if extrapolation is required in the 2nd, 3rd, or 4th dimensions.

NOTE: Watch out when input magnitude is either a longitude or a bearing,
as a discontinuity in the results appears below the first input (inputs below
the first are considered as being after the last).

LEFT:
- Improve accuracy of differential computation by doing before and after point,
  instead of only after.
- For Hermite interpolation, there exists an explicit expression for the
  differential computation.

SLOW:
- The linear_precompute interpolation method requires precomputing slopes at
  construction time (slower) but the computation of differentials is faster. The
  spline interpolation method requires computing the spline second differentials
  at construction time may times (slow), but also once at execution time (very slow). 
  Hermite also requires precomputation.

*/

namespace math {

	class table4V_diff;

// CLASS F_TABLE4V
// ===============
// ===============

class MATH_API f_table4V : public pred4v {
private:
	friend class table4V_diff_prec;
	friend class table4V_diff_real;
	/**< classes that have access to private members */

	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */
	math::vec1* _points4;
	/**< pointer to object of class vec1 (ordered magnitude vector of
	size l defining the table 1st inputs). */
	math::vec1* _points3;
	/**< pointer to object of class vec1 (ordered magnitude vector of
	size m defining the table 2nd inputs). */
	math::vec1* _points2;
	/**< pointer to object of class vec1 (ordered magnitude vector of
	size n defining the table 3rd inputs). */
	math::vec1* _points1;
	/**< pointer to object of class vec1 (ordered magnitude vector of
	size q defining the table 4th inputs). */
	math::vec4* _VVValues;
	/**< pointer to object of class vec4 (size l vector to pointers to
	size m vectors to pointers to size n vectors to pointers to size q vectors
	of vec1 objects) containing the table outputs */
	bool _del_flag;
	/**< flag that indicates if the _points4, _points3, _points2, _points1, and
	_VVValues vectors	shall be deleted by the destructor (1) or not (0) */

	math::logic::INTERP_MODE _interp_mode;
	/**< interpolation method name */
	math::interp* _interp;
	/**< pointer to interpolation method */
	bool _equi4;
	/**< true if _points4 vector equispaced, false otherwise */
	bool _equi3;
	/**< true if _points3 vector equispaced, false otherwise */
	bool _equi2;
	/**< true if _points2 vector equispaced, false otherwise */
	bool _equi1;
	/**< true if _points1 vector equispaced, false otherwise */
	double _points4_diff;
	/**< difference between consecutive _points4 */
	double _points3_diff;
	/**< difference between consecutive _points3 */
	double _points2_diff;
	/**< difference between consecutive _points2 */
	double _points1_diff;
	/**< difference between consecutive _points1 */

	math::table4V_diff* _functor_diff;
	/**< pointer to functor that solves differentials */
	math::vec4* _slopes_d4;
	/**< pointer to object of class vec4 (size l vector to pointer to
	size m vectors to pointers to size n vectors to pointers to size q vectors
	of vec1 objects) containing the differentials (slopes) with respect to the
	1st independent	magnitude. Derivate class coefficient always used. */
	math::vec4* _slopes_d3;
	/**< pointer to object of class vec4 (size l vector to pointer to
	size m vectors to pointers to size n vectors to pointers to size q vectors
	of vec1 objects) containing the differentials (slopes) with respect to the
	2nd independent	magnitude. Derivate class coefficient always used. */
	math::vec4* _slopes_d2;
	/**< pointer to object of class vec4 (size l vector to pointer to
	size m vectors to pointers to size n vectors to pointers to size q vectors
	of vec1 objects) containing the differentials (slopes) with respect to the
	3rd independent	magnitude. Derivate class coefficient always used. */
	math::vec4* _slopes_d1;
	/**< pointer to object of class vec4 (size l vector to pointer to
	size m vectors to pointers to size n vectors to pointers to size q vectors
	of vec1 objects) containing the differentials (slopes) with respect to the
	4th independent	magnitude. Derivate class coefficient always used. */
	math::hermite4v* _herm;
	/**< pointer to Hermite coefficients (if applicable) */

	math::range_checker* _checker4;
	/**< object controlling out of range for _points4, inactive by default */
	math::range_checker* _checker3;
	/**< object controlling out of range for _points3, active by default */
	math::range_checker* _checker2;
	/**< object controlling out of range for _points2, active by default */
	math::range_checker* _checker1;
	/**< object controlling out of range for _points1, active by default */
	math::pos_finder* _finder4;
	/**< object determining position of input with respect to _points4 */
	math::pos_finder* _finder3;
	/**< object determining position of input with respect to _points3 */
	math::pos_finder* _finder2;
	/**< object determining position of input with respect to _points2 */
	math::pos_finder* _finder1;
	/**< object determining position of input with respect to _points1 */

	void fill_up_slopes_lagrange_first_precompute();
	void fill_up_slopes_hermite();
	/**< fills up the _slopes_d4, _slopes_d3, _slopes_d2, and _slopes_d1 attributes based on
	_points4, _points3, _points2, _points1, and _VVValues */
	void fill_up_slopes_aux4_lagrange_first_precompute(math::vec4& slopes_d4,
													   unsigned short pos, 
													   unsigned short siz2,
													   unsigned short siz3,
													   unsigned short siz4);
	void fill_up_slopes_aux4_hermite(math::vec4& slopes_d4,
									 unsigned short index3,
									 unsigned short index2,
									 unsigned short index1);
	/**< fills up order "pos" of the _slopes_d4 4Dmatrix of differentials with respect to the
	first independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */
	void fill_up_slopes_aux3_lagrange_first_precompute(math::vec4& slopes_d3,
													   unsigned short pos, 
													   unsigned short siz2,
													   unsigned short siz3,
													   unsigned short siz4);
	void fill_up_slopes_aux3_hermite(math::vec4& slopes_d3,
									 unsigned short index4,
									 unsigned short index2,
									 unsigned short index1);
	/**< fills up order "pos" of the _slopes_d3 4Dmatrix of differentials with respect
	to the second independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */
	void fill_up_slopes_aux2_lagrange_first_precompute(math::vec4& slopes_d2,
													   unsigned short pos, 
													   unsigned short siz2,
													   unsigned short siz3,
													   unsigned short siz4);
	void fill_up_slopes_aux2_hermite(math::vec4& slopes_d2,
									 unsigned short index4, 
									 unsigned short index3,
									 unsigned short index1);
	/**< fills up order "pos" of the _slopes_d2 4Dmatrix of differentials with respect
	to the third independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */
	void fill_up_slopes_aux1_lagrange_first_precompute(math::vec4& slopes_d1,
													   unsigned short pos, 
													   unsigned short siz2,
													   unsigned short siz3,
													   unsigned short siz4);
	void fill_up_slopes_aux1_hermite(math::vec4& slopes_d1,
									 unsigned short index4, 
									 unsigned short index3,
									 unsigned short index2);
	/**< fills up order "pos" of the _slopes_d1 4Dmatrix of differentials with respect
	to the fourth independent magnitude with a matrix of size "siz2" by "siz3" by "siz4" */

	f_table4V();
	/**< empty constructor not implemented */
	f_table4V& operator=(const f_table4V& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_table4V(math::vec1* points4,
				math::vec1* points3,
				math::vec1* points2,
				math::vec1* points1,
				math::vec4* VVValues,
				math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on pointer to size l vector points4, pointer	to a
	size m vector points3, pointer to a size n vector points2, pointer to a
	size q vector point3, and pointer to a size l vector of size m vectors of size n
	vectors of size q vectors VVValues - deleted by destructor. */
	f_table4V(math::vec1& points4,
				math::vec1& points3,
				math::vec1& points2,
				math::vec1& points1,
				math::vec4& VVValues,
				math::logic::INTERP_MODE interp_mode = math::logic::lagrange_first_precompute);
	/**< constructor based on reference to size l vector points4, reference	to a
	size m vector points3, reference to a size n vector points2, reference to a
	size q vector points1, and reference to	a size l vector of size m vectors of 
	size n vectos of size q vectors VValues - not deleted by destructor. */
	f_table4V(const f_table4V&);
	/**< copy constructor */
	virtual ~f_table4V();
	/**< destructor */
	virtual f_table4V* clone() const;
	/**< cloner */
	bool operator==(const pred4v& op2) const;
	bool operator==(const f_table4V& op2) const;
	/**< overloaded operator == (equal) */
	bool operator!=(const f_table4V& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

	const math::vec1& get_points4() const {return *_points4;}
	/**< get the _points4 object as constant reference to read */
	const math::vec1& get_points3() const {return *_points3;}
	/**< get the _points3 object as constant reference to read */
	const math::vec1& get_points2() const {return *_points2;}
	/**< get the _points2 object as constant reference to read */
	const math::vec1& get_points1() const {return *_points1;}
	/**< get the _points1 object as constant reference to read */
	const math::vec4& get_VVValues() const {return *_VVValues;}
	/**< get the _VVValues object as constant reference to read */
	const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred4v above */
	const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred4v above */
	const math::logic::INTERP_MODE& get_interp_mode() const {return _interp_mode;}
	/**< get the interpolation mode to read */
	const math::interp& get_interp() const {return *_interp;}
	/**< get constant reference to interpolation method */
	const math::vec4& get_slopes_d4() const {return *_slopes_d4;}
	/**< return _slopes_d4 matrix to read. Only required to take a look at the 
	_slopes_d4 matrix in the unitary tests */
	const math::vec4& get_slopes_d3() const {return *_slopes_d3;}
	/**< return _slopes_d3 matrix to read. Only required to take a look at the 
	_slopes_d3 matrix in the unitary tests */
	const math::vec4& get_slopes_d2() const {return *_slopes_d2;}
	/**< return _slopes_d2 matrix to read. Only required to take a look at the 
	_slopes_d2 matrix in the unitary tests */
	const math::vec4& get_slopes_d1() const {return *_slopes_d1;}
	/**< return _slopes_d1 matrix to read. Only required to take a look at the 
	_slopes_d1 matrix in the unitary tests */
	bool get_equi4() const {return _equi4;}
	/**< returns true if _points4 vector equispaced, false otherwise */
	bool get_equi3() const {return _equi3;}
	/**< returns true if _points3 vector equispaced, false otherwise */
	bool get_equi2() const {return _equi2;}
	/**< returns true if _points2 vector equispaced, false otherwise */
	bool get_equi1() const {return _equi1;}
	/**< returns true if _points1 vector equispaced, false otherwise */
	const math::pos_finder& get_finder4() const {return *_finder4;}
	/**< returns constant reference to finder object */
	const math::pos_finder& get_finder3() const {return *_finder3;}
	/**< returns constant reference to finder object */
	const math::pos_finder& get_finder2() const {return *_finder2;}
	/**< returns constant reference to finder object */
	const math::pos_finder& get_finder1() const {return *_finder1;}
	/**< returns constant reference to finder object */
	const double& get_points4_diff() const {return _points4_diff;}
	/**< returns difference between consecutive _points4, 0. if not equispaced */
	const double& get_points3_diff() const {return _points3_diff;}
	/**< returns difference between consecutive _points3, 0. if not equispaced */
	const double& get_points2_diff() const {return _points2_diff;}
	/**< returns difference between consecutive _points2, 0. if not equispaced */
	const double& get_points1_diff() const {return _points1_diff;}
	/**< returns difference between consecutive _points1, 0. if not equispaced */
    const math::range_checker& get_checker1() const {return *_checker1;}
    /**< returns points1 out of range checker */
    const math::range_checker& get_checker2() const {return *_checker2;}
    /**< returns points2 out of range checker */
    const math::range_checker& get_checker3() const {return *_checker3;}
    /**< returns points3 out of range checker */
    const math::range_checker& get_checker4() const {return *_checker4;}
    /**< returns points4 out of range checker */    
    
	void activate_checker4();
	void deactivate_checker4();
	/**< Activates or deactives the out of range verification for _points4, which
	is inactive by default */
	void activate_checker3();
	void deactivate_checker3();
	/**< Activates or deactives the out of range verification for _points3, which
	is active by default */
	void activate_checker2();
	void deactivate_checker2();
	/**< Activates or deactives the out of range verification for _points2, which
	is active by default */
	void activate_checker1();
	void deactivate_checker1();
	/**< Activates or deactives the out of range verification for _points1, which
	is active by default */

	int compute_pos4(const double& input4) const;
	/**< Returns first position within _points4 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	int compute_pos3(const double& input3) const;
	/**< Returns first position within _points3 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	int compute_pos2(const double& input2) const;
	/**< Returns first position within _points2 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	int compute_pos1(const double& input1) const;
	/**< Returns first position within _points1 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	math::ratio* compute_ratio4(const double& input4,
									  const int& pos4) const;
	/**< Returns ratio of input4 with respect to the two points4 identified by pos4
	and pos4+1. */
	math::ratio* compute_ratio3(const double& input3,
									  const int& pos3) const;
	/**< Returns ratio of input3 with respect to the two points3 identified by pos3
	and pos3+1. */
	math::ratio* compute_ratio2(const double& input2,
									  const int& pos2) const;
	/**< Returns ratio of input2 with respect to the two points2 identified by pos2
	and pos2+1. */
	math::ratio* compute_ratio1(const double& input1,
									  const int& pos1) const;
	/**< Returns ratio of input1 with respect to the two points1 identified by pos1
	and pos1+1. */
	double compute_value(const int& pos4,
					   const int& pos3,
					   const int& pos2,
					   const int& pos1,
					   const math::ratio& ratio4,
					   const math::ratio& ratio3,
					   const math::ratio& ratio2,
					   const math::ratio& ratio1) const;
	/**< Fills up result magnitude by interpolating based on the input positions
	and ratios */
	double compute_diff(const int& pos4,
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
					  const double& input1_dt) const;
	/**< Fills up result differential based on positions and ratios */

    double value(const double& input4,
			   const double& input3,
			   const double& input2,
			   const double& input1) const;
	/**< see virtual function of class predicate4v above. */
	double d_dt(const double& input4,
			  const double& input3,
			  const double& input2,
			  const double& input1,
			  const double& input4_dt,
			  const double& input3_dt,
			  const double& input2_dt,
			  const double& input1_dt) const;
	/**< see virtual function of class pred4v above. */
}; // closes class f_table4V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABLE4V_DIFF
// ==================
// ==================

class MATH_API table4V_diff {
public:
	virtual void compute_diff(double& result,
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
							  const double& input1_dt) const = 0;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class table4V_diff

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABLE4V_DIFF_PREC
// =======================
// =======================

class MATH_API table4V_diff_prec: public table4V_diff {
private:
	math::f_table4V* _pred;
	/**< weak pointer to f_table4V class */
public:
	table4V_diff_prec(math::f_table4V& pred);
	/**< constructor based on equispaced four dimensional table */
	~table4V_diff_prec();
	/**< destructor */
	void compute_diff(double& result,
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
					  const double& input1_dt) const;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class table4V_diff_prec

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS TABLE4V_DIFF_REAL
// =======================
// =======================

class MATH_API table4V_diff_real: public table4V_diff {
private:
	math::f_table4V* _pred;
	/**< weak pointer to f_table4V class */
public:
	table4V_diff_real(math::f_table4V& pred);
	/**< constructor based on equispaced four dimensional table */
	~table4V_diff_real();
	/**< destructor */
	void compute_diff(double& result,
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
					  const double& input1_dt) const;
	/**< fill up result differential based on positions, inputs, and its partial
	differentials with time */
}; // closes class table4V_diff_real

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

} // closes namespace math

#endif





