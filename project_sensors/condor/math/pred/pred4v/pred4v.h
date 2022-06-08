

#ifndef MATH_PRED4V_H
#define MATH_PRED4V_H

#include "../../math.h"
#include "../pred.h"

/*
This file contains the base class pred4v for those function predicates that have
four independent magnitudes. The derivate classes are in other files in this
same folder.It is structured as a virtual base class (pred4v) with
its derived implementations.
- f_table4V containing (l) 1st input magnitudes, (m) 2nd input magnitudes,
  (n) 3rd input magnitudes, (q) 4th input magnitudes, and (l x m x n x q) output
  magnitudes representing the function value at each point. Different
  interpolation methods allowed.

They all share the methods value(x,y,z,v) and d_dt(x,y,z,v,dx_dt,dy_dt,dz_dt,dv_dt).
Their results are:
- f_table4V
	value(x,y,z,v)				  = employs linear interpolation between points.
	d_dt(x,y,z,v,
	     dx_dt,dy_dt,dz_dt,dv_dt) = d_dx * dx_dt + d_dy * dy_dt + d_dz * dz_dt +
									d_dv * dv_dt, where d_dx, d_dy, d_dz, and d_dv
									are the partial	differentials computed by linear
									interpolation
	Valid for any value of x (which generally represents altitude) but error if
	y, z, or v (which generally represent longitude, latitude, and time) outside
	intervals defined by _points.

SLOW FUNCTIONS:
- Determining the position within the _points1, _points2, _points3, and _points4
  vectors of the input values is expensive for f_table4V, as four complete searches
  are required. It is much cheaper for f_table4Veq, as the points are equispaced.
  Nothing can be done to accelerate this.
- The position is determined independently when calling the value and d_dt
  methods, which is redundant. It may be possible to modify the interface
  so the search is only executed once if both methods are called.

ACCURACY:
- No method presents any degradation in accuracy

LEFT:
- The f_table4V and f_table4Veq predicates return an error if called with an input
  outside the predefined interval for their 2nd, 3rd, and 4th independent variables.
  This may cause problems when iterating, in which the input may go outside the
  normal range.

MODE OF USE:
Employ the constructors to create instances of the classes:
- f_table4V predicate(points1,points2,points3,points4Values) creates an object
  predicate of class f_table4Veq, with points1 being an l size equispaced
  increasing order vector of independent magnitudes, points2 being a different
  m size equispaced increasing order vector of independent magnitudes, points3
  being a different n size equispaced increasing order vector of independent
  magnitudes, points4 being a different q size equispaced incresing order vector
  of independent magnituds, and Values an l size vector of m size vectors of
  size n vectors of size q vectors with the resulting magnitudes an each point.
  Linear interpolation is employed between points.

The functions value() and d_dt() can then be directly used, although it is more
common for these objects to be the first input parameters of the fun constructors.

*/

namespace math {

// CLASS PRED4V
// ============
// ============

class MATH_API pred4v : public pred{
public:
	virtual double value(const double& input1,
					   const double& input2,
					   const double& input3,
					   const double& input4) const = 0;
	/**< evaluates the function at the reference magnitudes input1, input2, 
	input3, and input4, and	writes the result at the reference magnitude result.
	Only the magnitude value is inserted into result, the units are assummed to
	be OK. */
	virtual double d_dt(const double& input1, const double& input2,
					  const double& input3, const double& input4,
					  const double& input1_dt, const double& input2_dt,
					  const double& input3_dt, const double& input4_dt) const = 0;
	/**< evaluates the function differential with time at the reference
	magnitudes input1, input2, input3, and input4 and their differentials with time
	input1_dt, input2_dt, input3_dt, and input4_dt and writes the result at the
	reference magnitude result. */
	virtual const logic::PRED_NAME& get_name() const = 0;
	/**< see virtual function of class pred above */
	virtual const std::string& get_st_name() const = 0;
	/**< see virtual function of class pred above */
	virtual ~pred4v() {}
    /**< destructor */
	virtual bool operator==(const pred4v& op2) const = 0;
	/**< overloaded operator == (equal) */
	inline virtual bool operator!=(const pred4v& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */
}; // closes pred4v class

} // closes namespace math

#endif





