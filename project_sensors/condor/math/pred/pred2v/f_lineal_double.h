#ifndef MATH_F_LINEAL_DOUBLE_H
#define MATH_F_LINEAL_DOUBLE_H

#include "../../math.h"
#include "pred2v.h"

/*
Derivate class of pred2v containing f_lineal_double objects. Refer to pred2v
for more info.
*/

namespace math {

// CLASS F_LINEAL_DOUBLE
// =====================
// =====================

class MATH_API f_lineal_double : public pred2v {
private:
	double _f0;
	/**< polynomial coefficient of grade 0 (x) and 0 (y) [no units] */
	double _f1x;
	/**< polynomial coefficient of grade 1 (x) and 0 (y) [no units] */
	double _f1y;
	/**< polynomial coefficient of grade 0 (x) and 1 (y) [no units] */
	double _f1xy;
	/**< polynomial coefficient of grade 1 (x) and 1 (y) [no units] */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_lineal_double();
	/**< empty constructor not implemented */
	f_lineal_double& operator=(const f_lineal_double& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_lineal_double(const double f0,
					const double f1x,
					const double f1y, 
					const double f1xy);	
	/**< constructor based on lineal function coefficients for input and
	output in standard units */
	f_lineal_double(const f_lineal_double&);
	/**< copy constructor */
	inline const double& get_f0() const {return _f0;}
	/**< get the _f0 value to read */
	inline const double& get_f1x() const {return _f1x;}
	/**< get the _f1 value to read */
	inline const double& get_f1y() const {return _f1y;}
	/**< get the _f1y value to read */
	inline const double& get_f1xy() const {return _f1xy;}
	/**< get the _f1xy value to read */
    double value(const double& input1,
                 const double& input2) const;
	/**< see virtual function of class pred2v above. */
    double d_dt(const double& input1,
                const double& input2,
                const double& input1_dt,
                const double& input2_dt) const;
	/**< see virtual function of class pred2v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred2v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred2v above */
	virtual f_lineal_double* clone() const;
	/**< cloner */
	bool operator==(const pred2v& op2) const;
	bool operator==(const f_lineal_double& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_lineal_double& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_lineal_double

} // closes namespace math

#endif



