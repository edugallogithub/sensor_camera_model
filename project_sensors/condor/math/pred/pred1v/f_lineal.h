#ifndef MATH_F_LINEAL_H
#define MATH_F_LINEAL_H

#include "../../math.h"
#include "pred1v.h"

/*
Derivate class of pred1v containing f_lineal objects. Refer to pred1v
for more info.
*/

namespace math {

// CLASS F_LINEAL
// ==============
// ==============

class MATH_API f_lineal : public pred1v {
private:
	double _f0;	
	/**< polynomial coefficient of grade 0 [no units] */
	double _f1;
	/**< polynomial coefficient of grade 1 [no units] */
	bool _identity;
	/**< true if _f0 == 0. and _f1 == 1. (identity function) */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_lineal();
	/**< empty constructor not implemented */
	f_lineal& operator=(const f_lineal& op2);
	/**< overloaded operator = (assignment) not implemented */
	void fill_identity();
	/**< sets the _identity flag */
public:
	f_lineal(const double f0,
			 const double f1);
	/**< constructor based on lineal function coefficients for input and
	output in standard units */
	f_lineal(const double x1,
			 const double x2,
			 const double y1,
			 const double y2);
	/**< constructor based on input variable values x1 and x2 and output
	variable values y1 and y2 in standard units */
	f_lineal(const f_lineal&);
	/**< copy constructor */
	inline const double& get_f0() const {return _f0;}
	/**< get the _f0 value to read */
	inline const double& get_f1() const {return _f1;}
	/**< get the _f1 value to read */
    double value(const double& input) const;
	/**< see virtual function of class pred1v above. */
    double d_dt(const double& input, const double& input_dt) const;
	/**< see virtual function of class pred1v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred1v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred1v above */
	virtual f_lineal* clone() const;
	/**< cloner */
	bool operator==(const pred1v& op2) const;
	bool operator==(const f_lineal& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_lineal& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */
	bool is_identity() const {return _identity;}
	/**< fun1 functions with f_lineal as predicate, f0 == 0 and f1 ==1, are
	fully equivalent to a magnitude as they	always return the input magnitude,
	but they are difficult to detect from the outside. This method returns 
	false except for that specific case. Only overloaded where it applies. */

}; // closes class f_lineal

} // closes namespace math

#endif



