#ifndef MATH_F_CONSTANT_H
#define MATH_F_CONSTANT_H

#include "../../math.h"
#include "pred0v.h"

/*
Derivate class of pred0v containing f_constant objects. Refer to pred0v
for more info.
*/

namespace math {

// CLASS F_CONSTANT
// ================
// ================

class MATH_API f_constant : public pred0v {
private:
	double _f0; 
	/**< constant value of the function, in standard units */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */
	
	f_constant();
	/**< empty constructor not implemented */
	f_constant& operator=(const f_constant& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_constant(const double& f0);
	/**< constructor based on constant value for output in standard units */
	f_constant(const f_constant&);
	/**< copy constructor */
	inline const double& get_f0() const	{return _f0;}
	/**< get the _f0 value to read */
	inline double value () const {return _f0;}
	/**< see virtual function of class pred0v above. */
	inline double d_dt() const	{return 0.;}
	/**< see virtual function of class pred0v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred0v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred0v above */
	virtual f_constant* clone() const;
	/**< cloner */
	bool operator==(const pred0v& op2) const;
	bool operator==(const f_constant& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_constant& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_constant

} // closes namespace math

#endif



