#ifndef MATH_F_NULL_H
#define MATH_F_NULL_H

#include "../../math.h"
#include "pred0v.h"

/*
Derivate class of pred0v containing f_null objects. Refer to pred0v
for more info.
*/

namespace math {

// CLASS F_NULL
// ============
// ============

class MATH_API f_null : public pred0v {
private:
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_null& operator=(const f_null& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_null();
	/**< constructor */
	f_null(const f_null&);
	/**< copy constructor */
	inline double value () const {return 0.;}
	/**< see virtual function of class pred0v above. */
	inline double d_dt() const {return 0.;}
	/**< see virtual function of class pred0v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred0v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred0v above */
	virtual f_null* clone() const;
	/**< cloner */
	bool operator==(const pred0v& op2) const;
	inline bool operator==(const f_null& op2) const 
	{return true;}
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_null& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_null

} // closes namespace math

#endif



