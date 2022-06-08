#ifndef MATH_F_CUBIC_H
#define MATH_F_CUBIC_H

#include "../../math.h"
#include "pred1v.h"

/*
Derivate class of pred1v containing f_cubic objects. Refer to pred1v
for more info.
*/

namespace math {

// CLASS F_CUBIC
// =============
// =============

class MATH_API f_cubic : public pred1v {
private:
	double _f0;	/**< polynomial coefficient of grade 0 [no units] */
	double _f1;	/**< polynomial coefficient of grade 1 [no units] */
	double _f2;	/**< polynomial coefficient of grade 2 [no units] */
	double _f3;	/**< polynomial coefficient of grade 3 [no units] */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_cubic();
	/**< empty constructor not implemented */
	f_cubic& operator=(const f_cubic& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_cubic(const double f0,
			const double f1,
			const double f2,
			const double f3);
	/**< constructor based on cubic function coefficients for input
	and output in standard units */
	f_cubic(const f_cubic&);
	/**< copy constructor */
	inline const double& get_f0() const {return _f0;}
	/**< get the _f0 value to read */
	inline const double& get_f1() const {return _f1;}
	/**< get the _f1 value to read */
	inline const double& get_f2() const {return _f2;}
	/**< get the _f2 value to read */
	inline const double& get_f3() const {return _f3;}
	/**< get the _f3 value to read */

    double value(const double& input) const;
	/**< see virtual function of class pred1v above. */
    double d_dt(const double& input, const double& input_dt) const;
	/**< see virtual function of class pred1v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred1v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred1v above */
	virtual f_cubic* clone() const;
	/**< cloner */
	bool operator==(const pred1v& op2) const;
	bool operator==(const f_cubic& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_cubic& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_cubic

} // closes namespace math

#endif



