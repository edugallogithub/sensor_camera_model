#ifndef MATH_F_PARABOLIC_H
#define MATH_F_PARABOLIC_H

#include "../../math.h"
#include "pred1v.h"

/*
Derivate class of pred1v containing f_parabolic objects. Refer to pred1v
for more info.
*/

namespace math {

// CLASS F_PARABOLIC
// =================
// =================

class MATH_API f_parabolic : public pred1v {
private:
	double _f0;	/**< polynomial coefficient of grade 0 [no units] */
	double _f1;	/**< polynomial coefficient of grade 1 [no units] */
	double _f2;	/**< polynomial coefficient of grade 2 [no units] */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_parabolic();
	/**< empty constructor not implemented */
	f_parabolic& operator=(const f_parabolic& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_parabolic(const double f0,
				const double f1,
				const double f2);
	/**< constructor based on parabolic function coefficients for input
	and output in standard units */
	f_parabolic(const f_parabolic&);
	/**< copy constructor */
	inline const double& get_f0() const {return _f0;}
	/**< get the _f0 value to read */
	inline const double& get_f1() const {return _f1;}
	/**< get the _f1 value to read */
	inline const double& get_f2() const {return _f2;}
	/**< get the _f2 value to read */
    double value(const double& input) const;
	/**< see virtual function of class pred1v above. */
    double d_dt(const double& input, const double& input_dt) const;
	/**< see virtual function of class pred1v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred1v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred1v above */
	virtual f_parabolic* clone() const;
	/**< cloner */
	bool operator==(const pred1v& op2) const;
	bool operator==(const f_parabolic& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_parabolic& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_parabolic

} // closes namespace math

#endif



