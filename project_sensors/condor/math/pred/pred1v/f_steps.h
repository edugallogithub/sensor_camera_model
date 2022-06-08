#ifndef MATH_F_STEPS_H
#define MATH_F_STEPS_H

#include "../../math.h"
#include "pred1v.h"

/*
Derivate class of pred1v containing f_steps objects. Refer to pred1v
for more info.
*/

namespace math {

// CLASS F_STEPS
// =============
// =============

class MATH_API f_steps : public pred1v {
private:
	math::vec1* _points;
	/**< pointer to object of class vec1 (ordered magnitude vector of size
	n+1 defining n intervals). */
	math::vec1* _values;
	/**< pointer to object of class vec1 (ordered magnitude vector of size
	n containing the value in each interval). */
	bool _del_flag;
	/**< flag that indicates if the _points and _values vectors shall be deleted
	by the destructor (true) or not (false) */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_steps();
	/**< empty constructor not implemented */
	f_steps& operator=(const f_steps& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_steps(math::vec1* points,
			math::vec1* values);
	/**< constructor based on pointers - deleted by destructor */
	f_steps(math::vec1& points,
			math::vec1& values);
	/**< constructor based on references - not deleted by destructor */
	f_steps(const f_steps&);
	/**< copy constructor */
	virtual ~f_steps();
	/**< destructor */
    double value(const double& input) const;
	/**< see virtual function of class pred1v above. */
    double d_dt(const double& input, const double& input_dt) const;
	/**< see virtual function of class pred1v above. */
	inline const math::vec1& get_points() const {return *_points;}
	/**< get the _points object as constant reference to read */
	inline const math::vec1& get_values() const {return *_values;}
	/**< get the _values object as constant reference to read */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred1v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred1v above */
	virtual f_steps* clone() const;
	/**< cloner */
	bool operator==(const pred1v& op2) const;
	bool operator==(const f_steps& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_steps& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes class f_steps

} // closes namespace math

#endif



