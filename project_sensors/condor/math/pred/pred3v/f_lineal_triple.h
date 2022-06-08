#ifndef MATH_F_LINEAL_TRIPLE_H
#define MATH_F_LINEAL_TRIPLE_H

#include "../../math.h"
#include "pred3v.h"

/*
Derivate class of pred3v containing f_lineal_triple objects. Refer to pred3v
for more info.
*/

namespace math {

// CLASS F_LINEAL_TRIPLE
// =====================
// =====================

class MATH_API f_lineal_triple : public pred3v {
private:
	double _f0;
	/**< polynomial coefficient of grade 0 (x), 0 (y), and 0 (z) [no units] */
	double _f1x;
	/**< polynomial coefficient of grade 1 (x), 0 (y), and 0 (z) [no units] */
	double _f1y;
	/**< polynomial coefficient of grade 0 (x), 1 (y), and 0 (z) [no units] */
	double _f1z;
	/**< polynomial coefficient of grade 0 (x), 0 (y), and 1 (z) [no units] */
	double _f1xy;
	/**< polynomial coefficient of grade 1 (x), 1 (y), and 0 (z) [no units] */
	double _f1xz;
	/**< polynomial coefficient of grade 1 (x), 0 (y), and 1 (z) [no units] */
	double _f1yz;
	/**< polynomial coefficient of grade 0 (x), 1 (y), and 1 (z) [no units] */
	double _f1xyz;
	/**< polynomial coefficient of grade 1 (x), 1 (y), and 1 (z) [no units] */
	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	f_lineal_triple();
	/**< empty constructor not implemented */
	f_lineal_triple& operator=(const f_lineal_triple& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
	f_lineal_triple(const double f0,
					const double f1x,
					const double f1y,
					const double f1z,
					const double f1xy,
					const double f1xz,
					const double f1yz,
					const double f1xyz);
	/**< constructor based on lineal function coefficients for input and
	output in standard units *//**< constructor */
	f_lineal_triple(const f_lineal_triple&);
	/**< copy constructor */
	inline const double& get_f0() const {return _f0;}
	/**< get the _f0 value to read */
	inline const double& get_f1x() const {return _f1x;}
	/**< get the _f1x value to read */
	inline const double& get_f1y() const {return _f1y;}
	/**< get the _f1y value to read */
	inline const double& get_f1z() const {return _f1z;}
	/**< get the _f1z value to read */
	inline const double& get_f1xy() const {return _f1xy;}
	/**< get the _f1xy value to read */
	inline const double& get_f1xz() const {return _f1xz;}
	/**< get the _f1xz value to read */
	inline const double& get_f1yz() const {return _f1yz;}
	/**< get the _f1yz value to read */
	inline const double& get_f1xyz() const {return _f1xyz;}

	/**< get the _f1xyz value to read */
    double value(const double& input1,
			   const double& input2,
			   const double& input3) const;
	/**< see virtual function of class pred3v above. */
    double d_dt(const double& input1, const double& input2, const double& input3,
                const double& input1_dt, const double& input2_dt, const double& input3_dt) const;
	/**< see virtual function of class pred3v above. */
	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function ofvalue class pred3v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred3v above */
	virtual f_lineal_triple* clone() const;
	/**< cloner */
	bool operator==(const pred3v& op2) const;
	bool operator==(const f_lineal_triple& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_lineal_triple& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

}; // closes f_lineal_triple class

} // closes namespace math

#endif





