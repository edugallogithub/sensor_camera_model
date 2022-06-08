

#ifndef MATH_F_TABLE1V_SPL_H
#define MATH_F_TABLE1V_SPL_H

#include "../../math.h"
#include "pred1v.h"

/*
Derivate class of pred1v containing f_table1V_spl objects. Refer to pred1v
for more info. AVOID USE IF POSSIBLE. OTHER TABLES PREFERRED. RESULTS NOT
GUARANTEED.
*/

namespace math {

// CLASS F_TABLE1V_SPL
// ===================
// ===================

class MATH_API f_table1V_spl : public pred1v {
private:
	std::vector<double>* _points1;
	/**< pointer to ordered vector of size n containing the table inputs */
	std::vector<double>* _values;
	/**< pointer to vector of size n containing the table outputs */
	std::vector<double>* _diffdiff;
	/**< pointer to vector of size n containing the second differentials of the 
	functions. */
	std::vector<double>* _dist;
	/**< pointer to vector of doubles of size n-1 containing the distances
	between each pair of consecutive inputs. */

	static const math::logic::PRED_NAME _name;
	/**< predicate name */
	static const std::string _st_name;
	/**< predicate name string */

	void complete_spline(const std::vector<double>& points,
						 const std::vector<double>& values);
	/**< complements the constructor based on two vectors of the same size (n), 
	one containing the inputs and another the outputs of the function. */

	f_table1V_spl();
	/**< empty constructor not implemented */
	f_table1V_spl& operator=(const f_table1V_spl& op2);
	/**< overloaded operator = (assignment) not implemented */
public:
    f_table1V_spl(std::vector<double>* points,
				  std::vector<double>* values);
	/**< constructor based on pointers to the vectors representing the 
	table inputs and outputs. Assumes standard units. */
	f_table1V_spl(const f_table1V_spl&);
	/**< copy constructor */
	virtual ~f_table1V_spl();
	/**< destructor */

	virtual f_table1V_spl* clone() const;
	/**< cloner */
	bool operator==(const pred1v& op2) const;
	bool operator==(const f_table1V_spl& op2) const;
	/**< overloaded operator == (equal) */
	inline bool operator!=(const f_table1V_spl& op2) const
	{return !(*this == op2);}
	/**< overloaded operator != (not equal) */

	inline const logic::PRED_NAME& get_name() const {return _name;}
	/**< see virtual function of class pred1v above */
	inline const std::string& get_st_name() const {return _st_name;}
	/**< see virtual function of class pred1v above */
	inline const std::vector<double>& get_points1() const {return *_points1;}
	/**< get the _points object as constant reference to read */
	inline const std::vector<double>& get_values() const {return *_values;}
	/**< get the _values object as constant reference to read */

	int compute_pos1(const double& input) const;
	/**< Returns first position within _points1 vector that shall be employed when
	interpolating to obtain the result corresponding to the input magnitude. */
	math::ratio* compute_ratio1(const double& input1,
									  const int& pos1) const; 
	/**< Returns ratio of input1 with respect to the two points1 identified by pos1 and
	pos1+1. */
	void compute_value(double& result,
					   const int& pos1, 
					   const math::ratio& ratio1) const;
	/**< Fills up result magnitude by interpolating based on the input position and 
	ratio */
	void compute_diff(double& result,
					  const int& pos1, 
					  const double& input1,
					  const double& input1_dt) const;
	/**< Fills up result differential based on position and ratio */

    double value(const double& input) const;
	/**< see virtual function of class predicate1v above. */
    double d_dt(const double& input, const double& input_dt) const;
	/**< see virtual function of class pred1v above. */

}; // closes class f_table1V_spl

} // closes namespace math

#endif




