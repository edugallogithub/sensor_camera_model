#include "f_null.h"

// CLASS F_NULL
// ============
// ============

const math::logic::PRED_NAME math::f_null::_name = math::logic::f_null;
/* predicate name */

const std::string math::f_null::_st_name = "f_null";
/* predicate name string */

math::f_null::f_null() {}
/* constructor */

math::f_null::f_null(const math::f_null&) {}
/* copy constructor */

math::f_null* math::f_null::clone() const {
	return new f_null(*this);
}
/* cloner */

bool math::f_null::operator==(const math::pred0v& op2) const {
	return (op2.get_name() == math::logic::f_null) ?
		(*this == static_cast<const f_null&>(op2)) : false;	
}
/* overloaded operator == (equal) */


