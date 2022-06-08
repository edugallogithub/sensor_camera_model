#ifndef MATH_TEST_PRED
#define MATH_TEST_PRED

#include "math_test.h"
#include <jail/unit_test.h>
#include <math/classifiers.h>

/*
This file test the generic functions and related predicates.
The strategy is to manually cover all possible cases.

Simple instances of each class are created and then the value and d_dt
functions are tested at a few key points. The units of both the output and
input variables are then changed and the results verified. The input
variables are also checked to see if their identity is recognized.

Special emphasis is placed in guaranteeing a proper interpolation
between magnitudes of type longitude or bearing.

*/

namespace math {
namespace test {
	
class Tpre_fun: public ::jail::unit_test {
public:
	Tpre_fun(jail::counter&);
	/**< constructor based on counter */
	void run();	
	/**< execute tests and write results on console */

	void test0v_f_null();			// test f_null, pred0v and fun0
	void test0v_f_constant();		// test f_constant, pred0v and fun0
	void test1v_f_lineal();			// test f_lineal, pred1v and fun1
	void test1v_f_parabolic();		// test f_parabolic, pred1v and fun1
	void test1v_f_cubic();			// test f_cubic, pred1v and fun1
	void test1v_f_steps();			// test f_steps, pred1v and fun1
	void test1v_f_table1V();		// test f_table1V, pred1v and fun1
	void test1v_f_table1V_spl();	// test f_table1Vspl, pred1v and fun1
	void test1v_f_table1Veq();		// test f_table1Veq, pred1v and fun1
	void test2v_f_lineal_double();	// test f_lineal_double, pred2v and fun2
	void test2v_f_table2V();		// test f_table2V, pred2v and fun2
	void test2v_f_table2Veq();		// test f_table2Veq, pred2v and fun2
	void test2v_f_tabular2V();		// test f_tabular2V, pred2v and fun2
	void test3v_f_lineal_triple();	// test f_lineal_triple, pred3v and fun3
	void test3v_f_table3V();		// test f_table3V, pred3v and fun3
	void test3v_f_table3Veq();		// test f_table3Veq, pred3v and fun3
	void test3v_f_tabular3V();		// test f_tabular3V, pred3v and fun3
	void test4v_f_table4V();		// test f_table4V, pred4v and fun4
	void test4v_f_table4Veq();		// test f_table4Veq, pred4v and fun4
};

}; // closes namespace test
}; // closes namespace math
#endif

