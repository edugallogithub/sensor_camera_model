#ifndef ATT_TEST_SCREW
#define ATT_TEST_SCREW

#include "ang_test.h"
#include <jail/unit_test.h>

/* This test manually checks the relationship between the screw, the dual quaternion, and the
 * rotation vector and translation. Change input values to see the result.
 * Reaches two important conclusions:
 * 1 --> Never go from screw to rot vector and translation directly as expression does not work for small rotations.
 * 2 --> Do not use acos when obtaining the screw from the dual quaternion. Instead, treat it as rotation
 *       vector with the appropriate formulas.
 */

namespace ang {
namespace test {
	
class Tscrew: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tscrew(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
	void test();

};

} // closes namespace test
} // closes namespace ang

#endif

