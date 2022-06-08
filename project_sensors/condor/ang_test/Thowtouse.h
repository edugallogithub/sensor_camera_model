#ifndef ATT_TEST_HOWTOUSE
#define ATT_TEST_HOWTOUSE

#include "ang_test.h"
#include <jail/unit_test.h>

/*
 * This test is a practical example of a camera situation, and is indeed the
 * MOST INTERESTING and intended to act as a reference for how to employ the
 * attitude classes.
 */

namespace ang {
namespace test {
	
class Thowtouse: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Thowtouse(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
    void test_se3_practical();  // rotation plus translation practical application
};

} // closes namespace test
} // closes namespace ang

#endif

