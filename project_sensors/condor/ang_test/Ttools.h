#ifndef ATT_TEST_TOOLS
#define ATT_TEST_TOOLS

#include "ang_test.h"
#include <jail/unit_test.h>

/*
 * This file contains tests to verify the proper behavior some methods of the tools class.
 */

namespace ang {
namespace test {
	
class Ttools: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Ttools(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;	

    /**< specific tests */
    void test_skew();
};

} // closes namespace test
} // closes namespace ang

#endif

