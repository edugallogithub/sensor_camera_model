#ifndef EAGLE_TEST_EAGLE_ZONES
#define EAGLE_TEST_EAGLE_ZONES

#include "eagle_test.h"
#include "jail/unit_test.h"

/*
This file creates sample images from different zones in the US.
*/

namespace eagle {
namespace test {
	
class Teagle_zones: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Teagle_zones(jail::counter&);
	/**< execute tests and write results on console */
	void run(int* argc, char **argv);

};

} // closes namespace test
} // closes namespace eagle


#endif

