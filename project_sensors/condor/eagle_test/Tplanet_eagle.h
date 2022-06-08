#ifndef EAGLE_TEST_PLANET_EAGLE
#define EAGLE_TEST_PLANET_EAGLE

#include "eagle_test.h"
#include "eagle/planet_eagle.h"
#include "eagle/earth_eagle_att.h"

#include "jail/unit_test.h"

/*
This file contains tests to verify the proper behavior of the files in the "planet_eagle" library.
*/

namespace eagle {
namespace test {

class Tplanet_eagle: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tplanet_eagle(jail::counter&);
	/**< execute tests and write results on console */
	void run(int* argc, char **argv); //  override;

	/**< simple test to verify that it indeed runs */
    void test1(eagle::earth_eagle& Oearth_eagle, eagle::planet_eagle& Oplanet_eagle);


};

} // closes namespace test

} // closes namespace eagle

#endif

