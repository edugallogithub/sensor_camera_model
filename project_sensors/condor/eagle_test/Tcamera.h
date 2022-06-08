#ifndef EAGLE_TEST_CAMERA
#define EAGLE_TEST_CAMERA

#include "eagle_test.h"
#include "eagle/earth_eagle_att.h"
#include "eagle/earth_eagle_final.h"
#include "eagle/earth_eagle_osg.h"
#include "jail/unit_test.h"


namespace eagle {
namespace test {
	
class Tcamera: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tcamera(jail::counter&);
	/**< execute tests and write results on console */
	void run() override;

    static void test_camera();

};

} // closes namespace test

} // closes namespace eagle

#endif

