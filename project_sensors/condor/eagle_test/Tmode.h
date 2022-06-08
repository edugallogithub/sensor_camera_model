#ifndef EAGLE_TEST_MODE
#define EAGLE_TEST_MODE

#include "eagle_test.h"
#include "eagle/earth_eagle_att.h"
#include "eagle/earth_eagle_final.h"
#include "eagle/earth_eagle_osg.h"
#include "jail/unit_test.h"



namespace eagle {
namespace test {
	
class Tmode: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Tmode(jail::counter&);
	/**< execute tests and write results on console */
	void run(int* argc, char **argv); //  override;

	/**< simple test to verify that it indeed runs */
    static void test1(eagle::earth_eagle& Oeagle);
    static void test2(eagle::earth_eagle& Oeagle);


};

} // closes namespace test
} // closes namespace eagle


#endif

