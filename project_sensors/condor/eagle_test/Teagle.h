#ifndef EAGLE_TEST_EAGLE
#define EAGLE_TEST_EAGLE

#include "eagle_test.h"
#include "eagle/earth_eagle_att.h"
#include "eagle/earth_eagle_final.h"
#include "eagle/earth_eagle_osg.h"
#include "jail/unit_test.h"

/*
This file contains tests to verify the proper behavior of the files in the "eagle" library.
*/

namespace eagle {
namespace test {
	
class Teagle: public ::jail::unit_test {
public:
	/**< constructor based on counter */
	explicit Teagle(jail::counter&);
	/**< execute tests and write results on console */
	void run(int* argc, char **argv); //  override;

	/**< simple test to verify that it indeed runs */
    static void test_simple(eagle::earth_eagle& Oeagle);

    /**< simple test to visually verify that images really have altitude information */
    static void test_altitude(eagle::earth_eagle& Oeagle);

	/**< compares the different transformations within earth eagle executing them both with the osg
	 * and attitude processes, ensuring they are identical */
	void test_osg_versus_att(eagle::earth_eagle_att& Oeagle_att, eagle::earth_eagle_osg& Oeagle_osg, eagle::earth_eagle_final& Oeagle_fin);

    /**< compares the time taken to employ earth eagle with both the osg and attitude transformations,
     * so the fastest one can be set as the default */
	void test_speed_osg_versus_att(eagle::earth_eagle_att& Oeagle_att, eagle::earth_eagle_osg& Oeagle_osg);

    /**< processes the same images with different converters */
    static void test_save_to_stream(eagle::earth_eagle& Oeagle1, eagle::earth_eagle& Oeagle2, eagle::earth_eagle& Oeagle3);

    /**< compares the time taken to obtain a monochrome image in opencv with the different converting methods from osg */
    static void test_speed_converter(eagle::earth_eagle& Oeagle1, eagle::earth_eagle& Oeagle2, eagle::earth_eagle& Oeagle3);

    /**< covers the influence of the transformation between body and camera frames */
    void test_G_bc(eagle::earth_eagle& Oeagle, const std::string& st_case);
};

} // closes namespace test
} // closes namespace eagle


#endif

