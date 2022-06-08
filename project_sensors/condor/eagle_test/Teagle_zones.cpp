#include "Teagle_zones.h"
#include "eagle/converter.h"
#include "eagle/earth_eagle_final.h"
#include "ang/rotate/euler.h"
#include "ang/transform/speu_rodrigues.h"
#include "env/geo.h"
#include "env/coord.h"
#include "acft/sens/camera.h"
#include "acft/acft/iner.h"
#include "math/logic/seeder.h"

eagle::test::Teagle_zones::Teagle_zones(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void eagle::test::Teagle_zones::run(int* argc, char **argv) {
	::jail::unit_test::run();

    env::geo_mix Ogeo(env::logic::zone_default);

    auto PG_bc = new ang::speu_rodrigues(ang::rodrigues({90.0 * math::constant::D2R(), 0., 0.}), {0., 0., 0.});
    sens::camera_test Ocam1(1024, 768, 37.9233, PG_bc, false);

    eagle::earth_eagle_final Oeagle1(Ogeo, Ocam1, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "zones", "arcgisonline.earth", argc, argv);

    env::geodetic_coord Ogisu;
    ang::euler Oeuler_nedbfs(0., 0., 0.);

    // Everglades
    double lambda_rad = 279.145 * math::constant::D2R();
    double phi_rad    = 25.9   * math::constant::D2R();
    double h_m        = +1500.0;

    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "everglades_low", 1.0);
    h_m        = +2500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "everglades_mid", 1.0);
    h_m        = +3500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "everglades_high", 1.0);

    // Farm
    lambda_rad = 271.0336 * math::constant::D2R();
    phi_rad    = 38.9628   * math::constant::D2R();
    h_m        = +1500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "farm_low", 1.0);
    h_m        = +2500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "farm_mid", 1.0);
    h_m        = +3500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "farm_high", 1.0);

    // Forest
    lambda_rad = 278.3411 * math::constant::D2R();
    phi_rad    = 37.5319   * math::constant::D2R();
    h_m        = +1500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "forest_low", 1.0);
    h_m        = +2500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "forest_mid", 1.0);
    h_m        = +3500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "forest_high", 1.0);

    // Desert
    lambda_rad = 247.1631 * math::constant::D2R();
    phi_rad    = 31.9719   * math::constant::D2R();
    h_m        = +1500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "desert_low", 1.0);
    h_m        = +2500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "desert_mid", 1.0);
    h_m        = +3500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "desert_high", 1.0);

    // Mountains
    lambda_rad = 239.028 * math::constant::D2R();
    phi_rad    = 48.4469   * math::constant::D2R();
    h_m        = +3000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "mountains_low", 1.0);
    h_m        = +4000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "mountains_mid", 1.0);
    h_m        = +5000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "mountains_high", 1.0);

    // Urban
    lambda_rad = 242.0544 * math::constant::D2R();
    phi_rad    = 33.8042   * math::constant::D2R();
    h_m        = +4500.0;
    h_m        = +1500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "urban_low", 1.0);
    h_m        = +2500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "urban_mid", 1.0);
    h_m        = +3500.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle1.process_single_input(Ogisu, Oeuler_nedbfs, "urban_high", 1.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    math::seeder Oseeder(1);
    acft::iner Oiner;
    sens::camera Ocam2(Oseeder, 512, 384, 49.2255, Oiner, true, false);
    eagle::earth_eagle_final Oeagle2(Ogeo, Ocam2, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "zones", "arcgisonline.earth", argc, argv);

    // desert generic
    lambda_rad = 247.1631 * math::constant::D2R();
    phi_rad    = 31.9719   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "desert_generic", 1.0);

    // desert specific
    lambda_rad = -112.2861 * math::constant::D2R();
    phi_rad    = 31.9158   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "desert_specific", 1.0);

    // farm generic
    lambda_rad = -87.8711 * math::constant::D2R();
    phi_rad    = 38.8610   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "farm_generic", 1.0);

    // farm specific
    lambda_rad = -87.5379 * math::constant::D2R();
    phi_rad    = 38.6877   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "farm_specific", 1.0);

    // forest generic
    lambda_rad = -72.4497 * math::constant::D2R();
    phi_rad    = 43.3431   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "forest_generic", 1.0);

    // forest specific
    lambda_rad = -72.4922 * math::constant::D2R();
    phi_rad    = 43.3771   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "forest_specific", 1.0);

    // mix generic
    lambda_rad = -88.9973 * math::constant::D2R();
    phi_rad    = 34.6792   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "mix_generic", 1.0);

    // mix specific
    lambda_rad = -89.0792 * math::constant::D2R();
    phi_rad    = 34.7103   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "mix_specific", 1.0);

    // prairie generic
    lambda_rad = -80.8985 * math::constant::D2R();
    phi_rad    = 25.9174   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "prairie_generic", 1.0);

    // prairie specific
    lambda_rad = -80.6685 * math::constant::D2R();
    phi_rad    = 25.9083   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "prairie_specific", 1.0);

    // urban generic
    lambda_rad = -118.2002 * math::constant::D2R();
    phi_rad    = 33.9244   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "urban_generic", 1.0);

    // urban specific
    lambda_rad = -118.2196 * math::constant::D2R();
    phi_rad    = 34.0075   * math::constant::D2R();
    h_m        = +2700.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    Oeagle2.process_single_input(Ogisu, Oeuler_nedbfs, "urban_specific", 1.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

	finished();
}
/* execute tests and write results on console */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////











