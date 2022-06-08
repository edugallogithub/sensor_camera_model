#include "Tplanet_eagle.h"

#include "eagle/converter.h"
#include "math/logic/seeder.h"
#include "ang/transform/speu_rodrigues.h"
#include "ang/rotate/euler.h"
#include "env/geo.h"
#include "env/coord.h"
#include "acft/sens/camera.h"
#include "acft/acft/iner.h"

eagle::test::Tplanet_eagle::Tplanet_eagle(jail::counter& Ocounter)
        : ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void eagle::test::Tplanet_eagle::run(int* argc, char **argv) {
    ::jail::unit_test::run();

    env::geo_mix Ogeo(env::logic::zone_default);
    auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{0., 0., 0.}), {0., 0., 0.});
    sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);

    eagle::planet_eagle    Oplanet_eagle(      Ocam, eagle::logic::image_jpg, false, eagle::logic::converter_osg_sav_opencv_no, 0.8, "planet_eagle", "WA_Everett_merged.earth", argc, argv);
    eagle::earth_eagle_att Oearth_eagle (Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "planet_eagle", "WA_Everett_merged.earth", argc, argv);

    test1(Oearth_eagle, Oplanet_eagle);

    finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void eagle::test::Tplanet_eagle::test1(eagle::earth_eagle& Oearth_eagle, eagle::planet_eagle& Oplanet_eagle) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad = -45.0 * d2r;
    double theta05_rad =  5.0 * d2r;
    double theta10_rad = 10.0 * d2r;
    double theta15_rad = 15.0 * d2r;
    double theta20_rad = 20.0 * d2r;
    ang::euler Oeuler_nedcam;

    Oeuler_nedcam.set_rad(psi_rad, 0.0, 0.0);
    Oearth_eagle.process_single_input(Ogisu, Oeuler_nedcam, "Enothing", 1.0);
    Oplanet_eagle.process_single_input(Ogisu, psi_rad, 0.0, "Pnothing");

    Oeuler_nedcam.set_rad(psi_rad, 0.0, theta05_rad);
    Oearth_eagle.process_single_input(Ogisu, Oeuler_nedcam, "Eroll05", 1.0);
    Oeuler_nedcam.set_rad(psi_rad, 0.0, theta10_rad);
    Oearth_eagle.process_single_input(Ogisu, Oeuler_nedcam, "Eroll10", 1.0);
    Oeuler_nedcam.set_rad(psi_rad, 0.0, theta15_rad);
    Oearth_eagle.process_single_input(Ogisu, Oeuler_nedcam, "Eroll15", 1.0);
    Oeuler_nedcam.set_rad(psi_rad, 0.0, theta20_rad);
    Oearth_eagle.process_single_input(Ogisu, Oeuler_nedcam, "Eroll20", 1.0);

    Oplanet_eagle.process_single_input(Ogisu, psi_rad, theta05_rad, "Ppitch05");
    Oplanet_eagle.process_single_input(Ogisu, psi_rad, theta10_rad, "Ppitch10");
    Oplanet_eagle.process_single_input(Ogisu, psi_rad, theta15_rad, "Ppitch15");
    Oplanet_eagle.process_single_input(Ogisu, psi_rad, theta20_rad, "Ppitch20");

} // closes test_simple

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

