#include "Tmode.h"
#include "math/logic/timer.h"
#include "math/logic/seeder.h"
#include "ang/rotate/euler.h"
#include "eagle/converter.h"
#include "env/geo.h"
#include "env/coord.h"
#include "acft/sens/camera.h"
#include "acft/acft/iner.h"

eagle::test::Tmode::Tmode(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void eagle::test::Tmode::run(int* argc, char **argv) {
	::jail::unit_test::run();

    env::geo_mix Ogeo(env::logic::zone_default);

    double lambda_rad = -2.13426900;
    double phi_rad    = +0.83643580;
    double h_m        = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);

    double d2r = math::constant::D2R();
    double psi_rad    = 8.0 * d2r;
    double theta_rad  = 10.0 * d2r;
    double xi_rad     = 4.0 * d2r;
    ang::euler euler_nedcam(psi_rad, theta_rad, xi_rad);

    {
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
//        eagle::earth_eagle_final Oeagle(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_vis_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_VIS_NOT", "WA_Everett.earth", argc, argv);
//
//        Oeagle.process_single_input(Ogisu, euler_nedcam, "test_vis_not", 1.0);
//        std::cout << "H1H1H1H1" << std::endl;
//        cv::Mat Omat;
//        Oeagle.obtain_image(Omat);
    }
//    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
    {
        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
        eagle::earth_eagle_final Oeagle(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_SAV_NOT", "WA_Everett.earth", argc, argv);

        Oeagle.process_single_input(Ogisu, euler_nedcam, "test_sav_not", 1.0);
        std::cout << "H2H2H2H2" << std::endl;
        cv::Mat Omat;
        Oeagle.obtain_image(Omat,"test_sav_not");
    }
    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
    {
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
//        eagle::earth_eagle_final Oeagle(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_str_opencv_vis, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_VIS_VIS", "WA_Everett.earth", true, argc, argv);

//        Oeagle.process_single_input(Ogisu, euler_nedcam, "test_vis_vis", 1.0);
//        std::cout << "H3H3H3H3" << std::endl;
//        cv::Mat Omat;
//        Oeagle.obtain_image(Omat);

//        cv::imshow("Viewer", Omat);
//        cv::waitKey(1000); // waits the minimum amount [1 msec] so image shows
    }
    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
    {
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
//        eagle::earth_eagle_final Oeagle(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_str_opencv_sav, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_VIS_SAV", "WA_Everett.earth", true, argc, argv);
//
//        Oeagle.process_single_input(Ogisu, euler_nedcam, "test_vis_sav", 1.0);
//        std::cout << "H4H4H4H4" << std::endl;
//        cv::Mat Omat;
//        Oeagle.obtain_image(Omat);
//
//        cv::imshow("Viewer", Omat);
//        cv::waitKey(1000); // waits the minimum amount [1 msec] so image shows
    }
    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;





	finished();
}
/* execute tests and write results on console */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Tmode::test1(eagle::earth_eagle& Oeagle) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad    = +0.83643580;
    double h_m        = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad    = 8.0 * d2r;
    double theta_rad  = 10.0 * d2r;
    double xi_rad     = 4.0 * d2r;
    ang::euler euler_nedcam(psi_rad, theta_rad, xi_rad);

    Oeagle.process_single_input(Ogisu, euler_nedcam, "test_simple", 1.0);
} // closes test_simple



void eagle::test::Tmode::test2(eagle::earth_eagle& Oeagle) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad    = +0.83643580;
    double h_m        = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad    = 8.0 * d2r;
    double theta_rad  = 10.0 * d2r;
    double xi_rad     = 4.0 * d2r;
    ang::euler euler_nedcam(psi_rad, theta_rad, xi_rad);

    Oeagle.process_single_input(Ogisu, euler_nedcam, "test_simple", 1.0);
} // closes test_simple




///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////









