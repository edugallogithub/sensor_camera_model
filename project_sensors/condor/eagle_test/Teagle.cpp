#include "Teagle.h"
#include "math/logic/timer.h"
#include "math/logic/seeder.h"
#include "ang/rotate/euler.h"
#include "eagle/converter.h"
#include "env/geo.h"
#include "env/coord.h"
#include "acft/sens/camera.h"
#include "acft/acft/iner.h"

eagle::test::Teagle::Teagle(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void eagle::test::Teagle::run(int* argc, char **argv) {
	::jail::unit_test::run();

    env::geo_mix Ogeo(env::logic::zone_default);
    math::seeder Oseeder(1);


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

    {
//        // just check that eagle final works, this test does not save any image
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
//        eagle::earth_eagle_final Oeagle_Everett(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett", "WA_Everett.earth", argc, argv);
//        test_simple(Oeagle_Everett);
    }

    {
        // visually check that altitude is employed when using databases (both individual tiffs and merged tiffs)
        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
        sens::camera_test Ocam(1024, 768, 25.0, PG_bc, false);
        eagle::earth_eagle_final Oeagle_Rainier(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Rainier", "WA_Rainier.earth", argc, argv);
        test_altitude(Oeagle_Rainier);
//        eagle::earth_eagle_final Oeagle_Rainier_merged(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Rainier_merged", "WA_Rainier_merged.earth", argc, argv);
//        test_altitude(Oeagle_Rainier_merged);
    }

    {
//        // Numerically and visually compare the earth eagle when running with my attitude transformation and those within OpenSceneGraph
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues({90.0 * math::constant::D2R() + 0.5, -0.2, 0.3}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);   //sens::camera Ocam(1024, 768, 15.0, Oiner, Gq_bc);
//        eagle::earth_eagle_att   Oeagle_att(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett", "WA_Everett.earth", argc, argv);
//        eagle::earth_eagle_osg   Oeagle_osg(      Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett", "WA_Everett.earth", argc, argv);
//        eagle::earth_eagle_final Oeagle_fin(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett", "WA_Everett.earth", argc, argv);
//        test_osg_versus_att(Oeagle_att, Oeagle_osg, Oeagle_fin);
//        test_speed_osg_versus_att(Oeagle_att, Oeagle_osg); // maintain commented as it takes time

        // visually compare the results of the saved images --> they should be identical
//        Oeagle_osg.process_text_file("test_Everett.csv", "osg_");
//        Oeagle_att.process_text_file("test_Everett.csv", "att_");

//        eagle::earth_eagle_att   Oeagle_att_merged(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, 0.8, "Everett_merged", "WA_Everett_merged.earth", argc, argv);
//        eagle::earth_eagle_osg   Oeagle_osg_merged(      Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, 0.8, "Everett_merged", "WA_Everett_merged.earth", argc, argv);
//        eagle::earth_eagle_final Oeagle_fin_merged(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, 0.8, "Everett_merged", "WA_Everett_merged.earth", argc, argv);
//        test_osg_versus_att(Oeagle_att_merged, Oeagle_osg_merged, Oeagle_fin_merged);
    }

    {
//        // compare different options to visualize and save image in both osg and opencv
//        auto PG_bc = new ang::speu_rodrigues(ang::rodrigues(ang::euler{90., 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam(1024, 768, 15.0, PG_bc, false);
//        eagle::earth_eagle_final Oeagle1(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_sav,        math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "osg_to_opencv", "WA_Everett_merged.earth", argc, argv);
//        eagle::earth_eagle_final Oeagle2(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_sav_stream, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "osg_to_opencv", "WA_Everett_merged.earth", argc, argv);
//        eagle::earth_eagle_final Oeagle3(Ogeo, Ocam, eagle::logic::image_jpg, eagle::logic::converter_osg_vis_opencv_sav,        math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "osg_to_opencv", "WA_Everett_merged.earth", argc, argv);

//        test_save_to_stream(Oeagle1, Oeagle2, Oeagle3);
//        test_speed_converter(Oeagle1, Oeagle2, Oeagle3);
    }

    {
//        // verify the variations with the transformation between body and camera
//        auto PG_bc1 = new ang::speu_rodrigues(ang::rodrigues({90.0 * math::constant::D2R(), 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam1 (1024, 768, 15.0, PG_bc1, false);
//        eagle::earth_eagle_final Oeagle_Everett1(Ogeo, Ocam1, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett1, "test_Gbc_01");
//
//        auto PG_bc2 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam2 (1024, 768, 15.0, PG_bc2, false);
//        eagle::earth_eagle_final Oeagle_Everett2(Ogeo, Ocam2, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett2, "test_Gbc_02");
//
//        auto PG_bc3 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 5 * math::constant::D2R(), 0.}), {0., 0., 0.});
//        sens::camera_test Ocam3 (1024, 768, 15.0, PG_bc3, false);
//        eagle::earth_eagle_final Oeagle_Everett3(Ogeo, Ocam3, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett3, "test_Gbc_03");
//
//        auto PG_bc4 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 5 * math::constant::D2R(), 2 * math::constant::D2R()}), {0., 0., 0.});
//        sens::camera_test Ocam4 (1024, 768, 15.0, PG_bc4, false);
//        eagle::earth_eagle_final Oeagle_Everett4(Ogeo, Ocam4, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett4, "test_Gbc_04");
//
//        auto PG_bc5 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 0., 0.}), {0., 0., 0.});
//        sens::camera_test Ocam5 (1024, 768, 15.0, PG_bc5, false);
//        eagle::earth_eagle_final Oeagle_Everett5(Ogeo, Ocam5, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett5, "test_Gbc_05");
//
//        auto PG_bc6 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 0., 0.}), {1000., 0., 0.}); // 100 m from body to camera along x body
//        sens::camera_test Ocam6 (1024, 768, 15.0, PG_bc6, false);
//        eagle::earth_eagle_final Oeagle_Everett6(Ogeo, Ocam6, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett6, "test_Gbc_06");
//
//        auto PG_bc7 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 0., 0.}), {1000., 100., 0.}); // 100 m from body to camera along x body, and also y body
//        sens::camera_test Ocam7 (1024, 768, 15.0, PG_bc7, false);
//        eagle::earth_eagle_final Oeagle_Everett7(Ogeo, Ocam7, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett7, "test_Gbc_07");
//
//        auto PG_bc8 = new ang::speu_rodrigues(ang::rodrigues({135.0 * math::constant::D2R(), 0., 0.}), {100., 100., -500.});
//        sens::camera_test Ocam8 (1024, 768, 15.0, PG_bc8, false);
//        eagle::earth_eagle_final Oeagle_Everett8(Ogeo, Ocam8, eagle::logic::image_jpg, eagle::logic::converter_osg_sav_opencv_no, math::logic::folder_hard_disk, math::logic::folder_hard_disk, 0.8, "Everett_Gbc", "WA_Everett.earth", argc, argv);
//        test_G_bc(Oeagle_Everett8, "test_Gbc_08");
    }

	finished();
}
/* execute tests and write results on console */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Teagle::test_simple(eagle::earth_eagle& Oeagle) {
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

void eagle::test::Teagle::test_altitude(eagle::earth_eagle& Oeagle) {
    double d2r = math::constant::D2R();
    env::geodetic_coord Ogisu;
    ang::euler Oeuler_nedcam;
    double lambda_rad, phi_rad, h_m, psi_deg, theta_deg, xi_deg;

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.847 * d2r; //+46.847 * d2r;;
    h_m        = +60000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  =  0.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_001", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.847* d2r; //+46.847 * d2r;;
    h_m        = +50000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  =  0.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_002", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.847* d2r; //+46.847 * d2r;;
    h_m        = +40000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  =  0.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_003", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.847 * d2r; //+46.847 * d2r;;
    h_m        = +30000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  =  0.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_004", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.82 * d2r; //+46.847 * d2r;;
    h_m        = +26000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 10.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_005", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.795 * d2r; //+46.847 * d2r;;
    h_m        = +23000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 20.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_006", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.77 * d2r; //+46.847 * d2r;;
    h_m        = +20000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 30.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_007", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.745 * d2r; //+46.847 * d2r;;
    h_m        = +17000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 40.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_008", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.72 * d2r; //+46.847 * d2r;;
    h_m        = +14000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 50.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_009", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.695 * d2r; //+46.847 * d2r;;
    h_m        = +11000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 60.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_010", 1.0);

    lambda_rad = -121.747 * d2r;
    phi_rad    =  46.67 * d2r; //+46.847 * d2r;;
    h_m        = +8000.0;
    Ogisu.set(lambda_rad, phi_rad, h_m);
    psi_deg    =  0.0;
    theta_deg  = 70.0;
    xi_deg     =  0.0;
    Oeuler_nedcam.set_rad(psi_deg*d2r, theta_deg*d2r, xi_deg*d2r);
    Oeagle.process_single_input(Ogisu, Oeuler_nedcam, "test_rainier_011", 1.0);

    int STOP = 8;
} // closes test_altitude

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Teagle::test_osg_versus_att(eagle::earth_eagle_att& Oeagle_att, eagle::earth_eagle_osg& Oeagle_osg, eagle::earth_eagle_final& Oeagle_fin) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;
    double psi_rad = 38.0 * d2r;
    double theta_rad = 14.0 * d2r;
    double xi_rad = 9.0 * d2r;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    ang::euler Oeuler_nedcam(psi_rad, theta_rad, xi_rad);

    Oeagle_osg.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_osg", 1.0);
    Oeagle_att.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_att", 1.0);
    Oeagle_fin.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_fin", 1.0); // this one verify visually

    ang::homogeneous H_enuned = Oeagle_att.get_H_enuned();
    ang::homogeneous H_enuned2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_enuned());

    check("H_enuned  (0,0):     ", H_enuned()(0,0), H_enuned2()(0,0), 1e-12);
    check("H_enuned  (0,1):     ", H_enuned()(0,1), H_enuned2()(0,1), 1e-12);
    check("H_enuned  (0,2):     ", H_enuned()(0,2), H_enuned2()(0,2), 1e-12);
    check("H_enuned  (0,3):     ", H_enuned()(0,3), H_enuned2()(0,3), 1e-6);
    check("H_enuned  (1,0):     ", H_enuned()(1,0), H_enuned2()(1,0), 1e-12);
    check("H_enuned  (1,1):     ", H_enuned()(1,1), H_enuned2()(1,1), 1e-12);
    check("H_enuned  (1,2):     ", H_enuned()(1,2), H_enuned2()(1,2), 1e-12);
    check("H_enuned  (1,3):     ", H_enuned()(1,3), H_enuned2()(1,3), 1e-6);
    check("H_enuned  (2,0):     ", H_enuned()(2,0), H_enuned2()(2,0), 1e-12);
    check("H_enuned  (2,1):     ", H_enuned()(2,1), H_enuned2()(2,1), 1e-12);
    check("H_enuned  (2,2):     ", H_enuned()(2,2), H_enuned2()(2,2), 1e-12);
    check("H_enuned  (2,3):     ", H_enuned()(2,3), H_enuned2()(2,3), 1e-6);
    check("H_enuned  (3,0):     ", H_enuned()(3,0), H_enuned2()(3,0), 1e-12);
    check("H_enuned  (3,1):     ", H_enuned()(3,1), H_enuned2()(3,1), 1e-12);
    check("H_enuned  (3,2):     ", H_enuned()(3,2), H_enuned2()(3,2), 1e-12);
    check("H_enuned  (3,3):     ", H_enuned()(3,3), H_enuned2()(3,3), 1e-12);

    ang::homogeneous H_ecefned = Oeagle_att.get_H_ecefned();
    ang::homogeneous H_ecefned2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_ecefned());

    check("H_ecefned  (0,0):    ", H_ecefned()(0,0), H_ecefned2()(0,0), 1e-12);
    check("H_ecefned  (0,1):    ", H_ecefned()(0,1), H_ecefned2()(0,1), 1e-12);
    check("H_ecefned  (0,2):    ", H_ecefned()(0,2), H_ecefned2()(0,2), 1e-12);
    check("H_ecefned  (0,3):    ", H_ecefned()(0,3), H_ecefned2()(0,3), 1e-6);
    check("H_ecefned  (1,0):    ", H_ecefned()(1,0), H_ecefned2()(1,0), 1e-12);
    check("H_ecefned  (1,1):    ", H_ecefned()(1,1), H_ecefned2()(1,1), 1e-12);
    check("H_ecefned  (1,2):    ", H_ecefned()(1,2), H_ecefned2()(1,2), 1e-12);
    check("H_ecefned  (1,3):    ", H_ecefned()(1,3), H_ecefned2()(1,3), 1e-6);
    check("H_ecefned  (2,0):    ", H_ecefned()(2,0), H_ecefned2()(2,0), 1e-12);
    check("H_ecefned  (2,1):    ", H_ecefned()(2,1), H_ecefned2()(2,1), 1e-12);
    check("H_ecefned  (2,2):    ", H_ecefned()(2,2), H_ecefned2()(2,2), 1e-12);
    check("H_ecefned  (2,3):    ", H_ecefned()(2,3), H_ecefned2()(2,3), 1e-6);
    check("H_ecefned  (3,0):    ", H_ecefned()(3,0), H_ecefned2()(3,0), 1e-12);
    check("H_ecefned  (3,1):    ", H_ecefned()(3,1), H_ecefned2()(3,1), 1e-12);
    check("H_ecefned  (3,2):    ", H_ecefned()(3,2), H_ecefned2()(3,2), 1e-12);
    check("H_ecefned  (3,3):    ", H_ecefned()(3,3), H_ecefned2()(3,3), 1e-12);

    ang::homogeneous H_nedbfs = Oeagle_att.get_H_nedbfs();
    ang::homogeneous H_nedbfs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_nedbfs());

    check("H_nedbfs  (0,0):     ", H_nedbfs()(0,0), H_nedbfs2()(0,0), 1e-12);
    check("H_nedbfs  (0,1):     ", H_nedbfs()(0,1), H_nedbfs2()(0,1), 1e-12);
    check("H_nedbfs  (0,2):     ", H_nedbfs()(0,2), H_nedbfs2()(0,2), 1e-12);
    check("H_nedbfs  (0,3):     ", H_nedbfs()(0,3), H_nedbfs2()(0,3), 1e-6);
    check("H_nedbfs  (1,0):     ", H_nedbfs()(1,0), H_nedbfs2()(1,0), 1e-12);
    check("H_nedbfs  (1,1):     ", H_nedbfs()(1,1), H_nedbfs2()(1,1), 1e-12);
    check("H_nedbfs  (1,2):     ", H_nedbfs()(1,2), H_nedbfs2()(1,2), 1e-12);
    check("H_nedbfs  (1,3):     ", H_nedbfs()(1,3), H_nedbfs2()(1,3), 1e-6);
    check("H_nedbfs  (2,0):     ", H_nedbfs()(2,0), H_nedbfs2()(2,0), 1e-12);
    check("H_nedbfs  (2,1):     ", H_nedbfs()(2,1), H_nedbfs2()(2,1), 1e-12);
    check("H_nedbfs  (2,2):     ", H_nedbfs()(2,2), H_nedbfs2()(2,2), 1e-12);
    check("H_nedbfs  (2,3):     ", H_nedbfs()(2,3), H_nedbfs2()(2,3), 1e-6);
    check("H_nedbfs  (3,0):     ", H_nedbfs()(3,0), H_nedbfs2()(3,0), 1e-12);
    check("H_nedbfs  (3,1):     ", H_nedbfs()(3,1), H_nedbfs2()(3,1), 1e-12);
    check("H_nedbfs  (3,2):     ", H_nedbfs()(3,2), H_nedbfs2()(3,2), 1e-12);
    check("H_nedbfs  (3,3):     ", H_nedbfs()(3,3), H_nedbfs2()(3,3), 1e-12);

    ang::homogeneous H_bfscrs = Oeagle_att.get_H_bfscrs();
    ang::homogeneous H_bfscrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_bfscrs());

    check("H_bfscrs  (0,0):     ", H_bfscrs()(0,0), H_bfscrs2()(0,0), 1e-12);
    check("H_bfscrs  (0,1):     ", H_bfscrs()(0,1), H_bfscrs2()(0,1), 1e-12);
    check("H_bfscrs  (0,2):     ", H_bfscrs()(0,2), H_bfscrs2()(0,2), 1e-12);
    check("H_bfscrs  (0,3):     ", H_bfscrs()(0,3), H_bfscrs2()(0,3), 1e-6);
    check("H_bfscrs  (1,0):     ", H_bfscrs()(1,0), H_bfscrs2()(1,0), 1e-12);
    check("H_bfscrs  (1,1):     ", H_bfscrs()(1,1), H_bfscrs2()(1,1), 1e-12);
    check("H_bfscrs  (1,2):     ", H_bfscrs()(1,2), H_bfscrs2()(1,2), 1e-12);
    check("H_bfscrs  (1,3):     ", H_bfscrs()(1,3), H_bfscrs2()(1,3), 1e-6);
    check("H_bfscrs  (2,0):     ", H_bfscrs()(2,0), H_bfscrs2()(2,0), 1e-12);
    check("H_bfscrs  (2,1):     ", H_bfscrs()(2,1), H_bfscrs2()(2,1), 1e-12);
    check("H_bfscrs  (2,2):     ", H_bfscrs()(2,2), H_bfscrs2()(2,2), 1e-12);
    check("H_bfscrs  (2,3):     ", H_bfscrs()(2,3), H_bfscrs2()(2,3), 1e-6);
    check("H_bfscrs  (3,0):     ", H_bfscrs()(3,0), H_bfscrs2()(3,0), 1e-12);
    check("H_bfscrs  (3,1):     ", H_bfscrs()(3,1), H_bfscrs2()(3,1), 1e-12);
    check("H_bfscrs  (3,2):     ", H_bfscrs()(3,2), H_bfscrs2()(3,2), 1e-12);
    check("H_bfscrs  (3,3):     ", H_bfscrs()(3,3), H_bfscrs2()(3,3), 1e-12);

    ang::homogeneous H_nedcrs = Oeagle_att.get_H_nedcrs();
    ang::homogeneous H_nedcrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_nedcrs());

    check("H_nedcrs  (0,0):     ", H_nedcrs()(0,0), H_nedcrs2()(0,0), 1e-12);
    check("H_nedcrs  (0,1):     ", H_nedcrs()(0,1), H_nedcrs2()(0,1), 1e-12);
    check("H_nedcrs  (0,2):     ", H_nedcrs()(0,2), H_nedcrs2()(0,2), 1e-12);
    check("H_nedcrs  (0,3):     ", H_nedcrs()(0,3), H_nedcrs2()(0,3), 1e-6);
    check("H_nedcrs  (1,0):     ", H_nedcrs()(1,0), H_nedcrs2()(1,0), 1e-12);
    check("H_nedcrs  (1,1):     ", H_nedcrs()(1,1), H_nedcrs2()(1,1), 1e-12);
    check("H_nedcrs  (1,2):     ", H_nedcrs()(1,2), H_nedcrs2()(1,2), 1e-12);
    check("H_nedcrs  (1,3):     ", H_nedcrs()(1,3), H_nedcrs2()(1,3), 1e-6);
    check("H_nedcrs  (2,0):     ", H_nedcrs()(2,0), H_nedcrs2()(2,0), 1e-12);
    check("H_nedcrs  (2,1):     ", H_nedcrs()(2,1), H_nedcrs2()(2,1), 1e-12);
    check("H_nedcrs  (2,2):     ", H_nedcrs()(2,2), H_nedcrs2()(2,2), 1e-12);
    check("H_nedcrs  (2,3):     ", H_nedcrs()(2,3), H_nedcrs2()(2,3), 1e-6);
    check("H_nedcrs  (3,0):     ", H_nedcrs()(3,0), H_nedcrs2()(3,0), 1e-12);
    check("H_nedcrs  (3,1):     ", H_nedcrs()(3,1), H_nedcrs2()(3,1), 1e-12);
    check("H_nedcrs  (3,2):     ", H_nedcrs()(3,2), H_nedcrs2()(3,2), 1e-12);
    check("H_nedcrs  (3,3):     ", H_nedcrs()(3,3), H_nedcrs2()(3,3), 1e-12);

    ang::homogeneous H_ecefcrs = Oeagle_att.get_H_ecefcrs();
    ang::homogeneous H_ecefcrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_ecefcrs());

    check("H_ecefcrs  (0,0):    ", H_ecefcrs()(0,0), H_ecefcrs2()(0,0), 1e-12);
    check("H_ecefcrs  (0,1):    ", H_ecefcrs()(0,1), H_ecefcrs2()(0,1), 1e-12);
    check("H_ecefcrs  (0,2):    ", H_ecefcrs()(0,2), H_ecefcrs2()(0,2), 1e-12);
    check("H_ecefcrs  (0,3):    ", H_ecefcrs()(0,3), H_ecefcrs2()(0,3), 1e-6);
    check("H_ecefcrs  (1,0):    ", H_ecefcrs()(1,0), H_ecefcrs2()(1,0), 1e-12);
    check("H_ecefcrs  (1,1):    ", H_ecefcrs()(1,1), H_ecefcrs2()(1,1), 1e-12);
    check("H_ecefcrs  (1,2):    ", H_ecefcrs()(1,2), H_ecefcrs2()(1,2), 1e-12);
    check("H_ecefcrs  (1,3):    ", H_ecefcrs()(1,3), H_ecefcrs2()(1,3), 1e-6);
    check("H_ecefcrs  (2,0):    ", H_ecefcrs()(2,0), H_ecefcrs2()(2,0), 1e-12);
    check("H_ecefcrs  (2,1):    ", H_ecefcrs()(2,1), H_ecefcrs2()(2,1), 1e-12);
    check("H_ecefcrs  (2,2):    ", H_ecefcrs()(2,2), H_ecefcrs2()(2,2), 1e-12);
    check("H_ecefcrs  (2,3):    ", H_ecefcrs()(2,3), H_ecefcrs2()(2,3), 1e-6);
    check("H_ecefcrs  (3,0):    ", H_ecefcrs()(3,0), H_ecefcrs2()(3,0), 1e-12);
    check("H_ecefcrs  (3,1):    ", H_ecefcrs()(3,1), H_ecefcrs2()(3,1), 1e-12);
    check("H_ecefcrs  (3,2):    ", H_ecefcrs()(3,2), H_ecefcrs2()(3,2), 1e-12);
    check("H_ecefcrs  (3,3):    ", H_ecefcrs()(3,3), H_ecefcrs2()(3,3), 1e-12);

    ang::homogeneous H_crsmrs = Oeagle_att.get_H_crsmrs();
    ang::homogeneous H_crsmrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_crsmrs());

    check("H_crsmrs  (0,0):     ", H_crsmrs()(0,0), H_crsmrs2()(0,0), 1e-12);
    check("H_crsmrs  (0,1):     ", H_crsmrs()(0,1), H_crsmrs2()(0,1), 1e-12);
    check("H_crsmrs  (0,2):     ", H_crsmrs()(0,2), H_crsmrs2()(0,2), 1e-12);
    check("H_crsmrs  (0,3):     ", H_crsmrs()(0,3), H_crsmrs2()(0,3), 1e-6);
    check("H_crsmrs  (1,0):     ", H_crsmrs()(1,0), H_crsmrs2()(1,0), 1e-12);
    check("H_crsmrs  (1,1):     ", H_crsmrs()(1,1), H_crsmrs2()(1,1), 1e-12);
    check("H_crsmrs  (1,2):     ", H_crsmrs()(1,2), H_crsmrs2()(1,2), 1e-12);
    check("H_crsmrs  (1,3):     ", H_crsmrs()(1,3), H_crsmrs2()(1,3), 1e-12);
    check("H_crsmrs  (2,0):     ", H_crsmrs()(2,0), H_crsmrs2()(2,0), 1e-12);
    check("H_crsmrs  (2,1):     ", H_crsmrs()(2,1), H_crsmrs2()(2,1), 1e-12);
    check("H_crsmrs  (2,2):     ", H_crsmrs()(2,2), H_crsmrs2()(2,2), 1e-12);
    check("H_crsmrs  (2,3):     ", H_crsmrs()(2,3), H_crsmrs2()(2,3), 1e-6);
    check("H_crsmrs  (3,0):     ", H_crsmrs()(3,0), H_crsmrs2()(3,0), 1e-12);
    check("H_crsmrs  (3,1):     ", H_crsmrs()(3,1), H_crsmrs2()(3,1), 1e-12);
    check("H_crsmrs  (3,2):     ", H_crsmrs()(3,2), H_crsmrs2()(3,2), 1e-12);
    check("H_crsmrs  (3,3):     ", H_crsmrs()(3,3), H_crsmrs2()(3,3), 1e-12);

    ang::homogeneous H_ecefmrs = Oeagle_att.get_H_ecefmrs();
    ang::homogeneous H_ecefmrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_ecefmrs());

    check("H_ecefmrs  (0,0):    ", H_ecefmrs()(0,0), H_ecefmrs2()(0,0), 1e-12);
    check("H_ecefmrs  (0,1):    ", H_ecefmrs()(0,1), H_ecefmrs2()(0,1), 1e-12);
    check("H_ecefmrs  (0,2):    ", H_ecefmrs()(0,2), H_ecefmrs2()(0,2), 1e-12);
    check("H_ecefmrs  (0,3):    ", H_ecefmrs()(0,3), H_ecefmrs2()(0,3), 1e-6);
    check("H_ecefmrs  (1,0):    ", H_ecefmrs()(1,0), H_ecefmrs2()(1,0), 1e-12);
    check("H_ecefmrs  (1,1):    ", H_ecefmrs()(1,1), H_ecefmrs2()(1,1), 1e-12);
    check("H_ecefmrs  (1,2):    ", H_ecefmrs()(1,2), H_ecefmrs2()(1,2), 1e-12);
    check("H_ecefmrs  (1,3):    ", H_ecefmrs()(1,3), H_ecefmrs2()(1,3), 1e-6);
    check("H_ecefmrs  (2,0):    ", H_ecefmrs()(2,0), H_ecefmrs2()(2,0), 1e-12);
    check("H_ecefmrs  (2,1):    ", H_ecefmrs()(2,1), H_ecefmrs2()(2,1), 1e-12);
    check("H_ecefmrs  (2,2):    ", H_ecefmrs()(2,2), H_ecefmrs2()(2,2), 1e-12);
    check("H_ecefmrs  (2,3):    ", H_ecefmrs()(2,3), H_ecefmrs2()(2,3), 1e-6);
    check("H_ecefmrs  (3,0):    ", H_ecefmrs()(3,0), H_ecefmrs2()(3,0), 1e-12);
    check("H_ecefmrs  (3,1):    ", H_ecefmrs()(3,1), H_ecefmrs2()(3,1), 1e-12);
    check("H_ecefmrs  (3,2):    ", H_ecefmrs()(3,2), H_ecefmrs2()(3,2), 1e-12);
    check("H_ecefmrs  (3,3):    ", H_ecefmrs()(3,3), H_ecefmrs2()(3,3), 1e-12);

   // ang::homogeneous H_mrsrrs = Oeagle_att.get_H_mrsrrs();
   // ang::homogeneous H_mrsrrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_mrsrrs());

    //check("H_mrsrrs  (0,0):     ", H_mrsrrs(0,0), H_mrsrrs2(0,0), 1e-12);
    //check("H_mrsrrs  (0,1):     ", H_mrsrrs(0,1), H_mrsrrs2(0,1), 1e-12);
    //check("H_mrsrrs  (0,2):     ", H_mrsrrs(0,2), H_mrsrrs2(0,2), 1e-12);
    //check("H_mrsrrs  (0,3):     ", H_mrsrrs(0,3), H_mrsrrs2(0,3), 1e-6);
    //check("H_mrsrrs  (1,0):     ", H_mrsrrs(1,0), H_mrsrrs2(1,0), 1e-12);
    //check("H_mrsrrs  (1,1):     ", H_mrsrrs(1,1), H_mrsrrs2(1,1), 1e-12);
    //check("H_mrsrrs  (1,2):     ", H_mrsrrs(1,2), H_mrsrrs2(1,2), 1e-12);
    //check("H_mrsrrs  (1,3):     ", H_mrsrrs(1,3), H_mrsrrs2(1,3), 1e-6);
    //check("H_mrsrrs  (2,0):     ", H_mrsrrs(2,0), H_mrsrrs2(2,0), 1e-12);
    //check("H_mrsrrs  (2,1):     ", H_mrsrrs(2,1), H_mrsrrs2(2,1), 1e-12);
    //check("H_mrsrrs  (2,2):     ", H_mrsrrs(2,2), H_mrsrrs2(2,2), 1e-12);
    //check("H_mrsrrs  (2,3):     ", H_mrsrrs(2,3), H_mrsrrs2(2,3), 1e-6);
    //check("H_mrsrrs  (3,0):     ", H_mrsrrs(3,0), H_mrsrrs2(3,0), 1e-12);
    //check("H_mrsrrs  (3,1):     ", H_mrsrrs(3,1), H_mrsrrs2(3,1), 1e-12);
    //check("H_mrsrrs  (3,2):     ", H_mrsrrs(3,2), H_mrsrrs2(3,2), 1e-12);
    //check("H_mrsrrs  (3,3):     ", H_mrsrrs(3,3), H_mrsrrs2(3,3), 1e-12);

    ang::homogeneous H_ecefrrs = Oeagle_att.get_H_ecefrrs();
    ang::homogeneous H_ecefrrs2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_ecefrrs());

    check("H_ecefrrs  (0,0):    ", H_ecefrrs()(0,0), H_ecefrrs2()(0,0), 1e-12);
    check("H_ecefrrs  (0,1):    ", H_ecefrrs()(0,1), H_ecefrrs2()(0,1), 1e-12);
    check("H_ecefrrs  (0,2):    ", H_ecefrrs()(0,2), H_ecefrrs2()(0,2), 1e-12);
    check("H_ecefrrs  (0,3):    ", H_ecefrrs()(0,3), H_ecefrrs2()(0,3), 1e-6);
    check("H_ecefrrs  (1,0):    ", H_ecefrrs()(1,0), H_ecefrrs2()(1,0), 1e-12);
    check("H_ecefrrs  (1,1):    ", H_ecefrrs()(1,1), H_ecefrrs2()(1,1), 1e-12);
    check("H_ecefrrs  (1,2):    ", H_ecefrrs()(1,2), H_ecefrrs2()(1,2), 1e-12);
    check("H_ecefrrs  (1,3):    ", H_ecefrrs()(1,3), H_ecefrrs2()(1,3), 1e-6);
    check("H_ecefrrs  (2,0):    ", H_ecefrrs()(2,0), H_ecefrrs2()(2,0), 1e-12);
    check("H_ecefrrs  (2,1):    ", H_ecefrrs()(2,1), H_ecefrrs2()(2,1), 1e-12);
    check("H_ecefrrs  (2,2):    ", H_ecefrrs()(2,2), H_ecefrrs2()(2,2), 1e-12);
    check("H_ecefrrs  (2,3):    ", H_ecefrrs()(2,3), H_ecefrrs2()(2,3), 1e-6);
    check("H_ecefrrs  (3,0):    ", H_ecefrrs()(3,0), H_ecefrrs2()(3,0), 1e-12);
    check("H_ecefrrs  (3,1):    ", H_ecefrrs()(3,1), H_ecefrrs2()(3,1), 1e-12);
    check("H_ecefrrs  (3,2):    ", H_ecefrrs()(3,2), H_ecefrrs2()(3,2), 1e-12);
    check("H_ecefrrs  (3,3):    ", H_ecefrrs()(3,3), H_ecefrrs2()(3,3), 1e-12);

    ang::homogeneous H_rrsecef = Oeagle_att.get_H_rrsecef();
    ang::homogeneous H_rrsecef2 = eagle::earth_eagle::convert_SE3_osg_2_att(Oeagle_osg.get_S_rrsecef());

    check("H_rrsecef  (0,0):    ", H_rrsecef()(0,0), H_rrsecef2()(0,0), 1e-12);
    check("H_rrsecef  (0,1):    ", H_rrsecef()(0,1), H_rrsecef2()(0,1), 1e-12);
    check("H_rrsecef  (0,2):    ", H_rrsecef()(0,2), H_rrsecef2()(0,2), 1e-12);
    check("H_rrsecef  (0,3):    ", H_rrsecef()(0,3), H_rrsecef2()(0,3), 1e-6);
    check("H_rrsecef  (1,0):    ", H_rrsecef()(1,0), H_rrsecef2()(1,0), 1e-12);
    check("H_rrsecef  (1,1):    ", H_rrsecef()(1,1), H_rrsecef2()(1,1), 1e-12);
    check("H_rrsecef  (1,2):    ", H_rrsecef()(1,2), H_rrsecef2()(1,2), 1e-12);
    check("H_rrsecef  (1,3):    ", H_rrsecef()(1,3), H_rrsecef2()(1,3), 1e-6);
    check("H_rrsecef  (2,0):    ", H_rrsecef()(2,0), H_rrsecef2()(2,0), 1e-12);
    check("H_rrsecef  (2,1):    ", H_rrsecef()(2,1), H_rrsecef2()(2,1), 1e-12);
    check("H_rrsecef  (2,2):    ", H_rrsecef()(2,2), H_rrsecef2()(2,2), 1e-12);
    check("H_rrsecef  (2,3):    ", H_rrsecef()(2,3), H_rrsecef2()(2,3), 1e-6);
    check("H_rrsecef  (3,0):    ", H_rrsecef()(3,0), H_rrsecef2()(3,0), 1e-12);
    check("H_rrsecef  (3,1):    ", H_rrsecef()(3,1), H_rrsecef2()(3,1), 1e-12);
    check("H_rrsecef  (3,2):    ", H_rrsecef()(3,2), H_rrsecef2()(3,2), 1e-12);
    check("H_rrsecef  (3,3):    ", H_rrsecef()(3,3), H_rrsecef2()(3,3), 1e-12);

} // closes test_compare

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////v

void eagle::test::Teagle::test_speed_osg_versus_att(eagle::earth_eagle_att& Oeagle_att, eagle::earth_eagle_osg& Oeagle_osg) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad = 38.0 * d2r;
    double theta_rad = 14.0 * d2r;
    double xi_rad = 9.0 * d2r;
    ang::euler Oeuler_nedcam(psi_rad, theta_rad, xi_rad);

    math::timer Otimer;
    Oeagle_osg.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_osg", 1.0); // ensure cache
    Otimer.start();
    for (unsigned int i = 0; i != 100; i++) {
        Oeagle_osg.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_osg", 1.0);
    }
    Otimer.stop();
    std::cout << "Time required to process 100 images with osg transformations: " << Otimer.measure_nanosec() * 1e-9 << " [sec]." << std::endl;

    Oeagle_att.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_att", 1.0);
    Otimer.start();
    for (unsigned int i = 0; i != 100; i++) {
        Oeagle_att.process_single_input(Ogisu, Oeuler_nedcam, "test_earth_eagle_att", 1.0);
    }
    Otimer.stop();
    std::cout << "Time required to process 100 images with att transformations: " << Otimer.measure_nanosec() * 1e-9 << " [sec]." << std::endl;

    std::cout << std::endl << std::endl;

} // closes test_speed

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Teagle::test_save_to_stream(eagle::earth_eagle& Oeagle1, eagle::earth_eagle& Oeagle2, eagle::earth_eagle& Oeagle3) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad = -18.0 * d2r;
    double theta_rad = -12.0 * d2r;
    double xi_rad = -6.0 * d2r;
    ang::euler Oeuler_nedcam(psi_rad, theta_rad, xi_rad);

    Oeagle1.process_single_input(Ogisu, Oeuler_nedcam, "converter_file", 1.0);
    cv::Mat Omat1;
    Oeagle1.obtain_image(Omat1,"converter_file");
    cv::imshow("Converter file", Omat1);
    cv::waitKey(0);

    Oeagle2.process_single_input(Ogisu, Oeuler_nedcam, "converter_file_stream", 1.0);
    cv::Mat Omat2;
    Oeagle2.obtain_image(Omat2,"converter_file_stream");
    cv::namedWindow("Converter file stream", CV_WINDOW_AUTOSIZE);
    cv::imshow("Converter file stream", Omat2);
    cv::waitKey(0);

    //////////////////////////////////////////////////////////
    // Lo que sigue falla porque cuando en earth_eagle_ang::process_single_input se hace _Peagle->frame(),
    // internamente deberia llamar automaticamente al operador () del converter de turno, pero no se hace.
    // En el ordenador antiguo (antes de Feb 2020) sí lo hacía. Al parecer esto no es relevante para los
    // dos casos anteriores pero sí aquí.
    //////////////////////////////////////////////////////////
    Oeagle3.process_single_input(Ogisu, Oeuler_nedcam, "converter_stream", 1.0);
    cv::Mat Omat3;
    Oeagle3.obtain_image(Omat3,"converter_stream");
    cv::namedWindow("Converter stream", CV_WINDOW_AUTOSIZE);
    cv::imshow("Converter stream", Omat3);
    cv::waitKey(0);

    cv::destroyWindow("Converter file stream");
    cv::destroyWindow("Converter stream");
} // closes test_save_to_stream

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Teagle::test_speed_converter(eagle::earth_eagle& Oeagle1, eagle::earth_eagle& Oeagle2, eagle::earth_eagle& Oeagle3) {
    double d2r = math::constant::D2R();
    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);
    double psi_rad = -18.0 * d2r;
    double theta_rad = -12.0 * d2r;
    double xi_rad = -6.0 * d2r;
    ang::euler Oeuler_nedcam(psi_rad, theta_rad, xi_rad);

    math::timer Otimer;

    Oeagle1.process_single_input(Ogisu, Oeuler_nedcam, "converter_file", 1.0); // ensure cache
    cv::Mat Omat1;
    Otimer.start();
    for (unsigned int i = 0; i != 100; i++) {
        Oeagle1.process_single_input(Ogisu, Oeuler_nedcam, "converter_file", 1.0);
        Oeagle1.obtain_image(Omat1,"converter_file");
    }
    Otimer.stop();
    std::cout << "Time required to process 100 images while also saving RGB file: " << Otimer.measure_nanosec() * 1e-9 << " [sec]." << std::endl;

    Oeagle2.process_single_input(Ogisu, Oeuler_nedcam, "converter_file_stream", 1.0);
    cv::Mat Omat2;
    Otimer.start();
    for (unsigned int i = 0; i != 100; i++) {
        Oeagle2.process_single_input(Ogisu, Oeuler_nedcam, "converter_file_stream", 1.0);
        Oeagle2.obtain_image(Omat2,"converter_file_stream");
    }
    Otimer.stop();
    std::cout << "Time required to process 100 images with also saving RGB file and employing stream: " << Otimer.measure_nanosec() * 1e-9 << " [sec]." << std::endl;

    Oeagle3.process_single_input(Ogisu, Oeuler_nedcam, "converter_stream", 1.0);
    cv::Mat Omat3;
    Otimer.start();
    for (unsigned int i = 0; i != 100; i++) {
        Oeagle3.process_single_input(Ogisu, Oeuler_nedcam, "converter_stream", 1.0);
        Oeagle3.obtain_image(Omat3,"converter_stream");
    }
    Otimer.stop();
    std::cout << "Time required to process 100 images without saving RGB file: " << Otimer.measure_nanosec() * 1e-9 << " [sec]." << std::endl;

    std::cout << std::endl << std::endl;

} // closes test_speed

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void eagle::test::Teagle::test_G_bc(eagle::earth_eagle& Oeagle, const std::string& st_case) {
    double d2r = math::constant::D2R();

    double lambda_rad = -2.13426900;
    double phi_rad = +0.83643580;
    double h_m = +1000.0;

 //   double lambda_rad = -121.747 * d2r;
   // double phi_rad    =  46.847 * d2r; //+46.847 * d2r;;
   // double h_m        = +50000.0;
    env::geodetic_coord Ogisu(lambda_rad, phi_rad, h_m);

    double psi_deg    =  0.0;
    double theta_deg  =  0.0;
    double xi_deg     =  0.0;
    ang::euler Oeuler_nedbfs(psi_deg * d2r, theta_deg * d2r, xi_deg * d2r);

    Oeagle.process_single_input(Ogisu, Oeuler_nedbfs, st_case, 1.0);

}
/* covers the influence of the transformation between body and camera frames */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////



