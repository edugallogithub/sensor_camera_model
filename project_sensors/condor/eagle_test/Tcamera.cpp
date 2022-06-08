#include "Tcamera.h"
#include "ang/rotate/euler.h"
#include "env/coord.h"
#include "env/geo.h"

eagle::test::Tcamera::Tcamera(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void eagle::test::Tcamera::run() {
	::jail::unit_test::run();

    /*
     * This is a simple test but clarifies many issues. Given point coordinates, camera
     * coordinates, and rotation, it gives you how that point is viewed from the camera.
     */

    test_camera();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void eagle::test::Tcamera::test_camera() {
    env::geo_mix Ogeo(env::logic::zone_default);

    env::geodetic_coord xC_gdt_rad_m(-2.08228593, 0.82349135, 580.72165881);
    double NC_m = Ogeo.radius_vert(xC_gdt_rad_m.get_phi_rad());
    env::cartesian_coord xC_car_m = Ogeo.geodetic2cartesian(xC_gdt_rad_m, NC_m);
    Eigen::Vector3d T_ece_m = xC_car_m(); // vector from Earth to Camera viewed in Earth

    ang::euler euler_en(xC_gdt_rad_m.get_lambda_rad(), - xC_gdt_rad_m.get_phi_rad() - math::constant::PIHALF(), 0.0); //TODO ////// CANDIDATE TO ITS OWN METHOD //////////////////////////////////////
    ang::dcm R_en(euler_en); // rotation from Earth to camera NED

    ang::euler euler_nc(-0.30543262, +1.13446401, 0.0); // Camera X looks right, Camera Y looks down, Camera Z looks ahead
    ang::dcm R_nc(euler_nc); // rotation from NED Camera to Camera

    ang::dcm R_ec = R_en * R_nc; // rotation from Earth to Camera

    env::geodetic_coord x1_gdt_rad_m(-2.08235622, 0.82362685, 355.24871982);
    double N1_m = Ogeo.radius_vert(x1_gdt_rad_m.get_phi_rad());
    env::cartesian_coord x1_car_m = Ogeo.geodetic2cartesian(x1_gdt_rad_m, N1_m);
    Eigen::Vector3d T_e1e_m = x1_car_m(); // vector from Earth to 1 viewed in Earth

    env::geodetic_coord x2_gdt_rad_m(-2.08234285, 0.82362978, 355.24871982);
    double N2_m = Ogeo.radius_vert(x2_gdt_rad_m.get_phi_rad());
    env::cartesian_coord x2_car_m = Ogeo.geodetic2cartesian(x2_gdt_rad_m, N2_m);
    Eigen::Vector3d T_e2e_m = x2_car_m(); // vector from Earth to 2 viewed in Earth

    env::geodetic_coord x3_gdt_rad_m(-2.08263378, 0.82424469, 355.55352353);
    double N3_m = Ogeo.radius_vert(x3_gdt_rad_m.get_phi_rad());
    env::cartesian_coord x3_car_m = Ogeo.geodetic2cartesian(x3_gdt_rad_m, N3_m);
    Eigen::Vector3d T_e3e_m = x3_car_m(); // vector from Earth to 3 viewed in Earth

    env::geodetic_coord x4_gdt_rad_m(-2.08264716, 0.82424175, 355.55352353);
    double N4_m = Ogeo.radius_vert(x4_gdt_rad_m.get_phi_rad());
    env::cartesian_coord x4_car_m = Ogeo.geodetic2cartesian(x4_gdt_rad_m, N4_m);
    Eigen::Vector3d T_e4e_m = x4_car_m(); // vector from Earth to 4 viewed in Earth

    Eigen::Vector3d T_c1e_m = T_e1e_m - T_ece_m; // vector from Camera to 1 viewed in Earth
    Eigen::Vector3d T_c2e_m = T_e2e_m - T_ece_m; // vector from Camera to 2 viewed in Earth
    Eigen::Vector3d T_c3e_m = T_e3e_m - T_ece_m; // vector from Camera to 3 viewed in Earth
    Eigen::Vector3d T_c4e_m = T_e4e_m - T_ece_m; // vector from Camera to 4 viewed in Earth

    Eigen::Vector3d T_c1c_m = R_ec / T_c1e_m; // vector from Camera to 1 viewed in Camera
    Eigen::Vector3d T_c2c_m = R_ec / T_c2e_m; // vector from Camera to 2 viewed in Camera
    Eigen::Vector3d T_c3c_m = R_ec / T_c3e_m; // vector from Camera to 3 viewed in Camera
    Eigen::Vector3d T_c4c_m = R_ec / T_c4e_m; // vector from Camera to 4 viewed in Camera

    std::cout << "T_e1e_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_e1e_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_e1e_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_e1e_m(2) << std::endl;
    std::cout << "T_ece_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_ece_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_ece_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_ece_m(2) << std::endl;
    std::cout << "T_c1e_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1e_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1e_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1e_m(2) << std::endl;
    std::cout << "T_c1c_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1c_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1c_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1c_m(2) << std::endl;
    std::cout << "T_c2c_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c2c_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c2c_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c2c_m(2) << std::endl;
    std::cout << "T_c3c_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c3c_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c3c_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c3c_m(2) << std::endl;
    std::cout << "T_c4c_m   " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c4c_m(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c4c_m(1)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c4c_m(2) << std::endl;

    // Datos: Image 1024 x 768 px, vfov = 49.2255 deg, 17e-6 m/px --> f = 19 mm = 1117.647 px
    double f_px = 1117.647;

    Eigen::Vector2d T_c1c_px; T_c1c_px << T_c1c_m(0) * f_px / T_c1c_m(2), T_c1c_m(1) * f_px / T_c1c_m(2);
    Eigen::Vector2d T_c2c_px; T_c2c_px << T_c2c_m(0) * f_px / T_c2c_m(2), T_c2c_m(1) * f_px / T_c2c_m(2);
    Eigen::Vector2d T_c3c_px; T_c3c_px << T_c3c_m(0) * f_px / T_c3c_m(2), T_c3c_m(1) * f_px / T_c3c_m(2);
    Eigen::Vector2d T_c4c_px; T_c4c_px << T_c4c_m(0) * f_px / T_c4c_m(2), T_c4c_m(1) * f_px / T_c4c_m(2);

    std::cout << "T_c1c_px  " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1c_px(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c1c_px(1) << std::endl;
    std::cout << "T_c2c_px  " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c2c_px(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c2c_px(1) << std::endl;
    std::cout << "T_c3c_px  " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c3c_px(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c3c_px(1) << std::endl;
    std::cout << "T_c4c_px  " << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c4c_px(0)
                              << std::fixed << std::showpos << std::setw(15) << std::setprecision(3) << T_c4c_px(1) << std::endl;


} // closes test_camera

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



