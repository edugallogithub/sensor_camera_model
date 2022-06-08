#include "earth_eagle_final.h"
#include "math/logic/share.h"
#include "math/logic/seeder.h"
#include "math/logic/log_writer_mgr.h"
#include "ang/tools.h"
#include "ang/rotate/euler.h"
#include "ang/transform/speu_rodrigues.h"
#include "ang/transform/homogeneous.h"
#include "env/coord.h"
#include "env/geo.h"
#include "acft/sens/camera.h"

#include <boost/filesystem.hpp>

// CLASS EARTH EAGLE_FINAL
// =======================
// =======================

eagle::earth_eagle_final::earth_eagle_final(const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                                           eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename, int* Xargc, char **Xargv)
: earth_eagle(Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, Xargc, Xargv, true), _Pgeo(&Ogeo) {
    /**< transformation that converts axes (1,2,3) in (2,1,-3), and viceversa (it is symmetric) */
    ang::euler euler_switch_rad(math::constant::PIHALF(), 0.0, math::constant::PI());
    ang::dcm R_switch(euler_switch_rad);
    _H_switch.set(R_switch, Eigen::Vector3d::Zero());
}
/* constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed,
 * enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
 * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
 * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */

void eagle::earth_eagle_final::process_text_file(const std::string& st_folder, bool flag) {
    if (flag == false) {return;}

    this->empty_images_folder();

    std::string st_case_guid = st_folder.substr(0,2);
    boost::filesystem::path path_outputs(_st_output_nav);
    boost::filesystem::path path_folder1("nav");
    boost::filesystem::path path_folder2(st_case_guid);
    boost::filesystem::path path_folder3(st_folder);
    boost::filesystem::path path_folder(path_outputs / path_folder1 / path_folder2 / path_folder3);
    boost::filesystem::path path_file("images_truth.txt");
    std::string st_input_trj_file = (path_folder / path_file).string();

    math::log_writer_mgr::get_log().write(std::string("Eagle:          Start creating images from case " + st_folder + "."));

    // open input file with trajectory 6 dof info
    std::ifstream Oinput;
    Oinput.open(st_input_trj_file);
    assert(Oinput);

    int id; // frame number
    double t_sec; // frame time
    env::geodetic_coord x_gdt_b_rad_m;
    ang::euler euler_nb;
    double ratio; // mass ratio

    // generate 1st image several times without saving it to avoid fuzzy images
    std::string line_str;
    std::getline(Oinput, line_str);  // read header line

    int i = -1;
    std::string st_i;

    while (Oinput.eof() == false) {
        Oinput >> id >> t_sec >> x_gdt_b_rad_m >> euler_nb >> ratio; // obtain time, geodetic position (PHI before LAMBDA), and attitude
        i++; // go to next image
        math::digits::to_string_5digits(i, st_i); // generate 5 digit string with image number
        this->process_single_input(x_gdt_b_rad_m, euler_nb, st_i, ratio);
        math::log_writer_mgr::get_log().write(std::string("Eagle:          Created image # " + st_i + "."));
    }
    Oinput.close();
    math::log_writer_mgr::get_log().write(std::string("Eagle:          Note that image # " + st_i + " has not been created (but message appears)."));
}
/* processes all images contained in the input text file (containing body frame position and attitude for the different images) based on
 * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
 * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
 * starting with 0. If flag is false, method does nothing. Note that file contains the body pose, not the camera (body to camera transformation
 * automatically added). */

void eagle::earth_eagle_final::process_text_file_ROZAS(const std::string& st_folder, bool flag) {
    if (flag == false) {return;}

    this->empty_images_folder();

    boost::filesystem::path p0(_st_output_nav);
    std::string st_file = "image_ground_truth.txt";
    std::string st_input = (p0 / st_file).string();

    std::cout << st_input << std::endl;

    // open input file with trajectory 6 dof info
    std::ifstream Oinput;
    Oinput.open(st_input);
    assert(Oinput);

    int id; // frame number
    double gt_time, ros_time, alo_time;
    double phi_rad, lambda_rad, h_m, q0, q1, q2, q3;

    double t_sec; // frame time
    env::geodetic_coord gisu;
    ang::euler euler_nedbfs;

    std::string line_str;
    std::getline(Oinput, line_str);  // read header line

    std::string st_i;
    while (Oinput.eof() == false) {
        Oinput >> id >> gt_time >> ros_time >> alo_time >> phi_rad >> lambda_rad >> h_m >> q0 >> q1 >> q2 >> q3;

        gisu.set(lambda_rad, phi_rad, h_m);
        ang::rodrigues rodrigues_nedbfs(-q0, q1, q2, q3);
        euler_nedbfs = rodrigues_nedbfs;
        math::digits::to_string_5digits(id, st_i); // generate 5 digit string with image number
        this->process_single_input(gisu, euler_nedbfs, st_i, 1.0); //////////////////////////////////////////////////// TODO hay que leer el ratio, no 1.0
        math::log_writer_mgr::get_log().write(std::string("Eagle:          Created image # " + st_i + "."));
    }
    Oinput.close();
}
/* same as process_text_file but reads different file format generated in real flights */

bool eagle::earth_eagle_final::process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) {
    this->set_camera(x_gdt_b_rad_m, euler_nb, ratio); // pass position and attitude to camera
    return this->process_image(st_i);
}
/* generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
 * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
 * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
 * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */

void eagle::earth_eagle_final::set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) {
    // transformation between Earth Centered Earth Fixed (ECEF - E) and North-East-Down (NED - N)
    ang::euler euler_en_rad = x_gdt_b_rad_m.obtain_euler_en();
    ang::dcm R_en(euler_en_rad);
    Eigen::Vector3d x_ene_m = _Pgeo->geodetic2cartesian(x_gdt_b_rad_m, _Pgeo->radius_vert(x_gdt_b_rad_m.get_phi_rad()))();
    ang::homogeneous H_en(R_en, x_ene_m);

    // transformation between North-East-Down (NED - N) and Body Fixed System (BFS - B)
    ang::dcm R_nb(euler_nb);
    ang::homogeneous H_nb(R_nb, Eigen::Vector3d::Zero());

    // transformation between Body Fixed System (BFS) and Rotated Modified Camera Reference System (RRS)
    // need to remove 90 [deg] of yaw (needs to be done here)
    const Eigen::Vector3d& T_bcb_m = _Pcam->get_Gq_bc_truth(ratio).get_T();
    ang::euler euler_bc(_Pcam->get_Gq_bc_truth(ratio).get_rodrigues());
    euler_bc.set_yaw_rad(euler_bc.get_yaw_rad() - math::constant::PIHALF());
    ang::tools::correct_yaw_rad(euler_bc.get_yaw_rad());
    ang::dcm R_bc(euler_bc); // note that I call it b but it is not as it has been rotated
    ang::homogeneous H_bc(R_bc, T_bcb_m);
    ang::homogeneous H_br = H_bc * _H_switch;

    _Peagle->getCamera()->setViewMatrix(eagle::earth_eagle::convert_SE3_att_2_osg((H_en * H_nb * H_br).inverse()));
}
/* set eagle body frame pose (takes body frame to camera frame transformation from camera) */


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////



