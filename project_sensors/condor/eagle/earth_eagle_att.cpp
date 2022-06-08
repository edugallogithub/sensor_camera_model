#include "earth_eagle_att.h"
#include "math/logic/share.h"
#include "math/logic/seeder.h"
#include "math/logic/log_writer_mgr.h"
#include "ang/rotate/euler.h"
#include "env/coord.h"
#include "env/geo.h"
#include "acft/sens/camera.h"

#include <boost/filesystem.hpp>

// CLASS EARTH EAGLE_ATT
// =====================
// =====================

eagle::earth_eagle_att::earth_eagle_att(const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                                           eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename, int* Xargc, char **Xargv)
: earth_eagle(Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, Xargc, Xargv, true), _Pgeo(&Ogeo) {

    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    // THIS CLASS IS NO LONGER VALID -- IT NEEDS WORK TO HAVE THE SAME BEHAVIOR AS EARTH_EAGLE_FINAL
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////

    Eigen::Vector3d v3zero = Eigen::Vector3d::Zero();
    ang::euler euler_enuned_rad(math::constant::PIHALF(), 0.0, math::constant::PI());
    ang::dcm R_enuned(euler_enuned_rad);
    _H_enuned.set(R_enuned, v3zero); // ENU-NED (note that by coincidence enu-ned equals ned-enu)
    _H_crsmrs = _H_enuned; // CRS-MRS (modified camera reference system, the same but with different orientation)
    //double rot_rad = 0.; // extra possible 90 [deg] camera rotation. It is zero if _flag_rotate is true.
    //if (flag_rotate == false) {rot_rad = math::constant::PIHALF();}
    //ang::euler euler_mrsrrs_rad(rot_rad, 0.0, 0.0);
    //ang::dcm R_mrsrrs(euler_mrsrrs_rad);
    //_H_mrsrrs.set(R_mrsrrs, v3zero); // MRS-RRS (rotated camera reference system, the same but rotated pi/2 with respect to MRS)

    _H_bfscrs = _Pcam->get_Gq_bc_truth(1.0); /////////////////////////// IT SHOULD BE BC INSTEAD OF RC
    //ang::dcm R_bfscrs(_Pcam->get_euler_bc());
    //_H_bfscrs.set(R_bfscrs, v3zero); // BFS-CRS (camera reference system)

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    // TODO: I need to remove 90 deg from the input angles before the transformation, exactly the same as in the earth_eagle_final constructor
    // With those changes, it should work because it worked before
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    throw std::runtime_error("NOT READY YET");
}
/* constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed,
 * enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
 * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
 * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */

void eagle::earth_eagle_att::process_text_file(const std::string& st_folder, bool flag) {
    if (flag == false) {return;}

    this->empty_images_folder();

    std::string st_case_guid = st_folder.substr(0,2);
    boost::filesystem::path path_outputs_nav(_st_output_nav);
    boost::filesystem::path path_folder1("nav");
    boost::filesystem::path path_folder2(st_case_guid);
    boost::filesystem::path path_folder3(st_folder);
    boost::filesystem::path path_folder_nav(path_outputs_nav / path_folder1 / path_folder2 / path_folder3);
    boost::filesystem::path path_file("images_truth.txt");
    std::string st_input_trj_file = (path_folder_nav / path_file).string();

    math::log_writer_mgr::get_log().write(std::string("Eagle:          Start creating images from case " + st_folder + "."));

    // open input file with trajectory 6 dof info
    std::ifstream Oinput;
    Oinput.open(st_input_trj_file);
    assert(Oinput);

    int id; // frame number
    double t_sec; // frame time
    env::geodetic_coord gisu;
    ang::euler euler_nedbfs;
    double ratio;

    // generate 1st image several times without saving it to avoid fuzzy images
    std::string line_str;
    std::getline(Oinput, line_str);  // read header line

    unsigned int i = 0;
    std::string st_i;
    Oinput >> id >> t_sec >> gisu >> euler_nedbfs >> ratio; // obtain time, geodetic position, and attitude
    do {
        math::log_writer_mgr::get_log().write(std::string("Eagle:          Created image # " + std::to_string(i) + "."));
        math::digits::to_string_5digits(i, st_i); // generate 5 digit string with image number

        this->process_single_input(gisu, euler_nedbfs, st_i, ratio);

        i++; // go to next image

        // do not put this line at start of loop as file not finished and while conditions forces repetition of last image
        Oinput >> id >> t_sec >> gisu >> euler_nedbfs; // obtain time, geodetic position, and attitude
    } while (Oinput.eof() == false);

    Oinput.close();
}
/* processes all images contained in the input text file (containing body frame position and attitude for the different images) based on
 * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
 * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
 * starting with 0. If flag is false, method does nothing. Note that file contains the body pose, not the camera (body to camera transformation
 * automatically added). */

bool eagle::earth_eagle_att::process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) {
    this->set_camera(x_gdt_b_rad_m, euler_nb, ratio); // pass position and attitude to camera
    return this->process_image(st_i);
}
/* generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
 * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
 * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
 * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */

void eagle::earth_eagle_att::set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) {
    // ECEF-NED
    ang::euler euler_ecefned_rad = x_gdt_b_rad_m.obtain_euler_en();
    ang::dcm R_ecefned(euler_ecefned_rad);
    Eigen::Vector3d x_ecefned_ecef = _Pgeo->geodetic2cartesian(x_gdt_b_rad_m, _Pgeo->radius_vert(x_gdt_b_rad_m.get_phi_rad()))();
    _H_ecefned.set(R_ecefned, x_ecefned_ecef);

    // ECEF-ENU
    //_H_ecefenu = _H_ecefned * _H_enuned; // not necessary but do not delete just in case

    // NED-BFS (body fixed system)
    ang::dcm R_nedbfs(euler_nb);
    _H_nedbfs.set(R_nedbfs, Eigen::Vector3d::Zero());

    // NED-CRS (camera reference system)
    _H_nedcrs = _H_nedbfs * _H_bfscrs;

    // ECEF-CRS
    _H_ecefcrs = _H_ecefned * _H_nedcrs;

    // ECEF-MRS
    _H_ecefmrs = _H_ecefcrs * _H_crsmrs;

    // ECEF-RRS
    _H_ecefrrs = _H_ecefmrs; // * _H_mrsrrs;

    // RRS-ECEF
    _H_rrsecef = _H_ecefrrs.inverse();

    _Peagle->getCamera()->setViewMatrix(eagle::earth_eagle::convert_SE3_att_2_osg(_H_rrsecef));
}
/* set eagle body frame pose (takes body frame to camera frame transformation from camera) */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////



