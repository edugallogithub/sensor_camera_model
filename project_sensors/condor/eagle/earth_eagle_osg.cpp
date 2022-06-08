#include "earth_eagle_osg.h"
#include "math/logic/share.h"
#include "math/logic/seeder.h"
#include "math/logic/constant.h"
#include "math/logic/log_writer_mgr.h"
#include "ang/rotate/euler.h"
#include "ang/transform/speu_rodrigues.h"
#include "env/coord.h"
#include "acft/sens/camera.h"

#include <boost/filesystem.hpp>

// CLASS EARTH EAGLE_OSG
// =====================
// =====================

eagle::earth_eagle_osg::earth_eagle_osg(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale,
                                           const std::string& st_images_name, const std::string& st_map_filename, int* Xargc, char **Xargv)
: earth_eagle(Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, Xargc, Xargv, true) {

    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    // THIS CLASS IS NO LONGER VALID -- IT NEEDS WORK TO HAVE THE SAME BEHAVIOR AS EARTH_EAGLE_FINAL
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////

    // constant homogeneous rotations
    _S_enuned.makeRotate(osg::PI, osg::X_AXIS, 0.0, osg::Y_AXIS, osg::PI_2, osg::Z_AXIS); // ENU-NED
    _S_crsmrs = _S_enuned; // CRS-MRS
    //double rot_rad = 0.; // extra possible 90 [deg] camera rotation. It is zero if _flag_rotate is true.
    //if (flag_rotate == false) {rot_rad = math::constant::PIHALF();}
    //_S_mrsrrs.makeRotate(0.0, osg::X_AXIS, 0.0, osg::Y_AXIS, 0.0, osg::Z_AXIS); // MRS-RRS

    ang::euler euler_bc(_Pcam->get_Gq_bc_truth(1.0).get_rodrigues());
    _S_bfscrs.makeRotate(euler_bc.get_bank_rad(), osg::X_AXIS,
                         euler_bc.get_pitch_rad(), osg::Y_AXIS,
                         euler_bc.get_yaw_rad(), osg::Z_AXIS); // BFS-CRS

    //_S_bfscrs.makeRotate(_Pcam->get_euler_bc().get_bank_rad(), osg::X_AXIS,
    //                     _Pcam->get_euler_bc().get_pitch_rad(), osg::Y_AXIS,
    //                     _Pcam->get_euler_bc().get_yaw_rad(), osg::Z_AXIS); // BFS-CRS

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
/* constructor based on eagle camera, enumerator indicating the image format employed, enumerator indicating
 * file conversion option (from osg to opencv), name of folder where images are to be saved (if applicable), name of file containing
 * the Earth maps, weak pointer to number of arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing
 * arguments in executable (AGAINST RULE) */

void eagle::earth_eagle_osg::process_text_file(const std::string& st_folder, bool flag) {
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

bool eagle::earth_eagle_osg::process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) {
    this->set_camera(x_gdt_b_rad_m, euler_nb, ratio); // pass position and attitude to camera
    return this->process_image(st_i);
}
/* generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
 * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
 * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
 * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */

void eagle::earth_eagle_osg::set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) {
    _Pellipsoid->computeLocalToWorldTransformFromLatLongHeight(x_gdt_b_rad_m.get_phi_rad(), x_gdt_b_rad_m.get_lambda_rad(), x_gdt_b_rad_m.get_h_m(), _S_ecefenu); // ECEF-ENU
    _S_ecefned = _S_enuned * _S_ecefenu; // ECEF-NED
    _S_nedbfs.makeRotate(euler_nb.get_bank_rad(), osg::X_AXIS, euler_nb.get_pitch_rad(), osg::Y_AXIS, euler_nb.get_yaw_rad(), osg::Z_AXIS); // NED-BFS
    _S_nedcrs =  _S_bfscrs * _S_nedbfs; // NED-CRS
    _S_ecefcrs = _S_nedcrs * _S_ecefned; // ECEF-CRS
    _S_ecefmrs = _S_crsmrs * _S_ecefcrs; // ECEF-MRS
    _S_ecefrrs = _S_ecefmrs; //_S_mrsrrs * _S_ecefmrs; // ECEF-RRS
    _S_rrsecef.invert(_S_ecefrrs); // RRS-ECEF
    _Peagle->getCamera()->setViewMatrix(_S_rrsecef);
}
/* set eagle body frame pose (takes body frame to camera frame transformation from camera) */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////



