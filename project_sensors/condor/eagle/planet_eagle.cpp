#include "planet_eagle.h"
#include "math/logic/share.h"
#include "math/logic/seeder.h"
#include "math/logic/log_writer_mgr.h"
#include "acft/sens/camera.h"

#include <boost/filesystem.hpp>
#include <osgViewer/Viewer>
#include <osgEarthAnnotation/FeatureNode>

/*
 * THIS CLASS WORKS AS INTENDED BUT IS OF NO USE TO THIS PROJECT, AS IT WAS INTENDED FOR CLARIFICATION PURPOSES ONLY.
 * AMAZINGLY, IF THE BELOW CODE IS DECOMMENTED, THE PROCESS TEXT FILE FUNCTION NOT HERE BUT IN THE EARTH EAGLE FINAL
 * CLASS ALWAYS CRASHES RIGHT BEFORE THE END WITHOUT GENERATING THE LAST IMAGE. SO KEEP IT COMMENTED FOR HISTORIC REASONS,
 * OR TAKE IT AWAY, NOT REALLY THAT IMPORTANT.
 */

// CLASS PLANET VIEWER
// ===================
// ===================

eagle::planet_eagle::planet_eagle(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, bool flag_rotate,
                                     eagle::logic::CONVERTER_ID converter_id, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename, int* Xargc, char **Xargv)
: _image_id(image_id), _Pcam(&Ocam), _LOD_scale(LOD_scale) {

//    if (flag_rotate == true) {
//        throw std::runtime_error("flag_rotate needs to be false.");
//    }
//
//    boost::filesystem::path path_eagle("eagle";
//    boost::filesystem::path path_maps("maps");
//    boost::filesystem::path path_earth_files("earth_files");
//
//    // complete path of each image excluding number and extension */
//    boost::filesystem::path path_outputs_folder(math::share::get_location(output_folder));
//    boost::filesystem::path path_images_filename(st_images_name);
//    boost::filesystem::path path_images_filename_(st_images_name + "_");
//    boost::filesystem::create_directory(path_outputs_folder / path_eagle/ path_images_filename);
//    _st_image_name_file = (path_outputs_folder / path_eagle / path_images_filename / path_images_filename_).string();
//
//    // complete path of earth maps
//    boost::filesystem::path path_config(math::share::condor_inputs_prefix);
//    boost::filesystem::path path_map_filename(st_map_filename);
//    std::string st_earth_filepath = (path_config / path_eagle / path_maps / path_earth_files / path_map_filename).string();
//
//    // arguments to be passed to osgViewer Viewer
//    // pointer to integer containing number of arguments passed to osgViewer Viewer, which are the same of the executable plus one added here
//    int* argc = Xargc;
//    // pointer to (array of C-string pointers), which is the same as weak pointer to (pointer to C-string pointers) containing
//    // the arguments passed to osgViewer Viewer, which are the same same of the executable plus one added here.
//    char** argv = Xargv;
//    argv[1] = new char[st_earth_filepath.size() + 1];
//    std::strcpy(argv[1], st_earth_filepath.c_str());
//    ++(*argc);
//    osg::ArgumentParser Oarguments(argc, argv); // arguments for eagle initialization
//
//    // eagle initialized with arguments so it can capture path of earth map
//    _Peagle = new osgViewer::Viewer(Oarguments);
//
//    // add Earth map to eagle
//    _Pnode = osgDB::readNodeFiles(Oarguments);
//    _Peagle->setSceneData(_Pnode.get());
//
//    // file extension
//    switch (image_id) {
//        case eagle::logic::image_jpg:
//            _st_extensiondot = ".jpg";
//            _st_extension    = "jpg";
//            break;
//        case eagle::logic::image_png:
//            _st_extensiondot = ".png";
//            _st_extension    = "png";
//            break;
//        case eagle::logic::image_bmp:
//            _st_extensiondot = ".bmp";
//            _st_extension    = "bmp";
//            break;
//        default:
//            throw std::runtime_error("Invalid image format option.");
//    }
//
//    // converter (from osg to opencv)
//    switch (converter_id) {
//        case eagle::logic::converter_osg_vis_opencv_no:
//            _Pconverter = new eagle::converter_osg_vis_opencv_no();
//            break;
//        case eagle::logic::converter_osg_sav_opencv_no:
//            _Pconverter = new eagle::converter_osg_sav_opencv_no();
//            break;
//        case eagle::logic::converter_osg_sav_opencv_vis:
//            _Pconverter = new eagle::converter_osg_sav_opencv_vis();
//            break;
//        case eagle::logic::converter_osg_sav_opencv_sav:
//            _Pconverter = new eagle::converter_osg_sav_opencv_sav();
//            break;
//        case eagle::logic::converter_osg_sav_opencv_vis_stream:
//            _Pconverter = new eagle::converter_osg_sav_opencv_vis_stream();
//            break;
//        case eagle::logic::converter_osg_sav_opencv_sav_stream:
//            _Pconverter = new eagle::converter_osg_sav_opencv_sav_stream();
//            break;
//        case eagle::logic::converter_osg_vis_opencv_vis:
//            _Pconverter = new eagle::converter_osg_vis_opencv_vis(_st_extension);
//            break;
//        case eagle::logic::converter_osg_vis_opencv_sav:
//            _Pconverter = new eagle::converter_osg_vis_opencv_sav(_st_extension);
//            break;
//        default:
//            throw std::runtime_error("Invalid converter option.");
//    }
//
//    // add screen capture to eagle
//    _Pcapture = new osgViewer::ScreenCaptureHandler(_Pconverter);
//    _Peagle->addEventHandler(_Pcapture);
//
//    // add manipulator to eagle
//    _Pmanip = new osgEarth::Util::EarthManipulator();
//    _Peagle->setCameraManipulator(_Pmanip.get());
//
//    // new set projection matrix (????????)
//    double aspect_ratio = (double) (_Pcam->get_width_px()) / (double) (_Pcam->get_height_px());
//    double near_z = 1;
//    double far_z = 10000;
//    _Peagle->getCamera()->setProjectionMatrixAsPerspective(_Pcam->get_vfov_deg(), aspect_ratio, near_z, far_z);
//    _Peagle->getCamera()->setProjectionResizePolicy(osg::Camera::FIXED);
//    _Peagle->getCamera()->setLODScale(_LOD_scale);
//    _Peagle->setUpViewInWindow(0, 0, _Pcam->get_width_px(), _Pcam->get_height_px()); // deprecated but replacement below does not work
//    //_Peagle->apply(new osgViewer::SingleWindow(0, 0, width_px, height_px));
//
//    // realize eagle
//    _Peagle->realize();

}
/* constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed, flag indicating
 * camera rotation, enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
 * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
 * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */

void eagle::planet_eagle::process_text_file(const std::string& st_folder) {
//    std::string st_case_guid = st_folder.substr(0,2);
//    boost::filesystem::path path_outputs(math::share::get_location(output_folder));
//    boost::filesystem::path path_folder1("nav");
//    boost::filesystem::path path_folder2(st_case_guid);
//    boost::filesystem::path path_folder3(st_folder);
//    boost::filesystem::path path_folder(path_outputs / path_folder1 / path_folder2 / path_folder3);
//    boost::filesystem::path path_file("images_truth.txt");
//    std::string st_input_trj_file = (path_folder / path_file).string();
//
//    math::log_writer_mgr::get_log().write(std::string("Eagle:          Start creating images from case " + st_folder + "."));
//
//    // open input file with trajectory 6 dof info
//    std::ifstream Oinput;
//    Oinput.open(st_input_trj_file);
//    assert(Oinput);
//
//    int id; // frame number
//    double t_sec; // frame time
//    env::geodetic_coord gisu;
//    ang::euler euler_nedcam;
//
//    // generate 1st image several times without saving it to avoid fuzzy images
//    std::string line_str;
//    std::getline(Oinput, line_str);  // read header line
//
//    unsigned int i = 0;
//    std::string st_i;
//    Oinput >> id >> t_sec >> gisu >> euler_nedcam; // obtain time, geodetic position, and attitude
//    do {
//        math::log_writer_mgr::get_log().write(std::string("Eagle:          Created image # " + std::to_string(i) + "."));
//        math::digits::to_string_5digits(i, st_i); // generate 5 digit string with image number
//
//        this->process_single_input(gisu, euler_nedcam.get_yaw_rad(), euler_nedcam.get_pitch_rad(), st_i);
//
//        i++; // go to next image
//
//        // do not put this line at start of loop as file not finished and while conditions forces repetition of last image
//        Oinput >> id >> t_sec >> gisu >> euler_nedcam; // obtain time, geodetic position, and attitude
//    } while (Oinput.eof() == false);
//
//    Oinput.close();
}
/* processes all images contained in the input text file (containing camera position and attitude for the different images) based on
 * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
 * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images. */

void eagle::planet_eagle::process_single_input(const env::geodetic_coord& Ogisu, const double& psi_rad, const double& theta_rad, const std::string& st_i) {
//    this->set_camera(Ogisu, psi_rad, theta_rad); // pass position and attitude to camera
//    this->process_image(st_i);
}
/* generates an OpenSceneGraph image based on the input camera position and attitude, visualizing the image in RGB and then following
 * the converter policy, either saving it to RGB file or stream). */

void eagle::planet_eagle::process_image(const std::string& st_i) {
//    _Peagle->frame(); // call frame once so the viewpoint is rendered at least once
//
//    _st_filename_mono = _st_image_name_file + st_i + _st_extensiondot; // generate complete image name
//    _st_filename_rgb  = _st_image_name_file + st_i + "_rgb" +  _st_extensiondot; // generate complete image name
//
//    // Wait for all requests to complete.  Not doing this is why some of your screenshots are lower resolution.
//    // This could take a long time depending on the scene, so if you have time constraints you could wait until the
//    // requests are no longer in progress or some time threshold has been reached if you'd like.
//    while (_Peagle->getDatabasePager()->getRequestsInProgress() == true) {
//        //OE_NOTICE << "Waiting on " << _st_filename_rgb << std::endl;
//        _Peagle->frame();
//    }
//
//    // capture the frame
//    _Pcapture->setFramesToCapture(1);
//    _Pcapture->captureNextFrame(*_Peagle);
//    _Pconverter->set_filename(_st_filename_rgb, _st_filename_mono, _image_id);
//    _Peagle->frame(); // saves the image
//    _Pcapture->stopCapture();
}
/* internal processing of image previously loaded into viewer. Inputs are image number. */

void eagle::planet_eagle::set_camera(const env::geodetic_coord& Ogisu, const double& psi_rad, const double& theta_rad) {
//    if (theta_rad < 0) {
//        throw std::runtime_error("It does not work with negative angles.");
//    }
//
//    double lambda_deg = Ogisu.get_lambda_rad() * math::constant::R2D();
//    double phi_deg    = Ogisu.get_phi_rad() * math::constant::R2D();
//    double h_m        = Ogisu.get_h_m();
//    double psi_deg    = psi_rad * math::constant::R2D();
//    double theta_deg  = theta_rad * math::constant::R2D();
//
//    osgEarth::Viewpoint Oviewpoint("Everett", lambda_deg, phi_deg, h_m, psi_deg - 90.0, theta_deg - 90.0, 0);
//
//    // set the viewpoint
//    _Pmanip->setViewpoint(Oviewpoint);
}
/* set eagle camera position and orientation. */

void eagle::planet_eagle::obtain_image(cv::Mat& Omat) const {
//    _Pconverter->obtain_image(Omat);
}
/* generates an OpenCV monochrome image based on the previously saved OpenSceneGraph RGB image according to the converter policy, filling
 * up the input object. */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////





