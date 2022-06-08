#include "earth_eagle.h"
#include "earth_eagle_final.h"
#include "earth_eagle_att.h"
#include "earth_eagle_osg.h"
#include "converter.h"
#include "math/logic/share.h"
#include "math/logic/seeder.h"
#include "math/logic/log_writer_mgr.h"
#include "ang/rotate/euler.h"
#include "env/coord.h"
#include "acft/sens/camera.h"
#include <boost/filesystem.hpp>
#include <string>

// CLASS EARTH EAGLE
// =================
// =================

eagle::earth_eagle::earth_eagle(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale,
                                   const std::string& st_images_name, const std::string& st_map_filename, int* Xargc, char **Xargv, bool flag)
: _image_id(image_id), _Pcam(&Ocam), _LOD_scale(LOD_scale), _st_output_nav(math::share::get_location(folder_nav)), _st_output_eagle(math::share::get_location(folder_eagle)) {

    boost::filesystem::path path_eagle("eagle");

    boost::filesystem::path path_maps("maps");
    boost::filesystem::path path_earth_files("earth_files");

    std::string st_case_guid = st_images_name.substr(0,2);
    boost::filesystem::path path_guid(st_case_guid);

    std::string st_zone = st_images_name.substr(3,2);
    boost::filesystem::path path_zone(st_zone);

    // complete path of each image excluding number and extension (folder emptied in function process text file, not here) */
    boost::filesystem::path path_outputs_nav(_st_output_nav);
    boost::filesystem::path path_outputs_eagle(_st_output_eagle);
    boost::filesystem::path path_images_filename(st_images_name);
    boost::filesystem::path path_images_filename_(st_images_name + "_");

    boost::filesystem::create_directory(path_outputs_eagle);
    boost::filesystem::create_directory(path_outputs_eagle / path_eagle);
    boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid);
    boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid / path_zone);
    boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename);
    _st_folder_name     = (path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename).string();
    _st_image_name_file = (path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename / path_images_filename_).string();

    // complete path of earth maps
    boost::filesystem::path path_config(math::share::condor_input);
    boost::filesystem::path path_map_filename(st_map_filename);
    std::string st_earth_filepath = (path_config / path_eagle / path_maps / path_earth_files / path_map_filename).string();

    // file extension
    switch (image_id) {
        case eagle::logic::image_jpg:
            _st_extensiondot = ".jpg";
            _st_extension = "jpg";
            break;
        case eagle::logic::image_png:
            _st_extensiondot = ".png";
            _st_extension = "png";
            break;
        case eagle::logic::image_bmp:
            _st_extensiondot = ".bmp";
            _st_extension = "bmp";
            break;
        default:
            throw std::runtime_error("Invalid image format option.");
    }

    if (flag == true) {
        // arguments to be passed to osgViewer Viewer
        // pointer to integer containing number of arguments passed to osgViewer Viewer, which are the same of the executable plus one added here
        int *argc = Xargc;
        // pointer to (array of C-string pointers), which is the same as weak pointer to (pointer to C-string pointers) containing
        // the arguments passed to osgViewer Viewer, which are the same same of the executable plus one added here.
        char **argv = Xargv;
        argv[1] = new char[st_earth_filepath.size() + 1];
        std::strcpy(argv[1], st_earth_filepath.c_str());
        ++(*argc);
        osg::ArgumentParser Oarguments(argc, argv); // arguments for eagle initialization

        //std::cout << argv[1] << std::endl; /////////////////////

        // eagle initialized with arguments so it can capture path of earth map
        _Peagle = new osgViewer::Viewer(Oarguments);

        // add Earth map to eagle (TAKES A REALLY LONG TIME, ENSURE IT IS ONLY DONE ONCE) !!!!!!!!!
        _Pnode = osgDB::readNodeFiles(Oarguments);

        // the next line was indicated by Jason from Pelican, and achieves lighter images. Without it, images look too dark.
        _Pnode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);

        _Peagle->setSceneData(_Pnode.get());
        _Pellipsoid = osgEarth::MapNode::get(_Pnode)->getTerrain()->getSRS()->getEllipsoid(); // get ellipsoid

        // converter (from osg to opencv)
        switch (converter_id) {
            case eagle::logic::converter_osg_vis_opencv_no:
                _Pconverter = new eagle::converter_osg_vis_opencv_no();
                break;
            case eagle::logic::converter_osg_sav_opencv_no:
                _Pconverter = new eagle::converter_osg_sav_opencv_no();
                break;
            case eagle::logic::converter_osg_sav_opencv_vis:
                _Pconverter = new eagle::converter_osg_sav_opencv_vis();
                break;
            case eagle::logic::converter_osg_sav_opencv_sav:
                _Pconverter = new eagle::converter_osg_sav_opencv_sav();
                break;
            case eagle::logic::converter_osg_sav_opencv_vis_stream:
                _Pconverter = new eagle::converter_osg_sav_opencv_vis_stream();
                break;
            case eagle::logic::converter_osg_sav_opencv_sav_stream:
                _Pconverter = new eagle::converter_osg_sav_opencv_sav_stream();
                break;
            case eagle::logic::converter_osg_str_opencv_no:
                _Pconverter = new eagle::converter_osg_str_opencv_no(_st_extension);
                break;
            case eagle::logic::converter_osg_str_opencv_vis:
                _Pconverter = new eagle::converter_osg_str_opencv_vis(_st_extension);
                break;
            case eagle::logic::converter_osg_str_opencv_sav:
                _Pconverter = new eagle::converter_osg_str_opencv_sav(_st_extension);
                break;
            case eagle::logic::converter_osg_savstr_opencv_vis:
                _Pconverter = new eagle::converter_osg_savstr_opencv_vis(_st_extension);
                break;
            case eagle::logic::converter_osg_savstr_opencv_sav:
                _Pconverter = new eagle::converter_osg_savstr_opencv_sav(_st_extension);
                break;
            default:
                throw std::runtime_error("Invalid converter option.");
        }

        // add screen capture to eagle
        _Pcapture = new osgViewer::ScreenCaptureHandler(_Pconverter);
        _Peagle->addEventHandler(_Pcapture);

        // tell the database pager to not modify the unref settings
        _Peagle->getDatabasePager()->setUnrefImageDataAfterApplyPolicy(true, false);
        // seems to be something about discarding small info in the images, but I do not see any difference
        _Peagle->getCamera()->setSmallFeatureCullingPixelSize(-1.0f);

        // new set projection matrix (????????)
        double aspect_ratio = (double) (_Pcam->get_width_px()) / (double) (_Pcam->get_height_px());
        double near_z = 1;
        double far_z = 10000;
        _Peagle->getCamera()->setProjectionMatrixAsPerspective(_Pcam->get_vfov_deg(), aspect_ratio, near_z, far_z);
        _Peagle->getCamera()->setProjectionResizePolicy(osg::Camera::FIXED);
        _Peagle->getCamera()->setLODScale(_LOD_scale); // default is 1, smaller more detail but slower
        _Peagle->setUpViewInWindow(0, 0, (int)(_Pcam->get_width_px()), (int)(_Pcam->get_height_px())); // deprecated but replacement below does not work
        //_Peagle->apply(new osgViewer::SingleWindow(0, 0, width_px, height_px));

        // thread-safe initialization of the OSG wrapper manager. Calling this here prevents the "unsupported wrapper" messages from OSG
        osgDB::Registry::instance()->getObjectWrapperManager()->findWrapper("osg::Image");

        // realize eagle
        _Peagle->realize();

        // DO NOT DELETE, as it is required for get_terrain_altitude methods, which at this time is not employed
        //auto Pquery = new osgEarth::Annotation::ElevationQuery(osgEarth::MapNode::get(_Pnode)->getMap());
        //Pquery->setMaxTilesToCache(10);
        //Pquery->setFallBackOnNoData(false);
    }
}
/* constructor based on eagle camera, enumerator indicating the image format employed, enumerator indicating file conversion option (from
 * osg to opencv), LOD scale (lower better - aim for 0.5-0.8), name of folder where images are to be saved (if applicable), name of file
 * containing the Earth maps, weak pointer to number of arguments in executable (AGAINST RULE), and weak pointer to array of C-string
 * pointers containing arguments in executable (AGAINST RULE). The flag indicates whether (true) to load all class attributes, or (false) to only load those not
 * related with Open Scene Graph, which applies to the blind derivate class. */

//eagle::earth_eagle::earth_eagle(const earth_eagle& O)
//: _image_id(O._image_id), _Pcam(O._Pcam), _LOD_scale(O._LOD_scale), _st_output_nav(O._st_output_nav), _st_output_eagle(O._st_output_eagle),
//_st_folder_name(O._st_folder_name), _st_image_name_file(O._st_image_name_file), _st_extensiondot(O._st_extensiondot),
//_st_extension(O._st_extension), _st_filename_rgb(O._st_filename_rgb), _st_filename_mono(O._st_filename_mono) {
//
//
//
//
//
//
//        boost::filesystem::path path_eagle("eagle");
//
//        boost::filesystem::path path_maps("maps");
//        boost::filesystem::path path_earth_files("earth_files");
//
//        std::string st_case_guid = st_images_name.substr(0,2);
//        boost::filesystem::path path_guid(st_case_guid);
//
//        std::string st_zone = st_images_name.substr(3,2);
//        boost::filesystem::path path_zone(st_zone);
//
//        // complete path of each image excluding number and extension (folder emptied in function process text file, not here) */
//        boost::filesystem::path path_outputs_nav(_st_output_nav);
//        boost::filesystem::path path_outputs_eagle(_st_output_eagle);
//        boost::filesystem::path path_images_filename(st_images_name);
//        boost::filesystem::path path_images_filename_(st_images_name + "_");
//
//        boost::filesystem::create_directory(path_outputs_eagle);
//        boost::filesystem::create_directory(path_outputs_eagle / path_eagle);
//        boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid);
//        boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid / path_zone);
//        boost::filesystem::create_directory(path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename);
//        _st_folder_name     = (path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename).string();
//        _st_image_name_file = (path_outputs_eagle / path_eagle / path_guid / path_zone / path_images_filename / path_images_filename_).string();
//
//        // complete path of earth maps
//        boost::filesystem::path path_config(math::share::condor_input);
//        boost::filesystem::path path_map_filename(st_map_filename);
//        std::string st_earth_filepath = (path_config / path_eagle / path_maps / path_earth_files / path_map_filename).string();
//
//        // file extension
//        switch (image_id) {
//            case eagle::logic::image_jpg:
//                _st_extensiondot = ".jpg";
//                _st_extension = "jpg";
//                break;
//            case eagle::logic::image_png:
//                _st_extensiondot = ".png";
//                _st_extension = "png";
//                break;
//            case eagle::logic::image_bmp:
//                _st_extensiondot = ".bmp";
//                _st_extension = "bmp";
//                break;
//            default:
//                throw std::runtime_error("Invalid image format option.");
//        }
//
//        if (flag == true) {
//            // arguments to be passed to osgViewer Viewer
//            // pointer to integer containing number of arguments passed to osgViewer Viewer, which are the same of the executable plus one added here
//            int *argc = Xargc;
//            // pointer to (array of C-string pointers), which is the same as weak pointer to (pointer to C-string pointers) containing
//            // the arguments passed to osgViewer Viewer, which are the same same of the executable plus one added here.
//            char **argv = Xargv;
//            argv[1] = new char[st_earth_filepath.size() + 1];
//            std::strcpy(argv[1], st_earth_filepath.c_str());
//            ++(*argc);
//            osg::ArgumentParser Oarguments(argc, argv); // arguments for eagle initialization
//
//            //std::cout << argv[1] << std::endl; /////////////////////
//
//            // eagle initialized with arguments so it can capture path of earth map
//            _Peagle = new osgViewer::Viewer(Oarguments);
//
//            // add Earth map to eagle (TAKES A REALLY LONG TIME, ENSURE IT IS ONLY DONE ONCE) !!!!!!!!!
//            _Pnode = osgDB::readNodeFiles(Oarguments);
//
//            // the next line was indicated by Jason from Pelican, and achieves lighter images. Without it, images look too dark.
//            _Pnode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
//
//            _Peagle->setSceneData(_Pnode.get());
//            _Pellipsoid = osgEarth::MapNode::get(_Pnode)->getTerrain()->getSRS()->getEllipsoid(); // get ellipsoid
//
//            // converter (from osg to opencv)
//            switch (converter_id) {
//                case eagle::logic::converter_osg_vis_opencv_no:
//                    _Pconverter = new eagle::converter_osg_vis_opencv_no();
//                    break;
//                case eagle::logic::converter_osg_sav_opencv_no:
//                    _Pconverter = new eagle::converter_osg_sav_opencv_no();
//                    break;
//                case eagle::logic::converter_osg_sav_opencv_vis:
//                    _Pconverter = new eagle::converter_osg_sav_opencv_vis();
//                    break;
//                case eagle::logic::converter_osg_sav_opencv_sav:
//                    _Pconverter = new eagle::converter_osg_sav_opencv_sav();
//                    break;
//                case eagle::logic::converter_osg_sav_opencv_vis_stream:
//                    _Pconverter = new eagle::converter_osg_sav_opencv_vis_stream();
//                    break;
//                case eagle::logic::converter_osg_sav_opencv_sav_stream:
//                    _Pconverter = new eagle::converter_osg_sav_opencv_sav_stream();
//                    break;
//                case eagle::logic::converter_osg_str_opencv_no:
//                    _Pconverter = new eagle::converter_osg_str_opencv_no(_st_extension);
//                    break;
//                case eagle::logic::converter_osg_str_opencv_vis:
//                    _Pconverter = new eagle::converter_osg_str_opencv_vis(_st_extension);
//                    break;
//                case eagle::logic::converter_osg_str_opencv_sav:
//                    _Pconverter = new eagle::converter_osg_str_opencv_sav(_st_extension);
//                    break;
//                case eagle::logic::converter_osg_savstr_opencv_vis:
//                    _Pconverter = new eagle::converter_osg_savstr_opencv_vis(_st_extension);
//                    break;
//                case eagle::logic::converter_osg_savstr_opencv_sav:
//                    _Pconverter = new eagle::converter_osg_savstr_opencv_sav(_st_extension);
//                    break;
//                default:
//                    throw std::runtime_error("Invalid converter option.");
//            }
//
//            // add screen capture to eagle
//            _Pcapture = new osgViewer::ScreenCaptureHandler(_Pconverter);
//            _Peagle->addEventHandler(_Pcapture);
//
//            // tell the database pager to not modify the unref settings
//            _Peagle->getDatabasePager()->setUnrefImageDataAfterApplyPolicy(true, false);
//            // seems to be something about discarding small info in the images, but I do not see any difference
//            _Peagle->getCamera()->setSmallFeatureCullingPixelSize(-1.0f);
//
//            // new set projection matrix (????????)
//            double aspect_ratio = (double) (_Pcam->get_width_px()) / (double) (_Pcam->get_height_px());
//            double near_z = 1;
//            double far_z = 10000;
//            _Peagle->getCamera()->setProjectionMatrixAsPerspective(_Pcam->get_vfov_deg(), aspect_ratio, near_z, far_z);
//            _Peagle->getCamera()->setProjectionResizePolicy(osg::Camera::FIXED);
//            _Peagle->getCamera()->setLODScale(_LOD_scale); // default is 1, smaller more detail but slower
//            _Peagle->setUpViewInWindow(0, 0, (int)(_Pcam->get_width_px()), (int)(_Pcam->get_height_px())); // deprecated but replacement below does not work
//            //_Peagle->apply(new osgViewer::SingleWindow(0, 0, width_px, height_px));
//
//            // thread-safe initialization of the OSG wrapper manager. Calling this here prevents the "unsupported wrapper" messages from OSG
//            osgDB::Registry::instance()->getObjectWrapperManager()->findWrapper("osg::Image");
//
//            // realize eagle
//            _Peagle->realize();
//
//            // DO NOT DELETE, as it is required for get_terrain_altitude methods, which at this time is not employed
//            //auto Pquery = new osgEarth::Annotation::ElevationQuery(osgEarth::MapNode::get(_Pnode)->getMap());
//            //Pquery->setMaxTilesToCache(10);
//            //Pquery->setFallBackOnNoData(false);
//        }
//
//}
/* copy constructor */

void eagle::earth_eagle::obtain_image(cv::Mat& Omat, const std::string& st_i) const {
    _Pconverter->obtain_image(Omat, st_i);
}
/* generates an OpenCV monochrome image based on the previously saved OpenSceneGraph RGB image according to the converter policy, filling
 * up the input object. */

bool eagle::earth_eagle::process_image(const std::string& st_i) {
    _Peagle->frame(); // call frame once so the viewpoint is rendered at least once

    _st_filename_mono = _st_image_name_file + st_i + _st_extensiondot; // generate complete image name
    _st_filename_rgb  = _st_image_name_file + st_i + "_rgb" +  _st_extensiondot; // generate complete image name

    if (boost::filesystem::exists(_st_filename_rgb) == true) {
        //std::cout << _st_filename_rgb << std::endl;
        return true;
    }

    // Wait for all requests to complete.  Not doing this is why some of your screenshots are lower resolution.
    // This could take a long time depending on the scene, so if you have time constraints you could wait until the
    // requests are no longer in progress or some time threshold has been reached if you'd like.
    while (_Peagle->getDatabasePager()->getRequestsInProgress() == true) {
        //OE_NOTICE << "Waiting on " << _st_filename_rgb << std::endl;
        _Peagle->frame();
    }

    // capture the frame
    _Pcapture->setFramesToCapture(1);
    _Pcapture->captureNextFrame(*_Peagle);
    _Pconverter->set_filename(st_i, _st_filename_rgb, _st_filename_mono, _image_id);
    _Peagle->frame(); // saves the image
    _Pcapture->stopCapture();
    return false;
}
/* internal processing of image previously loaded into viewer. Input is image number. Returns true if image already exists
 * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */

ang::homogeneous eagle::earth_eagle::convert_SE3_osg_2_att(const osg::Matrixd& O) {
    Eigen::Matrix3d RE_inv;
    RE_inv << O(0, 0), O(0, 1), O(0, 2), O(1, 0), O(1, 1), O(1, 2), O(2, 0), O(2, 1), O(2, 2);
    ang::dcm R_inv(RE_inv);

    Eigen::Vector3d x;
    x << O(3, 0), O(3, 1), O(3, 2);

    return {R_inv.inverse(), x};
}
/* convert osg Matrixd to ang::homogeneous representation of SE3 */

osg::Matrixd eagle::earth_eagle::convert_SE3_att_2_osg(const ang::homogeneous& H) {
    ang::dcm R_inv = H.get_inverse_dcm();
    Eigen::Vector3d x = H.get_T();
    return osg::Matrixd(R_inv()(0,0), R_inv()(0,1), R_inv()(0,2), 0.,
                        R_inv()(1,0), R_inv()(1,1), R_inv()(1,2), 0.,
                        R_inv()(2,0), R_inv()(2,1), R_inv()(2,2), 0.,
                        x(0),       x(1),       x(2),       1.);
}
/* convert ang::homogeneous representation of SE3 to osg Matrixd */

//double eagle::earth_eagle::get_terrain_altitude(const osgEarth::Annotation::Terrain& Oterrain, osgEarth::Annotation::ElevationQuery& Oquery, const double &lambda_deg, const double &phi_deg) const {
//    double query_resolution = 0;
//    double out_hamsl = 0.0;
//    double out_resolution = 0.0;
//    osgEarth::Annotation::GeoPoint mapPoint;
//    mapPoint.set(Oterrain.getSRS(), osg::Vec3d(lambda_deg, phi_deg, 0), osgEarth::AltitudeMode::ALTMODE_ABSOLUTE);
//    bool ok = Oquery.getElevation(mapPoint, out_hamsl, query_resolution, &out_resolution);
//    return out_hamsl;
//}
/* return terrain altitude (KEY CAPABILITY ALTHOUGH I DO NOT USE IT) */

void eagle::earth_eagle::fix_problems(const std::string& st_folder) {
    std::string st_case_guid = st_folder.substr(0,2);
    std::string st_zone = st_folder.substr(3,2);
    boost::filesystem::path path_outputs_nav(_st_output_nav);
    boost::filesystem::path path_folder1("nav");
    boost::filesystem::path path_folder2(st_case_guid);
    boost::filesystem::path path_folder3(st_zone);
    boost::filesystem::path path_folder4(st_folder);
    boost::filesystem::path path_file("images_truth.txt");
    std::string st_input_trj_file = (path_outputs_nav / path_folder1 / path_folder2 / path_folder3 / path_folder4 / path_file).string();

    std::ifstream Ostream;
    Ostream.open(st_input_trj_file);
    assert(Ostream);

    int id;
    std::string st_line, st_id, st_i;
    std::getline(Ostream, st_line); // first line contains text

    boost::filesystem::path path_filename_rgb;

    while (std::getline(Ostream, st_line)) {
        st_id = st_line.substr(0,6);
        id    = std::stoi(st_id);
        math::digits::to_string_5digits(id, st_i);
        path_filename_rgb = _st_image_name_file + st_i + "_rgb" + _st_extensiondot; // NOLINT(performance-inefficient-string-concatenation)

        if (boost::filesystem::exists(path_filename_rgb) == false) {
            std::string st_phi_rad    = st_line.substr(18,14);
            std::string st_lambda_rad = st_line.substr(32,14);
            std::string st_h_m        = st_line.substr(46,14);
            std::string st_psi_rad    = st_line.substr(60,14);
            std::string st_theta_rad  = st_line.substr(74,14);
            std::string st_xi_rad     = st_line.substr(88,14);
            std::string st_ratio      = st_line.substr(102,14);
            math::log_writer_mgr::get_log().write(std::string("Eagle:          Recreated image # " + st_i + "."));
            env::geodetic_coord x_gdt_b_rad_m(std::stod(st_lambda_rad), std::stod(st_phi_rad), std::stod(st_h_m));
            ang::euler Oeuler_nb(std::stod(st_psi_rad), std::stod(st_theta_rad), std::stod(st_xi_rad));
            double ratio = std::stod(st_ratio);

            //std::cout << id << "   " << ratio << std::endl; //////////////////////////
            //throw std::runtime_error("ABCDE");
            this->process_single_input(x_gdt_b_rad_m, Oeuler_nb, st_i, ratio);
            //this->color2mono(st_i);
        }
    }
    Ostream.close();
}
/* sometimes for unknown reasons an image is not created when it should, either when generating them one by one from the navigator
 * or when generating them all together from the "images_truth.txt" text file. This function reviews the "images_truth.txt" file ensuring that there
 * is an image for every line in the file, recreating those that for some reason were not created before. Fortunately the number of
 * missing images is extremely low, but there is often at least one in every trajectory. */

void eagle::earth_eagle::color2mono(const std::string& st_folder) {
    math::log_writer_mgr::get_log().write(std::string("Eagle:          Starts converting images to monochrome."));
    boost::filesystem::path path_outputs_folder(_st_output_eagle);
    boost::filesystem::path path_eagle("eagle");
    std::string st_case_guid = st_folder.substr(0,2);
    boost::filesystem::path path_guid(st_case_guid);
    boost::filesystem::path path_folder(path_outputs_folder / path_eagle / path_guid / st_folder);

    std::string st_filename_input;
    //std::string st_filename_output;
    //std::string st_ext;
    //std::string st_rgb = "_rgb";
    //int l;
    cv::Mat Omat;

    int count = 0;
    for(auto & it : boost::filesystem::directory_iterator(path_folder)){
        std::cout << count << std::endl; /////////////////////////////
        st_filename_input = it.path().filename().string();
        //l = st_filename_input.find(st_rgb);
        //st_filename_output = st_filename_input.substr(0, l) + ".jpg";
        Omat = cv::imread((path_folder / st_filename_input).string(), cv::IMREAD_GRAYSCALE);
        cv::imwrite((path_folder / st_filename_input).string(), Omat);
        //boost::filesystem::remove(path_folder / st_filename_input);
        count++;
    }
    math::log_writer_mgr::get_log().write(std::string("Eagle:          Finishes converting images to monochrome."));
}
/**< for some reason I can not make the option converter_osg_vis_opencv_sav work as it should to generate black and white images, and
 * I am stuck with the converter_osg_sav_opencv_no that generates color images, which take up extra space. This function opens all
 * color images in the input folder, saves them in black and white, and deletes the color copy. */

void eagle::earth_eagle::remove_color(math::logic::FOLDER output_folder, const std::string& st_folder) {
    std::string st_case_guid = st_folder.substr(0,2);
    std::string st_case_zone = st_folder.substr(3,2);

    boost::filesystem::path path_folder(math::share::get_location(output_folder));
    boost::filesystem::path path_eagle("eagle");
    boost::filesystem::path path_guid(st_case_guid);
    boost::filesystem::path path_zone_in(st_case_zone);
    boost::filesystem::path path_zone_out(st_case_zone + "BW");
    boost::filesystem::path path_seed(st_folder);

    boost::filesystem::path path_full_in(path_folder / path_eagle / path_guid / path_zone_in / path_seed);
    boost::filesystem::path path_full_out(path_folder / path_eagle / path_guid / path_zone_out / path_seed);

    boost::filesystem::create_directory((path_folder / path_eagle / path_guid / path_zone_out).string());
    boost::filesystem::create_directory((path_folder / path_eagle / path_guid / path_zone_out / path_seed).string());

    //std::cout << path_full_in.string() << std::endl;
    //std::cout << path_full_out.string() << std::endl;

    int count = 0;
    std::string st_filename_input, st_filename_output;
    cv::Mat Omat;

    std::vector<int> Vcompression;
    Vcompression.push_back(CV_IMWRITE_JPEG_QUALITY);
    Vcompression.push_back(95); // 0-100 --> 95 default, 100 maximum quality and biggest size
    //Vcompression.push_back(CV_IMWRITE_PNG_COMPRESSION);
    //Vcompression.push_back(9); // 0-9 --> 3 default, all lossless, 9 smallest size, slowest decoding speed

    for(auto & it : boost::filesystem::directory_iterator(path_full_in)){
        std::cout << count << std::endl; /////////////////////////////
        st_filename_input = it.path().filename().string();

        // comment one or the other for outputs in jpg or png
        st_filename_output = st_filename_input;
        //st_filename_output = st_filename_input.substr(0, st_filename_input.size() - 3) + "png";

        //std::cout << "AAA  " << st_filename_input << std::endl;
        //std::cout << "BBB  " << st_filename_noext << std::endl;

        Omat = cv::imread((path_full_in / st_filename_input).string(), cv::IMREAD_GRAYSCALE);
        cv::imwrite((path_full_out / st_filename_output).string(), Omat, Vcompression);
        count++;
    }
}
/* converts all images of a single folder to black and white */

std::string eagle::earth_eagle::obtain_folder(const int& case_guid, const env::logic::ZONE_ID& zone_id, const int& seed_order, const double& turb_factor,
                          const env::logic::OFFSETS_ID& offsets_id, const env::logic::WIND_ID& wind_id, const double& t_sec_gpsloss) {
    std::string st_case_guid = std::to_string(case_guid);
    switch (st_case_guid.size()) {
        case 1:
            st_case_guid.insert(0, "0");
            break;
        case 2:
            break;
        default:
            throw std::runtime_error("Incorrect case guidance choice.");
            break;
    }

    std::string st_zone = eagle::earth_eagle::zone2string(zone_id);
    std::string st_seeds = math::seeder::seed2string(seed_order);

    std::string st_turb_factor = std::to_string((int)(turb_factor * 10));
    switch (st_turb_factor.size()) {
        case 1:
            st_turb_factor.insert(0, "0");
            break;
        case 2:
            break;
        default:
            throw std::runtime_error("Incorrect factor choice.");
            break;
    }
    std::string st_t_sec_gpsloss = std::to_string((int)(t_sec_gpsloss));
    switch (st_t_sec_gpsloss.size()) {
        case 1:
            st_t_sec_gpsloss.insert(0, "000");
            break;
        case 2:
            st_t_sec_gpsloss.insert(0, "00");
            break;
        case 3:
            st_t_sec_gpsloss.insert(0, "0");
            break;
        case 4:
            break;
        default:
            throw std::runtime_error("Incorrect GPS loss choice.");
            break;
    }
    std::string st_offsets = std::to_string(offsets_id);
    switch (st_offsets.size()) {
        case 1:
            st_offsets.insert(0, "0");
            break;
        case 2:
            break;
        default:
            throw std::runtime_error("Incorrect offsets choice.");
            break;
    }
    std::string st_wind = std::to_string(wind_id);
    switch (st_wind.size()) {
        case 1:
            st_wind.insert(0, "0");
            break;
        case 2:
            break;
        default:
            throw std::runtime_error("Incorrect wind choice.");
            break;
    }

    return  st_case_guid + "_" + st_zone + "_" + st_seeds + "_" + st_turb_factor + "_" + st_offsets + "_" + st_wind + "_" + st_t_sec_gpsloss;
    //return  st_case_guid + "_" + st_seeds + "_" + st_turb_factor + "_" + st_offsets + "_" + st_wind + "_" + st_t_sec_gpsloss;
}
/* returns a string describing the execution (to be used as a folder name) based on the guidance case, the location, the seed order, the turbulence factor,
 * the offsets ID, the wind ID, and the time the GPS signal is lost. This function has nothing to do with the Eagle viewer itself, but is
 * located here for convenience. */

eagle::earth_eagle* eagle::earth_eagle::create(eagle::logic::EAGLE_TYPE eagle_type, const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                                  eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename,
                                  int* argc, char **argv) {
    switch (eagle_type) {
        case eagle::logic::eagle_final:
            return new eagle::earth_eagle_final(Ogeo, Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, argc, argv);
        case eagle::logic::eagle_att:
            return new eagle::earth_eagle_att(Ogeo, Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, argc, argv);
        case eagle::logic::eagle_osg:
            return new eagle::earth_eagle_osg(Ocam, image_id, converter_id, folder_nav, folder_eagle, LOD_scale, st_images_name, st_map_filename, argc, argv);
        case eagle::logic::eagle_blind:
            return new eagle::earth_eagle_blind(Ocam, image_id, folder_nav, folder_eagle, st_images_name, st_map_filename, argc, argv);
        case eagle::logic::eagle_size:
        default:
            throw std::runtime_error("Eagle model not available");
    }
}
/* create eagle pointer based on eagle type and rest of constructor attributes */

void eagle::earth_eagle::empty_images_folder() const {
    boost::filesystem::remove_all(_st_folder_name);
    boost::filesystem::create_directory(_st_folder_name);
    //boost::filesystem::remove_all(_st_image_name_file); // Edu Gallo replaced this line (which did nothing) by the two above on 15-May-2022
}
/* empty folder where images are located */

std::string eagle::earth_eagle::zone2string(env::logic::ZONE_ID zone_id) {
    switch (zone_id) {
        case env::logic::zone_default:
            return "DF";
        case env::logic::zone_desert:
            return "DS";
        case env::logic::zone_urban:
            return "UR";
        case env::logic::zone_everglades:
            return "EV";
        case env::logic::zone_forest:
            return "FR";
        case env::logic::zone_farm:
            return "FM";
        case env::logic::zone_mix:
            return "MX";
        case env::logic::zone_mexico:
            return "ME";
        case env::logic::zone_colombia:
            return "CO";
        case env::logic::zone_argentina:
            return "AR";
        case env::logic::zone_canada:
            return "CA";
        case env::logic::zone_rozas:
            return "RZ";
        case env::logic::zone_wisconsin:
            return "99";
        case env::logic::zone_pntt:
            return "PN";
        case env::logic::zone_river:
            return "RI";
        default:
            throw std::runtime_error("Location not available");
    }
}
/* assign a string to the location */

env::logic::ZONE_ID eagle::earth_eagle::string2zone(const std::string& st_zone) {
    if (st_zone == "DF") {return env::logic::zone_default;}
    else if (st_zone == "DS") {return env::logic::zone_desert;}
    else if (st_zone == "UR") {return env::logic::zone_urban;}
    else if (st_zone == "EV") {return env::logic::zone_everglades;}
    else if (st_zone == "FR") {return env::logic::zone_forest;}
    else if (st_zone == "FM") {return env::logic::zone_farm;}
    else if (st_zone == "MX") {return env::logic::zone_mix;}
    else if (st_zone == "ME") {return env::logic::zone_mexico;}
    else if (st_zone == "CO") {return env::logic::zone_colombia;}
    else if (st_zone == "AR") {return env::logic::zone_argentina;}
    else if (st_zone == "CA") {return env::logic::zone_canada;}
    else if (st_zone == "RZ") {return env::logic::zone_rozas;}
    else if (st_zone == "WS") {return env::logic::zone_wisconsin;}
    else if (st_zone == "PN") {return env::logic::zone_pntt;}
    else if (st_zone == "RI") {return env::logic::zone_river;}
    else {
        throw std::runtime_error("Location not available");
    }
}
/* returns location based on string */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////





