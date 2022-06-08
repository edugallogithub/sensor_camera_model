#ifndef EARTH_EAGLE
#define EARTH_EAGLE

#include "eagle.h"
#include "logic.h"
#include "converter.h"
#include "math/logic/logic.h"
#include "env/logic.h"
#include <string>
#include <osgEarthAnnotation/FeatureNode>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <opencv2/highgui/highgui.hpp>

/* This file contains the earth_eagle class, intended to generate Earth images obtained from satellite
 * data as if they were taken from a certain position and orientation.
 *
 * The homogeneous transformations required to prepare the eagle can be performed with the Open
 * Scene Graph methods or with my own attitude methods. They are fully equivalent. For that reason, there
 * exist two derivative classes that are fully equivalent. It seems (I made a test) that the attitude
 * transformations are slightly faster, so those should be used.
 *
 * Although I have tried lots of techniques (check converter class) to accelerate the transition from when
 * the image is generated in OpenSceneGraph to where the image is employed in OpenCV, none of them really
 * work as intended because OpenSceneGraph internally opens a new thread to generate the image and I can not
 * control when it is completed. Because of that, I use the converter_osg_sav_opencv_no, which generates
 * the image, then saves it as a text file, avoid using its obtain_image method, and instead manually load
 * the image from OpenCV after a hard coded delay to give time for the file to be properly filled up.
 * This is a huge waste of resources.
 */

/* TWO ISSUES LEFT TO SORT OUT:
 * 1. Rendition of altitude (works in readymap but not in arcgsonline)
 * 2. Cache options in email - ask Jesus
 *
 * QUESTIONS FOR JASON:
 * 1. Rendition of altitude if I can not sort it out
 * 2. I want to avoid converting it to jpg or png, even if I do not save it. Objective is processing
 *    image in opencv. I would like to pass the pixel values directly (kind of matrix), not converting
 *    it to an image format, which takes time and may involve compression (in case of jpg). Is the native
 *    format of osg equal to that of opencv for images? WHAT DO YOU RECOMMEND?? TIME IS CRITICAL.
 * 3. is there a way to avoid visualizing the images? And the eagle itself?
  */

namespace ang {
    class euler;
    class homogeneous;
}
namespace env {
    class geodetic_coord;
    class geo;
}
namespace acft {
    class iner;
}
namespace sens {
    class camera_base;
}

namespace eagle {

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS EARTH_EAGLE
// =================
// =================

class EAGLE_API earth_eagle {
protected:
    /**< enumerator describing the image format employed */
    eagle::logic::IMAGE_ID _image_id;
    /**< weak pointer to eagle camera */
    const sens::camera_base* const _Pcam;
    /**< LOD scale (lower better - aim for 0.5-0.8) */
    float _LOD_scale;

    /**< folder where output files are located */
    std::string _st_output_nav;
    /**< folder where output images are located */
    std::string _st_output_eagle;
    /**< complete path of folder where images are located */
    std::string _st_folder_name;
    /**< complete path of each image excluding number and extension */
    std::string _st_image_name_file;
    /**< file extension starting with dot */
    std::string _st_extensiondot;
    /**< file extension starting without dot */
    std::string _st_extension;
    /**< full name of latest RGB image from OpenSceneGraph */
    std::string _st_filename_rgb;
    /**< full name of latest monochrome image from OpenCV */
    std::string _st_filename_mono;

    /**< pointer to osg Earth ellipsoid */
    osg::ref_ptr<const osg::EllipsoidModel> _Pellipsoid;
    /**< pointer to image converter */
    osg::ref_ptr<eagle::converter> _Pconverter;
    /**< pointer to screen capture */
    osg::ref_ptr<osgViewer::ScreenCaptureHandler> _Pcapture;
    /**< pointer to node containing Earth map */
    osg::ref_ptr<osg::Node> _Pnode;
    /**< smart pointer to eagle */
    osg::ref_ptr<osgViewer::Viewer> _Peagle;

    /**< set eagle body frame pose (takes body frame to camera frame transformation from camera) */
    virtual void set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) = 0;
    /**< internal processing of image previously loaded into viewer. Input is image number. Returns true if image already exists
     * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */
    virtual bool process_image(const std::string& st_i);

    /**< return terrain altitude (KEY CAPABILITY ALTHOUGH I DO NOT USE IT) */
    //virtual double get_terrain_altitude(const osgEarth::Annotation::Terrain& Oterrain, osgEarth::Annotation::ElevationQuery& Oquery, const double &lambda_deg, const double &phi_deg) const;
public:
    /**< default constructor */
    earth_eagle() = delete;
    /**< constructor based on eagle camera, enumerator indicating the image format employed, enumerator indicating file conversion option (from
     * osg to opencv), LOD scale (lower better - aim for 0.5-0.8), name of folder where images are to be saved (if applicable), name of file containing
     * the Earth maps, weak pointer to number of arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing
     * arguments in executable (AGAINST RULE). The flag indicates whether (true) to load all class attributes, or (false) to only load those not
     * related with Open Scene Graph, which applies to the blind derivate class. */
    earth_eagle(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale,
                 const std::string& st_images_name, const std::string& st_map_filename, int* argc, char **argv, bool flag);
    /**< copy constructor */
    earth_eagle(const earth_eagle&) = delete;
    /**< move constructor */
    earth_eagle(earth_eagle&&) = delete;
    /**< destructor */
    virtual ~earth_eagle() = default;
    /**< copy assignment */
    earth_eagle& operator=(const earth_eagle&) = delete;
    /**< move assignment */
    earth_eagle& operator=(earth_eagle&&) = delete;

    /**< processes all images contained in the input text file (containing body frame position and attitude for the different images) based on
     * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
     * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
     * starting with 0. If flag is false, method does nothing. Note that file contains the body pose, not the camera (body to camera transformation
     * automatically added). */
    virtual void process_text_file(const std::string& st_folder, bool flag) = 0;
    /**< generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
     * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
     * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
     * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */
    virtual bool process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) = 0;
    /**< generates an OpenCV monochrome image based on the previously saved OpenSceneGraph RGB image according to the converter policy, filling
     * up the input object. */
    virtual void obtain_image(cv::Mat& Omat, const std::string& st_i) const;
    /**< sometimes for unknown reasons an image is not created when it should, either when generating them one by one from the navigator
     * or when generating them all together from the "images_truth.txt" text file. This function reviews the "images_truth.txt" file ensuring that there
     * is an image for every line in the file, recreating those that for some reason were not created before. Fortunately the number of
     * missing images is extremely low, but there is often at least one in every trajectory. */
    virtual void fix_problems(const std::string& st_folder);
    /**< for some reason I can not make the option converter_osg_vis_opencv_sav work as it should to generate black and white images, and
     * I am stuck with the converter_osg_sav_opencv_no that generates color images, which take up extra space. This function opens all
     * color images in the input folder, saves them in black and white, and deletes the color copy. */
    virtual void color2mono(const std::string& st_folder);
    /**< converts all images of a single folder to black and white */
    virtual void remove_color(math::logic::FOLDER output_folder, const std::string& st_folder);

    /**< empty folder where images are located */
    virtual void empty_images_folder() const;

    /**< convert osg Matrixd to ang::homogeneous representation of SE3 */
    static ang::homogeneous convert_SE3_osg_2_att(const osg::Matrixd& O);
    /**< convert ang::homogeneous representation of SE3 to osg Matrixd */
    static osg::Matrixd convert_SE3_att_2_osg(const ang::homogeneous& H);
    /**< assign a string to the location */
    static std::string zone2string(env::logic::ZONE_ID);
    /**< returns location based on string */
    static env::logic::ZONE_ID string2zone(const std::string& st_zone);

    /**< returns a string describing the execution (to be used as a folder name) based on the guidance case, the location, the seed order, the turbulence factor,
     * the offsets ID, the wind ID, and the time the GPS signal is lost. This function has nothing to do with the Eagle viewer itself, but is
     * located here for convenience. */
    static std::string obtain_folder(const int& case_guid, const env::logic::ZONE_ID& zone_id, const int& seed_order, const double& turb_factor,
                                     const env::logic::OFFSETS_ID& offsets_id, const env::logic::WIND_ID& wind_id, const double& t_sec_gpsloss);
    /**< create eagle pointer based on eagle type and rest of constructor attributes */
    static eagle::earth_eagle* create(eagle::logic::EAGLE_TYPE eagle_type, const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                                      eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename,
                                      int* argc, char **argv);
    /**< get converter */
    const eagle::converter& get_converter() const {return *_Pconverter;}
    /**< get complete path of each image excluding number and extension */
    const std::string& get_st_image_name_file() const {return _st_image_name_file;}
}; // closes class earth_eagle

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS EARTH_EAGLE_BLIND
// =======================
// =======================

/* This is just the dummy that does not do anything, but is useful to replace any of the other derivate eagle
 * classes when no images are required. */

class EAGLE_API earth_eagle_blind : public earth_eagle {
private:
    /**< set eagle body frame pose (takes body frame to camera frame transformation from camera). In this class it does not do anything. */
    void set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) override {}
public:
    /**< default constructor */
    earth_eagle_blind() = delete;
    /**< constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed, flag indicating
     * camera rotation, enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
     * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
     * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */
    earth_eagle_blind(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,  math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const std::string& st_images_name, const std::string& st_map_filename, int* argc, char** argv)
            : earth_eagle(Ocam, image_id, eagle::logic::converter_size, folder_nav, folder_eagle, 0., st_images_name, st_map_filename, argc, argv, false) {}
    /**< copy constructor */
    earth_eagle_blind(const earth_eagle_blind&) = delete;
    /**< move constructor */
    earth_eagle_blind(earth_eagle_blind&&) = delete;
    /**< destructor */
    ~earth_eagle_blind() override = default;
    /**< copy assignment */
    earth_eagle_blind& operator=(const earth_eagle_blind&) = delete;
    /**< move assignment */
    earth_eagle_blind& operator=(earth_eagle_blind&&) = delete;

    /**< empty folder where images are located. In this case it deletes nothing as it is intended to use already created images. */
    void empty_images_folder() const {}
    /**< processes all images contained in the input text file (containing camera position and attitude for the different images) based on
     * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
     * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
     * starting with 0. In this class it does not do anything, no matter what flag contains. Note that file contains the body pose, not the
     * camera (body to camera transformation automatically added).*/
    void process_text_file(const std::string& st_folder, bool flag) override {}
    /**< generates an OpenSceneGraph image based on the input body position and attitude, visualizing the image in RGB and then following
     * the converter policy, either saving it to RGB file or stream). In this class it does not do anything. Note that inputs position
     * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already
     * exists with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created.*/
    bool process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) override {return false;}
    /**< for some reason I can not make the option converter_osg_vis_opencv_sav work as it should to generate black and white images, and
     * I am stuck with the converter_osg_sav_opencv_no that generates color images, which take up extra space. This function opens all
     * color images in the input folder, saves them in black and white, and deletes the color copy. */
    void color2mono(const std::string& st_folder) override {}
    /**< sometimes for unknown reasons an image is not created when it should, either when generating them one by one from the navigator
     * or when generating them all together from the "images_truth.txt" text file. This function reviews the "images_truth.txt" file ensuring that there
     * is an image for every line in the file, recreating those that for some reason were not created before. Fortunately the number of
     * missing images is extremely low, but there is often at least one in every trajectory. */
    void fix_problems(const std::string& st_folder) override {}
}; // closes class earth_eagle_blind

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace eagle

#endif

