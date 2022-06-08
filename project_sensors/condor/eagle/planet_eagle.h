#ifndef PLANET_EAGLE
#define PLANET_EAGLE

#include "eagle.h"
#include "earth_eagle.h"
#include "converter.h"
#include <string>
#include <osgEarthUtil/EarthManipulator>

/* The difference between planet_eagle and earth_eagle is that planet_eagle employs the attitude
 * representation of osgEarth, not osg itself, which is simpler but lacks bank angles.
 * THIS FUNCTION DOES NOT ACCEPT ROLL ANGLES AND PITCH ONES NEED TO BE POSITIVE.
 * THIS FUNCTION DOES NOT ACCEPT A TRUE FLAG_ROTATE IN THE CONSTRUCTOR
 */

namespace  eagle {

// CLASS PLANET_EAGLE
// ==================
// ==================

class EAGLE_API planet_eagle {
protected:
    /**< enumerator describing the image format employed */
    eagle::logic::IMAGE_ID _image_id;
    /**< weak pointer to eagle camera */
    const sens::camera_base* _Pcam;
    /**< LOD scale (lower better - aim for 0.5-0.8) */
    float _LOD_scale;

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

    /**< pointer to image converter */
    osg::ref_ptr<eagle::converter> _Pconverter;
    /**< pointer to screen capture */
    osg::ref_ptr<osgViewer::ScreenCaptureHandler> _Pcapture;
    /**< pointer to node containing Earth map */
    osg::ref_ptr<osg::Node> _Pnode;
    /**< smart pointer to eagle */
    osg::ref_ptr<osgViewer::Viewer> _Peagle;
    /**< smart pointer to Earth manipulator */
    osg::ref_ptr<osgEarth::Util::EarthManipulator> _Pmanip;

    /**< set eagle camera position and orientation. */
    virtual void set_camera(const env::geodetic_coord& Ogisu, const double& psi_rad, const double& theta_rad);
    /**< internal processing of image previously loaded into viewer. Inputs are image number plus add-on string. */
    virtual void process_image(const std::string& st_i);
public:
    /**< default constructor */
    planet_eagle() = delete;
    /**< constructor based on eagle camera, enumerator indicating the image format employed, enumerator indicating file conversion option (from
     * osg to opencv), LOD scale (lower better - aim for 0.5-0.8), name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number
     * of arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */
    planet_eagle(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, bool flag_rotate,
                       eagle::logic::CONVERTER_ID converter_id, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename,
                       int* argc, char **argv);

    /**< copy constructor */
    planet_eagle(const planet_eagle&) = delete;
    /**< move constructor */
    planet_eagle(planet_eagle&&) = delete;
    /**< destructor */
    virtual ~planet_eagle() = default;
    /**< copy assignment */
    planet_eagle& operator=(const planet_eagle&) = delete;
    /**< move assignment */
    planet_eagle& operator=(planet_eagle&&) = delete;

    /**< processes all images contained in the input text file (containing camera position and attitude for the different images) based on
     * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
     * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images. */
    virtual void process_text_file(const std::string& st_folder);
    /**< generates an OpenSceneGraph RGB image based on the input camera position and attitude, visualizing the image in RGB and then following
     * the converter policy, either saving it to RGB file or stream). */
    virtual void process_single_input(const env::geodetic_coord& Ogisu, const double& psi_rad, const double& theta_rad, const std::string& st_i);
    /**< generates an OpenCV monochrome image based on the previously saved OpenSceneGraph RGB image according to the converter policy, filling
     * up the input object. */
    virtual void obtain_image(cv::Mat& Omat) const;

}; // closes class planet_eagle

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace viewer

#endif

