#ifndef EAGLE_CONVERTER
#define EAGLE_CONVERTER

#include "eagle.h"
#include "logic.h"

#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <opencv2/highgui/highgui.hpp>

/* This file contains different ways of taking an RGB image from OpenSceneGraph and obtaining
 * a monochrome OpenCV Mat object. The process involves two steps (in this order):
 * 1. operator() works with the OpenSceneGraph RGB image (always visualizing it, and then
 *    saving it in file, saving it in local stream, or not saving it).
 * 2. obtain_image() may load the image into monochrome OpenCV (from file or stream), and may
 *    then save it into file.
 *
 * The different derivate classes execute this process in different ways:
 *  1. converter_osg_vis_opencv_no         --> OpenSceneGraph (color): visualizes image
 *                                             OpenCV (b & w):         no action
 *  2. converter_osg_sav_opencv_no         --> OpenSceneGraph (color): visualizes image and saves it in text file
 *                                             OpenCV (b & w):         no action
 *  3. converter_osg_sav_opencv_vis        --> OpenSceneGraph (color): visualizes image and saves it in text file
 *                                             OpenCV (b & w):         reads text file to visualize image
 *  4. converter_osg_sav_opencv_sav        --> OpenSceneGraph (color): visualizes image and saves it in text file
 *                                             OpenCV (b & w):         reads text file to visualize and save image in another text file
 *  5. converter_osg_sav_opencv_vis_stream --> OpenSceneGraph (color): visualizes image and saves it in text file
 *                                             OpenCV (b & w):         reads text file through stream to visualize image
 *  6. converter_osg_sav_opencv_sav_stream --> OpenSceneGraph (color): visualizes image and saves it in text file
 *                                             OpenCV (b & w):         reads text file through stream to visualize and save image in another text file
 *  7. converter_osg_str_opencv_no         --> OpenSceneGraph (color): visualizes image and internally stores it in stream
 *                                             OpenCV (b & w):         no action
 *  8. converter_osg_str_opencv_vis        --> OpenSceneGraph (color): visualizes image and internally stores it in stream
 *                                             OpenCV (b & w):         loads stream to visualize image
 *  9. converter_osg_str_opencv_sav        --> OpenSceneGraph (color): visualizes image and internally stores it in stream
 *                                             OpenCV (b & w):         loads stream to visualize and save image in text file
 * 10. converter_osg_savstr_opencv_vis     --> OpenSceneGraph (color): visualizes image, saves it in text file, and internally stores it in stream
 *                                             OpenCV (b & w):         loads stream to visualize image
 * 11. converter_osg_savstr_opencv_sav     --> OpenSceneGraph (color): visualizes image, saves it in text file, and internally stores it in stream
 *                                             OpenCV (b & w):         loads stream to visualize and save image in text file
 */

namespace eagle {

// CLASS CONVERTER
// ===============
// ===============

class EAGLE_API converter: public osgViewer::ScreenCaptureHandler::CaptureOperation {
protected:
    /**< five digit code with image number */
    std::string _st_i;
    /**< full path of RGB image to be saved to file */
    std::string _filename_rgb;
    /**< full path of monochrome image to be saved to file */
    std::string _filename_mono;
    /**< enumerator specifying format of image to be saved to file */
    eagle::logic::IMAGE_ID _image_id;
    /**< vector of integers containing opencv image saving compression choices */
    std::vector<int> _Vcompress;
public:
    /**< constructor */
    converter();
    /**< process image in OpenSceneGraph */
    virtual void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) = 0;
    /**< passes image to OpenCV */
    virtual void obtain_image(cv::Mat& Omat, const std::string& st_i) = 0;
    /**< set the name of the images to be saved */
    virtual void set_filename(const std::string& st_i, const std::string& filename_rgb, const std::string& filename_mono, eagle::logic::IMAGE_ID image_id);
}; // closes class converter

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 01 CONVERTER_OSG_VIS_OPENCV_NO
// ====================================
// ====================================

/* Visualizes images in OpenSceneGraph, no action in OpenCV */
class EAGLE_API converter_osg_vis_opencv_no : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_vis_opencv_no

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 02 CONVERTER_OSG_SAV_OPENCV_NO
// ====================================
// ====================================

/* Visualizes and saves images in OpenSceneGraph, no action in OpenCV */
class EAGLE_API converter_osg_sav_opencv_no : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_sav_opencv_no

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 03 CONVERTER_OSG_SAV_OPENCV_VIS
// =====================================
// =====================================

/* Visualizes and saves images in OpenSceneGraph, visualizes images in OpenCV (reading previous file) */
class EAGLE_API converter_osg_sav_opencv_vis : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_sav_opencv_vis

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 04 CONVERTER_OSG_SAV_OPENCV_SAV
// =====================================
// =====================================

/* Visualizes and saves images in OpenSceneGraph, visualizes and saves images in OpenCV (reading previous file) */
class EAGLE_API converter_osg_sav_opencv_sav : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_sav_opencv_sav

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 05 CONVERTER_OSG_SAV_OPENCV_VIS_STREAM
// ============================================
// ============================================

/**< Visualizes and saves images in OpenSceneGraph, visualizes images in OpenCV (reading stream of previous file) */
class EAGLE_API converter_osg_sav_opencv_vis_stream : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_sav_opencv_vis_stream

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 06 CONVERTER_OSG_SAV_OPENCV_SAV_STREAM
// ============================================
// ============================================

/**< Visualizes and saves images in OpenSceneGraph, visualizes and saves images in OpenCV (reading stream of previous file) */
class EAGLE_API converter_osg_sav_opencv_sav_stream : public converter {
public:
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_sav_opencv_sav_stream

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 07 CONVERTER_OSG_STR_OPENCV_NO
// ====================================
// ====================================

/**< Visualizes images in OpenSceneGraph, saves it to stream, no action in OpenCV */
class EAGLE_API converter_osg_str_opencv_no : public converter {
private:
    /**< manages reader and writing from stream */
    osg::ref_ptr<osgDB::ReaderWriter> _Preaderwriter;
    /**< buffer where image is temporarily located */
    std::stringstream _Ostream_img;
public:
    /**< default constructor */
    converter_osg_str_opencv_no() = delete;
    /**< constructor */
    explicit converter_osg_str_opencv_no(const std::string& st_extension);
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_str_opencv_no

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 08 CONVERTER_OSG_STR_OPENCV_VIS
// =====================================
// =====================================

/**< Visualizes images in OpenSceneGraph, saves it to stream, and visualizes it in OpenCV */
class EAGLE_API converter_osg_str_opencv_vis : public converter {
private:
    /**< manages reader and writing from stream */
    osg::ref_ptr<osgDB::ReaderWriter> _Preaderwriter;
    /**< buffer where image is temporarily located */
    std::stringstream _Ostream_img;
public:
    /**< default constructor */
    converter_osg_str_opencv_vis() = delete;
    /**< constructor */
    explicit converter_osg_str_opencv_vis(const std::string& st_extension);
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_str_opencv_vis

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 09 CONVERTER_OSG_STR_OPENCV_SAV
// =====================================
// =====================================

/**< Visualizes images in OpenSceneGraph, saves it to stream, and visualizes and saves it in OpenCV */
class EAGLE_API converter_osg_str_opencv_sav : public converter {
private:
    /**< manages reader and writing from stream */
    osg::ref_ptr<osgDB::ReaderWriter> _Preaderwriter;
    /**< buffer where image is temporarily located */
    std::stringstream _Ostream_img;
public:
    /**< default constructor */
    converter_osg_str_opencv_sav() = delete;
    /**< constructor */
    explicit converter_osg_str_opencv_sav(const std::string& st_extension);
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_str_opencv_sav

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 10 CONVERTER_OSG_SAVSTR_OPENCV_VIS
// ========================================
// ========================================

/**< Visualizes and saves images in OpenSceneGraph, saves it to stream, and visualizes it in OpenCV */
class EAGLE_API converter_osg_savstr_opencv_vis : public converter {
private:
    /**< manages reader and writing from stream */
    osg::ref_ptr<osgDB::ReaderWriter> _Preaderwriter;
    /**< buffer where image is temporarily located */
    std::stringstream _Ostream_img;
public:
    /**< default constructor */
    converter_osg_savstr_opencv_vis() = delete;
    /**< constructor */
    explicit converter_osg_savstr_opencv_vis(const std::string& st_extension);
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_savosg_str_opencv_vis

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 11 CONVERTER_OSG_SAVSTR_OPENCV_SAV
// ========================================
// ========================================

/**< Visualizes and saves images in OpenSceneGraph, saves it to stream, and visualizes and saves it in OpenCV */
class EAGLE_API converter_osg_savstr_opencv_sav : public converter {
private:
    /**< manages reader and writing from stream */
    osg::ref_ptr<osgDB::ReaderWriter> _Preaderwriter;
    /**< buffer where image is temporarily located */
    std::stringstream _Ostream_img;
public:
    /**< default constructor */
    converter_osg_savstr_opencv_sav() = delete;
    /**< constructor */
    explicit converter_osg_savstr_opencv_sav(const std::string& st_extension);
    /**< process image in OpenSceneGraph */
    void operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) override;
    /**< passes image to OpenCV */
    void obtain_image(cv::Mat& Omat, const std::string& st_i) override;
}; // closes class converter_osg_savstr_opencv_sav

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace eagle

#endif






