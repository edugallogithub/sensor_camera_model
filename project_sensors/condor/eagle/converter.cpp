#include "converter.h"
#include <stdio.h>
#include <osgDB/WriteFile>

using namespace std;

// CLASS CONVERTER
// ===============
// ===============

eagle::converter::converter()
: _st_i("99999"), _filename_rgb("screenshot_rgb.jpg"), _filename_mono("screenshot.jpg"), _image_id(eagle::logic::image_size) {
}
/* constructor */

void eagle::converter::set_filename(const std::string& st_i, const std::string& filename_rgb, const std::string& filename_mono, eagle::logic::IMAGE_ID image_id) {
    _st_i          = st_i;
    _filename_rgb  = filename_rgb;
    _filename_mono = filename_mono;
    _image_id      = image_id;
    switch (_image_id) {
        case eagle::logic::image_jpg:
            _Vcompress = std::vector<int>(2);
            _Vcompress[0] = cv::IMWRITE_JPEG_QUALITY;
            _Vcompress[1] = 100;
            break;
        case eagle::logic::image_png:
            _Vcompress = std::vector<int>(2);
            _Vcompress[0] = cv::IMWRITE_PNG_COMPRESSION;
            _Vcompress[1] = 0;
            break;
        case eagle::logic::image_bmp:
            _Vcompress = std::vector<int>(0);
            break;
        case eagle::logic::image_size:
        default:
            throw std::runtime_error("Invalid image format.");
    }
}
/* set the name of the images to be saved */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 01 CONVERTER_OSG_VIS_OPENCV_NO
// ====================================
// ====================================

void eagle::converter_osg_vis_opencv_no::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_vis_opencv_no::obtain_image(cv::Mat& Omat, const std::string& st_i) {
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 02 CONVERTER_OSG_SAV_OPENCV_NO
// ====================================
// ====================================

void eagle::converter_osg_sav_opencv_no::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_sav_opencv_no::obtain_image(cv::Mat& Omat, const std::string& st_i) {
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 03 CONVERTER_OSG_SAV_OPENCV_VIS
// =====================================
// =====================================

void eagle::converter_osg_sav_opencv_vis::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    //osg::setNotifyLevel(osg::ALWAYS); //////////////////////////////////////

    // saves RGB image into _filename_rgb file location
    bool flag = osgDB::writeImageFile(Oimage, _filename_rgb);
    std::cout << _st_i << "  " << flag << std::endl; /////////////////////////
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_sav_opencv_vis::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }
    Omat = cv::imread(_filename_rgb, cv::IMREAD_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 04 CONVERTER_OSG_SAV_OPENCV_SAV
// =====================================
// =====================================

void eagle::converter_osg_sav_opencv_sav::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_sav_opencv_sav::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }
    Omat = cv::imread(_filename_rgb, cv::IMREAD_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    // saves monochrome image into _filename_mono location
    cv::imwrite(_filename_mono, Omat, _Vcompress);
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 05 CONVERTER_OSG_SAV_OPENCV_VIS_STREAM
// ============================================
// ============================================

void eagle::converter_osg_sav_opencv_vis_stream::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_sav_opencv_vis_stream::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    // file stream based on _filename_rgb file
    std::ifstream Ostream_img(_filename_rgb, std::ios::binary);

    // obtain stream length
    Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    int buffer_length = Ostream_img.tellg(); // length of stream equals current position
    Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream

    char* Pbuffer = new char[buffer_length]; // char array of proper size
    Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 06 CONVERTER_OSG_SAV_OPENCV_SAV_STREAM
// ============================================
// ============================================

void eagle::converter_osg_sav_opencv_sav_stream::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_sav_opencv_sav_stream::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    // file stram based on _filename_rgb file
    std::ifstream Ostream_img(_filename_rgb, std::ios::binary);

    // obtain stream length
    Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    int buffer_length = Ostream_img.tellg(); // length of stream equals current position
    Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream

    char* Pbuffer = new char[buffer_length]; // char array of proper size
    Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;

    // saves monochrome image into _filename_mono location
    cv::imwrite(_filename_mono, Omat, _Vcompress);
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 07 CONVERTER_OSG_STR_OPENCV_NO
// ====================================
// ====================================

eagle::converter_osg_str_opencv_no::converter_osg_str_opencv_no(const std::string& st_extension)
        : _Preaderwriter(osgDB::Registry::instance()->getReaderWriterForExtension(st_extension)) {
}
/* constructor */

void eagle::converter_osg_str_opencv_no::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // resets stream contents as otherwise it just adds to previous image
    _Ostream_img.str("");

    // saves RGB image into stream (not file)
    _Preaderwriter->writeImage(Oimage, _Ostream_img);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_str_opencv_no::obtain_image(cv::Mat& Omat, const std::string& st_i) {
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 08 CONVERTER_OSG_STR_OPENCV_VIS
// =====================================
// =====================================

eagle::converter_osg_str_opencv_vis::converter_osg_str_opencv_vis(const std::string& st_extension)
: _Preaderwriter(osgDB::Registry::instance()->getReaderWriterForExtension(st_extension)) {
}
/* constructor */

void eagle::converter_osg_str_opencv_vis::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // resets stream contents as otherwise it just adds to previous image
    _Ostream_img.str("");

    // saves RGB image into stream (not file)
    _Preaderwriter->writeImage(Oimage, _Ostream_img);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_str_opencv_vis::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    std::cout << std::endl << " 2nd - obtain vis_vis action begin" << std::endl << std::endl;
    int buffer_length;

    // obtain stream length
    _Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    buffer_length = _Ostream_img.tellg(); // length of stream equals current position
    if (_Ostream_img.fail() == true) { throw std::runtime_error("Failed to read image."); }

    // show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    _Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream
    char *Pbuffer = new char[buffer_length]; // char array of proper size
    _Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 09 CONVERTER_OSG_STR_OPENCV_SAV
// =====================================
// =====================================

eagle::converter_osg_str_opencv_sav::converter_osg_str_opencv_sav(const std::string& st_extension)
        : _Preaderwriter(osgDB::Registry::instance()->getReaderWriterForExtension(st_extension)) {
}
/* constructor */

void eagle::converter_osg_str_opencv_sav::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // resets stream contents as otherwise it just adds to previous image
    _Ostream_img.str("");

    // saves RGB image into stream (not file)
    _Preaderwriter->writeImage(Oimage, _Ostream_img);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_str_opencv_sav::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    int buffer_length;
    // obtain stream length
    _Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    buffer_length = _Ostream_img.tellg(); // length of stream equals current position
    if (_Ostream_img.fail() == true) { throw std::runtime_error("Failed to read image."); }

    // show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    _Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream
    char *Pbuffer = new char[buffer_length]; // char array of proper size
    _Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;

    // saves monochrome image into _filename_mono location
    cv::imwrite(_filename_mono, Omat, _Vcompress);
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 10 CONVERTER_OSG_SAVSTR_OPENCV_VIS
// ========================================
// ========================================

eagle::converter_osg_savstr_opencv_vis::converter_osg_savstr_opencv_vis(const std::string& st_extension)
        : _Preaderwriter(osgDB::Registry::instance()->getReaderWriterForExtension(st_extension)) {
}
/* constructor */

void eagle::converter_osg_savstr_opencv_vis::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);

    // resets stream contents as otherwise it just adds to previous image
    _Ostream_img.str("");

    // saves RGB image into stream (not file)
    _Preaderwriter->writeImage(Oimage, _Ostream_img);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_savstr_opencv_vis::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    int buffer_length;
    // obtain stream length
    _Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    buffer_length = _Ostream_img.tellg(); // length of stream equals current position
    if (_Ostream_img.fail() == true) { throw std::runtime_error("Failed to read image."); }

    //  show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    _Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream
    char *Pbuffer = new char[buffer_length]; // char array of proper size
    _Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS 11 CONVERTER_OSG_SAVSTR_OPENCV_SAV
// ========================================
// ========================================

eagle::converter_osg_savstr_opencv_sav::converter_osg_savstr_opencv_sav(const std::string& st_extension)
        : _Preaderwriter(osgDB::Registry::instance()->getReaderWriterForExtension(st_extension)) {
}
/* constructor */

void eagle::converter_osg_savstr_opencv_sav::operator()(const osg::Image& Oimage, const unsigned int context_id_not_used) {
    // saves RGB image into _filename_rgb file location
    osgDB::writeImageFile(Oimage, _filename_rgb);

    // resets stream contents as otherwise it just adds to previous image
    _Ostream_img.str("");

    // saves RGB image into stream (not file)
    _Preaderwriter->writeImage(Oimage, _Ostream_img);
}
/* process image in OpenSceneGraph */

void eagle::converter_osg_savstr_opencv_sav::obtain_image(cv::Mat& Omat, const std::string& st_i) {
    if (_st_i.compare(st_i) != 0) {
        throw std::runtime_error("Image number does not coincide.");
    }

    int buffer_length;
    // obtain stream length
    _Ostream_img.seekg(0, std::ios::end); // goes to end of stream
    buffer_length = _Ostream_img.tellg(); // length of stream equals current position
    if (_Ostream_img.fail() == true) { throw std::runtime_error("Failed to read image."); }

     // show specific values on console
    //std::cout << static_cast<ushort>(Pbuffer[0]) << std::endl;
    //std::cout << static_cast<ushort>(Pbuffer[buffer_length-1]) << std::endl;

    _Ostream_img.seekg(0, std::ios::beg); // goes to begin of stream
    char *Pbuffer = new char[buffer_length]; // char array of proper size
    _Ostream_img.read(Pbuffer, buffer_length); // fill up buffer from stream

    // loads monochrome image from buffer
    Omat = cv::imdecode(cv::Mat(1, buffer_length, CV_8UC1, Pbuffer), CV_LOAD_IMAGE_GRAYSCALE);
    if (Omat.empty() == true) {
        throw std::runtime_error("Loaded image is empty.");
    }
    delete[] Pbuffer;

    // saves monochrome image into _filename_mono location
    cv::imwrite(_filename_mono, Omat, _Vcompress);
}
/* passes image to OpenCV */

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////








