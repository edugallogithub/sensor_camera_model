#ifndef EARTH_EAGLE_OSG
#define EARTH_EAGLE_OSG

#include "eagle.h"
#include "earth_eagle.h"

namespace acft {
    class iner;
}
namespace  eagle {

// CLASS EARTH_EAGLE_OSG
// =====================
// =====================

class EAGLE_API earth_eagle_osg : public earth_eagle {
private:
    /**< transformation between East-North-Up (ENU) and North-East-Down (NED) */
    osg::Matrixd _S_enuned;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and East-North-Up (ENU) */
    osg::Matrixd _S_ecefenu;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and North-East-Down (NED) */
    osg::Matrixd _S_ecefned;
    /**< transformation between North-East-Down (NED) and Body Fixed Syxtem (BFS) */
    osg::Matrixd _S_nedbfs;
    /**< transformation between Body Fixed System (BFS) and Camera Reference System (CRS) */
    osg::Matrixd _S_bfscrs;
    /**< transformation between North-East-Down (NED) and Camera Reference System (CRS) */
    osg::Matrixd _S_nedcrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Camera Reference System (CRS) */
    osg::Matrixd _S_ecefcrs;
    /**< transformation between Camera Reference System (CRS) and Modified Camera Reference System (MRS) */
    osg::Matrixd _S_crsmrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Modified Camera Reference System (MRS) */
    osg::Matrixd _S_ecefmrs;
    /**< transformation between Modified Camera Reference system (MRS) and Rotated Modified Camera Reference System (RRS) */
    //osg::Matrixd _S_mrsrrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Rotated Modified Camera Reference System (RRS) */
    osg::Matrixd _S_ecefrrs;
    /**< transformation between Rotated Modified Camera Reference System (RRS) and Earth Centered Earth Fixed (ECEF) */
    osg::Matrixd _S_rrsecef;

    /**< set eagle body frame pose (takes body frame to camera frame transformation from camera) */
    void set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) override;
public:
    /**< default constructor */
    earth_eagle_osg() = delete;
    /**< constructor based on eagle camera, enumerator indicating the image format employed, enumerator indicating
     * file conversion option (from osg to opencv), name of folder where images are to be saved (if applicable), name of file containing
     * the Earth maps, weak pointer to number of arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing
     * arguments in executable (AGAINST RULE) */
    earth_eagle_osg(const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id, eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale,
                     const std::string& st_images_name, const std::string& st_map_filename, int* argc, char **argv);
    /**< copy constructor */
    earth_eagle_osg(const earth_eagle_osg&) = delete;
    /**< move constructor */
    earth_eagle_osg(earth_eagle_osg&&) = delete;
    /**< destructor */
    ~earth_eagle_osg() override = default;
    /**< copy assignment */
    earth_eagle_osg& operator=(const earth_eagle_osg&) = delete;
    /**< move assignment */
    earth_eagle_osg& operator=(earth_eagle_osg&&) = delete;

    /**< processes all images contained in the input text file (containing body frame position and attitude for the different images) based on
     * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
     * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
     * starting with 0. If flag is false, method does nothing. Note that file contains the body pose, not the camera (body to camera transformation
     * automatically added). */
    void process_text_file(const std::string& st_folder, bool flag) override;
    /**< generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
     * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
     * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
     * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */
    bool process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) override;

    /**< Getters intended for testing */
    const osg::Matrixd& get_S_enuned() const {return _S_enuned;}
    const osg::Matrixd& get_S_ecefenu() const {return _S_ecefenu;}
    const osg::Matrixd& get_S_ecefned() const {return _S_ecefned;}
    const osg::Matrixd& get_S_nedbfs() const {return _S_nedbfs;}
    const osg::Matrixd& get_S_bfscrs() const {return _S_bfscrs;}
    const osg::Matrixd& get_S_nedcrs() const {return _S_nedcrs;}
    const osg::Matrixd& get_S_ecefcrs() const {return _S_ecefcrs;}
    const osg::Matrixd& get_S_crsmrs() const {return _S_crsmrs;}
    const osg::Matrixd& get_S_ecefmrs() const {return _S_ecefmrs;}
    //const osg::Matrixd& get_S_mrsrrs() const {return _S_mrsrrs;}
    const osg::Matrixd& get_S_ecefrrs() const {return _S_ecefrrs;}
    const osg::Matrixd& get_S_rrsecef() const {return _S_rrsecef;}

}; // closes class earth_eagle_osg

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace eagle

#endif

