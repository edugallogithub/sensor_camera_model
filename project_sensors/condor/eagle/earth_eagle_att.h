#ifndef EARTH_EAGLE_ATT
#define EARTH_EAGLE_ATT

#include "eagle.h"
#include "earth_eagle.h"
#include "ang/transform/homogeneous.h"

namespace acft {
    class iner;
}
namespace  eagle {

// CLASS EARTH_EAGLE_ATT
// =====================
// =====================

class EAGLE_API earth_eagle_att : public earth_eagle {
private:
    /**< weak pointer to ellipsoidal Earth object */
    const env::geo* const _Pgeo;

    /**< transformation between East-North-Up (ENU) and North-East-Down (NED) */
    ang::homogeneous _H_enuned;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and East-North-Up (ENU) */
    //ang::homogeneous _H_ecefenu; // not necessary but do not delete just in case
    /**< transformation between Earth Centered Earth Fixed (ECEF) and North-East-Down (NED) */
    ang::homogeneous _H_ecefned;
    /**< transformation between North-East-Down (NED) and Body Fixed System (BFS) */
    ang::homogeneous _H_nedbfs;
    /**< transformation between Body Fixed System (BFS) and Camera Reference System (CRS) */
    ang::homogeneous _H_bfscrs;
    /**< transformation between North-East-Down (NED) and Camera Reference System (CRS) */
    ang::homogeneous _H_nedcrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Camera Reference System (CRS) */
    ang::homogeneous _H_ecefcrs;
    /**< transformation between Camera Reference System (CRS) and Modified Camera Reference System (MRS) */
    ang::homogeneous _H_crsmrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Modified Camera Reference System (MRS) */
    ang::homogeneous _H_ecefmrs;
    /**< transformation between Modified Camera Reference system (MRS) and Rotated Modified Camera Reference System (RRS) */
    //ang::homogeneous _H_mrsrrs;
    /**< transformation between Earth Centered Earth Fixed (ECEF) and Rotated Modified Camera Reference System (RRS) */
    ang::homogeneous _H_ecefrrs;
    /**< transformation between Rotated Modified Camera Reference System (RRS) and Earth Centered Earth Fixed (ECEF) */
    ang::homogeneous _H_rrsecef;

    /**< set eagle body frame pose (takes body frame to camera frame transformation from camera) */
    void set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) override;
public:
    /**< default constructor */
    earth_eagle_att() = delete;
    /**< constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed,
     * enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
     * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
     * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */
    earth_eagle_att(const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                     eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename,
                     int* argc, char **argv);
    /**< copy constructor */
    earth_eagle_att(const earth_eagle_att&) = delete;
    /**< move constructor */
    earth_eagle_att(earth_eagle_att&&) = delete;
    /**< destructor */
    ~earth_eagle_att() override = default;
    /**< copy assignment */
    earth_eagle_att& operator=(const earth_eagle_att&) = delete;
    /**< move assignment */
    earth_eagle_att& operator=(earth_eagle_att&&) = delete;

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
    const ang::homogeneous& get_H_enuned() const {return _H_enuned;}
    //const ang::homogeneous& get_H_ecefenu() const {return _H_ecefenu;} // not necessary but do not delete just in case
    const ang::homogeneous& get_H_ecefned() const {return _H_ecefned;}
    const ang::homogeneous& get_H_nedbfs() const {return _H_nedbfs;}
    const ang::homogeneous& get_H_bfscrs() const {return _H_bfscrs;}
    const ang::homogeneous& get_H_nedcrs() const {return _H_nedcrs;}
    const ang::homogeneous& get_H_ecefcrs() const {return _H_ecefcrs;}
    const ang::homogeneous& get_H_crsmrs() const {return _H_crsmrs;}
    const ang::homogeneous& get_H_ecefmrs() const {return _H_ecefmrs;}
    //const ang::homogeneous& get_H_mrsrrs() const {return _H_mrsrrs;}
    const ang::homogeneous& get_H_ecefrrs() const {return _H_ecefrrs;}
    const ang::homogeneous& get_H_rrsecef() const {return _H_rrsecef;}

}; // closes class earth_eagle_att

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace eagle

#endif

