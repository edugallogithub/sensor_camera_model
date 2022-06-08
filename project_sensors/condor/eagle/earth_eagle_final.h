#ifndef EARTH_EAGLE_FINAL
#define EARTH_EAGLE_FINAL

#include "eagle.h"
#include "earth_eagle.h"
#include "ang/transform/homogeneous.h"

namespace acft {
    class iner;
}
namespace eagle {

// CLASS EARTH_EAGLE_FINAL
// =======================
// =======================

/* This is the derivative class to employ. Equal to earth_eagle_att, but with less attributes,
 * hence more direct, but more difficult to check intermediate results.
 *
 * The reference systems considered are the following:
 * - Earth Centered Earth Fixed (ECEF)
 * - North-East-Down (NED)
 * - Body Fixed System (BFS)
 * - Camera Reference System (CRS)
 * - Modified Camera Reference system (MRS) --> to employ OpenSceneGraph format
 * - Rotated Modified Camera Reference System (RRS) --> possible 90 [deg] extra if camera rotated
 */

class EAGLE_API earth_eagle_final : public earth_eagle {
private:
    /**< weak pointer to ellipsoidal Earth object */
    const env::geo* const _Pgeo;

    /**< transformation that converts axes (1,2,3) in (2,1,-3), and viceversa (it is symmetric) */
    ang::homogeneous _H_switch;

    /**< set eagle body frame pose (takes body frame to camera frame transformation from camera) */
    void set_camera(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const double& ratio) override;
public:
    /**< default constructor */
    earth_eagle_final() = delete;
    /**< constructor based on Earth ellipsoidal model, eagle camera, enumerator indicating the image format employed,
     * enumerator indicating file conversion option (from osg to opencv), LOD scale (lower better - aim for 0.5-0.8),
     * name of folder where images are to be saved (if applicable), name of file containing the Earth maps, weak pointer to number of
     * arguments in executable (AGAINST RULE), and weak pointer to array of C-string pointers containing arguments in executable (AGAINST RULE) */
    earth_eagle_final(const env::geo& Ogeo, const sens::camera_base& Ocam, eagle::logic::IMAGE_ID image_id,
                     eagle::logic::CONVERTER_ID converter_id, math::logic::FOLDER folder_nav, math::logic::FOLDER folder_eagle, const float& LOD_scale, const std::string& st_images_name, const std::string& st_map_filename,
                      int* argc, char **argv);
    /**< copy constructor */
    earth_eagle_final(const earth_eagle_final&) = delete;
    /**< move constructor */
    earth_eagle_final(earth_eagle_final&&) = delete;
    /**< destructor */
    ~earth_eagle_final() override = default;
    /**< copy assignment */
    earth_eagle_final& operator=(const earth_eagle_final&) = delete;
    /**< move assignment */
    earth_eagle_final& operator=(earth_eagle_final&&) = delete;

    /**< processes all images contained in the input text file (containing body frame position and attitude for the different images) based on
     * OpenSceneGraph, without the conversion to OpenCV. Should be employed with a file (not stream) converter, so images are visualized in RGB
     * and then saved into files. If other converter is employed, images will be visualized but not saved. Automatically numbers images
     * starting with 0. If flag is false, method does nothing. Note that file contains the body pose, not the camera (body to camera transformation
     * automatically added). */
    void process_text_file(const std::string& st_folder, bool flag) override;
    /**< same as process_text_file but reads different file format generated in real flights */
    void process_text_file_ROZAS(const std::string& st_folder, bool flag);

    /**< generates an OpenSceneGraph RGB image based on the input body frame pose (takes body frame to camera frame transformation from camera),
     * visualizing the image in RGB and then following the converter policy, either saving it to RGB file or stream). Note that inputs position
     * and attitude are those of the body, not the camera (body to camera transformation automatically added). Returns true if image already exists
     * with 100% confidence (it already existed), false if it is creating it new and hence may take a while to be created. */
    bool process_single_input(const env::geodetic_coord& x_gdt_b_rad_m, const ang::euler& euler_nb, const std::string& st_i, const double& ratio) override;


}; // closes class earth_eagle_final

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace eagle

#endif

