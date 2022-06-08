#ifndef ACFT_CAMERA
#define ACFT_CAMERA

#include "../acft.h"
#include "ang/rotate/euler.h"
#include "ang/transform/speu_rodrigues.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

/*
 * f = 16 [mm] -->  fov = 57.0924 [deg]
 * f = 19 [mm] -->  fov = 49.2255 [deg]
 * f = 24 [mm] -->  fov = 39.8680 [deg]
 * f = 30 [mm] -->  fov = 32.3584 [deg]
 *
 * Tested sizes:
 *  768 x 1024 [px]
 * 1056 x 1408 [px]
 *  576 x  768 [px]
 *
 */


namespace math {
    class seeder;
}
namespace ang {
    class speu_rodrigues;
}
namespace acft {
    class iner;
}

namespace sens {

// CLASS CAMERA_BASE
// =================
// =================

class ACFT_API camera_base {
protected:
    /**< number of horizontal pixels (number of columns) */
    unsigned int _width_px;
    /**< number of vertical pixels (number of rows) */
    unsigned int _height_px;
    /**< field of view */
    double _vfov_deg;
    /**< focal length [pixels] */
    double _f_px;
    /**< horizontal principal point location (from left side) */
    double _cx_px;
    /**< vertical principal point location (from top) */
    double _cy_px;
public:
    /**< default constructor */
    camera_base() = delete;
    /**< constructor based on images width and height, together with vertical field of view. */
    camera_base(const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg);
    /**< copy constructor */
    camera_base(const camera_base&) = delete;
    /**< move constructor */
    camera_base(camera_base&&) = delete;
    /**< destructor */
    virtual ~camera_base() = default;
    /**< copy assignment */
    camera_base& operator=(const camera_base&) = delete;
    /**< move assignment */
    camera_base& operator=(camera_base&&) = delete;

    /**< projection from image (pixels) to unit sphere camera coordinates */
    Eigen::Vector3d convert_img2uscrs(const Eigen::Vector2d& x_img_px) const;
    /**< projection from image (pixels) to unit sphere camera coordinates */
    Eigen::Vector3d convert_img2uscrs(const double& x_img_pxi, const double& x_img_pxii) const;
    /**< projection from image (pixels) to unit plane camera coordinates */
    Eigen::Vector3d convert_img2upcrs(const Eigen::Vector2d& x_img_px) const;
    /**< projection from image (pixels) to unit plane camera coordinates */
    Eigen::Vector3d convert_img2upcrs(const double& x_img_pxi, const double& x_img_pxii) const;

    /**< transforms cartesian coordinates from camera reference system (CRS) to image reference system (IMG).
    Unit sphere or unit plane camera reference systems (USCRS or UPCRS) are also acceptable inputs. */
    Eigen::Vector2d convert_crs2img(const Eigen::Vector3d& x_crs_m) const;
    /**< transforms cartesian coordinates from unit plane (3rd coord == 1) camera reference system (UPCRS) to image reference system (IMG) */
    Eigen::Vector2d convert_upcrs2img(const Eigen::Vector3d& x_crs_up) const;

    /**< projection from unit sphere camera reference system (USCRS) to unit plane camera reference system (UPCRS) */
    static Eigen::Vector3d convert_uscrs2upcrs(const Eigen::Vector3d& x_crs_us);
    /**< projection from unit plane camera reference system (UPCRS) to unit sphere camera reference system (USCRS) */
    static Eigen::Vector3d convert_upcrs2uscrs(const Eigen::Vector3d& x_crs_up);

    /**< returns true if input pixel coordinates are inside frame and not closer than input boundary to frame border */
    bool is_in_frame(const Eigen::Vector2i& x_img_px, int boundary_px) const;
    /**< returns true if input pixel coordinates are inside frame and not closer than input boundary to frame border */
    bool is_in_frame(const double&x_img_pxi, const double& x_img_pxii, int boundary_px) const;
    /**< returns true if input pixel coordinates (at input pyramid level) are inside frame and not closer than input boundary to frame border */
    bool is_in_frame(const Eigen::Vector2i& x_img_px, int boundary_px, int level) const;
    /**< returns true if input pixel coordinates (at input pyramid level) are inside frame and not closer than input boundary to frame border */
    bool is_in_frame(const double& x_img_pxi, const double& x_img_pxii, int boundary_px, int level) const;

    /**< get effective focal length based on pyramid level */
    double get_f_px_pyramid(const int& pyramid_level) const {return _f_px / sens::camera_base::pyramid_scale(pyramid_level);}
    /**< return the scale corresponding to a given image pyramid level --> 1 - 2 - 4 - 8 - 16 */
    static int pyramid_scale(const int& pyramid_level) {return (1 << pyramid_level);}

    /**< transform input image (pixels) from level 0 to input level */
    static Eigen::Vector2d pyramid_scale_from_level0(const Eigen::Vector2d& x_img_px, const int& pyramid_level);
    /**< transform input pixel distance from level 0 to input level */
    static double pyramid_scale_from_level0(const double& px, const int& pyramid_level);
    /**< transform input image (pixels) to level 0 from input level */
    static Eigen::Vector2d pyramid_scale_to_level0(const Eigen::Vector2d& x_img_px, const int& pyramid_level);
    /**< transform input pixel distance to level 0 from input level */
    static double pyramid_scale_to_level0(const double& px, const int& pyramid_level);

    /**< get number of horizontal pixels (number of columns) */
    const unsigned int& get_width_px() const {return _width_px;}
    /**< get number of vertical pixels (number of rows) */
    const unsigned int& get_height_px() const {return _height_px;}
    /**< get field of view */
    const double& get_vfov_deg() const {return _vfov_deg;}
    /**< get focal length [pixels] */
    const double& get_f_px() const {return _f_px;}
    /**< get horizontal principal point location (from left side) */
    const double& get_cx_px() const {return _cx_px;}
    /**< get vertical principal point location (from top) */
    const double& get_cy_px() const {return _cy_px;}

    /**< get true transformation from body to camera based on mass ratio */
    virtual ang::speu_rodrigues get_Gq_bc_truth(const double& ratio) const = 0;
    /**< get estimated translation from body to camera based on mass ratio */
    virtual ang::speu_rodrigues get_Gq_bc_est(const double& ratio) const = 0;

    /**< describe camera in stream */
    virtual void create_text(std::ostream& Ostream) const = 0;
}; // closes class camera_base

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS CAMERA
// ============
// ============

class ACFT_API camera : public camera_base {
private:
    /**< weak pointer to aircraft inertia */
    const acft::iner* const _Piner;

    /**< standard deviation in each direction when measuring location of camera frame */
    double _sigma_T_rcb_m;
    /**< standard deviation of each real Euler angle (not error) when fixing camera to aircraft */
    double _sigma_ypb_bc_deg;
    /**< standard deviation of each Euler angles when measuring orientation of camera frame */
    double _sigma_euler_bc_deg;

    /**< true translation from aircraft reference (wing trailing edge) "r" to camera "c" viewed in body "b" */
    Eigen::Array3d _Trcb_m_truth;
    /**< estimated translation from aircraft reference (wing trailing edge) "r" to camera "c" viewed in body "b" */
    Eigen::Array3d _Trcb_m_est;

    /**< true euler angles from body to camera frames */
    ang::euler _euler_bc_truth;
    /**< estimated euler angles from body to camera frames */
    ang::euler _euler_bc_est;

    /**< true transformation from aircraft reference (wing trailing edge) "r" to camera "c" */
    ang::speu_rodrigues _Gq_rc_truth;
    /**< estimated transformation from aircraft reference (wing trailing edge) "r" to camera "c" */
    ang::speu_rodrigues _Gq_rc_est;
public:
    /**< default constructor */
    camera() = delete;
    /**< constructor based on seeder, images width and height, vertical field of view, and aircraft inertial properties.
     * Flag rotate is true if image is rotated 90 [deg] (x right, y down), false otherwise (x up, y right).
     * Flag console true shows summary in console, false does not. */
    camera(math::seeder& Oseeder, const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg, const acft::iner& Oiner, bool flag_rotate, bool flag_console);
    /**< copy constructor */
    camera(const camera&) = delete;
    /**< move constructor */
    camera(camera&&) = delete;
    /**< destructor */
    ~camera() override = default;
    /**< copy assignment */
    camera& operator=(const camera&) = delete;
    /**< move assignment */
    camera& operator=(camera&&) = delete;

    /**< get standard deviation in each direction when measuring location of camera */
    const double& get_sigma_T_rcb_m() const {return _sigma_T_rcb_m;}
    /**< get standard deviation of each real Euler angle (not error) when fixing camera to aircraft */
    const double& get_sigma_ypb_bc_deg() const {return _sigma_ypb_bc_deg;}
    /**< get standard deviation of each Euler angles when measuring orientation of camera frame */
    const double& get_sigma_euler_bc_deg() const {return _sigma_euler_bc_deg;}

    /**< get true transformation from body to camera based on mass ratio */
    ang::speu_rodrigues get_Gq_bc_truth(const double& ratio) const override;
    /**< get estimated translation from body to camera based on mass ratio */
    ang::speu_rodrigues get_Gq_bc_est(const double& ratio) const override;
    /**< describe camera in stream */
    void create_text(std::ostream& Ostream) const override;
}; // closes class camera

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS CAMERA_TEST
// =================
// =================

class ACFT_API camera_test : public camera_base {
private:
    /**< true and estimated (they are the same) transformation between body and camera */
    const ang::speu_rodrigues* const _PGq_bc;
public:
    /**< default constructor */
    camera_test() = delete;
    /**< constructor based on images width and height, together with vertical field of view. Also includes transformation
     * between body and camera frames. Flag console true shows summary in console, false does not. */
    camera_test(const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg, ang::speu_rodrigues* PGq_bc, bool flag_console);
    /**< copy constructor */
    camera_test(const camera_test&) = delete;
    /**< move constructor */
    camera_test(camera_test&&) = delete;
    /**< destructor */
    ~camera_test() override;
    /**< copy assignment */
    camera_test& operator=(const camera_test&) = delete;
    /**< move assignment */
    camera_test& operator=(camera_test&&) = delete;

    /**< get true transformation from body to camera based on mass ratio */
    ang::speu_rodrigues get_Gq_bc_truth(const double& ratio_NOT_USED) const override;
    /**< get estimated translation from body to camera based on mass ratio */
    ang::speu_rodrigues get_Gq_bc_est(const double& ratio_NOT_USED) const override;
    /**< describe camera in stream */
    void create_text(std::ostream& Ostream) const override;
}; // closes class camera_test

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

} // closes namespace sens

#endif







