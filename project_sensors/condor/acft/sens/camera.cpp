#include "camera.h"
#include "../acft/iner.h"
#include "math/logic/constant.h"
#include "math/logic/seeder.h"
#include "ang/rotate/rodrigues.h"
#include "ang/transform/speu_rodrigues.h"
#include <iostream>
#include <random>

// CLASS CAMERA_BASE
// =================
// =================

sens::camera_base::camera_base(const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg)
: _width_px(width_px), _height_px(height_px), _vfov_deg(vfov_deg), _f_px(height_px / (2.0 * tan(0.5 * vfov_deg * math::constant::D2R()))),
  _cx_px(0.5 * width_px - 0.5), _cy_px(0.5 * height_px - 0.5) {
}
/* constructor based on images width and height, together with vertical field of view. */

Eigen::Vector3d sens::camera_base::convert_img2uscrs(const Eigen::Vector2d& x_img_px) const {
    return this->convert_img2uscrs(x_img_px[0], x_img_px[1]);
}
/* projection from image (pixels) to unit sphere camera coordinates */

Eigen::Vector3d sens::camera_base::convert_img2uscrs(const double& x_img_pxi, const double& x_img_pxii) const {
    Eigen::Vector3d xyz;
    xyz[0] = (x_img_pxi - _cx_px) / _f_px;
    xyz[1] = (x_img_pxii - _cy_px) / _f_px;
    xyz[2] = 1.0;
    return xyz.normalized();
}
/* projection from image (pixels) to unit sphere camera coordinates */

Eigen::Vector3d sens::camera_base::convert_img2upcrs(const Eigen::Vector2d& x_img_px) const {
    return this->convert_img2upcrs(x_img_px[0], x_img_px[1]);
}
/* projection from image (pixels) to unit plane camera coordinates */

Eigen::Vector3d sens::camera_base::convert_img2upcrs(const double& x_img_pxi, const double& x_img_pxii) const {
    return {(x_img_pxi - _cx_px) / _f_px, (x_img_pxii - _cy_px) / _f_px, 1.0};
}
/* projection from image (pixels) to unit plane camera coordinates */

Eigen::Vector2d sens::camera_base::convert_crs2img(const Eigen::Vector3d& x_crs_m) const {
    return this->convert_upcrs2img(Eigen::Vector3d(x_crs_m(0)/x_crs_m(2), x_crs_m(1)/x_crs_m(2), 1.0));
}
/* transforms cartesian coordinates from camera reference system (CRS) to image reference system (IMG).
Unit sphere or unit plane camera reference systems (USCRS or UPCRS) are also acceptable inputs. */

Eigen::Vector2d sens::camera_base::convert_upcrs2img(const Eigen::Vector3d& x_crs_up) const {
    // third input coordinate is always one and not employed
    return {_f_px * x_crs_up(0) + _cx_px, _f_px * x_crs_up(1) + _cy_px};
}
/* transforms cartesian coordinates from unit plane (3rd coord == 1) camera reference system (UPCRS) to image reference system (IMG) */

Eigen::Vector3d sens::camera_base::convert_uscrs2upcrs(const Eigen::Vector3d& x_crs_us) {
    return x_crs_us / x_crs_us(2);
}
/* projection from unit sphere camera reference system (USCRS) to unit plane camera reference system (UPCRS) */

Eigen::Vector3d sens::camera_base::convert_upcrs2uscrs(const Eigen::Vector3d& x_crs_up) {
    return x_crs_up.normalized();
}
/* projection from unit plane camera reference system (UPCRS) to unit sphere camera reference system (USCRS) */

bool sens::camera_base::is_in_frame(const Eigen::Vector2i& x_img_px, int boundary_px) const {
    if ((x_img_px[0] >= boundary_px) && (x_img_px[0] < (_width_px  - boundary_px)) &&
        (x_img_px[1] >= boundary_px) && (x_img_px[1] < (_height_px - boundary_px))) {
        return true;
    }
    return false;
}
/* returns true if input pixel coordinates are inside frame and not closer than input boundary to frame border */

bool sens::camera_base::is_in_frame(const double& x_img_pxi, const double& x_img_pxii, int boundary_px) const {
    if ((x_img_pxi  >= boundary_px) && (x_img_pxi  < (_width_px  - boundary_px)) &&
        (x_img_pxii >= boundary_px) && (x_img_pxii < (_height_px - boundary_px))) {
        return true;
    }
    return false;
}
/* returns true if input pixel coordinates are inside frame and not closer than input boundary to frame border */

bool sens::camera_base::is_in_frame(const Eigen::Vector2i& x_img_px, int boundary_px, int level) const {
    if ((x_img_px[0] >= boundary_px) && (x_img_px[0] < _width_px  / sens::camera::pyramid_scale(level) - boundary_px) &&
        (x_img_px[1] >= boundary_px) && (x_img_px[1] < _height_px / sens::camera::pyramid_scale(level) - boundary_px)) {
        return true;
    }
    return false;
}
/* returns true if input pixel coordinates (at input pyramid level) are inside frame and not closer than input boundary to frame border */

bool sens::camera_base::is_in_frame(const double& x_img_pxi, const double& x_img_pxii, int boundary_px, int level) const {
    if ((x_img_pxi  >= boundary_px) && (x_img_pxi  < _width_px  / sens::camera::pyramid_scale(level) - boundary_px) &&
        (x_img_pxii >= boundary_px) && (x_img_pxii < _height_px / sens::camera::pyramid_scale(level) - boundary_px)) {
        return true;
    }
    return false;
}
/* returns true if input pixel coordinates (at input pyramid level) are inside frame and not closer than input boundary to frame border */

Eigen::Vector2d sens::camera_base::pyramid_scale_from_level0(const Eigen::Vector2d& x_img_px, const int& level) {
    return x_img_px / sens::camera_base::pyramid_scale(level);
}
/* transform input image (pixels) from level 0 to input level */

double sens::camera_base::pyramid_scale_from_level0(const double& px, const int& level) {
    return px / sens::camera_base::pyramid_scale(level);
}
/* transform input pixel distance from level 0 to input level */

Eigen::Vector2d sens::camera_base::pyramid_scale_to_level0(const Eigen::Vector2d& x_img_px, const int& level) {
    return x_img_px * sens::camera_base::pyramid_scale(level);
}
/* transform input image (pixels) to level 0 from input level */

double sens::camera_base::pyramid_scale_to_level0(const double& px, const int& level) {
    return px * sens::camera_base::pyramid_scale(level);
}
/* transform input pixel distance to level 0 from input level */

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS CAMERA
// ============
// ============

sens::camera::camera(math::seeder& Oseeder, const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg, const acft::iner& Oiner, bool flag_rotate, bool flag_console)
: camera_base(width_px, height_px, vfov_deg), _Piner(&Oiner), _sigma_T_rcb_m(0.002), _sigma_ypb_bc_deg(0.1), _sigma_euler_bc_deg(0.01) {

    // seed generator
    std::ranlux24_base Ogen(Oseeder.provide_seed(math::seeder::seeder_cam));

    // standard normal distribution
    std::normal_distribution<double> Odist(0.,1.);

    // camera is located 50 [cm] forward and 12 [cm] below wing trailing edge
    _Trcb_m_truth << 0.5, 0.0, 0.12;

    // estimated camera position with respect to reference (wing) in body
    _Trcb_m_est << _Trcb_m_truth(0) + _sigma_T_rcb_m * Odist(Ogen),
                   _Trcb_m_truth(1) + _sigma_T_rcb_m * Odist(Ogen),
                   _Trcb_m_truth(2) + _sigma_T_rcb_m * Odist(Ogen);
    //_Trcb_m_est = _Trcb_m_truth; ///////////////////////////////////////////////////////////////////////////////////////

    // true Euler angles from body to camera
    if (flag_rotate == true) {
        _euler_bc_truth.set_yaw_rad((90.0 + _sigma_ypb_bc_deg * Odist(Ogen)) * math::constant::D2R());
    }
    else {
        _euler_bc_truth.set_yaw_rad(_sigma_ypb_bc_deg * Odist(Ogen) * math::constant::D2R());
    }
    _euler_bc_truth.set_pitch_rad((0.0 + _sigma_ypb_bc_deg * Odist(Ogen)) * math::constant::D2R());
    _euler_bc_truth.set_bank_rad ((0.0 + _sigma_ypb_bc_deg * Odist(Ogen)) * math::constant::D2R());

    // estimated Euler angles from body to camera
    _euler_bc_est.set_yaw_rad  (_euler_bc_truth.get_yaw_rad()   + _sigma_euler_bc_deg * Odist(Ogen) * math::constant::D2R());
    _euler_bc_est.set_pitch_rad(_euler_bc_truth.get_pitch_rad() + _sigma_euler_bc_deg * Odist(Ogen) * math::constant::D2R());
    _euler_bc_est.set_bank_rad (_euler_bc_truth.get_bank_rad()  + _sigma_euler_bc_deg * Odist(Ogen) * math::constant::D2R());
    //_euler_bc_est = _euler_bc_truth; /////////////////////////////////////////////////////////////////////////////////////

    // true and estimated transformations from aircraft reference (wing trailing edge) "r" to camera "c" */
    _Gq_rc_truth.set(ang::rodrigues(_euler_bc_truth), _Trcb_m_truth);
    _Gq_rc_est.set(ang::rodrigues(_euler_bc_est), _Trcb_m_est);

    if (flag_console == true) {
        this->create_text(std::cout);
    }
}
/* constructor based on seeder, images width and height, vertical field of view, and aircraft inertial properties.
 * Flag rotate is true if image is rotated 90 [deg] (x right, y down), false otherwise (x up, y right).
 * Flag console true shows summary in console, false does not. */

ang::speu_rodrigues sens::camera::get_Gq_bc_truth(const double& ratio) const {
    return ang::speu_rodrigues(_Gq_rc_truth.get_rodrigues(), _Trcb_m_truth - _Piner->get_Trbb_m(ratio));
}
/* get true transformation from body to camera based on mass ratio */

ang::speu_rodrigues sens::camera::get_Gq_bc_est(const double& ratio) const {
    return ang::speu_rodrigues(_Gq_rc_est.get_rodrigues(), _Trcb_m_est - _Piner->get_Trbb_m(ratio));
}
/* get estimated translation from body to camera based on mass ratio */

void sens::camera::create_text(std::ostream& Ostream) const {
    Ostream << std::endl << "CAMERA:" << std::endl << std::endl;
    Ostream << "Trcb truth [m]:      "
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_truth(0)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_truth(1)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_truth(2) << std::endl;
    Ostream << "Trcb est [m]:        "
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_est(0)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_est(1)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _Trcb_m_est(2) << std::endl;
    Ostream << "euler bc truth [deg]:"
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_truth.get_yaw_rad()   * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_truth.get_pitch_rad() * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_truth.get_bank_rad()  * math::constant::R2D() << std::endl;
    Ostream << "euler bc est [deg]:  "
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_est.get_yaw_rad()   * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_est.get_pitch_rad() * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _euler_bc_est.get_bank_rad()  * math::constant::R2D() << std::endl;
}
/* describe platform model in stream */

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS CAMERA_TEST
// =================
// =================

sens::camera_test::camera_test(const unsigned int& width_px, const unsigned int& height_px, const double& vfov_deg, ang::speu_rodrigues* PGq_bc, bool flag_console)
: camera_base(width_px, height_px, vfov_deg), _PGq_bc(PGq_bc) {
    if (flag_console == true) {
        this->create_text(std::cout);
    }
}
/* constructor based on images width and height, together with vertical field of view. Also includes transformation
 * between body and camera frames. Flag console true shows summary in console, false does not. */

sens::camera_test::~camera_test() {
    delete _PGq_bc;
}
/* destructor */

ang::speu_rodrigues sens::camera_test::get_Gq_bc_truth(const double& ratio_NOT_USED) const {
    return *_PGq_bc;
}
/* get true transformation from body to camera based on mass ratio */

ang::speu_rodrigues sens::camera_test::get_Gq_bc_est(const double& ratio_NOT_USED) const {
    return *_PGq_bc;
}
/* get estimated translation from body to camera based on mass ratio */

void sens::camera_test::create_text(std::ostream& Ostream) const {
    Ostream << std::endl << "CAMERA:" << std::endl << std::endl;
    Ostream << "Tbcb truth & est [m]:      "
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _PGq_bc->get_T()(0)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _PGq_bc->get_T()(1)
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << _PGq_bc->get_T()(2) << std::endl;
    ang::euler euler_bc(_PGq_bc->get_rodrigues());
    Ostream << "euler bc truth & est [deg]:"
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << euler_bc.get_yaw_rad()   * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << euler_bc.get_pitch_rad() * math::constant::R2D()
            << std::fixed << std::setw(9) << std::setprecision(4) << std::showpos << euler_bc.get_bank_rad()  * math::constant::R2D() << std::endl;
}
/* describe platform model in stream */

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



















