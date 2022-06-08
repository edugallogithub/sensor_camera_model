#ifndef ACFT_TRJ_VIS_OUT
#define ACFT_TRJ_VIS_OUT

#include "../acft.h"
#include "../logic.h"
#include "trj.h"
#include "ang/transform/speu_rodrigues.h"
#include "ang/transform/trfv.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

namespace st {

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS IMAGE INPUT STATE ST_VIS_OUT
// ==================================
// ==================================

class ACFT_API st_vis_out {
private:
    /**< image id (starts in 0 and grows) */
    int _id;
    /**< time */
    double _t_sec;

    /**< center of gravity (not camera center) geodetic coordinates */
    env::geodetic_coord _x_gdt_b_rad_m;
    /**< unit quaternion from NED to body */
    ang::rodrigues _q_nb;
    /**< transformation from Earth to body */
    ang::speu_rodrigues _Gq_eb;
    /**< translation from initial world frame to body viewed in initial world */
    Eigen::Vector3d _x_w0bw0_m;
    /**< transformation from initial world to camera */
    ang::speu_rodrigues _Gq_w0c;
    /**< transform vector between previous (reference) and current body transformations */
    ang::trfv _tau_b_Deltarc;
    /**< transform vector between previous (reference) and current camera transformations */
    ang::trfv _tau_c_Deltarc;
    /**< NED absolute velocity (computed from differences, so noisy) */
    Eigen::Vector3d _v_n_mps;
    /**< NED absolute velocity (computed from differences, so noisy --> need to low pass filter) */
    Eigen::Vector3d _v_n_mps_smooth;
    /**< NED low frequency wind field */
    Eigen::Vector3d _vlf_n_mps;

    /**< transformation from initial world to camera (image align sparse output, only useful for plotting) */
    ang::speu_rodrigues _Gq_w0c_img_align_sparse;
    /**< transform vector between previous (reference) and current body transformations (image align sparse output, only useful for plotting) */
    ang::trfv _tau_b_Deltarc_img_align_sparse;

    /**< smooth (low pass filter) geometric altitude obtained by navigation filter */
    double _h_m_nav_out_smooth;
    /**< smooth (low pass filter) geometric altitude obtained by vision */
    double _h_m_vis_out_smooth;
    /**< difference between vision and inertial navigation smoothed estimated altitudes */
    double _h_m_vis_nav_smooth_diff;
    /**< time derivative of difference between vision and inertial navigation smoothed estimated altitudes computed over 100 images */
    double _hdot_vis_nav_diff_msec_100img_before;

    /**< smooth (low pass filter) body pitch obtained by navigation filter */
    double _theta_rad_nav_out_smooth;
    /**< smooth (low pass filter) body pitch obtained by vision */
    double _theta_rad_vis_out_smooth;
    /**< difference between vision and inertial navigation smoothed estimated body pitch */
    double _theta_deg_vis_nav_smooth_diff;

    /**< smooth (low pass filter) body bank obtained by navigation filter */
    double _xi_rad_nav_out_smooth;
    /**< smooth (low pass filter) body bank obtained by vision */
    double _xi_rad_vis_out_smooth;
    /**< difference between vision and inertial navigation smoothed estimated body bank */
    double _xi_deg_vis_nav_smooth_diff;

    /**< algorithm current stage (1st frame, 2nd frame, default frame, relocalizing) */
    vis::logic::HANDLING_STAGE _handling_stage;

    /**< center of gravity geodetic coordinates as visual sensor input */
    env::geodetic_coord _x_gdt_b_rad_m_sensed;
    /**< center of gravity cartesian coordinates as visual sensor input */
    env::cartesian_coord _x_car_b_m_sensed;
    /**< NED absolute velocity as visual sensor input */
    Eigen::Vector3d _v_n_mps_sensed;
    /**< ECEF absolute velocity as visual sensor input */
    Eigen::Vector3d _v_e_mps_sensed;
    /**> squared standard deviation of the NED absolute velocity with respect to its mean or average */
    Eigen::Vector3d _v_n_mps2_var;
    /**> squared standard deviation of the ECEF absolute velocity with respect to its mean or average */
    Eigen::Vector3d _v_e_mps2_var;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_vis_out() = default;
    /**< copy constructor */
    st_vis_out(const st_vis_out&) = default;
    /**< move constructor */
    st_vis_out(st_vis_out&&) = default;
    /**< destructor */
    ~st_vis_out() = default;
    /**< copy assignment */
    st_vis_out& operator=(const st_vis_out&) = default;
    /**< move assignment */
    st_vis_out& operator=(st_vis_out&&) = default;

    /**< get image id (starts in 0 and grows) */
    int& get_id() {return _id;}
    const int& get_id() const {return _id;}
    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}

    /**< get center of gravity (not camera center) geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m() {return _x_gdt_b_rad_m;}
    const env::geodetic_coord& get_x_gdt_b_rad_m() const {return _x_gdt_b_rad_m;}
    /**< get unit quaternion from NED to body */
    ang::rodrigues& get_q_nb() {return _q_nb;}
    const ang::rodrigues& get_q_nb() const {return _q_nb;}
    /**< get transformation from Earth to body */
    ang::speu_rodrigues& get_Gq_eb() {return _Gq_eb;}
    const ang::speu_rodrigues& get_Gq_eb() const {return _Gq_eb;}
    /**< get translation from initial world frame to body viewed in initial world */
    Eigen::Vector3d& get_x_w0bw0_m() {return _x_w0bw0_m;}
    const Eigen::Vector3d& get_x_w0bw0_m() const {return _x_w0bw0_m;}
    /**< get transformation from initial world to camera */
    ang::speu_rodrigues& get_Gq_w0c() {return _Gq_w0c;}
    const ang::speu_rodrigues& get_Gq_w0c() const {return _Gq_w0c;}
    /**< get transform vector between previous (reference) and current body transformations */
    ang::trfv& get_tau_b_Deltarc() {return _tau_b_Deltarc;}
    const ang::trfv& get_tau_b_Deltarc() const {return _tau_b_Deltarc;}
    /**< get transform vector between previous (reference) and current camera transformations */
    ang::trfv& get_tau_c_Deltarc() {return _tau_c_Deltarc;}
    const ang::trfv& get_tau_c_Deltarc() const {return _tau_c_Deltarc;}
    /**< get NED absolute velocity */
    Eigen::Vector3d& get_v_n_mps() {return _v_n_mps;}
    const Eigen::Vector3d& get_v_n_mps() const {return _v_n_mps;}
    /**< get smoothed NED absolute velocity */
    Eigen::Vector3d& get_v_n_mps_smooth() {return _v_n_mps_smooth;}
    const Eigen::Vector3d& get_v_n_mps_smooth() const {return _v_n_mps_smooth;}
    /**< get NED low frequency wind field */
    Eigen::Vector3d& get_vlf_n_mps() {return _vlf_n_mps;}
    const Eigen::Vector3d& get_vlf_n_mps() const {return _vlf_n_mps;}

    /**< get transformation from initial world to camera (image align sparse output, only useful for plotting) */
    ang::speu_rodrigues& get_Gq_w0c_img_align_sparse() {return _Gq_w0c_img_align_sparse;}
    const ang::speu_rodrigues& get_Gq_w0c_img_align_sparse() const {return _Gq_w0c_img_align_sparse;}
    /**< get transform vector between previous (reference) and current body transformations (image align sparse output, only useful for plotting) */
    ang::trfv& get_tau_b_Deltarc_img_align_sparse() {return _tau_b_Deltarc_img_align_sparse;}
    const ang::trfv& get_tau_b_Deltarc_img_align_sparse() const {return _tau_b_Deltarc_img_align_sparse;}

    /**< get smooth (low pass filter) geometric altitude obtained by navigation filter */
    double& get_h_m_nav_out_smooth() {return _h_m_nav_out_smooth;}
    const double& get_h_m_nav_out_smooth() const {return _h_m_nav_out_smooth;}
    /**< get smooth (low pass filter) geometric altitude obtained by vision */
    double& get_h_m_vis_out_smooth() {return _h_m_vis_out_smooth;}
    const double& get_h_m_vis_out_smooth() const {return _h_m_vis_out_smooth;}
    /**< get difference between vision and inertial navigation smoothed estimated altitudes */
    double& get_h_m_vis_nav_smooth_diff() {return _h_m_vis_nav_smooth_diff;}
    const double& get_h_m_vis_nav_smooth_diff() const {return _h_m_vis_nav_smooth_diff;}
    /**< get time derivative of difference between vision and inertial navigation smoothed estimated altitudes computed over 100 images */
    double& get_hdot_vis_nav_diff_msec_100img_before() {return _hdot_vis_nav_diff_msec_100img_before;}
    const double& get_hdot_vis_nav_diff_msec_100img_before() const {return _hdot_vis_nav_diff_msec_100img_before;}

    /**< get smooth (low pass filter) body pitch obtained by navigation filter */
    double& get_theta_rad_nav_out_smooth() {return _theta_rad_nav_out_smooth;}
    const double& get_theta_rad_nav_out_smooth() const {return _theta_rad_nav_out_smooth;}
    /**< get smooth (low pass filter) body pitch obtained by vision */
    double& get_theta_rad_vis_out_smooth() {return _theta_rad_vis_out_smooth;}
    const double& get_theta_rad_vis_out_smooth() const {return _theta_rad_vis_out_smooth;}
    /**< get difference between vision and inertial navigation smoothed estimated body pitch */
    double& get_theta_deg_vis_nav_smooth_diff() {return _theta_deg_vis_nav_smooth_diff;}
    const double& get_theta_deg_vis_nav_smooth_diff() const {return _theta_deg_vis_nav_smooth_diff;}

    /**< get smooth (low pass filter) body bank obtained by navigation filter */
    double& get_xi_rad_nav_out_smooth() {return _xi_rad_nav_out_smooth;}
    const double& get_xi_rad_nav_out_smooth() const {return _xi_rad_nav_out_smooth;}
    /**< get smooth (low pass filter) body bank obtained by vision */
    double& get_xi_rad_vis_out_smooth() {return _xi_rad_vis_out_smooth;}
    const double& get_xi_rad_vis_out_smooth() const {return _xi_rad_vis_out_smooth;}
    /**< get difference between vision and inertial navigation smoothed estimated body bank */
    double& get_xi_deg_vis_nav_smooth_diff() {return _xi_deg_vis_nav_smooth_diff;}
    const double& get_xi_deg_vis_nav_smooth_diff() const {return _xi_deg_vis_nav_smooth_diff;}

    /**< get algorithm current stage (1st frame, 2nd frame, default frame, relocalizing) */
    vis::logic::HANDLING_STAGE& get_handling_stage() {return _handling_stage;}
    const vis::logic::HANDLING_STAGE& get_handling_stage() const {return _handling_stage;}

    /**< get center of gravity geodetic coordinates as visual sensor input */
    env::geodetic_coord& get_x_gdt_b_rad_m_sensed() {return _x_gdt_b_rad_m_sensed;}
    const env::geodetic_coord& get_x_gdt_b_rad_m_sensed() const {return _x_gdt_b_rad_m_sensed;}
    /**< get center of gravity cartesian coordinates as visual sensor input */
    env::cartesian_coord& get_x_car_b_m_sensed() {return _x_car_b_m_sensed;}
    const env::cartesian_coord& get_x_car_b_m_sensed() const {return _x_car_b_m_sensed;}
    /**< get NED absolute velocity as visual sensor input */
    Eigen::Vector3d& get_v_n_mps_sensed() {return _v_n_mps_sensed;}
    const Eigen::Vector3d& get_v_n_mps_sensed() const {return _v_n_mps_sensed;}
    /**< get ECEF absolute velocity as visual sensor input */
    Eigen::Vector3d& get_v_e_mps_sensed() {return _v_e_mps_sensed;}
    const Eigen::Vector3d& get_v_e_mps_sensed() const {return _v_e_mps_sensed;}
    /**> get squared standard deviation of the NED absolute velocity with respect to its mean or average */
    Eigen::Vector3d& get_v_n_mps2_var() {return _v_n_mps2_var;}
    const Eigen::Vector3d& get_v_n_mps2_var() const {return _v_n_mps2_var;}
    /**> get squared standard deviation of the ECEF absolute velocity with respect to its mean or average */
    Eigen::Vector3d& get_v_e_mps2_var() {return _v_e_mps2_var;}
    const Eigen::Vector3d& get_v_e_mps2_var() const {return _v_e_mps2_var;}
}; // closes class st_vis_out

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS IMAGE INPUT TRAJECTORY TRJ_VIS_OUT
// =======================================
// =======================================

class ACFT_API trj_vis_out : public trj{
private:
    /**< vector of sensor input states */
    std::vector<st::st_vis_out,Eigen::aligned_allocator<st::st_vis_out>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_vis_out() = delete;
    /**< constructor based on time separation between consecutive images, size of trajectory,
     * time separation between consecutive truth samples, and number of operations */
    trj_vis_out(const double& Deltat_sec, const unsigned int& nel, const double& Deltat_sec_truth, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, (int)(Deltat_sec / Deltat_sec_truth)), _Vst(nel) {
        if (std::remainder(_Deltat_sec, Deltat_sec_truth) > math::constant::EPS()) {throw std::runtime_error("Time separation between samples does not match.");}
    }
    /**< copy constructor */
    trj_vis_out(const trj_vis_out&) = delete;
    /**< move constructor */
    trj_vis_out(trj_vis_out&&) = delete;
    /**< destructor */
    ~trj_vis_out() override = default;
    /**< copy assignment */
    trj_vis_out& operator=(const trj_vis_out&) = delete;
    /**< move assignment */
    trj_vis_out& operator=(trj_vis_out&&) = delete;

    /**< get vector of input image states */
    std::vector<st::st_vis_out,Eigen::aligned_allocator<st::st_vis_out>>& operator()() {return _Vst;}
    const std::vector<st::st_vis_out,Eigen::aligned_allocator<st::st_vis_out>>& operator()() const {return _Vst;}
    /**< get vector of input image states */
    std::vector<st::st_vis_out,Eigen::aligned_allocator<st::st_vis_out>>& get() {return _Vst;}
    const std::vector<st::st_vis_out,Eigen::aligned_allocator<st::st_vis_out>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_vis_out

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif
