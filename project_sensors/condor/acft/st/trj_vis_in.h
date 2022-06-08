#ifndef ACFT_TRJ_VIS_IN
#define ACFT_TRJ_VIS_IN

#include "../acft.h"
#include "trj.h"
#include "env/coord.h"
#include "ang/rotate/rotv.h"
#include "ang/rotate/euler.h"
#include "ang/rotate/rodrigues.h"
#include "ang/transform/speu_rodrigues.h"
#include "ang/transform/trfv.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

namespace st {

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS IMAGE INPUT STATE ST_VIS_IN
// =================================
// =================================

class ACFT_API st_vis_in {
private:
    /**< image id (starts in 0 and grows) */
    int _id;
    /**< time */
    double _t_sec;

    /**< center of gravity (not camera center) true geodetic coordinates */
    env::geodetic_coord _x_gdt_b_rad_m_truth;
    /**< true unit quaternion from NED to body */
    ang::rodrigues _q_nb_truth;
    /**< true transformation from Earth to body */
    ang::speu_rodrigues _Gq_eb_truth;
    /**< true NED absolute velocity */
    Eigen::Vector3d _v_n_mps_truth;
    /**< true mass ratio */
    double _ratio_truth;
    /**< true transform vector between previous (reference) and current image body poses */
    ang::trfv _tau_b_Deltarc_truth;

    /**< center of gravity (not camera center) estimated geodetic coordinates */
    env::geodetic_coord _x_gdt_b_rad_m_est;
    /**< estimated unit quaternion from NED to body */
    ang::rodrigues _q_nb_est;
    /**< estimated transformation from Earth to body */
    ang::speu_rodrigues _Gq_eb_est;
    /**< estimated NED absolute velocity */
    Eigen::Vector3d _v_n_mps_est;
    /**< estimated mass ratio */
    double _ratio_est;
    /**< estimated transform vector between previous (reference) and current image body poses */
    ang::trfv _tau_b_Deltarc_est;

    /**< estimated airspeed vector in body */
    Eigen::Vector3d _vtas_b_mps_est;
    /**< estimated rotation vector from previous to current frame viewed in previous frame based on
     * estimated unit quaternions of both previous and current frames */
    ang::rotv _rotv_bprev_Deltarc_est;
    /**< estimated translation due to airspeed exclusively (no wind) between previous and current
     * frame viewed in previous NED frame */
    Eigen::Vector3d _Delta_xtas_nprev_m_est;
    /**< estimated translation due to airspeed exclusively (no wind) between previous and current
     * frame viewed in previous body frame */
    Eigen::Vector3d _Delta_xtas_bprev_m_est;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_vis_in() = default;
    /**< copy constructor */
    st_vis_in(const st_vis_in&) = default;
    /**< move constructor */
    st_vis_in(st_vis_in&&) = default;
    /**< destructor */
    ~st_vis_in() = default;
    /**< copy assignment */
    st_vis_in& operator=(const st_vis_in&) = default;
    /**< move assignment */
    st_vis_in& operator=(st_vis_in&&) = default;

    /**< get image id (starts in 0 and grows) */
    int& get_id() {return _id;}
    const int& get_id() const {return _id;}
    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}

    /**< get center of gravity (not camera center) true geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m_truth() {return _x_gdt_b_rad_m_truth;}
    const env::geodetic_coord& get_x_gdt_b_rad_m_truth() const {return _x_gdt_b_rad_m_truth;}
    /**< get true unit quaternion from NED to body */
    ang::rodrigues& get_q_nb_truth() {return _q_nb_truth;}
    const ang::rodrigues& get_q_nb_truth() const {return _q_nb_truth;}
    /**< get true transformation from Earth to body */
    ang::speu_rodrigues& get_Gq_eb_truth() {return _Gq_eb_truth;}
    const ang::speu_rodrigues& get_Gq_eb_truth() const {return _Gq_eb_truth;}
    /**< get true NED absolute velocity */
    Eigen::Vector3d& get_v_n_mps_truth() {return _v_n_mps_truth;}
    const Eigen::Vector3d& get_v_n_mps_truth() const {return _v_n_mps_truth;}
    /**< get true mass ratio */
    double& get_ratio_truth() {return _ratio_truth;}
    const double& get_ratio_truth() const {return _ratio_truth;}
    /**< get true transform vector between previous (reference) and current image body poses */
    ang::trfv& get_tau_b_Deltarc_truth() {return _tau_b_Deltarc_truth;}
    const ang::trfv& get_tau_b_Deltarc_truth() const {return _tau_b_Deltarc_truth;}

    /**< get center of gravity (not camera center) estimated geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m_est() {return _x_gdt_b_rad_m_est;}
    const env::geodetic_coord& get_x_gdt_b_rad_m_est() const {return _x_gdt_b_rad_m_est;}
    /**< get estimated unit quaternion from NED to body */
    ang::rodrigues& get_q_nb_est() {return _q_nb_est;}
    const ang::rodrigues& get_q_nb_est() const {return _q_nb_est;}
    /**< get estimated transformation from Earth to body */
    ang::speu_rodrigues& get_Gq_eb_est() {return _Gq_eb_est;}
    const ang::speu_rodrigues& get_Gq_eb_est() const {return _Gq_eb_est;}
    /**< get estimated NED absolute velocity */
    Eigen::Vector3d& get_v_n_mps_est() {return _v_n_mps_est;}
    const Eigen::Vector3d& get_v_n_mps_est() const {return _v_n_mps_est;}
    /**< get estimated mass ratio */
    double& get_ratio_est() {return _ratio_est;}
    const double& get_ratio_est() const {return _ratio_est;}
    /**< get estimated transform vector between previous (reference) and current image body poses */
    ang::trfv& get_tau_b_Deltarc_est() {return _tau_b_Deltarc_est;}
    const ang::trfv& get_tau_b_Deltarc_est() const {return _tau_b_Deltarc_est;}

    /**< get estimated airspeed vector in body */
    Eigen::Vector3d& get_vtas_b_mps_est() {return _vtas_b_mps_est;}
    const Eigen::Vector3d& get_vtas_b_mps_est() const {return _vtas_b_mps_est;}
    /**< get estimated rotation vector from previous to current frame viewed in previous frame based on
     * estimated unit quaternions of both previous and current frames */
    ang::rotv& get_rotv_bprev_Deltarc_est() {return _rotv_bprev_Deltarc_est;}
    const ang::rotv& get_rotv_bprev_Deltarc_est() const {return _rotv_bprev_Deltarc_est;}
    /**< estimated translation due to airspeed exclusively (no wind) between previous and current
     * frame viewed in previous NED frame */
    Eigen::Vector3d& get_Delta_xtas_nprev_m_est() {return _Delta_xtas_nprev_m_est;}
    const Eigen::Vector3d& get_Delta_xtas_nprev_m_est() const {return _Delta_xtas_nprev_m_est;}
    /**< get estimated translation due to airspeed exclusively (no wind) between previous and current
     * frame viewed in previous body frame */
    Eigen::Vector3d& get_Delta_xtas_bprev_m_est() {return _Delta_xtas_bprev_m_est;}
    const Eigen::Vector3d& get_Delta_xtas_bprev_m_est() const {return _Delta_xtas_bprev_m_est;}
}; // closes class st_vis_in

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS IMAGE INPUT TRAJECTORY TRJ_VIS_IN
// =======================================
// =======================================

class ACFT_API trj_vis_in : public trj{
private:
    /**< vector of sensor input states */
    std::vector<st::st_vis_in,Eigen::aligned_allocator<st::st_vis_in>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_vis_in() = delete;
    /**< constructor based on time separation between consecutive images, size of trajectory,
     * time separation between consecutive truth samples, and number of operations */
    trj_vis_in(const double& Deltat_sec, const unsigned int& nel, const double& Deltat_sec_truth, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, (int)(Deltat_sec / Deltat_sec_truth)), _Vst(nel) {
        if (std::remainder(_Deltat_sec, Deltat_sec_truth) > math::constant::EPS()) {throw std::runtime_error("Time separation between samples does not match.");}
    }
    /**< copy constructor */
    trj_vis_in(const trj_vis_in&) = delete;
    /**< move constructor */
    trj_vis_in(trj_vis_in&&) = delete;
    /**< destructor */
    ~trj_vis_in() override = default;
    /**< copy assignment */
    trj_vis_in& operator=(const trj_vis_in&) = delete;
    /**< move assignment */
    trj_vis_in& operator=(trj_vis_in&&) = delete;

    /**< get vector of input image states */
    std::vector<st::st_vis_in,Eigen::aligned_allocator<st::st_vis_in>>& operator()() {return _Vst;}
    const std::vector<st::st_vis_in,Eigen::aligned_allocator<st::st_vis_in>>& operator()() const {return _Vst;}
    /**< get vector of input image states */
    std::vector<st::st_vis_in,Eigen::aligned_allocator<st::st_vis_in>>& get() {return _Vst;}
    const std::vector<st::st_vis_in,Eigen::aligned_allocator<st::st_vis_in>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_vis_in

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif
