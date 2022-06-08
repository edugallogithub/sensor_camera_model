#ifndef ACFT_TRJ_TRUTH
#define ACFT_TRJ_TRUTH

#include "../acft.h"
#include "trj.h"
#include "ang/rotate/rodrigues.h"
#include "ang/rotate/so3_tangent.h"
#include "ang/transform/dual.h"
#include "ang/transform/se3_tangent.h"
#include "env/coord.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <iostream>

namespace st {

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS STATE ST_TRUTH
// ====================
// ====================

class ACFT_API st_truth {
private:
    /**< time */
    double _t_sec;
    /**< geodetic coordinates */
    env::geodetic_coord _x_gdt_b_rad_m;
    /**< BFS absolute velocity */
    Eigen::Vector3d _v_b_mps;
    /**< BFS to NED quaternion */
    ang::rodrigues _q_nb;
    /**< angular speed of BFS over IRS expressed in BFS */
    ang::so3_tangent _w_ibb_rps;
    /**< aircraft mass */
    double _m_kg;
    /**< twist from ECEF to BFS */
    ang::se3_tangent _xi_ebb_mrps;
    /**< unit dual quaternion from ECEF to BFS */
    ang::dual _z_eb;
    /**< control parameters (throttle, elevator, ailerons, rudder) */
    Eigen::Array4d _delta_control;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_truth() = default;
    /**< copy constructor */
    st_truth(const st_truth&) = default;
    /**< move constructor */
    st_truth(st_truth&&) = default;
    /**< destructor */
    ~st_truth() = default;
    /**< copy assignment */
    st_truth& operator=(const st_truth&) = default;
    /**< move assignment */
    st_truth& operator=(st_truth&&) = default;

    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}
    /**< get geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m() {return _x_gdt_b_rad_m;}
    const env::geodetic_coord& get_x_gdt_b_rad_m() const {return _x_gdt_b_rad_m;}
    /**< get BFS absolute velocity */
    Eigen::Vector3d& get_v_b_mps() {return _v_b_mps;}
    const Eigen::Vector3d& get_v_b_mps() const {return _v_b_mps;}
    /**< get BFS to NED quaternion */
    ang::rodrigues& get_q_nb() {return _q_nb;}
    const ang::rodrigues& get_q_nb() const {return _q_nb;}
    /**< get angular speed of BFS over IRS expressed in BFS */
    ang::so3_tangent& get_w_ibb_rps() {return _w_ibb_rps;}
    const ang::so3_tangent& get_w_ibb_rps() const {return _w_ibb_rps;}
    /**< get aircraft mass */
    double& get_m_kg() {return _m_kg;}
    const double& get_m_kg() const {return _m_kg;}
    /**< get twist from ECEF to BFS */
    ang::se3_tangent& get_xi_ebb_mrps() {return _xi_ebb_mrps;}
    const ang::se3_tangent& get_xi_ebb_mrps() const {return _xi_ebb_mrps;}
    /**< get unit dual quaternion from ECEF to BFS */
    ang::dual& get_z_eb() {return _z_eb;}
    const ang::dual& get_z_eb() const {return _z_eb;}
    /**< get control parameters (throttle, elevator, ailerons, rudder) */
    Eigen::Array4d& get_delta_control() {return _delta_control;}
    const Eigen::Array4d& get_delta_control() const {return _delta_control;}    
}; // closes class st_truth

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS STATE DIFFERENTIAL ST_DIFF
// ================================
// ================================

class ACFT_API st_diff {
private:
    /**< time derivative of geodetic coordinates */
    Eigen::Array3d _dx_gdt_rad_m_dt;
    /**< time derivative of BFS absolute velocity */
    Eigen::Vector3d _dv_b_mps_dt;
    /**< time derivative of BFS to NED quaternion */
    ang::quat _dq_nb_dt;
    /**< angular speed from BFS to NED in BFS (perturbation for quaternion integration) */
    ang::so3_tangent _w_nbb_rps;
    /**< time derivative of angular speed of BFS over IRS expressed in BFS */
    Eigen::Vector3d _dw_ibb_rps_dt;
    /**< time derivative of aircraft mass */
    double _dm_kg_dt;
    /**< twist from ECEF to BFS in BFS (perturbation for dual quaternion integration) */
    ang::se3_tangent _xi_ebb_mrps;
    /**< time derivative of twist from ECEF to BFS viewed in BFS */
    Eigen::Vector6d _dxi_ebb_mrps_dt;
    /**< time derivative control parameters (throttle, elevator, ailerons, rudder) */
    Eigen::Array4d _ddelta_control_dt;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_diff() = default;
    /**< copy constructor */
    st_diff(const st_diff&) = delete;
    /**< move constructor */
    st_diff(st_diff&&) = delete;
    /**< destructor */
    ~st_diff() = default;
    /**< copy assignment */
    st_diff& operator=(const st_diff&) = delete;
    /**< move assignment */
    st_diff& operator=(st_diff&&) = delete;

    /**< get time derivative of geodetic coordinates */
    Eigen::Array3d& get_dx_gdt_rad_m_dt() {return _dx_gdt_rad_m_dt;}
    const Eigen::Array3d& get_dx_gdt_rad_m_dt() const {return _dx_gdt_rad_m_dt;}
    /**< get time derivative of BFS absolute velocity */
    Eigen::Vector3d& get_dv_b_mps_dt() {return _dv_b_mps_dt;}
    const Eigen::Vector3d& get_dv_b_mps_dt() const {return _dv_b_mps_dt;}
    /**< get time derivative of BFS to NED quaternion */
    ang::quat& get_dq_nb_dt() {return _dq_nb_dt;}
    const ang::quat& get_dq_nb_dt() const {return _dq_nb_dt;}
    /**< get angular speed from BFS to NED in BFS (perturbation for quaternion integration) */
    ang::so3_tangent& get_w_nbb_rps() {return _w_nbb_rps;}
    const ang::so3_tangent& get_w_nbb_rps() const {return _w_nbb_rps;}
    /**< get time derivative of angular speed of BFS over IRS expressed in BFS */
    Eigen::Vector3d& get_dw_ibb_rps_dt() {return _dw_ibb_rps_dt;}
    const Eigen::Vector3d& get_dw_ibb_rps_dt() const {return _dw_ibb_rps_dt;}
    /**< get time derivative of aircraft mass */
    double& get_dm_kg_dt() {return _dm_kg_dt;}
    const double& get_dm_kg_dt() const {return _dm_kg_dt;}
    /**< get twist from ECEF to BFS in BFS (perturbation for dual quaternion integration) */
    ang::se3_tangent& get_xi_ebb_mrps() {return _xi_ebb_mrps;}
    const ang::se3_tangent& get_xi_ebb_mrps() const {return _xi_ebb_mrps;}
    /**< get time derivative of twist from ECEF to BFS viewed in BFS */
    Eigen::Vector6d& get_dxi_ebb_mrps_dt() {return _dxi_ebb_mrps_dt;}
    const Eigen::Vector6d& get_dxi_ebb_mrps_dt() const {return _dxi_ebb_mrps_dt;}
    /**< get time derivative control parameters (throttle, elevator, ailerons, rudder) */
    Eigen::Array4d& get_ddelta_control_dt() {return _ddelta_control_dt;}
    const Eigen::Array4d& get_ddelta_control_dt() const {return _ddelta_control_dt;}
}; // closes class st_diff

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS TRUTH TRAJECTORY TRJ_TRUTH
// ================================
// ================================

class ACFT_API trj_truth : public trj {
private:
    /**< vector of truth states */
    std::vector<st::st_truth,Eigen::aligned_allocator<st::st_truth>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_truth() = delete;
    /**< constructor based on time separation between consecutive samples and number of operations */
    trj_truth(const double& Deltat_sec, const unsigned int& nel, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, 1), _Vst(nel) {
    }
    /**< copy constructor */
    trj_truth(const trj_truth&) = delete;
    /**< move constructor */
    trj_truth(trj_truth&&) = delete;
    /**< destructor */
    ~trj_truth() = default;
    /**< copy assignment */
    trj_truth& operator=(const trj_truth&) = delete;
    /**< move assignment */
    trj_truth& operator=(trj_truth&&) = delete;

    /**< get vector of truth states */
	std::vector<st::st_truth,Eigen::aligned_allocator<st::st_truth>>& operator()() {return _Vst;}
    const std::vector<st::st_truth,Eigen::aligned_allocator<st::st_truth>>& operator()() const {return _Vst;}
    /**< get vector of truth states */
    std::vector<st::st_truth,Eigen::aligned_allocator<st::st_truth>>& get() {return _Vst;}
    const std::vector<st::st_truth,Eigen::aligned_allocator<st::st_truth>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_truth

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif















