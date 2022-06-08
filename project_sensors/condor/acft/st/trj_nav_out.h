#ifndef ACFT_TRJ_NAV_OUT
#define ACFT_TRJ_NAV_OUT

#include "../acft.h"
#include "trj.h"
#include "ang/rotate/rodrigues.h"
#include "ang/rotate/euler.h"
#include "ang/rotate/so3_tangent.h"
#include "ang/transform/dual.h"
#include "ang/transform/se3_tangent.h"
#include "env/coord.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

namespace st {

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS NAVIGATION OUTPUT STATE ST_NAV_OUT
// ========================================
// ========================================

class ACFT_API st_nav_out {
private:
    /**< time */
    double _t_sec;

    /**< ===== ===== ===== Variables computed in air data filter ===== ===== ===== */
    /**< true airspeed */
    double _vtas_mps;
    /**< Euler angles from WFS to BFS */
    ang::euler _euler_wb;
    /**< atmospheric temperature */
    double _T_degK;
    /**< pressure altitude */
    double _Hp_m;
    /**< rate of climb */
    double _roc_mps;
    /**< temperature offset at mean sea level */
    double _DeltaT_degK;
    /**< true airspeed differential with time (only available from some air data filters) */
    double _vtas_dot_mps2;
    /**< euler angles from WFS to BFS differential with time (only available from some air data filters) */
    Eigen::Vector3d _euler_wb_dot_rps;

    /**< ===== ===== ===== Variables computed in attitude filter ===== ===== ===== */
    /**< BFS to NED quaternion */
    ang::rodrigues _q_nb;
    /**< ECEF to BFS quaternion */
    ang::rodrigues _q_eb;
    /**< angular velocity of BFS with respect to NED viewed in BFS */
    ang::so3_tangent _w_nbb_rps;
    /**< angular velocity of BFS with respect to ECEF viewed in BFS */
    ang::so3_tangent _w_ebb_rps;
    /**< angular velocity of BFS with respect to ECEF viewed in ECEF */
    ang::so3_tangent _w_ebe_rps;
    /**< gyroscope full error (everything except white noise)s */
    Eigen::Vector3d _E_gyr_rps;
    /**< magnetometer full error (everything except white noise) */
    Eigen::Vector3d _E_mag_nT;
    /**< magnetic field error (model minus real values) */
    Eigen::Vector3d _B_n_nT_dev;
    /**< body yaw angle covariance */
    double _covar_psi_rad2;
    /**< body pitch angle covariance */
    double _covar_theta_rad2;
    /**< body roll angle covariance */
    double _covar_xi_rad2;

    /**< ===== ===== ===== Variables computed in gps or position filter ===== ===== ===== */
    /**< geodetic coordinates of center of mass (body)*/
    env::geodetic_coord _x_gdt_b_rad_m;
    /**< cartesian coordinates of center of mass (body) */
    env::cartesian_coord _x_car_b_m;
    /**< NED absolute velocity */
    Eigen::Vector3d _v_n_mps;
    /**< ECEF absolute velocity */
    Eigen::Vector3d _v_e_mps;
    /**< BFS absolute velocity */
    Eigen::Vector3d _v_b_mps;
    /**< specific force */
    Eigen::Vector3d _f_ibb_mps2;
    /**< accelerometer full error (everything except white noise) */
    Eigen::Vector3d _E_acc_mps2;
    /**< pressure offset at mean sea level */
    double _Deltap_pa;
    /**< NED low frequency wind field */
    Eigen::Vector3d _vlf_n_mps;

    /**< ==== ===== ===== Fusion variables ===== ===== ===== */
    /**< unit dual quaternion from ECEF to BFS */
    ang::dual _z_eb;
    /**< twist from ECEF to BFS */
    ang::se3_tangent _xi_ebb_mrps;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_nav_out() = default;
    /**< copy constructor */
    st_nav_out(const st_nav_out&) = default;
    /**< move constructor */
    st_nav_out(st_nav_out&&) = default;
    /**< destructor */
    ~st_nav_out() = default;
    /**< copy assignment */
    st_nav_out& operator=(const st_nav_out&) = default;
    /**< move assignment */
    st_nav_out& operator=(st_nav_out&&) = default;

    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}

    /**< ===== ===== ===== Variables computed in air data filter ===== ===== ===== */
    /**< get true airspeed */
    double& get_vtas_mps() {return _vtas_mps;}
    const double& get_vtas_mps() const {return _vtas_mps;}
    /**< get Euler angles from WFS to BFS */
    ang::euler& get_euler_wb() {return _euler_wb;}
    const ang::euler& get_euler_wb() const {return _euler_wb;}
    /**< get atmospheric temperature */
    double& get_T_degK() {return _T_degK;}
    const double& get_T_degK() const {return _T_degK;}
    /**< get pressure altitude */
    double& get_Hp_m() {return _Hp_m;}
    const double& get_Hp_m() const {return _Hp_m;}
    /**< get rate of climb */
    double& get_roc_mps() {return _roc_mps;}
    const double& get_roc_mps() const {return _roc_mps;}
    /**< get temperature offset at mean sea level */
    double& get_DeltaT_degK() {return _DeltaT_degK;}
    const double& get_DeltaT_degK() const {return _DeltaT_degK;}
    /**< get true airspeed differential with time (only available from some air data filters) */
    double& get_vtas_dot_mps2() {return _vtas_dot_mps2;}
    const double& get_vtas_dot_mps2() const {return _vtas_dot_mps2;}
    /**< euler angles from WFS to BFS differential with time (only available from some air data filters) */
    Eigen::Vector3d& get_euler_wb_dot_rps() {return _euler_wb_dot_rps;}
    const Eigen::Vector3d& get_euler_wb_dot_rps() const {return _euler_wb_dot_rps;}

    /**< ===== ===== ===== Variables computed in attitude filter ===== ===== ===== */
    /**< get BFS to NED quaternion */
    ang::rodrigues& get_q_nb() {return _q_nb;}
    const ang::rodrigues& get_q_nb() const {return _q_nb;}
    /**< get ECEF to BFS quaternion */
    ang::rodrigues& get_q_eb() {return _q_eb;}
    const ang::rodrigues& get_q_eb() const {return _q_eb;}
    /**< get angular velocity of BFS with respect to NED viewed in BFS */
    ang::so3_tangent& get_w_nbb_rps() {return _w_nbb_rps;}
    const ang::so3_tangent& get_w_nbb_rps() const {return _w_nbb_rps;}
    /**< get angular velocity of BFS with respect to ECEF viewed in BFS */
    ang::so3_tangent& get_w_ebb_rps() {return _w_ebb_rps;}
    const ang::so3_tangent& get_w_ebb_rps() const {return _w_ebb_rps;}
    /**< get angular velocity of BFS with respect to ECEF viewed in ECEF */
    ang::so3_tangent& get_w_ebe_rps() {return _w_ebe_rps;}
    const ang::so3_tangent& get_w_ebe_rps() const {return _w_ebe_rps;}
    /**< get gyroscope full error (everything except white noise) */
    Eigen::Vector3d& get_E_gyr_rps() {return _E_gyr_rps;}
    const Eigen::Vector3d& get_E_gyr_rps() const {return _E_gyr_rps;}
    /**< get magnetometer full error (everything except white noise) */
    Eigen::Vector3d& get_E_mag_nT() {return _E_mag_nT;}
    const Eigen::Vector3d& get_E_mag_nT() const {return _E_mag_nT;}
    /**< get magnetic field error (model minus real values) */
    Eigen::Vector3d& get_B_n_nT_dev() {return _B_n_nT_dev;}
    const Eigen::Vector3d& get_B_n_nT_dev() const {return _B_n_nT_dev;}
    /**< body yaw angle covariance */
    double& get_covar_psi_rad2() {return _covar_psi_rad2;}
    const double& get_covar_psi_rad2() const {return _covar_psi_rad2;}
    /**< body pitch angle covariance */
    double& get_covar_theta_rad2() {return _covar_theta_rad2;}
    const double& get_covar_theta_rad2() const {return _covar_theta_rad2;}
    /**< body roll angle covariance */
    double& get_covar_xi_rad2() {return _covar_xi_rad2;}
    const double& get_covar_xi_rad2() const {return _covar_xi_rad2;}

    /**< ===== ===== ===== Variables computed in gps or position filter ===== ===== ===== */
    /**< get geodetic coordinates of center of mass (b) */
    env::geodetic_coord& get_x_gdt_b_rad_m() {return _x_gdt_b_rad_m;}
    const env::geodetic_coord& get_x_gdt_b_rad_m() const {return _x_gdt_b_rad_m;}
    /**< get cartesian coordinates of center of mass (b) */
    env::cartesian_coord& get_x_car_b_m() {return _x_car_b_m;}
    const env::cartesian_coord& get_x_car_b_m() const {return _x_car_b_m;}
    /**< get NED absolute velocity */
    Eigen::Vector3d& get_v_n_mps() {return _v_n_mps;}
    const Eigen::Vector3d& get_v_n_mps() const {return _v_n_mps;}
    /**< get ECEF absolute velocity */
    Eigen::Vector3d& get_v_e_mps() {return _v_e_mps;}
    const Eigen::Vector3d& get_v_e_mps() const {return _v_e_mps;}
    /**< get BFS absolute velocity */
    Eigen::Vector3d& get_v_b_mps() {return _v_b_mps;}
    const Eigen::Vector3d& get_v_b_mps() const {return _v_b_mps;}
    /**< get specific force */
    Eigen::Vector3d& get_f_ibb_mps2() {return _f_ibb_mps2;}
    const Eigen::Vector3d& get_f_ibb_mps2() const {return _f_ibb_mps2;}
    /**< get accelerometer full error (everything except white noise) */
    Eigen::Vector3d& get_E_acc_mps2() {return _E_acc_mps2;}
    const Eigen::Vector3d& get_E_acc_mps2() const {return _E_acc_mps2;}
    /**< get pressure offset at mean sea level */
    double& get_Deltap_pa() {return _Deltap_pa;}
    const double& get_Deltap_pa() const {return _Deltap_pa;}
    /**< get NED low frequency wind field */
    Eigen::Vector3d& get_vlf_n_mps() {return _vlf_n_mps;}
    const Eigen::Vector3d& get_vlf_n_mps() const {return _vlf_n_mps;}

    /**< ==== ===== ===== Fusion variables ===== ===== ===== */
    /**< get unit dual quaternion from ECEF to BFS */
    ang::dual& get_z_eb() {return _z_eb;}
    const ang::dual& get_z_eb() const {return _z_eb;}
    /**< get twist from ECEF to BFS */
    ang::se3_tangent& get_xi_ebb_mrps() {return _xi_ebb_mrps;}
    const ang::se3_tangent& get_xi_ebb_mrps() const {return _xi_ebb_mrps;}
}; // closes class st_nav_out

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS NAVIGATION OUTPUT TRAJECTORY TRJ_NAV_OUT
// ==============================================
// ==============================================

class ACFT_API trj_nav_out : public trj {
private:
    /**< vector of navigation states */
    std::vector<st::st_nav_out,Eigen::aligned_allocator<st::st_nav_out>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_nav_out() = delete;
    /**< constructor based on time separation between consecutive navigation samples, size of trajectory,
    time separation between consecutive truth samples, and number of operations */
    trj_nav_out(const double& Deltat_sec, const unsigned int& nel, const double& Deltat_sec_truth, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, (int)(Deltat_sec / Deltat_sec_truth)), _Vst(nel)  {
        if (std::remainder(_Deltat_sec, Deltat_sec_truth) > math::constant::EPS()) {throw std::runtime_error("Time separation between samples does not match.");}
    }
    /**< copy constructor */
    trj_nav_out(const trj_nav_out&) = delete;
    /**< move constructor */
    trj_nav_out(trj_nav_out&&) = delete;
    /**< destructor */
    ~trj_nav_out() = default;
    /**< copy assignment */
    trj_nav_out& operator=(const trj_nav_out&) = delete;
    /**< move assignment */
    trj_nav_out& operator=(trj_nav_out&&) = delete;

    /**< get vector of navigation states */
    std::vector<st::st_nav_out,Eigen::aligned_allocator<st::st_nav_out>>& operator()() {return _Vst;}
    const std::vector<st::st_nav_out,Eigen::aligned_allocator<st::st_nav_out>>& operator()() const {return _Vst;}
    /**< get vector of navigation states */
    std::vector<st::st_nav_out,Eigen::aligned_allocator<st::st_nav_out>>& get() {return _Vst;}
    const std::vector<st::st_nav_out,Eigen::aligned_allocator<st::st_nav_out>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_nav_out

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif















