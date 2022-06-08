#ifndef ACFT_TRJ_NAV_IN
#define ACFT_TRJ_NAV_IN

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

// CLASS NAVIGATION INPUT STATE ST_NAV_IN
// ======================================
// ======================================

class ACFT_API st_nav_in {
private:
    /**< time */
    double _t_sec;
    /**< geodetic coordinates */
    env::geodetic_coord _x_gdt_b_rad_m;
    /**< cartesian coordinates */
    env::cartesian_coord _x_car_b_m;
    /**< NED to BFS quaternion */
    ang::rodrigues _q_nb;
    /**< ECEF to BFS quaternion */
    ang::rodrigues _q_eb;
    /**< unit dual quaternion from ECEF to BFS */
    ang::dual _z_eb;
    /**< angular velocity of BFS with respect to NED viewed in BFS. */
    ang::so3_tangent _w_nbb_rps;
    /**< angular velocity of BFS with respect to ECEF viewed in BFS. */
    ang::so3_tangent _w_ebb_rps;
    /**< angular velocity of BFS with respect to ECEF viewed in ECEF. */
    ang::so3_tangent _w_ebe_rps;
    /**< NED absolute velocity. */
    Eigen::Vector3d _v_n_mps;
    /**< ECEF absolute velocity. */
    Eigen::Vector3d _v_e_mps;
    /**< twist from ECEF to BFS */
    ang::se3_tangent _xi_ebb_mrps;
    /**< true airspeed */
    double _vtas_mps;
    /**< Euler angles from WFS to BFS */
    ang::euler _euler_wb;
    /**< atmospheric temperature */
    double _T_degK;
    /**< pressure altitude */
    double _Hp_m;
    /**< specific force */
    Eigen::Vector3d _f_ibb_mps2;
    /**< radius of curvature of prime vertical N */
    double _N_m;
    /**< radius of curvature of meridian M */
    double _M_m;
    /**< magnetic field error (model minus real values) */
    Eigen::Vector3d _B_n_nT_dev;
    /**< NED low frequency wind field */
    Eigen::Vector3d _vlf_n_mps;
    /**< temperature offset at mean sea level */
    double _DeltaT_degK;
    /**< pressure offset at mean sea level */
    double _Deltap_pa;
    /**< horizontal air distance */
    double _dist_air_hor_m;
    /**< horizontal ground distance */
    double _dist_grd_hor_m;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_nav_in() = default;
    /**< copy constructor */
    st_nav_in(const st_nav_in&) = default;
    /**< move constructor */
    st_nav_in(st_nav_in&&) = default;
    /**< destructor */
    ~st_nav_in() = default;
    /**< copy assignment */
    st_nav_in& operator=(const st_nav_in&) = default;
    /**< move assignment */
    st_nav_in& operator=(st_nav_in&&) = default;

    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}
    /**< get geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m() {return _x_gdt_b_rad_m;}
    const env::geodetic_coord& get_x_gdt_b_rad_m() const {return _x_gdt_b_rad_m;}
    /**< get cartesian coordinates */
    env::cartesian_coord& get_x_car_b_m() {return _x_car_b_m;}
    const env::cartesian_coord& get_x_car_b_m() const {return _x_car_b_m;}
    /**< get NED to BFS quaternion */
    ang::rodrigues& get_q_nb() {return _q_nb;}
    const ang::rodrigues& get_q_nb() const {return _q_nb;}
    /**< get ECEF to BFS quaternion */
    ang::rodrigues& get_q_eb() {return _q_eb;}
    const ang::rodrigues& get_q_eb() const {return _q_eb;}
    /**< get unit dual quaternion from ECEF to BFS */
    ang::dual& get_z_eb() {return _z_eb;}
    const ang::dual& get_z_eb() const {return _z_eb;}
    /**< get angular velocity of BFS with respect to NED viewed in BFS. */
    ang::so3_tangent& get_w_nbb_rps() {return _w_nbb_rps;}
    const ang::so3_tangent& get_w_nbb_rps() const {return _w_nbb_rps;}
    /**< get angular velocity of BFS with respect to ECEF viewed in BFS. */
    ang::so3_tangent& get_w_ebb_rps() {return _w_ebb_rps;}
    const ang::so3_tangent& get_w_ebb_rps() const {return _w_ebb_rps;}
    /**< get angular velocity of BFS with respect to ECEF viewed in ECEF. */
    ang::so3_tangent& get_w_ebe_rps() {return _w_ebe_rps;}
    const ang::so3_tangent& get_w_ebe_rps() const {return _w_ebe_rps;}
    /**< get NED absolute velocity. */
    Eigen::Vector3d& get_v_n_mps() {return _v_n_mps;}
    const Eigen::Vector3d& get_v_n_mps() const {return _v_n_mps;}
    /**< get ECEF absolute velocity. */
    Eigen::Vector3d& get_v_e_mps() {return _v_e_mps;}
    const Eigen::Vector3d& get_v_e_mps() const {return _v_e_mps;}
    /**< get twist from ECEF to BFS */
    ang::se3_tangent& get_xi_ebb_mrps() {return _xi_ebb_mrps;}
    const ang::se3_tangent& get_xi_ebb_mrps() const {return _xi_ebb_mrps;}
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
    /**< get specific force */
    Eigen::Vector3d& get_f_ibb_mps2() {return _f_ibb_mps2;}
    const Eigen::Vector3d& get_f_ibb_mps2() const {return _f_ibb_mps2;}
    /**< get radius of curvature of prime vertical N */
    double& get_N_m() {return _N_m;}
    const double& get_N_m() const {return _N_m;}
    /**< get radius of curvature of meridian M */
    double& get_M_m() {return _M_m;}
    const double& get_M_m() const {return _M_m;}
    /**< get magnetic field error (model minus real values) */
    Eigen::Vector3d& get_B_n_nT_dev() {return _B_n_nT_dev;}
    const Eigen::Vector3d& get_B_n_nT_dev() const {return _B_n_nT_dev;}
    /**< get NED low frequency wind field */
    Eigen::Vector3d& get_vlf_n_mps() {return _vlf_n_mps;}
    const Eigen::Vector3d& get_vlf_n_mps() const {return _vlf_n_mps;}
    /**< get temperature offset at mean sea level */
    double& get_DeltaT_degK() {return _DeltaT_degK;}
    const double& get_DeltaT_degK() const {return _DeltaT_degK;}
    /**< get pressure offset at mean sea level */
    double& get_Deltap_pa() {return _Deltap_pa;}
    const double& get_Deltap_pa() const {return _Deltap_pa;}
    /**< get horizontal air distance */
    double& get_dist_air_hor_m() {return _dist_air_hor_m;}
    const double& get_dist_air_hor_m() const {return _dist_air_hor_m;}
    /**< get horizontal ground distance */
    double& get_dist_grd_hor_m() {return _dist_grd_hor_m;}
    const double& get_dist_grd_hor_m() const {return _dist_grd_hor_m;}
}; // closes class st_nav_in

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS NAVIGATION INPUT TRAJECTORY TRJ_NAV_IN
// ============================================
// ============================================

class ACFT_API trj_nav_in : public trj {
private:
    /**< vector of navigation states */
    std::vector<st::st_nav_in,Eigen::aligned_allocator<st::st_nav_in>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_nav_in() = delete;
    /**< constructor based on time separation between consecutive navigation samples, size of trajectory,
    time separation between consecutive truth samples, and number of operations */
    trj_nav_in(const double& Deltat_sec, const unsigned int& nel, const double& Deltat_sec_truth, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, (int)(Deltat_sec / Deltat_sec_truth)), _Vst(nel) {
        if (std::remainder(_Deltat_sec, Deltat_sec_truth) > math::constant::EPS()) {throw "Time separation between samples does not match.";}
    }
    /**< copy constructor */
    trj_nav_in(const trj_nav_in&) = delete;
    /**< move constructor */
    trj_nav_in(trj_nav_in&&) = delete;
    /**< destructor */
    ~trj_nav_in() = default;
    /**< copy assignment */
    trj_nav_in& operator=(const trj_nav_in&) = delete;
    /**< move assignment */
    trj_nav_in& operator=(trj_nav_in&&) = delete;

    /**< get vector of navigation states */
    std::vector<st::st_nav_in,Eigen::aligned_allocator<st::st_nav_in>>& operator()() {return _Vst;}
    const std::vector<st::st_nav_in,Eigen::aligned_allocator<st::st_nav_in>>& operator()() const {return _Vst;}
    /**< get vector of navigation states */
    std::vector<st::st_nav_in,Eigen::aligned_allocator<st::st_nav_in>>& get() {return _Vst;}
    const std::vector<st::st_nav_in,Eigen::aligned_allocator<st::st_nav_in>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_nav_in

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif















