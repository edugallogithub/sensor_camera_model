#ifndef ACFT_TRJ_IMG_TRUTH
#define ACFT_TRJ_IMG_TRUTH

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

// CLASS IMAGE INPUT STATE ST_IMG_TRUTH
// ====================================
// ====================================

class ACFT_API st_img_truth {
private:
    /**< image id (starts in 0 and grows) */
    int _id;
    /**< time */
    double _t_sec;
    /**< center of gravity (not camera center) geodetic coordinates (to generate images). */
    env::geodetic_coord _x_gdt_b_rad_m;
    /**< Euler angles from NED to BFS (to generate images) */
    ang::euler _euler_nb;
    /**< mass ratio */
    double _ratio;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    st_img_truth() = default;
    /**< copy constructor */
    st_img_truth(const st_img_truth&) = default;
    /**< move constructor */
    st_img_truth(st_img_truth&&) = default;
    /**< destructor */
    ~st_img_truth() = default;
    /**< copy assignment */
    st_img_truth& operator=(const st_img_truth&) = default;
    /**< move assignment */
    st_img_truth& operator=(st_img_truth&&) = default;

    /**< get image id (starts in 0 and grows) */
    int& get_id() {return _id;}
    const int& get_id() const {return _id;}
    /**< get time */
    double& get_t_sec() {return _t_sec;}
    const double& get_t_sec() const {return _t_sec;}
    /**< get center of gravity (not camera center) geodetic coordinates */
    env::geodetic_coord& get_x_gdt_b_rad_m() {return _x_gdt_b_rad_m;}
    const env::geodetic_coord& get_x_gdt_b_rad_m() const {return _x_gdt_b_rad_m;}
    /**< get Euler angles from NED to BFS */
    ang::euler& get_euler_nb() {return _euler_nb;}
    const ang::euler& get_euler_nb() const {return _euler_nb;}
    /**< get mass ratio */
    double& get_ratio() {return _ratio;}
    const double& get_ratio() const {return _ratio;}
}; // closes class st_img_truth

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS IMAGE INPUT TRAJECTORY TRJ_IMG_TRUTH
// ==========================================
// ==========================================

class ACFT_API trj_img_truth : public trj{
private:
    /**< vector of sensor input states */
    std::vector<st::st_img_truth,Eigen::aligned_allocator<st::st_img_truth>> _Vst;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    trj_img_truth() = delete;
    /**< constructor based on time separation between consecutive images, size of trajectory,
     * time separation between consecutive truth samples, and number of operations */
    trj_img_truth(const double& Deltat_sec, const unsigned int& nel, const double& Deltat_sec_truth, const unsigned short& nel_op)
        : trj(Deltat_sec, nel, nel_op, (int)(Deltat_sec / Deltat_sec_truth)), _Vst(nel) {
        if (std::remainder(_Deltat_sec, Deltat_sec_truth) > math::constant::EPS()) {throw std::runtime_error("Time separation between samples does not match.");}
    }
    /**< copy constructor */
    trj_img_truth(const trj_img_truth&) = delete;
    /**< move constructor */
    trj_img_truth(trj_img_truth&&) = delete;
    /**< destructor */
    ~trj_img_truth() override = default;
    /**< copy assignment */
    trj_img_truth& operator=(const trj_img_truth&) = delete;
    /**< move assignment */
    trj_img_truth& operator=(trj_img_truth&&) = delete;

    /**< get vector of input image states */
    std::vector<st::st_img_truth,Eigen::aligned_allocator<st::st_img_truth>>& operator()() {return _Vst;}
    const std::vector<st::st_img_truth,Eigen::aligned_allocator<st::st_img_truth>>& operator()() const {return _Vst;}
    /**< get vector of input image states */
    std::vector<st::st_img_truth,Eigen::aligned_allocator<st::st_img_truth>>& get() {return _Vst;}
    const std::vector<st::st_img_truth,Eigen::aligned_allocator<st::st_img_truth>>& get() const {return _Vst;}

    /**< resize trajectory to new size, leaving only the first members. Number of operations does not change. */
    void resize_st(const unsigned int& nel) override {
        _Vst.resize(nel);
        _nel = nel;
    }
}; // closes class trj_img_truth

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace st

#endif
