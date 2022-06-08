#ifndef ACFT_INIT_ERROR
#define ACFT_INIT_ERROR

#include "../acft.h"
#include "../logic.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace math {
    class seeder;
}
namespace ang {
    class rodrigues;
}
namespace env {
    class earth;
}
namespace st {
    class sti;
    class st_sens_out;
    class error_gen_att;
    class error_gen_triple;
    class error_gen_single;

// CLASS INIT_ERROR
// ================
// ================

class ACFT_API init_error {
private:
    /**< pointer to initial Euler angles error generator */
    st::error_gen_att* _Piniterr_eul;
    /**< pointer to initial accelerometer error generator */
    st::error_gen_triple* _Piniterr_acc;
    /**< pointer to initial gyroscope error generator */
    st::error_gen_triple* _Piniterr_gyr;
    /**< pointer to initial magnetometer error generator */
    st::error_gen_triple* _Piniterr_mag;
    /**< pointer to initial NED magnetic field error generator */
    st::error_gen_triple* _Piniterr_mgn;
    /**< pointer to initial height over terrain error generator */
    st::error_gen_single* _Piniterr_hgh;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**< default constructor */
    init_error() = delete;
    /**< constructor based on seed order and identifiers */
    init_error(math::seeder& Oseeder,
               st::logic::INITEUL_ID initeul_id, st::logic::INITACC_ID initacc_id, st::logic::INITGYR_ID initgyr_id,
               st::logic::INITMAG_ID initmag_id, st::logic::INITMGN_ID initmgn_id, st::logic::INITHGH_ID inithgh_id);
    /**< copy constructor */
    init_error(const init_error&) = delete;
    /**< move constructor */
    init_error(init_error&&) = delete;
    /**< destructor */
    ~init_error();
    /**< copy assignment */
    init_error& operator=(const init_error&) = delete;
    /**< move assignment */
    init_error& operator=(init_error&&) = delete;

    /**< fills up the initial estimated aircraft attitude based on the initial conditions,
     * and writes results in input stream */
    void eval_eul(const st::sti& Osti, ang::rodrigues& q_nb_init, std::ostream& Ostream) const;
    /**< modifies the input rodrigues parameters with the initialization error */
    void eval_eul(ang::rodrigues& q_nb_init) const;

    /**< fills up the initial estimated full accelerometer error based on the initial conditions,
     * and writes the results in input stream */
    void eval_acc(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Eacc_init_std_mps2, Eigen::Vector3d& Eacc_init_mps2, std::ostream& Ostream) const;
    /**< fills up the initial estimated full gyroscope error based on the initial conditions,
     * and writes the results in input stream */
    void eval_gyr(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Egyr_init_std_rps, Eigen::Vector3d& Egyr_init_rps, std::ostream& Ostream) const;
    /**< fills up the initial estimated full magnetometer error based on the initial conditions,
     * and writes the results in input stream */
    void eval_mag(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Emag_init_std_nT, Eigen::Vector3d& Emag_init_nT, std::ostream& Ostream) const;
    /**< fills up the initial estimated difference of the magnetic field between the models are reality
     * based on the initial conditions, and writes results in input stream */
    void eval_mgn(const st::sti& Osti, const env::earth& Oearth, Eigen::Vector3d& Berror_init_std_nT, Eigen::Vector3d& Berror_init_nT, std::ostream& Ostream) const;

    /**< returns the initial error in the estimated height over the terrain */
    double eval_hgh() const;

}; // closes class init_error

}; // closes namespace st

#endif


