#include "init_error.h"
#include "ang/rotate/euler.h"
#include "env/earth.h"
#include "env/geo.h"
#include "env/mag.h"
#include "env/coord.h"
#include "error_gen_triple.h"
#include "error_gen_single.h"
#include "error_gen_att.h"
#include "math/logic/seeder.h"
#include "acft/st/sti.h"
#include "acft/st/trj_sens_out.h"

// CLASS INIT_ERROR
// ================
// ================

st::init_error::init_error(math::seeder& Oseeder,
                            st::logic::INITEUL_ID initeul_id, st::logic::INITACC_ID initacc_id, st::logic::INITGYR_ID initgyr_id,
                            st::logic::INITMAG_ID initmag_id, st::logic::INITMGN_ID initmgn_id, st::logic::INITHGH_ID inithgh_id)
: _Piniterr_eul(nullptr), _Piniterr_acc(nullptr), _Piniterr_gyr(nullptr), _Piniterr_mag(nullptr), _Piniterr_mgn(nullptr), _Piniterr_hgh(nullptr) {
    _Piniterr_eul = st::error_gen_att::create_init_eul_error_generator(initeul_id, Oseeder.provide_seed(math::seeder::seeder_initeul));
    _Piniterr_acc = st::error_gen_triple::create_init_acc_error_generator(initacc_id, Oseeder.provide_seed(math::seeder::seeder_initacc));
    _Piniterr_gyr = st::error_gen_triple::create_init_gyr_error_generator(initgyr_id, Oseeder.provide_seed(math::seeder::seeder_initgyr));
    _Piniterr_mag = st::error_gen_triple::create_init_mag_error_generator(initmag_id, Oseeder.provide_seed(math::seeder::seeder_initmag));
    _Piniterr_mgn = st::error_gen_triple::create_init_mgn_error_generator(initmgn_id, Oseeder.provide_seed(math::seeder::seeder_initmgn));
    _Piniterr_hgh = st::error_gen_single::create_init_height_error_generator(inithgh_id, Oseeder.provide_seed(math::seeder::seeder_initheight));
}
/* constructor based on seed order, specific accelerometer, gyroscope, magnetometer, pressure, temperature, airspeed, angle of attack, angle of sideslip */

st::init_error::~init_error() {
    delete _Piniterr_eul;
    delete _Piniterr_acc;
    delete _Piniterr_gyr;
    delete _Piniterr_mag;
    delete _Piniterr_mgn;
    delete _Piniterr_hgh;
}
/* destructor */

void st::init_error::eval_eul(const st::sti& Osti, ang::rodrigues& q_nb_init, std::ostream& Ostream) const {
    double q_std_sigma_deg = _Piniterr_eul->get_sigma_deg();
    ang::rodrigues q_add(_Piniterr_eul->eval());
    ang::euler euler_nb_zero(Osti.get_psi_rad(), Osti.get_theta_rad(), Osti.get_xi_rad());
    ang::rodrigues q_nb_zero(euler_nb_zero);
    q_nb_init = q_add * q_nb_zero;
    ang::euler euler_nb_init(q_nb_init);

    Ostream << std::endl;
    Ostream << "Euler truth [deg]:      " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Osti.get_psi_rad() * math::constant::R2D()
            << std::setw(9) << Osti.get_theta_rad() * math::constant::R2D()
            << std::setw(9) << Osti.get_xi_rad() * math::constant::R2D()
            << std::endl;
    Ostream << "Euler filter [deg]:     " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << euler_nb_init.get_yaw_rad() * math::constant::R2D()
            << std::setw(9) << euler_nb_init.get_pitch_rad() * math::constant::R2D()
            << std::setw(9) << euler_nb_init.get_bank_rad() * math::constant::R2D()
            << std::endl;
    Ostream << "Euler filter std [deg]: " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << q_std_sigma_deg
            << std::endl;
}
/* fills up the initial estimated aircraft attitude based on the initial conditions
 * and writes results in input stream */

void st::init_error::eval_eul(ang::rodrigues& q_nb_init) const {
    q_nb_init = _Piniterr_eul->eval() * q_nb_init;
}
/* modifies the input rodrigues parameters with the initialization error */

void st::init_error::eval_acc(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Eacc_init_std_mps2, Eigen::Vector3d& Eacc_init_mps2, std::ostream& Ostream) const {
    Eigen::Array3d init_acc_std        = Eigen::Array3d::Ones() * _Piniterr_acc->get_sigma();
    Eigen::Array3d init_acc_multiply   = Eigen::Array3d::Ones() + _Piniterr_acc->eval().array();
    Eacc_init_std_mps2 = (Ost_sens_out_ini.get_E_acc_mps2().array() * init_acc_std).matrix();
    Eacc_init_mps2     = (Ost_sens_out_ini.get_E_acc_mps2().array() * init_acc_multiply).matrix();

    // TODO: The initial error Eacc_init_rps is properly computed, but not its standard deviation Eacc_init_std_rps, as it makes use
    // of the realization Ost_sens_out_init.get_E_acc_mps2(), while it should employ fixed terms such as the bias offset and the random walk
    // I found this on Nov 2021 but do not correct it as all results are based on existing method
    // This is a serious error with implications for the position filter initialization

    Ostream << std::endl;
    Ostream << "Eacc init truth [mps2]: " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Ost_sens_out_ini.get_E_acc_mps2()(0)
            << std::setw(9) << Ost_sens_out_ini.get_E_acc_mps2()(1)
            << std::setw(9) << Ost_sens_out_ini.get_E_acc_mps2()(2)
            << std::endl;
    Ostream << "Eacc init filter [mps2]:" << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Eacc_init_mps2(0)
            << std::setw(9) << Eacc_init_mps2(1)
            << std::setw(9) << Eacc_init_mps2(2)
            << std::endl;
    Ostream << "Eacc init std [mps2]:   " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Eacc_init_std_mps2(0)
            << std::setw(9) << Eacc_init_std_mps2(1)
            << std::setw(9) << Eacc_init_std_mps2(2)
            << std::endl;
}
/* fills up the initial estimated full accelerometer error based on the initial conditions,
 * and writes the results in input stream */

void st::init_error::eval_gyr(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Egyr_init_std_rps, Eigen::Vector3d& Egyr_init_rps, std::ostream& Ostream) const {
    Eigen::Array3d init_gyr_std       = Eigen::Array3d::Ones() * _Piniterr_gyr->get_sigma();
    Eigen::Array3d init_gyr_multiply  = Eigen::Array3d::Ones() + _Piniterr_gyr->eval().array();
    Egyr_init_std_rps = (Ost_sens_out_ini.get_E_gyr_rps().array() * init_gyr_std).matrix();
    Egyr_init_rps     = (Ost_sens_out_ini.get_E_gyr_rps().array() * init_gyr_multiply).matrix();

    // TODO: The initial error Egyr_init_rps is properly computed, but not its standard deviation Egyr_init_std_rps, as it makes use
    // of the realization Ost_sens_out_init.get_E_gyr_rps(), while it should employ fixed terms such as the bias offset and the random walk
    // I found this on Nov 2021 but do not correct it as all results are based on existing method
    // This is a serious error with implications for the attitude filter initialization

    //std::cout << "AAAAAAAAAAAAA" << std::fixed << std::setprecision(8) << std::showpos << init_gyr_std << std::endl;
    //std::cout << "BBBBBBBBBBBBB" << std::fixed << std::setprecision(8) << std::showpos << init_gyr_multiply << std::endl;
    //std::cout << "CCCCCCCCCCCCC" << std::fixed << std::setprecision(8) << std::showpos << Egyr_init_std_rps << std::endl;
    //std::cout << "DDDDDDDDDDDDD" << std::fixed << std::setprecision(8) << std::showpos << Egyr_init_rps << std::endl;
    //std::cout << "EEEEEEEEEEEEE" << std::fixed << std::setprecision(8) << std::showpos << Ost_sens_out_ini.get_E_gyr_rps() << std::endl;

    Ostream << std::endl;
    Ostream << "Egyr init truth [dps]:  " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Ost_sens_out_ini.get_E_gyr_rps()(0) * math::constant::R2D()
            << std::setw(9) << Ost_sens_out_ini.get_E_gyr_rps()(1) * math::constant::R2D()
            << std::setw(9) << Ost_sens_out_ini.get_E_gyr_rps()(2) * math::constant::R2D()
            << std::endl;
    Ostream << "Egyr init filter [dps]: " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Egyr_init_rps(0) * math::constant::R2D()
            << std::setw(9) << Egyr_init_rps(1) * math::constant::R2D()
            << std::setw(9) << Egyr_init_rps(2) * math::constant::R2D()
            << std::endl;
    Ostream << "Egyr init std [dps]:    " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Egyr_init_std_rps(0) * math::constant::R2D()
            << std::setw(9) << Egyr_init_std_rps(1) * math::constant::R2D()
            << std::setw(9) << Egyr_init_std_rps(2) * math::constant::R2D()
            << std::endl;
}
/* fills up the initial estimated full gyroscope error based on the initial conditions,
 * and writes the results in input stream */

void st::init_error::eval_mag(const st::st_sens_out& Ost_sens_out_ini, Eigen::Vector3d& Emag_init_std_nT, Eigen::Vector3d& Emag_init_nT, std::ostream& Ostream) const {
    Eigen::Array3d init_mag_std      = Eigen::Array3d::Ones() * _Piniterr_mag->get_sigma();
    Eigen::Array3d init_mag_multiply = Eigen::Array3d::Ones() + _Piniterr_mag->eval().array();
    Emag_init_std_nT = (Ost_sens_out_ini.get_E_mag_nT().array() * init_mag_std).matrix();
    Emag_init_nT     = (Ost_sens_out_ini.get_E_mag_nT().array() * init_mag_multiply).matrix();

    // TODO: The initial error Emag_init_rps is properly computed, but not its standard deviation Emag_init_std_rps, as it makes use
    // of the realization Ost_sens_out_init.get_E_mag_nT(), while it should employ fixed terms such as the bias offset and the hard iron magnetism
    // I found this on Nov 2021 but do not correct it as all results are based on existing method
    // This is a serious error with implications for the attitude filter initialization

    Ostream << std::endl;
    Ostream << "Emag init truth [nT]:   " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Ost_sens_out_ini.get_E_mag_nT()(0)
            << std::setw(9) << Ost_sens_out_ini.get_E_mag_nT()(1)
            << std::setw(9) << Ost_sens_out_ini.get_E_mag_nT()(2)
            << std::endl;
    Ostream << "Emag init filter [nT]:  " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Emag_init_nT(0)
            << std::setw(9) << Emag_init_nT(1)
            << std::setw(9) << Emag_init_nT(2)
            << std::endl;
    Ostream << "Emag init std [nT]:     " << std::fixed << std::setprecision(3) << std::showpos
            << std::setw(9) << Emag_init_std_nT(0)
            << std::setw(9) << Emag_init_std_nT(1)
            << std::setw(9) << Emag_init_std_nT(2)
            << std::endl;
}
/* fills up the initial estimated full magnetometer error based on the initial conditions,
 * and writes the results in input stream */

void st::init_error::eval_mgn(const st::sti& Osti, const env::earth& Oearth, Eigen::Vector3d& Berror_init_std_nT, Eigen::Vector3d& Berror_init_nT, std::ostream& Ostream) const {
    env::geodetic_coord x_gdt_rad_m_init(Osti.get_lambda_rad(), Osti.get_phi_rad(), Osti.get_h_m());
    Eigen::Vector3d B_n_nT_init_truth = Oearth.get_geo().get_mag().compute_B_n_nT_truth(Osti.get_t_sec(), x_gdt_rad_m_init);
    Eigen::Vector3d B_n_nT_init_model = Oearth.get_geo().get_mag().compute_B_n_nT_model(Osti.get_t_sec(), x_gdt_rad_m_init);
    Eigen::Vector3d B_n_nT_init_diff  = B_n_nT_init_model - B_n_nT_init_truth;
    Eigen::Array3d init_Berror_std      = Eigen::Array3d::Ones() * _Piniterr_mgn->get_sigma();
    Eigen::Array3d init_Berror_multiply = Eigen::Array3d::Ones() + _Piniterr_mgn->eval().array();
    Berror_init_std_nT = (B_n_nT_init_diff.array() * init_Berror_std).matrix();
    Berror_init_nT     = (B_n_nT_init_diff.array() * init_Berror_multiply).matrix();

    // TODO: The initial error Berror_init_nT is properly computed, but not its standard deviation Berror_init_std_nT, as it makes use
    // of the realization B_n_nT_init_diff), while it should employ fixed terms such as the mag class sigma_B_nT_north = 138.0 nT,
    // sigma_B_nT_east = 89.0 nT, and sigma_B_nT_down = 165.0 nT. I found this on Nov 2021 but do not correct it as all results are
    // based on existing method. This is a serious error with implications for the attitude filter initialization.

    Ostream << std::endl;
    Ostream << "B_n init truth [nT]:    " << std::fixed << std::setprecision(1) << std::showpos
            << std::setw(9) << B_n_nT_init_truth(0)
            << std::setw(9) << B_n_nT_init_truth(1)
            << std::setw(9) << B_n_nT_init_truth(2)
            << std::endl;
    Ostream << "B_n init model [nT]:    " << std::fixed << std::setprecision(1) << std::showpos
            << std::setw(9) << B_n_nT_init_model(0)
            << std::setw(9) << B_n_nT_init_model(1)
            << std::setw(9) << B_n_nT_init_model(2)
            << std::endl;
    Ostream << "B_n init diff truth [nT]:" << std::fixed << std::setprecision(2) << std::showpos
            << std::setw(8) << B_n_nT_init_diff(0)
            << std::setw(9) << B_n_nT_init_diff(1)
            << std::setw(9) << B_n_nT_init_diff(2)
            << std::endl;
    Ostream << "B_n init diff filter [nT]:" << std::fixed << std::setprecision(2) << std::showpos
            << std::setw(7) << Berror_init_nT(0)
            << std::setw(9) << Berror_init_nT(1)
            << std::setw(9) << Berror_init_nT(2)
            << std::endl;
    Ostream << "B_n init diff std [nT]: " << std::fixed << std::setprecision(2) << std::showpos
            << std::setw(9) << Berror_init_std_nT(0)
            << std::setw(9) << Berror_init_std_nT(1)
            << std::setw(9) << Berror_init_std_nT(2)
            << std::endl;
}
/* fills up the initial estimated difference of the magnetic field between the models are reality
 * based on the initial conditions, and writes results in input stream */

double st::init_error::eval_hgh() const {
    return _Piniterr_hgh->eval();
}
/* returns the initial error in the estimated height over the terrain, and writes results in input stream */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////







