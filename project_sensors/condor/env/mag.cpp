#include "mag.h"
#include "ang/rotate/euler.h"
#include "env/coord.h"
#include <iostream>
#include <iomanip>

// CLASS MAG
// =========
// =========

ang::euler env::mag::compute_euler_mag(const Eigen::Vector3d& mag_n_nT) {
    return ang::euler(ang::euler::obtain_yaw_forward(mag_n_nT), ang::euler::obtain_pitch_forward(mag_n_nT), 0.0);
}
/* returns the magnetic Euler angles (declination, minus inclination, 0) based on the magnetic field */

env::mag* env::mag::create_mag(env::logic::ZONE_ID zone_id, env::logic::REALISM_ID realism_id, std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen) {
    double sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down;
    double d2r = math::constant::D2R();
    switch (realism_id) {
        case env::logic::realism_grav_yes_magn_yes:
        case env::logic::realism_grav_no_magn_yes:
            sigma_B_nT_north = 138.0;
            sigma_B_nT_east  = 89.0;
            sigma_B_nT_down  = 165.0;
            break;
        case env::logic::realism_grav_yes_magn_no:
        case env::logic::realism_grav_no_magn_no:
            sigma_B_nT_north = 0.;
            sigma_B_nT_east  = 0.;
            sigma_B_nT_down  = 0.;
            break;
        case env::logic::realism_size:
        default:
            throw std::runtime_error("Realism option not available");
    }

    env::mag* Pmag = nullptr;
    switch (zone_id) {
        case env::logic::zone_default:
            Pmag = new env::mag_constant(Eigen::Vector3d(17478.0, -1531.5, 52134.2), dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break;
        case env::logic::zone_desert: {
            env::geodetic_coord x_gdt_rad_m_ul(+248.0 * d2r, +33.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+249.0 * d2r, +33.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+248.0 * d2r, +32.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+249.0 * d2r, +32.0 * d2r, +3000);
            Eigen::Vector3d B_n_nT_ul(23708.8, 4205.8, 40325.4);
            Eigen::Vector3d B_n_nT_ur(23672.4, 4064.6, 40535.5);
            Eigen::Vector3d B_n_nT_ll(24050.6, 4199.7, 39432.9);
            Eigen::Vector3d B_n_nT_lr(24017.7, 4060.8, 39642.7);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_urban: {
            env::geodetic_coord x_gdt_rad_m_ul(+241.0 * d2r, +34.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+242.0 * d2r, +34.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+241.0 * d2r, +33.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+242.0 * d2r, +33.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(23618.6, 5015.2, 39617.6);
            Eigen::Vector3d B_n_nT_ur(23583.6, 4920.5, 39854.7);
            Eigen::Vector3d B_n_nT_ll(23940.3, 4996.0, 38740.2);
            Eigen::Vector3d B_n_nT_lr(23910.0, 4903.5, 38976.2);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_everglades: {
            env::geodetic_coord x_gdt_rad_m_ul(+279.0 * d2r, +26.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+280.0 * d2r, +26.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+279.0 * d2r, +25.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+280.0 * d2r, +25.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(25077.3, -2795.3, 35743.1);
            Eigen::Vector3d B_n_nT_ur(25069.0, -3065.0, 35586.0);
            Eigen::Vector3d B_n_nT_ll(25314.9, -2759.2, 34794.9);
            Eigen::Vector3d B_n_nT_lr(25302.5, -3030.2, 34639.6);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_forest: {
            env::geodetic_coord x_gdt_rad_m_ul(+287.0 * d2r, +44.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+288.0 * d2r, +44.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+287.0 * d2r, +43.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+288.0 * d2r, +43.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(18731.8, -4671.3, 49143.0);
            Eigen::Vector3d B_n_nT_ur(18812.9, -4846.3, 48922.8);
            Eigen::Vector3d B_n_nT_ll(19218.1, -4701.3, 48456.7);
            Eigen::Vector3d B_n_nT_lr(19297.3, -4878.9, 48232.0);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_farm: {
            env::geodetic_coord x_gdt_rad_m_ul(+272.0 * d2r, +39.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+273.0 * d2r, +39.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+272.0 * d2r, +38.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+273.0 * d2r, +38.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(20504.5, -1182.5, 47737.6);
            Eigen::Vector3d B_n_nT_ur(20506.9, -1449.3, 47679.8);
            Eigen::Vector3d B_n_nT_ll(20950.1, -1165.7, 46967.4);
            Eigen::Vector3d B_n_nT_lr(20951.2, -1434.3, 46909.8);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_mix: {
            env::geodetic_coord x_gdt_rad_m_ul(+271.0 * d2r, +35.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+272.0 * d2r, +35.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+271.0 * d2r, +34.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+272.0 * d2r, +34.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(22211.8,  -837.9, 44581.8);
            Eigen::Vector3d B_n_nT_ur(22204.8, -1110.6, 44539.6);
            Eigen::Vector3d B_n_nT_ll(22603.9,  -816.6, 43736.2);
            Eigen::Vector3d B_n_nT_lr(22595.3, -1090.1, 43694.4);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_mexico: {
            env::geodetic_coord x_gdt_rad_m_ul(+259.5 * d2r, +20.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+260.5 * d2r, +20.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+259.5 * d2r, +19.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+260.5 * d2r, +19.5 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(27056.8, 	2357.6,	30307.8);
            Eigen::Vector3d B_n_nT_ur(27017.5, 	2161.0,	30432.4);
            Eigen::Vector3d B_n_nT_ll(27277.5, 	2384.1,	29273.4);
            Eigen::Vector3d B_n_nT_lr(27236.4, 	2190.9, 29399.1);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_colombia: {
            env::geodetic_coord x_gdt_rad_m_ul(+289.5 * d2r, +5.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+290.5 * d2r, +5.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+289.5 * d2r, +4.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+290.5 * d2r, +4.5 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(26615.5, -4728.6, 14064.6);
            Eigen::Vector3d B_n_nT_ur(26571.5, -5023.6,	13799.1);
            Eigen::Vector3d B_n_nT_ll(26558.1, -4673.3,	13133.1);
            Eigen::Vector3d B_n_nT_lr(26509.5, -4971.9,	12869.5);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_argentina: {
            env::geodetic_coord x_gdt_rad_m_ul(+299.5 * d2r, -24.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+300.5 * d2r, -24.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+299.5 * d2r, -25.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+300.5 * d2r, -25.5 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(19120.4, -4205.9, -10532.6);
            Eigen::Vector3d B_n_nT_ur(18926.5, -4447.3,	-10769.9);
            Eigen::Vector3d B_n_nT_ll(18905.6, -4026.8,	-10953.1);
            Eigen::Vector3d B_n_nT_lr(18711.0, -4268.6,	-11178.7);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_canada: {
            env::geodetic_coord x_gdt_rad_m_ul(+254.5 * d2r, +52.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+255.5 * d2r, +52.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+254.5 * d2r, +51.5 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+255.5 * d2r, +51.5 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(14038.7, 2354.1, 55209.6);
            Eigen::Vector3d B_n_nT_ur(13961.6, 2164.1, 55301.2);
            Eigen::Vector3d B_n_nT_ll(14619.9, 2428.6, 54753.5);
            Eigen::Vector3d B_n_nT_lr(14543.3, 2236.5, 54851.0);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }

        case env::logic::zone_rozas:
            Pmag = new env::mag_constant(Eigen::Vector3d(24175.8, -749.4, 38813.1), dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break;
        case env::logic::zone_wisconsin: {
            //env::geodetic_coord x_gdt_rad_m_ul(-90.0, +45.0, +3000);
            //env::geodetic_coord x_gdt_rad_m_ur(-89.0, +45.0, +3000);
            //env::geodetic_coord x_gdt_rad_m_ll(-90.0, +44.0, +3000);
            //env::geodetic_coord x_gdt_rad_m_lr(-89.0, +44.0, +3000);
            env::geodetic_coord x_gdt_rad_m_ul(+270.0 * d2r, +45.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+271.0 * d2r, +45.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+270.0 * d2r, +44.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+271.0 * d2r, +44.0 * d2r, +3000);
            Eigen::Vector3d B_n_nT_ul(17455.7,  -764.6, 52192.7);
            Eigen::Vector3d B_n_nT_ur(17449.5, -1019.6, 52158.5);
            Eigen::Vector3d B_n_nT_ll(17983.9,  -739.1, 51571.5);
            Eigen::Vector3d B_n_nT_lr(17977.2,  -997.0, 51538.4);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        case env::logic::zone_pntt:
            Pmag = new env::mag_constant(Eigen::Vector3d(20677.1, -346.8, 47696.0), dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break;
        case env::logic::zone_river: {
            env::geodetic_coord x_gdt_rad_m_ul(+271.0 * d2r, +35.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ur(+272.0 * d2r, +35.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_ll(+271.0 * d2r, +34.0 * d2r, +3000);
            env::geodetic_coord x_gdt_rad_m_lr(+272.0 * d2r, +34.0 * d2r, +3000);

            Eigen::Vector3d B_n_nT_ul(22211.8,  -837.9, 44581.8);
            Eigen::Vector3d B_n_nT_ur(22204.8, -1110.6, 44539.6);
            Eigen::Vector3d B_n_nT_ll(22603.9,  -816.6, 43736.2);
            Eigen::Vector3d B_n_nT_lr(22595.3, -1090.1, 43694.4);
            Pmag = new env::mag_linear(B_n_nT_ul, x_gdt_rad_m_ul,
                                       B_n_nT_ur, x_gdt_rad_m_ur,
                                       B_n_nT_ll, x_gdt_rad_m_ll,
                                       B_n_nT_lr, x_gdt_rad_m_lr,
                                       dist_normal, gen, sigma_B_nT_north, sigma_B_nT_east, sigma_B_nT_down);
            break; }
        default:
            throw std::runtime_error("Location not available");
    }
    return Pmag;
}
/* return pointer to magnetic object based on location, realism ID, normal distribution, and seed generator. */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS MAG_CONSTANT
// ==================
// ==================

env::mag_constant::mag_constant(const Eigen::Vector3d& B_n_nT_model)
: _B_n_nT_model(B_n_nT_model), _B_n_nT_truth(B_n_nT_model) {
}
/* constructor based on NED magnetic field. Model and real magnetic fields both coincide with input. */

env::mag_constant::mag_constant(const Eigen::Vector3d& B_n_nT_model, std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen, const double& sigma_B_nT_north, const double& sigma_B_nT_east, const double& sigma_B_nT_down)
: _B_n_nT_model(B_n_nT_model) {
    _B_n_nT_truth = _B_n_nT_model + Eigen::Vector3d(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down* dist_normal(gen));
}
/* constructor based on NED magnetic field, normal distribution, seed generator, and three NED standard deviations. */

Eigen::Vector3d env::mag_constant::compute_B_n_nT_model(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    return _B_n_nT_model;
}
/* computes the MODEL Earth magnetic field based on time and geodetic coordinates */

Eigen::Vector3d env::mag_constant::compute_B_n_nT_truth(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    return _B_n_nT_truth;
}
/* computes the REAL Earth magnetic field based on time and geodetic coordinates */

void env::mag_constant::create_text(std::ostream& Ostream, const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    Ostream << std::endl << "MAGNETISM (CONSTANT) ERROR [model-truth-diff]:" << std::endl << std::endl;
    Ostream << "North [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(0)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_truth(0)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(0) - _B_n_nT_truth(0)
            << std::endl;
    Ostream << "East  [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(1)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_truth(1)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(1) - _B_n_nT_truth(1)
            << std::endl;
    Ostream << "Down  [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(2)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_truth(2)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << _B_n_nT_model(2) - _B_n_nT_truth(2)
            << std::endl;
}
/* describe magnetic model and realism in stream */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS MAG_LINEAR
// ================
// ================

env::mag_linear::mag_linear(const Eigen::Vector3d& B_n_nT_model_ul, const env::geodetic_coord& x_gdt_rad_m_ul,
                            const Eigen::Vector3d& B_n_nT_model_ur, const env::geodetic_coord& x_gdt_rad_m_ur,
                            const Eigen::Vector3d& B_n_nT_model_ll, const env::geodetic_coord& x_gdt_rad_m_ll,
                            const Eigen::Vector3d& B_n_nT_model_lr, const env::geodetic_coord& x_gdt_rad_m_lr)
: _B_n_nT_model_ll(B_n_nT_model_ll),
  _B_n_nT_model_lr(B_n_nT_model_lr),
  _B_n_nT_model_ul_ll(B_n_nT_model_ul - B_n_nT_model_ll),
  _B_n_nT_model_ur_lr(B_n_nT_model_ur - B_n_nT_model_lr),
  _phi_rad_low(x_gdt_rad_m_ll.get_phi_rad()),
  _lambda_rad_left(x_gdt_rad_m_ll.get_lambda_rad()),
  _Delta_phi_rad(x_gdt_rad_m_ul.get_phi_rad() - x_gdt_rad_m_ll.get_phi_rad()),
  _Delta_lambda_rad(x_gdt_rad_m_lr.get_lambda_rad() - x_gdt_rad_m_ll.get_lambda_rad()),
  _B_n_nT_truth_ll(B_n_nT_model_ll),
  _B_n_nT_truth_lr(B_n_nT_model_lr),
  _B_n_nT_truth_ul_ll(B_n_nT_model_ul - B_n_nT_model_ll),
  _B_n_nT_truth_ur_lr(B_n_nT_model_ur - B_n_nT_model_lr) {
    if (std::fabs(x_gdt_rad_m_ul.get_lambda_rad() - x_gdt_rad_m_ll.get_lambda_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong longitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ur.get_lambda_rad() - x_gdt_rad_m_lr.get_lambda_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong longitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ul.get_phi_rad() - x_gdt_rad_m_ur.get_phi_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong latitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ll.get_phi_rad() - x_gdt_rad_m_lr.get_phi_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong latitudes in magnetic model.");
    }
//    std::cout << "MAGNETIC FIELD ERROR (model - truth) [nT]: " << std::endl << std::fixed << std::setprecision(1) << std::showpos
//              << std::setw(7) << 0
//              << std::setw(7) << 0
//              << std::setw(7) << 0
//              << std::endl << std::endl;
}
/* constructor based on NED magnetic field and geodetic coordinates at four points forming square. Model and real magnetic
 * fields both coincide with each other and input. Altitude neglected. Linear interpolation in between and also outside square */

env::mag_linear::mag_linear(const Eigen::Vector3d& B_n_nT_model_ul, const env::geodetic_coord& x_gdt_rad_m_ul,
                            const Eigen::Vector3d& B_n_nT_model_ur, const env::geodetic_coord& x_gdt_rad_m_ur,
                            const Eigen::Vector3d& B_n_nT_model_ll, const env::geodetic_coord& x_gdt_rad_m_ll,
                            const Eigen::Vector3d& B_n_nT_model_lr, const env::geodetic_coord& x_gdt_rad_m_lr,
                            std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen,
                            const double& sigma_B_nT_north, const double& sigma_B_nT_east, const double& sigma_B_nT_down)
: _B_n_nT_model_ll(B_n_nT_model_ll),
  _B_n_nT_model_lr(B_n_nT_model_lr),
  _B_n_nT_model_ul_ll(B_n_nT_model_ul - B_n_nT_model_ll),
  _B_n_nT_model_ur_lr(B_n_nT_model_ur - B_n_nT_model_lr),
  _phi_rad_low(x_gdt_rad_m_ll.get_phi_rad()),
  _lambda_rad_left(x_gdt_rad_m_ll.get_lambda_rad()),
  _Delta_phi_rad(x_gdt_rad_m_ul.get_phi_rad() - x_gdt_rad_m_ll.get_phi_rad()),
  _Delta_lambda_rad(x_gdt_rad_m_lr.get_lambda_rad() - x_gdt_rad_m_ll.get_lambda_rad())
{
    if (std::fabs(x_gdt_rad_m_ul.get_lambda_rad() - x_gdt_rad_m_ll.get_lambda_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong longitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ur.get_lambda_rad() - x_gdt_rad_m_lr.get_lambda_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong longitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ul.get_phi_rad() - x_gdt_rad_m_ur.get_phi_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong latitudes in magnetic model.");
    }
    if (std::fabs(x_gdt_rad_m_ll.get_phi_rad() - x_gdt_rad_m_lr.get_phi_rad()) > math::constant::EPS()) {
        throw std::runtime_error("Wrong latitudes in magnetic model.");
    }
    Eigen::Vector3d error_B_nT(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down* dist_normal(gen));

//    std::cout << "MAGNETIC FIELD ERROR (model - truth) [nT]: " << std::endl << std::fixed << std::setprecision(1) << std::showpos
//              << std::setw(7) << - error_B_nT(0)
//              << std::setw(7) << - error_B_nT(1)
//              << std::setw(7) << - error_B_nT(2)
//              << std::endl << std::endl;

    _B_n_nT_truth_ll    = B_n_nT_model_ll + error_B_nT;
    _B_n_nT_truth_lr    = B_n_nT_model_lr + error_B_nT;
    _B_n_nT_truth_ul_ll = B_n_nT_model_ul + error_B_nT - _B_n_nT_truth_ll;
    _B_n_nT_truth_ur_lr = B_n_nT_model_ur + error_B_nT - _B_n_nT_truth_lr;

    // PREVIOUSLY DIFFERENT ERROR IN EACH CORNER, NOW THE SAME
    //_B_n_nT_truth_ll    = B_n_nT_model_ll + Eigen::Vector3d(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down * dist_normal(gen));
    //_B_n_nT_truth_lr    = B_n_nT_model_lr + Eigen::Vector3d(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down * dist_normal(gen));
    //_B_n_nT_truth_ul_ll = B_n_nT_model_ul + Eigen::Vector3d(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down * dist_normal(gen)) - _B_n_nT_truth_ll;
    //_B_n_nT_truth_ur_lr = B_n_nT_model_ur + Eigen::Vector3d(sigma_B_nT_north * dist_normal(gen), sigma_B_nT_east * dist_normal(gen), sigma_B_nT_down * dist_normal(gen)) - _B_n_nT_truth_lr;
}
/* constructor based on NED magnetic field and geodetic coordinates at four points forming square, normal distribution,
 * seed generator, and three NED standard deviations. Altitude neglected. Linear interpolation in between and also outside square */

Eigen::Vector3d env::mag_linear::compute_B_n_nT_model(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    double ratio1 = (x_gdt_rad_m.get_phi_rad() - _phi_rad_low) / _Delta_phi_rad;
    double ratio2 = (x_gdt_rad_m.get_lambda_rad() - _lambda_rad_left) / _Delta_lambda_rad;
    return (1.0 - ratio2) * (_B_n_nT_model_ll + ratio1 * _B_n_nT_model_ul_ll) + ratio2 * (_B_n_nT_model_lr + ratio1 * _B_n_nT_model_ur_lr);
}
/* computes the MODEL Earth magnetic field based on time and geodetic coordinates */

Eigen::Vector3d env::mag_linear::compute_B_n_nT_truth(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    double ratio1 = (x_gdt_rad_m.get_phi_rad() - _phi_rad_low) / _Delta_phi_rad;
    double ratio2 = (x_gdt_rad_m.get_lambda_rad() - _lambda_rad_left) / _Delta_lambda_rad;
    return (1.0 - ratio2) * (_B_n_nT_truth_ll + ratio1 * _B_n_nT_truth_ul_ll) + ratio2 * (_B_n_nT_truth_lr + ratio1 * _B_n_nT_truth_ur_lr);
}
/* computes the REAL Earth magnetic field based on time and geodetic coordinates */

void env::mag_linear::create_text(std::ostream& Ostream, const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const {
    Eigen::Vector3d B_n_nT_model = this->compute_B_n_nT_model(t_sec, x_gdt_rad_m);
    Eigen::Vector3d B_n_nT_truth = this->compute_B_n_nT_truth(t_sec, x_gdt_rad_m);

    Ostream << std::endl << "MAGNETISM (LINEAR) ERROR [model-truth-diff]:" << std::endl << std::endl;
    Ostream << "North [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(0)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_truth(0)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(0) - B_n_nT_truth(0)
            << std::endl;
    Ostream << "East  [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(1)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_truth(1)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(1) - B_n_nT_truth(1)
            << std::endl;
    Ostream << "Down  [nT]: "
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(2)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_truth(2)
            << std::fixed << std::setw(9) << std::setprecision(1) << std::showpos << B_n_nT_model(2) - B_n_nT_truth(2)
            << std::endl;
}
/* describe magnetic model and realism in stream */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////








