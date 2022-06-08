#include "sti.h"
#include "trj_truth.h"
#include "../guid/guid.h"
#include "ang/tools.h"
#include "ang/rotate/euler.h"
#include "env/earth.h"
#include "env/speed.h"
#include "env/geo.h"
#include "env/offsets.h"
#include "env/wind.h"
#include "env/turb.h"
#include "env/atm.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS INITIAL CONDITIONS STI
// ============================
// ============================

st::sti::sti(const sti& Osti)
: _t_sec(Osti._t_sec), _m_kg(Osti._m_kg), _vtas_mps(Osti._vtas_mps),
  _lambda_deg(Osti._lambda_deg), _lambda_rad(Osti._lambda_rad), _phi_deg(Osti._phi_deg),   _phi_rad(Osti._phi_rad), _h_m(Osti._h_m), _Hp_m(Osti._Hp_m),
  _theta_deg(Osti._theta_deg),   _theta_rad(Osti._theta_rad),   _xi_deg(Osti._xi_deg),     _xi_rad(Osti._xi_rad),
  _psi_deg(Osti._psi_deg),       _psi_rad(Osti._psi_rad),       _chi_deg(Osti._chi_deg),   _chi_rad(Osti._chi_rad),
  _alpha_deg(Osti._alpha_deg),   _alpha_rad(Osti._alpha_rad),   _beta_deg(Osti._beta_deg), _beta_rad(Osti._beta_rad),
  _w_nbb_rps(Osti._w_nbb_rps), _hground_m(Osti._hground_m),
  _delta_control(Osti._delta_control) {
}
/* copy constructor */

void st::sti::complete_misc(const double& t_sec, const double& m_kg, const double& vtas_mps, const double& alpha_deg, const double& beta_deg, const ang::so3_tangent& w_nbb_rps, const double& hground_m) {
    _t_sec      = t_sec;
    _m_kg       = m_kg;
    _vtas_mps   = vtas_mps;
    _alpha_deg  = alpha_deg;
    _beta_deg   = beta_deg;
    _w_nbb_rps  = w_nbb_rps;
    _hground_m  = hground_m;

    _alpha_rad  = _alpha_deg * math::constant::D2R();
    _beta_rad   = _beta_deg * math::constant::D2R();
}
/* add time, mass, true airspeed, angle of attack, angle of sideslip, aircraft rotation speed, and ground altitude */

void st::sti::complete_control(const double& delta_thr, const double& delta_elv, const double& delta_ail, const double& delta_rud) {
    _delta_control << delta_thr, delta_elv, delta_ail, delta_rud;
}
/* add control parameters (throttle, elevator, ailerons, rudder) */

st::sti* st::sti::create_sti(const control::guid& Oguid, env::logic::ZONE_ID zone_id, const double& t_sec_init, const unsigned short& case_guid, const unsigned short& seed) {
    switch (Oguid.get_sti_id()) {
        case st::logic::sti_h3000_tas25_psi00: {
            auto Psti = new st::sti_h_psi();
            Psti->complete_misc(t_sec_init, 19.5, 25.0, 3.0, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_h(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, h_m
            Psti->complete_psi_theta_xi(0.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_h3000_tas30_psi00: {
            auto Psti = new st::sti_h_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_h(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, h_m
            Psti->complete_psi_theta_xi(0.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_h3000_tas30_psi30: {
            auto Psti = new st::sti_h_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_h(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, h_m
            Psti->complete_psi_theta_xi(30.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_psi00: {
            auto Psti = new st::sti_Hp_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_psi_theta_xi(0.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_psi90: {
            auto Psti = new st::sti_Hp_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_psi_theta_xi(90.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_chi00: {
            auto Psti = new st::sti_Hp_chi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_chi_theta_xi(0.0, 0.0, 0.0); // chi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_chi90: {
            auto Psti = new st::sti_Hp_chi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_chi_theta_xi(90.0, 0.0, 0.0); // chi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_chiXX: {
            auto Psti = new st::sti_Hp_chi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_chi_theta_xi(Oguid.get().front()->get_guid_val(control::logic::cntr_AIL), 0.0, 0.0); // chi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_Hp3000_tas30_psiXX: {
            auto Psti = new st::sti_Hp_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 3000.0); // lambda_deg, phi_deg, Hp_m
            Psti->complete_psi_theta_xi(Oguid.get().front()->get_guid_val(control::logic::cntr_AIL), 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_HpXX_tasXX_chiXX: {
            auto Psti = new st::sti_Hp_chi();
            Psti->complete_misc(t_sec_init, 19.5, Oguid.get().front()->get_guid_val(control::logic::cntr_THR), -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), Oguid.get().front()->get_guid_val(control::logic::cntr_ELV)); // lambda_deg, phi_deg, Hp_m
            Psti->complete_chi_theta_xi(Oguid.get().front()->get_guid_val(control::logic::cntr_AIL), 0.0, 0.0); // chi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_HpXX_tasXX_psiXX: {
            auto Psti = new st::sti_Hp_psi();
            Psti->complete_misc(t_sec_init, 19.5, Oguid.get().front()->get_guid_val(control::logic::cntr_THR), -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_Hp(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), Oguid.get().front()->get_guid_val(control::logic::cntr_ELV)); // lambda_deg, phi_deg, Hp_m
            Psti->complete_psi_theta_xi(Oguid.get().front()->get_guid_val(control::logic::cntr_AIL), 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_rozas: {
            auto Psti = new st::sti_h_psi();
            Psti->complete_misc(t_sec_init, 19.5, 30.0, -0.2, 0.0, ang::so3_tangent(Eigen::Vector3d::Zero()), st::sti::get_hground_m_init(zone_id, case_guid, seed)); // t_sec, m_kg, vtas_mps, alpha_deg, beta_deg, w_nbb_rps (3), hground_m
            Psti->complete_lambda_phi_h(st::sti::get_lambda_deg_init(zone_id, case_guid, seed), st::sti::get_phi_deg_init(zone_id, case_guid, seed), 1300.0); // lambda_deg, phi_deg, h_m (1300.0, 2300.0, 3300.0)
            Psti->complete_psi_theta_xi(40.0, 0.0, 0.0); // psi_deg, theta_deg, xi_deg
            Psti->complete_control(0.3, 0.0, 0.0, 0.0); // control parameters (throttle, elevator, ailerons, rudder)
            return Psti;
        }
        case st::logic::sti_size:
        default:
            throw std::runtime_error("Initial conditions case not available");
    }
}
/* creates initial conditions based on guidance and initial time */

double st::sti::get_lambda_deg_init(env::logic::ZONE_ID zone_id, const unsigned short& case_guid, const unsigned short& seed) {
    switch (zone_id) {
        case env::logic::zone_default:
            return 270.0;
        case env::logic::zone_desert: {
            if (case_guid == 1) {
                if ((seed == 7) || (seed == 8) || (seed == 30) || (seed == 47) || (seed == 77)) {
                    return 248.037484;
                }
                else if (seed == 92) {
                    return 247.997637;
                }
                else {
                    return 248.001185;
                }

            } else {
                return 248.001185;
            }
        }
        case env::logic::zone_urban:
            return 241.799731;
        case env::logic::zone_everglades:
            return 279.088834;
        case env::logic::zone_forest: {
            if (case_guid == 1) {
                if ((seed == 13) || (seed == 44) || (seed == 66)) {
                    return 287.507253;
                }
                else if ((seed == 23) || (seed == 26) || (seed == 76)) {
                    return 286.391731;
                }
                else if ((seed == 68) || (seed == 72)) {
                    return 288.589289;
                }
                else {
                    return 287.490805;
                }
            }
            else {
                return 287.490805;
            }
        }
        case env::logic::zone_farm: {
            if (case_guid == 1) {
                if ((seed == 15) || (seed == 16) || (seed == 19) || (seed == 31) || (seed == 50) || (seed == 79) || (seed == 87)) {
                    return 273.029625;
                }
                else if (seed == 42) {
                    return 272.091038;
                }
                else if ((seed == 59) || (seed == 91)) {
                    return 273.043256;
                }
                else if (seed == 77) {
                    return 272.132580;
                }
                else {
                    return 272.122371;
                }
            }
            else {
                return 272.122371;
            }
        }
        case env::logic::zone_mix: {
            if ((case_guid == 1) && (seed == 4)) {
                return 271.023307;
            }
            else {
                return 270.984538;
            }
        }
        case env::logic::zone_mexico:
            return 260.0;
        case env::logic::zone_colombia:
            return 290.0;
        case env::logic::zone_argentina:
            return 300.0;
        case env::logic::zone_canada:
            return 255.0;
        case env::logic::zone_rozas:
            return 352.514219;
        case env::logic::zone_wisconsin:
            return 270.0;
        case env::logic::zone_pntt:
            return 268.883747;
        case env::logic::zone_river:
            return 271.908733;
        default:
            throw std::runtime_error("Zone not available");
    }
}
/* get initial longitude based on location */

double st::sti::get_phi_deg_init(env::logic::ZONE_ID zone_id, const unsigned short& case_guid, const unsigned short& seed) {
    switch (zone_id) {
        case env::logic::zone_default:
            return 45.0;
        case env::logic::zone_desert: {
            if (case_guid == 1) {
                if ((seed == 7) || (seed == 8) || (seed == 30) || (seed == 47) || (seed == 77)) {
                    return 32.385922;
                }
                else if (seed == 92) {
                    return 32.081079;
                }
                else {
                    return 32.157903;
                }

            } else {
                return 32.157903;
            }
        }
        case env::logic::zone_urban:
            return 33.924426;
        case env::logic::zone_everglades:
            return 25.855172;
        case env::logic::zone_forest: {
            if (case_guid == 1) {
                if ((seed == 13) || (seed == 44) || (seed == 66)) {
                    return 42.748124;
                }
                else if ((seed == 23) || (seed == 26) || (seed == 76)) {
                    return 43.332651;
                }
                else if ((seed == 68) || (seed == 72)) {
                    return 42.911662;
                }
                else {
                    return 43.354486;
                }
            }
            else {
                return 43.354486;
            }
        }
        case env::logic::zone_farm: {
            if (case_guid == 1) {
                if ((seed == 15) || (seed == 16) || (seed == 19) || (seed == 31) || (seed == 50) || (seed == 79) || (seed == 87)) {
                    return 38.245024;
                }
                else if (seed == 42) {
                    return 38.859908;
                }
                else if ((seed == 59) || (seed == 91)) {
                    return 38.085929;
                }
                else if (seed == 77) {
                    return 38.817388;
                }
                else {
                    return 38.865625;
                }
            }
            else {
                return 38.865625;
            }
        }
        case env::logic::zone_mix: {
            if ((case_guid == 1) && (seed == 4)) {
                return 34.846473;
            }
            else {
                return 34.720636;
            }
        }
        case env::logic::zone_mexico:
            return 20.0;
        case env::logic::zone_colombia:
            return 5.0;
        case env::logic::zone_argentina:
            return -25.0;
        case env::logic::zone_canada:
            return 52.0;
        case env::logic::zone_rozas:
            return 43.094976;
        case env::logic::zone_wisconsin:
            return 45.0;
        case env::logic::zone_pntt:
            return 38.731938;
        case env::logic::zone_river:
            return 35.012125;
        default:
            throw std::runtime_error("Zone not available");
    }
}
/* get initial latitude based on location */

double st::sti::get_hground_m_init(env::logic::ZONE_ID zone_id, const unsigned short& case_guid, const unsigned short& seed) {
    switch (zone_id) {
        case env::logic::zone_default:
            return 0.0; // not valid for images
        case env::logic::zone_desert: {
            if (case_guid == 1) {
                if ((seed == 7) || (seed == 8) || (seed == 30) || (seed == 47) || (seed == 77)) {
                    return 517.0;
                }
                else if (seed == 92) {
                    return 691.0;
                }
                else {
                    return 661.0;
                }
            }
            else {
                return 661.0;
            }
        }
        case env::logic::zone_urban:
            return 26.0;
        case env::logic::zone_everglades:
            return 10.0;
        case env::logic::zone_forest: {
            if (case_guid == 1) {
                if ((seed == 13) || (seed == 44) || (seed == 66)) {
                    return 88.0;
                }
                else if ((seed == 23) || (seed == 26) || (seed == 76)) {
                    return 97.0;
                }
                else if ((seed == 68) || (seed == 72)) {
                    return 78.0;
                }
                else {
                    return 200.0;
                }
            }
            else {
                return 200.0;
            }
        }
        case env::logic::zone_farm: {
            if (case_guid == 1) {
                if ((seed == 15) || (seed == 16) || (seed == 19) || (seed == 31) || (seed == 50) || (seed == 79) || (seed == 87)) {
                    return 155.0;
                }
                else if (seed == 42) {
                    return 138.0;
                }
                else if ((seed == 59) || (seed == 91)) {
                    return 135.0;
                }
                else if (seed == 77) {
                    return 136.0;
                }
                else {
                    return 144.0;
                }
            }
            else {
                return 144.0;
            }
        }
        case env::logic::zone_mix: {
            if ((case_guid == 1) && (seed == 4)) {
                return 158.0;
            }
            return 133.0;
        }
        case env::logic::zone_mexico:
            return 2914.0;
        case env::logic::zone_colombia:
            return 111.9;
        case env::logic::zone_argentina:
            return 115.2;
        case env::logic::zone_canada:
            return 561.7;
        case env::logic::zone_rozas:
            return 426.0;
        case env::logic::zone_wisconsin:
            return 0; // not valid for images
        case env::logic::zone_pntt:
            return 210.0;
        case env::logic::zone_river:
            return 221.0;
        default:
            throw std::runtime_error("Zone not available");
    }
}
/* get initial altitude based on location (only for guidance initial altitude in some scenarios) */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS INITIAL CONDITIONS STI_h_psi
// ==================================
// ==================================

st::sti_h_psi::sti_h_psi(const sti_h_psi& Osti)
: sti(*this) {
}
/* copy constructor */

void st::sti_h_psi::complete_lambda_phi_h(const double& lambda_deg, const double& phi_deg, const double& h_m) {
    _lambda_deg = lambda_deg;
    _phi_deg    = phi_deg;
    _h_m        = h_m;
    _Hp_m       = std::nan("");

    ang::tools::correct_longitude_deg(_lambda_deg);

    _lambda_rad = _lambda_deg * math::constant::D2R();
    _phi_rad    = _phi_deg * math::constant::D2R();
}
/* add geodetic position (with geometric altitude) */

void st::sti_h_psi::complete_psi_theta_xi(const double& psi_deg, const double& theta_deg, const double& xi_deg) {
    _psi_deg    = psi_deg;
    _theta_deg  = theta_deg;
    _xi_deg     = xi_deg;
    ang::tools::correct_yaw_deg(_psi_deg);

    _psi_rad    = _psi_deg * math::constant::D2R();
    _theta_rad  = _theta_deg * math::constant::D2R();
    _xi_rad     = _xi_deg * math::constant::D2R();

    _chi_deg    = std::nan("");
    _chi_rad    = std::nan("");
}
/* add body yaw, body pitch, body roll */

void st::sti_h_psi::complete_and_fill_up_st_truth(st::st_truth& Ost_init, const env::earth& Oearth) {
    Ost_init.get_t_sec() = _t_sec; // time

    Ost_init.get_x_gdt_b_rad_m().get_lambda_rad() = _lambda_rad; // geodetic coordinates
    Ost_init.get_x_gdt_b_rad_m().get_phi_rad() = _phi_rad;
    Ost_init.get_x_gdt_b_rad_m().get_h_m() = _h_m;

    ang::euler euler_nedbfs_rad(_psi_rad, _theta_rad, _xi_rad);
    Ost_init.get_q_nb() = euler_nedbfs_rad; // quaternion

    Ost_init.get_m_kg() = _m_kg; // mass

    Ost_init.get_delta_control() << _delta_control; // control parameters

    double N_m = Oearth.get_geo().radius_vert(_phi_rad);
    double M_m = Oearth.get_geo().radius_mer(_phi_rad, N_m);

    double H_m = Oearth.get_geo().htoH(_h_m, _phi_rad);
    double DeltaT_degK = Oearth.get_offsets().compute_DeltaT_degK(_t_sec, _lambda_rad, _phi_rad);
    double Deltap_pa = Oearth.get_offsets().compute_Deltap_pa(_t_sec, _lambda_rad, _phi_rad);
    env::atm Oatm(DeltaT_degK, Deltap_pa);
    _Hp_m = Oatm.H2Hp(H_m, false);

    ang::euler euler_wb(-_beta_rad, _alpha_rad, 0.);
    Eigen::Vector3d vtas_b_mps = env::speed::vtas2vtasbfs(_vtas_mps, euler_wb);

    Eigen::Vector3d vlf_n_mps = Oearth.get_wind().compute_wind_ned(_t_sec, Ost_init.get_x_gdt_b_rad_m());
    Eigen::Vector3d vlf_b_mps = Ost_init.get_q_nb() / vlf_n_mps;

    double hground_m = 0.; ////////////////////////////////////////////
    double height_m  = _h_m - hground_m;
    Eigen::Vector3d vhf_b_mps = Oearth.get_turb().compute_wind_high_frequency_twomils_bfs(_t_sec, height_m, ang::euler::obtain_yaw_forward(vlf_n_mps), Ost_init.get_q_nb());

    Ost_init.get_v_b_mps() = vtas_b_mps + vlf_b_mps + vhf_b_mps; // speed
    Eigen::Vector3d v_n_mps = Ost_init.get_q_nb() * Ost_init.get_v_b_mps();

    _chi_rad   = ang::euler::obtain_yaw_forward(v_n_mps);
    _chi_deg   = _chi_rad * math::constant::R2D();
    _gamma_rad = ang::euler::obtain_pitch_forward(v_n_mps);
    _gamma_deg = _gamma_rad * math::constant::R2D();

    ang::so3_tangent w_enn_rps = Oearth.get_geo().compute_wenn_rps(v_n_mps, N_m, M_m, _phi_rad, _h_m);
    Ost_init.get_w_ibb_rps() = _w_nbb_rps() + Ost_init.get_q_nb() / (w_enn_rps() + Oearth.get_geo().compute_wien_rps(_phi_rad)());
    Ost_init.get_xi_ebb_mrps().set_vi(Ost_init.get_v_b_mps());
    Ost_init.get_xi_ebb_mrps().set_w(ang::so3_tangent(_w_nbb_rps() + Ost_init.get_q_nb() / w_enn_rps()));

    ang::rodrigues q_en(Ost_init.get_x_gdt_b_rad_m().obtain_euler_en());
    Ost_init.get_z_eb() = ang::dual(q_en * Ost_init.get_q_nb(), Oearth.get_geo().geodetic2cartesian(Ost_init.get_x_gdt_b_rad_m(), N_m)());

    double v_mps   = Ost_init.get_v_b_mps().norm();
    double vlf_mps = vlf_n_mps.norm();
    double vhf_mps = vhf_b_mps.norm();

    std::cout << std::endl << "INITIAL CONDITIONS:" << std::endl << std::endl
              << "t [sec]:          " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _t_sec << std::endl
              << "m [kg]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _m_kg << std::endl
              << "h [m]:            " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _h_m << std::endl
              << "Hp [m]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _Hp_m << std::endl
              << "vtas [mps]        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _vtas_mps << std::endl
              << "vlf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vlf_mps << std::endl
              << "vhf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vhf_mps << std::endl
              << "v [mps]           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_mps << std::endl
              << "DeltaT [degK]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << DeltaT_degK << std::endl
              << "psi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _psi_deg << std::endl
              << "theta [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _theta_deg << std::endl
              << "xi [deg]:         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _xi_deg << std::endl
              << "chi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _chi_deg << std::endl
              << "gamma [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _gamma_deg << std::endl
              << "alpha [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _alpha_deg << std::endl
              << "beta [deg]:       " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _beta_deg << std::endl
              << "delta_thr [-]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(0) << std::endl
              << "delta_elv [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(1) << std::endl
              << "delta_ail [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(2) << std::endl
              << "delta_rud [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(3) << std::endl
              << "vtas b i [mps]:   " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[0] << std::endl
              << "vtas b ii [mps]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[1] << std::endl
              << "vtas b iii [mps]: " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[2] << std::endl
              << "v b i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[0] << std::endl
              << "v b ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[1] << std::endl
              << "v b iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[2] << std::endl
              << "v n i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[0] << std::endl
              << "v n ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[1] << std::endl
              << "v n iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[2] << std::endl
              << std::endl;
}
/* fill up state vector with initial conditions */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS INITIAL CONDITIONS STI_h_chi
// ===================================
// ===================================

st::sti_h_chi::sti_h_chi(const sti_h_chi& Osti)
: sti(*this) {
}
/* copy constructor */

void st::sti_h_chi::complete_lambda_phi_h(const double& lambda_deg, const double& phi_deg, const double& h_m) {
    _lambda_deg = lambda_deg;
    _phi_deg    = phi_deg;
    _h_m        = h_m;
    _Hp_m       = std::nan("");

    ang::tools::correct_longitude_deg(_lambda_deg);

    _lambda_rad = _lambda_deg * math::constant::D2R();
    _phi_rad    = _phi_deg * math::constant::D2R();
}
/* add geodetic position */

void st::sti_h_chi::complete_chi_theta_xi(const double& chi_deg, const double& theta_deg, const double& xi_deg) {
    _chi_deg    = chi_deg;
    _theta_deg  = theta_deg;
    _xi_deg     = xi_deg;

    ang::tools::correct_yaw_deg(_chi_deg);
    _chi_rad    = _chi_deg * math::constant::D2R();
    _theta_rad  = _theta_deg * math::constant::D2R();
    _xi_rad     = _xi_deg * math::constant::D2R();

    _psi_deg    = std::nan("");
    _psi_rad    = std::nan("");
}
/* add body yaw, body pitch, body roll */

void st::sti_h_chi::complete_and_fill_up_st_truth(st::st_truth& Ost_init, const env::earth& Oearth) {
    Ost_init.get_t_sec() = _t_sec; // time

    Ost_init.get_m_kg() = _m_kg; // mass

    Ost_init.get_delta_control() << _delta_control; // control parameters

    double H_m = Oearth.get_geo().htoH(_h_m, _phi_rad);
    double DeltaT_degK = Oearth.get_offsets().compute_DeltaT_degK(_t_sec, _lambda_rad, _phi_rad);
    double Deltap_pa = Oearth.get_offsets().compute_Deltap_pa(_t_sec, _lambda_rad, _phi_rad);
    env::atm Oatm(DeltaT_degK, Deltap_pa);
    _Hp_m = Oatm.H2Hp(H_m, false);

    Ost_init.get_x_gdt_b_rad_m().get_lambda_rad() = _lambda_rad; // geodetic coordinates
    Ost_init.get_x_gdt_b_rad_m().get_phi_rad() = _phi_rad;
    Ost_init.get_x_gdt_b_rad_m().get_h_m() = _h_m;

    double N_m = Oearth.get_geo().radius_vert(_phi_rad);
    double M_m = Oearth.get_geo().radius_mer(_phi_rad, N_m);

    ang::euler euler_wb_rad(-_beta_rad, _alpha_rad, 0.);
    Eigen::Vector3d vtas_b_mps = env::speed::vtas2vtasbfs(_vtas_mps, euler_wb_rad);

    Eigen::Vector3d vlf_n_mps = Oearth.get_wind().compute_wind_ned(_t_sec, Ost_init.get_x_gdt_b_rad_m());

    // Note that at high altitudes q_nb does not matter when obtaining turbulence. THIS IS ONLY VALID FOR HIGH ALTITUDES TODO
    double hground_m = 0.; ////////////////////////////////////////////
    double height_m  = _h_m - hground_m;
    ang::euler euler_nb_dummy(0., 0., 0.);
    ang::rodrigues q_nb_dummy(euler_nb_dummy);
    Eigen::Vector3d vhf_b_mps = Oearth.get_turb().compute_wind_high_frequency_twomils_bfs(_t_sec, height_m, ang::euler::obtain_yaw_forward(vlf_n_mps), q_nb_dummy);

    st::sti_h_chi::psi_u Func;
    std::vector<double> par(12);
    par[0]  = _chi_rad;
    par[1]  = _theta_rad;
    par[2]  = _xi_rad;
    par[3]  = vtas_b_mps(0);
    par[4]  = vtas_b_mps(1);
    par[5]  = vtas_b_mps(2);
    par[6]  = vlf_n_mps(0);
    par[7]  = vlf_n_mps(1);
    par[8]  = vlf_n_mps(2);
    par[9]  = vhf_b_mps(0);
    par[10] = vhf_b_mps(1);
    par[11] = vhf_b_mps(2);
    double psi1_rad = _chi_rad - math::constant::PI() / 18;
    double psi2_rad = _chi_rad + math::constant::PI() / 18;
    _psi_rad = Func.find_zero_secant(par, 0., psi1_rad, psi2_rad, 1e-10);
    ang::tools::correct_yaw_rad(_psi_rad);
    _psi_deg = _psi_rad * math::constant::R2D();

    ang::euler euler_nedbfs_rad(_psi_rad, _theta_rad, _xi_rad);
    Ost_init.get_q_nb() = euler_nedbfs_rad; // quaternion

    Eigen::Vector3d vlf_b_mps = Ost_init.get_q_nb() / vlf_n_mps;
    Ost_init.get_v_b_mps() = vtas_b_mps + vlf_b_mps + vhf_b_mps; // speed
    Eigen::Vector3d v_n_mps = Ost_init.get_q_nb() * Ost_init.get_v_b_mps();

    ang::so3_tangent w_enn_rps = Oearth.get_geo().compute_wenn_rps(v_n_mps, N_m, M_m, _phi_rad, _h_m);
    Ost_init.get_w_ibb_rps() = _w_nbb_rps() + Ost_init.get_q_nb() / (w_enn_rps() + Oearth.get_geo().compute_wien_rps(_phi_rad)());
    Ost_init.get_xi_ebb_mrps().set_vi(Ost_init.get_v_b_mps());
    Ost_init.get_xi_ebb_mrps().set_w(ang::so3_tangent(_w_nbb_rps() + Ost_init.get_q_nb() / w_enn_rps()));

    ang::rodrigues q_en(Ost_init.get_x_gdt_b_rad_m().obtain_euler_en());
    Ost_init.get_z_eb() = ang::dual(q_en * Ost_init.get_q_nb(), Oearth.get_geo().geodetic2cartesian(Ost_init.get_x_gdt_b_rad_m(), N_m)());

    _gamma_rad = ang::euler::obtain_pitch_forward(v_n_mps);
    _gamma_deg = _gamma_rad * math::constant::R2D();

    double v_mps   = Ost_init.get_v_b_mps().norm();
    double vlf_mps = vlf_n_mps.norm();
    double vhf_mps = vhf_b_mps.norm();

    std::cout << std::endl << "INITIAL CONDITIONS:" << std::endl << std::endl
              << "t [sec]:          " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _t_sec << std::endl
              << "m [kg]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _m_kg << std::endl
              << "h [m]:            " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _h_m << std::endl
              << "Hp [m]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _Hp_m << std::endl
              << "vtas [mps]        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _vtas_mps << std::endl
              << "vlf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vlf_mps << std::endl
              << "vhf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vhf_mps << std::endl
              << "v [mps]           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_mps << std::endl
              << "DeltaT [degK]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << DeltaT_degK << std::endl
              << "psi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _psi_deg << std::endl
              << "theta [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _theta_deg << std::endl
              << "xi [deg]:         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _xi_deg << std::endl
              << "chi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _chi_deg << std::endl
              << "gamma [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _gamma_deg << std::endl
              << "alpha [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _alpha_deg << std::endl
              << "beta [deg]:       " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _beta_deg << std::endl
              << "delta_thr [-]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(0) << std::endl
              << "delta_elv [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(1) << std::endl
              << "delta_ail [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(2) << std::endl
              << "delta_rud [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(3) << std::endl
              << "vtas b i [mps]:   " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[0] << std::endl
              << "vtas b ii [mps]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[1] << std::endl
              << "vtas b iii [mps]: " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[2] << std::endl
              << "v b i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[0] << std::endl
              << "v b ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[1] << std::endl
              << "v b iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[2] << std::endl
              << "v n i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[0] << std::endl
              << "v n ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[1] << std::endl
              << "v n iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[2] << std::endl
              << std::endl;
}
/* fill up state vector with initial conditions */

double st::sti_h_chi::compute_psi_diff(const double& psi_rad, const double& chi_rad, const double& theta_rad, const double& xi_rad,
                                       const double& vtas_bi_mps, const double& vtas_bii_mps, const double& vtas_biii_mps,
                                       const double& vwind_ni_mps, const double& vwind_nii_mps, const double& vwind_niii_mps,
                                       const double& vturb_bi_mps, const double& vturb_bii_mps, const double& vturb_biii_mps) {
    Eigen::Vector3d vtas_b_mps(vtas_bi_mps, vtas_bii_mps, vtas_biii_mps);
    Eigen::Vector3d vwind_n_mps(vwind_ni_mps, vwind_nii_mps, vwind_niii_mps);
    Eigen::Vector3d vturb_b_mps(vturb_bi_mps, vturb_bii_mps, vturb_biii_mps);
    ang::euler euler_nb(psi_rad, theta_rad, xi_rad);
    ang::dcm R_nb(euler_nb);
    Eigen::Vector3d vtas_n_mps = R_nb * vtas_b_mps;
    Eigen::Vector3d v_n_mps = vtas_n_mps + vwind_n_mps + R_nb * vturb_b_mps;
    double new_chi_rad = ang::euler::obtain_yaw_forward(v_n_mps);
    return ang::tools::angle_diff_rad(new_chi_rad, chi_rad);
}
/* function that returns difference between input body bearing psi and what can be computed from rest of inputs.
 * Difference should be zero. */

double st::sti_h_chi::psi_u::exec(const double& psi_rad, const std::vector<double>& par) {
    // par[0]  = chi_rad
    // par[1]  = theta_rad
    // par[2]  = xi_rad
    // par[3]  = vtas_bi_mps
    // par[4]  = vtas_bii_mps
    // par[5]  = vtas_biii_mps
    // par[6]  = vwind_ni_mps
    // par[7]  = vwind_nii_mps
    // par[8]  = vwind_niii_mps
    // par[9]  = vturb_bi_mps
    // par[10] = vturb_bii_mps
    // par[11] = vturb_biii_mps
    return st::sti_h_chi::compute_psi_diff(psi_rad, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11]);
}
/* encapsulates minimization function to obtain body bearing psi */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS INITIAL CONDITIONS STI_Hp_chi
// ===================================
// ===================================

st::sti_Hp_chi::sti_Hp_chi(const sti_Hp_chi& Osti)
: sti(*this) {
}
/* copy constructor */

void st::sti_Hp_chi::complete_lambda_phi_Hp(const double& lambda_deg, const double& phi_deg, const double& Hp_m) {
    _lambda_deg = lambda_deg;
    _phi_deg    = phi_deg;
    _Hp_m       = Hp_m;
    _h_m        = std::nan("");

    ang::tools::correct_longitude_deg(_lambda_deg);

    _lambda_rad = _lambda_deg * math::constant::D2R();
    _phi_rad    = _phi_deg * math::constant::D2R();
}
/* add pressure geodetic position (based on pressure altitude) */

void st::sti_Hp_chi::complete_chi_theta_xi(const double& chi_deg, const double& theta_deg, const double& xi_deg) {
    _chi_deg    = chi_deg;
    _theta_deg  = theta_deg;
    _xi_deg     = xi_deg;

    ang::tools::correct_yaw_deg(_chi_deg);
    _chi_rad    = _chi_deg * math::constant::D2R();
    _theta_rad  = _theta_deg * math::constant::D2R();
    _xi_rad     = _xi_deg * math::constant::D2R();

    _psi_deg    = std::nan("");
    _psi_rad    = std::nan("");
}
/* add body yaw, body pitch, body roll */

void st::sti_Hp_chi::complete_and_fill_up_st_truth(st::st_truth& Ost_init, const env::earth& Oearth) {
    Ost_init.get_t_sec() = _t_sec; // time

    Ost_init.get_m_kg() = _m_kg; // mass

    Ost_init.get_delta_control() << _delta_control; // control parameters

    double DeltaT_degK = Oearth.get_offsets().compute_DeltaT_degK(_t_sec, _lambda_rad, _phi_rad);
    double Deltap_pa = Oearth.get_offsets().compute_Deltap_pa(_t_sec, _lambda_rad, _phi_rad);
    env::atm Oatm(DeltaT_degK, Deltap_pa);
    double H_m = Oatm.Hp2H(_Hp_m);
    double Hp_m_bias = Oatm.H2Hp(H_m, false); // only purpose is to initialize function H2Hp for later use, do not remove
    _h_m = Oearth.get_geo().H2h(H_m, _phi_rad);

    Ost_init.get_x_gdt_b_rad_m().get_lambda_rad() = _lambda_rad; // geodetic coordinates
    Ost_init.get_x_gdt_b_rad_m().get_phi_rad() = _phi_rad;
    Ost_init.get_x_gdt_b_rad_m().get_h_m() = _h_m;

    double N_m = Oearth.get_geo().radius_vert(_phi_rad);
    double M_m = Oearth.get_geo().radius_mer(_phi_rad, N_m);

    ang::euler euler_wb_rad(-_beta_rad, _alpha_rad, 0.);
    Eigen::Vector3d vtas_b_mps = env::speed::vtas2vtasbfs(_vtas_mps, euler_wb_rad);

    Eigen::Vector3d vlf_n_mps = Oearth.get_wind().compute_wind_ned(_t_sec, Ost_init.get_x_gdt_b_rad_m());

    // Note that at high altitudes q_nb does not matter when obtaining turbulence. THIS IS ONLY VALID FOR HIGH ALTITUDES TODO
    double hground_m = 0.; ////////////////////////////////////////////
    double height_m  = _h_m - hground_m;
    ang::euler euler_nb_dummy(0., 0., 0.);
    ang::rodrigues q_nb_dummy(euler_nb_dummy);
    Eigen::Vector3d vhf_b_mps = Oearth.get_turb().compute_wind_high_frequency_twomils_bfs(_t_sec, height_m, ang::euler::obtain_yaw_forward(vlf_n_mps), q_nb_dummy);

    st::sti_Hp_chi::psi_u Func;
    std::vector<double> par(12);
    par[0]  = _chi_rad;
    par[1]  = _theta_rad;
    par[2]  = _xi_rad;
    par[3]  = vtas_b_mps(0);
    par[4]  = vtas_b_mps(1);
    par[5]  = vtas_b_mps(2);
    par[6]  = vlf_n_mps(0);
    par[7]  = vlf_n_mps(1);
    par[8]  = vlf_n_mps(2);
    par[9]  = vhf_b_mps(0);
    par[10] = vhf_b_mps(1);
    par[11] = vhf_b_mps(2);

    double psi1_rad = _chi_rad - math::constant::PI() / 18;
    double psi2_rad = _chi_rad + math::constant::PI() / 18;
    _psi_rad = Func.find_zero_secant(par, 0., psi1_rad, psi2_rad, 1e-10);
    ang::tools::correct_yaw_rad(_psi_rad);
    _psi_deg = _psi_rad * math::constant::R2D();

    ang::euler euler_nedbfs_rad(_psi_rad, _theta_rad, _xi_rad);
    Ost_init.get_q_nb() = euler_nedbfs_rad; // quaternion

    Eigen::Vector3d vlf_b_mps = Ost_init.get_q_nb() / vlf_n_mps;
    Ost_init.get_v_b_mps() = vtas_b_mps + vlf_b_mps + vhf_b_mps; // speed
    Eigen::Vector3d v_n_mps = Ost_init.get_q_nb() * Ost_init.get_v_b_mps();

    ang::so3_tangent w_enn_rps = Oearth.get_geo().compute_wenn_rps(v_n_mps, N_m, M_m, _phi_rad, _h_m);
    Ost_init.get_w_ibb_rps() = _w_nbb_rps() + Ost_init.get_q_nb() / (w_enn_rps() + Oearth.get_geo().compute_wien_rps(_phi_rad)());
    Ost_init.get_xi_ebb_mrps().set_vi(Ost_init.get_v_b_mps());
    Ost_init.get_xi_ebb_mrps().set_w(ang::so3_tangent(_w_nbb_rps() + Ost_init.get_q_nb() / w_enn_rps()));

    ang::rodrigues q_en(Ost_init.get_x_gdt_b_rad_m().obtain_euler_en());
    Ost_init.get_z_eb() = ang::dual(q_en * Ost_init.get_q_nb(), Oearth.get_geo().geodetic2cartesian(Ost_init.get_x_gdt_b_rad_m(), N_m)());

    _gamma_rad = ang::euler::obtain_pitch_forward(v_n_mps);
    _gamma_deg = _gamma_rad * math::constant::R2D();

    double v_mps   = Ost_init.get_v_b_mps().norm();
    double vlf_mps = vlf_n_mps.norm();
    double vhf_mps = vhf_b_mps.norm();

    std::cout << std::endl << "INITIAL CONDITIONS:" << std::endl << std::endl
              << "t [sec]:          " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _t_sec << std::endl
              << "m [kg]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _m_kg << std::endl
              << "h [m]:            " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _h_m << std::endl
              << "Hp [m]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _Hp_m << std::endl
              << "vtas [mps]        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _vtas_mps << std::endl
              << "vlf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vlf_mps << std::endl
              << "vhf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vhf_mps << std::endl
              << "v [mps]           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_mps << std::endl
              << "DeltaT [degK]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << DeltaT_degK << std::endl
              << "psi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _psi_deg << std::endl
              << "theta [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _theta_deg << std::endl
              << "xi [deg]:         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _xi_deg << std::endl
              << "chi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _chi_deg << std::endl
              << "gamma [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _gamma_deg << std::endl
              << "alpha [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _alpha_deg << std::endl
              << "beta [deg]:       " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _beta_deg << std::endl
              << "delta_thr [-]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(0) << std::endl
              << "delta_elv [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(1) << std::endl
              << "delta_ail [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(2) << std::endl
              << "delta_rud [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(3) << std::endl
              << "vtas b i [mps]:   " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[0] << std::endl
              << "vtas b ii [mps]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[1] << std::endl
              << "vtas b iii [mps]: " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[2] << std::endl
              << "v b i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[0] << std::endl
              << "v b ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[1] << std::endl
              << "v b iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[2] << std::endl
              << "v n i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[0] << std::endl
              << "v n ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[1] << std::endl
              << "v n iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[2] << std::endl
              << std::endl;
}
/* fill up state vector with initial conditions */

double st::sti_Hp_chi::compute_psi_diff(const double& psi_rad, const double& chi_rad, const double& theta_rad, const double& xi_rad,
                                        const double& vtas_bi_mps, const double& vtas_bii_mps, const double& vtas_biii_mps,
                                        const double& vwind_ni_mps, const double& vwind_nii_mps, const double& vwind_niii_mps,
                                        const double& vturb_bi_mps, const double& vturb_bii_mps, const double& vturb_biii_mps) {
    Eigen::Vector3d vtas_b_mps(vtas_bi_mps, vtas_bii_mps, vtas_biii_mps);
    Eigen::Vector3d vwind_n_mps(vwind_ni_mps, vwind_nii_mps, vwind_niii_mps);
    Eigen::Vector3d vturb_b_mps(vturb_bi_mps, vturb_bii_mps, vturb_biii_mps);
    ang::euler euler_nb(psi_rad, theta_rad, xi_rad);
    ang::dcm R_nb(euler_nb);
    Eigen::Vector3d vtas_n_mps = R_nb * vtas_b_mps;
    Eigen::Vector3d v_n_mps = vtas_n_mps + vwind_n_mps + R_nb * vturb_b_mps;
    double new_chi_rad = ang::euler::obtain_yaw_forward(v_n_mps);
    return ang::tools::angle_diff_rad(new_chi_rad, chi_rad);
}
/* function that returns difference between input body bearing psi and what can be computed from rest of inputs.
 * Difference should be zero. */

double st::sti_Hp_chi::psi_u::exec(const double& psi_rad, const std::vector<double>& par) {
    // par[0]  = chi_rad
    // par[1]  = theta_rad
    // par[2]  = xi_rad
    // par[3]  = vtas_bi_mps
    // par[4]  = vtas_bii_mps
    // par[5]  = vtas_biii_mps
    // par[6]  = vwind_ni_mps
    // par[7]  = vwind_nii_mps
    // par[8]  = vwind_niii_mps
    // par[9]  = vturb_bi_mps
    // par[10] = vturb_bii_mps
    // par[11] = vturb_biii_mps
    return st::sti_Hp_chi::compute_psi_diff(psi_rad, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11]);
}
/* encapsulates minimization function to obtain body bearing psi */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS INITIAL CONDITIONS STI_Hp_psi
// ===================================
// ===================================

st::sti_Hp_psi::sti_Hp_psi(const sti_Hp_psi& Osti)
: sti(*this) {
}
/* copy constructor */

void st::sti_Hp_psi::complete_lambda_phi_Hp(const double& lambda_deg, const double& phi_deg, const double& Hp_m) {
    _lambda_deg = lambda_deg;
    _phi_deg    = phi_deg;
    _Hp_m       = Hp_m;
    _h_m        = std::nan("");

    ang::tools::correct_longitude_deg(_lambda_deg);

    _lambda_rad = _lambda_deg * math::constant::D2R();
    _phi_rad    = _phi_deg * math::constant::D2R();
}
/* add pressure geodetic position (based on pressure altitude) */

void st::sti_Hp_psi::complete_psi_theta_xi(const double& psi_deg, const double& theta_deg, const double& xi_deg) {
    _psi_deg    = psi_deg;
    _theta_deg  = theta_deg;
    _xi_deg     = xi_deg;
    ang::tools::correct_yaw_deg(_psi_deg);

    _psi_rad    = _psi_deg * math::constant::D2R();
    _theta_rad  = _theta_deg * math::constant::D2R();
    _xi_rad     = _xi_deg * math::constant::D2R();

    _chi_deg    = std::nan("");
    _chi_rad    = std::nan("");
}
/* add body yaw, body pitch, body roll */

void st::sti_Hp_psi::complete_and_fill_up_st_truth(st::st_truth& Ost_init, const env::earth& Oearth) {
    Ost_init.get_t_sec() = _t_sec; // time

    Ost_init.get_m_kg() = _m_kg; // mass

    Ost_init.get_delta_control() << _delta_control; // control parameters

    ang::euler euler_nedbfs_rad(_psi_rad, _theta_rad, _xi_rad);
    Ost_init.get_q_nb() = euler_nedbfs_rad; // quaternion

    double DeltaT_degK = Oearth.get_offsets().compute_DeltaT_degK(_t_sec, _lambda_rad, _phi_rad);
    double Deltap_pa = Oearth.get_offsets().compute_Deltap_pa(_t_sec, _lambda_rad, _phi_rad);
    env::atm Oatm(DeltaT_degK, Deltap_pa);
    double H_m = Oatm.Hp2H(_Hp_m);
    double Hp_m_bias = Oatm.H2Hp(H_m, false); // only purpose is to initialize function H2Hp for later use, do not remove
    _h_m = Oearth.get_geo().H2h(H_m, _phi_rad);

    Ost_init.get_x_gdt_b_rad_m().get_lambda_rad() = _lambda_rad; // geodetic coordinates
    Ost_init.get_x_gdt_b_rad_m().get_phi_rad() = _phi_rad;
    Ost_init.get_x_gdt_b_rad_m().get_h_m() = _h_m;

    double N_m = Oearth.get_geo().radius_vert(_phi_rad);
    double M_m = Oearth.get_geo().radius_mer(_phi_rad, N_m);

    ang::euler euler_wb_rad(-_beta_rad, _alpha_rad, 0.);
    Eigen::Vector3d vtas_b_mps = env::speed::vtas2vtasbfs(_vtas_mps, euler_wb_rad);

    Eigen::Vector3d vlf_n_mps = Oearth.get_wind().compute_wind_ned(_t_sec, Ost_init.get_x_gdt_b_rad_m());
    Eigen::Vector3d vlf_b_mps = Ost_init.get_q_nb() / vlf_n_mps;

    double hground_m = 0.; ////////////////////////////////////////////
    double height_m  = _h_m - hground_m;
    Eigen::Vector3d vhf_b_mps = Oearth.get_turb().compute_wind_high_frequency_twomils_bfs(_t_sec, height_m, ang::euler::obtain_yaw_forward(vlf_n_mps), Ost_init.get_q_nb());

    Ost_init.get_v_b_mps() = vtas_b_mps + vlf_b_mps + vhf_b_mps; // speed
    Eigen::Vector3d v_n_mps = Ost_init.get_q_nb() * Ost_init.get_v_b_mps();

    _chi_rad   = ang::euler::obtain_yaw_forward(v_n_mps);
    _chi_deg   = _chi_rad * math::constant::R2D();
    _gamma_rad = ang::euler::obtain_pitch_forward(v_n_mps);
    _gamma_deg = _gamma_rad * math::constant::R2D();

    ang::so3_tangent w_enn_rps = Oearth.get_geo().compute_wenn_rps(v_n_mps, N_m, M_m, _phi_rad, _h_m);
    Ost_init.get_w_ibb_rps() = _w_nbb_rps() + Ost_init.get_q_nb() / (w_enn_rps() + Oearth.get_geo().compute_wien_rps(_phi_rad)());
    Ost_init.get_xi_ebb_mrps().set_vi(Ost_init.get_v_b_mps());
    Ost_init.get_xi_ebb_mrps().set_w(ang::so3_tangent(_w_nbb_rps() + Ost_init.get_q_nb() / w_enn_rps()));

    ang::rodrigues q_en(Ost_init.get_x_gdt_b_rad_m().obtain_euler_en());
    Ost_init.get_z_eb() = ang::dual(q_en * Ost_init.get_q_nb(), Oearth.get_geo().geodetic2cartesian(Ost_init.get_x_gdt_b_rad_m(), N_m)());

    double v_mps   = Ost_init.get_v_b_mps().norm();
    double vlf_mps = vlf_n_mps.norm();
    double vhf_mps = vhf_b_mps.norm();

    std::cout << std::endl << "INITIAL CONDITIONS:" << std::endl << std::endl
              << "t [sec]:          " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _t_sec << std::endl
              << "m [kg]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _m_kg << std::endl
              << "h [m]:            " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _h_m << std::endl
              << "Hp [m]:           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _Hp_m << std::endl
              << "vtas [mps]        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _vtas_mps << std::endl
              << "vlf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vlf_mps << std::endl
              << "vhf [mps]         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vhf_mps << std::endl
              << "v [mps]           " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_mps << std::endl
              << "DeltaT [degK]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << DeltaT_degK << std::endl
              << "psi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _psi_deg << std::endl
              << "theta [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _theta_deg << std::endl
              << "xi [deg]:         " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _xi_deg << std::endl
              << "chi [deg]:        " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _chi_deg << std::endl
              << "gamma [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _gamma_deg << std::endl
              << "alpha [deg]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _alpha_deg << std::endl
              << "beta [deg]:       " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _beta_deg << std::endl
              << "delta_thr [-]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(0) << std::endl
              << "delta_elv [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(1) << std::endl
              << "delta_ail [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(2) << std::endl
              << "delta_rud [deg]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << _delta_control(3) << std::endl
              << "vtas b i [mps]:   " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[0] << std::endl
              << "vtas b ii [mps]:  " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[1] << std::endl
              << "vtas b iii [mps]: " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << vtas_b_mps[2] << std::endl
              << "v b i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[0] << std::endl
              << "v b ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[1] << std::endl
              << "v b iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << Ost_init.get_v_b_mps()[2] << std::endl
              << "v n i [mps]:      " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[0] << std::endl
              << "v n ii [mps]:     " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[1] << std::endl
              << "v n iii [mps]:    " << std::fixed << std::setw(12) << std::setprecision(3) << std::showpos << v_n_mps[2] << std::endl
              << std::endl;
}
/* fill up state vector with initial conditions */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




