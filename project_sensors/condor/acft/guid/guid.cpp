#include "guid.h"
#include "trg.h"
#include "../st/sti.h"
#include "math/logic/seeder.h"
#include "ang/tools.h"
#include <iostream>
#include <iomanip>

// CLASS GUIDANCE_OPERATION
// ========================
// ========================

control::guid_op::guid_op(control::logic::THR_ID guid_thr_id,
                         const double& guid_thr_val,    control::logic::ELV_ID guid_elv_id,
                         const double& guid_elv_val,    control::logic::AIL_ID guid_ail_id,
                         const double& guid_ail_val,    control::logic::RUD_ID guid_rud_id,
                         const double& guid_rud_val,    control::logic::TRG_ID guid_trg_id,
                         const double& guid_trg_val)
: _guid_thr_id(guid_thr_id), _guid_elv_id(guid_elv_id), _guid_ail_id(guid_ail_id), _guid_rud_id(guid_rud_id),
  _guid_val(guid_thr_val, guid_elv_val, guid_ail_val, guid_rud_val),
  _guid_trg_id(guid_trg_id), _guid_trg_val(guid_trg_val),
  _Pguid_trg_functor(nullptr) {
    switch (_guid_trg_id) {
        case control::logic::trgg_t_sec:
            _Pguid_trg_functor = new control::trg_t_sec(*this);
            break;
        case control::logic::trgg_Deltat_sec:
            _Pguid_trg_functor = new control::trg_Deltat_sec(*this);
            break;
        case control::logic::trgg_h_m:
            _Pguid_trg_functor = new control::trg_h_m(*this);
            break;
        case control::logic::trgg_Hp_m:
            _Pguid_trg_functor = new control::trg_Hp_m(*this);
            break;
        case control::logic::trgg_gamma_deg:
            _Pguid_trg_functor = new control::trg_gamma_deg(*this);
            break;
        case control::logic::trgg_chi_deg:
            _Pguid_trg_functor = new control::trg_chi_deg(*this);
            break;
        case control::logic::trgg_psi_deg:
            _Pguid_trg_functor = new control::trg_psi_deg(*this);
            break;
        case control::logic::trgg_id_size:
        default:
            throw std::runtime_error("Trigger case not available");
    }
}
/* constructor based on guidance identificators */

control::guid_op::~guid_op() {
    delete _Pguid_trg_functor;
}
/* destructor */

void control::guid_op::init_trg(const st::st_nav_out& Ost_nav_out) {
    _Pguid_trg_functor->initialize(Ost_nav_out);
}
/* initialize trigger so it later jumps when it switches sign */

bool control::guid_op::eval_trg(const st::st_nav_out& Ost_nav_out) const {
    return _Pguid_trg_functor->eval(Ost_nav_out);
}
/* evaluate trigger, returning true when it jumps */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS GUIDANCE
// ==============
// ==============

control::guid::guid(const unsigned short& nel_op, const double& t_sec_end, const int& seed, const unsigned short& case_guid)
: _sti_id(st::logic::sti_size), _nel_op(nel_op), _Vguid_op(nel_op, nullptr), _t_sec_end(t_sec_end), _case_guid(case_guid),
  _gen(seed), _dist_uniform_bearing(-179, 180), _dist_normal(0.,1.), _dist_uniform_time(10,50) {
}
/* constructor based on number of operations (user responsibility to fill it up afterwards) and final time.
 * Seed not employed*/

control::guid::~guid() {
    for (unsigned short i = 0; i != _nel_op; ++i) {
        delete _Vguid_op[i];
    }
}
/* destructor */

void control::guid::add_op(const unsigned short& index, control::guid_op* Pop) {
    delete _Vguid_op[index];
    _Vguid_op[index] = Pop;
}
/* adds operation at position indicated by index */

using namespace control::logic;

control::guid* control::guid::create_guid(math::seeder& Oseeder, const unsigned short& case_guid, env::logic::ZONE_ID zone_id, const double& t_sec_gpsloss, bool flag_console) {
    switch (case_guid) {

        // ==================================================================
        // ===== ===== ===== Random Navigation Trajectories ===== ===== =====
        // ==================================================================

        case 1: { // Benchmark Scenario --> change of bearing, airspeed, and pressure altitude


            double t_sec_end = 3800.0;
            auto Pguid = new control::guid(6, t_sec_end, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;

            double beta_deg = 0.0;

            double chi_deg_ini = Pguid->_dist_uniform_bearing(Pguid->_gen);
            double chi_deg_end_eff, chi_deg_end_trg, angle_diff_deg, xi_deg;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_end_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_end_eff, chi_deg_ini);
            } while (std::fabs(angle_diff_deg) < 10.0);
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg          = +10.0;
                chi_deg_end_trg = chi_deg_end_eff - 5.;
            }
            else {
                xi_deg          = -10.0;
                chi_deg_end_trg = chi_deg_end_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_end_trg);

            double Deltat_sec_tas;
            do {
                Deltat_sec_tas = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas < 150.0);

            double vtas_mps_ini;
            do {
                vtas_mps_ini = 29.0 + 1.5 * Pguid->_dist_normal(Pguid->_gen);
            } while ((vtas_mps_ini > 34.0) || (vtas_mps_ini < 24.0));
            double vtas_mps_end, vtas_mps_diff;
            do {
                vtas_mps_diff = 1.5 * Pguid->_dist_normal(Pguid->_gen);
                vtas_mps_end  = vtas_mps_ini + vtas_mps_diff;

            } while ((std::fabs(vtas_mps_diff) < 0.5) || (vtas_mps_end > 34.0) || (vtas_mps_end < 24.0));

            double Deltat_sec_Hp;
            do {
                Deltat_sec_Hp = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while ((Deltat_sec_Hp < 150.0) || ((Deltat_sec_tas + Deltat_sec_Hp) > 2500.0));
            double Hp_m_ini = 2700.0 + 200.0 * Pguid->_dist_normal(Pguid->_gen);
            double Hp_m_end, Hp_m_diff, gammaTAS_deg;
            do { // ensure that there is a pressure altitude change of at least 100 [m]
                Hp_m_diff = 300.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (std::fabs(Hp_m_diff) < 100.0);
            Hp_m_end = Hp_m_ini + Hp_m_diff;
            // path angle depends on initial and final pressure altitudes
            if (Hp_m_diff > 0) {gammaTAS_deg = +2.0;}
            else               {gammaTAS_deg = -2.0;}

            double Xt_sec_gpsloss = 100.0;
            double t_sec_turn;
            do {
                t_sec_turn = Xt_sec_gpsloss + 30.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (t_sec_turn < (Xt_sec_gpsloss + 15.0));

            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_ini,     rud_beta_deg, beta_deg, trgg_t_sec,        t_sec_turn));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_xi_deg,  xi_deg,          rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_end_trg));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_tas));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_Hp));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_gammaTAS_deg, gammaTAS_deg, ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Hp_m,       Hp_m_end));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_end,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_end));
            if (flag_console == true) {
                Pguid->create_text(std::cout);
            }

            return Pguid;



//            double t_sec_end = 3800.0;
//            auto Pguid = new control::guid(6, t_sec_end, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
//            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
//
//            double beta_deg = 0.0;
//
//            // Turn magnitude
//            double chi_deg_ini = Pguid->_dist_uniform_bearing(Pguid->_gen);
//            double chi_deg_end_eff, chi_deg_end_trg, angle_diff_deg, xi_deg;
//            do { // ensure that there is a direction change of at least 10 [deg]
//                chi_deg_end_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
//                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_end_eff, chi_deg_ini);
//            } while (std::fabs(angle_diff_deg) < 10.0);
//            // turn direction depends on initial and final bearings
//            // turning end 5 [deg] earlier for better control
//            if (angle_diff_deg > 0) {
//                xi_deg          = +10.0;
//                chi_deg_end_trg = chi_deg_end_eff - 5.;
//            }
//            else {
//                xi_deg          = -10.0;
//                chi_deg_end_trg = chi_deg_end_eff + 5.;
//            }
//            ang::tools::correct_yaw_deg(chi_deg_end_trg);
//
//            // Airspeed change timing and magnitude
//            double Deltat_sec_tas;
//            do {
//                Deltat_sec_tas = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
//            } while (Deltat_sec_tas < 150.0);
//
//            double vtas_mps_ini;
//            do {
//                vtas_mps_ini = 29.0 + 1.5 * Pguid->_dist_normal(Pguid->_gen);
//            } while ((vtas_mps_ini > 34.0) || (vtas_mps_ini < 24.0));
//            double vtas_mps_end, vtas_mps_diff;
//            do {
//                vtas_mps_diff = 1.5 * Pguid->_dist_normal(Pguid->_gen);
//                vtas_mps_end  = vtas_mps_ini + vtas_mps_diff;
//
//            } while ((std::fabs(vtas_mps_diff) < 0.5) || (vtas_mps_end > 34.0) || (vtas_mps_end < 24.0));
//
//            // Altitude change timing and magnitude
//            double Deltat_sec_Hp;
//            do {
//                Deltat_sec_Hp = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
//            } while ((Deltat_sec_Hp < 150.0) || ((Deltat_sec_tas + Deltat_sec_Hp) > 2500.0));
//            double Hp_m_ini = 2700.0 + 200.0 * Pguid->_dist_normal(Pguid->_gen);
//            double Hp_m_end, Hp_m_diff, gammaTAS_deg;
//            do { // ensure that there is a pressure altitude change of at least 100 [m]
//                Hp_m_diff = 300.0 * Pguid->_dist_normal(Pguid->_gen);
//            } while (std::fabs(Hp_m_diff) < 100.0);
//            Hp_m_end = Hp_m_ini + Hp_m_diff;
//            // path angle depends on initial and final pressure altitudes
//            if (Hp_m_diff > 0) {gammaTAS_deg = +2.0;}
//            else               {gammaTAS_deg = -2.0;}
//
//            // Turn timing
//            double t_sec_turn;
//            do {
//                t_sec_turn = t_sec_gpsloss + 30.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen);
//            } while (t_sec_turn < (t_sec_gpsloss + 15.0));
//
//            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_ini,     rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_turn));
//            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_xi_deg,  xi_deg,          rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_end_trg));
//            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_tas));
//            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_Hp));
//            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_gammaTAS_deg, gammaTAS_deg, ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Hp_m,       Hp_m_end));
//            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_end,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_end));
//            if (flag_console == true) {
//                Pguid->create_text(std::cout);
//            }
//
//            return Pguid;




        }

        case 2: { // same as 1 but shorter (1000 sec instead of 3800 sec). Employed for VNSE Sensitivity Analysis
            double t_sec_end = 1000.0; // TODO original #1 has 3800.0
            auto Pguid = new control::guid(6, t_sec_end, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;

            double beta_deg = 0.0;
            double chi_deg_ini = Pguid->_dist_uniform_bearing(Pguid->_gen);
            double chi_deg_end_eff, chi_deg_end_trg, angle_diff_deg, xi_deg;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_end_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);  // TODO this is the original of #1
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_end_eff, chi_deg_ini);
            } while (std::fabs(angle_diff_deg) < 10.0);
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                xi_deg          = +10.0; // TODO this is the original of #1
                //xi_deg         = +5.0; // TODO XI_05
                //xi_deg         = +15.0; // TODO XI_15
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                chi_deg_end_trg = chi_deg_end_eff - 5.;
            }
            else {
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                xi_deg          = -10.0; // TODO this is the original of #1
                //xi_deg          = -5.0; // TODO XI_05
                //xi_deg          = -15.0; // TODO XI_15
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                chi_deg_end_trg = chi_deg_end_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_end_trg);

            double Deltat_sec_tas;
            do {
                Deltat_sec_tas = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas < 150.0);

            double vtas_mps_ini;
            do {
                vtas_mps_ini = 29.0 + 1.5 * Pguid->_dist_normal(Pguid->_gen);
            } while ((vtas_mps_ini > 34.0) || (vtas_mps_ini < 24.0));
            double vtas_mps_end, vtas_mps_diff;
            do {
                vtas_mps_diff = 1.5 * Pguid->_dist_normal(Pguid->_gen);
                vtas_mps_end  = vtas_mps_ini + vtas_mps_diff;

            } while ((std::fabs(vtas_mps_diff) < 0.5) || (vtas_mps_end > 34.0) || (vtas_mps_end < 24.0));

            double Deltat_sec_Hp;
            do {
                Deltat_sec_Hp = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while ((Deltat_sec_Hp < 150.0) || ((Deltat_sec_tas + Deltat_sec_Hp) > 2500.0));
            double Hp_m_end, Hp_m_diff, gammaTAS_deg;
            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////
            double Hp_m_ini = 2700.0 + 200.0 * Pguid->_dist_normal(Pguid->_gen); // TODO this is the original of #1
            //double Hp_m_ini = 2200.0 + 150.0 * Pguid->_dist_normal(Pguid->_gen); // TODO H_2200_150
            //double Hp_m_ini = 1700.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen); // TODO H_1700_100
            //double Hp_m_ini = 1200.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen); // TODO H_1200_050
            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////
            do { // ensure that there is a pressure altitude change of at least 100 [m]
                Hp_m_diff = 300.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (std::fabs(Hp_m_diff) < 100.0);
            Hp_m_end = Hp_m_ini + Hp_m_diff;
            // path angle depends on initial and final pressure altitudes
            if (Hp_m_diff > 0) {gammaTAS_deg = +2.0;}
            else               {gammaTAS_deg = -2.0;}

            double t_sec_turn;
            do {
                t_sec_turn = t_sec_gpsloss + 30.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (t_sec_turn < (t_sec_gpsloss + 15.0));

            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_ini,     rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_turn));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_xi_deg,  xi_deg,          rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_end_trg));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_tas));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_ini,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_Hp));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_gammaTAS_deg, gammaTAS_deg, ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_Hp_m,       Hp_m_end));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, vtas_mps_end, elv_Hp_m,         Hp_m_end,     ail_chi_deg, chi_deg_end_eff, rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_end));
            if (flag_console == true) {
                Pguid->create_text(std::cout);
            }

            return Pguid;
        }

        case 3: { // Alternate Scenario --> eight bearing changes at constant airspeed and pressure altitude
            // dummy variables are necessary to maintain commonality with intent #1

            double t_sec_end = 500.0;
            auto Pguid = new control::guid(17, t_sec_end, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;

            double beta_deg = 0.0;
            double angle_diff_deg;

            double chi_deg_0 = Pguid->_dist_uniform_bearing(Pguid->_gen);

            double chi_deg_1_eff, chi_deg_1_trg, xi_deg_1;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_1_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_1_eff, chi_deg_0);
            } while (std::fabs(angle_diff_deg) < 10.0);
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg_1      = +10.0;
                chi_deg_1_trg = chi_deg_1_eff - 5.;
            }
            else {
                xi_deg_1      = -10.0;
                chi_deg_1_trg = chi_deg_1_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_1_trg);

            double Deltat_sec_tas_dummy;
            do {
                Deltat_sec_tas_dummy = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas_dummy < 150.0);

            double vtas_mps_ini;
            do {
                vtas_mps_ini = 29.0 + 1.5 * Pguid->_dist_normal(Pguid->_gen);
            } while ((vtas_mps_ini > 34.0) || (vtas_mps_ini < 24.0));
            double vtas_mps_end_dummy, vtas_mps_diff_dummy;
            do {
                vtas_mps_diff_dummy = 1.5 * Pguid->_dist_normal(Pguid->_gen);
                vtas_mps_end_dummy  = vtas_mps_ini + vtas_mps_diff_dummy;

            } while ((std::fabs(vtas_mps_diff_dummy) < 0.5) || (vtas_mps_end_dummy > 34.0) || (vtas_mps_end_dummy < 24.0));

            double Deltat_sec_Hp_dummy;
            do {
                Deltat_sec_Hp_dummy = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while ((Deltat_sec_Hp_dummy < 150.0) || ((Deltat_sec_tas_dummy + Deltat_sec_Hp_dummy) > 2500.0));

            double Hp_m_ini = 2700.0 + 200.0 * Pguid->_dist_normal(Pguid->_gen);
            double Hp_m_end_dummy, Hp_m_diff_dummy;
            do { // ensure that there is a pressure altitude change of at least 100 [m]
                Hp_m_diff_dummy = 300.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (std::fabs(Hp_m_diff_dummy) < 100.0);
            Hp_m_end_dummy = Hp_m_ini + Hp_m_diff_dummy;

            double t_sec_turn_1;
            do {
                t_sec_turn_1 = t_sec_gpsloss + 30.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (t_sec_turn_1 < (t_sec_gpsloss + 15.0));

            // up this point I have the same calls to distributions as guidance intent #1
            double Deltat_sec_turn_2 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_3 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_4 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_5 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_6 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_7 = Pguid->_dist_uniform_time(Pguid->_gen);
            double Deltat_sec_turn_8 = Pguid->_dist_uniform_time(Pguid->_gen);

            double chi_deg_2_eff, chi_deg_2_trg, xi_deg_2;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_2_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_2_eff, chi_deg_1_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_2      = +10.0;
                chi_deg_2_trg = chi_deg_2_eff - 5.;
            }
            else {
                xi_deg_2      = -10.0;
                chi_deg_2_trg = chi_deg_2_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_2_trg);

            double chi_deg_3_eff, chi_deg_3_trg, xi_deg_3;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_3_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_3_eff, chi_deg_2_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_3      = +10.0;
                chi_deg_3_trg = chi_deg_3_eff - 5.;
            }
            else {
                xi_deg_3      = -10.0;
                chi_deg_3_trg = chi_deg_3_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_3_trg);

            double chi_deg_4_eff, chi_deg_4_trg, xi_deg_4;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_4_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_4_eff, chi_deg_3_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_4      = +10.0;
                chi_deg_4_trg = chi_deg_4_eff - 5.;
            }
            else {
                xi_deg_4      = -10.0;
                chi_deg_4_trg = chi_deg_4_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_4_trg);

            double chi_deg_5_eff, chi_deg_5_trg, xi_deg_5;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_5_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_5_eff, chi_deg_4_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_5      = +10.0;
                chi_deg_5_trg = chi_deg_5_eff - 5.;
            }
            else {
                xi_deg_5      = -10.0;
                chi_deg_5_trg = chi_deg_5_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_5_trg);

            double chi_deg_6_eff, chi_deg_6_trg, xi_deg_6;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_6_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_6_eff, chi_deg_5_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_6      = +10.0;
                chi_deg_6_trg = chi_deg_6_eff - 5.;
            }
            else {
                xi_deg_6      = -10.0;
                chi_deg_6_trg = chi_deg_6_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_6_trg);

            double chi_deg_7_eff, chi_deg_7_trg, xi_deg_7;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_7_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_7_eff, chi_deg_6_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_7      = +10.0;
                chi_deg_7_trg = chi_deg_7_eff - 5.;
            }
            else {
                xi_deg_7      = -10.0;
                chi_deg_7_trg = chi_deg_7_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_7_trg);

            double chi_deg_8_eff, chi_deg_8_trg, xi_deg_8;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_8_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_8_eff, chi_deg_7_eff);
            } while (std::fabs(angle_diff_deg) < 10.0);
            if (angle_diff_deg > 0) {
                xi_deg_8      = +10.0;
                chi_deg_8_trg = chi_deg_8_eff - 5.;
            }
            else {
                xi_deg_8      = -10.0;
                chi_deg_8_trg = chi_deg_8_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_8_trg);

            Pguid->add_op(0,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_0,     rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_turn_1));
            Pguid->add_op(1,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_1,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_1_trg));
            Pguid->add_op(2,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_1_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_2));
            Pguid->add_op(3,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_2,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_2_trg));
            Pguid->add_op(4,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_2_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_3));
            Pguid->add_op(5,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_3,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_3_trg));
            Pguid->add_op(6,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_3_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_4));
            Pguid->add_op(7,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_4,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_4_trg));
            Pguid->add_op(8,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_4_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_5));
            Pguid->add_op(9,  new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_5,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_5_trg));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_5_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_6));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_6,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_6_trg));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_6_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_7));
            Pguid->add_op(13, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_7,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_7_trg));
            Pguid->add_op(14, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_7_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_8));
            Pguid->add_op(15, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_xi_deg,  xi_deg_8,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_8_trg));
            Pguid->add_op(16, new control::guid_op(thr_vtas_mps, vtas_mps_ini, elv_Hp_m, Hp_m_ini, ail_chi_deg, chi_deg_8_eff, rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_end));

            if (flag_console == true) {
                Pguid->create_text(std::cout);
            }

            return Pguid;
        }

        case 4: { // Intermediate Scenario -> Turn + speed change + turn + altitude change + turn + speed change + turn + altitude change
            // dummy variables are necessary to maintain commonality with intent #1

            double t_sec_end = 2000.0;
            auto Pguid = new control::guid(15, t_sec_end, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;

            double beta_deg = 0.0;
            double angle_diff_deg;

            // 1st turn magnitude
            double chi_deg_0 = Pguid->_dist_uniform_bearing(Pguid->_gen);
            double chi_deg_1_eff, chi_deg_1_trg, xi_deg_1;
            do { // ensure that there is a direction change of at least 10 [deg]
                chi_deg_1_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_1_eff, chi_deg_0);
            } while (std::fabs(angle_diff_deg) < 10.0);
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg_1      = +10.0;
                chi_deg_1_trg = chi_deg_1_eff - 5.;
            }
            else {
                xi_deg_1      = -10.0;
                chi_deg_1_trg = chi_deg_1_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_1_trg);

            // 1st airspeed change dummy timing and magnitude
            double Deltat_sec_tas_dummy;
            do {
                Deltat_sec_tas_dummy = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas_dummy < 150.0);

            double vtas_mps_0;
            do {
                vtas_mps_0 = 29.0 + 1.5 * Pguid->_dist_normal(Pguid->_gen);
            } while ((vtas_mps_0 > 34.0) || (vtas_mps_0 < 24.0));
            double vtas_mps_1, vtas_mps_1_diff;
            do {
                vtas_mps_1_diff = 1.5 * Pguid->_dist_normal(Pguid->_gen);
                vtas_mps_1  = vtas_mps_0 + vtas_mps_1_diff;
            } while ((std::fabs(vtas_mps_1_diff) < 0.5) || (vtas_mps_1 > 34.0) || (vtas_mps_1 < 24.0));

            // 1st altitude change dummy timing and magnitude
            double Deltat_sec_Hp_dummy;
            do {
                Deltat_sec_Hp_dummy = 500.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while ((Deltat_sec_Hp_dummy < 150.0) || ((Deltat_sec_tas_dummy + Deltat_sec_Hp_dummy) > 2500.0));
            double Hp_m_0 = 2700.0 + 200.0 * Pguid->_dist_normal(Pguid->_gen);
            double Hp_m_1, Hp_m_1_diff, gammaTAS_1_deg;
            do { // ensure that there is a pressure altitude change of at least 100 [m]
                Hp_m_1_diff = 300.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (std::fabs(Hp_m_1_diff) < 100.0);
            Hp_m_1 = Hp_m_0 + Hp_m_1_diff;
            // path angle depends on initial and final pressure altitudes
            if (Hp_m_1_diff > 0) {gammaTAS_1_deg = +2.0;}
            else                 {gammaTAS_1_deg = -2.0;}

            // Timings for all maneuvers
            double t_sec_turn_1, Deltat_sec_tas_1, Deltat_sec_turn_2, Deltat_sec_Hp_1, Deltat_sec_turn_3, Deltat_sec_tas_2, Deltat_sec_turn_4, Deltat_sec_Hp_2;
            do {
                t_sec_turn_1 = t_sec_gpsloss + 30.0 + 50.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (t_sec_turn_1 < (t_sec_gpsloss + 15.0));
            do {
                Deltat_sec_tas_1 = 200.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas_1 < 50.0);
            do {
                Deltat_sec_turn_2 = 200.0 + 10.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_turn_2 < 50.0);
            do {
                Deltat_sec_Hp_1 = 200.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_Hp_1 < 50.0);
            do {
                Deltat_sec_turn_3 = 200.0 + 10.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_turn_3 < 50.0);
            do {
                Deltat_sec_tas_2 = 200.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_tas_2 < 50.0);
            do {
                Deltat_sec_turn_4 = 200.0 + 10.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_turn_4 < 50.0);
            do {
                Deltat_sec_Hp_2 = 200.0 + 100.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (Deltat_sec_Hp_2 < 50.0);

            // 2nd turn magnitude
            double chi_deg_2_eff, chi_deg_2_trg, xi_deg_2;
            do { // ensure that there is a direction change of at least 60 [deg] but no more than 120 [deg]
                chi_deg_2_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_2_eff, chi_deg_1_eff);
            } while ((std::fabs(angle_diff_deg) < 60.0) || (std::fabs(angle_diff_deg) > 120.0));
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg_2      = +10.0;
                chi_deg_2_trg = chi_deg_2_eff - 5.;
            }
            else {
                xi_deg_2      = -10.0;
                chi_deg_2_trg = chi_deg_2_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_2_trg);

            // 3rd turn magnitude
            double chi_deg_3_eff, chi_deg_3_trg, xi_deg_3;
            do { // ensure that there is a direction change of at least 60 [deg] but no more than 120 [deg]
                chi_deg_3_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_3_eff, chi_deg_2_eff);
            } while ((std::fabs(angle_diff_deg) < 60.0) || (std::fabs(angle_diff_deg) > 120.0));
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg_3      = +10.0;
                chi_deg_3_trg = chi_deg_3_eff - 5.;
            }
            else {
                xi_deg_3      = -10.0;
                chi_deg_3_trg = chi_deg_3_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_3_trg);

            // 4thd turn magnitude
            double chi_deg_4_eff, chi_deg_4_trg, xi_deg_4;
            do { // ensure that there is a direction change of at least 60 [deg] but no more than 120 [deg]
                chi_deg_4_eff = Pguid->_dist_uniform_bearing(Pguid->_gen);
                angle_diff_deg = ang::tools::angle_diff_deg(chi_deg_4_eff, chi_deg_3_eff);
            } while ((std::fabs(angle_diff_deg) < 60.0) || (std::fabs(angle_diff_deg) > 120.0));
            // turn direction depends on initial and final bearings
            // turning end 5 [deg] earlier for better control
            if (angle_diff_deg > 0) {
                xi_deg_4      = +10.0;
                chi_deg_4_trg = chi_deg_4_eff - 5.;
            }
            else {
                xi_deg_4      = -10.0;
                chi_deg_4_trg = chi_deg_4_eff + 5.;
            }
            ang::tools::correct_yaw_deg(chi_deg_3_trg);

            // 2nd airspeed change magnitude
            double vtas_mps_2, vtas_mps_2_diff;
            do {
                vtas_mps_2_diff = 1.5 * Pguid->_dist_normal(Pguid->_gen);
                vtas_mps_2  = vtas_mps_1 + vtas_mps_2_diff;
            } while ((std::fabs(vtas_mps_2_diff) < 0.5) || (vtas_mps_2 > 34.0) || (vtas_mps_2 < 24.0));

            // 2nd altitude change and magnitude
            double Hp_m_2, Hp_m_2_diff, gammaTAS_2_deg;
            do { // ensure that there is a pressure altitude change of at least 100 [m]
                Hp_m_2_diff = 300.0 * Pguid->_dist_normal(Pguid->_gen);
            } while (std::fabs(Hp_m_2_diff) < 100.0);
            Hp_m_2 = Hp_m_1 + Hp_m_2_diff;
            // path angle depends on initial and final pressure altitudes
            if (Hp_m_2_diff > 0) {gammaTAS_2_deg = +2.0;}
            else                 {gammaTAS_2_deg = -2.0;}

            Pguid->add_op( 0, new control::guid_op(thr_vtas_mps, vtas_mps_0, elv_Hp_m,         Hp_m_0,         ail_chi_deg, chi_deg_0,     rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_turn_1));
            Pguid->add_op( 1, new control::guid_op(thr_vtas_mps, vtas_mps_0, elv_Hp_m,         Hp_m_0,         ail_xi_deg,  xi_deg_1,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_1_trg));
            Pguid->add_op( 2, new control::guid_op(thr_vtas_mps, vtas_mps_0, elv_Hp_m,         Hp_m_0,         ail_chi_deg, chi_deg_1_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_tas_1));
            Pguid->add_op( 3, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_0,         ail_chi_deg, chi_deg_1_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_2));
            Pguid->add_op( 4, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_0,         ail_xi_deg,  xi_deg_2,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_2_trg));
            Pguid->add_op( 5, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_0,         ail_chi_deg, chi_deg_2_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_Hp_1));
            Pguid->add_op( 6, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_gammaTAS_deg, gammaTAS_1_deg, ail_chi_deg, chi_deg_2_eff, rud_beta_deg, beta_deg, trgg_Hp_m,       Hp_m_1));
            Pguid->add_op( 7, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_1,         ail_chi_deg, chi_deg_2_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_3));
            Pguid->add_op( 8, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_1,         ail_xi_deg,  xi_deg_3,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_3_trg));
            Pguid->add_op( 9, new control::guid_op(thr_vtas_mps, vtas_mps_1, elv_Hp_m,         Hp_m_1,         ail_chi_deg, chi_deg_3_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_tas_2));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, vtas_mps_2, elv_Hp_m,         Hp_m_1,         ail_chi_deg, chi_deg_3_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_turn_4));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, vtas_mps_2, elv_Hp_m,         Hp_m_1,         ail_xi_deg,  xi_deg_4,      rud_beta_deg, beta_deg, trgg_chi_deg,    chi_deg_4_trg));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, vtas_mps_2, elv_Hp_m,         Hp_m_1,         ail_chi_deg, chi_deg_4_eff, rud_beta_deg, beta_deg, trgg_Deltat_sec, Deltat_sec_Hp_2));
            Pguid->add_op(13, new control::guid_op(thr_vtas_mps, vtas_mps_2, elv_gammaTAS_deg, gammaTAS_2_deg, ail_chi_deg, chi_deg_4_eff, rud_beta_deg, beta_deg, trgg_Hp_m,       Hp_m_2));
            Pguid->add_op(14, new control::guid_op(thr_vtas_mps, vtas_mps_2, elv_Hp_m,         Hp_m_2,         ail_chi_deg, chi_deg_4_eff, rud_beta_deg, beta_deg, trgg_t_sec,      t_sec_end));

            if (flag_console == true) {
                Pguid->create_text(std::cout);
            }

            return Pguid;
        }

        // ========================================================
        // ===== ===== ===== Control fine tuning ===== ===== =====
        // ========================================================

        case 51: { // turn initiation based on absolute speed bearing
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg,  30.0, rud_beta_deg, 0., trgg_t_sec,   50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,   10.0, rud_beta_deg, 0., trgg_t_sec,  150.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,  -10.0, rud_beta_deg, 0., trgg_t_sec,  250.0));
            return Pguid;
        }
        case 52: { // turn initiation based on body heading
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_psiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec, 50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_t_sec, 150.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_t_sec, 250.0));
            return Pguid;
        }
        case 53: { // speed changes
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,   50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 34.0, elv_Hp_m, 3000.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,  150.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 26.0, elv_Hp_m, 3000.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,  250.0));
            return Pguid;
        }
        case 54: { // body pitch changes
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,      3000.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,   50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,   +4.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,  150.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,   -4.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,  250.0));
            return Pguid;
        }
        case 55: { // turn conclusion based on absolute speed bearing (note trigger 5 deg less than next op)
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg,  30.0, rud_beta_deg, 0., trgg_t_sec,    50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,   10.0, rud_beta_deg, 0., trgg_chi_deg, 115.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, 120.0, rud_beta_deg, 0., trgg_t_sec,   250.0));
            return Pguid;
        }
        case 56: { // turn conclusion based on body heading (note trigger 5 deg less than next op)
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_psiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg,  30.0, rud_beta_deg, 0., trgg_t_sec,    50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,   10.0, rud_beta_deg, 0., trgg_psi_deg, 115.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg, 120.0, rud_beta_deg, 0., trgg_t_sec,   250.0));
            return Pguid;
        }
        case 57: { // altitude changes changes (pressure altitude)
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,         3000.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,   50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg,   +2.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_Hp_m,  3100.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,         3100.0, ail_chi_deg, 30.0, rud_beta_deg, 0., trgg_t_sec,  250.0));
            return Pguid;
        }
        case 58: {
            auto Pguid = new control::guid(3, 250.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 29.0, elv_Hp_m, 2700.0, ail_chi_deg,  30.0, rud_beta_deg, 0., trgg_t_sec,    50.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 29.0, elv_Hp_m, 2700.0, ail_xi_deg,   20.0, rud_beta_deg, 0., trgg_chi_deg, 215.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 29.0, elv_Hp_m, 2700.0, ail_chi_deg, 220.0, rud_beta_deg, 0., trgg_t_sec,   250.0));
            return Pguid;
        }

        // ========================================================
        // ===== ===== ===== PNTT Trajectories ===== ===== =====
        // ========================================================

        case 99: {
            auto Pguid = new control::guid(1, 5000.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_HpXX_tasXX_chiXX;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 29.0, elv_Hp_m, 2700.0, ail_chi_deg,  30.0, rud_beta_deg, 0., trgg_t_sec,    5000.0));
            return Pguid;
        }

        // ========================================================
        // ===== ===== ===== Control Trajectories ===== ===== =====
        // ========================================================

        case 151: { // primary loops - accelerate and decelerate in vtas +- 2 [mps] steps - starts at 3000 - 25
            // Calm: max errors are vtas [mps] +2.003, theta [deg] -0.591, xi [deg] +0.652, beta [deg] -0.208 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -2.4883, theta [deg] -1.226, xi [deg] -0.838, beta [deg] +1.091
            // NONE  mean-std error vtas [mps] +1.43e-3 +1.19e-1, alpha [deg] +2.18e-3 +1.15e-1, beta [deg] -1.63e-3 +1.32e-1, Hp [m] +8.59e-1 +4.63e0
            // XACT  mean-std error rtov [deg] +5.82e-1 +3.39e-1, w_nbb [dps] +6.61e-2 +3.08e-2, bias_gyr [dps] +9.61e-3 +1.41e-2, bias_mag [nT] +2.42e+1 +4.50e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -2.472, theta [deg] -1.275, xi [deg] -0.781, beta [deg] +1.2351
            // mean-std error rtov [deg] +5.35e-1 +3.87e-1, w_nbb [dps] +6.59e-2 +2.84e-2, bias_gyr [dps] +8.63e-3 +1.07e-2, bias_mag [nT] +2.44e+1 +1.90e+0
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas25_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 25.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 27.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 29.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 31.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 33.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 31.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 29.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 27.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 25.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 152: { // primary loops - accelerate and decelerate in vtas +- 4 [mps] steps - starts at 3000 - 25
            // Calm: max errors are vtas [mps] +4.004, theta [deg] -1.058, xi [deg] +1.264, beta [deg] +0.398 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -4.279, theta [deg] -1.637, xi [deg] +1.179, beta [deg] +1.229
            // NONE  mean-std error vtas [mps] +9.25e-4 +1.42e-1, alpha [deg] +4.99e-3 +1.18e-1, beta [deg] -2.29e-3 +1.33e-1, Hp [m] +9.40e-1 +4.98e0
            // XACT  mean-std error rtov [deg] +6.05e-1 +3.78e-1, w_nbb [dps] +6.80e-2 +3.30e-2, bias_gyr [dps] +1.38e-2 +1.69e-2, bias_mag [nT] +2.42e+1 +5.60e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -4.533, theta [deg] -1.709, xi [deg] +1.253, beta [deg] +1.267
            // mean-std error rtov [deg] +6.17e-1 +3.61e-1, w_nbb [dps] +6.73e-2 +2.90e-2, bias_gyr [dps] +1.24e-2 +1.26e-2, bias_mag [nT] +2.38e+1 +1.98e+0
            auto Pguid = new control::guid(5, 160.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas25_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 25.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 29.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 33.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 29.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 25.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 153: { // primary loops - pull up and push down in theta +- 2 [deg] steps - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -0.422, theta [deg] +2.002, xi [deg] -0.104, beta [deg] -0.109 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -0.972, theta [deg] -2.446, xi [deg] -0.546, beta [deg] -0.973 (seed 1)
            // NONE  mean-std error vtas [mps] +1.52e-3 +1.04e-1, alpha [deg] +2.15e-3 +1.15e-1, beta [deg] -1.65e-3 +1.29e-1, Hp [m] +2.32e-2 +5.63e0
            // XACT  mean-std error rtov [deg] +7.69e-1 +3.55e-1, w_nbb [dps] +6.79e-2 +3.50e-2, bias_gyr [dps] +1.05e-2 +1.17e-2, bias_mag [nT] +1.10e+1 +2.05e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -0.990, theta [deg] -2.445, xi [deg] -0.583, beta [deg] -1.018 (seed 1)
            // mean-std error rtov [deg] +5.79e-1 +3.31e-1, w_nbb [dps] +6.81e-2 +3.54e-2, bias_gyr [dps] +8.75e-3 +1.11e-2, bias_mag [nT] +1.13e+1 +1.30e+0
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -4.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 154: { // primary loops - pull up and push down in theta +- 4 [deg] steps - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -0.838, theta [deg] -4.001, xi [deg] -0.208, beta [deg] -0.109 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -1.257, theta [deg] -4.308, xi [deg] -0.588, beta [deg] -1.038 (seed 1)
            // NONE  mean-std error vtas [mps] +1.12e-3 +1.09e-1, alpha [deg] +4.96e-3 +1.22e-1, beta [deg] -2.30e-3 +1.29e-1, Hp [m] +3.98e-1 +5.91e0
            // XACT  mean-std error rtov [deg] +8.43e-1 +3.83e-1, w_nbb [dps] +7.16e-2 +4.77e-2, bias_gyr [dps] +1.55e-2 +1.31e-2, bias_mag [nT] +1.16e+1 +2.47e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -1.299, theta [deg] -4.318, xi [deg] -0.639, beta [deg] -1.082 (seed 1)
            // mean-std error rtov [deg] +6.79e-1 +3.51e-1, w_nbb [dps] +7.18e-2 +4.92e-2, bias_gyr [dps] +1.31e-2 +1.30e-2, bias_mag [nT] +1.16e+1 +9.41e-1
            auto Pguid = new control::guid(5, 160.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -4.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 155: { // primary loops - roll in and roll out in xi +- 5 [deg] steps - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -0.054, theta [deg] +0.220, xi [deg] +5.003, beta [deg] +1.262 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -0.876, theta [deg] +0.839, xi [deg] +5.170, beta [deg] -2.093 (seed 1)
            // NONE  mean-std error vtas [mps] +1.99e-3 +1.05e-1, alpha [deg] +1.62e-3 +1.13e-1, beta [deg] -2.20e-3 +1.32e-1, Hp [m] -1.92e-1 +4.42e0
            // XACT  mean-std error rtov [deg] +3.70e-1 +3.63e-1, w_nbb [dps] +6.75e-2 +3.52e-2, bias_gyr [dps] +6.88e-3 +1.10e-2, bias_mag [nT] +5.55e+0 +3.50e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -0.931, theta [deg] +0.847, xi [deg] +5.190, beta [deg] -2.184 (seed 1)
            // mean-std error rtov [deg] +3.41e-1 +3.19e-1, w_nbb [dps] +6.79e-2 +3.55e-2, bias_gyr [dps] +6.13e-3 +1.01e-2, bias_mag [nT] +5.02e+0 +3.18e+0
            auto Pguid = new control::guid(13, 400.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  15.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(9,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -15.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 156: { // primary loops - roll in and roll out in xi +- 10 [deg] steps - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -0.193, theta [deg] +0.938, xi [deg] +10.010, beta [deg] +2.517 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] +0.884, theta [deg] +1.440, xi [deg] +10.161, beta [deg] -3.142 (seed 1)
            // NONE  mean-std error vtas [mps] +1.99e-3 +1.05e-1, alpha [deg] +1.62e-3 +1.14e-1, beta [deg] -2.20e-3 +1.37e-1, Hp [m] -1.95e-1 +4.62e0
            // XACT  mean-std error rtov [deg] +3.84e-1 +3.84e-1, w_nbb [dps] +7.09e-2 +5.42e-2, bias_gyr [dps] +6.78e-3 +1.11e-2, bias_mag [nT] +5.34e+0 +3.95e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -0.914, theta [deg] +1.439, xi [deg] +10.190, beta [deg] -3.267 (seed 1)
            // mean-std error rtov [deg] +3.73e-1 +3.39e-1, w_nbb [dps] +7.14e-2 +5.50e-2, bias_gyr [dps] +6.09e-3 +1.03e-2, bias_mag [nT] +4.32e+0 +3.31e+0
            auto Pguid = new control::guid(13, 400.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  20.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  30.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  20.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -20.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(9,  new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -30.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -20.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 157: { // primary loops - simultaneous acceleration - deceleration and push-pull maneuvers - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -2.003, theta [deg] -2.001, xi [deg] +0.635, beta [deg] -0.210 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -2.481, theta [deg] -2.445, xi [deg] -0.809, beta [deg] +1.057 (seed 1)
            // NONE  mean-std error vtas [mps] +1.52e-3 +1.11e-1, alpha [deg] +2.15e-3 +1.14e-1, beta [deg] -1.65e-3 +1.29e-1, Hp [m] +4.71e-2 +5.53e0
            // XACT  mean-std error rtov [deg] +6.79e-1 +3.53e-1, w_nbb [dps] +6.75e-2 +3.41e-2, bias_gyr [dps] +9.85e-3 +1.17e-2, bias_mag [nT] +1.20e+1 +1.98e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -2.462, theta [deg] -2.454, xi [deg] -0.777, beta [deg] +1.097 (seed 1)
            // mean-std error rtov [deg] +4.98e-1 +3.26e-1, w_nbb [dps] +6.77e-2 +3.44e-2, bias_gyr [dps] +8.19e-3 +1.12e-2, bias_mag [nT] +1.24e+1 +9.73e-1
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg,  2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 28.0, elv_theta_deg,  2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg, -2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 28.0, elv_theta_deg, -2.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,  0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 158: { // primary loops - simultaneous acceleration - deceleration and roll in - roll out maneuvers - starts at 3000 - 30
            // Calm: max errors are vtas [mps] +2.004, theta [deg] -0.585, xi [deg] +5.008, beta [deg] +1.487 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -2.480, theta [deg] -1.189, xi [deg] +5.380, beta [deg] +1.873 (seed 1)
            // NONE  mean-std error vtas [mps] +1.52e-3 +1.15e-1, alpha [deg] +2.16e-3 +1.12e-1, beta [deg] -1.66e-3 +1.32e-1, Hp [m] +9.99e-3 +5.02e+0
            // XACT  mean-std error rtov [deg] +5.47e-1 +4.07e-1, w_nbb [dps] +6.82e-2 +3.78e-2, bias_gyr [dps] +9.39e-3 +1.27e-2, bias_mag [nT] +8.13e+0 +3.41e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -2.567, theta [deg] -1.246, xi [deg] +5.394, beta [deg] +1.912 (seed 1)
            // mean-std error rtov [deg] +6.27e-1 +3.46e-1, w_nbb [dps] +6.85e-2 +3.82e-2, bias_gyr [dps] +9.07e-3 +1.19e-2, bias_mag [nT] +8.95e+0 +2.26e+0
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg, 0.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 34.0, elv_theta_deg, 0.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg, 0.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 28.0, elv_theta_deg, 0.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 26.0, elv_theta_deg, 0.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 28.0, elv_theta_deg, 0.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, 0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 159: { // primary loops - simultaneous push-pull and roll in - roll out maneuvers - starts at 3000 - 30
            // Calm: max errors are vtas [mps] +0.427, theta [deg] +2.002, xi [deg] +5.006, beta [deg] +1.253 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] +0.941, theta [deg] -2.446, xi [deg] +5.268, beta [deg] +1.564 (seed 1)
            // NONE  mean-std error vtas [mps] +1.52e-3 +1.04e-1, alpha [deg] +2.16e-3 +1.15e-1, beta [deg] -1.66e-3 +1.31e-1, Hp [m] +2.40e-2 +5.74e+0
            // XACT  mean-std error rtov [deg] +4.15e-1 +3.42e-1, w_nbb [dps] +6.89e-2 +4.17e-2, bias_gyr [dps] +7.79e-3 +1.22e-2, bias_mag [nT] +7.81e+0 +3.34e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] +0.990, theta [deg] -2.454, xi [deg] +5.286, beta [deg] +1.609 (seed 1)
            // mean-std error rtov [deg] +3.63e-1 +3.25e-1, w_nbb [dps] +6.93e-2 +4.25e-2, bias_gyr [dps] +6.59e-3 +1.14e-2, bias_mag [nT] +8.37e+0 +1.78e+0
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  2.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  2.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -2.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -4.0, ail_xi_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg, -2.0, ail_xi_deg,  -5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 160: { // primary loops - simultaneous acceleration-deceleration, push-pull and roll in - roll out maneuvers - starts at 3000 - 30
            // Calm: max errors are vtas [mps] -2.009, theta [deg] +2.002, xi [deg] -5.005, beta [deg] -1.444 (no turbulence) (GEO SPH, DeltaT = 10)
            // Open: max errors are vtas [mps] -2.562, theta [deg] -2.421, xi [deg] -5.105, beta [deg] -1.672 (seed 1)
            // NONE  mean-std error vtas [mps] +1.52e-3 +1.16e-1, alpha [deg] +2.11e-3 +1.15e-1, beta [deg] -1.64e-3 +1.29e-1, Hp [m] -7.12e-1 +5.10e+0
            // XACT  mean-std error rtov [deg] +4.21e-1 +3.59e-1, w_nbb [dps] +6.86e-2 +4.12e-2, bias_gyr [dps] +8.18e-3 +1.23e-2, bias_mag [nT] +7.16e+0 +3.52e+0
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -2.562, theta [deg] -2.424, xi [deg] -5.151, beta [deg] -1.748 (seed 1)
            // mean-std error rtov [deg] +4.08e-1 +3.17e-1, w_nbb [dps] +6.90e-2 +4.18e-2, bias_gyr [dps] +7.33e-3 +1.15e-2, bias_mag [nT] +6.44e+0 +3.22e+0
            auto Pguid = new control::guid(9, 280.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  0.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg,  2.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg,  2.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 32.0, elv_theta_deg,  2.0, ail_xi_deg,  15.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,  10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 28.0, elv_theta_deg,  2.0, ail_xi_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_theta_deg,  4.0, ail_xi_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 161: { // secondary loops - evaluation of Hp target in straight constant speed level flight - starts at 3000 - 30
            // Open: max errors are vtas [mps] +9.42e-1 theta [deg] +7.54e-1 xi [deg] -4.56e-1 beta [deg] -1.01e+0 Hp [m] +7.75e+0
            //       mean-std error vtas [mps] +2.03e-3 +1.03e-1, alpha [deg] +5.47e-3 +1.13e-1, beta [deg] -3.39e-3 +1.30e-1, Hp [m] +1.92e-2 +8.80e-1
            //       mean-std error rtov [deg] +7.25e-1 +3.27e-1, w_nbb [dps] +6.61e-2 +2.83e-2, bias_gyr [dps] +1.19e-2 +1.27e-2, bias_mag [nT] +1.07e+1 +9.67e-1
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] +9.39e-1 theta [deg]-1.39e+0 xi [deg] +5.45e-1 beta [deg] -1.06e+0 Hp [m] +8.64e+0
            // mean-std error rtov [deg] +6.96e-1 +3.47e-1, w_nbb [dps] +7.95e-2 +3.71e-2, bias_gyr [dps] +1.17e-2 +1.27e-2, bias_mag [nT] +1.10e+1 +1.21e+0
            auto Pguid = new control::guid(1, 200.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid); // 30.0 // 200.0
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_t_sec, 200.0));
            return Pguid;
        }
        case 162: { // secondary loops - evaluation of Hp target in constant absolute bearing constant speed level flight - starts at 3000 - 30
            // Open: max errors are vtas [mps] -0.968 theta [deg] -0.512 xi [deg] -0.704 beta [deg] -0.975 Hp [m] +1.10 chi [deg] -0.197
            //       mean-std error vtas [mps] +2.17e-3 +1.05e-1, alpha [deg] +5.46e-3 +1.16e-1, beta [deg] -3.39e-3 +1.29e-1, Hp [m] +1.27e-2 +4.880e-1
            //       mean-std error rtov [deg] +6.93e-1 +3.06e-1, w_nbb [dps] +6.38e-2 +2.71e-2, bias_gyr [dps] +1.15e-2 +1.19e-2, bias_mag [nT] +7.62e+0 +9.65e-1
            // Closed (no wind, zero weather, default constant magnetism, init time -10):
            // max error vtas [mps] -0.940 theta [deg] -1.394 xi [deg] -0.795 beta [deg] -1.069 Hp [m] +8.702 chi [deg] -0.457
            // mean-std error rtov [deg] +5.73e-1 +3.20e-1, w_nbb [dps] +7.93e-2 +3.70e-2, bias_gyr [dps] +9.73e-3 +1.10e-2, bias_mag [nT] +7.82e+0 +1.04e+0
            auto Pguid = new control::guid(1, 200.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, 0.0, rud_beta_deg, 0., trgg_t_sec, 200.0));
            return Pguid;
        }
        case 163: { // secondary loops - altitude changes (both Hp and h with CONNECTIONS using gammaTAS) - starts at 3000 - 30
            // IF I DO NOT WANT CONNECTIONS, I NEED TO FIND A WAY TO DEACTIVATE THE INTEGRAL CONTROL
            // PROBLEM WITH CONNECTIONS IS OVERSHOOTING AS THE AIRCRAFT ARRIVES TO THE TARGET AND ONLY REACTS ONCE IT HAS REACHED IT, WITH INERTIA CARRYING IT OVER THE VALUE
            // Open: max errors are vtas [mps] +1.489, theta [deg] +8.156, xi [deg] -0.570, beta [deg] -1.035, Hp [m] +1.105, h [m] -0.611, gammaTAS [deg] +1.703 (seed 1)
            auto Pguid = new control::guid(9, 360.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,         3000.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_t_sec,        40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg,    2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Hp_m,       3050.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,         3050.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec,   30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg,   -2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Hp_m,       3000.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m,         3000.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec,   30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg,   -2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_h_m,        2950.0));
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_h_m,          2950.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec,   30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg,    2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_h_m,        3000.0));
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_h_m,          3000.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            return Pguid;
            }
        case 164: { // secondary loops - gammaTAS changes - starts at 3000 - 30
            // maximum errors are vtas [mps] -1.006, theta [deg] +2.828, xi [deg] -0.553, beta [deg] -1.044, gammaTAS [deg] +2.763 (seed 1)
            auto Pguid = new control::guid(5, 160.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 4.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 2.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 4.0, ail_xi_deg, 0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 165: { // secondary loops - bearing changes (chi, chiTAS, and psi with CONNECTIONS using xi) - starts at 3000 - 30
            // IF I DO NOT WANT CONNECTIONS, I NEED TO FIND A WAY TO DEACTIVATE THE INTEGRAL CONTROL
            // PROBLEM WITH CONNECTIONS IS OVERSHOOTING AS THE AIRCRAFT ARRIVES TO THE TARGET AND ONLY REACTS ONCE IT HAS REACHED IT, WITH INERTIA CARRYING IT OVER THE VALUE
            // maximum errors are vtas [mps] -0.976, theta [deg] +0.597, xi [deg] -11.639, beta [deg] +3.283,
            // Hp [m] +1.102, chi [deg] -1.132, psi [deg] -4.548, chiTAS [deg] -1.344 (seed 1)
            auto Pguid = new control::guid(13, 300.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg,      0.0, rud_beta_deg, 0., trgg_t_sec,       40.0));
            Pguid->add_op(1,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,      10.0, rud_beta_deg, 0., trgg_chi_deg,     45.0));
            Pguid->add_op(2,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg,     45.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            Pguid->add_op(3,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,     -10.0, rud_beta_deg, 0., trgg_chi_deg,    -25.0));
            Pguid->add_op(4,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg,    -25.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            Pguid->add_op(5,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,     -10.0, rud_beta_deg, 0., trgg_chi_deg,    -50.0)); // no chiTAs trigger
            Pguid->add_op(6,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chiTAS_deg, -50.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            Pguid->add_op(7,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,      10.0, rud_beta_deg, 0., trgg_chi_deg,    -10.0)); // no chiTAS trigger
            Pguid->add_op(8,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chiTAS_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            Pguid->add_op(9,  new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,      10.0, rud_beta_deg, 0., trgg_psi_deg,     30.0));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg,     30.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_xi_deg,     -10.0, rud_beta_deg, 0., trgg_psi_deg,    -35.0));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg,    -35.0, rud_beta_deg, 0., trgg_Deltat_sec,  30.0));
            return Pguid;
        }
        case 166: { // secondary loops - muTAS changes - starts at 3000 - 30
            // maximum errors are vtas [mps] -0.886, theta [deg] +1.147, xi [deg] +22.518, beta [deg] +4.554, gammaTAS [deg] +1.128, muTAS [deg] -20.023 (seed 1)
            auto Pguid = new control::guid(5, 160.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_muTAS_deg,   0.0, rud_beta_deg, 0., trgg_t_sec,      40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_muTAS_deg,  15.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_muTAS_deg,   5.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_muTAS_deg, -15.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_gammaTAS_deg, 0.0, ail_muTAS_deg,   0.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }
        case 167: { // secondary loops - bearing changes (chi, chiTAS, and psi with CONNECTIONS using muTAS) - starts at 3000 - 30
            // IF I DO NOT WANT CONNECTIONS, I NEED TO FIND A WAY TO DEACTIVATE THE INTEGRAL CONTROL
            // PROBLEM WITH CONNECTIONS IS OVERSHOOTING AS THE AIRCRAFT ARRIVES TO THE TARGET AND ONLY REACTS ONCE IT HAS REACHED IT, WITH INERTIA CARRYING IT OVER THE VALUE
            // maximum errors are vtas [mps] -0.972, theta [deg] +0.559, xi [deg] +15.616, beta [deg] -3.867,
            // Hp [m] +1.102, chi [deg] +1.166, psi [deg] +3.724, muTAS [deg] +10.619, chiTAS [deg] +2.682 (seeds1)
            auto Pguid = new control::guid(13, 300.0, Oseeder.provide_seed(math::seeder::seeder_guid), case_guid);
            Pguid->_sti_id = st::logic::sti_h3000_tas30_psi00;
            Pguid->add_op(0, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, 0.0, rud_beta_deg, 0., trgg_t_sec, 40.0));
            Pguid->add_op(1, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, 10.0, rud_beta_deg, 0., trgg_chi_deg, 45.0));
            Pguid->add_op(2, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, 45.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(3, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, -10.0, rud_beta_deg, 0., trgg_chi_deg, -25.0));
            Pguid->add_op(4, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chi_deg, -25.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(5, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, -10.0, rud_beta_deg, 0., trgg_chi_deg, -50.0)); // no chiTAs trigger
            Pguid->add_op(6, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chiTAS_deg, -50.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(7, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, 10.0, rud_beta_deg, 0., trgg_chi_deg, -10.0)); // no chiTAS trigger
            Pguid->add_op(8, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_chiTAS_deg, -10.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(9, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, 10.0, rud_beta_deg, 0., trgg_psi_deg, 30.0));
            Pguid->add_op(10, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg, 30.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            Pguid->add_op(11, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_muTAS_deg, -10.0, rud_beta_deg, 0., trgg_psi_deg, -35.0));
            Pguid->add_op(12, new control::guid_op(thr_vtas_mps, 30.0, elv_Hp_m, 3000.0, ail_psi_deg, -35.0, rud_beta_deg, 0., trgg_Deltat_sec, 30.0));
            return Pguid;
        }

        default:
            throw std::runtime_error("Guidance case not available");
    }
}
/* returns pointer to guidance objectives based on seeder and case identifier. */

void control::guid::create_text(std::ostream& Ostream) const {
    switch (_case_guid) {
        case 1:
            Ostream << std::endl << "INTENT:" << std::endl << std::endl;
            Ostream << "Initial  chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn       t [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[1]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Final    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Initial vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_THR) << std::endl;
            Ostream << "Speed Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_trg_val()                     << std::endl;
            Ostream << "Final   vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_val(control::logic::cntr_THR) << std::endl;
            Ostream << "Initial   Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_ELV) << std::endl;
            Ostream << "Alt   Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[3]->get_guid_trg_val()                     << std::endl;
            Ostream << "Alt gammaTAS [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[4]->get_guid_val(control::logic::cntr_ELV) << std::endl;
            Ostream << "Final     Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_val(control::logic::cntr_ELV) << std::endl;
            break;
        case 2:
            Ostream << std::endl << "INTENT:" << std::endl << std::endl;
            Ostream << "Initial  chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn       t [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[1]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Final    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Initial vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_THR) << std::endl;
            Ostream << "Initial   Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_ELV) << std::endl;
            break;
        case 3:
            Ostream << std::endl << "INTENT:" << std::endl << std::endl;
            Ostream << "Initial    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #1      t [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #1     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[1]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #1    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "vtas           [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_THR) << std::endl;
            Ostream << "Hp               [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_ELV) << std::endl;
            Ostream << "Turn #2 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #2     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[3]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #2    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[4]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #3 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[4]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #3     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #3    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[6]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #4 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[6]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #4     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[7]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #4    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[8]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #5 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[8]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #5     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[9]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #5    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[10]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #6 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[10]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #6     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[11]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #6    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[12]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #7 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[12]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #7     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[13]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #7    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[14]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #8 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[14]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #8     xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[15]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #8    chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[16]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            break;
        case 4:
            Ostream << std::endl << "INTENT:" << std::endl << std::endl;
            Ostream << "Initial     chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Initial    vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_THR) << std::endl;
            Ostream << "Initial      Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_val(control::logic::cntr_ELV) << std::endl;

            Ostream << "Turn #1       t [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[0]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #1      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[1]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #1     chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_val(control::logic::cntr_AIL) << std::endl;

            Ostream << "Speed #1 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[2]->get_guid_trg_val()                     << std::endl;
            Ostream << "Speed #1   vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[3]->get_guid_val(control::logic::cntr_THR) << std::endl;

            Ostream << "Turn #2  Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[3]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #2      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[4]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #2     chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_val(control::logic::cntr_AIL) << std::endl;

            Ostream << "Alt #1   Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[5]->get_guid_trg_val()                     << std::endl;
            Ostream << "Alt #1       Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[7]->get_guid_val(control::logic::cntr_ELV) << std::endl;

            Ostream << "Turn #3  Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[7]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #3      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[8]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #3     chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[9]->get_guid_val(control::logic::cntr_AIL) << std::endl;

            Ostream << "Speed #2 Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[9]->get_guid_trg_val()                     << std::endl;
            Ostream << "Speed #2   vtas [mps]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[10]->get_guid_val(control::logic::cntr_THR) << std::endl;

            Ostream << "Turn #4  Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[10]->get_guid_trg_val()                     << std::endl;
            Ostream << "Turn #4      xi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[11]->get_guid_val(control::logic::cntr_AIL) << std::endl;
            Ostream << "Turn #4     chi [deg]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[12]->get_guid_val(control::logic::cntr_AIL) << std::endl;

            Ostream << "Alt #2   Deltat [sec]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[12]->get_guid_trg_val()                     << std::endl;
            Ostream << "Alt #2       Hp   [m]: " << std::fixed << std::setw(7)  << std::setprecision(1) << std::showpos << _Vguid_op[14]->get_guid_val(control::logic::cntr_ELV) << std::endl;

            break;
        default:
            std::cout << "Guidance not described in text file." << std::endl;
            break;
    }
}
/* describe guidance in stream */

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
















