#include "geo.h"
#include "mag.h"
#include "coord.h"
#include "math/logic/share.h"
#include "math/vec/vec1.h"
#include "math/vec/vec2.h"
#include "math/pred/pred2v/f_table2V.h"
#include "ang/tools.h"
#include "ang/rotate/rotv.h"
#include "ang/rotate/euler.h"

#include <boost/filesystem.hpp>

// CLASS GEO
// =========
// =========

const double	env::geo::_g0_mps2 = 9.80665;
/* standard acceleration of free fall */
const double	env::geo::_GM_m3ps2 = 3986004.418e8;
/* Earth gravitational constant (atmosphere included)[m3 / s2] */
const double	env::geo::_w_irsecefecefiii_rps = 7.292115e-5;
/* Earth rotation velocity */
const double	env::geo::_a_m          = 6378137.0;
/* ellipsoid semimajor */
const double	env::geo::_f            = 1 / 298.257223563;
/* ellipsoid flattening [-] */
const double	env::geo::_C20          = -0.484166774985e-3;
/* EGM96 2nd degree zonal harmonical */
const double    env::geo::_gcMSLE_mps2  = 9.7803253359;
/* theoretical WGS84 gravity value at Equator */
const double    env::geo::_gcMSLP_mps2  = 9.8321849378;
/* theoretical WGS84 gravity value at poles */

const std::string env::geo::_sgeom = "/H2h_32_geom.txt";
const std::string env::geo::_sgeop = "/H2h_32_geop.txt";
const std::string env::geo::_slat = "/H2h_32_lat.txt";
const std::string env::geo::_ssgeom = "/htoH_32_geom.txt";
const std::string env::geo::_ssgeop = "/htoH_32_geop.txt";
const std::string env::geo::_sslat = "/htoH_32_lat.txt";
/* names required internally */

void env::geo::set_ready() {
    _b_m = _a_m * (1. - _f);
    _a2 = pow(_a_m, 2);
    _b2 = pow(_b_m, 2);
    _e2 = _f * (2. - _f);
    _e = sqrt(_e2);
    _g2 = (_a2 - _b2) / _b2;
    _g = sqrt(_g2);
    _m = (pow(_w_irsecefecefiii_rps,2) * pow(_a_m,2) * _b_m) / _GM_m3ps2;
    _k = (_b_m * _gcMSLP_mps2) / (_a_m * _gcMSLE_mps2) - 1.0;
}
/* updates all private attributes based on semimajor and flattening */

env::geo::geo(env::logic::ZONE_ID zone_id)
: _gen(1), _dist_normal(0.,1.), _dist_uniform(-179, 180),
  _sigma_gc_mps2(0.), _sigma_gc_deg(0.), _Delta_gc_mps2(0.), _Delta_gc_rad(0.), _Delta_gc_bearing_rad(0.),
  _R_gc(Eigen::Matrix3d::Identity()) {
    // compute other class constants
    this->set_ready();

    boost::filesystem::path p0(math::share::condor_input);
    boost::filesystem::path p1("geopotential");
    std::string file = (p0 / p1).string();

    // H->h conversion
    double H_init, H_step;
    unsigned short H_size;
    std::ifstream H_mystream((file + _sgeop).c_str());
    H_mystream >> H_init >> H_step >> H_size;
    H_mystream.close();
    auto PH_values = new math::vec1(H_size);
    for (unsigned short i = 0; i != H_size; ++i) {
        PH_values->set(i, H_init + i * H_step);
    }
    double phi_init, phi_step;
    unsigned short phi_size;
    std::ifstream phi_mystream((file + _slat).c_str());
    phi_mystream >> phi_init >> phi_step >> phi_size;
    phi_mystream.close();
    auto Pphi_values = new math::vec1(phi_size);
    for (unsigned short i = 0; i != phi_size; ++i) {
        Pphi_values->set(i, phi_init + i * phi_step);
    }
    auto Ph_values = new math::vec2(file + _sgeom);
    _H2h = new math::f_table2V(Pphi_values, PH_values, Ph_values, math::logic::lagrange_first);

    // h->H conversion
    double h_init, h_step;
    unsigned short h_size;
    std::ifstream h_mystream((file + _ssgeom).c_str());
    h_mystream >> h_init >> h_step >> h_size;
    h_mystream.close();
    auto Qh_values = new math::vec1(h_size);
    for (unsigned short i = 0; i != h_size; ++i) {
        Qh_values->set(i, h_init + i * h_step);
    }
    double varphi_init, varphi_step;
    unsigned short varphi_size;
    std::ifstream varphi_mystream((file + _sslat).c_str());
    varphi_mystream >> varphi_init >> varphi_step >> varphi_size;
    varphi_mystream.close();
    auto Qphi_values = new math::vec1(varphi_size);
    for (unsigned short i = 0; i != varphi_size; ++i) {
        Qphi_values->set(i, varphi_init + i * varphi_step);
    }
    auto QH_values = new math::vec2(file + _ssgeop);
    _htoH = new math::f_table2V(Qphi_values, Qh_values, QH_values, math::logic::lagrange_first);

    _Pmag = env::mag::create_mag(zone_id, env::logic::realism_grav_no_magn_no, _dist_normal, _gen);
}
/* constructor based on location. Model and real gravitation and magnetic fields coincide. */

env::geo::geo(env::logic::ZONE_ID zone_id, env::logic::REALISM_ID realism_id, const int& seed)
: _gen(seed), _dist_normal(0.,1.), _dist_uniform(-179, 180) {
    // compute other class constants
    this->set_ready();

    boost::filesystem::path p0(math::share::condor_input);
    boost::filesystem::path p1("geopotential");
    std::string file = (p0 / p1).string();

    // H->h conversion
    double H_init, H_step;
    unsigned short H_size;
    std::ifstream H_mystream((file + _sgeop).c_str());
    H_mystream >> H_init >> H_step >> H_size;
    H_mystream.close();
    auto PH_values = new math::vec1(H_size);
    for (unsigned short i = 0; i != H_size; ++i) {
        PH_values->set(i, H_init + i * H_step);
    }
    double phi_init, phi_step;
    unsigned short phi_size;
    std::ifstream phi_mystream((file + _slat).c_str());
    phi_mystream >> phi_init >> phi_step >> phi_size;
    phi_mystream.close();
    auto Pphi_values = new math::vec1(phi_size);
    for (unsigned short i = 0; i != phi_size; ++i) {
        Pphi_values->set(i, phi_init + i * phi_step);
    }
    auto Ph_values = new math::vec2(file + _sgeom);
    _H2h = new math::f_table2V(Pphi_values, PH_values, Ph_values, math::logic::lagrange_first);

    // h->H conversion
    double h_init, h_step;
    unsigned short h_size;
    std::ifstream h_mystream((file + _ssgeom).c_str());
    h_mystream >> h_init >> h_step >> h_size;
    h_mystream.close();
    auto Qh_values = new math::vec1(h_size);
    for (unsigned short i = 0; i != h_size; ++i) {
        Qh_values->set(i, h_init + i * h_step);
    }
    double varphi_init, varphi_step;
    unsigned short varphi_size;
    std::ifstream varphi_mystream((file + _sslat).c_str());
    varphi_mystream >> varphi_init >> varphi_step >> varphi_size;
    varphi_mystream.close();
    auto Qphi_values = new math::vec1(varphi_size);
    for (unsigned short i = 0; i != varphi_size; ++i) {
        Qphi_values->set(i, varphi_init + i * varphi_step);
    }
    auto QH_values = new math::vec2(file + _ssgeop);
    _htoH = new math::f_table2V(Qphi_values, Qh_values, QH_values, math::logic::lagrange_first);

    switch (realism_id) {
        case env::logic::realism_grav_yes_magn_yes:
        case env::logic::realism_grav_yes_magn_no:
            _sigma_gc_mps2 = 0.0001;
            _sigma_gc_deg  = 0.0028;
            break;
        case env::logic::realism_grav_no_magn_yes:
        case env::logic::realism_grav_no_magn_no:
            _sigma_gc_mps2 = 0.;
            _sigma_gc_deg  = 0.;
            break;
        case env::logic::realism_size:
        default:
            throw std::runtime_error("Realism option not available");
    }

    _Delta_gc_mps2 = _sigma_gc_mps2 * _dist_normal(_gen);
    _Delta_gc_rad  = _sigma_gc_deg  * _dist_normal(_gen) * math::constant::D2R();
    _Delta_gc_bearing_rad = _dist_uniform(_gen) * math::constant::D2R();

    _R_gc = ang::rotv(Eigen::Vector3d(cos(_Delta_gc_bearing_rad), sin(_Delta_gc_bearing_rad), 0.), _Delta_gc_rad);

    _Pmag = env::mag::create_mag(zone_id, realism_id, _dist_normal, _gen); // IMPORTANT: Pass distribution to magnetism after using it (for repeatability purposes)
}
/* constructor based on location, realism enumerator, and realism seed. Seed only applies to realistic (not model) gravitation and magnetism. */

env::geo::~geo() {
    delete _Pmag;
    delete _H2h;
    delete _htoH;
}
/* destructor */

/* ===== ===== ===== Coordinate Change ===== ===== ===== */
/* ===================================================== */

env::cartesian_coord env::geo::geocentric2cartesian(const env::geocentric_coord& sisu) const {
	return {sisu.get_r_m() * sin(sisu.get_theta_rad()) * cos(sisu.get_lambda_rad()),
		sisu.get_r_m() * sin(sisu.get_theta_rad()) * sin(sisu.get_lambda_rad()),
		sisu.get_r_m() * cos(sisu.get_theta_rad())};
}
/* transforms geocentric coordinates into cartesian ones. */

env::geocentric_coord env::geo::cartesian2geocentric(const env::cartesian_coord& cisu) const {
    double lambda_rad = atan2(cisu.get_x2_m(), cisu.get_x1_m());
    ang::tools::correct_longitude_rad(lambda_rad);
	return {atan2(sqrt(pow(cisu.get_x1_m(), 2.) + pow(cisu.get_x2_m(), 2.)), cisu.get_x3_m()),
		lambda_rad,
		sqrt(pow(cisu.get_x1_m(), 2.) + pow(cisu.get_x2_m(), 2.) + pow(cisu.get_x3_m(), 2.))};
}
/* transforms cartesian coordinates into geocentric ones. */

env::cartesian_coord env::geo::geodetic2cartesian(const env::geodetic_coord& gisu, const double& N_m) const {
    return {(gisu.get_h_m() + N_m) * cos(gisu.get_phi_rad()) * cos(gisu.get_lambda_rad()),
            (gisu.get_h_m() + N_m) * cos(gisu.get_phi_rad()) * sin(gisu.get_lambda_rad()),
            (gisu.get_h_m() + (N_m * (1. - _e2))) * sin(gisu.get_phi_rad())};
}
/* transforms geodetic coordinates into cartesian ones. Requires radius of curvature of meridian. */

env::geocentric_coord env::geo::geodetic2geocentric(const env::geodetic_coord& gisu, const double& N_m) const {
    double sinphi = sin(gisu.get_phi_rad());
    double cosphi = cos(gisu.get_phi_rad());
    return {atan2((N_m + gisu.get_h_m()) * cosphi, (N_m * (1 - _e2) + gisu.get_h_m()) * sinphi),
            gisu.get_lambda_rad(),
            sqrt(pow((N_m + gisu.get_h_m()),2) * pow(cosphi,2) + pow((N_m * (1 - _e2) + gisu.get_h_m()),2) * pow(sinphi,2))};
}
/* transforms geodetic coordinates into geocentric ones. Requires radius of curvature of meridian. */

env::geodetic_coord env::geo::geocentriccartesian2geodetic(const env::geocentric_coord& sisu, const env::cartesian_coord& cisu) const {
    double r = sqrt(pow(cisu.get_x1_m(), 2.) + pow(cisu.get_x2_m(), 2.));
    double phi0_rad = -(sisu.get_theta_rad() - math::constant::PIHALF());
    double x3 = cisu.get_x3_m();
    double N_m = this->radius_vert(phi0_rad);
    double phi_rad, delta = 1.0, h_m = 0;
    do {
        phi_rad = atan2(x3, r * (1 - (N_m * _e2) / (N_m + h_m)));
        N_m = this->radius_vert(phi_rad);
        if (fabs(phi_rad) > math::constant::PI() / 3) {
            h_m = x3 / sin(phi_rad) - N_m * (1 - _e2);
        }
        else {
            h_m = r / cos(phi_rad) - N_m;
        }
        delta = fabs(phi_rad - phi0_rad);
        phi0_rad = phi_rad;
    } while (delta > 1e-14);

    // I fix here the longitude that comes from the geocentric coordinates so it is [0,360) instead of (-180,180].
    // Right thing to do would be to fix generation of geocetric coordintes instead of here but I am afraid
    // other things such as gravity may change.
    double lambda_rad = sisu.get_lambda_rad();
    ang::tools::correct_longitude_rad(lambda_rad);
    return {lambda_rad, phi_rad, h_m};
}
/* transforms geocentric plus cartesian coordinates into geodetic ones. */

Eigen::Matrix3d env::geo::jacobian_cartesian_geodetic(const env::geodetic_coord& gisu, const double& N_m) const {
    double sinlambda = sin(gisu.get_lambda_rad());
    double coslambda = cos(gisu.get_lambda_rad());
    double sinphi    = sin(gisu.get_phi_rad());
    double cosphi    = cos(gisu.get_phi_rad());
    double factor    = (N_m + gisu.get_h_m());

    Eigen::Matrix3d res;

    res(0,0) = - factor * cosphi * sinlambda;
    res(0,1) = - factor * sinphi * coslambda;
    res(0,2) = + cosphi * coslambda;
    res(1,0) = + factor * cosphi * coslambda;
    res(1,1) = - factor * sinphi * sinlambda;
    res(1,2) = + cosphi * sinlambda;
    res(2,0) = 0.;
    res(2,1) = + (N_m * (1. - _e2) + gisu.get_h_m()) * cosphi;
    res(2,2) = + sinphi;

    return res;
}
/* computes the jacobian of the cartesian coordinates with respect to the geodetic ones */

/* ===== ===== ===== Radii of Curvature ===== ===== ===== */
/* ====================================================== */

double env::geo::radius_vert(const double& phi_rad) const {
    return _a_m / sqrt(1 - _e2 * pow(sin(phi_rad), 2));
}
/* returns the radius of curvature of prime vertical [m] based on latitude */

double env::geo::radius_mer(const double& phi_rad, const double& N_m) const {
    return N_m / (1 + _g2 * pow(cos(phi_rad), 2));
}
/* returns the radius of curvature of meridian [m] based on latitude */

/* ===== ===== ===== Geodetic <--> Geopotential Altitude Conversion ===== ===== ===== */
/* ================================================================================== */

double env::geo::H2h(const double& H_m, const double& phi_rad) const {
    return _H2h->value(fabs(phi_rad), H_m);
}
/* returns geodetic altitude [m] based on geopotential altitude and latitude */

double env::geo::htoH(const double& h_m, const double& phi_rad) const {
    return _htoH->value(fabs(phi_rad), h_m);
}
/* returns geopotential altitude [m] based on geodetic altitude and latitude */

/* ===== ===== ===== Gravitation, Centrifugal, and Gravity ===== ===== ===== */
/* ========================================================================= */

Eigen::Vector3d env::geo::compute_gravitation_n(const env::geodetic_coord& x_gdt_rad_m, const double& N_m) const {
    // transform geodetic// coordinates into geocentric ones
    env::geocentric_coord x_gct_rad_m = this->geodetic2geocentric(x_gdt_rad_m, N_m);
    // obtain gravitation in spherical coordinates (SRS)
    double sin0 = sin(x_gct_rad_m.get_theta_rad());
    double cos0 = cos(x_gct_rad_m.get_theta_rad());
    Eigen::Vector3d g_srs_mps2(-_GM_m3ps2 / pow(x_gct_rad_m.get_r_m(),2) * pow(_a_m / x_gct_rad_m.get_r_m(), 2) * _C20 * 3.0 * sqrt(5.0) * sin0 * cos0,
                               0.0,
                               -_GM_m3ps2 / pow(x_gct_rad_m.get_r_m(),2) * (1 + 3 * pow(_a_m / x_gct_rad_m.get_r_m(),2) * _C20 * (sqrt(5.) / 2.) * (3 * pow(cos0,2) - 1)));
    // Euler angles to rotate from SRS to NED
    ang::euler euler_srsned_rad(0.0, -x_gct_rad_m.get_theta_rad() - x_gdt_rad_m.get_phi_rad() - math::constant::PIHALF(), 0.0);
    // rotate gravitation to NED
    Eigen::Vector3d g_ned_mps2 = euler_srsned_rad / g_srs_mps2;
    return g_ned_mps2;
}
/* computes the gravitation acceleration viewed in NED based on geodetic position and radius of curvature of prime vertical.
 * Based on ellipsoidal EGM96 model, which is not orthogonal to ellipsoid. */

Eigen::Vector3d env::geo::compute_centrifugal_n(const env::geodetic_coord& x_gdt_rad_m, const double& N_m) const {
//	// transform geodetic coordinates into geocentric ones
//	env::geocentric_coord x_gct_rad_m = this->geodetic2geocentric(x_gdt_rad_m, N_m);
//	// obtain centrifugal acceleration in speherical coordinate(SRS)
//	double sin0 = sin(x_gct_rad_m.get_theta_rad());
//	double cos0 = cos(x_gct_rad_m.get_theta_rad());
//	double temp = pow(_w_irsecefecefiii_rps, 2) * x_gct_rad_m.get_r_m() * sin0;
//    Eigen::Vector3d f_srs_mps2(cos0 * temp, 0.0, sin0 * temp);
//	// Euler angles to rotate from SRS to NED
//	ang::euler euler_srsned_rad(0.0, -x_gct_rad_m.get_theta_rad() - x_gdt_rad_m.get_phi_rad() - math::constant::PIHALF(), 0.0);
//	// rotate centrifugal acceleration to NED
//    Eigen::Vector3d f_ned_mps2 = euler_srsned_rad / f_srs_mps2;
//	return f_ned_mps2;
    return - this->compute_transport_n(N_m, x_gdt_rad_m.get_h_m(), x_gdt_rad_m.get_phi_rad());
}
/* computes the centrifugal acceleration viewed in NED based on geodetic position and radius of curvature of prime vertical.
 * Employs true centrifugal effect. */

Eigen::Vector3d env::geo::compute_gravity_n_model(const env::geodetic_coord& x_gdt_rad_m) const {
    double sinphi = sin(x_gdt_rad_m.get_phi_rad());
    double pow2sinphi = pow(sinphi, 2);

    double gcMSL_mps2 = _gcMSLE_mps2 * (1 + _k * pow2sinphi) / sqrt(1 - _e2 * pow2sinphi);
    double gc_mps2 = gcMSL_mps2 * (1 - 2.0 / _a_m * (1.0 + _f + _m - 2.0 * _f * pow2sinphi) * x_gdt_rad_m.get_h_m() + 3 / pow(_a_m,2) * pow(x_gdt_rad_m.get_h_m(), 2));
    return Eigen::Vector3d{0., 0., gc_mps2};

    ///////////////////////////////////////////////////////////////////////////////
    // MAKE THIS FORMULA MORE EFFICIENT BY COMBINING FACTORS
    ///////////////////////////////////////////////////////////////////////////////
}
/* computes the MODEL gravity acceleration viewed in NED based on geodetic position.
 * Based on ellipsoidal WGS84 model, which is orthogonal to ellipsoid. Does not add the
 * gravitation plus the centrifugal effect but employs the Somigliana gravity formula,
 * which is always orthogonal to the ellipsoid. */

Eigen::Vector3d env::geo::compute_gravity_n_truth(const env::geodetic_coord& x_gdt_rad_m) const {
    Eigen::Vector3d gc_n_mps2_model = this->compute_gravity_n_model(x_gdt_rad_m);
    double gc_mps2_norm = gc_n_mps2_model.norm();
    return _R_gc * gc_n_mps2_model* (1.0 + _Delta_gc_mps2 / gc_mps2_norm);
}
/* computes the REAL gravity acceleration viewed in NED based on geodetic position.
* It first calls the compute_gravity_ned_model and then modifies it to add stochastic realism. */

Eigen::Vector3d env::geo::compute_gravity_n_truth(const Eigen::Vector3d& gc_n_mps2_model) const {
    double gc_mps2_norm = gc_n_mps2_model.norm();
    return _R_gc * gc_n_mps2_model* (1.0 + _Delta_gc_mps2 / gc_mps2_norm);
}
/* modifies the input MODEL gravity acceleration viewed in NED adding stochastic realism. */

void env::geo::create_text(std::ostream& Ostream, const env::geodetic_coord& x_gdt_rad_m) const {
    Eigen::Vector3d gc_n_mps2_model = this->compute_gravity_n_model(x_gdt_rad_m);
    Eigen::Vector3d gc_n_mps2_truth = this->compute_gravity_n_truth(gc_n_mps2_model);

    Ostream << std::endl << "GRAVITY ERROR [model-truth-diff]:" << std::endl << std::endl;
    Ostream << "North [mps2]: "
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(0)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_truth(0)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(0) - gc_n_mps2_truth(0)
            << std::endl;
    Ostream << "East  [mps2]: "
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(1)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_truth(1)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(1) - gc_n_mps2_truth(1)
            << std::endl;
    Ostream << "Down  [mps2]: "
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(2)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_truth(2)
            << std::fixed << std::setw(10) << std::setprecision(6) << std::showpos << gc_n_mps2_model(2) - gc_n_mps2_truth(2)
            << std::endl;
}
/* describe gravity model and realism in stream */

/* ===== ===== ===== Conversion between Ground Speed and Geodetic Coordinates Differentials ===== ===== ===== */
/* ========================================================================================================== */

Eigen::Vector3d env::geo::vned_to_xgdtdot(const Eigen::Vector3d& v_n_mps, const env::geodetic_coord& x_gdt_rad_m, const double& N_m, const double& M_m) {
    return {v_n_mps(1) / ((N_m + x_gdt_rad_m.get_h_m()) * cos(x_gdt_rad_m.get_phi_rad())), v_n_mps(0) / (M_m + x_gdt_rad_m.get_h_m()), - v_n_mps(2)};
}
/* computes the time derivative of the geodetic coordinates based on the absolute speed in NED and the radii of curvature */

Eigen::Vector3d env::geo::xgdtdot_to_vned(const Eigen::Vector3d& xdot_gdt_rad_m, const env::geodetic_coord& x_gdt_rad_m, const double& N_m, const double& M_m) {
    return {xdot_gdt_rad_m(1) * (M_m + x_gdt_rad_m.get_h_m()), xdot_gdt_rad_m(0) * (N_m + x_gdt_rad_m.get_h_m()) * cos(x_gdt_rad_m.get_phi_rad()),- xdot_gdt_rad_m(2)};
}
/* computes the absolute speed in NED based o the time derivative of the geodetic coordinates and the radii of curvature */

Eigen::Matrix3d env::geo::jac_xgdtdot_wrt_vned(const env::geodetic_coord& x_gdt_rad_m, const double& N_m, const double& M_m) {
    Eigen::Matrix3d res = Eigen::Matrix3d::Zero();

    res(0,1) = 1. / ((N_m + x_gdt_rad_m.get_h_m()) * cos(x_gdt_rad_m.get_phi_rad()));
    res(1,0) = 1. / (M_m + x_gdt_rad_m.get_h_m());
    res(2,2) = -1.;

    return res;
}
/* computes jacobian of the time derivative of the geodetic coordinates with respect to the absolute velocity in NED */

/* ===== ===== ===== Earth Angular Velocity ===== ===== ===== */
/* ========================================================== */

ang::so3_tangent env::geo::compute_wiee_rps() const {
    return {0., 0., _w_irsecefecefiii_rps};
}
/* computes the Earth angular velocity [rps] (ECEF with respect to IRS) viewed in ECEF */

ang::so3_tangent env::geo::compute_wien_rps(const double& phi_rad) const {
	return {_w_irsecefecefiii_rps * cos(phi_rad), 0.0, -_w_irsecefecefiii_rps * sin(phi_rad)};
}
/* computes the Earth angular velocity [rps] (ECEF with respect to IRS) viewed in NED */

Eigen::Matrix3d env::geo::jac_wiee_wrt_ve() const {
    return Eigen::Matrix3d::Zero();
}
/* computes the jacobian of the Earth angular velocity viewed in ECEF with respect to the ground velocity viewed in ECEF */

Eigen::Matrix3d env::geo::jac_wien_wrt_vn() const {
    return Eigen::Matrix3d::Zero();
}
/* computes the jacobian of the Earth angular velocity viewed in NED with respect to the ground velocity viewed in NED */

/* ===== ===== ===== Motion Angular Velocity ===== ===== ===== */
/* =========================================================== */

ang::so3_tangent env::geo::compute_wene_rps(const Eigen::Vector3d& v_n_mps, const double& N_m, const double& M_m, const double& lambda_rad, const double& phi_rad, const double& h_m) const {
    return {v_n_mps(0) * sin(lambda_rad) / (M_m + h_m), - v_n_mps(0) * cos(lambda_rad) / (M_m + h_m), v_n_mps(1) / cos(phi_rad) / (N_m + h_m)};
}
/* computes the motion angular velocity [rps] (NED with respect to ECEF) viewed in ECEF. */

ang::so3_tangent env::geo::compute_wene_rps(const double& lambdadot_rps, const double& phidot_rps, const double& lambda_rad) const {
    return {phidot_rps * sin(lambda_rad), - phidot_rps * cos(lambda_rad), lambdadot_rps};
}
/*< computes the motion angular velocity [rps] (NED with respect to ECEF) viewed in ECEF. */

ang::so3_tangent env::geo::compute_wene_rps(const double& N_m, const double& M_m, const double& lambda_rad, const double& phi_rad, const double& h_m, const Eigen::Vector3d& v_e_mps) const {
    double sinlambda = sin(lambda_rad);
    double coslambda = cos(lambda_rad);
    double sinphi    = sin(phi_rad);
    double cosphi    = cos(phi_rad);
    double v_i_n     = - sinphi * coslambda * v_e_mps(0) - sinphi * sinlambda * v_e_mps(1) + cosphi * v_e_mps(2);

    return {+ v_i_n * sinlambda / (M_m + h_m),
            - v_i_n * coslambda / (M_m + h_m),
            + (- sinlambda * v_e_mps(0) + coslambda * v_e_mps(1)) / (N_m + h_m) / cosphi};
}
/* computes the motion angular velocity [rps] (NED with respect to ECEF) viewed in ECEF */

ang::so3_tangent env::geo::compute_wenn_rps(const Eigen::Vector3d& v_n_mps, const double& N_m, const double& M_m, const double& phi_rad, const double& h_m) const {
    return {v_n_mps(1) / (N_m + h_m), -v_n_mps(0) / (M_m + h_m), -v_n_mps(1) * tan(phi_rad) / (N_m + h_m)};
}
/* computes the motion angular velocity [rps] (NED with respect to ECEF) viewed in NED. */

ang::so3_tangent env::geo::compute_wenn_rps(const double& lambdadot_rps, const double& phidot_rps, const double& phi_rad) const {
    return {lambdadot_rps * cos(phi_rad), - phidot_rps, - lambdadot_rps * sin(phi_rad)};
}
/* computes the motion angular velocity [rps] (NED with respect to ECEF) viewed in NED. */

Eigen::Matrix3d env::geo::jac_wene_wrt_ve(const double& N_m, const double& M_m, const double& lambda_rad, const double& phi_rad, const double& h_m) const {
    double sinlambda = sin(lambda_rad);
    double coslambda = cos(lambda_rad);
    double sinphi    = sin(phi_rad);
    double cosphi    = cos(phi_rad);
    double factor0   = + sinlambda / (M_m + h_m);
    double factor1   = - coslambda / (M_m + h_m);
    double factor2   = 1. / cosphi / (N_m + h_m);

    Eigen::Matrix3d res;
    res(0,0) = factor0 * (- sinphi * coslambda);
    res(0,1) = factor0 * (- sinphi * sinlambda);
    res(0,2) = factor0 * (+ cosphi);
    res(1,0) = factor1 * (- sinphi * coslambda);
    res(1,1) = factor1 * (- sinphi * sinlambda);
    res(1,2) = factor1 * (+ cosphi);
    res(2,0) = factor2 * (- sinlambda);
    res(2,1) = factor2 * (+ coslambda);
    res(2,2) = 0.;
    return res;
}
/* computes the jacobian of the motion angular velocity viewed in ECEF with respect to the ground velocity viewed in ECEF */

Eigen::Matrix3d env::geo::jac_wenn_wrt_vn(const double& N_m, const double& M_m, const double& phi_rad, const double& h_m) const {
    Eigen::Matrix3d res = Eigen::Matrix3d::Zero();
    res(0,1) = 1. / (N_m + h_m);
    res(1,0) = - 1. / (M_m + h_m);
    res(2,1) = - tan(phi_rad) / (N_m + h_m);
    return res;
}
/* computes the jacobian of the motion angular velocity viewed in NED with respect to the ground velocity viewed in NED */

/* ===== ===== ===== Coriolis Acceleration ===== ===== ===== */
/* ========================================================= */

Eigen::Vector3d env::geo::compute_coriolis_e(const Eigen::Vector3d& v_e_mps, const ang::so3_tangent& w_iee_rps) const {
    return w_iee_rps().cross(v_e_mps) * 2.0;
}
/* computes the Coriolis acceleration [mps2] viewed in ECEF */

Eigen::Vector3d env::geo::compute_coriolis_e(const Eigen::Vector3d& v_e_mps) const {
    return {- 2. * _w_irsecefecefiii_rps * v_e_mps(1), 2.0 * _w_irsecefecefiii_rps * v_e_mps(0), 0.};
}
/* computes the Coriolis acceleration [mps2] viewed in ECEF */

Eigen::Vector3d env::geo::compute_coriolis_n(const Eigen::Vector3d& v_n_mps, const ang::so3_tangent& w_ien_rps) const {
    return w_ien_rps().cross(v_n_mps) * 2.0;
}
/* computes the Coriolis acceleration [mps2] viewed in NED */

Eigen::Vector3d env::geo::compute_coriolis_n(const Eigen::Vector3d& v_n_mps, const double& phi_rad) const {
    return this->compute_wien_rps(phi_rad)().cross(v_n_mps) * 2.0;
}
/* computes the Coriolis acceleration [mps2] viewed in NED */

Eigen::Matrix3d env::geo::jac_coriolis_e_wrt_ve(const ang::so3_tangent& w_iee_rps) const {
    return ang::tools::skew3(w_iee_rps()) * 2.0;
}
/* computes the jacobian of the Coriolis acceleration viewed in ECEF with respect to the ground speed viewed in ECEF */

Eigen::Matrix3d env::geo::jac_coriolis_e_wrt_ve() const {
    Eigen::Matrix3d res = Eigen::Matrix3d::Zero();
    res(0,1) = - 2.0 * _w_irsecefecefiii_rps;
    res(1,0) = + 2.0 * _w_irsecefecefiii_rps;
    return res;
}
/* computes the jacobian of the Coriolis acceleration viewed in eCEF with respect to the ground speed viewed in ECEF */

Eigen::Matrix3d env::geo::jac_coriolis_n_wrt_vn(const ang::so3_tangent& w_ien_rps) const {
    return ang::tools::skew3(w_ien_rps()) * 2.0;
}
/* computes the jacobian of the Coriolis acceleration viewed in NED with respect to the ground speed viewed in NED */

Eigen::Matrix3d env::geo::jac_coriolis_n_wrt_vn(const double& phi_rad) const {
    return ang::tools::skew3(this->compute_wien_rps(phi_rad)()) * 2.0;
}
/* computes the jacobian of the Coriolis acceleration viewed in NED with respect to the ground speed viewed in NED */

/* ===== ===== ===== Transport Acceleration ===== ===== ===== */
/* ========================================================== */

Eigen::Vector3d env::geo::compute_transport_e(const double& N_m, const double& h_m, const double& lambda_rad, const double& phi_rad) const {
    double factor = pow(_w_irsecefecefiii_rps,2) * (N_m + h_m);
    double cosphi = cos(phi_rad);
    return {- factor * cosphi * cos(lambda_rad), - factor * cosphi * sin(lambda_rad), 0.};
}
/* computes the transport acceleration [mps2] viewed in ECEF */

Eigen::Vector3d env::geo::compute_transport_n(const double& N_m, const double& h_m, const double& phi_rad) const {
    double factor = pow(_w_irsecefecefiii_rps,2) * (N_m + h_m);
    double cosphi = cos(phi_rad);
    return {factor * sin(phi_rad) * cosphi, 0.0, factor * pow(cosphi,2)};
}
/* computes the transport acceleration [mps2] viewed in NED */

Eigen::Matrix3d env::geo::jac_transport_e_wrt_ve() const {
    return Eigen::Matrix3d::Zero();
}
/* computes the jacobian of the transport acceleration viewed in ECEF with respect to the ground velocity viewed in ECEF */

Eigen::Matrix3d env::geo::jac_transport_n_wrt_vn() const {
    return Eigen::Matrix3d::Zero();
}
/* computes the jacobian of the transport acceleration viewed in NED with respect to the ground velocity viewed in NED */

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// CLASS GEO_MIX
// =============
// =============

const double env::geo_mix::_Re_m = 6356766.0;
/* Earth radius */

env::geo_mix::geo_mix(env::logic::ZONE_ID zone_id)
: geo(zone_id) {
}
/* constructor based on location. Model and real gravitation and magnetic fields coincide. */

env::geo_mix::geo_mix(env::logic::ZONE_ID zone_id, env::logic::REALISM_ID realism_id, const int& seed)
: geo(zone_id, realism_id, seed) {
}
/* constructor based on location, realism enumerator, and realism seed. Seed only applies to realistic (not model) gravitation and magnetism. */

/* ===== ===== ===== Geodetic <--> Geopotential Altitude Conversion ===== ===== ===== */
/* ================================================================================== */

double env::geo_mix::H2h(const double& H_m, const double& phi_rad) const {
    return _Re_m * H_m / (_Re_m - H_m);
}
/* returns geodetic altitude [m] based on geopotential altitude and latitude */

double env::geo_mix::htoH(const double& h_m, const double& phi_rad) const {
    return _Re_m * h_m / (_Re_m + h_m);
}
/* returns geopotential altitude [m] based on geodetic altitude and latitude */

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////












