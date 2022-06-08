#include "error_gen_triple.h"
#include "math/logic/constant.h"
#include <iostream>

// CLASS ERROR_GENERATOR_TRIPLE
// ============================
// ============================

st::error_gen_triple::error_gen_triple(const double& sigma, const int& seed)
: _sigma(sigma), _gen(seed), _dist(0.,1.) {
    _res = _sigma * Eigen::Vector3d(_dist(_gen), _dist(_gen),_dist(_gen));
}
/* constructor based on white noise and seed */

const Eigen::Vector3d& st::error_gen_triple::eval() const {
    return _res;
}
/* return new sensor measurement */

st::error_gen_triple* st::error_gen_triple::create_init_acc_error_generator(st::logic::INITACC_ID initacc_id, const int& seed) {
    st::error_gen_triple* Pres = nullptr;
    switch(initacc_id) {
        case st::logic::initacc_id_zero:
            Pres = new st::error_gen_triple(0., seed);
            break;
        case st::logic::initacc_id_base:
            Pres = new st::error_gen_triple(0.01, seed); // 1% in each component
            break;
        case st::logic::initacc_id_worse:
            Pres = new st::error_gen_triple(0.02, seed); // 2% in each component
            break;
        case st::logic::initacc_id_worst:
            Pres = new st::error_gen_triple(0.05, seed); // 5% in each component
            break;
        case st::logic::initacc_size:
            throw std::runtime_error("Initial accelerometer error generator model not available");
        default:
            throw std::runtime_error("Initial accelerometer error generator model not available");
    }
    return Pres;
}
/* create initial accelerometer error generator */

st::error_gen_triple* st::error_gen_triple::create_init_gyr_error_generator(st::logic::INITGYR_ID initgyr_id, const int& seed) {
    st::error_gen_triple* Pres = nullptr;
    switch(initgyr_id) {
        case st::logic::initgyr_id_zero:
            Pres = new st::error_gen_triple(0., seed);
            break;
        case st::logic::initgyr_id_base:
            Pres = new st::error_gen_triple(0.01, seed); // 1% in each component
            break;
        case st::logic::initgyr_id_worse:
            Pres = new st::error_gen_triple(0.02, seed); // 2% in each component
            break;
        case st::logic::initgyr_id_worst:
            Pres = new st::error_gen_triple(0.05, seed); // 5% in each component
            break;
        case st::logic::initgyr_size:
            throw std::runtime_error("Initial gyroscope error generator model not available");
        default:
            throw std::runtime_error("Initial gyroscope error generator model not available");
    }
    return Pres;
}
/* create initial gyroscope error generator */

st::error_gen_triple* st::error_gen_triple::create_init_mag_error_generator(st::logic::INITMAG_ID initmag_id, const int& seed) {
    st::error_gen_triple* Pres = nullptr;
    switch(initmag_id) {
        case st::logic::initmag_id_zero:
            Pres = new st::error_gen_triple(0., seed);
            break;
        case st::logic::initmag_id_base:
            Pres = new st::error_gen_triple(0.01, seed); // 1% in each component
            break;
        case st::logic::initmag_id_worse:
            Pres = new st::error_gen_triple(0.02, seed); // 2% in each component
            break;
        case st::logic::initmag_id_worst:
            Pres = new st::error_gen_triple(0.05, seed); // 5% in each component
            break;
        case st::logic::initmag_size:
            throw std::runtime_error("Initial magnetometer error generator model not available");
        default:
            throw std::runtime_error("Initial magnetometer error generator model not available");
    }
    return Pres;
}
/* create initial magnetometer error generator */

st::error_gen_triple* st::error_gen_triple::create_init_mgn_error_generator(st::logic::INITMGN_ID initmgn_id, const int& seed) {
    st::error_gen_triple* Pres = nullptr;
    switch(initmgn_id) {
        case st::logic::initmgn_id_zero:
            Pres = new st::error_gen_triple(0., seed);
            break;
        case st::logic::initmgn_id_base:
            Pres = new st::error_gen_triple(0.1, seed); // 10% of difference between real and model magnetic field in each component
            break;
        case st::logic::initmgn_id_worse:
            Pres = new st::error_gen_triple(0.2, seed); // 20% of difference between real and model magnetic field in each component
            break;
        case st::logic::initmgn_id_worst:
            Pres = new st::error_gen_triple(0.5, seed); // 50% of difference between real and model magnetic field in each component
            break;
        case st::logic::initmgn_size:
            throw std::runtime_error("Initial NED magnetic field error generator model not available");
        default:
            throw std::runtime_error("Initial NED magnetic field error generator model not available");
    }
    return Pres;
}
/* create initial NED magnetic field error generator */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////










