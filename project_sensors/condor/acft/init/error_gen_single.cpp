#include "error_gen_single.h"
#include "math/logic/constant.h"
#include <iostream>

// CLASS ERROR_GENERATOR_SINGLE
// ============================
// ============================

st::error_gen_single::error_gen_single(const double& sigma, const int& seed)
: _sigma(sigma), _gen(seed), _dist(0.,1.) {
    _res = _sigma * _dist(_gen);
}
/* constructor based on white noise and seed */

const double& st::error_gen_single::eval() const {
    return _res;
}
/* return new sensor measurement */

st::error_gen_single* st::error_gen_single::create_init_height_error_generator(st::logic::INITHGH_ID initheight_id, const int &seed) {
    st::error_gen_single* Pres = nullptr;
    switch(initheight_id) {
        case st::logic::inithgh_id_zero:
            Pres = new st::error_gen_single(0., seed);
            break;
        case st::logic::inithgh_id_base:
            Pres = new st::error_gen_single(50.0, seed); // 50 [m]
            break;
        case st::logic::inithgh_id_worse:
            Pres = new st::error_gen_single(100.0, seed); // 100 [m]
            break;
        case st::logic::inithgh_size:
            throw std::runtime_error("Initial height error generator model not available");
        default:
            throw std::runtime_error("Initial height error generator model not available");
    }
    return Pres;
}
/* create initial height over terrain error generator */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////










