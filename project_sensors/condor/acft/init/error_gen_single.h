#ifndef ACFT_ERROR_GEN_SINGLE
#define ACFT_ERROR_GEN_SINGLE

#include "../acft.h"
#include "../logic.h"

#include <random>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace st {

// CLASS ERROR_GENERATOR_SINGLE
// ============================
// ============================

class ACFT_API error_gen_single {
private:
    /**< white noise [unit] */
    double _sigma;
    /**< seed generator */
    std::ranlux24_base _gen;
    /**< standard normal distribution */
    std::normal_distribution<double> _dist;

    /**< evaluation result, so it can be retrieved multiple times */
    double _res;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**< default constructor */
    error_gen_single() = delete;
    /**< constructor based on white noise and seed */
    error_gen_single(const double& sigma, const int& seed);
    /**< copy constructor */
    error_gen_single(const error_gen_single&) = delete;
    /**< move constructor */
    error_gen_single(error_gen_single&&) = delete;
    /**< destructor */
    ~error_gen_single() = default;
    /**< copy assignment */
    error_gen_single& operator=(const error_gen_single&) = delete;
    /**< move assignment */
    error_gen_single& operator=(error_gen_single&&) = delete;

    /**< return sensor measurement */
    const double& eval() const;

    /**< create initial height over terrain error generator */
    static st::error_gen_single* create_init_height_error_generator(st::logic::INITHGH_ID, const int& seed);

    /**< get white noise [unit] to read */
    const double& get_sigma() const {return _sigma;}
}; // closes class error_gen_triple

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}; // closes namespace st

#endif
