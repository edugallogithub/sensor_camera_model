#ifndef MATH_LOW_PASS_H
#define MATH_LOW_PASS_H

#include "../math.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace math {

/*
 * Implements the discrete time low pass filter by Finn Haugen
 */

// CLASS LOW_PASS_SINGLE
// =====================
// =====================

class MATH_API low_pass_single {
private:
    /**< filter time constant */
    double _Tf_sec;
    /**< filter time sample */
    double _Deltat_sec;
    /**< filter ratio */
    double _alpha;
    /**< filter ratio factor */
    double _beta;
    /**< filter state */
    double _y_prev;
public:
    /**< default constructor */
    low_pass_single() = delete;
    /**< constructor based on filter time constant and filter time sample (constant at least 5 times sample) */
    low_pass_single(const double& Tf_sec, const double& Deltat_sec);
    /**< copy constructor */
    low_pass_single(const low_pass_single&) = default;
    /**< move constructor */
    low_pass_single(low_pass_single&&) = default;
    /**< destructor */
    ~low_pass_single() = default;
    /**< copy assignment */
    low_pass_single& operator=(const low_pass_single&) = default;
    /**< move assignment */
    low_pass_single& operator=(low_pass_single&&) = default;

    /**< update filter characteristics without modifying its state */
    void update(const double& Tf_sec, const double& Deltat_sec);
    /**< initialize filter (first evaluation) */
    void init(const double& u0);
    /**< evaluate filter */
    double eval(const double& u);
    /**< evaluate filter without modifying its state */
    double eval_without_modify(const double& u) const;
    /**< temporary for showing results */
    double XXeval(const double& u);

    /**< get filter time constant to read */
    const double& get_Tf_sec() const {return _Tf_sec;}
    /**< get filter time sample to read */
    const double& get_Deltat_sec() const {return _Deltat_sec;}
    /**< get filter ratio to read */
    const double& get_alpha() const {return _alpha;}
    /**< get filter ratio factor to read */
    const double& get_beta() const {return _beta;}
    /**< get filter state to read */
    const double& get_y_prev() const {return _y_prev;}
}; // closes class low_pass_single

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// CLASS LOW_PASS_TRIPLE
// =====================
// =====================

class MATH_API low_pass_triple {
private:
    /**< filter time constant */
    double _Tf_sec;
    /**< filter time sample */
    double _Deltat_sec;
    /**< filter ratio */
    double _alpha;
    /**< filter ratio factor */
    double _beta;
    /**< filter state */
    Eigen::Vector3d _y_prev;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    low_pass_triple() = delete;
    /**< constructor based on filter time constant and filter time sample (constant at least 5 times sample) */
    low_pass_triple(const double& Tf_sec, const double& Deltat_sec);
    /**< copy constructor */
    low_pass_triple(const low_pass_triple&) = default;
    /**< move constructor */
    low_pass_triple(low_pass_triple&&) = default;
    /**< destructor */
    ~low_pass_triple() = default;
    /**< copy assignment */
    low_pass_triple& operator=(const low_pass_triple&) = default;
    /**< move assignment */
    low_pass_triple& operator=(low_pass_triple&&) = default;

    /**< update filter characteristics without modifying its state */
    void update(const double& Tf_sec, const double& Deltat_sec);
    /**< initialize filter (first evaluation) */
    void init(const Eigen::Vector3d& u0);
    /**< evaluate filter */
    Eigen::Vector3d eval(const Eigen::Vector3d& u);
    /**< evaluate filter without modifying its state */
    Eigen::Vector3d eval_without_modify(const Eigen::Vector3d& u) const;
    /**< temporary for showing results */
    Eigen::Vector3d XXeval(const Eigen::Vector3d& u);

    /**< get filter time constant to read */
    const double& get_Tf_sec() const {return _Tf_sec;}
    /**< get filter time sample to read */
    const double& get_Deltat_sec() const {return _Deltat_sec;}
    /**< get filter ratio to read */
    const double& get_alpha() const {return _alpha;}
    /**< get filter ratio factor to read */
    const double& get_beta() const {return _beta;}
    /**< get filter state to read */
    const Eigen::Vector3d& get_y_prev() const {return _y_prev;}
}; // closes class low_pass_triple

} // closes namespace math

#endif

