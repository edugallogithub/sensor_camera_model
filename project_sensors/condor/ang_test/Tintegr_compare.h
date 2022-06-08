#ifndef ATT_TEST_INTEGR_COMPARE
#define ATT_TEST_INTEGR_COMPARE

#include "ang_test.h"
#include "ang/quat.h"
#include "ang/rotate/so3_tangent.h"
#include "ang/rotate/rodrigues.h"
#include "ang/transform/speu_rodrigues.h"
#include <jail/unit_test.h>
#include <Eigen/Core>
#include <vector>

/*
 * This text compares the similarity for the different representations of integrating
 * in Euclidean space and then normalizing versus integrating in the tangent space by
 * means of the plus operator.
 * It also checks the similarity of integrating in the global or tangent spaces.
 */

namespace ang {
namespace test {

class Tintegr_compare: public ::jail::unit_test {
public:
    /**< constructor based on counter */
    explicit Tintegr_compare(jail::counter&);
    /**< execute tests and write results on console */
    void run() override;

    /**< specific tests */
    void test_so3_local();
    void test_so3_global();
    void test_so3_local_vs_global();
    void test_se3_local();
    void test_se3_global();
    void test_se3_local_vs_global();

};

} // closes namespace test
} // closes namespace ang

#endif

