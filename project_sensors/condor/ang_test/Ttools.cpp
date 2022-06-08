#include "Ttools.h"

#include "ang/auxiliary.h"
#include "ang/rotate/euler.h"
#include "ang/rotate/rodrigues.h"
#include "ang/tools.h"

#include <iostream>

ang::test::Ttools::Ttools(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Ttools::run() {
	::jail::unit_test::run();

    test_skew();            std::cout << std::endl << std::endl;

    finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Ttools::test_skew() {
    std::cout << "===== ===== Test SKEW ===== ===== " << std::endl;
    std::cout << "===== ===== ========= ===== ===== " << std::endl;

    Eigen::Vector3d a(1.0, 2.0, 3.0);
    Eigen::Vector3d b(7.0, 5.0, 17.0);
    Eigen::Vector3d c = a.cross(b);
    Eigen::Vector3d c2 = ang::tools::skew3(a) * b;
    Eigen::Vector3d c3 = ang::tools::right_skew3(b) * a;
    Eigen::Vector3d c4 = - b.cross(a);
    Eigen::Vector3d c5 = - ang::tools::skew3(b) * a;
    Eigen::Vector3d a2 = ang::tools::skew3_inverse(ang::tools::skew3(a));
    Eigen::Vector3d b2 = ang::tools::right_skew3_inverse(ang::tools::right_skew3(b));

    check("skew3 x:                ", c(0), c2(0), 1e-12);
    check("skew3 y:                ", c(1), c2(1), 1e-12);
    check("skew3 z:                ", c(2), c2(2), 1e-12);
    check("right_skew3 x:          ", c(0), c3(0), 1e-12);
    check("right_skew3 y:          ", c(1), c3(1), 1e-12);
    check("right_skew3 z:          ", c(2), c3(2), 1e-12);
    check("skew3 x:                ", c(0), c4(0), 1e-12);
    check("skew3 y:                ", c(1), c4(1), 1e-12);
    check("skew3 z:                ", c(2), c4(2), 1e-12);
    check("right_skew3 x:          ", c(0), c5(0), 1e-12);
    check("right_skew3 y:          ", c(1), c5(1), 1e-12);
    check("right_skew3 z:          ", c(2), c5(2), 1e-12);
    check("skew3_inverse x:        ", a(0), a2(0), 1e-12);
    check("skew3_inverse y:        ", a(1), a2(1), 1e-12);
    check("skew3_inverse z:        ", a(2), a2(2), 1e-12);
    check("right_skew3_inverse x:  ", b(0), b2(0), 1e-12);
    check("right_skew3_inverse y:  ", b(1), b2(1), 1e-12);
    check("right_skew3_inverse z:  ", b(2), b2(2), 1e-12);

    ang::quat A(1.0, 2.0, 3.0, 4.0);
    ang::quat B(7.0, 5.0, 17.0, 51.0);
    ang::quat C = A * B;
    ang::quat C2 = ang::tools::skew4(A) * B;
    ang::quat C3 = ang::tools::right_skew4(B) * A;
    ang::quat A2 = ang::tools::skew4_inverse(ang::tools::skew4(A));
    ang::quat B2 = ang::tools::right_skew4_inverse(ang::tools::right_skew4(B));

    check("skew4 0:                ", C(0), C2(0), 1e-12);
    check("skew4 x:                ", C(1), C2(1), 1e-12);
    check("skew4 y:                ", C(2), C2(2), 1e-12);
    check("skew4 z:                ", C(3), C2(3), 1e-12);
    check("right_skew4 0:          ", C(0), C3(0), 1e-12);
    check("right_skew4 x:          ", C(1), C3(1), 1e-12);
    check("right_skew4 y:          ", C(2), C3(2), 1e-12);
    check("right_skew4 z:          ", C(3), C3(3), 1e-12);
    check("skew4_inverse 0:        ", A(0), A2(0), 1e-12);
    check("skew4_inverse x:        ", A(1), A2(1), 1e-12);
    check("skew4_inverse y:        ", A(2), A2(2), 1e-12);
    check("skew4_inverse z:        ", A(3), A2(3), 1e-12);
    check("right_skew4_inverse 0:  ", B(0), B2(0), 1e-12);
    check("right_skew4_inverse x:  ", B(1), B2(1), 1e-12);
    check("right_skew4_inverse y:  ", B(2), B2(2), 1e-12);
    check("right_skew4_inverse z:  ", B(3), B2(3), 1e-12);

    Eigen::Vector3d V(2.0, -31.0, 47.0);
    ang::quat QV; QV << 0.0, V;
    ang::quat D = A * QV;
    ang::quat D2 = ang::tools::skew43(A) * V;
    ang::quat D3 = ang::tools::right_skew43(V) * A;
    ang::quat D4 = ang::tools::skew4(A) * QV;
    ang::quat D5 = ang::tools::right_skew4(QV) * A;

    check("skew43 0:                ", D(0), D2(0), 1e-12);
    check("skew43 x:                ", D(1), D2(1), 1e-12);
    check("skew43 y:                ", D(2), D2(2), 1e-12);
    check("skew43 z:                ", D(3), D2(3), 1e-12);
    check("right_skew43 0:          ", D(0), D3(0), 1e-12);
    check("right_skew43 x:          ", D(1), D3(1), 1e-12);
    check("right_skew43 y:          ", D(2), D3(2), 1e-12);
    check("right_skew43 z:          ", D(3), D3(3), 1e-12);
    check("skew4 0:                 ", D(0), D4(0), 1e-12);
    check("skew4 x:                 ", D(1), D4(1), 1e-12);
    check("skew4 y:                 ", D(2), D4(2), 1e-12);
    check("skew4 z:                 ", D(3), D4(3), 1e-12);
    check("right_skew4 0:           ", D(0), D5(0), 1e-12);
    check("right_skew4 x:           ", D(1), D5(1), 1e-12);
    check("right_skew4 y:           ", D(2), D5(2), 1e-12);
    check("right_skew4 z:           ", D(3), D5(3), 1e-12);

} // closes test_skew

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////




































