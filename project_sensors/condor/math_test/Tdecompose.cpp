#include "Tdecompose.h"
#include "math/logic/constant.h"
#include "math/math/linear_fit_lsq.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

math::test::Tdecompose::Tdecompose(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tdecompose::run() {
	::jail::unit_test::run();

    test_PartialPivLU();
    test_FullPivLU();
    test_LDLT();
    test_BDCSVD();
    test_matrix_exp();
    test_linear_fit_lsq();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_PartialPivLU() {
    Eigen::Matrix3d A; A << 1., 2., 3., 4., 2., 6., 3., 8., 5.;
    Eigen::Vector3d b(4., 2., 7.);

    Eigen::PartialPivLU<Eigen::Matrix3d> Oplu(A);

    Eigen::Matrix3d LU  = Oplu.matrixLU();
    Eigen::Matrix3d L   = LU.template triangularView<Eigen::UnitLower>(); // unit lower triangular
    Eigen::Matrix3d U   = LU.template triangularView<Eigen::Upper>(); // upper triangular
    Eigen::Matrix3d PI  = Oplu.permutationP();
    Eigen::Matrix3d P   = PI.inverse(); // permutation matrix

    Eigen::Matrix3d A2 = P * L * U; // manual reconstruction of A
    Eigen::Matrix3d A3 = Oplu.reconstructedMatrix(); // reconstructed A, for debug purposes

    std::cout << "===== PartialPivLU Decomposition ====="  << std::endl;
    std::cout << "======================================"  << std::endl;
    std::cout << "A:     " << std::endl << A  << std::endl << std::endl;
    std::cout << "L:     " << std::endl << L  << std::endl << std::endl;
    std::cout << "U:     " << std::endl << U  << std::endl << std::endl;
    std::cout << "LU:    " << std::endl << LU << std::endl << std::endl;
    std::cout << "P:     " << std::endl << P  << std::endl << std::endl;
    std::cout << "A2:    " << std::endl << A2 << std::endl << std::endl;
    std::cout << "A3:    " << std::endl << A3 << std::endl << std::endl;

    double det           = A.determinant(); // determinant
    Eigen::Matrix3d AI2  = A.inverse(); // inverse
    Eigen::Matrix3d AI   = Oplu.inverse(); // inverse
    Eigen::Matrix3d I1   = A * AI;
    Eigen::Matrix3d I2   = AI * A;

    std::cout << "det:   " << std::endl << det << std::endl << std::endl;
    std::cout << "AI2:   " << std::endl << AI2 << std::endl << std::endl;
    std::cout << "AI:    " << std::endl << AI  << std::endl << std::endl;
    std::cout << "I1:    " << std::endl << I1  << std::endl << std::endl;
    std::cout << "I2:    " << std::endl << I2  << std::endl << std::endl;

    Eigen::Vector3d x  = Oplu.solve(b);
    Eigen::Vector3d b2 = A * x;

    std::cout << "b:     " << std::endl << b.transpose()  << std::endl << std::endl;
    std::cout << "x:     " << std::endl << x.transpose()  << std::endl << std::endl;
    std::cout << "b2:    " << std::endl << b2.transpose() << std::endl << std::endl;

    check("PartialPivLU A==A2 0 - 0 ", A(0,0), A2(0,0), 1e-12);
    check("PartialPivLU A==A2 0 - 1 ", A(0,1), A2(0,1), 1e-12);
    check("PartialPivLU A==A2 0 - 2 ", A(0,2), A2(0,2), 1e-12);
    check("PartialPivLU A==A2 1 - 0 ", A(1,0), A2(1,0), 1e-12);
    check("PartialPivLU A==A2 1 - 1 ", A(1,1), A2(1,1), 1e-12);
    check("PartialPivLU A==A2 1 - 2 ", A(1,2), A2(1,2), 1e-12);
    check("PartialPivLU A==A2 2 - 0 ", A(2,0), A2(2,0), 1e-12);
    check("PartialPivLU A==A2 2 - 1 ", A(2,1), A2(2,1), 1e-12);
    check("PartialPivLU A==A2 2 - 2 ", A(2,2), A2(2,2), 1e-12);

    check("PartialPivLU A==A3 0 - 0 ", A(0,0), A3(0,0), 1e-12);
    check("PartialPivLU A==A3 0 - 1 ", A(0,1), A3(0,1), 1e-12);
    check("PartialPivLU A==A3 0 - 2 ", A(0,2), A3(0,2), 1e-12);
    check("PartialPivLU A==A3 1 - 0 ", A(1,0), A3(1,0), 1e-12);
    check("PartialPivLU A==A3 1 - 1 ", A(1,1), A3(1,1), 1e-12);
    check("PartialPivLU A==A3 1 - 2 ", A(1,2), A3(1,2), 1e-12);
    check("PartialPivLU A==A3 2 - 0 ", A(2,0), A3(2,0), 1e-12);
    check("PartialPivLU A==A3 2 - 1 ", A(2,1), A3(2,1), 1e-12);
    check("PartialPivLU A==A3 2 - 2 ", A(2,2), A3(2,2), 1e-12);

    check("PartialPivLU AI==AI2 0 - 0 ", AI(0,0), AI2(0,0), 1e-12);
    check("PartialPivLU AI==AI2 0 - 1 ", AI(0,1), AI2(0,1), 1e-12);
    check("PartialPivLU AI==AI2 0 - 2 ", AI(0,2), AI2(0,2), 1e-12);
    check("PartialPivLU AI==AI2 1 - 0 ", AI(1,0), AI2(1,0), 1e-12);
    check("PartialPivLU AI==AI2 1 - 1 ", AI(1,1), AI2(1,1), 1e-12);
    check("PartialPivLU AI==AI2 1 - 2 ", AI(1,2), AI2(1,2), 1e-12);
    check("PartialPivLU AI==AI2 2 - 0 ", AI(2,0), AI2(2,0), 1e-12);
    check("PartialPivLU AI==AI2 2 - 1 ", AI(2,1), AI2(2,1), 1e-12);
    check("PartialPivLU AI==AI2 2 - 2 ", AI(2,2), AI2(2,2), 1e-12);

    check("PartialPivLU I1 0 - 0 ", I1(0,0), 1.0, 1e-12);
    check("PartialPivLU I1 0 - 1 ", I1(0,1), 0.0, 1e-12);
    check("PartialPivLU I1 0 - 2 ", I1(0,2), 0.0, 1e-12);
    check("PartialPivLU I1 1 - 0 ", I1(1,0), 0.0, 1e-12);
    check("PartialPivLU I1 1 - 1 ", I1(1,1), 1.0, 1e-12);
    check("PartialPivLU I1 1 - 2 ", I1(1,2), 0.0, 1e-12);
    check("PartialPivLU I1 2 - 0 ", I1(2,0), 0.0, 1e-12);
    check("PartialPivLU I1 2 - 1 ", I1(2,1), 0.0, 1e-12);
    check("PartialPivLU I1 2 - 2 ", I1(2,2), 1.0, 1e-12);

    check("PartialPivLU I2 0 - 0 ", I2(0,0), 1.0, 1e-12);
    check("PartialPivLU I2 0 - 1 ", I2(0,1), 0.0, 1e-12);
    check("PartialPivLU I2 0 - 2 ", I2(0,2), 0.0, 1e-12);
    check("PartialPivLU I2 1 - 0 ", I2(1,0), 0.0, 1e-12);
    check("PartialPivLU I2 1 - 1 ", I2(1,1), 1.0, 1e-12);
    check("PartialPivLU I2 1 - 2 ", I2(1,2), 0.0, 1e-12);
    check("PartialPivLU I2 2 - 0 ", I2(2,0), 0.0, 1e-12);
    check("PartialPivLU I2 2 - 1 ", I2(2,1), 0.0, 1e-12);
    check("PartialPivLU I2 2 - 2 ", I2(2,2), 1.0, 1e-12);

    check("PartialPivLU b==b2 0     ", b(0), b2(0), 1e-12);
    check("PartialPivLU b==b2 1     ", b(1), b2(1), 1e-12);
    check("PartialPivLU b==b2 2     ", b(2), b2(2), 1e-12);
} // closes test_PartialPivLU

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_FullPivLU() {

    Eigen::Matrix3d A; A << 1., 2., 3., 4., 2., 6., 3., 8., 5.;
    Eigen::Vector3d b(4., 2., 7.);

    Eigen::FullPivLU<Eigen::Matrix3d> Oflu(A);

    Eigen::Matrix3d LU  = Oflu.matrixLU();
    Eigen::Matrix3d L   = LU.template triangularView<Eigen::UnitLower>(); // unit lower triangular
    Eigen::Matrix3d U   = LU.template triangularView<Eigen::Upper>(); // upper triangular
    Eigen::Matrix3d PI  = Oflu.permutationP();
    Eigen::Matrix3d P   = PI.inverse(); // permutation matrix
    Eigen::Matrix3d QI  = Oflu.permutationQ();
    Eigen::Matrix3d Q   = QI.inverse(); // permutation matrix

    Eigen::Matrix3d A2 = P * L * U * Q; // manual reconstruction of A
    Eigen::Matrix3d A3 = Oflu.reconstructedMatrix(); // reconstructed A, for debug purposes

    std::cout << "===== FullPivLU Decomposition ====="     << std::endl;
    std::cout << "==================================="     << std::endl;
    std::cout << "A:     " << std::endl << A  << std::endl << std::endl;
    std::cout << "L:     " << std::endl << L  << std::endl << std::endl;
    std::cout << "U:     " << std::endl << U  << std::endl << std::endl;
    std::cout << "LU:    " << std::endl << LU << std::endl << std::endl;
    std::cout << "P:     " << std::endl << P  << std::endl << std::endl;
    std::cout << "Q:     " << std::endl << Q  << std::endl << std::endl;
    std::cout << "A2:    " << std::endl << A2 << std::endl << std::endl;
    std::cout << "A3:    " << std::endl << A3 << std::endl << std::endl;

    double det           = A.determinant(); // determinant
    Eigen::Matrix3d AI2  = A.inverse(); // inverse
    Eigen::Matrix3d AI   = Oflu.inverse(); // inverse
    Eigen::Matrix3d I1   = A * AI;
    Eigen::Matrix3d I2   = AI * A;
    bool flag_invert     = Oflu.isInvertible();
    long rank            = Oflu.rank();

    std::cout << "det:         " << std::endl << det         << std::endl << std::endl;
    std::cout << "AI2:         " << std::endl << AI2         << std::endl << std::endl;
    std::cout << "AI:          " << std::endl << AI          << std::endl << std::endl;
    std::cout << "I1:          " << std::endl << I1          << std::endl << std::endl;
    std::cout << "I2:          " << std::endl << I2          << std::endl << std::endl;
    std::cout << "flag_invert: " << std::endl << flag_invert << std::endl << std::endl;
    std::cout << "rank:        " << std::endl << rank        << std::endl << std::endl;

    Eigen::Vector3d x  = Oflu.solve(b);
    Eigen::Vector3d b2 = A * x;

    std::cout << "b:     " << std::endl << b.transpose()  << std::endl << std::endl;
    std::cout << "x:     " << std::endl << x.transpose()  << std::endl << std::endl;
    std::cout << "b2:    " << std::endl << b2.transpose() << std::endl << std::endl;

    check("FullPivLU A==A2 0 - 0 ", A(0,0), A2(0,0), 1e-12);
    check("FullPivLU A==A2 0 - 1 ", A(0,1), A2(0,1), 1e-12);
    check("FullPivLU A==A2 0 - 2 ", A(0,2), A2(0,2), 1e-12);
    check("FullPivLU A==A2 1 - 0 ", A(1,0), A2(1,0), 1e-12);
    check("FullPivLU A==A2 1 - 1 ", A(1,1), A2(1,1), 1e-12);
    check("FullPivLU A==A2 1 - 2 ", A(1,2), A2(1,2), 1e-12);
    check("FullPivLU A==A2 2 - 0 ", A(2,0), A2(2,0), 1e-12);
    check("FullPivLU A==A2 2 - 1 ", A(2,1), A2(2,1), 1e-12);
    check("FullPivLU A==A2 2 - 2 ", A(2,2), A2(2,2), 1e-12);

    check("FullPivLU A==A3 0 - 0 ", A(0,0), A3(0,0), 1e-12);
    check("FullPivLU A==A3 0 - 1 ", A(0,1), A3(0,1), 1e-12);
    check("FullPivLU A==A3 0 - 2 ", A(0,2), A3(0,2), 1e-12);
    check("FullPivLU A==A3 1 - 0 ", A(1,0), A3(1,0), 1e-12);
    check("FullPivLU A==A3 1 - 1 ", A(1,1), A3(1,1), 1e-12);
    check("FullPivLU A==A3 1 - 2 ", A(1,2), A3(1,2), 1e-12);
    check("FullPivLU A==A3 2 - 0 ", A(2,0), A3(2,0), 1e-12);
    check("FullPivLU A==A3 2 - 1 ", A(2,1), A3(2,1), 1e-12);
    check("FullPivLU A==A3 2 - 2 ", A(2,2), A3(2,2), 1e-12);

    check("FullPivLU AI==AI2 0 - 0 ", AI(0,0), AI2(0,0), 1e-12);
    check("FullPivLU AI==AI2 0 - 1 ", AI(0,1), AI2(0,1), 1e-12);
    check("FullPivLU AI==AI2 0 - 2 ", AI(0,2), AI2(0,2), 1e-12);
    check("FullPivLU AI==AI2 1 - 0 ", AI(1,0), AI2(1,0), 1e-12);
    check("FullPivLU AI==AI2 1 - 1 ", AI(1,1), AI2(1,1), 1e-12);
    check("FullPivLU AI==AI2 1 - 2 ", AI(1,2), AI2(1,2), 1e-12);
    check("FullPivLU AI==AI2 2 - 0 ", AI(2,0), AI2(2,0), 1e-12);
    check("FullPivLU AI==AI2 2 - 1 ", AI(2,1), AI2(2,1), 1e-12);
    check("FullPivLU AI==AI2 2 - 2 ", AI(2,2), AI2(2,2), 1e-12);

    check("FullPivLU I1 0 - 0 ", I1(0,0), 1.0, 1e-12);
    check("FullPivLU I1 0 - 1 ", I1(0,1), 0.0, 1e-12);
    check("FullPivLU I1 0 - 2 ", I1(0,2), 0.0, 1e-12);
    check("FullPivLU I1 1 - 0 ", I1(1,0), 0.0, 1e-12);
    check("FullPivLU I1 1 - 1 ", I1(1,1), 1.0, 1e-12);
    check("FullPivLU I1 1 - 2 ", I1(1,2), 0.0, 1e-12);
    check("FullPivLU I1 2 - 0 ", I1(2,0), 0.0, 1e-12);
    check("FullPivLU I1 2 - 1 ", I1(2,1), 0.0, 1e-12);
    check("FullPivLU I1 2 - 2 ", I1(2,2), 1.0, 1e-12);

    check("FullPivLU I2 0 - 0 ", I2(0,0), 1.0, 1e-12);
    check("FullPivLU I2 0 - 1 ", I2(0,1), 0.0, 1e-12);
    check("FullPivLU I2 0 - 2 ", I2(0,2), 0.0, 1e-12);
    check("FullPivLU I2 1 - 0 ", I2(1,0), 0.0, 1e-12);
    check("FullPivLU I2 1 - 1 ", I2(1,1), 1.0, 1e-12);
    check("FullPivLU I2 1 - 2 ", I2(1,2), 0.0, 1e-12);
    check("FullPivLU I2 2 - 0 ", I2(2,0), 0.0, 1e-12);
    check("FullPivLU I2 2 - 1 ", I2(2,1), 0.0, 1e-12);
    check("FullPivLU I2 2 - 2 ", I2(2,2), 1.0, 1e-12);

    check("FullPivLU b==b2 0     ", b(0), b2(0), 1e-12);
    check("FullPivLU b==b2 1     ", b(1), b2(1), 1e-12);
    check("FullPivLU b==b2 2     ", b(2), b2(2), 1e-12);

} // closes test_FullPivLU

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_LDLT() {
    Eigen::Matrix3d A; A << 3., 0., 0., 0., 4., 0., 0., 0., 5.;
    Eigen::Vector3d b(4., 2., 7.);

    Eigen::LDLT<Eigen::Matrix3d> Oldlt(A); // for LDLT A needs to be positive or negative semidefinite
    Eigen::ComputationInfo Ainfo = Oldlt.info(); // 0 if OK, 1 if not OK

    Eigen::Matrix3d I;  I.setIdentity();

    Eigen::Matrix3d L  = Oldlt.matrixL(); // lower triangular with ones in diagonal
    Eigen::Matrix3d U  = Oldlt.matrixU(); // transpose of L
    Eigen::Matrix3d D  = Oldlt.vectorD().real().asDiagonal(); // diagonal matrix D
    Eigen::Matrix3d P  = Oldlt.transpositionsP() * I; // permutations matrix P
    Eigen::Matrix3d PT = Oldlt.transpositionsP().transpose() * I; // permutations matrix P transpose

    Eigen::Matrix3d LDU  = L * D * U; // manual reconstruction of LDLT
    Eigen::Matrix3d A2 = PT * LDU * P; // manual reconstruction of A
    Eigen::Matrix3d A3 = Oldlt.reconstructedMatrix(); // reconstructed A, for debug purposes

    Eigen::Vector3d x  = Oldlt.solve(b);
    Eigen::Vector3d b2 = A * x;

    std::cout << "===== LDLT Decomposition =====" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << "A:     " << std::endl << A              << std::endl << std::endl;
    std::cout << "Ainfo: " << std::endl << Ainfo          << std::endl << std::endl;
    std::cout << "L:     " << std::endl << L              << std::endl << std::endl;
    std::cout << "D:     " << std::endl << D              << std::endl << std::endl;
    std::cout << "U:     " << std::endl << U              << std::endl << std::endl;
    std::cout << "P:     " << std::endl << P              << std::endl << std::endl;
    std::cout << "PT:    " << std::endl << PT             << std::endl << std::endl;
    std::cout << "LDU:   " << std::endl << LDU            << std::endl << std::endl;
    std::cout << "A2:    " << std::endl << A2             << std::endl << std::endl;
    std::cout << "A3:    " << std::endl << A3             << std::endl << std::endl;
    std::cout << "b:     " << std::endl << b.transpose()  << std::endl << std::endl;
    std::cout << "x:     " << std::endl << x.transpose()  << std::endl << std::endl;
    std::cout << "b2:    " << std::endl << b2.transpose() << std::endl << std::endl;

    check("LDLT A==A2 0 - 0 ", A(0,0), A2(0,0), 1e-12);
    check("LDLT A==A2 0 - 1 ", A(0,1), A2(0,1), 1e-12);
    check("LDLT A==A2 0 - 2 ", A(0,2), A2(0,2), 1e-12);
    check("LDLT A==A2 1 - 0 ", A(1,0), A2(1,0), 1e-12);
    check("LDLT A==A2 1 - 1 ", A(1,1), A2(1,1), 1e-12);
    check("LDLT A==A2 1 - 2 ", A(1,2), A2(1,2), 1e-12);
    check("LDLT A==A2 2 - 0 ", A(2,0), A2(2,0), 1e-12);
    check("LDLT A==A2 2 - 1 ", A(2,1), A2(2,1), 1e-12);
    check("LDLT A==A2 2 - 2 ", A(2,2), A2(2,2), 1e-12);

    check("LDLT A==A3 0 - 0 ", A(0,0), A3(0,0), 1e-12);
    check("LDLT A==A3 0 - 1 ", A(0,1), A3(0,1), 1e-12);
    check("LDLT A==A3 0 - 2 ", A(0,2), A3(0,2), 1e-12);
    check("LDLT A==A3 1 - 0 ", A(1,0), A3(1,0), 1e-12);
    check("LDLT A==A3 1 - 1 ", A(1,1), A3(1,1), 1e-12);
    check("LDLT A==A3 1 - 2 ", A(1,2), A3(1,2), 1e-12);
    check("LDLT A==A3 2 - 0 ", A(2,0), A3(2,0), 1e-12);
    check("LDLT A==A3 2 - 1 ", A(2,1), A3(2,1), 1e-12);
    check("LDLT A==A3 2 - 2 ", A(2,2), A3(2,2), 1e-12);

    check("LDLT b==b2 0     ", b(0), b2(0), 1e-12);
    check("LDLT b==b2 1     ", b(1), b2(1), 1e-12);
    check("LDLT b==b2 2     ", b(2), b2(2), 1e-12);
} // closes test_LDLT

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_BDCSVD() {
    std::cout << "===== BDCSVD Decomposition =====" << std::endl;
    std::cout << "================================" << std::endl;

    /* Executing solve(), this is, solving least squares, requires Unitary matrixes to
     * be computed, although thin ones suffice. Thin ones are also available for
     * dynamic size matrices, use Full Unitaries in case of fixed size. Thin Unitaries
     * is MUCH QUICKER.
     */

    Eigen::MatrixXf A1  = Eigen::MatrixXf::Random(3,2);
    Eigen::VectorXf b1  = Eigen::VectorXf::Random(3);
    Eigen::Vector2f x1a = A1.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b1);
    Eigen::Vector2f x1b = A1.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b1);
    //Eigen::Vector2f x1e = A1.bdcSvd().solve(b1); // DOES NOT WORK: Unitary matrices required
    Eigen::BDCSVD<Eigen::MatrixXf> Obdcsvd1c(A1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::BDCSVD<Eigen::MatrixXf> Obdcsvd1d(A1, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Eigen::BDCSVD<Eigen::MatrixXf> Obdcsvd1f(A1); // DOES NOT WORK: Unitary matrices required
    Eigen::Vector2f x1c = Obdcsvd1c.solve(b1);
    Eigen::Vector2f x1d = Obdcsvd1d.solve(b1);
    //Eigen::Vector2f x1f = Obdcsvd1f.solve(b1);

    std::cout << "A1:    " << std::endl << A1              << std::endl << std::endl;
    std::cout << "b1:    " << b1.transpose()  << std::endl << std::endl;
    std::cout << "x1a:   " << x1a.transpose() << std::endl << std::endl;
    std::cout << "x1b:   " << x1b.transpose() << std::endl << std::endl;
    std::cout << "x1c:   " << x1c.transpose() << std::endl << std::endl;
    std::cout << "x1d:   " << x1d.transpose() << std::endl << std::endl;

    // ThinU and ThinV only available with dynamic sizes
    Eigen::Matrix<float,3,2> A2; A2 << 1, 2, 3, 1.5, 1, 2.7;
    Eigen::Vector3f          b2; b2 << 1.2, 1.4, 1.6;
    //Eigen::Vector2f x2a = A2.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b2); // DOES NOT WORK: Full Unitary for fixed size matrices
    Eigen::Vector2f x2b = A2.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b2);
    //Eigen::Vector2f x2e = A2.bdcSvd().solve(b2); // DOES NOT WORK: Unitary matrices required
    //Eigen::BDCSVD<Eigen::Matrix<float,3,2>> Obdcsvd2c(A2, Eigen::ComputeThinU | Eigen::ComputeThinV);  // DOES NOT WORK: Full Unitary for fixed size matrices
    Eigen::BDCSVD<Eigen::Matrix<float,3,2>> Obdcsvd2d(A2, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Eigen::BDCSVD<Eigen::Matrix<float,3,2>> Obdcsvd2f(A2); // DOES NOT WORK: Unitary matrices required
    Eigen::Vector2f x2d = Obdcsvd2d.solve(b2);
    //Eigen::Vector2f x2f = Obdcsvd2f.solve(b1);

    std::cout << "A2:    " << A2              << std::endl << std::endl;
    std::cout << "b2:    " << b2.transpose()  << std::endl << std::endl;
    std::cout << "x2b:   " << x2b.transpose() << std::endl << std::endl;
    std::cout << "x2d:   " << x2d.transpose() << std::endl << std::endl;

    Eigen::Matrix<double,10,3> A;
    A << 1, 2, 3, 2, 2, 3, 1, 3, 3, 1, 2, 4, 7, 2, 4, 7, 4, 4, 5, 4, 1, 1, 7, 2, 1, 7, 5, 6, 7, 1;
    Eigen::Matrix<double,10,1> b;
    b << 32.3, 35.9, 37.2, 37.8, 62.3, 71.9, 46.1, 50.8, 69.2, 65.2;

    Eigen::BDCSVD<Eigen::Matrix<double,10,3>> Obdcsvd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Vector3d x  = Obdcsvd.solve(b);

    std::cout << "===== BDCSVD Decomposition =====" << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << "A:     " << std::endl << A              << std::endl << std::endl;
    std::cout << "b:     " << std::endl << b.transpose()  << std::endl << std::endl;
    std::cout << "x:     " << std::endl << x.transpose()  << std::endl << std::endl;
} // closes test_BDCSVD

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_matrix_exp() {
    double pi = math::constant::PI();
    Eigen::Matrix3d A;
    A << 0, -pi/4, 0, pi/4, 0, 0, 0, 0, 0;

    Eigen::Matrix3d Aexp = A.exp();
    Eigen::Matrix3d A2 = Aexp.log();

    std::cout << "===== Matrix Exponential =====" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << "A:    " << std::endl << A    << std::endl << std::endl;
    std::cout << "Aexp: " << std::endl << Aexp << std::endl << std::endl;
    std::cout << "A2:   " << std::endl << A2   << std::endl << std::endl;

    check("Matrix exp and log 0 - 0 ", A(0,0), A2(0,0), 1e-12);
    check("Matrix exp and log 0 - 1 ", A(0,1), A2(0,1), 1e-12);
    check("Matrix exp and log 0 - 2 ", A(0,2), A2(0,2), 1e-12);
    check("Matrix exp and log 1 - 0 ", A(1,0), A2(1,0), 1e-12);
    check("Matrix exp and log 1 - 1 ", A(1,1), A2(1,1), 1e-12);
    check("Matrix exp and log 1 - 2 ", A(1,2), A2(1,2), 1e-12);
    check("Matrix exp and log 2 - 0 ", A(2,0), A2(2,0), 1e-12);
    check("Matrix exp and log 2 - 1 ", A(2,1), A2(2,1), 1e-12);
    check("Matrix exp and log 2 - 2 ", A(2,2), A2(2,2), 1e-12);

    Eigen::Array33d Bexp = A.array().exp();
    Eigen::Array33d A3   = Bexp.log();

    std::cout << "===== Matrix Coefficient Exponential =====" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Aexp: " << std::endl << Bexp << std::endl << std::endl;
    std::cout << "A3:   " << std::endl << A3   << std::endl << std::endl;

    check("Matrix exp and log 0 - 0 ", A(0,0), A3(0,0), 1e-12);
    check("Matrix exp and log 0 - 1 ", A(0,1), A3(0,1), 1e-12);
    check("Matrix exp and log 0 - 2 ", A(0,2), A3(0,2), 1e-12);
    check("Matrix exp and log 1 - 0 ", A(1,0), A3(1,0), 1e-12);
    check("Matrix exp and log 1 - 1 ", A(1,1), A3(1,1), 1e-12);
    check("Matrix exp and log 1 - 2 ", A(1,2), A3(1,2), 1e-12);
    check("Matrix exp and log 2 - 0 ", A(2,0), A3(2,0), 1e-12);
    check("Matrix exp and log 2 - 1 ", A(2,1), A3(2,1), 1e-12);
    check("Matrix exp and log 2 - 2 ", A(2,2), A3(2,2), 1e-12);
} // closes test_matrix_exp

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdecompose::test_linear_fit_lsq() {

    Eigen::MatrixXd in(11,1);
    in << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    Eigen::MatrixXd out(11,1);
    out << 2.1, 2.4, 1.7, 1.9, 3.0, 2.9, 0.8, 1.1, 1.7, 1.8, 2.3;

    double err_mean, err_std;
    Eigen::VectorXd err(11);
    Eigen::Vector2d x;

    math::linear_fit_lsq::compute(x, err_mean, err_std, err, in, out);


    } // closes test_linear_fit_lsq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////