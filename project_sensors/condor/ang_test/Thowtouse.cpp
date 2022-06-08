#include "Thowtouse.h"

#include "../ang/rotate/dcm.h"
#include "../ang/rotate/rodrigues.h"
#include "../ang/rotate/euler.h"
#include "../ang/rotate/rotv.h"
#include "../ang/transform/speu_rodrigues.h"
#include "../ang/transform/trfv.h"
#include <iostream>

ang::test::Thowtouse::Thowtouse(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void ang::test::Thowtouse::run() {
	::jail::unit_test::run();

    test_se3_practical();               std::cout << std::endl << std::endl;

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void ang::test::Thowtouse::test_se3_practical() {
    double d2r = math::constant::D2R();

    // The input frame (WRS) has orientation North(1), East(2), Down(3).
    // Let's assume that the aircraft is moving northwards.
    // The camera is originally (CRS1) oriented pointing down, x pointing rightwards, and y facing backwards,
    // this is, East(1), South(2), Down(3). The Euler angles are hence (90, 0, 0), only yaw.
    // The camera is originally (CRS1) located at an altitude of 150 [m], 1 [m] to the
    // North, and 25 [m] to the East --> x_wc1w =(1,25,-150). The translation vector (vector from CRS1
    // origin to WRS origin expressed in CRS1) x_c1wc1 = (-25, 1, 150), or x_wc1w = (1,25,-150).
    //
    // The input (WRS) vector of (5,8,-12) gets thus transformed (CRS1) into (-17 , -4, 138).
    //
    // -> If the camera has no rotation and a translation (vector from CRS2 origin to CRS1 origin expressed in CRS2 origin) of
    //    (-10, -3, -40), the resulting frame (CRS2) has the same Euler angles (90, 0, 0) and translation (vector from
    //    CRS2 origin to WRS origin expressed in CRS2) of (-35, -2, 110).
    //
    // -> If the camera has no translation and a rotation (Euler angles in CRS1) of (0, 10, 0), only pitch,
    //    the resulting frame (CRS3) has the Euler angles (WRS to CRS3) of (90, 10, 0)
    //    and the same translation (vector from CRS3 origin to WRS origin expressed in CRS3) of  (-50.66, 1, 143.38)
    //
    // -> If the camera has a translation (vector from CRS4 origin to CRS1 origin expressed in CRS4 origin) of
    //    (-2.90215, -3, -41.1288) and a rotation (Euler angles in CRS1) of (0, 10, 0), only pitch, the resulting frame (CRS4)
    //    has the Euler angles (WRS to CRS4) of (90, 10, 0) and translation (vector from CRS4 origin to WRS origin
    //    expressed in CRS4) of (-53.5696, -2, 102.251).

    // inputs
    ang::rodrigues q_wc1(ang::euler(90 * d2r, 0 * d2r, 0 * d2r));       // camera rotation based on Euler angles (WRS to CRS1)
    ang::rodrigues q_c1w = q_wc1.inverse();                 // rotation CRS to WRS
    Eigen::Vector3d x_c1wc1(-25.0, 1.0, 150.0);       // camera translation (CRS1 to WRS in CRS1)
    Eigen::Vector3d x_wc1w = - (q_wc1 * x_c1wc1);     // camera translation (WRS to CRS1 in WRS)

    ang::rodrigues q_c1c2(ang::euler(0 * d2r, 0 * d2r, 0 * d2r));       // camera rotation based on Euler angles (CRS1 to CRS2)
    ang::rodrigues q_c2c1 = q_c1c2.inverse();               // rotation CRS2 to CRS1
    Eigen::Vector3d x_c2c1c2(-10., -3., -40.);        // camera translation (CRS2 to CRS1 in CRS2)
    Eigen::Vector3d x_c1c2c1 = - (q_c1c2 * x_c2c1c2); // camera translation (CRS1 to CRS2 in CRS1)

    ang::rodrigues q_c1c3(ang::euler(0 * d2r, 10 * d2r, 0 * d2r));      // camera rotation based on Euler angles (CRS1 to CRS3)
    ang::rodrigues q_c3c1 = q_c1c3.inverse();               // rotation CRS3 to CRS1
    Eigen::Vector3d x_c3c1c3(0., 0., 0.);             // camera translation (CRS3 to CRS1 in CRS3)
    Eigen::Vector3d x_c1c3c1 = - (q_c1c3 * x_c3c1c3); // camera translation (CRS1 to CRS3 in CRS1)

    ang::rodrigues q_c1c4(ang::euler(0 * d2r, 10 * d2r, 0 * d2r));      // camera rotation based on Euler angles (CRS1 to CRS4)
    ang::rodrigues q_c4c1 = q_c1c4.inverse();               // rotation CRS4 to CRS1
    Eigen::Vector3d x_c4c1c4(-2.90215, -3, -41.1288); // camera translation (CRS4 to CRS1 in CRS4)
    Eigen::Vector3d x_c1c4c1 = - (q_c1c4 * x_c4c1c4); // camera translation (CRS1 to CRS4 in CRS1)

    // transformations (rotation plus translation)
    ang::speu_rodrigues Gq_wc1(q_wc1,   x_wc1w);        ang::speu_rodrigues Gq_c1w(q_c1w,   x_c1wc1);
    ang::speu_rodrigues Gq_c1c2(q_c1c2, x_c1c2c1);      ang::speu_rodrigues Gq_c2c1(q_c2c1, x_c2c1c2);
    ang::speu_rodrigues Gq_c1c3(q_c1c3, x_c1c3c1);      ang::speu_rodrigues Gq_c3c1(q_c3c1, x_c3c1c3);
    ang::speu_rodrigues Gq_c1c4(q_c1c4, x_c1c4c1);      ang::speu_rodrigues Gq_c4c1(q_c4c1, x_c4c1c4);

    std::cout << "Gq_wc1 Euler angles:        " << ang::euler(Gq_wc1.get_rodrigues()) << std::endl;
    std::cout << "Gq_wc1 translation (WRS):   " << Gq_wc1.get_T() << std::endl;
    std::cout << "Gq_c1w translation (CRS1):  " << Gq_wc1.get_inverse_T() << std::endl;
    std::cout << "Gq_c1w Euler angles:        " << ang::euler(Gq_c1w.get_rodrigues()) << std::endl;
    std::cout << "Gq_c1w translation (CRS1):  " << Gq_c1w.get_T() << std::endl;
    std::cout << "Gq_wc1 translation (WRS):   " << Gq_c1w.get_inverse_T() << std::endl << std::endl;

    std::cout << "Gq_c1c2 Euler angles:       " << ang::euler(Gq_c1c2.get_rodrigues()) << std::endl;
    std::cout << "Gq_c1c2 translation (CRS1): " << Gq_c1c2.get_T() << std::endl;
    std::cout << "Gq_c2c1 translation (CRS2): " << Gq_c1c2.get_inverse_T() << std::endl;
    std::cout << "Gq_c2c1 Euler angles:       " << ang::euler(Gq_c2c1.get_rodrigues()) << std::endl;
    std::cout << "Gq_c2c1 translation (CRS2): " << Gq_c2c1.get_T() << std::endl;
    std::cout << "Gq_c1c2 translation (CRS1): " << Gq_c2c1.get_inverse_T() << std::endl << std::endl;

    std::cout << "Gq_c1c3 Euler angles:       " << ang::euler(Gq_c1c3.get_rodrigues()) << std::endl;
    std::cout << "Gq_c1c3 translation (CRS1): " << Gq_c1c3.get_T() << std::endl;
    std::cout << "Gq_c3c1 translation (CRS3): " << Gq_c1c3.get_inverse_T() << std::endl;
    std::cout << "Gq_c3c1 Euler angles:       " << ang::euler(Gq_c3c1.get_rodrigues()) << std::endl;
    std::cout << "Gq_c3c1 translation (CRS3): " << Gq_c3c1.get_T() << std::endl;
    std::cout << "Gq_c1c3 translation (CRS1): " << Gq_c3c1.get_inverse_T() << std::endl << std::endl;

    std::cout << "Gq_c1c4 Euler angles:       " << ang::euler(Gq_c1c4.get_rodrigues()) << std::endl;
    std::cout << "Gq_c1c4 translation (CRS1): " << Gq_c1c4.get_T() << std::endl;
    std::cout << "Gq_c4c1 translation (CRS4): " << Gq_c1c4.get_inverse_T() << std::endl;
    std::cout << "Gq_c4c1 Euler angles:       " << ang::euler(Gq_c4c1.get_rodrigues()) << std::endl;
    std::cout << "Gq_c4c1 translation (CRS4): " << Gq_c4c1.get_T() << std::endl;
    std::cout << "Gq_c1c4 translation (CRS1): " << Gq_c4c1.get_inverse_T() << std::endl << std::endl;

    // ===== ===== ===== WAYS TO OBTAIN SAME RESULT IN TRANSFORMATION WRS->CRS1 ===== ===== =====
    Eigen::Vector3d v_w(5.0, 8.0, -12.0); // input vector (WRS)
    Eigen::Vector3d Av_c1 = q_c1w * v_w + x_c1wc1;
    Eigen::Vector3d Bv_c1 = Gq_c1w * v_w;
    Eigen::Vector3d Cv_c1 = q_c1w.inverse() / v_w + x_c1wc1;
    Eigen::Vector3d Dv_c1 = Gq_c1w.inverse() / v_w;
    Eigen::Vector3d Ev_c1 = Gq_wc1 / v_w;

    std::cout << "Av_c1:    " << Av_c1 << std::endl;
    std::cout << "Bv_c1:    " << Bv_c1 << std::endl;
    std::cout << "Cv_c1:    " << Cv_c1 << std::endl;
    std::cout << "Dv_c1:    " << Dv_c1 << std::endl;
    std::cout << "Ev_c1:    " << Ev_c1 << std::endl << std::endl;

    // ===== ===== ===== 4 WAYS TO OBTAIN SAME RESULT IN TRANSFORMATION CRS1->WRS ===== ===== =====
    Eigen::Vector3d Av_w = q_c1w / Av_c1 - q_c1w / x_c1wc1;
    Eigen::Vector3d Bv_w = Gq_c1w / Bv_c1;
    Eigen::Vector3d Cv_w = q_c1w.inverse() * Cv_c1 - q_c1w.inverse() * x_c1wc1;
    Eigen::Vector3d Dv_w = Gq_c1w.inverse() * Dv_c1;
    Eigen::Vector3d Ev_w = Gq_wc1 * Ev_c1;
    std::cout << "v_w:     " << v_w << std::endl;
    std::cout << "Av_w:    " << Av_w << std::endl;
    std::cout << "Bv_w:    " << Bv_w << std::endl;
    std::cout << "Cv_w:    " << Cv_w << std::endl;
    std::cout << "Dv_w:    " << Dv_w << std::endl;
    std::cout << "Ev_w:    " << Ev_w << std::endl << std::endl;

    // ===== ===== ===== WAYS TO OBTAIN THE SAME RESULT WHEN CONCATENATING TRANSFORMATIONS ===== ===== =====
    ang::speu_rodrigues Gqa_wc2 = Gq_wc1 * Gq_c1c2;
    ang::speu_rodrigues Gqb_wc2(q_wc1 * q_c1c2, x_wc1w + q_wc1 * x_c1c2c1);
    ang::speu_rodrigues Gqa_c2w = Gq_c2c1 * Gq_c1w;
    ang::speu_rodrigues Gqb_c2w(q_c2c1 * q_c1w, x_c2c1c2 + q_c2c1 * x_c1wc1);

    std::cout << "Translation only -> Ta_wc2 Euler angles:       " << ang::euler(Gqa_wc2.get_rodrigues()) << std::endl;
    std::cout << "Translation only -> Ta_wc2 translation (WRS):  " << Gqa_wc2.get_T() << std::endl;
    std::cout << "Translation only -> Ta_c2w translation (CRS2): " << Gqa_wc2.get_inverse_T() << std::endl;
    std::cout << "Translation only -> Tb_wc2 Euler angles:       " << ang::euler(Gqb_wc2.get_rodrigues()) << std::endl;
    std::cout << "Translation only -> Tb_wc2 translation (WRS):  " << Gqb_wc2.get_T() << std::endl;
    std::cout << "Translation only -> Tb_c2w translation (CRS2): " << Gqb_wc2.get_inverse_T() << std::endl;
    std::cout << "Translation only -> Ta_c2w Euler angles:       " << ang::euler(Gqa_c2w.get_rodrigues()) << std::endl;
    std::cout << "Translation only -> Ta_c2w translation (CRS2): " << Gqa_c2w.get_T() << std::endl;
    std::cout << "Translation only -> Ta_wc2 translation (WRS):  " << Gqa_c2w.get_inverse_T() << std::endl;
    std::cout << "Translation only -> Tb_c2w Euler angles:       " << ang::euler(Gqb_c2w.get_rodrigues()) << std::endl;
    std::cout << "Translation only -> Tb_c2w translation (CRS2): " << Gqb_c2w.get_T() << std::endl;
    std::cout << "Translation only -> Tb_wc2 translation (WRS):  " << Gqb_c2w.get_inverse_T() << std::endl << std::endl;

    // ===== ===== ===== WAYS TO OBTAIN THE SAME RESULT WHEN CONCATENATING TRANSFORMATIONS ===== ===== =====
    ang::speu_rodrigues Gqa_wc3 = Gq_wc1 * Gq_c1c3;
    ang::speu_rodrigues Gqb_wc3(q_wc1 * q_c1c3, x_wc1w + q_wc1 * x_c1c3c1);
    ang::speu_rodrigues Gqa_c3w = Gq_c3c1 * Gq_c1w;
    ang::speu_rodrigues Gqb_c3w(q_c3c1 * q_c1w, x_c3c1c3 + q_c3c1 * x_c1wc1);

    std::cout << "Rotation only -> Ta_wc3 Euler angles:       " << ang::euler(Gqa_wc3.get_rodrigues()) << std::endl;
    std::cout << "Rotation only -> Ta_wc3 translation (WRS):  " << Gqa_wc3.get_T() << std::endl;
    std::cout << "Rotation only -> Ta_c3w translation (CRS3): " << Gqa_wc3.get_inverse_T() << std::endl;
    std::cout << "Rotation only -> Tb_wc3 Euler angles:       " << ang::euler(Gqb_wc3.get_rodrigues()) << std::endl;
    std::cout << "Rotation only -> Tb_wc3 translation (WRS):  " << Gqb_wc3.get_T() << std::endl;
    std::cout << "Rotation only -> Tb_c3w translation (CRS3): " << Gqb_wc3.get_inverse_T() << std::endl;
    std::cout << "Rotation only -> Ta_c3w Euler angles:       " << ang::euler(Gqa_c3w.get_rodrigues()) << std::endl;
    std::cout << "Rotation only -> Ta_c3w translation (CRS3): " << Gqa_c3w.get_T() << std::endl;
    std::cout << "Rotation only -> Ta_wc3 translation (WRS):  " << Gqa_c3w.get_inverse_T() << std::endl;
    std::cout << "Rotation only -> Tb_c3w Euler angles:       " << ang::euler(Gqb_c3w.get_rodrigues()) << std::endl;
    std::cout << "Rotation only -> Tb_c3w translation (CRS3): " << Gqb_c3w.get_T() << std::endl;
    std::cout << "Rotation only -> Tb_wc3 translation (WRS):  " << Gqb_c3w.get_inverse_T() << std::endl << std::endl;

    // ===== ===== ===== WAYS TO OBTAIN THE SAME RESULT WHEN CONCATENATING TRANSFORMATIONS ===== ===== =====
    ang::speu_rodrigues Gqa_wc4 = Gq_wc1 * Gq_c1c4;
    ang::speu_rodrigues Gqb_wc4(q_wc1 * q_c1c4, x_wc1w + q_wc1 * x_c1c4c1);
    ang::speu_rodrigues Gqa_c4w = Gq_c4c1 * Gq_c1w;
    ang::speu_rodrigues Gqb_c4w(q_c4c1 * q_c1w, x_c4c1c4 + q_c4c1 * x_c1wc1);

    std::cout << "Translation and rotation -> Ta_wc4 Euler angles:       " << ang::euler(Gqa_wc4.get_rodrigues()) << std::endl;
    std::cout << "Translation and rotation -> Ta_wc4 translation (WRS):  " << Gqa_wc4.get_T() << std::endl;
    std::cout << "Translation and rotation -> Ta_c4w translation (CRS4): " << Gqa_wc4.get_inverse_T() << std::endl;
    std::cout << "Translation and rotation -> Tb_wc4 Euler angles:       " << ang::euler(Gqb_wc4.get_rodrigues()) << std::endl;
    std::cout << "Translation and rotation -> Tb_wc4 translation (WRS):  " << Gqb_wc4.get_T() << std::endl;
    std::cout << "Translation and rotation -> Tb_c4w translation (CRS4): " << Gqb_wc4.get_inverse_T() << std::endl;
    std::cout << "Translation and rotation -> Ta_c4w Euler angles:       " << ang::euler(Gqa_c4w.get_rodrigues()) << std::endl;
    std::cout << "Translation and rotation -> Ta_c4w translation (CRS4): " << Gqa_c4w.get_T() << std::endl;
    std::cout << "Translation and rotation -> Ta_wc4 translation (WRS):  " << Gqa_c4w.get_inverse_T() << std::endl;
    std::cout << "Translation and rotation -> Tb_c4w Euler angles:       " << ang::euler(Gqb_c4w.get_rodrigues()) << std::endl;
    std::cout << "Translation and rotation -> Tb_c4w translation (CRS4): " << Gqb_c4w.get_T() << std::endl;
    std::cout << "Translation and rotation -> Tb_wc4 translation (WRS):  " << Gqb_c4w.get_inverse_T() << std::endl << std::endl;

    // ===== ===== ===== 2 WAYS TO USE THE TWIST ===== ===== =====
    // The negative of the twist is the same as the twist of the inverse transformation
    ang::trfv taua_c4w(Gqa_c4w);
    std::cout << "Ta_c4w_twist Euler angles:       " << ang::euler(taua_c4w.get_rotv()) << std::endl;
    std::cout << "Ta_c4w_twist translation (CRS1): " << taua_c4w.get_T().transpose() << std::endl;
    std::cout << "Ta_c4w_twist translation (WRS):  " << taua_c4w.get_inverse_T().transpose() << std::endl;
    std::cout << std::endl;

} // closes test_se3_practical_w2c

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////






































