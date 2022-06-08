#include "Tlow_pass.h"

math::test::Tlow_pass::Tlow_pass(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tlow_pass::run() {
	::jail::unit_test::run();

	test1();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void math::test::Tlow_pass::test1() {

    std::vector<double> u(21);

    u[0] = 0;
    u[1] = 10;
    u[2] = 20;
    u[3] = 30;
    u[4] = 40;
    u[5] = 50;
    u[6] = 60;
    u[7] = 70;
    u[8] = 80;
    u[9] = 90;
    u[10] = 100;
    u[11] = 110;
    u[12] = 120;
    u[13] = 130;
    u[14] = 140;
    u[15] = 150;
    u[16] = 160;
    u[17] = 170;
    u[18] = 180;
    u[19] = 190;
    u[20] = 200;

    math::low_pass_single Tlp0(0., 1.);
    math::low_pass_single Tlp1(1., 1.);
    math::low_pass_single Tlp2(2., 1.);
    math::low_pass_single Tlp4(4., 1.);
    math::low_pass_single Tlp8(8., 1.);
    math::low_pass_single Tlp16(16., 1.);
    math::low_pass_single Tlp32(32., 1.);
    std::vector<double> y0(21);
    std::vector<double> y1(21);
    std::vector<double> y2(21);
    std::vector<double> y4(21);
    std::vector<double> y8(21);
    std::vector<double> y16(21);
    std::vector<double> y32(21);

    Tlp0.init(u[0]);            Tlp1.init(u[0]);            Tlp2.init(u[0]);            Tlp4.init(u[0]);            Tlp8.init(u[0]);            Tlp16.init(u[0]);               Tlp32.init(u[0]);
    y0[0]  = Tlp0.eval(u[0]);   y1[0]  = Tlp1.eval(u[0]);   y2[0]  = Tlp2.eval(u[0]);   y4[0]  = Tlp4.eval(u[0]);   y8[0]  = Tlp8.eval(u[0]);   y16[0]  = Tlp16.eval(u[0]);     y32[0]  = Tlp32.eval(u[0]);
    y0[1]  = Tlp0.eval(u[1]);   y1[1]  = Tlp1.eval(u[1]);   y2[1]  = Tlp2.eval(u[1]);   y4[1]  = Tlp4.eval(u[1]);   y8[1]  = Tlp8.eval(u[1]);   y16[1]  = Tlp16.eval(u[1]);     y32[1]  = Tlp32.eval(u[1]);
    y0[2]  = Tlp0.eval(u[2]);   y1[2]  = Tlp1.eval(u[2]);   y2[2]  = Tlp2.eval(u[2]);   y4[2]  = Tlp4.eval(u[2]);   y8[2]  = Tlp8.eval(u[2]);   y16[2]  = Tlp16.eval(u[2]);     y32[2]  = Tlp32.eval(u[2]);
    y0[3]  = Tlp0.eval(u[3]);   y1[3]  = Tlp1.eval(u[3]);   y2[3]  = Tlp2.eval(u[3]);   y4[3]  = Tlp4.eval(u[3]);   y8[3]  = Tlp8.eval(u[3]);   y16[3]  = Tlp16.eval(u[3]);     y32[3]  = Tlp32.eval(u[3]);
    y0[4]  = Tlp0.eval(u[4]);   y1[4]  = Tlp1.eval(u[4]);   y2[4]  = Tlp2.eval(u[4]);   y4[4]  = Tlp4.eval(u[4]);   y8[4]  = Tlp8.eval(u[4]);   y16[4]  = Tlp16.eval(u[4]);     y32[4]  = Tlp32.eval(u[4]);
    y0[5]  = Tlp0.eval(u[5]);   y1[5]  = Tlp1.eval(u[5]);   y2[5]  = Tlp2.eval(u[5]);   y4[5]  = Tlp4.eval(u[5]);   y8[5]  = Tlp8.eval(u[5]);   y16[5]  = Tlp16.eval(u[5]);     y32[5]  = Tlp32.eval(u[5]);
    y0[6]  = Tlp0.eval(u[6]);   y1[6]  = Tlp1.eval(u[6]);   y2[6]  = Tlp2.eval(u[6]);   y4[6]  = Tlp4.eval(u[6]);   y8[6]  = Tlp8.eval(u[6]);   y16[6]  = Tlp16.eval(u[6]);     y32[6]  = Tlp32.eval(u[6]);
    y0[7]  = Tlp0.eval(u[7]);   y1[7]  = Tlp1.eval(u[7]);   y2[7]  = Tlp2.eval(u[7]);   y4[7]  = Tlp4.eval(u[7]);   y8[7]  = Tlp8.eval(u[7]);   y16[7]  = Tlp16.eval(u[7]);     y32[7]  = Tlp32.eval(u[7]);
    y0[8]  = Tlp0.eval(u[8]);   y1[8]  = Tlp1.eval(u[8]);   y2[8]  = Tlp2.eval(u[8]);   y4[8]  = Tlp4.eval(u[8]);   y8[8]  = Tlp8.eval(u[8]);   y16[8]  = Tlp16.eval(u[8]);     y32[8]  = Tlp32.eval(u[8]);
    y0[9]  = Tlp0.eval(u[9]);   y1[9]  = Tlp1.eval(u[9]);   y2[9]  = Tlp2.eval(u[9]);   y4[9]  = Tlp4.eval(u[9]);   y8[9]  = Tlp8.eval(u[9]);   y16[9]  = Tlp16.eval(u[9]);     y32[9]  = Tlp32.eval(u[9]);
    y0[10] = Tlp0.eval(u[10]);  y1[10] = Tlp1.eval(u[10]);  y2[10] = Tlp2.eval(u[10]);  y4[10] = Tlp4.eval(u[10]);  y8[10] = Tlp8.eval(u[10]);  y16[10] = Tlp16.eval(u[10]);    y32[10] = Tlp32.eval(u[10]);
    y0[11] = Tlp0.eval(u[11]);  y1[11] = Tlp1.eval(u[11]);  y2[11] = Tlp2.eval(u[11]);  y4[11] = Tlp4.eval(u[11]);  y8[11] = Tlp8.eval(u[11]);  y16[11] = Tlp16.eval(u[11]);    y32[11] = Tlp32.eval(u[11]);
    y0[12] = Tlp0.eval(u[12]);  y1[12] = Tlp1.eval(u[12]);  y2[12] = Tlp2.eval(u[12]);  y4[12] = Tlp4.eval(u[12]);  y8[12] = Tlp8.eval(u[12]);  y16[12] = Tlp16.eval(u[12]);    y32[12] = Tlp32.eval(u[12]);
    y0[13] = Tlp0.eval(u[13]);  y1[13] = Tlp1.eval(u[13]);  y2[13] = Tlp2.eval(u[13]);  y4[13] = Tlp4.eval(u[13]);  y8[13] = Tlp8.eval(u[13]);  y16[13] = Tlp16.eval(u[13]);    y32[13] = Tlp32.eval(u[13]);
    y0[14] = Tlp0.eval(u[14]);  y1[14] = Tlp1.eval(u[14]);  y2[14] = Tlp2.eval(u[14]);  y4[14] = Tlp4.eval(u[14]);  y8[14] = Tlp8.eval(u[14]);  y16[14] = Tlp16.eval(u[14]);    y32[14] = Tlp32.eval(u[14]);
    y0[15] = Tlp0.eval(u[15]);  y1[15] = Tlp1.eval(u[15]);  y2[15] = Tlp2.eval(u[15]);  y4[15] = Tlp4.eval(u[15]);  y8[15] = Tlp8.eval(u[15]);  y16[15] = Tlp16.eval(u[15]);    y32[15] = Tlp32.eval(u[15]);
    y0[16] = Tlp0.eval(u[16]);  y1[16] = Tlp1.eval(u[16]);  y2[16] = Tlp2.eval(u[16]);  y4[16] = Tlp4.eval(u[16]);  y8[16] = Tlp8.eval(u[16]);  y16[16] = Tlp16.eval(u[16]);    y32[16] = Tlp32.eval(u[16]);
    y0[17] = Tlp0.eval(u[17]);  y1[17] = Tlp1.eval(u[17]);  y2[17] = Tlp2.eval(u[17]);  y4[17] = Tlp4.eval(u[17]);  y8[17] = Tlp8.eval(u[17]);  y16[17] = Tlp16.eval(u[17]);    y32[17] = Tlp32.eval(u[17]);
    y0[18] = Tlp0.eval(u[18]);  y1[18] = Tlp1.eval(u[18]);  y2[18] = Tlp2.eval(u[18]);  y4[18] = Tlp4.eval(u[18]);  y8[18] = Tlp8.eval(u[18]);  y16[18] = Tlp16.eval(u[18]);    y32[18] = Tlp32.eval(u[18]);
    y0[19] = Tlp0.eval(u[19]);  y1[19] = Tlp1.eval(u[19]);  y2[19] = Tlp2.eval(u[19]);  y4[19] = Tlp4.eval(u[19]);  y8[19] = Tlp8.eval(u[19]);  y16[19] = Tlp16.eval(u[19]);    y32[19] = Tlp32.eval(u[19]);
    y0[20] = Tlp0.eval(u[20]);  y1[20] = Tlp1.eval(u[20]);  y2[20] = Tlp2.eval(u[20]);  y4[20] = Tlp4.eval(u[20]);  y8[20] = Tlp8.eval(u[20]);  y16[20] = Tlp16.eval(u[20]);    y32[20] = Tlp32.eval(u[20]);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::setw(7) << "0:  " << std::setw(9) << u[0]  << std::setw(9) << y0[0]  << std::setw(9) << y1[0]  << std::setw(9) << y2[0]  << std::setw(9) << y4[0]  << std::setw(9) << y8[0]  << std::setw(9) << y16[0]  << std::setw(9) << y32[0]  << std::endl;
    std::cout << std::setw(7) << "1:  " << std::setw(9) << u[1]  << std::setw(9) << y0[1]  << std::setw(9) << y1[1]  << std::setw(9) << y2[1]  << std::setw(9) << y4[1]  << std::setw(9) << y8[1]  << std::setw(9) << y16[1]  << std::setw(9) << y32[1]  << std::endl;
    std::cout << std::setw(7) << "2:  " << std::setw(9) << u[2]  << std::setw(9) << y0[2]  << std::setw(9) << y1[2]  << std::setw(9) << y2[2]  << std::setw(9) << y4[2]  << std::setw(9) << y8[2]  << std::setw(9) << y16[2]  << std::setw(9) << y32[2]  << std::endl;
    std::cout << std::setw(7) << "3:  " << std::setw(9) << u[3]  << std::setw(9) << y0[3]  << std::setw(9) << y1[3]  << std::setw(9) << y2[3]  << std::setw(9) << y4[3]  << std::setw(9) << y8[3]  << std::setw(9) << y16[3]  << std::setw(9) << y32[3]  << std::endl;
    std::cout << std::setw(7) << "4:  " << std::setw(9) << u[4]  << std::setw(9) << y0[4]  << std::setw(9) << y1[4]  << std::setw(9) << y2[4]  << std::setw(9) << y4[4]  << std::setw(9) << y8[4]  << std::setw(9) << y16[4]  << std::setw(9) << y32[4]  << std::endl;
    std::cout << std::setw(7) << "5:  " << std::setw(9) << u[5]  << std::setw(9) << y0[5]  << std::setw(9) << y1[5]  << std::setw(9) << y2[5]  << std::setw(9) << y4[5]  << std::setw(9) << y8[5]  << std::setw(9) << y16[5]  << std::setw(9) << y32[5]  << std::endl;
    std::cout << std::setw(7) << "6:  " << std::setw(9) << u[6]  << std::setw(9) << y0[6]  << std::setw(9) << y1[6]  << std::setw(9) << y2[6]  << std::setw(9) << y4[6]  << std::setw(9) << y8[6]  << std::setw(9) << y16[6]  << std::setw(9) << y32[6]  << std::endl;
    std::cout << std::setw(7) << "7:  " << std::setw(9) << u[7]  << std::setw(9) << y0[7]  << std::setw(9) << y1[7]  << std::setw(9) << y2[7]  << std::setw(9) << y4[7]  << std::setw(9) << y8[7]  << std::setw(9) << y16[7]  << std::setw(9) << y32[7]  << std::endl;
    std::cout << std::setw(7) << "8:  " << std::setw(9) << u[8]  << std::setw(9) << y0[8]  << std::setw(9) << y1[8]  << std::setw(9) << y2[8]  << std::setw(9) << y4[8]  << std::setw(9) << y8[8]  << std::setw(9) << y16[8]  << std::setw(9) << y32[8]  << std::endl;
    std::cout << std::setw(7) << "9:  " << std::setw(9) << u[9]  << std::setw(9) << y0[9]  << std::setw(9) << y1[9]  << std::setw(9) << y2[9]  << std::setw(9) << y4[9]  << std::setw(9) << y8[9]  << std::setw(9) << y16[9]  << std::setw(9) << y32[9]  << std::endl;
    std::cout << std::setw(7) << "10: " << std::setw(9) << u[10] << std::setw(9) << y0[10] << std::setw(9) << y1[10] << std::setw(9) << y2[10] << std::setw(9) << y4[10] << std::setw(9) << y8[10] << std::setw(9) << y16[10] << std::setw(9) << y32[10] << std::endl;
    std::cout << std::setw(7) << "11: " << std::setw(9) << u[11] << std::setw(9) << y0[11] << std::setw(9) << y1[11] << std::setw(9) << y2[11] << std::setw(9) << y4[11] << std::setw(9) << y8[11] << std::setw(9) << y16[11] << std::setw(9) << y32[11] << std::endl;
    std::cout << std::setw(7) << "12: " << std::setw(9) << u[12] << std::setw(9) << y0[12] << std::setw(9) << y1[12] << std::setw(9) << y2[12] << std::setw(9) << y4[12] << std::setw(9) << y8[12] << std::setw(9) << y16[12] << std::setw(9) << y32[12] << std::endl;
    std::cout << std::setw(7) << "13: " << std::setw(9) << u[13] << std::setw(9) << y0[13] << std::setw(9) << y1[13] << std::setw(9) << y2[13] << std::setw(9) << y4[13] << std::setw(9) << y8[13] << std::setw(9) << y16[13] << std::setw(9) << y32[13] << std::endl;
    std::cout << std::setw(7) << "14: " << std::setw(9) << u[14] << std::setw(9) << y0[14] << std::setw(9) << y1[14] << std::setw(9) << y2[14] << std::setw(9) << y4[14] << std::setw(9) << y8[14] << std::setw(9) << y16[14] << std::setw(9) << y32[14] << std::endl;
    std::cout << std::setw(7) << "15: " << std::setw(9) << u[15] << std::setw(9) << y0[15] << std::setw(9) << y1[15] << std::setw(9) << y2[15] << std::setw(9) << y4[15] << std::setw(9) << y8[15] << std::setw(9) << y16[15] << std::setw(9) << y32[15] << std::endl;
    std::cout << std::setw(7) << "16: " << std::setw(9) << u[16] << std::setw(9) << y0[16] << std::setw(9) << y1[16] << std::setw(9) << y2[16] << std::setw(9) << y4[16] << std::setw(9) << y8[16] << std::setw(9) << y16[16] << std::setw(9) << y32[16] << std::endl;
    std::cout << std::setw(7) << "17: " << std::setw(9) << u[17] << std::setw(9) << y0[17] << std::setw(9) << y1[17] << std::setw(9) << y2[17] << std::setw(9) << y4[17] << std::setw(9) << y8[17] << std::setw(9) << y16[17] << std::setw(9) << y32[17] << std::endl;
    std::cout << std::setw(7) << "18: " << std::setw(9) << u[18] << std::setw(9) << y0[18] << std::setw(9) << y1[18] << std::setw(9) << y2[18] << std::setw(9) << y4[18] << std::setw(9) << y8[18] << std::setw(9) << y16[18] << std::setw(9) << y32[18] << std::endl;
    std::cout << std::setw(7) << "19: " << std::setw(9) << u[19] << std::setw(9) << y0[19] << std::setw(9) << y1[19] << std::setw(9) << y2[19] << std::setw(9) << y4[19] << std::setw(9) << y8[19] << std::setw(9) << y16[19] << std::setw(9) << y32[19] << std::endl;
    std::cout << std::setw(7) << "20: " << std::setw(9) << u[20] << std::setw(9) << y0[20] << std::setw(9) << y1[20] << std::setw(9) << y2[20] << std::setw(9) << y4[20] << std::setw(9) << y8[20] << std::setw(9) << y16[20] << std::setw(9) << y32[20] << std::endl;


} // closes test1

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////










