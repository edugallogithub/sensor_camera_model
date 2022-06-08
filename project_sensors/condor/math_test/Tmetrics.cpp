#include "Tmetrics.h"
#include "math/logic/share.h"

#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <list>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

math::test::Tmetrics::Tmetrics(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tmetrics::run() {
	::jail::unit_test::run();

    test_classic();
    test_robust();
    test_comparison_std_mad();
    //test_plot();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tmetrics::test_classic() {
    std::vector<int> Vi;
    std::vector<double> Vd;
    std::list<int> Li;
    std::list<double> Ld;
    std::vector<Eigen::Array3i> Ve(10);
    std::vector<Eigen::Array3d> Vf(10);

    Vi.push_back(1);    Vd.push_back(1.);   Li.push_back(1);    Ld.push_back(1.);   Ve[0] << 1, 10, 100;    Vf[0] << 1., 10., 100.;
    Vi.push_back(2);    Vd.push_back(2.);   Li.push_back(2);    Ld.push_back(2.);   Ve[1] << 2, 20, 200;    Vf[1] << 2., 20., 200.;
    Vi.push_back(3);    Vd.push_back(3.);   Li.push_back(3);    Ld.push_back(3.);   Ve[2] << 3, 30, 300;    Vf[2] << 3., 30., 300.;
    Vi.push_back(4);    Vd.push_back(4.);   Li.push_back(4);    Ld.push_back(4.);   Ve[3] << 4, 40, 400;    Vf[3] << 4., 40., 400.;
    Vi.push_back(5);    Vd.push_back(5.);   Li.push_back(5);    Ld.push_back(5.);   Ve[4] << 5, 50, 500;    Vf[4] << 5., 50., 500.;
    Vi.push_back(6);    Vd.push_back(6.);   Li.push_back(6);    Ld.push_back(6.);   Ve[5] << 6, 60, 600;    Vf[5] << 6., 60., 600.;
    Vi.push_back(7);    Vd.push_back(7.);   Li.push_back(7);    Ld.push_back(7.);   Ve[6] << 7, 70, 700;    Vf[6] << 7., 70., 700.;
    Vi.push_back(8);    Vd.push_back(8.);   Li.push_back(8);    Ld.push_back(8.);   Ve[7] << 8, 80, 800;    Vf[7] << 8., 80., 800.;
    Vi.push_back(9);    Vd.push_back(9.);   Li.push_back(9);    Ld.push_back(9.);   Ve[8] << 9, 90, 900;    Vf[8] << 9., 90., 900.;
    Vi.push_back(10);   Vd.push_back(10.);  Li.push_back(10);   Ld.push_back(10.);  Ve[9] << 10, 100, 1000; Vf[9] << 10., 100., 1000.;

    int    mean_i = (1  + 2  + 3  + 4  + 5  + 6  + 7  + 8  + 9  + 10)  / 10;
    double mean_d = (1. + 2. + 3. + 4. + 5. + 6. + 7. + 8. + 9. + 10.) / 10;
    int    mean_e2 = (10  + 20  + 30  + 40  + 50  + 60  + 70  + 80  + 90  + 100)  / 10;
    int    mean_e3 = (100  + 200  + 300  + 400  + 500  + 600  + 700  + 800  + 900  + 1000)  / 10;
    Eigen::Array3i mean_e(mean_i, mean_e2, mean_e3);
    Eigen::Array3d mean_f(mean_d, 10. * mean_d, 100. * mean_d);

    auto   std_i = (int)(sqrt(((1 - mean_i)*(1 - mean_i) + (2 - mean_i)*(2 - mean_i) + (3 - mean_i)*(3 - mean_i) + (4 - mean_i)*(4 - mean_i) + (5  - mean_i)*(5  - mean_i) +
                               (6 - mean_i)*(6 - mean_i) + (7 - mean_i)*(7 - mean_i) + (8 - mean_i)*(8 - mean_i) + (9 - mean_i)*(9 - mean_i) + (10 - mean_i)*(10 - mean_i)) / 10));
    double std_d = sqrt(((1. - mean_d)*(1. - mean_d) + (2. - mean_d)*(2. - mean_d) + (3. - mean_d)*(3. - mean_d) + (4. - mean_d)*(4. - mean_d) + (5.  - mean_d)*(5.  - mean_d) +
                         (6. - mean_d)*(6. - mean_d) + (7. - mean_d)*(7. - mean_d) + (8. - mean_d)*(8. - mean_d) + (9. - mean_d)*(9. - mean_d) + (10. - mean_d)*(10. - mean_d)) / 10.);
    auto   std_e2 = (int)(sqrt(((10 - mean_e2)*(10 - mean_e2) + (20 - mean_e2)*(20 - mean_e2) + (30 - mean_e2)*(30 - mean_e2) + (40 - mean_e2)*(40 - mean_e2) + (50  - mean_e2)*(50  - mean_e2) +
                                (60 - mean_e2)*(60 - mean_e2) + (70 - mean_e2)*(70 - mean_e2) + (80 - mean_e2)*(80 - mean_e2) + (90 - mean_e2)*(90 - mean_e2) + (100 - mean_e2)*(100 - mean_e2)) / 10));
    auto   std_e3 = (int)(sqrt(((100 - mean_e3)*(100 - mean_e3) + (200 - mean_e3)*(200 - mean_e3) + (300 - mean_e3)*(300 - mean_e3) + (400 - mean_e3)*(400 - mean_e3) + (500  - mean_e3)*(500  - mean_e3) +
                                (600 - mean_e3)*(600 - mean_e3) + (700 - mean_e3)*(700 - mean_e3) + (800 - mean_e3)*(800 - mean_e3) + (900 - mean_e3)*(900 - mean_e3) + (1000 - mean_e3)*(1000 - mean_e3)) / 10));
    Eigen::Array3i std_e(std_i, std_e2, std_e3);
    Eigen::Array3d std_f(std_d, 10. * std_d, 100. * std_d);
    
    auto   rms_i = (int)(sqrt((1*1 + 2*2 + 3*3 + 4*4 + 5*5 + 6*6 + 7*7 + 8*8 + 9*9 + 10*10) / 10));
    double rms_d = sqrt((1.*1. + 2.*2. + 3.*3. + 4.*4. + 5.*5. + 6.*6. + 7.*7. + 8.*8. + 9.*9. + 10.*10.) / 10.);
    auto   rms_e2 = (int)(sqrt((10*10 + 20*20 + 30*30 + 40*40 + 50*50 + 60*60 + 70*70 + 80*80 + 90*90 + 100*100) / 10));
    auto   rms_e3 = (int)(sqrt((100*100 + 200*200 + 300*300 + 400*400 + 500*500 + 600*600 + 700*700 + 800*800 + 900*900 + 1000*1000) / 10));
    Eigen::Array3i rms_e(rms_i, rms_e2, rms_e3);
    Eigen::Array3d rms_f(rms_d, 10. * rms_d, 100. * rms_d);

    int    smax_i = 10;
    double smax_d = 10.;
    int    max_i = 10;
    double max_d = 10.;
    int    min_i = 1;
    double min_d = 1.;
    
    int            mean1_Vi = math::mean(Vi.begin(), Vi.end());
    double         mean1_Vd = math::mean(Vd.begin(), Vd.end());
    Eigen::Array3i mean1_Ve = math::mean(Ve.begin(), Ve.end());
    Eigen::Array3d mean1_Vf = math::mean(Vf.begin(), Vf.end());

    int            mean2_Vi = math::mean(Vi);
    double         mean2_Vd = math::mean(Vd);
    int            mean2_Li = math::mean(Li);
    double         mean2_Ld = math::mean(Ld);
    Eigen::Array3i mean2_Ve = math::mean(Ve);
    Eigen::Array3d mean2_Vf = math::mean(Vf);

    int            mean3_Vi = math::mean_vector(Vi);
    double         mean3_Vd = math::mean_vector(Vd);
    Eigen::Array3i mean3_Ve = math::mean_vector(Ve);
    Eigen::Array3d mean3_Vf = math::mean_vector(Vf);

    check("mean1_Vi    ", mean_i,    mean1_Vi,    1e-15);
    check("mean1_Vd    ", mean_d,    mean1_Vd,    1e-15);
    check("mean1_Ve(0) ", mean_e(0), mean1_Ve(0), 1e-15);
    check("mean1_Ve(1) ", mean_e(1), mean1_Ve(1), 1e-14);
    check("mean1_Ve(2) ", mean_e(2), mean1_Ve(2), 1e-13);
    check("mean1_Vf(0) ", mean_f(0), mean1_Vf(0), 1e-15);
    check("mean1_Vf(1) ", mean_f(1), mean1_Vf(1), 1e-14);
    check("mean1_Vf(2) ", mean_f(2), mean1_Vf(2), 1e-13);

    check("mean2_Vi    ", mean_i,    mean2_Vi,    1e-15);
    check("mean2_Vd    ", mean_d,    mean2_Vd,    1e-15);
    check("mean2_Li    ", mean_i,    mean2_Li,    1e-15);
    check("mean2_Ld    ", mean_d,    mean2_Ld,    1e-15);
    check("mean2_Ve(0) ", mean_e(0), mean2_Ve(0), 1e-15);
    check("mean2_Ve(1) ", mean_e(1), mean2_Ve(1), 1e-14);
    check("mean2_Ve(2) ", mean_e(2), mean2_Ve(2), 1e-13);
    check("mean2_Vf(0) ", mean_f(0), mean2_Vf(0), 1e-15);
    check("mean2_Vf(1) ", mean_f(1), mean2_Vf(1), 1e-14);
    check("mean2_Vf(2) ", mean_f(2), mean2_Vf(2), 1e-13);

    check("mean3_Vi    ", mean_i,    mean3_Vi,    1e-15);
    check("mean3_Vd    ", mean_d,    mean3_Vd,    1e-15);
    check("mean3_Ve(0) ", mean_e(0), mean3_Ve(0), 1e-15);
    check("mean3_Ve(1) ", mean_e(1), mean3_Ve(1), 1e-14);
    check("mean3_Ve(2) ", mean_e(2), mean3_Ve(2), 1e-13);
    check("mean3_Vf(0) ", mean_f(0), mean3_Vf(0), 1e-15);
    check("mean3_Vf(1) ", mean_f(1), mean3_Vf(1), 1e-14);
    check("mean3_Vf(2) ", mean_f(2), mean3_Vf(2), 1e-13);

    int            std1_Vi = math::std(Vi.begin(), Vi.end(), mean1_Vi);
    double         std1_Vd = math::std(Vd.begin(), Vd.end(), mean1_Vd);
    Eigen::Array3i std1_Ve = math::std(Ve.begin(), Ve.end(), mean1_Ve);
    Eigen::Array3d std1_Vf = math::std(Vf.begin(), Vf.end(), mean1_Vf);

    int            std2_Vi = math::std(Vi, mean2_Vi);
    double         std2_Vd = math::std(Vd, mean2_Vd);
    int            std2_Li = math::std(Li, mean2_Li);
    double         std2_Ld = math::std(Ld, mean2_Ld);
    Eigen::Array3i std2_Ve = math::std(Ve, mean2_Ve);
    Eigen::Array3d std2_Vf = math::std(Vf, mean2_Vf);

    int            std3_Vi = math::std_vector(Vi, mean3_Vi);
    double         std3_Vd = math::std_vector(Vd, mean3_Vd);
    Eigen::Array3i std3_Ve = math::std_vector(Ve, mean3_Ve);
    Eigen::Array3d std3_Vf = math::std_vector(Vf, mean3_Vf);

    check("std1_Vi     ", std_i,    std1_Vi,    1e-15);
    check("std1_Vd     ", std_d,    std1_Vd,    1e-15);
    check("std1_Ve(0)  ", std_e(0), std1_Ve(0), 1e-15);
    check("std1_Ve(1)  ", std_e(1), std1_Ve(1), 1e-14);
    check("std1_Ve(2)  ", std_e(2), std1_Ve(2), 1e-13);
    check("std1_Vf(0)  ", std_f(0), std1_Vf(0), 1e-15);
    check("std1_Vf(1)  ", std_f(1), std1_Vf(1), 1e-14);
    check("std1_Vf(2)  ", std_f(2), std1_Vf(2), 1e-13);

    check("std2_Vi     ", std_i,    std2_Vi,    1e-15);
    check("std2_Vd     ", std_d,    std2_Vd,    1e-15);
    check("std2_Li     ", std_i,    std2_Li,    1e-15);
    check("std2_Ld     ", std_d,    std2_Ld,    1e-15);
    check("std2_Ve(0)  ", std_e(0), std2_Ve(0), 1e-15);
    check("std2_Ve(1)  ", std_e(1), std2_Ve(1), 1e-14);
    check("std2_Ve(2)  ", std_e(2), std2_Ve(2), 1e-13);
    check("std2_Vf(0)  ", std_f(0), std2_Vf(0), 1e-15);
    check("std2_Vf(1)  ", std_f(1), std2_Vf(1), 1e-14);
    check("std2_Vf(2)  ", std_f(2), std2_Vf(2), 1e-13);

    check("std3_Vi     ", std_i,    std3_Vi,    1e-15);
    check("std3_Vd     ", std_d,    std3_Vd,    1e-15);
    check("std3_Ve(0)  ", std_e(0), std3_Ve(0), 1e-15);
    check("std3_Ve(1)  ", std_e(1), std3_Ve(1), 1e-14);
    check("std3_Ve(2)  ", std_e(2), std3_Ve(2), 1e-13);
    check("std3_Vf(0)  ", std_f(0), std3_Vf(0), 1e-15);
    check("std3_Vf(1)  ", std_f(1), std3_Vf(1), 1e-14);
    check("std3_Vf(2)  ", std_f(2), std3_Vf(2), 1e-13);

    int            rms1_Vi = math::rms(Vi.begin(), Vi.end());
    double         rms1_Vd = math::rms(Vd.begin(), Vd.end());
    Eigen::Array3i rms1_Ve = math::rms(Ve.begin(), Ve.end());
    Eigen::Array3d rms1_Vf = math::rms(Vf.begin(), Vf.end());

    int            rms2_Vi = math::rms(Vi);
    double         rms2_Vd = math::rms(Vd);
    int            rms2_Li = math::rms(Li);
    double         rms2_Ld = math::rms(Ld);
    Eigen::Array3i rms2_Ve = math::rms(Ve);
    Eigen::Array3d rms2_Vf = math::rms(Vf);

    int            rms3_Vi = math::rms_vector(Vi);
    double         rms3_Vd = math::rms_vector(Vd);
    Eigen::Array3i rms3_Ve = math::rms_vector(Ve);
    Eigen::Array3d rms3_Vf = math::rms_vector(Vf);

    check("rms1_Vi     ", rms_i,    rms1_Vi,    1e-15);
    check("rms1_Vd     ", rms_d,    rms1_Vd,    1e-15);
    check("rms1_Ve(0)  ", rms_e(0), rms1_Ve(0), 1e-15);
    check("rms1_Ve(1)  ", rms_e(1), rms1_Ve(1), 1e-14);
    check("rms1_Ve(2)  ", rms_e(2), rms1_Ve(2), 1e-13);
    check("rms1_Vf(0)  ", rms_f(0), rms1_Vf(0), 1e-15);
    check("rms1_Vf(1)  ", rms_f(1), rms1_Vf(1), 1e-14);
    check("rms1_Vf(2)  ", rms_f(2), rms1_Vf(2), 1e-13);

    check("rms2_Vi     ", rms_i,    rms2_Vi,    1e-15);
    check("rms2_Vd     ", rms_d,    rms2_Vd,    1e-15);
    check("rms2_Li     ", rms_i,    rms2_Li,    1e-15);
    check("rms2_Ld     ", rms_d,    rms2_Ld,    1e-15);
    check("rms2_Ve(0)  ", rms_e(0), rms2_Ve(0), 1e-15);
    check("rms2_Ve(1)  ", rms_e(1), rms2_Ve(1), 1e-14);
    check("rms2_Ve(2)  ", rms_e(2), rms2_Ve(2), 1e-13);
    check("rms2_Vf(0)  ", rms_f(0), rms2_Vf(0), 1e-15);
    check("rms2_Vf(1)  ", rms_f(1), rms2_Vf(1), 1e-14);
    check("rms2_Vf(2)  ", rms_f(2), rms2_Vf(2), 1e-13);

    check("rms3_Vi     ", rms_i,    rms3_Vi,    1e-15);
    check("rms3_Vd     ", rms_d,    rms3_Vd,    1e-15);
    check("rms3_Ve(0)  ", rms_e(0), rms3_Ve(0), 1e-15);
    check("rms3_Ve(1)  ", rms_e(1), rms3_Ve(1), 1e-14);
    check("rms3_Ve(2)  ", rms_e(2), rms3_Ve(2), 1e-13);
    check("rms3_Vf(0)  ", rms_f(0), rms3_Vf(0), 1e-15);
    check("rms3_Vf(1)  ", rms_f(1), rms3_Vf(1), 1e-14);
    check("rms3_Vf(2)  ", rms_f(2), rms3_Vf(2), 1e-13);

    int    smax1_Vi = math::smax(Vi.begin(), Vi.end());
    double smax1_Vd = math::smax(Vd.begin(), Vd.end());
    int    smax1_Li = math::smax(Li.begin(), Li.end());
    double smax1_Ld = math::smax(Ld.begin(), Ld.end());
    int    smax2_Vi = math::smax(Vi);
    double smax2_Vd = math::smax(Vd);
    int    smax2_Li = math::smax(Li);
    double smax2_Ld = math::smax(Ld);
    int    smax3_Vi = math::smax_vector(Vi);
    double smax3_Vd = math::smax_vector(Vd);

    check("smax1_Vi  ", smax_i, smax1_Vi, 1e-15);
    check("smax1_Vd  ", smax_d, smax1_Vd, 1e-15);
    check("smax1_Li  ", smax_i, smax1_Li, 1e-15);
    check("smax1_Ld  ", smax_d, smax1_Ld, 1e-15);
    check("smax2_Vi  ", smax_i, smax2_Vi, 1e-15);
    check("smax2_Vd  ", smax_d, smax2_Vd, 1e-15);
    check("smax2_Li  ", smax_i, smax2_Li, 1e-15);
    check("smax2_Ld  ", smax_d, smax2_Ld, 1e-15);
    check("smax3_Vi  ", smax_i, smax3_Vi, 1e-15);
    check("smax3_Vd  ", smax_d, smax3_Vd, 1e-15);

    unsigned long pos1, pos2, pos3, pos4;
    int    max3_Vi = math::max_vector_pos(Vi, pos1);
    double max3_Vd = math::max_vector_pos(Vd, pos2);
    int    min3_Vi = math::min_vector_pos(Vi, pos3);
    double min3_Vd = math::min_vector_pos(Vd, pos4);

    check("max3_Vi  ", max_i, max3_Vi, 1e-15);
    check("max3_Vd  ", max_d, max3_Vd, 1e-15);
    check("min3_Vi  ", min_i, min3_Vi, 1e-15);
    check("min3_Vd  ", min_d, min3_Vd, 1e-15);
    check("max3_Vi  ", pos1, 9);
    check("max3_Vd  ", pos2, 9);
    check("min3_Vi  ", pos3, 0);
    check("min3_Vd  ", pos4, 0);
} // closes test_classic

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tmetrics::test_robust() {
    std::vector<int> Vi;
    std::vector<double> Vd;
    std::vector<int> Vi2;
    std::vector<double> Vd2;

    Vi.push_back(7);    Vd.push_back(7.);   Vi2.push_back(7);    Vd2.push_back(7.);
    Vi.push_back(3);    Vd.push_back(3.);   Vi2.push_back(3);    Vd2.push_back(3.);
    Vi.push_back(1);    Vd.push_back(1.);   Vi2.push_back(1);    Vd2.push_back(1.);
    Vi.push_back(4);    Vd.push_back(4.);   Vi2.push_back(4);    Vd2.push_back(4.);
    Vi.push_back(10);   Vd.push_back(10.);  Vi2.push_back(10);   Vd2.push_back(10.);
    Vi.push_back(9);    Vd.push_back(9.);   Vi2.push_back(9);    Vd2.push_back(9.);
    Vi.push_back(2);    Vd.push_back(2.);   Vi2.push_back(2);    Vd2.push_back(2.);
    Vi.push_back(5);    Vd.push_back(5.);   Vi2.push_back(5);    Vd2.push_back(5.);
    Vi.push_back(6);    Vd.push_back(6.);   Vi2.push_back(6);    Vd2.push_back(6.);
    Vi.push_back(11);   Vd.push_back(11.);  Vi2.push_back(11);   Vd2.push_back(11.);
    Vi.push_back(8);    Vd.push_back(8.);   Vi2.push_back(8);    Vd2.push_back(8.);

    int median_i1    = math::median(Vi);
    double median_d1 = math::median(Vd);
    int median_i2    = math::median_vector(Vi2);
    double median_d2 = math::median_vector(Vd2);

    check("median_i  ", median_i1, 6,      1e-15);
    check("median_i  ", median_i1, Vi[5],  1e-15);
    check("median_d  ", median_d1, 6.,     1e-15);
    check("median_d  ", median_d1, Vd[5],  1e-15);
    check("median_i  ", median_i2, 6,      1e-15);
    check("median_i  ", median_i2, Vi2[5], 1e-15);
    check("median_d  ", median_d2, 6.,     1e-15);
    check("median_d  ", median_d2, Vd2[5], 1e-15);

    int mad_i1    = math::mad(Vi, median_i1);
    double mad_d1 = math::mad(Vd, median_d1);
    int mad_i2    = math::mad_vector(Vi2, median_i2);
    double mad_d2 = math::mad_vector(Vd2, median_d2);

    check("mad_i  ", mad_i1, 3,      1e-15);
    check("mad_d  ", mad_d1, 3.,     1e-15);
    check("mad_i  ", mad_i2, 3,      1e-15);
    check("mad_d  ", mad_d2, 3.,     1e-15);

} // closes test_robust

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tmetrics::test_comparison_std_mad() {
    std::default_random_engine Ogenerator1;
    double mean0 = 10.0;
    double std0  = 3.0;
    std::normal_distribution<double> Odistribution1(mean0, std0);
    unsigned long nel = 1000;
    std::vector<double> Vd(nel);
    for (int i = 0; i < nel; ++i) {
        Vd[i] = Odistribution1(Ogenerator1);
    }
    double mean1 = math::mean(Vd);
    double std1  = math::std(Vd, mean1);

    Vd.push_back(250.0);
    double mean2 = math::mean(Vd);
    double std2  = math::std(Vd, mean1);

    double median = math::median(Vd);
    double mad    = 1.4826 * math::mad(Vd, median); // multiply so it can be compared with std

    std::cout << std::setw(14) << "mean"
              << std::setw(14) << "std"
              << std::setw(14) << "mean error"
              << std::setw(14) << "std error"
              << std::endl;
    std::cout << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean0 - mean0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std0 - std0
              << std::endl;
    std::cout << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean1
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std1
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean1 - mean0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std1 - std0
              << std::endl;
    std::cout << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean2
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std2
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mean2 - mean0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << std2 - std0
              << std::endl;
    std::cout << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << median
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mad
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << median - mean0
              << std::scientific << std::setw(14) << std::setprecision(5) << std::showpos << mad - std0
              << std::endl;

} // closes test_comparison_std_mad

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tmetrics::test_plot() {
    // generates a metric_type.txt file to be employed in thesis chapter 12 (scenarios) for the different types of metrics
    const unsigned int seed1 = 1, seed2 = 88, seed3 = 200;
    std::default_random_engine Ogenerator1(seed1);
    std::default_random_engine Ogenerator2(seed2);
    std::default_random_engine Ogenerator3(seed3);
    std::normal_distribution<double> Onormal(0.0,1.0);

    int nel = 1001;
    std::vector<double> Vmu1(nel,0.), Vstd1(nel,0.), Vstd2(nel,0.), Vstd3(nel,0.), Vcurve1(nel,0.), Vcurve2(nel,0.);
    std::vector<double> Vnobias_bounded(nel,0.), Vnobias_bounded1(nel,0.), Vnobias_bounded2(nel,0.);
    std::vector<double> Vbias_bounded(nel,0.), Vbias_bounded1(nel,0.), Vbias_bounded2(nel,0.);
    std::vector<double> Vnobias_drift(nel,0.), Vnobias_drift1(nel,0.), Vnobias_drift2(nel,0.);
    std::vector<double> Vbias_drift(nel,0.), Vbias_drift1(nel,0.), Vbias_drift2(nel,0.);
    std::vector<double> Vbounded(nel,0.), Vbounded1(nel,0.), Vbounded2(nel,0.);
    std::vector<double> Vdrift(nel,0.), Vdrift1(nel,0.), Vdrift2(nel,0.);

    for (int x=0; x != nel; ++x) {
        Vmu1[x] = 0.05 * Onormal(Ogenerator1);
    }

    double Ax1 = 100.;
    double Ax2 = 200.;
    double Ay0 = 1.0;
    double Ay1 = 2.5;
    double Ay2 = 3.0;
    double Ac = Ay0;
    double Aa = (Ay2 - 2. * Ay1 + Ay0) / (2 * Ax1 * Ax1);
    double Ab = (Ay1 - Ay0 - Aa * Ax1 * Ax1) / Ax1;
    for (int x=0; x != nel; ++x) {
        if (x <= Ax2) {Vstd1[x] = Aa * x * x + Ab * x + Ac;}
        else {Vstd1[x] = Vstd1[Ax2];}
    }

    double Bx1 = 250.;
    double Bx2 = 500.;
    double By0 = 0.0;
    double By1 = 0.4;
    double By2 = 1.0;
    double Bc = By0;
    double Ba = (By2 - 2. * By1 + By0) / (2 * Bx1 * Bx1);
    double Bb = (By1 - By0 - Ba * Bx1 * Bx1) / Bx1;
    for (int x=0; x != nel; ++x) {
        if (x <= Bx2) {Vcurve1[x] = Ba * x * x + Bb * x + Bc;}
        else {Vcurve1[x] = Vcurve1[Bx2];}
    }

    double Cx1 = 200.;
    double Cx2 = 400.;
    double Cy0 = 0.8;
    double Cy1 = 2.0;
    double Cy2 = 3.0;
    double Cc = Cy0;
    double Ca = (Cy2 - 2. * Cy1 + Cy0) / (2 * Cx1 * Cx1);
    double Cb = (Cy1 - Cy0 - Ca * Cx1 * Cx1) / Cx1;
    for (int x=0; x != nel; ++x) {
        if (x <= Cx2) {Vstd2[x] = Ca * x * x + Cb * x + Cc;}
        else {Vstd2[x] = Vstd2[Cx2];}
    }

    double Dx1 =  500.;
    double Dx2 = 1000.;
    double Dy0 = 0.;
    double Dy1 = -0.7;
    double Dy2 = -1.6;
    double Dc = Dy0;
    double Da = (Dy2 - 2. * Dy1 + Dy0) / (2 * Dx1 * Dx1);
    double Db = (Dy1 - Dy0 - Da * Dx1 * Dx1) / Dx1;
    for (int x=0; x != nel; ++x) {
        if (x <= Dx2) {Vstd3[x] = Da * x * x + Db * x + Dc;}
        else {Vstd3[x] = Vstd3[Dx2];}
    }

    double Ex1 =  500.;
    double Ex2 = 1000.;
    double Ey0 = 0.0;
    double Ey1 = 0.8;
    double Ey2 = 1.5;
    double Ec = Ey0;
    double Ea = (Ey2 - 2. * Ey1 + Ey0) / (2 * Ex1 * Ex1);
    double Eb = (Ey1 - Ey0 - Ea * Ex1 * Ex1) / Ex1;
    for (int x=0; x != nel; ++x) {
        if (x <= Ex2) { Vcurve2[x] = Ea * x * x + Eb * x + Ec; }
        else { Vcurve2[x] = Vcurve2[Ex2]; }
    }

    for (int x=0; x != nel; ++x) {
        Vnobias_bounded[x]  = Vmu1[x];
        Vnobias_bounded1[x] = Vnobias_bounded[x] + Vstd1[x];
        Vnobias_bounded2[x] = Vnobias_bounded[x] - Vstd1[x];
        Vbias_bounded[x]    = Vmu1[x] + Vcurve1[x] + 0.5;
        Vbias_bounded1[x]   = Vbias_bounded[x] + Vstd2[x];
        Vbias_bounded2[x]   = Vbias_bounded[x] - Vstd2[x];
        Vnobias_drift[x]    = Vmu1[x];
        Vnobias_drift1[x]   = Vnobias_drift[x] + 0.7 * Vstd1[x] + 2.0e-3 * x;
        Vnobias_drift2[x]   = Vnobias_drift[x] - 0.7 * Vstd1[x] - 2.0e-3 * x;
        Vbias_drift[x]      = Vmu1[x] + Vcurve1[x] + 0.5;
        Vbias_drift1[x]     = Vbias_drift[x] + 0.7 * Vstd1[x] + 2.0e-3 * x;
        Vbias_drift2[x]     = Vbias_drift[x] - 0.7 * Vstd1[x] - 2.0e-3 * x;
        Vbounded[x]         = Vmu1[x] + Vstd1[x];
        Vbounded1[x]        = Vbounded[x] + 0.7 * Vstd1[x];
        Vbounded2[x]        = Vbounded[x] - 0.7 * Vstd1[x];
        Vdrift[x]           = Vmu1[x] + 0.85 * Vstd1[x] - 0.3 + 1.3e-3 * x;
        Vdrift1[x]          = Vdrift[x] + 0.7 * Vstd1[x];
        Vdrift2[x]          = Vdrift[x] - 0.7 * Vstd1[x];
    }

    std::string st_folder = math::share::condor_output_hard_disk;
    boost::filesystem::path path_folder(st_folder);
    std::string st_env("env");
    std::string st_file("metric_type.txt");
    std::string st_file_complete = (path_folder / st_env / st_file).string();

    std::ofstream Oout;
    Oout.open(st_file_complete);
    for (int x=0; x != nel; ++x) {
        Oout << std::fixed << std::setw(4) << std::setprecision(0) << std::noshowpos << x
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_bounded[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_bounded1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_bounded2[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_bounded[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_bounded1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_bounded2[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_drift[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_drift1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vnobias_drift2[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_drift[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_drift1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbias_drift2[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbounded[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbounded1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vbounded2[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vdrift[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vdrift1[x]
             << std::fixed << std::setw(8) << std::setprecision(3) << std::showpos   << Vdrift2[x]
             << std::endl;
    }
    Oout.close();
} // closes test_plot

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////






































