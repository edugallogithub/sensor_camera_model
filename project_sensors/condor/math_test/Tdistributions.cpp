#include "Tdistributions.h"
#include "math/logic/timer.h"

#include <cstdlib>
#include <iostream>
#include <random>
#include <functional>

math::test::Tdistributions::Tdistributions(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tdistributions::run() {
	::jail::unit_test::run();

   test_normal_distribution();
   test_generator();
   test_normal_one_distribution_two_generators();
   test_random();
   // test_compare_generators();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdistributions::test_normal_distribution() {
    // This test proves that a given generator created without seed always produces the same sequence
    std::default_random_engine Ogenerator1;
    std::normal_distribution<double> Odistribution1(0.0,1.0);
    unsigned long nel = 10;
    std::vector<double> Ares(nel), Bres(nel), Cres(nel);
    for (int i=0; i < nel; ++i) {
        Ares[i] = 2.0 * Odistribution1(Ogenerator1);
        std::cout << Ares[i] << std::endl;
    }

    std::cout << "===================" << std::endl;

    std::default_random_engine Ogenerator2;
    std::normal_distribution<double> Odistribution2(0.0,2.0);
    for (int i=0; i < nel; ++i) {
        Bres[i] = Odistribution2(Ogenerator2);
        std::cout << Bres[i] << std::endl;
    }

    std::cout << "===================" << std::endl;

    std::default_random_engine Ogenerator3;
    std::normal_distribution<double> Odistribution3(0.0,2.0);
    auto Onormal = std::bind (Odistribution3, Ogenerator3);
    for (int i=0; i < nel; ++i) {
        Cres[i] = Onormal();
        std::cout << Cres[i] << std::endl;
    }

    check("Normal distribution 00a  ", Ares[0], Bres[0], 1e-12);
    check("Normal distribution 01a  ", Ares[1], Bres[1], 1e-12);
    check("Normal distribution 02a  ", Ares[2], Bres[2], 1e-12);
    check("Normal distribution 03a  ", Ares[3], Bres[3], 1e-12);
    check("Normal distribution 04a  ", Ares[4], Bres[4], 1e-12);
    check("Normal distribution 05a  ", Ares[5], Bres[5], 1e-12);
    check("Normal distribution 06a  ", Ares[6], Bres[6], 1e-12);
    check("Normal distribution 07a  ", Ares[7], Bres[7], 1e-12);
    check("Normal distribution 08a  ", Ares[8], Bres[8], 1e-12);
    check("Normal distribution 09a  ", Ares[9], Bres[9], 1e-12);
    check("Normal distribution 00b  ", Ares[0], Cres[0], 1e-12);
    check("Normal distribution 01b  ", Ares[1], Cres[1], 1e-12);
    check("Normal distribution 02b  ", Ares[2], Cres[2], 1e-12);
    check("Normal distribution 03b  ", Ares[3], Cres[3], 1e-12);
    check("Normal distribution 04b  ", Ares[4], Cres[4], 1e-12);
    check("Normal distribution 05b  ", Ares[5], Cres[5], 1e-12);
    check("Normal distribution 06b  ", Ares[6], Cres[6], 1e-12);
    check("Normal distribution 07b  ", Ares[7], Cres[7], 1e-12);
    check("Normal distribution 08b  ", Ares[8], Cres[8], 1e-12);
    check("Normal distribution 09b  ", Ares[9], Cres[9], 1e-12);
} // closes test_normal_distribution

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdistributions::test_generator() {
    // This test proves that a given generator created with the same seed always produces the same sequence
    // The problem now is how to generate the seeds
    unsigned seed1 = 432432;
    unsigned seed2 = 324324;

    std::default_random_engine Ogenerator1(seed1);
    std::normal_distribution<double> Odistribution1(0.0,1.0);
    unsigned long nel = 10;
    std::vector<double> Ares(nel), Bres(nel), Cres(nel);
    for (int i=0; i < nel; ++i) {
        Ares[i] = 2.0 * Odistribution1(Ogenerator1);
        std::cout << Ares[i] << std::endl;
    }

    std::cout << "===================" << std::endl;

    std::default_random_engine Ogenerator2(seed2);
    std::normal_distribution<double> Odistribution2(0.0,2.0);
    for (int i=0; i < nel; ++i) {
        Bres[i] = Odistribution2(Ogenerator2);
        std::cout << Bres[i] << std::endl;
    }

    std::cout << "===================" << std::endl;

    std::default_random_engine Ogenerator3(seed2);
    std::normal_distribution<double> Odistribution3(0.0,2.0);
    auto Onormal = std::bind (Odistribution3, Ogenerator3);
    for (int i=0; i < nel; ++i) {
        Cres[i] = Onormal();
        std::cout << Cres[i] << std::endl;
    }
    check("Normal distribution 00a  ", Bres[0], Cres[0], 1e-12);
    check("Normal distribution 01a  ", Bres[1], Cres[1], 1e-12);
    check("Normal distribution 02a  ", Bres[2], Cres[2], 1e-12);
    check("Normal distribution 03a  ", Bres[3], Cres[3], 1e-12);
    check("Normal distribution 04a  ", Bres[4], Cres[4], 1e-12);
    check("Normal distribution 05a  ", Bres[5], Cres[5], 1e-12);
    check("Normal distribution 06a  ", Bres[6], Cres[6], 1e-12);
    check("Normal distribution 07a  ", Bres[7], Cres[7], 1e-12);
    check("Normal distribution 08a  ", Bres[8], Cres[8], 1e-12);
    check("Normal distribution 09a  ", Bres[9], Cres[9], 1e-12);
} // closes test_generator

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdistributions::test_normal_one_distribution_two_generators() {
    // This text proves that all distributions of the same type are identical
    // and the results only depend on the generator, not the distribution.
    // However, distributions have memory, and their results do not only depend
    // on the generator but also on their previous history. Note below that
    // c4 is different than a4 or b4, although their inputs (created by Ogen11,
    // Ogen12, and Ogen13) are identical. The reason is that a4 and b4 are the
    // 4th execution of Odist, while c4 is only the 3rd execution of Odist1.

    unsigned int seed1 = 1; // Ogen11, Ogen12 and Ogen13 should be identical as created with same seed
    std::ranlux24_base Ogen11(seed1);
    std::ranlux24_base Ogen12(seed1);
    std::ranlux24_base Ogen13(seed1);
    unsigned int seed2 = 888;  // Ogen21, Ogen22 and Ogen23 should be identical as created with same seed
    std::ranlux24_base Ogen21(seed2);
    std::ranlux24_base Ogen22(seed2);
    std::ranlux24_base Ogen23(seed2);
    std::normal_distribution<double> Odist;; // Odist, Odist1, and Odist2 should be identical
    std::normal_distribution<double> Odist1;
    std::normal_distribution<double> Odist2;

    double a1 = Odist(Ogen11);
    double a2 = Odist(Ogen11);
    double a3 = Odist(Ogen21);
    double a4 = Odist(Ogen11);

    double b1 = Odist(Ogen12);
    double b2 = Odist(Ogen12);
    double b3 = Odist(Ogen22);
    double b4 = Odist(Ogen12);

    double c1 = Odist1(Ogen13);
    double c2 = Odist1(Ogen13);
    double c3 = Odist2(Ogen23);
    double c4 = Odist1(Ogen13);

    check("Generator equality    ", a1, b1, 1e-12);
    check("Generator equality    ", a2, b2, 1e-12);
    check("Generator equality    ", a3, b3, 1e-12);
    check("Generator equality    ", a4, b4, 1e-12);
    check("Distribution equality ", a1, c1, 1e-12);
    check("Distribution equality ", a2, c2, 1e-12);
    check("Distribution equality ", a3, c3, 1e-12);
    //check("Distribution equality ", a4, c4, 1e-12); // do not uncomment, they should be different (whole point of test)
} // closes test_normal_one_distribution_two_generators

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdistributions::test_random() {
    // This shows how to generate repeated seeds
    // It remains to be seen which generators are fastest
    unsigned seed = 1;
    std::default_random_engine Ogen(seed);
    std::uniform_int_distribution<> Odistr;
    for (int i=0; i < 10; ++i) {
        int number = Odistr(Ogen);
        std::cout << number << std::endl;
    }

    std::cout << "===================" << std::endl;
} // closes test_random

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template<typename T>
class Random {
private:
    T generator;
public:
    Random() : generator(std::chrono::high_resolution_clock::now().time_since_epoch().count()) {}

    int generate_integer(int begin, int end) {
        return std::uniform_int_distribution<int>(begin, end - 1)(generator);
    }
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tdistributions::test_compare_generators() {
    // This is a speed comparison test among the different generators. ONLY VALID IN RELEASE.
    // Best are ranlux24_base and ranlux48_base
    // Worst (by far) are ranlux24, ranlux48, knuth_b

    constexpr int n = 300000000;
    Random<std::default_random_engine>  Odefault_random_engine; // only one not portable (diff results in diff machines)
    Random<std::minstd_rand>            Ominstd_rand;
    Random<std::minstd_rand0>           Ominstd_rand0;
    Random<std::mt19937>                Omt19937;
    Random<std::mt19937_64>             Omt19937_64;
    Random<std::ranlux24_base>          Oranlux24_base;
    Random<std::ranlux48_base>          Oranlux48_base;
    //Random<std::ranlux24>               Oranlux24;
    //Random<std::ranlux48>               Oranlux48;
    //Random<std::knuth_b>                Oknuth_b;

    math::timer t;
    for (int j = 0; j < 3; ++j) {

        t.start();
        for (int i = 0; i < n; ++i) {(std::rand() % 10);}
        t.stop();
        std::cout << "rand                      " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Odefault_random_engine.generate_integer(0, 10);}
        t.stop();
        std::cout << "default_random_engine     " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Ominstd_rand.generate_integer(0, 10);}
        t.stop();
        std::cout << "minstd_rand               " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Ominstd_rand0.generate_integer(0, 10);}
        t.stop();
        std::cout << "minstd_rand0              " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Omt19937.generate_integer(0, 10);}
        t.stop();
        std::cout << "mt19937                   " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Omt19937_64.generate_integer(0, 10);}
        t.stop();
        std::cout << "mt19937_64                " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Oranlux24_base.generate_integer(0, 10);}
        t.stop();
        std::cout << "ranlux24_base             " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Oranlux48_base.generate_integer(0, 10);}
        t.stop();
        std::cout << "ranlux48_base             " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        /*t.start();
        for (int i = 0; i < n; ++i) {Oranlux24.generate_integer(0, 10);}
        t.stop();
        std::cout << "ranlux24                  " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Oranlux48.generate_integer(0, 10);}
        t.stop();
        std::cout << "ranlux48                  " << t.measure_nanosec() /n << " ns per rand" << std::endl;

        t.start();
        for (int i = 0; i < n; ++i) {Oknuth_b.generate_integer(0, 10);}
        t.stop();
        std::cout << "knuth_b                   " << t.measure_nanosec() /n << " ns per rand" << std::endl;*/
    }
} // closes test_compare_generators

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////






