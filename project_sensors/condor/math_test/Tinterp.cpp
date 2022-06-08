#include "Tinterp.h"
#include "math/vec/vec3.h"
#include "math/vec/algorithm.h"

math::test::Tinterp::Tinterp(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tinterp::run() {
	::jail::unit_test::run();

	test_search();
	test_spline();
    test_hermite_1d(math::logic::hermite_second);
    test_hermite_1d(math::logic::hermite_first);
	test_hermite_2d(math::logic::hermite_first);
	test_hermite_2d(math::logic::hermite_second);
	test_hermite_3d(math::logic::hermite_first);
	test_hermite_3d(math::logic::hermite_second);

	finished();
}
/* execute tests and write results on console */


/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void math::test::Tinterp::test_spline() {

	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// This script is MatLab is where the results for the tests are taken from
	// This test does exactly the same and hence should obtain the same results
	// >> points = [2., 2.5, 3.5, 5.5, 7.5, 10.5]';
	// >> values = [10.0, 20.0, 15.0, 5.0, 25.0, 50.0]';
	// >> Xpoints = [1. , 2., 2.5, 3.5, 4., 5.5, 7.5, 9.2, 10.5, 15.]';
	// >> cs = csape(points, values);
	// >> csdiff = fnder(cs);
	// >> fnval(cs, Xpoints) // Evaluates function
	// >> fnval(csdiff, Xpoints) // Evaluates function differentials
	// >> (fnval(cs, Xpoints+0.001) - fnval(cs, Xpoints-0.001)) / 0.002 // another way of computing differentials
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

	math::vec1 points(6);
	points.set(0, 2.0);
	points.set(1, 2.5);
	points.set(2, 3.5);
	points.set(3, 5.5);
	points.set(4, 7.5);
	points.set(5, 10.5);
	math::vec1 values(6);
	values.set(0, 10.);
	values.set(1, 20.);
	values.set(2, 15.);
	values.set(3, 5.);
	values.set(4, 25.);
	values.set(5, 50.);

	math::vec1 Xpoints(10);
	Xpoints.set(0, 1.0);
	Xpoints.set(1, 2.0);
	Xpoints.set(2, 2.5);
	Xpoints.set(3, 3.5);
	Xpoints.set(4, 4.0);
	Xpoints.set(5, 5.5);
	Xpoints.set(6, 7.5);
	Xpoints.set(7, 9.2);
	Xpoints.set(8, 10.5);
	Xpoints.set(9, 15.0);

	std::vector<double> Xvalues(10);
	math::f_table1V Otab(points, values, math::logic::spline);
	for (int i = 0; i != Xpoints.size1(); ++i) {
        Xvalues[i] = Otab.value(Xpoints[i]);
	}

	check("Spline-001 Evaluation        ", 	Xvalues[0],   -49.2524,     1e-4);
	check("Spline-002 Evaluation        ", 	Xvalues[1],   values[0],  1e-4);
	check("Spline-003 Evaluation        ", 	Xvalues[2],   values[1],  1e-4);
	check("Spline-004 Evaluation        ", 	Xvalues[3],   values[2],  1e-4);
	check("Spline-005 Evaluation        ", 	Xvalues[4],   10.0717,      1e-4);
	check("Spline-006 Evaluation        ", 	Xvalues[5],   values[3],  1e-4);
	check("Spline-007 Evaluation        ", 	Xvalues[6],   values[4],  1e-4);
	check("Spline-008 Evaluation        ", 	Xvalues[7],   44.8350,      1e-4);
	check("Spline-009 Evaluation        ", 	Xvalues[8],   values[5],  1e-4);
	check("Spline-010 Evaluation        ", 	Xvalues[9],   -102.2331,    1e-4);

	math::vec1 Xdiffs(10);
	double dum(1.);
	for (int i = 0; i != Xpoints.size1(); ++i) {
        Xdiffs[i] = Otab.d_dt(Xpoints[i], dum);
	}

	check("Spline-021 Differential      ", 	Xdiffs[0],   88.9589,      1e-1);
	check("Spline-022 Differential      ", 	Xdiffs[1],   31.9048,      1e-1);
	check("Spline-023 Differential      ", 	Xdiffs[2],   8.6849,       1e-1);
	check("Spline-024 Differential      ", 	Xdiffs[3],   -10.9191,     1e-1);
	check("Spline-025 Differential      ", 	Xdiffs[4],   -8.6551,      1e-1);
	check("Spline-026 Differential      ", 	Xdiffs[5],   3.1448,       1e-1);
	check("Spline-027 Differential      ", 	Xdiffs[6],   13.3399,      1e-1);
	check("Spline-028 Differential      ", 	Xdiffs[7],   8.4722,       1e-1);
	check("Spline-029 Differential      ", 	Xdiffs[8],   -1.4167,      1e-1);
	check("Spline-030 Differential      ", 	Xdiffs[9],   -76.9152,     1e-1);

	double stop_here = 88;

} // closes test_spline

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void math::test::Tinterp::test_search() {
	math::vec1 Ovec1(10);
	Ovec1.set(0, 10);
	Ovec1.set(1, 20);
	Ovec1.set(2, 30);
	Ovec1.set(3, 40);
	Ovec1.set(4, 50);
	Ovec1.set(5, 60);
	Ovec1.set(6, 70);
	Ovec1.set(7, 80);
	Ovec1.set(8, 90);
	Ovec1.set(9, 100);
	double diff = 10;
	
	std::vector<double> Oinput(22);
	Oinput[0] = -25;
	Oinput[1] = 0;
	Oinput[2] = 5;
	Oinput[3] = 10;
	Oinput[4] = 11;
	Oinput[5] = 19;
	Oinput[6] = 20;
	Oinput[7] = 24;
	Oinput[8] = 28;
	Oinput[9] = 30;
	Oinput[10] = 32;
	Oinput[11] = 40;
	Oinput[12] = 70;
	Oinput[13] = 74;
	Oinput[14] = 80;
	Oinput[15] = 87;
	Oinput[16] = 90;
	Oinput[17] = 99;
	Oinput[18] = 100;
	Oinput[19] = 101;
	Oinput[20] = 105;
	Oinput[21] = 125;

	math::interp_lagrange_first Ilin;
	math::interp_lagrange_second Iqua;
	math::interp_lagrange_third Icub;

	std::vector<int> X(22, -999);
	std::vector<int> Y(22, -999);
	std::vector<int> M2(22, -999);
	std::vector<int> M3(22, -999);
	std::vector<int> M4(22, -999);

	X[0]  = math::search_binary(Oinput[0], Ovec1);
	Y[0]  = math::search_equispaced(Oinput[0], Ovec1, diff);
	M2[0]  = Ilin.find_index(X[0], 10);  
	M3[0]  = Iqua.find_index(X[0], 10);  
	M4[0]  = Icub.find_index(X[0], 10);  

	X[1]  = math::search_binary(Oinput[1], Ovec1);
	Y[1]  = math::search_equispaced(Oinput[1], Ovec1, diff);
	M2[1]  = Ilin.find_index(X[1], 10);  
	M3[1]  = Iqua.find_index(X[1], 10);  
	M4[1]  = Icub.find_index(X[1], 10);  

	X[2]  = math::search_binary(Oinput[2], Ovec1);
	Y[2]  = math::search_equispaced(Oinput[2], Ovec1, diff);
	M2[2]  = Ilin.find_index(X[2], 10);  
	M3[2]  = Iqua.find_index(X[2], 10);  
	M4[2]  = Icub.find_index(X[2], 10);  

	X[3]  = math::search_binary(Oinput[3], Ovec1);
	Y[3]  = math::search_equispaced(Oinput[3], Ovec1, diff);
	M2[3]  = Ilin.find_index(X[3], 10);  
	M3[3]  = Iqua.find_index(X[3], 10);  
	M4[3]  = Icub.find_index(X[3], 10);  

	X[4]  = math::search_binary(Oinput[4], Ovec1);
	Y[4]  = math::search_equispaced(Oinput[4], Ovec1, diff);
	M2[4]  = Ilin.find_index(X[4], 10);  
	M3[4]  = Iqua.find_index(X[4], 10);  
	M4[4]  = Icub.find_index(X[4], 10);  

	X[5]  = math::search_binary(Oinput[5], Ovec1);
	Y[5]  = math::search_equispaced(Oinput[5], Ovec1, diff);
	M2[5]  = Ilin.find_index(X[5], 10);   
	M3[5]  = Iqua.find_index(X[5], 10);  
	M4[5]  = Icub.find_index(X[5], 10);  

	X[6]  = math::search_binary(Oinput[6], Ovec1);
	Y[6]  = math::search_equispaced(Oinput[6], Ovec1, diff);
	M2[6]  = Ilin.find_index(X[6], 10);  
	M3[6]  = Iqua.find_index(X[6], 10);  
	M4[6]  = Icub.find_index(X[6], 10);  

	X[7]  = math::search_binary(Oinput[7], Ovec1);
	Y[7]  = math::search_equispaced(Oinput[7], Ovec1, diff);
	M2[7]  = Ilin.find_index(X[7], 10);  
	M3[7]  = Iqua.find_index(X[7], 10);  
	M4[7]  = Icub.find_index(X[7], 10);  

	X[8]  = math::search_binary(Oinput[8], Ovec1);
	Y[8]  = math::search_equispaced(Oinput[8], Ovec1, diff);
	M2[8]  = Ilin.find_index(X[8], 10);  
	M3[8]  = Iqua.find_index(X[8], 10);  
	M4[8]  = Icub.find_index(X[8], 10);  

	X[9]  = math::search_binary(Oinput[9], Ovec1);
	Y[9]  = math::search_equispaced(Oinput[9], Ovec1, diff);
	M2[9]  = Ilin.find_index(X[9], 10);  
	M3[9]  = Iqua.find_index(X[9], 10);  
	M4[9]  = Icub.find_index(X[9], 10);  

	X[10] = math::search_binary(Oinput[10], Ovec1);
	Y[10] = math::search_equispaced(Oinput[10], Ovec1, diff);
	M2[10]  = Ilin.find_index(X[10], 10);  
	M3[10]  = Iqua.find_index(X[10], 10);  
	M4[10]  = Icub.find_index(X[10], 10);  

	X[11] = math::search_binary(Oinput[11], Ovec1);
	Y[11] = math::search_equispaced(Oinput[11], Ovec1, diff);
	M2[11]  = Ilin.find_index(X[11], 10);
	M3[11]  = Iqua.find_index(X[11], 10);  
	M4[11]  = Icub.find_index(X[11], 10);  

	X[12] = math::search_binary(Oinput[12], Ovec1);
	Y[12] = math::search_equispaced(Oinput[12], Ovec1, diff);
	M2[12]  = Ilin.find_index(X[12], 10);  
	M3[12]  = Iqua.find_index(X[12], 10);  
	M4[12]  = Icub.find_index(X[12], 10);  

	X[13] = math::search_binary(Oinput[13], Ovec1);
	Y[13] = math::search_equispaced(Oinput[13], Ovec1, diff);
	M2[13]  = Ilin.find_index(X[13], 10);  
	M3[13]  = Iqua.find_index(X[13], 10);  
	M4[13]  = Icub.find_index(X[13], 10);  

	X[14] = math::search_binary(Oinput[14], Ovec1);
	Y[14] = math::search_equispaced(Oinput[14], Ovec1, diff);
	M2[14]  = Ilin.find_index(X[14], 10);  
	M3[14]  = Iqua.find_index(X[14], 10);  
	M4[14]  = Icub.find_index(X[14], 10);  

	X[15] = math::search_binary(Oinput[15], Ovec1);
	Y[15] = math::search_equispaced(Oinput[15], Ovec1, diff);
	M2[15]  = Ilin.find_index(X[15], 10);  
	M3[15]  = Iqua.find_index(X[15], 10);  
	M4[15]  = Icub.find_index(X[15], 10);  

	X[16] = math::search_binary(Oinput[16], Ovec1);
	Y[16] = math::search_equispaced(Oinput[16], Ovec1, diff);
	M2[16]  = Ilin.find_index(X[16], 10);  
	M3[16]  = Iqua.find_index(X[16], 10);  
	M4[16]  = Icub.find_index(X[16], 10);  

	X[17] = math::search_binary(Oinput[17], Ovec1);
	Y[17] = math::search_equispaced(Oinput[17], Ovec1, diff);
	M2[17]  = Ilin.find_index(X[17], 10);  
	M3[17]  = Iqua.find_index(X[17], 10);  
	M4[17]  = Icub.find_index(X[17], 10);  

	X[18] = math::search_binary(Oinput[18], Ovec1);
	Y[18] = math::search_equispaced(Oinput[18], Ovec1, diff);
	M2[18]  = Ilin.find_index(X[18], 10);  
	M3[18]  = Iqua.find_index(X[18], 10);  
	M4[18]  = Icub.find_index(X[18], 10);  

	X[19] = math::search_binary(Oinput[19], Ovec1);
	Y[19] = math::search_equispaced(Oinput[19], Ovec1, diff);
	M2[19]  = Ilin.find_index(X[19], 10);  
	M3[19]  = Iqua.find_index(X[19], 10);  
	M4[19]  = Icub.find_index(X[19], 10);  

	X[20] = math::search_binary(Oinput[20], Ovec1);
	Y[20] = math::search_equispaced(Oinput[20], Ovec1, diff);
	M2[20]  = Ilin.find_index(X[20], 10);  
	M3[20]  = Iqua.find_index(X[20], 10);  
	M4[20]  = Icub.find_index(X[20], 10);  

	X[21] = math::search_binary(Oinput[21], Ovec1);
	Y[21] = math::search_equispaced(Oinput[21], Ovec1, diff);
	M2[21]  = Ilin.find_index(X[21], 10);  
	M3[21]  = Iqua.find_index(X[21], 10);  
	M4[21]  = Icub.find_index(X[21], 10);  

	check("Search-001 Equal searches         ",	X[0],	Y[0]);
	check("Search-002 Equal searches         ",	X[1],	Y[1]);
	check("Search-003 Equal searches         ",	X[2],	Y[2]);
	check("Search-004 Equal searches         ",	X[3],	Y[3]);
	check("Search-005 Equal searches         ",	X[4],	Y[4]);
	check("Search-006 Equal searches         ",	X[5],	Y[5]);
	check("Search-007 Equal searches         ",	X[6],	Y[6]);
	check("Search-008 Equal searches         ",	X[7],	Y[7]);
	check("Search-009 Equal searches         ",	X[8],	Y[8]);
	check("Search-010 Equal searches         ",	X[9],	Y[9]);
	check("Search-011 Equal searches         ",	X[10],	Y[10]);
	check("Search-012 Equal searches         ",	X[11],	Y[11]);
	check("Search-013 Equal searches         ",	X[12],	Y[12]);
	check("Search-014 Equal searches         ",	X[13],	Y[13]);
	check("Search-015 Equal searches         ",	X[14],	Y[14]);
	check("Search-016 Equal searches         ",	X[15],	Y[15]);
	check("Search-017 Equal searches         ",	X[16],	Y[16]);
	check("Search-018 Equal searches         ",	X[17],	Y[17]);
	check("Search-019 Equal searches         ",	X[18],	Y[18]);
	check("Search-020 Equal searches         ",	X[19],	Y[19]);
	check("Search-021 Equal searches         ",	X[20],	Y[20]);
	check("Search-022 Equal searches         ",	X[21],	Y[21]);

	check("Search-031 Find for linear       ",	M2[0],	X[0] + 1);
	check("Search-032 Find for quadratic    ",	M3[0],	X[0] + 1);
	check("Search-033 Find for cubic        ",	M4[0],	X[0] + 1);
	check("Search-034 Find for linear       ",	M2[1],	X[1] + 1);
	check("Search-035 Find for quadratic    ",	M3[1],	X[1] + 1);
	check("Search-036 Find for cubic        ",	M4[1],	X[1] + 1);
	check("Search-037 Find for linear       ",	M2[2],	X[2] + 1);
	check("Search-038 Find for quadratic    ",	M3[2],	X[2] + 1);
	check("Search-039 Find for cubic        ",	M4[2],	X[2] + 1);
	check("Search-040 Find for linear       ",	M2[3],	X[3]);
	check("Search-041 Find for quadratic    ",	M3[3],	X[3]);
	check("Search-042 Find for cubic        ",	M4[3],	X[3]);
	check("Search-043 Find for linear       ",	M2[4],	X[4]);
	check("Search-044 Find for quadratic    ",	M3[4],	X[4]);
	check("Search-045 Find for cubic        ",	M4[4],	X[4]);
	check("Search-046 Find for linear       ",	M2[5],	X[5]);
	check("Search-047 Find for quadratic    ",	M3[5],	X[5]);
	check("Search-048 Find for cubic        ",	M4[5],	X[5]);
	check("Search-049 Find for linear       ",	M2[6],	X[6]);
	check("Search-050 Find for quadratic    ",	M3[6],	X[6] - 1);
	check("Search-051 Find for cubic        ",	M4[6],	X[6] - 1);
	check("Search-052 Find for linear       ",	M2[7],	X[7]);
	check("Search-053 Find for quadratic    ",	M3[7],	X[7] - 1);
	check("Search-054 Find for cubic        ",	M4[7],	X[7] - 1);
	check("Search-055 Find for linear       ",	M2[8],	X[8]);
	check("Search-056 Find for quadratic    ",	M3[8],	X[8] - 1);
	check("Search-057 Find for cubic        ",	M4[8],	X[8] - 1);
	check("Search-058 Find for linear       ",	M2[9],	X[9]);
	check("Search-059 Find for quadratic    ",	M3[9],	X[9] - 1);
	check("Search-060 Find for cubic        ",	M4[9],	X[9] - 1);
	check("Search-061 Find for linear       ",	M2[10],	X[10]);
	check("Search-062 Find for quadratic    ",	M3[10],	X[10] - 1);
	check("Search-063 Find for cubic        ",	M4[10],	X[10] - 1);
	check("Search-064 Find for linear       ",	M2[11],	X[11]);
	check("Search-065 Find for quadratic    ",	M3[11],	X[11] - 1);
	check("Search-066 Find for cubic        ",	M4[11],	X[11] - 1);
	check("Search-067 Find for linear       ",	M2[12],	X[12]);
	check("Search-068 Find for quadratic    ",	M3[12],	X[12] - 1);
	check("Search-069 Find for cubic        ",	M4[12],	X[12] - 1);
	check("Search-070 Find for linear       ",	M2[13],	X[13]);
	check("Search-071 Find for quadratic    ",	M3[13],	X[13] - 1);
	check("Search-072 Find for cubic        ",	M4[13],	X[13] - 1);
	check("Search-073 Find for linear       ",	M2[14],	X[14]);
	check("Search-074 Find for quadratic    ",	M3[14],	X[14] - 1);
	check("Search-075 Find for cubic        ",	M4[14],	X[14] - 1);
	check("Search-076 Find for linear       ",	M2[15],	X[15]);
	check("Search-077 Find for quadratic    ",	M3[15],	X[15] - 1);
	check("Search-078 Find for cubic        ",	M4[15],	X[15] - 1);
	check("Search-079 Find for linear       ",	M2[16],	X[16]);
	check("Search-080 Find for quadratic    ",	M3[16],	X[16] - 1);
	check("Search-081 Find for cubic        ",	M4[16],	X[16] - 2);
	check("Search-082 Find for linear       ",	M2[17],	X[17]);
	check("Search-083 Find for quadratic    ",	M3[17],	X[17] - 1);
	check("Search-084 Find for cubic        ",	M4[17],	X[17] - 2);
	check("Search-085 Find for linear       ",	M2[18],	X[18] - 1);
	check("Search-086 Find for quadratic    ",	M3[18],	X[18] - 2);
	check("Search-087 Find for cubic        ",	M4[18],	X[18] - 3);
	check("Search-088 Find for linear       ",	M2[19],	X[19] - 1);
	check("Search-089 Find for quadratic    ",	M3[19],	X[19] - 2);
	check("Search-090 Find for cubic        ",	M4[19],	X[19] - 3);
	check("Search-090 Find for linear       ",	M2[20],	X[20] - 1);
	check("Search-091 Find for quadratic    ",	M3[20],	X[20] - 2);
	check("Search-092 Find for cubic        ",	M4[20],	X[20] - 3);
	check("Search-093 Find for linear       ",	M2[21],	X[21] - 1);
	check("Search-094 Find for quadratic    ",	M3[21],	X[21] - 2);
	check("Search-095 Find for cubic        ",	M4[21],	X[21] - 3);

	double stop_here = 0.;
} // closes test_search

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tinterp::test_hermite_1d(math::logic::INTERP_MODE interp_mode) {
	math::vec1 points(6);
	points.set(0, 2.0);
	points.set(1, 2.5);
	points.set(2, 3.5);
	points.set(3, 5.5);
	points.set(4, 7.5);
	points.set(5, 10.5);
	math::vec1 values(6);
	values.set(0, 10.);
	values.set(1, 20.);
	values.set(2, 15.);
	values.set(3, 5.);
	values.set(4, 25.);
	values.set(5, 50.);
	math::f_table1V Otab(points, values, interp_mode);

	math::vec1 Xpoints(10);
	Xpoints.set(0, 1.0);
	Xpoints.set(1, 2.0);
	Xpoints.set(2, 2.5);
	Xpoints.set(3, 3.5);
	Xpoints.set(4, 4.0);
	Xpoints.set(5, 5.5);
	Xpoints.set(6, 7.5);
	Xpoints.set(7, 9.2);
	Xpoints.set(8, 10.5);
	Xpoints.set(9, 15.0);

	std::vector<double> Xvalues(10);
	for (int i = 0; i != Xpoints.size1(); ++i) {
        Xvalues[i] = Otab.value(Xpoints[i]);
	}

	check("Hermite_1v-002 Evaluation        ", 	Xvalues[1],   values[0],  1e-10);
	check("Hermite_1v-003 Evaluation        ", 	Xvalues[2],   values[1],  1e-10);
	check("Hermite_1v-004 Evaluation        ", 	Xvalues[3],   values[2],  1e-10);
	check("Hermite_1v-006 Evaluation        ", 	Xvalues[5],   values[3],  1e-10);
	check("Hermite_1v-007 Evaluation        ", 	Xvalues[6],   values[4],  1e-10);
	check("Hermite_1v-009 Evaluation        ", 	Xvalues[8],   values[5],  1e-10);

	math::vec1 Xdiffs(10);
    double dum(1.);
	for (int i = 0; i != Xpoints.size1(); ++i) {
        Xdiffs[i] = Otab.d_dt(Xpoints[i], dum);
	}

	check("Hermite_1v-022 Differential      ", 	Xdiffs[1],   Otab.get_slopes()[0], 	1e-1);
	check("Hermite_1v-023 Differential      ", 	Xdiffs[2],   Otab.get_slopes()[1],     1e-1);
	check("Hermite_1v-024 Differential      ", 	Xdiffs[3],   Otab.get_slopes()[2],     1e-1);
	check("Hermite_1v-026 Differential      ", 	Xdiffs[5],   Otab.get_slopes()[3],     1e-1);
	check("Hermite_1v-027 Differential      ", 	Xdiffs[6],   Otab.get_slopes()[4],     1e-1);
	check("Hermite_1v-029 Differential      ", 	Xdiffs[8],   Otab.get_slopes()[5],     1e-1);
} // closes test_hermite_1d

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void math::test::Tinterp::test_hermite_2d(math::logic::INTERP_MODE interp_mode) {
	math::vec1 pointsE(6);
	pointsE.set(0, 2.0);
	pointsE.set(1, 2.5);
	pointsE.set(2, 3.5);
	pointsE.set(3, 5.5);
	pointsE.set(4, 7.5);
	pointsE.set(5, 10.5);
	math::vec1 pointsF(4);
	pointsF.set(0, 50.0);
	pointsF.set(1, 54.5);
	pointsF.set(2, 65.5);
	pointsF.set(3, 72.0);

	math::vec2 Values(4, 6);
	Values.set(0, 0, 10.);	Values.set(1, 0, 15.);	Values.set(2, 0, 12.);	Values.set(3, 0, 22.);
	Values.set(0, 1, 2.);	Values.set(1, 1, 12.);	Values.set(2, 1, 17.);	Values.set(3, 1, 12.);
	Values.set(0, 2, 84.);	Values.set(1, 2, 5.);	Values.set(2, 2, 34.);	Values.set(3, 2, 72.);
	Values.set(0, 3, 32.);	Values.set(1, 3, 55.);	Values.set(2, 3, 66.);	Values.set(3, 3, 41.);
	Values.set(0, 4, 51.);	Values.set(1, 4, 67.);	Values.set(2, 4, 5.);	Values.set(3, 4, 34.);
	Values.set(0, 5, 11.);	Values.set(1, 5, 72.);	Values.set(2, 5, 88.);	Values.set(3, 5, 64.);
	math::f_table2V Otab(pointsE, pointsF, Values, interp_mode);

	// Make simmetric copy to check results
	math::vec2 AValues(6, 4);
	AValues.set(0, 0, 10.);	AValues.set(1, 0, 2.);	AValues.set(2, 0, 84.);	
	AValues.set(3, 0, 32.);	AValues.set(4, 0, 51.);	AValues.set(5, 0, 11.);	
	AValues.set(0, 1, 15.);	AValues.set(1, 1, 12.);	AValues.set(2, 1, 5.);	
	AValues.set(3, 1, 55.);	AValues.set(4, 1, 67.);	AValues.set(5, 1, 72.);	
	AValues.set(0, 2, 12.);	AValues.set(1, 2, 17.);	AValues.set(2, 2, 34.);	
	AValues.set(3, 2, 66.);	AValues.set(4, 2, 5.);	AValues.set(5, 2, 88.);	
	AValues.set(0, 3, 22.);	AValues.set(1, 3, 12.);	AValues.set(2, 3, 72.);	
	AValues.set(3, 3, 41.);	AValues.set(4, 3, 34.);	AValues.set(5, 3, 64.);	
	math::f_table2V OAtab(pointsF, pointsE, AValues, interp_mode);

	math::vec1 XpointsE(10);
	XpointsE.set(0, 1.0);
	XpointsE.set(1, 2.0);
	XpointsE.set(2, 2.5);
	XpointsE.set(3, 3.5);
	XpointsE.set(4, 4.0);
	XpointsE.set(5, 5.5);
	XpointsE.set(6, 7.5);
	XpointsE.set(7, 9.2);
	XpointsE.set(8, 10.5);
	XpointsE.set(9, 15.0);

	math::vec1 XpointsF(7);
	XpointsF.set(0, 44.0);
	XpointsF.set(1, 50.0);
	XpointsF.set(2, 54.5);
	XpointsF.set(3, 59.2);
	XpointsF.set(4, 65.5);
	XpointsF.set(5, 72.0);
	XpointsF.set(6, 82.0);

	math::vec2 XValues(7, 10);

	for (int i = 0; i != XpointsF.size1(); ++i) {
		for (int j = 0; j != XpointsE.size1(); ++j) {
			//std::cout << "i= " << i << " - " << XpointsF[i] << std::endl;
			//std::cout << "j= " << j << " - " << XpointsE[j] << std::endl;
            XValues.get(i, j) = Otab.value(XpointsE[j], XpointsF[i]);
			//std::cout << XValues.get(i,j) << std::endl << std::endl;
		}
	}

	double a = Values.get(0, 0);
	double b = XValues.get(1, 1);

	// Check 2D table evaluation at the nodes
	check("Hermite_2v-001 Evaluation        ", 	XValues.get(1, 1),   Values.get(0, 0),     1e-4);
	check("Hermite_2v-002 Evaluation        ", 	XValues.get(2, 1),   Values.get(1, 0),     1e-4);
	check("Hermite_2v-003 Evaluation        ", 	XValues.get(4, 1),   Values.get(2, 0),     1e-4);
	check("Hermite_2v-004 Evaluation        ", 	XValues.get(5, 1),   Values.get(3, 0),     1e-4);
	check("Hermite_2v-011 Evaluation        ", 	XValues.get(1, 2),   Values.get(0, 1),     1e-4);
	check("Hermite_2v-012 Evaluation        ", 	XValues.get(2, 2),   Values.get(1, 1),     1e-4);
	check("Hermite_2v-013 Evaluation        ", 	XValues.get(4, 2),   Values.get(2, 1),     1e-4);
	check("Hermite_2v-014 Evaluation        ", 	XValues.get(5, 2),   Values.get(3, 1),     1e-4);
	check("Hermite_2v-021 Evaluation        ", 	XValues.get(1, 3),   Values.get(0, 2),     1e-4);
	check("Hermite_2v-022 Evaluation        ", 	XValues.get(2, 3),   Values.get(1, 2),     1e-4);
	check("Hermite_2v-023 Evaluation        ", 	XValues.get(4, 3),   Values.get(2, 2),     1e-4);
	check("Hermite_2v-024 Evaluation        ", 	XValues.get(5, 3),   Values.get(3, 2),     1e-4);
	check("Hermite_2v-031 Evaluation        ", 	XValues.get(1, 5),   Values.get(0, 3),     1e-4);
	check("Hermite_2v-032 Evaluation        ", 	XValues.get(2, 5),   Values.get(1, 3),     1e-4);
	check("Hermite_2v-033 Evaluation        ", 	XValues.get(4, 5),   Values.get(2, 3),     1e-4);
	check("Hermite_2v-034 Evaluation        ", 	XValues.get(5, 5),   Values.get(3, 3),     1e-4);
	check("Hermite_2v-041 Evaluation        ", 	XValues.get(1, 6),   Values.get(0, 4),     1e-4);
	check("Hermite_2v-042 Evaluation        ", 	XValues.get(2, 6),   Values.get(1, 4),     1e-4);
	check("Hermite_2v-043 Evaluation        ", 	XValues.get(4, 6),   Values.get(2, 4),     1e-4);
	check("Hermite_2v-044 Evaluation        ", 	XValues.get(5, 6),   Values.get(3, 4),     1e-4);
	check("Hermite_2v-051 Evaluation        ", 	XValues.get(1, 8),   Values.get(0, 5),     1e-4);
	check("Hermite_2v-052 Evaluation        ", 	XValues.get(2, 8),   Values.get(1, 5),     1e-4);
	check("Hermite_2v-053 Evaluation        ", 	XValues.get(4, 8),   Values.get(2, 5),     1e-4);
	check("Hermite_2v-054 Evaluation        ", 	XValues.get(5, 8),   Values.get(3, 5),     1e-4);

	math::vec2 XAValues(10, 7);
	for (int i = 0; i != XpointsF.size1(); ++i) {
		for (int j = 0; j != XpointsE.size1(); ++j) {
            XAValues.get(j, i) = OAtab.value(XpointsF[i], XpointsE[j]);
		}
	}

	// Check 2D simmetrical (or traspose) table evaluation at the nodes
	check("Hermite_2v-101 Evaluation trasp  ", 	XAValues.get(1, 1),   Values.get(0, 0),     1e-4);
	check("Hermite_2v-102 Evaluation trasp  ", 	XAValues.get(1, 2),   Values.get(1, 0),     1e-4);
	check("Hermite_2v-103 Evaluation trasp  ", 	XAValues.get(1, 4),   Values.get(2, 0),     1e-4);
	check("Hermite_2v-104 Evaluation trasp  ", 	XAValues.get(1, 5),   Values.get(3, 0),     1e-4);
	check("Hermite_2v-111 Evaluation trasp  ", 	XAValues.get(2, 1),   Values.get(0, 1),     1e-4);
	check("Hermite_2v-112 Evaluation trasp  ", 	XAValues.get(2, 2),   Values.get(1, 1),     1e-4);
	check("Hermite_2v-113 Evaluation trasp  ", 	XAValues.get(2, 4),   Values.get(2, 1),     1e-4);
	check("Hermite_2v-114 Evaluation trasp  ", 	XAValues.get(2, 5),   Values.get(3, 1),     1e-4);
	check("Hermite_2v-121 Evaluation trasp  ", 	XAValues.get(3, 1),   Values.get(0, 2),     1e-4);
	check("Hermite_2v-122 Evaluation trasp  ", 	XAValues.get(3, 2),   Values.get(1, 2),     1e-4);
	check("Hermite_2v-123 Evaluation trasp  ", 	XAValues.get(3, 4),   Values.get(2, 2),     1e-4);
	check("Hermite_2v-124 Evaluation trasp  ", 	XAValues.get(3, 5),   Values.get(3, 2),     1e-4);
	check("Hermite_2v-131 Evaluation trasp  ", 	XAValues.get(5, 1),   Values.get(0, 3),     1e-4);
	check("Hermite_2v-132 Evaluation trasp  ", 	XAValues.get(5, 2),   Values.get(1, 3),     1e-4);
	check("Hermite_2v-133 Evaluation trasp  ", 	XAValues.get(5, 4),   Values.get(2, 3),     1e-4);
	check("Hermite_2v-134 Evaluation trasp  ", 	XAValues.get(5, 5),   Values.get(3, 3),     1e-4);
	check("Hermite_2v-141 Evaluation trasp  ", 	XAValues.get(6, 1),   Values.get(0, 4),     1e-4);
	check("Hermite_2v-142 Evaluation trasp  ", 	XAValues.get(6, 2),   Values.get(1, 4),     1e-4);
	check("Hermite_2v-143 Evaluation trasp  ", 	XAValues.get(6, 4),   Values.get(2, 4),     1e-4);
	check("Hermite_2v-144 Evaluation trasp  ", 	XAValues.get(6, 5),   Values.get(3, 4),     1e-4);
	check("Hermite_2v-151 Evaluation trasp  ", 	XAValues.get(8, 1),   Values.get(0, 5),     1e-4);
	check("Hermite_2v-152 Evaluation trasp  ", 	XAValues.get(8, 2),   Values.get(1, 5),     1e-4);
	check("Hermite_2v-153 Evaluation trasp  ", 	XAValues.get(8, 4),   Values.get(2, 5),     1e-4);
	check("Hermite_2v-154 Evaluation trasp  ", 	XAValues.get(8, 5),   Values.get(3, 5),     1e-4);

	// Compare 2D table evaluation with its simmetrical both at the nodes and at other values
	check("Hermite_2v-201 Evaluation symm   ", 	XValues.get(0, 0),   XAValues.get(0, 0),     1e-4);
	check("Hermite_2v-202 Evaluation symm   ", 	XValues.get(1, 0),   XAValues.get(0, 1),     1e-4);
	check("Hermite_2v-203 Evaluation symm   ", 	XValues.get(2, 0),   XAValues.get(0, 2),     1e-4);
	check("Hermite_2v-204 Evaluation symm   ", 	XValues.get(3, 0),   XAValues.get(0, 3),     1e-4);
	check("Hermite_2v-205 Evaluation symm   ", 	XValues.get(4, 0),   XAValues.get(0, 4),     1e-4);
	check("Hermite_2v-206 Evaluation symm   ", 	XValues.get(5, 0),   XAValues.get(0, 5),     1e-4);
	check("Hermite_2v-207 Evaluation symm   ", 	XValues.get(6, 0),   XAValues.get(0, 6),     1e-4);
	
	check("Hermite_2v-211 Evaluation symm   ", 	XValues.get(0, 1),   XAValues.get(1, 0),     1e-4);
	check("Hermite_2v-212 Evaluation symm   ", 	XValues.get(1, 1),   XAValues.get(1, 1),     1e-4);
	check("Hermite_2v-213 Evaluation symm   ", 	XValues.get(2, 1),   XAValues.get(1, 2),     1e-4);
	check("Hermite_2v-214 Evaluation symm   ", 	XValues.get(3, 1),   XAValues.get(1, 3),     1e-4);
	check("Hermite_2v-215 Evaluation symm   ", 	XValues.get(4, 1),   XAValues.get(1, 4),     1e-4);
	check("Hermite_2v-216 Evaluation symm   ", 	XValues.get(5, 1),   XAValues.get(1, 5),     1e-4);
	check("Hermite_2v-217 Evaluation symm   ", 	XValues.get(6, 1),   XAValues.get(1, 6),     1e-4);

	check("Hermite_2v-221 Evaluation symm   ", 	XValues.get(0, 2),   XAValues.get(2, 0),     1e-4);
	check("Hermite_2v-222 Evaluation symm   ", 	XValues.get(1, 2),   XAValues.get(2, 1),     1e-4);
	check("Hermite_2v-223 Evaluation symm   ", 	XValues.get(2, 2),   XAValues.get(2, 2),     1e-4);
	check("Hermite_2v-224 Evaluation symm   ", 	XValues.get(3, 2),   XAValues.get(2, 3),     1e-4);
	check("Hermite_2v-225 Evaluation symm   ", 	XValues.get(4, 2),   XAValues.get(2, 4),     1e-4);
	check("Hermite_2v-226 Evaluation symm   ", 	XValues.get(5, 2),   XAValues.get(2, 5),     1e-4);
	check("Hermite_2v-227 Evaluation symm   ", 	XValues.get(6, 2),   XAValues.get(2, 6),     1e-4);

	check("Hermite_2v-231 Evaluation symm   ", 	XValues.get(0, 3),   XAValues.get(3, 0),     1e-4);
	check("Hermite_2v-232 Evaluation symm   ", 	XValues.get(1, 3),   XAValues.get(3, 1),     1e-4);
	check("Hermite_2v-233 Evaluation symm   ", 	XValues.get(2, 3),   XAValues.get(3, 2),     1e-4);
	check("Hermite_2v-234 Evaluation symm   ", 	XValues.get(3, 3),   XAValues.get(3, 3),     1e-4);
	check("Hermite_2v-231 Evaluation symm   ", 	XValues.get(4, 3),   XAValues.get(3, 4),     1e-4);
	check("Hermite_2v-232 Evaluation symm   ", 	XValues.get(5, 3),   XAValues.get(3, 5),     1e-4);
	check("Hermite_2v-233 Evaluation symm   ", 	XValues.get(6, 3),   XAValues.get(3, 6),     1e-4);

	check("Hermite_2v-241 Evaluation symm   ", 	XValues.get(0, 4),   XAValues.get(4, 0),     1e-4);
	check("Hermite_2v-242 Evaluation symm   ", 	XValues.get(1, 4),   XAValues.get(4, 1),     1e-4);
	check("Hermite_2v-243 Evaluation symm   ", 	XValues.get(2, 4),   XAValues.get(4, 2),     1e-4);
	check("Hermite_2v-244 Evaluation symm   ", 	XValues.get(3, 4),   XAValues.get(4, 3),     1e-4);
	check("Hermite_2v-245 Evaluation symm   ", 	XValues.get(4, 4),   XAValues.get(4, 4),     1e-4);
	check("Hermite_2v-246 Evaluation symm   ", 	XValues.get(5, 4),   XAValues.get(4, 5),     1e-4);
	check("Hermite_2v-247 Evaluation symm   ", 	XValues.get(6, 4),   XAValues.get(4, 6),     1e-4);

	check("Hermite_2v-251 Evaluation symm   ", 	XValues.get(0, 5),   XAValues.get(5, 0),     1e-4);
	check("Hermite_2v-252 Evaluation symm   ", 	XValues.get(1, 5),   XAValues.get(5, 1),     1e-4);
	check("Hermite_2v-253 Evaluation symm   ", 	XValues.get(2, 5),   XAValues.get(5, 2),     1e-4);
	check("Hermite_2v-254 Evaluation symm   ", 	XValues.get(3, 5),   XAValues.get(5, 3),     1e-4);
	check("Hermite_2v-255 Evaluation symm   ", 	XValues.get(4, 5),   XAValues.get(5, 4),     1e-4);
	check("Hermite_2v-256 Evaluation symm   ", 	XValues.get(5, 5),   XAValues.get(5, 5),     1e-4);
	check("Hermite_2v-257 Evaluation symm   ", 	XValues.get(6, 5),   XAValues.get(5, 6),     1e-4);

	check("Hermite_2v-261 Evaluation symm   ", 	XValues.get(0, 6),   XAValues.get(6, 0),     1e-4);
	check("Hermite_2v-262 Evaluation symm   ", 	XValues.get(1, 6),   XAValues.get(6, 1),     1e-4);
	check("Hermite_2v-263 Evaluation symm   ", 	XValues.get(2, 6),   XAValues.get(6, 2),     1e-4);
	check("Hermite_2v-264 Evaluation symm   ", 	XValues.get(3, 6),   XAValues.get(6, 3),     1e-4);
	check("Hermite_2v-265 Evaluation symm   ", 	XValues.get(4, 6),   XAValues.get(6, 4),     1e-4);
	check("Hermite_2v-266 Evaluation symm   ", 	XValues.get(5, 6),   XAValues.get(6, 5),     1e-4);
	check("Hermite_2v-267 Evaluation symm   ", 	XValues.get(6, 6),   XAValues.get(6, 6),     1e-4);

	check("Hermite_2v-271 Evaluation symm   ", 	XValues.get(0, 7),   XAValues.get(7, 0),     1e-4);
	check("Hermite_2v-272 Evaluation symm   ", 	XValues.get(1, 7),   XAValues.get(7, 1),     1e-4);
	check("Hermite_2v-273 Evaluation symm   ", 	XValues.get(2, 7),   XAValues.get(7, 2),     1e-4);
	check("Hermite_2v-274 Evaluation symm   ", 	XValues.get(3, 7),   XAValues.get(7, 3),     1e-4);
	check("Hermite_2v-275 Evaluation symm   ", 	XValues.get(4, 7),   XAValues.get(7, 4),     1e-4);
	check("Hermite_2v-276 Evaluation symm   ", 	XValues.get(5, 7),   XAValues.get(7, 5),     1e-4);
	check("Hermite_2v-277 Evaluation symm   ", 	XValues.get(6, 7),   XAValues.get(7, 6),     1e-4);

	check("Hermite_2v-281 Evaluation symm   ", 	XValues.get(0, 8),   XAValues.get(8, 0),     1e-4);
	check("Hermite_2v-282 Evaluation symm   ", 	XValues.get(1, 8),   XAValues.get(8, 1),     1e-4);
	check("Hermite_2v-283 Evaluation symm   ", 	XValues.get(2, 8),   XAValues.get(8, 2),     1e-4);
	check("Hermite_2v-284 Evaluation symm   ", 	XValues.get(3, 8),   XAValues.get(8, 3),     1e-4);
	check("Hermite_2v-285 Evaluation symm   ", 	XValues.get(4, 8),   XAValues.get(8, 4),     1e-4);
	check("Hermite_2v-286 Evaluation symm   ", 	XValues.get(5, 8),   XAValues.get(8, 5),     1e-4);
	check("Hermite_2v-287 Evaluation symm   ", 	XValues.get(6, 8),   XAValues.get(8, 6),     1e-4);

	check("Hermite_2v-291 Evaluation symm   ", 	XValues.get(0, 9),   XAValues.get(9, 0),     1e-4);
	check("Hermite_2v-292 Evaluation symm   ", 	XValues.get(1, 9),   XAValues.get(9, 1),     1e-4);
	check("Hermite_2v-293 Evaluation symm   ", 	XValues.get(2, 9),   XAValues.get(9, 2),     1e-4);
	check("Hermite_2v-294 Evaluation symm   ", 	XValues.get(3, 9),   XAValues.get(9, 3),     1e-4);
	check("Hermite_2v-295 Evaluation symm   ", 	XValues.get(4, 9),   XAValues.get(9, 4),     1e-4);
	check("Hermite_2v-296 Evaluation symm   ", 	XValues.get(5, 9),   XAValues.get(9, 5),     1e-4);
	check("Hermite_2v-297 Evaluation symm   ", 	XValues.get(6, 9),   XAValues.get(9, 6),     1e-4);

	math::vec2 XDiffs(7, 10);
	double dum1(1.), dum2(1.);

	for (int i = 0; i != XpointsF.size1(); ++i) {
		for (int j = 0; j != XpointsE.size1(); ++j) {
            XDiffs.get(i, j) = Otab.d_dt(XpointsE[j], XpointsF[i], dum1, dum2);
		}
	}


    Otab.get_slopes_d2().get(0, 0); // prueba 7///////////////////////

	// Check 2D table differentials at the nodes
	check("Hermite_2v-301 Differential      ", 	XDiffs.get(1, 1), Otab.get_slopes_d2().get(0, 0) + Otab.get_slopes_d1().get(0, 0),     1e-1);
	check("Hermite_2v-302 Differential      ", 	XDiffs.get(2, 1), Otab.get_slopes_d2().get(0, 1) + Otab.get_slopes_d1().get(1, 0),     1e-1);
	check("Hermite_2v-303 Differential      ", 	XDiffs.get(4, 1), Otab.get_slopes_d2().get(0, 2) + Otab.get_slopes_d1().get(2, 0),     1e-1);
	check("Hermite_2v-304 Differential      ", 	XDiffs.get(5, 1), Otab.get_slopes_d2().get(0, 3) + Otab.get_slopes_d1().get(3, 0),     1e-1);
	check("Hermite_2v-311 Differential      ", 	XDiffs.get(1, 2), Otab.get_slopes_d2().get(1, 0) + Otab.get_slopes_d1().get(0, 1),     1e-0);
	check("Hermite_2v-312 Differential      ", 	XDiffs.get(2, 2), Otab.get_slopes_d2().get(1, 1) + Otab.get_slopes_d1().get(1, 1),     1e-1);
	check("Hermite_2v-313 Differential      ", 	XDiffs.get(4, 2), Otab.get_slopes_d2().get(1, 2) + Otab.get_slopes_d1().get(2, 1),     1e-1);
	check("Hermite_2v-314 Differential      ", 	XDiffs.get(5, 2), Otab.get_slopes_d2().get(1, 3) + Otab.get_slopes_d1().get(3, 1),     1e-0);
	check("Hermite_2v-321 Differential      ", 	XDiffs.get(1, 3), Otab.get_slopes_d2().get(2, 0) + Otab.get_slopes_d1().get(0, 2),     1e-0);
	check("Hermite_2v-322 Differential      ", 	XDiffs.get(2, 3), Otab.get_slopes_d2().get(2, 1) + Otab.get_slopes_d1().get(1, 2),     1e-1);
	check("Hermite_2v-323 Differential      ", 	XDiffs.get(4, 3), Otab.get_slopes_d2().get(2, 2) + Otab.get_slopes_d1().get(2, 2),     1e-1);
	check("Hermite_2v-324 Differential      ", 	XDiffs.get(5, 3), Otab.get_slopes_d2().get(2, 3) + Otab.get_slopes_d1().get(3, 2),     1e-1);
	check("Hermite_2v-331 Differential      ", 	XDiffs.get(1, 5), Otab.get_slopes_d2().get(3, 0) + Otab.get_slopes_d1().get(0, 3),     1e-1);
	check("Hermite_2v-332 Differential      ", 	XDiffs.get(2, 5), Otab.get_slopes_d2().get(3, 1) + Otab.get_slopes_d1().get(1, 3),     1e-1);
	check("Hermite_2v-333 Differential      ", 	XDiffs.get(4, 5), Otab.get_slopes_d2().get(3, 2) + Otab.get_slopes_d1().get(2, 3),     1e-1);
	check("Hermite_2v-334 Differential      ", 	XDiffs.get(5, 5), Otab.get_slopes_d2().get(3, 3) + Otab.get_slopes_d1().get(3, 3),     1e-1);
	check("Hermite_2v-341 Differential      ", 	XDiffs.get(1, 6), Otab.get_slopes_d2().get(4, 0) + Otab.get_slopes_d1().get(0, 4),     1e-1);
	check("Hermite_2v-342 Differential      ", 	XDiffs.get(2, 6), Otab.get_slopes_d2().get(4, 1) + Otab.get_slopes_d1().get(1, 4),     1e-1);
	check("Hermite_2v-343 Differential      ", 	XDiffs.get(4, 6), Otab.get_slopes_d2().get(4, 2) + Otab.get_slopes_d1().get(2, 4),     1e-1);
	check("Hermite_2v-344 Differential      ", 	XDiffs.get(5, 6), Otab.get_slopes_d2().get(4, 3) + Otab.get_slopes_d1().get(3, 4),     1e-1);
	check("Hermite_2v-351 Differential      ", 	XDiffs.get(1, 8), Otab.get_slopes_d2().get(5, 0) + Otab.get_slopes_d1().get(0, 5),     1e-1);
	check("Hermite_2v-352 Differential      ", 	XDiffs.get(2, 8), Otab.get_slopes_d2().get(5, 1) + Otab.get_slopes_d1().get(1, 5),     1e-1);
	check("Hermite_2v-353 Differential      ", 	XDiffs.get(4, 8), Otab.get_slopes_d2().get(5, 2) + Otab.get_slopes_d1().get(2, 5),     1e-1);
	check("Hermite_2v-354 Differential      ", 	XDiffs.get(5, 8), Otab.get_slopes_d2().get(5, 3) + Otab.get_slopes_d1().get(3, 5),     1e-1);

	math::vec2 XADiffs(7, 10);

	for (int i = 0; i != XpointsF.size1(); ++i) {
		for (int j = 0; j != XpointsE.size1(); ++j) {
            XADiffs.get(i, j) = OAtab.d_dt(XpointsF[i], XpointsE[j], dum1, dum2);
		}
	}

	// Check 2D simmetrical (or traspose) table differentials at the nodes
	check("Hermite_2v-401 Differential trasp", 	XADiffs.get(1, 1), OAtab.get_slopes_d2().get(0, 0) + OAtab.get_slopes_d1().get(0, 0),     1e-1);
	check("Hermite_2v-402 Differential trasp", 	XADiffs.get(2, 1), OAtab.get_slopes_d2().get(1, 0) + OAtab.get_slopes_d1().get(0, 1),     1e-1);
	check("Hermite_2v-403 Differential trasp", 	XADiffs.get(4, 1), OAtab.get_slopes_d2().get(2, 0) + OAtab.get_slopes_d1().get(0, 2),     1e-1);
	check("Hermite_2v-404 Differential trasp", 	XADiffs.get(5, 1), OAtab.get_slopes_d2().get(3, 0) + OAtab.get_slopes_d1().get(0, 3),     1e-1);
	check("Hermite_2v-411 Differential trasp", 	XADiffs.get(1, 2), OAtab.get_slopes_d2().get(0, 1) + OAtab.get_slopes_d1().get(1, 0),     1e-0);
	check("Hermite_2v-412 Differential trasp", 	XADiffs.get(2, 2), OAtab.get_slopes_d2().get(1, 1) + OAtab.get_slopes_d1().get(1, 1),     1e-1);
	check("Hermite_2v-413 Differential trasp", 	XADiffs.get(4, 2), OAtab.get_slopes_d2().get(2, 1) + OAtab.get_slopes_d1().get(1, 2),     1e-1);
	check("Hermite_2v-414 Differential trasp", 	XADiffs.get(5, 2), OAtab.get_slopes_d2().get(3, 1) + OAtab.get_slopes_d1().get(1, 3),     1e-0);
	check("Hermite_2v-421 Differential trasp", 	XADiffs.get(1, 3), OAtab.get_slopes_d2().get(0, 2) + OAtab.get_slopes_d1().get(2, 0),     1e-0);
	check("Hermite_2v-422 Differential trasp", 	XADiffs.get(2, 3), OAtab.get_slopes_d2().get(1, 2) + OAtab.get_slopes_d1().get(2, 1),     1e-1);
	check("Hermite_2v-423 Differential trasp", 	XADiffs.get(4, 3), OAtab.get_slopes_d2().get(2, 2) + OAtab.get_slopes_d1().get(2, 2),     1e-1);
	check("Hermite_2v-424 Differential trasp", 	XADiffs.get(5, 3), OAtab.get_slopes_d2().get(3, 2) + OAtab.get_slopes_d1().get(2, 3),     1e-1);
	check("Hermite_2v-431 Differential trasp", 	XADiffs.get(1, 5), OAtab.get_slopes_d2().get(0, 3) + OAtab.get_slopes_d1().get(3, 0),     1e-1);
	check("Hermite_2v-432 Differential trasp", 	XADiffs.get(2, 5), OAtab.get_slopes_d2().get(1, 3) + OAtab.get_slopes_d1().get(3, 1),     1e-1);
	check("Hermite_2v-433 Differential trasp", 	XADiffs.get(4, 5), OAtab.get_slopes_d2().get(2, 3) + OAtab.get_slopes_d1().get(3, 2),     1e-1);
	check("Hermite_2v-434 Differential trasp", 	XADiffs.get(5, 5), OAtab.get_slopes_d2().get(3, 3) + OAtab.get_slopes_d1().get(3, 3),     1e-1);
	check("Hermite_2v-441 Differential trasp", 	XADiffs.get(1, 6), OAtab.get_slopes_d2().get(0, 4) + OAtab.get_slopes_d1().get(4, 0),     1e-1);
	check("Hermite_2v-442 Differential trasp", 	XADiffs.get(2, 6), OAtab.get_slopes_d2().get(1, 4) + OAtab.get_slopes_d1().get(4, 1),     1e-1);
	check("Hermite_2v-443 Differential trasp", 	XADiffs.get(4, 6), OAtab.get_slopes_d2().get(2, 4) + OAtab.get_slopes_d1().get(4, 2),     1e-1);
	check("Hermite_2v-444 Differential trasp", 	XADiffs.get(5, 6), OAtab.get_slopes_d2().get(3, 4) + OAtab.get_slopes_d1().get(4, 3),     1e-1);
	check("Hermite_2v-451 Differential trasp", 	XADiffs.get(1, 8), OAtab.get_slopes_d2().get(0, 5) + OAtab.get_slopes_d1().get(5, 0),     1e-1);
	check("Hermite_2v-452 Differential trasp", 	XADiffs.get(2, 8), OAtab.get_slopes_d2().get(1, 5) + OAtab.get_slopes_d1().get(5, 1),     1e-1);
	check("Hermite_2v-453 Differential trasp", 	XADiffs.get(4, 8), OAtab.get_slopes_d2().get(2, 5) + OAtab.get_slopes_d1().get(5, 2),     1e-1);
	check("Hermite_2v-454 Differential trasp", 	XADiffs.get(5, 8), OAtab.get_slopes_d2().get(3, 5) + OAtab.get_slopes_d1().get(5, 3),     1e-1);

	// Compare 2D table differentials with its simmetrical both at the nodes and at other values
	check("Hermite_2v-501 Differential symm ", 	XDiffs.get(0, 0),   XADiffs.get(0, 0),     1e-4);
	check("Hermite_2v-502 Differential symm ", 	XDiffs.get(0, 1),   XADiffs.get(0, 1),     1e-4);
	check("Hermite_2v-503 Differential symm ", 	XDiffs.get(0, 2),   XADiffs.get(0, 2),     1e-4);
	check("Hermite_2v-504 Differential symm ", 	XDiffs.get(0, 3),   XADiffs.get(0, 3),     1e-4);
	check("Hermite_2v-505 Differential symm ", 	XDiffs.get(0, 4),   XADiffs.get(0, 4),     1e-4);
	check("Hermite_2v-506 Differential symm ", 	XDiffs.get(0, 5),   XADiffs.get(0, 5),     1e-4);
	check("Hermite_2v-507 Differential symm ", 	XDiffs.get(0, 6),   XADiffs.get(0, 6),     1e-4);
	check("Hermite_2v-508 Differential symm ", 	XDiffs.get(0, 7),   XADiffs.get(0, 7),     1e-4);
	check("Hermite_2v-509 Differential symm ", 	XDiffs.get(0, 8),   XADiffs.get(0, 8),     1e-4);
	check("Hermite_2v-510 Differential symm ", 	XDiffs.get(0, 9),   XADiffs.get(0, 9),     1e-4);

	check("Hermite_2v-511 Differential symm ", 	XDiffs.get(1, 0),   XADiffs.get(1, 0),     1e-4);
	check("Hermite_2v-512 Differential symm ", 	XDiffs.get(1, 1),   XADiffs.get(1, 1),     1e-4);
	check("Hermite_2v-513 Differential symm ", 	XDiffs.get(1, 2),   XADiffs.get(1, 2),     1e-4);
	check("Hermite_2v-514 Differential symm ", 	XDiffs.get(1, 3),   XADiffs.get(1, 3),     1e-4);
	check("Hermite_2v-515 Differential symm ", 	XDiffs.get(1, 4),   XADiffs.get(1, 4),     1e-4);
	check("Hermite_2v-516 Differential symm ", 	XDiffs.get(1, 5),   XADiffs.get(1, 5),     1e-4);
	check("Hermite_2v-517 Differential symm ", 	XDiffs.get(1, 6),   XADiffs.get(1, 6),     1e-4);
	check("Hermite_2v-518 Differential symm ", 	XDiffs.get(1, 7),   XADiffs.get(1, 7),     1e-4);
	check("Hermite_2v-519 Differential symm ", 	XDiffs.get(1, 8),   XADiffs.get(1, 8),     1e-4);
	check("Hermite_2v-520 Differential symm ", 	XDiffs.get(1, 9),   XADiffs.get(1, 9),     1e-4);

	check("Hermite_2v-521 Differential symm ", 	XDiffs.get(2, 0),   XADiffs.get(2, 0),     1e-4);
	check("Hermite_2v-522 Differential symm ", 	XDiffs.get(2, 1),   XADiffs.get(2, 1),     1e-4);
	check("Hermite_2v-523 Differential symm ", 	XDiffs.get(2, 2),   XADiffs.get(2, 2),     1e-4);
	check("Hermite_2v-524 Differential symm ", 	XDiffs.get(2, 3),   XADiffs.get(2, 3),     1e-4);
	check("Hermite_2v-525 Differential symm ", 	XDiffs.get(2, 4),   XADiffs.get(2, 4),     1e-4);
	check("Hermite_2v-526 Differential symm ", 	XDiffs.get(2, 5),   XADiffs.get(2, 5),     1e-4);
	check("Hermite_2v-527 Differential symm ", 	XDiffs.get(2, 6),   XADiffs.get(2, 6),     1e-4);
	check("Hermite_2v-528 Differential symm ", 	XDiffs.get(2, 7),   XADiffs.get(2, 7),     1e-4);
	check("Hermite_2v-529 Differential symm ", 	XDiffs.get(2, 8),   XADiffs.get(2, 8),     1e-4);
	check("Hermite_2v-530 Differential symm ", 	XDiffs.get(2, 9),   XADiffs.get(2, 9),     1e-4);

	check("Hermite_2v-531 Differential symm ", 	XDiffs.get(3, 0),   XADiffs.get(3, 0),     1e-4);
	check("Hermite_2v-532 Differential symm ", 	XDiffs.get(3, 1),   XADiffs.get(3, 1),     1e-4);
	check("Hermite_2v-533 Differential symm ", 	XDiffs.get(3, 2),   XADiffs.get(3, 2),     1e-4);
	check("Hermite_2v-534 Differential symm ", 	XDiffs.get(3, 3),   XADiffs.get(3, 3),     1e-4);
	check("Hermite_2v-535 Differential symm ", 	XDiffs.get(3, 4),   XADiffs.get(3, 4),     1e-4);
	check("Hermite_2v-536 Differential symm ", 	XDiffs.get(3, 5),   XADiffs.get(3, 5),     1e-4);
	check("Hermite_2v-537 Differential symm ", 	XDiffs.get(3, 6),   XADiffs.get(3, 6),     1e-4);
	check("Hermite_2v-538 Differential symm ", 	XDiffs.get(3, 7),   XADiffs.get(3, 7),     1e-4);
	check("Hermite_2v-539 Differential symm ", 	XDiffs.get(3, 8),   XADiffs.get(3, 8),     1e-4);
	check("Hermite_2v-540 Differential symm ", 	XDiffs.get(3, 9),   XADiffs.get(3, 9),     1e-4);

	check("Hermite_2v-541 Differential symm ", 	XDiffs.get(4, 0),   XADiffs.get(4, 0),     1e-4);
	check("Hermite_2v-542 Differential symm ", 	XDiffs.get(4, 1),   XADiffs.get(4, 1),     1e-4);
	check("Hermite_2v-543 Differential symm ", 	XDiffs.get(4, 2),   XADiffs.get(4, 2),     1e-4);
	check("Hermite_2v-544 Differential symm ", 	XDiffs.get(4, 3),   XADiffs.get(4, 3),     1e-4);
	check("Hermite_2v-545 Differential symm ", 	XDiffs.get(4, 4),   XADiffs.get(4, 4),     1e-4);
	check("Hermite_2v-546 Differential symm ", 	XDiffs.get(4, 5),   XADiffs.get(4, 5),     1e-4);
	check("Hermite_2v-547 Differential symm ", 	XDiffs.get(4, 6),   XADiffs.get(4, 6),     1e-4);
	check("Hermite_2v-548 Differential symm ", 	XDiffs.get(4, 7),   XADiffs.get(4, 7),     1e-4);
	check("Hermite_2v-549 Differential symm ", 	XDiffs.get(4, 8),   XADiffs.get(4, 8),     1e-4);
	check("Hermite_2v-550 Differential symm ", 	XDiffs.get(4, 9),   XADiffs.get(4, 9),     1e-4);

	check("Hermite_2v-551 Differential symm ", 	XDiffs.get(5, 0),   XADiffs.get(5, 0),     1e-4);
	check("Hermite_2v-552 Differential symm ", 	XDiffs.get(5, 1),   XADiffs.get(5, 1),     1e-4);
	check("Hermite_2v-553 Differential symm ", 	XDiffs.get(5, 2),   XADiffs.get(5, 2),     1e-4);
	check("Hermite_2v-554 Differential symm ", 	XDiffs.get(5, 3),   XADiffs.get(5, 3),     1e-4);
	check("Hermite_2v-555 Differential symm ", 	XDiffs.get(5, 4),   XADiffs.get(5, 4),     1e-4);
	check("Hermite_2v-556 Differential symm ", 	XDiffs.get(5, 5),   XADiffs.get(5, 5),     1e-4);
	check("Hermite_2v-557 Differential symm ", 	XDiffs.get(5, 6),   XADiffs.get(5, 6),     1e-4);
	check("Hermite_2v-558 Differential symm ", 	XDiffs.get(5, 7),   XADiffs.get(5, 7),     1e-4);
	check("Hermite_2v-559 Differential symm ", 	XDiffs.get(5, 8),   XADiffs.get(5, 8),     1e-4);
	check("Hermite_2v-560 Differential symm ", 	XDiffs.get(5, 9),   XADiffs.get(5, 9),     1e-4);

	check("Hermite_2v-561 Differential symm ", 	XDiffs.get(6, 0),   XADiffs.get(6, 0),     1e-4);
	check("Hermite_2v-562 Differential symm ", 	XDiffs.get(6, 1),   XADiffs.get(6, 1),     1e-4);
	check("Hermite_2v-563 Differential symm ", 	XDiffs.get(6, 2),   XADiffs.get(6, 2),     1e-4);
	check("Hermite_2v-564 Differential symm ", 	XDiffs.get(6, 3),   XADiffs.get(6, 3),     1e-4);
	check("Hermite_2v-565 Differential symm ", 	XDiffs.get(6, 4),   XADiffs.get(6, 4),     1e-4);
	check("Hermite_2v-566 Differential symm ", 	XDiffs.get(6, 5),   XADiffs.get(6, 5),     1e-4);
	check("Hermite_2v-567 Differential symm ", 	XDiffs.get(6, 6),   XADiffs.get(6, 6),     1e-4);
	check("Hermite_2v-568 Differential symm ", 	XDiffs.get(6, 7),   XADiffs.get(6, 7),     1e-4);
	check("Hermite_2v-569 Differential symm ", 	XDiffs.get(6, 8),   XADiffs.get(6, 8),     1e-4);
	check("Hermite_2v-570 Differential symm ", 	XDiffs.get(6, 9),   XADiffs.get(6, 9),     1e-4);

	double stop_here = 88;

} // closes test_hermite_2d

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void math::test::Tinterp::test_hermite_3d(math::logic::INTERP_MODE interp_mode) {
	math::vec1 pointsE(3);
	pointsE.set(0, 5.);
	pointsE.set(1, 20.);
	pointsE.set(2, 50.);
	math::vec1 pointsF(6);
	pointsF.set(0, 2.0);
	pointsF.set(1, 2.5);
	pointsF.set(2, 3.5);
	pointsF.set(3, 5.5);
	pointsF.set(4, 7.5);
	pointsF.set(5, 10.5);
	math::vec1 pointsG(4);
	pointsG.set(0, 50.0);
	pointsG.set(1, 54.5);
	pointsG.set(2, 65.5);
	pointsG.set(3, 72.0);

	math::vec3 VValues(4, 6, 3);
	VValues.set(0, 0, 0, 10.);	VValues.set(1, 0, 0, 15.);	VValues.set(2, 0, 0, 12.);	VValues.set(3, 0, 0, 22.);
	VValues.set(0, 1, 0, 2.);	VValues.set(1, 1, 0, 12.);	VValues.set(2, 1, 0, 17.);	VValues.set(3, 1, 0, 12.);
	VValues.set(0, 2, 0, 84.);	VValues.set(1, 2, 0, 5.);	VValues.set(2, 2, 0, 34.);	VValues.set(3, 2, 0, 72.);
	VValues.set(0, 3, 0, 32.);	VValues.set(1, 3, 0, 55.);	VValues.set(2, 3, 0, 66.);	VValues.set(3, 3, 0, 41.);
	VValues.set(0, 4, 0, 51.);	VValues.set(1, 4, 0, 67.);	VValues.set(2, 4, 0, 5.);	VValues.set(3, 4, 0, 34.);
	VValues.set(0, 5, 0, 11.);	VValues.set(1, 5, 0, 72.);	VValues.set(2, 5, 0, 88.);	VValues.set(3, 5, 0, 64.);
	VValues.set(0, 0, 1, 60.);	VValues.set(1, 0, 1, 25.);	VValues.set(2, 0, 1, 2.);	VValues.set(3, 0, 1, 42.);
	VValues.set(0, 1, 1, 42.);	VValues.set(1, 1, 1, 2.);	VValues.set(2, 1, 1, 37.);	VValues.set(3, 1, 1, 24.);
	VValues.set(0, 2, 1, 4.);	VValues.set(1, 2, 1, 54.);	VValues.set(2, 2, 1, 4.);	VValues.set(3, 2, 1, 33.);
	VValues.set(0, 3, 1, 37.);	VValues.set(1, 3, 1, 25.);	VValues.set(2, 3, 1, 69.);	VValues.set(3, 3, 1, 4.);
	VValues.set(0, 4, 1, 54.);	VValues.set(1, 4, 1, 61.);	VValues.set(2, 4, 1, 42.);	VValues.set(3, 4, 1, 42.);
	VValues.set(0, 5, 1, 41.);	VValues.set(1, 5, 1, 22.);	VValues.set(2, 5, 1, 28.);	VValues.set(3, 5, 1, 32.);
	VValues.set(0, 0, 2, 78.);	VValues.set(1, 0, 2, 52.);	VValues.set(2, 0, 2, 72.);	VValues.set(3, 0, 2, 33.);
	VValues.set(0, 1, 2, 52.);	VValues.set(1, 1, 2, 12.);	VValues.set(2, 1, 2, 7.);	VValues.set(3, 1, 2, 82.);
	VValues.set(0, 2, 2, 41.);	VValues.set(1, 2, 2, 58);	VValues.set(2, 2, 2, 14.);	VValues.set(3, 2, 2, 22.);
	VValues.set(0, 3, 2, 72.);	VValues.set(1, 3, 2, 15.);	VValues.set(2, 3, 2, 13.);	VValues.set(3, 3, 2, 14.);
	VValues.set(0, 4, 2, 11.);	VValues.set(1, 4, 2, 31.);	VValues.set(2, 4, 2, 21.);	VValues.set(3, 4, 2, 1.);
	VValues.set(0, 5, 2, 1.);	VValues.set(1, 5, 2, 29.);	VValues.set(2, 5, 2, 78.);	VValues.set(3, 5, 2, 33.);
	math::f_table3V Otab(pointsE, pointsF, pointsG, VValues, interp_mode);

	// Make simmetric copy to check results
	math::vec1 ApointsE(pointsG);
	math::vec1 ApointsF(pointsF);
	math::vec1 ApointsG(pointsE);

	math::vec3 AVValues(3, 6, 4);
	AVValues.set(0, 0, 0, 10.);	AVValues.set(1, 0, 0, 60.);	AVValues.set(2, 0, 0, 78.);
	AVValues.set(0, 1, 0, 2.);	AVValues.set(1, 1, 0, 42.);	AVValues.set(2, 1, 0, 52.);
	AVValues.set(0, 2, 0, 84.);	AVValues.set(1, 2, 0, 4.);	AVValues.set(2, 2, 0, 41.);
	AVValues.set(0, 3, 0, 32.);	AVValues.set(1, 3, 0, 37.);	AVValues.set(2, 3, 0, 72.);
	AVValues.set(0, 4, 0, 51.);	AVValues.set(1, 4, 0, 54.);	AVValues.set(2, 4, 0, 11.);	
	AVValues.set(0, 5, 0, 11.);	AVValues.set(1, 5, 0, 41.);	AVValues.set(2, 5, 0, 1.);
	AVValues.set(0, 0, 1, 15.);	AVValues.set(1, 0, 1, 25.);	AVValues.set(2, 0, 1, 52.);
	AVValues.set(0, 1, 1, 12.);	AVValues.set(1, 1, 1, 2.);	AVValues.set(2, 1, 1, 12.);
	AVValues.set(0, 2, 1, 5.);	AVValues.set(1, 2, 1, 54.);	AVValues.set(2, 2, 1, 58.);
	AVValues.set(0, 3, 1, 55.);	AVValues.set(1, 3, 1, 25.);	AVValues.set(2, 3, 1, 15.);
	AVValues.set(0, 4, 1, 67.);	AVValues.set(1, 4, 1, 61.);	AVValues.set(2, 4, 1, 31.);	
	AVValues.set(0, 5, 1, 72.);	AVValues.set(1, 5, 1, 22.);	AVValues.set(2, 5, 1, 29.);
	AVValues.set(0, 0, 2, 12.);	AVValues.set(1, 0, 2, 2.);	AVValues.set(2, 0, 2, 72.);
	AVValues.set(0, 1, 2, 17.);	AVValues.set(1, 1, 2, 37.);	AVValues.set(2, 1, 2, 7.);
	AVValues.set(0, 2, 2, 34.);	AVValues.set(1, 2, 2, 4.);	AVValues.set(2, 2, 2, 14.);
	AVValues.set(0, 3, 2, 66.);	AVValues.set(1, 3, 2, 69.);	AVValues.set(2, 3, 2, 13.);
	AVValues.set(0, 4, 2, 5.);	AVValues.set(1, 4, 2, 42.);	AVValues.set(2, 4, 2, 21.);	
	AVValues.set(0, 5, 2, 88.);	AVValues.set(1, 5, 2, 28.);	AVValues.set(2, 5, 2, 78.);
	AVValues.set(0, 0, 3, 22.);	AVValues.set(1, 0, 3, 42.);	AVValues.set(2, 0, 3, 33.);
	AVValues.set(0, 1, 3, 12.);	AVValues.set(1, 1, 3, 24.);	AVValues.set(2, 1, 3, 82.);
	AVValues.set(0, 2, 3, 72.);	AVValues.set(1, 2, 3, 33.);	AVValues.set(2, 2, 3, 22.);
	AVValues.set(0, 3, 3, 41.);	AVValues.set(1, 3, 3, 4.);	AVValues.set(2, 3, 3, 14.);
	AVValues.set(0, 4, 3, 34.);	AVValues.set(1, 4, 3, 42.);	AVValues.set(2, 4, 3, 1.);	
	AVValues.set(0, 5, 3, 64.);	AVValues.set(1, 5, 3, 32.);	AVValues.set(2, 5, 3, 33.);
	math::f_table3V OAtab(ApointsE, ApointsF, ApointsG, AVValues, interp_mode);

	math::vec1 XpointsE(5);
	XpointsE.set(0, 5.0);
	XpointsE.set(1, 10.0);
	XpointsE.set(2, 20.0);
	XpointsE.set(3, 35.0);
	XpointsE.set(4, 50.0);
	
	math::vec1 XpointsF(9);
	XpointsF.set(0, 2.0);
	XpointsF.set(1, 2.2);
	XpointsF.set(2, 2.5);
	XpointsF.set(3, 3.5);
	XpointsF.set(4, 4.0);
	XpointsF.set(5, 5.5);
	XpointsF.set(6, 7.5);
	XpointsF.set(7, 9.2);
	XpointsF.set(8, 10.5);

	math::vec1 XpointsG(6);
	XpointsG.set(0, 50.0);
	XpointsG.set(1, 52.0);
	XpointsG.set(2, 54.5);
	XpointsG.set(3, 59.2);
	XpointsG.set(4, 65.5);
	XpointsG.set(5, 72.0);;

	math::vec3 XVValues(6, 9, 5);

	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
                XVValues.get(k, j, i) = Otab.value(XpointsE.get(i), XpointsF.get(j), XpointsG.get(k));
			}
		}
	}

	// Check 3D table evaluation at the nodes
	check("Hermite_3v-001 Evaluation        ", 	XVValues.get(0, 0, 0),   VValues.get(0, 0, 0),     1e-4);
	check("Hermite_3v-002 Evaluation        ", 	XVValues.get(2, 0, 0),   VValues.get(1, 0, 0),     1e-4);
	check("Hermite_3v-003 Evaluation        ", 	XVValues.get(4, 0, 0),   VValues.get(2, 0, 0),     1e-4);
	check("Hermite_3v-004 Evaluation        ", 	XVValues.get(5, 0, 0),   VValues.get(3, 0, 0),     1e-4);
	check("Hermite_3v-005 Evaluation        ", 	XVValues.get(0, 2, 0),   VValues.get(0, 1, 0),     1e-4);
	check("Hermite_3v-006 Evaluation        ", 	XVValues.get(2, 2, 0),   VValues.get(1, 1, 0),     1e-4);
	check("Hermite_3v-007 Evaluation        ", 	XVValues.get(4, 2, 0),   VValues.get(2, 1, 0),     1e-4);
	check("Hermite_3v-008 Evaluation        ", 	XVValues.get(5, 2, 0),   VValues.get(3, 1, 0),     1e-4);
	check("Hermite_3v-009 Evaluation        ", 	XVValues.get(0, 3, 0),   VValues.get(0, 2, 0),     1e-4);
	check("Hermite_3v-010 Evaluation        ", 	XVValues.get(2, 3, 0),   VValues.get(1, 2, 0),     1e-4);
	check("Hermite_3v-011 Evaluation        ", 	XVValues.get(4, 3, 0),   VValues.get(2, 2, 0),     1e-4);
	check("Hermite_3v-012 Evaluation        ", 	XVValues.get(5, 3, 0),   VValues.get(3, 2, 0),     1e-4);
	check("Hermite_3v-013 Evaluation        ", 	XVValues.get(0, 5, 0),   VValues.get(0, 3, 0),     1e-4);
	check("Hermite_3v-014 Evaluation        ", 	XVValues.get(2, 5, 0),   VValues.get(1, 3, 0),     1e-4);
	check("Hermite_3v-015 Evaluation        ", 	XVValues.get(4, 5, 0),   VValues.get(2, 3, 0),     1e-4);
	check("Hermite_3v-016 Evaluation        ", 	XVValues.get(5, 5, 0),   VValues.get(3, 3, 0),     1e-4);
	check("Hermite_3v-017 Evaluation        ", 	XVValues.get(0, 6, 0),   VValues.get(0, 4, 0),     1e-4);
	check("Hermite_3v-018 Evaluation        ", 	XVValues.get(2, 6, 0),   VValues.get(1, 4, 0),     1e-4);
	check("Hermite_3v-019 Evaluation        ", 	XVValues.get(4, 6, 0),   VValues.get(2, 4, 0),     1e-4);
	check("Hermite_3v-020 Evaluation        ", 	XVValues.get(5, 6, 0),   VValues.get(3, 4, 0),     1e-4);
	check("Hermite_3v-021 Evaluation        ", 	XVValues.get(0, 8, 0),   VValues.get(0, 5, 0),     1e-4);
	check("Hermite_3v-022 Evaluation        ", 	XVValues.get(2, 8, 0),   VValues.get(1, 5, 0),     1e-4);
	check("Hermite_3v-023 Evaluation        ", 	XVValues.get(4, 8, 0),   VValues.get(2, 5, 0),     1e-4);
	check("Hermite_3v-024 Evaluation        ", 	XVValues.get(5, 8, 0),   VValues.get(3, 5, 0),     1e-4);

	check("Hermite_3v-031 Evaluation        ", 	XVValues.get(0, 0, 2),   VValues.get(0, 0, 1),     1e-4);
	check("Hermite_3v-032 Evaluation        ", 	XVValues.get(2, 0, 2),   VValues.get(1, 0, 1),     1e-4);
	check("Hermite_3v-033 Evaluation        ", 	XVValues.get(4, 0, 2),   VValues.get(2, 0, 1),     1e-4);
	check("Hermite_3v-034 Evaluation        ", 	XVValues.get(5, 0, 2),   VValues.get(3, 0, 1),     1e-4);
	check("Hermite_3v-035 Evaluation        ", 	XVValues.get(0, 2, 2),   VValues.get(0, 1, 1),     1e-4);
	check("Hermite_3v-036 Evaluation        ", 	XVValues.get(2, 2, 2),   VValues.get(1, 1, 1),     1e-4);
	check("Hermite_3v-037 Evaluation        ", 	XVValues.get(4, 2, 2),   VValues.get(2, 1, 1),     1e-4);
	check("Hermite_3v-038 Evaluation        ", 	XVValues.get(5, 2, 2),   VValues.get(3, 1, 1),     1e-4);
	check("Hermite_3v-039 Evaluation        ", 	XVValues.get(0, 3, 2),   VValues.get(0, 2, 1),     1e-4);
	check("Hermite_3v-040 Evaluation        ", 	XVValues.get(2, 3, 2),   VValues.get(1, 2, 1),     1e-4);
	check("Hermite_3v-041 Evaluation        ", 	XVValues.get(4, 3, 2),   VValues.get(2, 2, 1),     1e-4);
	check("Hermite_3v-042 Evaluation        ", 	XVValues.get(5, 3, 2),   VValues.get(3, 2, 1),     1e-4);
	check("Hermite_3v-043 Evaluation        ", 	XVValues.get(0, 5, 2),   VValues.get(0, 3, 1),     1e-4);
	check("Hermite_3v-044 Evaluation        ", 	XVValues.get(2, 5, 2),   VValues.get(1, 3, 1),     1e-4);
	check("Hermite_3v-045 Evaluation        ", 	XVValues.get(4, 5, 2),   VValues.get(2, 3, 1),     1e-4);
	check("Hermite_3v-046 Evaluation        ", 	XVValues.get(5, 5, 2),   VValues.get(3, 3, 1),     1e-4);
	check("Hermite_3v-047 Evaluation        ", 	XVValues.get(0, 6, 2),   VValues.get(0, 4, 1),     1e-4);
	check("Hermite_3v-048 Evaluation        ", 	XVValues.get(2, 6, 2),   VValues.get(1, 4, 1),     1e-4);
	check("Hermite_3v-049 Evaluation        ", 	XVValues.get(4, 6, 2),   VValues.get(2, 4, 1),     1e-4);
	check("Hermite_3v-050 Evaluation        ", 	XVValues.get(5, 6, 2),   VValues.get(3, 4, 1),     1e-4);
	check("Hermite_3v-051 Evaluation        ", 	XVValues.get(0, 8, 2),   VValues.get(0, 5, 1),     1e-4);
	check("Hermite_3v-052 Evaluation        ", 	XVValues.get(2, 8, 2),   VValues.get(1, 5, 1),     1e-4);
	check("Hermite_3v-053 Evaluation        ", 	XVValues.get(4, 8, 2),   VValues.get(2, 5, 1),     1e-4);
	check("Hermite_3v-054 Evaluation        ", 	XVValues.get(5, 8, 2),   VValues.get(3, 5, 1),     1e-4);

	check("Hermite_3v-061 Evaluation        ", 	XVValues.get(0, 0, 4),   VValues.get(0, 0, 2),     1e-4);
	check("Hermite_3v-062 Evaluation        ", 	XVValues.get(2, 0, 4),   VValues.get(1, 0, 2),     1e-4);
	check("Hermite_3v-063 Evaluation        ", 	XVValues.get(4, 0, 4),   VValues.get(2, 0, 2),     1e-4);
	check("Hermite_3v-064 Evaluation        ", 	XVValues.get(5, 0, 4),   VValues.get(3, 0, 2),     1e-4);
	check("Hermite_3v-065 Evaluation        ", 	XVValues.get(0, 2, 4),   VValues.get(0, 1, 2),     1e-4);
	check("Hermite_3v-066 Evaluation        ", 	XVValues.get(2, 2, 4),   VValues.get(1, 1, 2),     1e-4);
	check("Hermite_3v-067 Evaluation        ", 	XVValues.get(4, 2, 4),   VValues.get(2, 1, 2),     1e-4);
	check("Hermite_3v-068 Evaluation        ", 	XVValues.get(5, 2, 4),   VValues.get(3, 1, 2),     1e-4);
	check("Hermite_3v-069 Evaluation        ", 	XVValues.get(0, 3, 4),   VValues.get(0, 2, 2),     1e-4);
	check("Hermite_3v-070 Evaluation        ", 	XVValues.get(2, 3, 4),   VValues.get(1, 2, 2),     1e-4);
	check("Hermite_3v-071 Evaluation        ", 	XVValues.get(4, 3, 4),   VValues.get(2, 2, 2),     1e-4);
	check("Hermite_3v-072 Evaluation        ", 	XVValues.get(5, 3, 4),   VValues.get(3, 2, 2),     1e-4);
	check("Hermite_3v-073 Evaluation        ", 	XVValues.get(0, 5, 4),   VValues.get(0, 3, 2),     1e-4);
	check("Hermite_3v-074 Evaluation        ", 	XVValues.get(2, 5, 4),   VValues.get(1, 3, 2),     1e-4);
	check("Hermite_3v-075 Evaluation        ", 	XVValues.get(4, 5, 4),   VValues.get(2, 3, 2),     1e-4);
	check("Hermite_3v-076 Evaluation        ", 	XVValues.get(5, 5, 4),   VValues.get(3, 3, 2),     1e-4);
	check("Hermite_3v-077 Evaluation        ", 	XVValues.get(0, 6, 4),   VValues.get(0, 4, 2),     1e-4);
	check("Hermite_3v-078 Evaluation        ", 	XVValues.get(2, 6, 4),   VValues.get(1, 4, 2),     1e-4);
	check("Hermite_3v-079 Evaluation        ", 	XVValues.get(4, 6, 4),   VValues.get(2, 4, 2),     1e-4);
	check("Hermite_3v-080 Evaluation        ", 	XVValues.get(5, 6, 4),   VValues.get(3, 4, 2),     1e-4);
	check("Hermite_3v-081 Evaluation        ", 	XVValues.get(0, 8, 4),   VValues.get(0, 5, 2),     1e-4);
	check("Hermite_3v-082 Evaluation        ", 	XVValues.get(2, 8, 4),   VValues.get(1, 5, 2),     1e-4);
	check("Hermite_3v-083 Evaluation        ", 	XVValues.get(4, 8, 4),   VValues.get(2, 5, 2),     1e-4);
	check("Hermite_3v-084 Evaluation        ", 	XVValues.get(5, 8, 4),   VValues.get(3, 5, 2),     1e-4);

	math::vec3 XAVValues(6, 9, 5);

	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
                XAVValues.get(k, j, i) = OAtab.value(XpointsG.get(k), XpointsF.get(j), XpointsE.get(i));
			}
		}
	}

	// Check 3D simmetrical (or traspose) table evaluation at the nodes
	check("Hermite_3v-101 Evaluation trasp  ", 	XAVValues.get(0, 0, 0),   VValues.get(0, 0, 0),     1e-4);
	check("Hermite_3v-102 Evaluation trasp  ", 	XAVValues.get(2, 0, 0),   VValues.get(1, 0, 0),     1e-4);
	check("Hermite_3v-103 Evaluation trasp  ", 	XAVValues.get(4, 0, 0),   VValues.get(2, 0, 0),     1e-4);
	check("Hermite_3v-104 Evaluation trasp  ", 	XAVValues.get(5, 0, 0),   VValues.get(3, 0, 0),     1e-4);
	check("Hermite_3v-105 Evaluation trasp  ", 	XAVValues.get(0, 2, 0),   VValues.get(0, 1, 0),     1e-4);
	check("Hermite_3v-106 Evaluation trasp  ", 	XAVValues.get(2, 2, 0),   VValues.get(1, 1, 0),     1e-4);
	check("Hermite_3v-107 Evaluation trasp  ", 	XAVValues.get(4, 2, 0),   VValues.get(2, 1, 0),     1e-4);
	check("Hermite_3v-108 Evaluation trasp  ", 	XAVValues.get(5, 2, 0),   VValues.get(3, 1, 0),     1e-4);
	check("Hermite_3v-109 Evaluation trasp  ", 	XAVValues.get(0, 3, 0),   VValues.get(0, 2, 0),     1e-4);
	check("Hermite_3v-110 Evaluation trasp  ", 	XAVValues.get(2, 3, 0),   VValues.get(1, 2, 0),     1e-4);
	check("Hermite_3v-111 Evaluation trasp  ", 	XAVValues.get(4, 3, 0),   VValues.get(2, 2, 0),     1e-4);
	check("Hermite_3v-112 Evaluation trasp  ", 	XAVValues.get(5, 3, 0),   VValues.get(3, 2, 0),     1e-4);
	check("Hermite_3v-113 Evaluation trasp  ", 	XAVValues.get(0, 5, 0),   VValues.get(0, 3, 0),     1e-4);
	check("Hermite_3v-114 Evaluation trasp  ", 	XAVValues.get(2, 5, 0),   VValues.get(1, 3, 0),     1e-4);
	check("Hermite_3v-115 Evaluation trasp  ", 	XAVValues.get(4, 5, 0),   VValues.get(2, 3, 0),     1e-4);
	check("Hermite_3v-116 Evaluation trasp  ", 	XAVValues.get(5, 5, 0),   VValues.get(3, 3, 0),     1e-4);
	check("Hermite_3v-117 Evaluation trasp  ", 	XAVValues.get(0, 6, 0),   VValues.get(0, 4, 0),     1e-4);
	check("Hermite_3v-118 Evaluation trasp  ", 	XAVValues.get(2, 6, 0),   VValues.get(1, 4, 0),     1e-4);
	check("Hermite_3v-119 Evaluation trasp  ", 	XAVValues.get(4, 6, 0),   VValues.get(2, 4, 0),     1e-4);
	check("Hermite_3v-120 Evaluation trasp  ", 	XAVValues.get(5, 6, 0),   VValues.get(3, 4, 0),     1e-4);
	check("Hermite_3v-121 Evaluation trasp  ", 	XAVValues.get(0, 8, 0),   VValues.get(0, 5, 0),     1e-4);
	check("Hermite_3v-122 Evaluation trasp  ", 	XAVValues.get(2, 8, 0),   VValues.get(1, 5, 0),     1e-4);
	check("Hermite_3v-123 Evaluation trasp  ", 	XAVValues.get(4, 8, 0),   VValues.get(2, 5, 0),     1e-4);
	check("Hermite_3v-124 Evaluation trasp  ", 	XAVValues.get(5, 8, 0),   VValues.get(3, 5, 0),     1e-4);

	check("Hermite_3v-131 Evaluation trasp  ", 	XAVValues.get(0, 0, 2),   VValues.get(0, 0, 1),     1e-4);
	check("Hermite_3v-132 Evaluation trasp  ", 	XAVValues.get(2, 0, 2),   VValues.get(1, 0, 1),     1e-4);
	check("Hermite_3v-133 Evaluation trasp  ", 	XAVValues.get(4, 0, 2),   VValues.get(2, 0, 1),     1e-4);
	check("Hermite_3v-134 Evaluation trasp  ", 	XAVValues.get(5, 0, 2),   VValues.get(3, 0, 1),     1e-4);
	check("Hermite_3v-135 Evaluation trasp  ", 	XAVValues.get(0, 2, 2),   VValues.get(0, 1, 1),     1e-4);
	check("Hermite_3v-136 Evaluation trasp  ", 	XAVValues.get(2, 2, 2),   VValues.get(1, 1, 1),     1e-4);
	check("Hermite_3v-137 Evaluation trasp  ", 	XAVValues.get(4, 2, 2),   VValues.get(2, 1, 1),     1e-4);
	check("Hermite_3v-138 Evaluation trasp  ", 	XAVValues.get(5, 2, 2),   VValues.get(3, 1, 1),     1e-4);
	check("Hermite_3v-139 Evaluation trasp  ", 	XAVValues.get(0, 3, 2),   VValues.get(0, 2, 1),     1e-4);
	check("Hermite_3v-140 Evaluation trasp  ", 	XAVValues.get(2, 3, 2),   VValues.get(1, 2, 1),     1e-4);
	check("Hermite_3v-141 Evaluation trasp  ", 	XAVValues.get(4, 3, 2),   VValues.get(2, 2, 1),     1e-4);
	check("Hermite_3v-142 Evaluation trasp  ", 	XAVValues.get(5, 3, 2),   VValues.get(3, 2, 1),     1e-4);
	check("Hermite_3v-143 Evaluation trasp  ", 	XAVValues.get(0, 5, 2),   VValues.get(0, 3, 1),     1e-4);
	check("Hermite_3v-144 Evaluation trasp  ", 	XAVValues.get(2, 5, 2),   VValues.get(1, 3, 1),     1e-4);
	check("Hermite_3v-145 Evaluation trasp  ", 	XAVValues.get(4, 5, 2),   VValues.get(2, 3, 1),     1e-4);
	check("Hermite_3v-146 Evaluation trasp  ", 	XAVValues.get(5, 5, 2),   VValues.get(3, 3, 1),     1e-4);
	check("Hermite_3v-147 Evaluation trasp  ", 	XAVValues.get(0, 6, 2),   VValues.get(0, 4, 1),     1e-4);
	check("Hermite_3v-148 Evaluation trasp  ", 	XAVValues.get(2, 6, 2),   VValues.get(1, 4, 1),     1e-4);
	check("Hermite_3v-149 Evaluation trasp  ", 	XAVValues.get(4, 6, 2),   VValues.get(2, 4, 1),     1e-4);
	check("Hermite_3v-150 Evaluation trasp  ", 	XAVValues.get(5, 6, 2),   VValues.get(3, 4, 1),     1e-4);
	check("Hermite_3v-151 Evaluation trasp  ", 	XAVValues.get(0, 8, 2),   VValues.get(0, 5, 1),     1e-4);
	check("Hermite_3v-152 Evaluation trasp  ", 	XAVValues.get(2, 8, 2),   VValues.get(1, 5, 1),     1e-4);
	check("Hermite_3v-153 Evaluation trasp  ", 	XAVValues.get(4, 8, 2),   VValues.get(2, 5, 1),     1e-4);
	check("Hermite_3v-154 Evaluation trasp  ", 	XAVValues.get(5, 8, 2),   VValues.get(3, 5, 1),     1e-4);

	check("Hermite_3v-161 Evaluation trasp  ", 	XAVValues.get(0, 0, 4),   VValues.get(0, 0, 2),     1e-4);
	check("Hermite_3v-162 Evaluation trasp  ", 	XAVValues.get(2, 0, 4),   VValues.get(1, 0, 2),     1e-4);
	check("Hermite_3v-163 Evaluation trasp  ", 	XAVValues.get(4, 0, 4),   VValues.get(2, 0, 2),     1e-4);
	check("Hermite_3v-164 Evaluation trasp  ", 	XAVValues.get(5, 0, 4),   VValues.get(3, 0, 2),     1e-4);
	check("Hermite_3v-165 Evaluation trasp  ", 	XAVValues.get(0, 2, 4),   VValues.get(0, 1, 2),     1e-4);
	check("Hermite_3v-166 Evaluation trasp  ", 	XAVValues.get(2, 2, 4),   VValues.get(1, 1, 2),     1e-4);
	check("Hermite_3v-167 Evaluation trasp  ", 	XAVValues.get(4, 2, 4),   VValues.get(2, 1, 2),     1e-4);
	check("Hermite_3v-168 Evaluation trasp  ", 	XAVValues.get(5, 2, 4),   VValues.get(3, 1, 2),     1e-4);
	check("Hermite_3v-169 Evaluation trasp  ", 	XAVValues.get(0, 3, 4),   VValues.get(0, 2, 2),     1e-4);
	check("Hermite_3v-170 Evaluation trasp  ", 	XAVValues.get(2, 3, 4),   VValues.get(1, 2, 2),     1e-4);
	check("Hermite_3v-171 Evaluation trasp  ", 	XAVValues.get(4, 3, 4),   VValues.get(2, 2, 2),     1e-4);
	check("Hermite_3v-172 Evaluation trasp  ", 	XAVValues.get(5, 3, 4),   VValues.get(3, 2, 2),     1e-4);
	check("Hermite_3v-173 Evaluation trasp  ", 	XAVValues.get(0, 5, 4),   VValues.get(0, 3, 2),     1e-4);
	check("Hermite_3v-174 Evaluation trasp  ", 	XAVValues.get(2, 5, 4),   VValues.get(1, 3, 2),     1e-4);
	check("Hermite_3v-175 Evaluation trasp  ", 	XAVValues.get(4, 5, 4),   VValues.get(2, 3, 2),     1e-4);
	check("Hermite_3v-176 Evaluation trasp  ", 	XAVValues.get(5, 5, 4),   VValues.get(3, 3, 2),     1e-4);
	check("Hermite_3v-177 Evaluation trasp  ", 	XAVValues.get(0, 6, 4),   VValues.get(0, 4, 2),     1e-4);
	check("Hermite_3v-178 Evaluation trasp  ", 	XAVValues.get(2, 6, 4),   VValues.get(1, 4, 2),     1e-4);
	check("Hermite_3v-179 Evaluation trasp  ", 	XAVValues.get(4, 6, 4),   VValues.get(2, 4, 2),     1e-4);
	check("Hermite_3v-180 Evaluation trasp  ", 	XAVValues.get(5, 6, 4),   VValues.get(3, 4, 2),     1e-4);
	check("Hermite_3v-181 Evaluation trasp  ", 	XAVValues.get(0, 8, 4),   VValues.get(0, 5, 2),     1e-4);
	check("Hermite_3v-182 Evaluation trasp  ", 	XAVValues.get(2, 8, 4),   VValues.get(1, 5, 2),     1e-4);
	check("Hermite_3v-183 Evaluation trasp  ", 	XAVValues.get(4, 8, 4),   VValues.get(2, 5, 2),     1e-4);
	check("Hermite_3v-184 Evaluation trasp  ", 	XAVValues.get(5, 8, 4),   VValues.get(3, 5, 2),     1e-4);

	// Compare 3D table evaluation with its simmetrical both at the nodes and at other values
	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
				check("Hermite_3v-300 Evaluation symm   ", 	XVValues.get(k, j, i),   XAVValues.get(k, j, i),     1e-4);		
			}
		}
	}

	math::vec3 XDDiffs(6, 9, 5);
	double dum1(1.), dum2(1.), dum3(1.);

	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
                XDDiffs.get(k, j, i) = Otab.d_dt(XpointsE.get(i), XpointsF.get(j), XpointsG.get(k), dum1, dum2, dum3);
			}
		}
	}

	// Check 3D table differentials at the nodes
	check("Hermite_3v-301 Differential      ", 	XDDiffs.get(0, 0, 0),   Otab.get_slopes_d1().get(0, 0, 0) + Otab.get_slopes_d2().get(0, 0, 0) + Otab.get_slopes_d3().get(0, 0, 0),    1e-1);
	check("Hermite_3v-302 Differential      ", 	XDDiffs.get(2, 0, 0),   Otab.get_slopes_d1().get(1, 0, 0) + Otab.get_slopes_d2().get(1, 0, 0) + Otab.get_slopes_d3().get(1, 0, 0),    1e-1);
	check("Hermite_3v-303 Differential      ", 	XDDiffs.get(4, 0, 0),   Otab.get_slopes_d1().get(2, 0, 0) + Otab.get_slopes_d2().get(2, 0, 0) + Otab.get_slopes_d3().get(2, 0, 0),    1e-1);
	check("Hermite_3v-304 Differential      ", 	XDDiffs.get(5, 0, 0),   Otab.get_slopes_d1().get(3, 0, 0) + Otab.get_slopes_d2().get(3, 0, 0) + Otab.get_slopes_d3().get(3, 0, 0),    1e-1);
	check("Hermite_3v-305 Differential      ", 	XDDiffs.get(0, 2, 0),   Otab.get_slopes_d1().get(0, 1, 0) + Otab.get_slopes_d2().get(0, 1, 0) + Otab.get_slopes_d3().get(0, 1, 0),    1e-1);
	check("Hermite_3v-306 Differential      ", 	XDDiffs.get(2, 2, 0),   Otab.get_slopes_d1().get(1, 1, 0) + Otab.get_slopes_d2().get(1, 1, 0) + Otab.get_slopes_d3().get(1, 1, 0),    1e-1);
	check("Hermite_3v-307 Differential      ", 	XDDiffs.get(4, 2, 0),   Otab.get_slopes_d1().get(2, 1, 0) + Otab.get_slopes_d2().get(2, 1, 0) + Otab.get_slopes_d3().get(2, 1, 0),    1e-1);
	check("Hermite_3v-308 Differential      ", 	XDDiffs.get(5, 2, 0),   Otab.get_slopes_d1().get(3, 1, 0) + Otab.get_slopes_d2().get(3, 1, 0) + Otab.get_slopes_d3().get(3, 1, 0),    1e-1);
	check("Hermite_3v-309 Differential      ", 	XDDiffs.get(0, 3, 0),   Otab.get_slopes_d1().get(0, 2, 0) + Otab.get_slopes_d2().get(0, 2, 0) + Otab.get_slopes_d3().get(0, 2, 0),    1e-1);
	check("Hermite_3v-310 Differential      ", 	XDDiffs.get(2, 3, 0),   Otab.get_slopes_d1().get(1, 2, 0) + Otab.get_slopes_d2().get(1, 2, 0) + Otab.get_slopes_d3().get(1, 2, 0),    1e-1);
	check("Hermite_3v-311 Differential      ", 	XDDiffs.get(4, 3, 0),   Otab.get_slopes_d1().get(2, 2, 0) + Otab.get_slopes_d2().get(2, 2, 0) + Otab.get_slopes_d3().get(2, 2, 0),    1e-1);
	check("Hermite_3v-312 Differential      ", 	XDDiffs.get(5, 3, 0),   Otab.get_slopes_d1().get(3, 2, 0) + Otab.get_slopes_d2().get(3, 2, 0) + Otab.get_slopes_d3().get(3, 2, 0),    1e-1);
	check("Hermite_3v-313 Differential      ", 	XDDiffs.get(0, 5, 0),   Otab.get_slopes_d1().get(0, 3, 0) + Otab.get_slopes_d2().get(0, 3, 0) + Otab.get_slopes_d3().get(0, 3, 0),    1e-1);
	check("Hermite_3v-314 Differential      ", 	XDDiffs.get(2, 5, 0),   Otab.get_slopes_d1().get(1, 3, 0) + Otab.get_slopes_d2().get(1, 3, 0) + Otab.get_slopes_d3().get(1, 3, 0),    1e-1);
	check("Hermite_3v-315 Differential      ", 	XDDiffs.get(4, 5, 0),   Otab.get_slopes_d1().get(2, 3, 0) + Otab.get_slopes_d2().get(2, 3, 0) + Otab.get_slopes_d3().get(2, 3, 0),	  1e-1);
	check("Hermite_3v-316 Differential      ", 	XDDiffs.get(5, 5, 0),   Otab.get_slopes_d1().get(3, 3, 0) + Otab.get_slopes_d2().get(3, 3, 0) + Otab.get_slopes_d3().get(3, 3, 0),    1e-1);
	check("Hermite_3v-317 Differential      ", 	XDDiffs.get(0, 6, 0),   Otab.get_slopes_d1().get(0, 4, 0) + Otab.get_slopes_d2().get(0, 4, 0) + Otab.get_slopes_d3().get(0, 4, 0),    1e-1);
	check("Hermite_3v-318 Differential      ", 	XDDiffs.get(2, 6, 0),   Otab.get_slopes_d1().get(1, 4, 0) + Otab.get_slopes_d2().get(1, 4, 0) + Otab.get_slopes_d3().get(1, 4, 0),    1e-1);
	check("Hermite_3v-319 Differential      ", 	XDDiffs.get(4, 6, 0),   Otab.get_slopes_d1().get(2, 4, 0) + Otab.get_slopes_d2().get(2, 4, 0) + Otab.get_slopes_d3().get(2, 4, 0),    1e-1);
	check("Hermite_3v-320 Differential      ", 	XDDiffs.get(5, 6, 0),   Otab.get_slopes_d1().get(3, 4, 0) + Otab.get_slopes_d2().get(3, 4, 0) + Otab.get_slopes_d3().get(3, 4, 0),    1e-1);
	check("Hermite_3v-321 Differential      ", 	XDDiffs.get(0, 8, 0),   Otab.get_slopes_d1().get(0, 5, 0) + Otab.get_slopes_d2().get(0, 5, 0) + Otab.get_slopes_d3().get(0, 5, 0),    1e-1);
	check("Hermite_3v-322 Differential      ", 	XDDiffs.get(2, 8, 0),   Otab.get_slopes_d1().get(1, 5, 0) + Otab.get_slopes_d2().get(1, 5, 0) + Otab.get_slopes_d3().get(1, 5, 0),    1e-1);
	check("Hermite_3v-323 Differential      ", 	XDDiffs.get(4, 8, 0),   Otab.get_slopes_d1().get(2, 5, 0) + Otab.get_slopes_d2().get(2, 5, 0) + Otab.get_slopes_d3().get(2, 5, 0),    1e-1);
	check("Hermite_3v-324 Differential      ", 	XDDiffs.get(5, 8, 0),   Otab.get_slopes_d1().get(3, 5, 0) + Otab.get_slopes_d2().get(3, 5, 0) + Otab.get_slopes_d3().get(3, 5, 0),    1e-1);

	check("Hermite_3v-331 Differential      ", 	XDDiffs.get(0, 0, 2),   Otab.get_slopes_d1().get(0, 0, 1) + Otab.get_slopes_d2().get(0, 0, 1) + Otab.get_slopes_d3().get(0, 0, 1),    1e-1);
	check("Hermite_3v-332 Differential      ", 	XDDiffs.get(2, 0, 2),   Otab.get_slopes_d1().get(1, 0, 1) + Otab.get_slopes_d2().get(1, 0, 1) + Otab.get_slopes_d3().get(1, 0, 1),    1e-1);
	check("Hermite_3v-333 Differential      ", 	XDDiffs.get(4, 0, 2),   Otab.get_slopes_d1().get(2, 0, 1) + Otab.get_slopes_d2().get(2, 0, 1) + Otab.get_slopes_d3().get(2, 0, 1),    1e-1);
	check("Hermite_3v-334 Differential      ", 	XDDiffs.get(5, 0, 2),   Otab.get_slopes_d1().get(3, 0, 1) + Otab.get_slopes_d2().get(3, 0, 1) + Otab.get_slopes_d3().get(3, 0, 1),    1e-1);
	check("Hermite_3v-335 Differential      ", 	XDDiffs.get(0, 2, 2),   Otab.get_slopes_d1().get(0, 1, 1) + Otab.get_slopes_d2().get(0, 1, 1) + Otab.get_slopes_d3().get(0, 1, 1),    1e-1);
	check("Hermite_3v-336 Differential      ", 	XDDiffs.get(2, 2, 2),   Otab.get_slopes_d1().get(1, 1, 1) + Otab.get_slopes_d2().get(1, 1, 1) + Otab.get_slopes_d3().get(1, 1, 1),    1e-1);
	check("Hermite_3v-337 Differential      ", 	XDDiffs.get(4, 2, 2),   Otab.get_slopes_d1().get(2, 1, 1) + Otab.get_slopes_d2().get(2, 1, 1) + Otab.get_slopes_d3().get(2, 1, 1),    1e-1);
	check("Hermite_3v-338 Differential      ", 	XDDiffs.get(5, 2, 2),   Otab.get_slopes_d1().get(3, 1, 1) + Otab.get_slopes_d2().get(3, 1, 1) + Otab.get_slopes_d3().get(3, 1, 1),    1e-1);
	check("Hermite_3v-339 Differential      ", 	XDDiffs.get(0, 3, 2),   Otab.get_slopes_d1().get(0, 2, 1) + Otab.get_slopes_d2().get(0, 2, 1) + Otab.get_slopes_d3().get(0, 2, 1),    1e-1);
	check("Hermite_3v-340 Differential      ", 	XDDiffs.get(2, 3, 2),   Otab.get_slopes_d1().get(1, 2, 1) + Otab.get_slopes_d2().get(1, 2, 1) + Otab.get_slopes_d3().get(1, 2, 1),    1e-1);
	check("Hermite_3v-341 Differential      ", 	XDDiffs.get(4, 3, 2),   Otab.get_slopes_d1().get(2, 2, 1) + Otab.get_slopes_d2().get(2, 2, 1) + Otab.get_slopes_d3().get(2, 2, 1),    1e-1);
	check("Hermite_3v-342 Differential      ", 	XDDiffs.get(5, 3, 2),   Otab.get_slopes_d1().get(3, 2, 1) + Otab.get_slopes_d2().get(3, 2, 1) + Otab.get_slopes_d3().get(3, 2, 1),    1e-1);
	check("Hermite_3v-343 Differential      ", 	XDDiffs.get(0, 5, 2),   Otab.get_slopes_d1().get(0, 3, 1) + Otab.get_slopes_d2().get(0, 3, 1) + Otab.get_slopes_d3().get(0, 3, 1),    1e-1);
	check("Hermite_3v-344 Differential      ", 	XDDiffs.get(2, 5, 2),   Otab.get_slopes_d1().get(1, 3, 1) + Otab.get_slopes_d2().get(1, 3, 1) + Otab.get_slopes_d3().get(1, 3, 1),    1e-1);
	check("Hermite_3v-345 Differential      ", 	XDDiffs.get(4, 5, 2),   Otab.get_slopes_d1().get(2, 3, 1) + Otab.get_slopes_d2().get(2, 3, 1) + Otab.get_slopes_d3().get(2, 3, 1),	  1e-1);
	check("Hermite_3v-346 Differential      ", 	XDDiffs.get(5, 5, 2),   Otab.get_slopes_d1().get(3, 3, 1) + Otab.get_slopes_d2().get(3, 3, 1) + Otab.get_slopes_d3().get(3, 3, 1),    1e-1);
	check("Hermite_3v-347 Differential      ", 	XDDiffs.get(0, 6, 2),   Otab.get_slopes_d1().get(0, 4, 1) + Otab.get_slopes_d2().get(0, 4, 1) + Otab.get_slopes_d3().get(0, 4, 1),    1e-1);
	check("Hermite_3v-348 Differential      ", 	XDDiffs.get(2, 6, 2),   Otab.get_slopes_d1().get(1, 4, 1) + Otab.get_slopes_d2().get(1, 4, 1) + Otab.get_slopes_d3().get(1, 4, 1),    1e-1);
	check("Hermite_3v-349 Differential      ", 	XDDiffs.get(4, 6, 2),   Otab.get_slopes_d1().get(2, 4, 1) + Otab.get_slopes_d2().get(2, 4, 1) + Otab.get_slopes_d3().get(2, 4, 1),    1e-1);
	check("Hermite_3v-350 Differential      ", 	XDDiffs.get(5, 6, 2),   Otab.get_slopes_d1().get(3, 4, 1) + Otab.get_slopes_d2().get(3, 4, 1) + Otab.get_slopes_d3().get(3, 4, 1),    1e-1);
	check("Hermite_3v-351 Differential      ", 	XDDiffs.get(0, 8, 2),   Otab.get_slopes_d1().get(0, 5, 1) + Otab.get_slopes_d2().get(0, 5, 1) + Otab.get_slopes_d3().get(0, 5, 1),    1e-1);
	check("Hermite_3v-352 Differential      ", 	XDDiffs.get(2, 8, 2),   Otab.get_slopes_d1().get(1, 5, 1) + Otab.get_slopes_d2().get(1, 5, 1) + Otab.get_slopes_d3().get(1, 5, 1),    1e-1);
	check("Hermite_3v-353 Differential      ", 	XDDiffs.get(4, 8, 2),   Otab.get_slopes_d1().get(2, 5, 1) + Otab.get_slopes_d2().get(2, 5, 1) + Otab.get_slopes_d3().get(2, 5, 1),    1e-1);
	check("Hermite_3v-354 Differential      ", 	XDDiffs.get(5, 8, 2),   Otab.get_slopes_d1().get(3, 5, 1) + Otab.get_slopes_d2().get(3, 5, 1) + Otab.get_slopes_d3().get(3, 5, 1),    1e-1);

	check("Hermite_3v-361 Differential      ", 	XDDiffs.get(0, 0, 4),   Otab.get_slopes_d1().get(0, 0, 2) + Otab.get_slopes_d2().get(0, 0, 2) + Otab.get_slopes_d3().get(0, 0, 2),    1e-1);
	check("Hermite_3v-362 Differential      ", 	XDDiffs.get(2, 0, 4),   Otab.get_slopes_d1().get(1, 0, 2) + Otab.get_slopes_d2().get(1, 0, 2) + Otab.get_slopes_d3().get(1, 0, 2),    1e-1);
	check("Hermite_3v-363 Differential      ", 	XDDiffs.get(4, 0, 4),   Otab.get_slopes_d1().get(2, 0, 2) + Otab.get_slopes_d2().get(2, 0, 2) + Otab.get_slopes_d3().get(2, 0, 2),    1e-1);
	check("Hermite_3v-364 Differential      ", 	XDDiffs.get(5, 0, 4),   Otab.get_slopes_d1().get(3, 0, 2) + Otab.get_slopes_d2().get(3, 0, 2) + Otab.get_slopes_d3().get(3, 0, 2),    1e-1);
	check("Hermite_3v-365 Differential      ", 	XDDiffs.get(0, 2, 4),   Otab.get_slopes_d1().get(0, 1, 2) + Otab.get_slopes_d2().get(0, 1, 2) + Otab.get_slopes_d3().get(0, 1, 2),    1e-1);
	check("Hermite_3v-366 Differential      ", 	XDDiffs.get(2, 2, 4),   Otab.get_slopes_d1().get(1, 1, 2) + Otab.get_slopes_d2().get(1, 1, 2) + Otab.get_slopes_d3().get(1, 1, 2),    1e-1);
	check("Hermite_3v-367 Differential      ", 	XDDiffs.get(4, 2, 4),   Otab.get_slopes_d1().get(2, 1, 2) + Otab.get_slopes_d2().get(2, 1, 2) + Otab.get_slopes_d3().get(2, 1, 2),    1e-1);
	check("Hermite_3v-368 Differential      ", 	XDDiffs.get(5, 2, 4),   Otab.get_slopes_d1().get(3, 1, 2) + Otab.get_slopes_d2().get(3, 1, 2) + Otab.get_slopes_d3().get(3, 1, 2),    1e-1);
	check("Hermite_3v-369 Differential      ", 	XDDiffs.get(0, 3, 4),   Otab.get_slopes_d1().get(0, 2, 2) + Otab.get_slopes_d2().get(0, 2, 2) + Otab.get_slopes_d3().get(0, 2, 2),    1e-1);
	check("Hermite_3v-370 Differential      ", 	XDDiffs.get(2, 3, 4),   Otab.get_slopes_d1().get(1, 2, 2) + Otab.get_slopes_d2().get(1, 2, 2) + Otab.get_slopes_d3().get(1, 2, 2),    1e-1);
	check("Hermite_3v-371 Differential      ", 	XDDiffs.get(4, 3, 4),   Otab.get_slopes_d1().get(2, 2, 2) + Otab.get_slopes_d2().get(2, 2, 2) + Otab.get_slopes_d3().get(2, 2, 2),    1e-1);
	check("Hermite_3v-372 Differential      ", 	XDDiffs.get(5, 3, 4),   Otab.get_slopes_d1().get(3, 2, 2) + Otab.get_slopes_d2().get(3, 2, 2) + Otab.get_slopes_d3().get(3, 2, 2),    1e-1);
	check("Hermite_3v-373 Differential      ", 	XDDiffs.get(0, 5, 4),   Otab.get_slopes_d1().get(0, 3, 2) + Otab.get_slopes_d2().get(0, 3, 2) + Otab.get_slopes_d3().get(0, 3, 2),    1e-1);
	check("Hermite_3v-374 Differential      ", 	XDDiffs.get(2, 5, 4),   Otab.get_slopes_d1().get(1, 3, 2) + Otab.get_slopes_d2().get(1, 3, 2) + Otab.get_slopes_d3().get(1, 3, 2),    1e-1);
	check("Hermite_3v-375 Differential      ", 	XDDiffs.get(4, 5, 4),   Otab.get_slopes_d1().get(2, 3, 2) + Otab.get_slopes_d2().get(2, 3, 2) + Otab.get_slopes_d3().get(2, 3, 2),	  1e-1);
	check("Hermite_3v-376 Differential      ", 	XDDiffs.get(5, 5, 4),   Otab.get_slopes_d1().get(3, 3, 2) + Otab.get_slopes_d2().get(3, 3, 2) + Otab.get_slopes_d3().get(3, 3, 2),    1e-1);
	check("Hermite_3v-377 Differential      ", 	XDDiffs.get(0, 6, 4),   Otab.get_slopes_d1().get(0, 4, 2) + Otab.get_slopes_d2().get(0, 4, 2) + Otab.get_slopes_d3().get(0, 4, 2),    1e-1);
	check("Hermite_3v-378 Differential      ", 	XDDiffs.get(2, 6, 4),   Otab.get_slopes_d1().get(1, 4, 2) + Otab.get_slopes_d2().get(1, 4, 2) + Otab.get_slopes_d3().get(1, 4, 2),    1e-1);
	check("Hermite_3v-379 Differential      ", 	XDDiffs.get(4, 6, 4),   Otab.get_slopes_d1().get(2, 4, 2) + Otab.get_slopes_d2().get(2, 4, 2) + Otab.get_slopes_d3().get(2, 4, 2),    1e-1);
	check("Hermite_3v-380 Differential      ", 	XDDiffs.get(5, 6, 4),   Otab.get_slopes_d1().get(3, 4, 2) + Otab.get_slopes_d2().get(3, 4, 2) + Otab.get_slopes_d3().get(3, 4, 2),    1e-1);
	check("Hermite_3v-381 Differential      ", 	XDDiffs.get(0, 8, 4),   Otab.get_slopes_d1().get(0, 5, 2) + Otab.get_slopes_d2().get(0, 5, 2) + Otab.get_slopes_d3().get(0, 5, 2),    1e-1);
	check("Hermite_3v-382 Differential      ", 	XDDiffs.get(2, 8, 4),   Otab.get_slopes_d1().get(1, 5, 2) + Otab.get_slopes_d2().get(1, 5, 2) + Otab.get_slopes_d3().get(1, 5, 2),    1e-1);
	check("Hermite_3v-383 Differential      ", 	XDDiffs.get(4, 8, 4),   Otab.get_slopes_d1().get(2, 5, 2) + Otab.get_slopes_d2().get(2, 5, 2) + Otab.get_slopes_d3().get(2, 5, 2),    1e-1);
	check("Hermite_3v-384 Differential      ", 	XDDiffs.get(5, 8, 4),   Otab.get_slopes_d1().get(3, 5, 2) + Otab.get_slopes_d2().get(3, 5, 2) + Otab.get_slopes_d3().get(3, 5, 2),    1e-1);

	math::vec3 XADDiffs(6, 9, 5);

	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
                XADDiffs.get(k, j, i) = OAtab.d_dt(XpointsG.get(k), XpointsF.get(j), XpointsE.get(i), dum1, dum2, dum3);
			}
		}
	}

	// Check 3D simmetrical (or traspose) table differentials at the nodes
	check("Hermite_3v-401 Differential trasp", 	XADDiffs.get(0, 0, 0),   Otab.get_slopes_d1().get(0, 0, 0) + Otab.get_slopes_d2().get(0, 0, 0) + Otab.get_slopes_d3().get(0, 0, 0),    1e-1);
	check("Hermite_3v-402 Differential trasp", 	XADDiffs.get(2, 0, 0),   Otab.get_slopes_d1().get(1, 0, 0) + Otab.get_slopes_d2().get(1, 0, 0) + Otab.get_slopes_d3().get(1, 0, 0),    1e-1);
	check("Hermite_3v-403 Differential trasp", 	XADDiffs.get(4, 0, 0),   Otab.get_slopes_d1().get(2, 0, 0) + Otab.get_slopes_d2().get(2, 0, 0) + Otab.get_slopes_d3().get(2, 0, 0),    1e-1);
	check("Hermite_3v-404 Differential trasp", 	XADDiffs.get(5, 0, 0),   Otab.get_slopes_d1().get(3, 0, 0) + Otab.get_slopes_d2().get(3, 0, 0) + Otab.get_slopes_d3().get(3, 0, 0),    1e-1);
	check("Hermite_3v-405 Differential trasp", 	XADDiffs.get(0, 2, 0),   Otab.get_slopes_d1().get(0, 1, 0) + Otab.get_slopes_d2().get(0, 1, 0) + Otab.get_slopes_d3().get(0, 1, 0),    1e-1);
	check("Hermite_3v-406 Differential trasp", 	XADDiffs.get(2, 2, 0),   Otab.get_slopes_d1().get(1, 1, 0) + Otab.get_slopes_d2().get(1, 1, 0) + Otab.get_slopes_d3().get(1, 1, 0),    1e-1);
	check("Hermite_3v-407 Differential trasp", 	XADDiffs.get(4, 2, 0),   Otab.get_slopes_d1().get(2, 1, 0) + Otab.get_slopes_d2().get(2, 1, 0) + Otab.get_slopes_d3().get(2, 1, 0),    1e-1);
	check("Hermite_3v-408 Differential trasp", 	XADDiffs.get(5, 2, 0),   Otab.get_slopes_d1().get(3, 1, 0) + Otab.get_slopes_d2().get(3, 1, 0) + Otab.get_slopes_d3().get(3, 1, 0),    1e-1);
	check("Hermite_3v-409 Differential trasp", 	XADDiffs.get(0, 3, 0),   Otab.get_slopes_d1().get(0, 2, 0) + Otab.get_slopes_d2().get(0, 2, 0) + Otab.get_slopes_d3().get(0, 2, 0),    1e-1);
	check("Hermite_3v-410 Differential trasp", 	XADDiffs.get(2, 3, 0),   Otab.get_slopes_d1().get(1, 2, 0) + Otab.get_slopes_d2().get(1, 2, 0) + Otab.get_slopes_d3().get(1, 2, 0),    1e-1);
	check("Hermite_3v-411 Differential trasp", 	XADDiffs.get(4, 3, 0),   Otab.get_slopes_d1().get(2, 2, 0) + Otab.get_slopes_d2().get(2, 2, 0) + Otab.get_slopes_d3().get(2, 2, 0),    1e-1);
	check("Hermite_3v-412 Differential trasp", 	XADDiffs.get(5, 3, 0),   Otab.get_slopes_d1().get(3, 2, 0) + Otab.get_slopes_d2().get(3, 2, 0) + Otab.get_slopes_d3().get(3, 2, 0),    1e-1);
	check("Hermite_3v-413 Differential trasp", 	XADDiffs.get(0, 5, 0),   Otab.get_slopes_d1().get(0, 3, 0) + Otab.get_slopes_d2().get(0, 3, 0) + Otab.get_slopes_d3().get(0, 3, 0),    1e-1);
	check("Hermite_3v-414 Differential trasp", 	XADDiffs.get(2, 5, 0),   Otab.get_slopes_d1().get(1, 3, 0) + Otab.get_slopes_d2().get(1, 3, 0) + Otab.get_slopes_d3().get(1, 3, 0),    1e-1);
	check("Hermite_3v-415 Differential trasp", 	XADDiffs.get(4, 5, 0),   Otab.get_slopes_d1().get(2, 3, 0) + Otab.get_slopes_d2().get(2, 3, 0) + Otab.get_slopes_d3().get(2, 3, 0),	  1e-1);
	check("Hermite_3v-416 Differential trasp", 	XADDiffs.get(5, 5, 0),   Otab.get_slopes_d1().get(3, 3, 0) + Otab.get_slopes_d2().get(3, 3, 0) + Otab.get_slopes_d3().get(3, 3, 0),    1e-1);
	check("Hermite_3v-417 Differential trasp", 	XADDiffs.get(0, 6, 0),   Otab.get_slopes_d1().get(0, 4, 0) + Otab.get_slopes_d2().get(0, 4, 0) + Otab.get_slopes_d3().get(0, 4, 0),    1e-1);
	check("Hermite_3v-418 Differential trasp", 	XADDiffs.get(2, 6, 0),   Otab.get_slopes_d1().get(1, 4, 0) + Otab.get_slopes_d2().get(1, 4, 0) + Otab.get_slopes_d3().get(1, 4, 0),    1e-1);
	check("Hermite_3v-419 Differential trasp", 	XADDiffs.get(4, 6, 0),   Otab.get_slopes_d1().get(2, 4, 0) + Otab.get_slopes_d2().get(2, 4, 0) + Otab.get_slopes_d3().get(2, 4, 0),    1e-1);
	check("Hermite_3v-420 Differential trasp", 	XADDiffs.get(5, 6, 0),   Otab.get_slopes_d1().get(3, 4, 0) + Otab.get_slopes_d2().get(3, 4, 0) + Otab.get_slopes_d3().get(3, 4, 0),    1e-1);
	check("Hermite_3v-421 Differential trasp", 	XADDiffs.get(0, 8, 0),   Otab.get_slopes_d1().get(0, 5, 0) + Otab.get_slopes_d2().get(0, 5, 0) + Otab.get_slopes_d3().get(0, 5, 0),    1e-1);
	check("Hermite_3v-422 Differential trasp", 	XADDiffs.get(2, 8, 0),   Otab.get_slopes_d1().get(1, 5, 0) + Otab.get_slopes_d2().get(1, 5, 0) + Otab.get_slopes_d3().get(1, 5, 0),    1e-1);
	check("Hermite_3v-423 Differential trasp", 	XADDiffs.get(4, 8, 0),   Otab.get_slopes_d1().get(2, 5, 0) + Otab.get_slopes_d2().get(2, 5, 0) + Otab.get_slopes_d3().get(2, 5, 0),    1e-1);
	check("Hermite_3v-424 Differential trasp", 	XADDiffs.get(5, 8, 0),   Otab.get_slopes_d1().get(3, 5, 0) + Otab.get_slopes_d2().get(3, 5, 0) + Otab.get_slopes_d3().get(3, 5, 0),    1e-1);

	check("Hermite_3v-431 Differential trasp", 	XADDiffs.get(0, 0, 2),   Otab.get_slopes_d1().get(0, 0, 1) + Otab.get_slopes_d2().get(0, 0, 1) + Otab.get_slopes_d3().get(0, 0, 1),    1e-1);
	check("Hermite_3v-432 Differential trasp", 	XADDiffs.get(2, 0, 2),   Otab.get_slopes_d1().get(1, 0, 1) + Otab.get_slopes_d2().get(1, 0, 1) + Otab.get_slopes_d3().get(1, 0, 1),    1e-1);
	check("Hermite_3v-433 Differential trasp", 	XADDiffs.get(4, 0, 2),   Otab.get_slopes_d1().get(2, 0, 1) + Otab.get_slopes_d2().get(2, 0, 1) + Otab.get_slopes_d3().get(2, 0, 1),    1e-1);
	check("Hermite_3v-434 Differential trasp", 	XADDiffs.get(5, 0, 2),   Otab.get_slopes_d1().get(3, 0, 1) + Otab.get_slopes_d2().get(3, 0, 1) + Otab.get_slopes_d3().get(3, 0, 1),    1e-1);
	check("Hermite_3v-435 Differential trasp", 	XADDiffs.get(0, 2, 2),   Otab.get_slopes_d1().get(0, 1, 1) + Otab.get_slopes_d2().get(0, 1, 1) + Otab.get_slopes_d3().get(0, 1, 1),    1e-1);
	check("Hermite_3v-436 Differential trasp", 	XADDiffs.get(2, 2, 2),   Otab.get_slopes_d1().get(1, 1, 1) + Otab.get_slopes_d2().get(1, 1, 1) + Otab.get_slopes_d3().get(1, 1, 1),    1e-1);
	check("Hermite_3v-437 Differential trasp", 	XADDiffs.get(4, 2, 2),   Otab.get_slopes_d1().get(2, 1, 1) + Otab.get_slopes_d2().get(2, 1, 1) + Otab.get_slopes_d3().get(2, 1, 1),    1e-1);
	check("Hermite_3v-438 Differential trasp", 	XADDiffs.get(5, 2, 2),   Otab.get_slopes_d1().get(3, 1, 1) + Otab.get_slopes_d2().get(3, 1, 1) + Otab.get_slopes_d3().get(3, 1, 1),    1e-1);
	check("Hermite_3v-439 Differential trasp", 	XADDiffs.get(0, 3, 2),   Otab.get_slopes_d1().get(0, 2, 1) + Otab.get_slopes_d2().get(0, 2, 1) + Otab.get_slopes_d3().get(0, 2, 1),    1e-1);
	check("Hermite_3v-440 Differential trasp", 	XADDiffs.get(2, 3, 2),   Otab.get_slopes_d1().get(1, 2, 1) + Otab.get_slopes_d2().get(1, 2, 1) + Otab.get_slopes_d3().get(1, 2, 1),    1e-1);
	check("Hermite_3v-441 Differential trasp", 	XADDiffs.get(4, 3, 2),   Otab.get_slopes_d1().get(2, 2, 1) + Otab.get_slopes_d2().get(2, 2, 1) + Otab.get_slopes_d3().get(2, 2, 1),    1e-1);
	check("Hermite_3v-442 Differential trasp", 	XADDiffs.get(5, 3, 2),   Otab.get_slopes_d1().get(3, 2, 1) + Otab.get_slopes_d2().get(3, 2, 1) + Otab.get_slopes_d3().get(3, 2, 1),    1e-1);
	check("Hermite_3v-443 Differential trasp", 	XADDiffs.get(0, 5, 2),   Otab.get_slopes_d1().get(0, 3, 1) + Otab.get_slopes_d2().get(0, 3, 1) + Otab.get_slopes_d3().get(0, 3, 1),    1e-1);
	check("Hermite_3v-444 Differential trasp", 	XADDiffs.get(2, 5, 2),   Otab.get_slopes_d1().get(1, 3, 1) + Otab.get_slopes_d2().get(1, 3, 1) + Otab.get_slopes_d3().get(1, 3, 1),    1e-1);
	check("Hermite_3v-445 Differential trasp", 	XADDiffs.get(4, 5, 2),   Otab.get_slopes_d1().get(2, 3, 1) + Otab.get_slopes_d2().get(2, 3, 1) + Otab.get_slopes_d3().get(2, 3, 1),	  1e-1);
	check("Hermite_3v-446 Differential trasp", 	XADDiffs.get(5, 5, 2),   Otab.get_slopes_d1().get(3, 3, 1) + Otab.get_slopes_d2().get(3, 3, 1) + Otab.get_slopes_d3().get(3, 3, 1),    1e-1);
	check("Hermite_3v-447 Differential trasp", 	XADDiffs.get(0, 6, 2),   Otab.get_slopes_d1().get(0, 4, 1) + Otab.get_slopes_d2().get(0, 4, 1) + Otab.get_slopes_d3().get(0, 4, 1),    1e-1);
	check("Hermite_3v-448 Differential trasp", 	XADDiffs.get(2, 6, 2),   Otab.get_slopes_d1().get(1, 4, 1) + Otab.get_slopes_d2().get(1, 4, 1) + Otab.get_slopes_d3().get(1, 4, 1),    1e-1);
	check("Hermite_3v-449 Differential trasp", 	XADDiffs.get(4, 6, 2),   Otab.get_slopes_d1().get(2, 4, 1) + Otab.get_slopes_d2().get(2, 4, 1) + Otab.get_slopes_d3().get(2, 4, 1),    1e-1);
	check("Hermite_3v-450 Differential trasp", 	XADDiffs.get(5, 6, 2),   Otab.get_slopes_d1().get(3, 4, 1) + Otab.get_slopes_d2().get(3, 4, 1) + Otab.get_slopes_d3().get(3, 4, 1),    1e-1);
	check("Hermite_3v-451 Differential trasp", 	XADDiffs.get(0, 8, 2),   Otab.get_slopes_d1().get(0, 5, 1) + Otab.get_slopes_d2().get(0, 5, 1) + Otab.get_slopes_d3().get(0, 5, 1),    1e-1);
	check("Hermite_3v-452 Differential trasp", 	XADDiffs.get(2, 8, 2),   Otab.get_slopes_d1().get(1, 5, 1) + Otab.get_slopes_d2().get(1, 5, 1) + Otab.get_slopes_d3().get(1, 5, 1),    1e-1);
	check("Hermite_3v-453 Differential trasp", 	XADDiffs.get(4, 8, 2),   Otab.get_slopes_d1().get(2, 5, 1) + Otab.get_slopes_d2().get(2, 5, 1) + Otab.get_slopes_d3().get(2, 5, 1),    1e-1);
	check("Hermite_3v-454 Differential trasp", 	XADDiffs.get(5, 8, 2),   Otab.get_slopes_d1().get(3, 5, 1) + Otab.get_slopes_d2().get(3, 5, 1) + Otab.get_slopes_d3().get(3, 5, 1),    1e-1);

	check("Hermite_3v-461 Differential trasp", 	XADDiffs.get(0, 0, 4),   Otab.get_slopes_d1().get(0, 0, 2) + Otab.get_slopes_d2().get(0, 0, 2) + Otab.get_slopes_d3().get(0, 0, 2),    1e-1);
	check("Hermite_3v-462 Differential trasp", 	XADDiffs.get(2, 0, 4),   Otab.get_slopes_d1().get(1, 0, 2) + Otab.get_slopes_d2().get(1, 0, 2) + Otab.get_slopes_d3().get(1, 0, 2),    1e-1);
	check("Hermite_3v-463 Differential trasp", 	XADDiffs.get(4, 0, 4),   Otab.get_slopes_d1().get(2, 0, 2) + Otab.get_slopes_d2().get(2, 0, 2) + Otab.get_slopes_d3().get(2, 0, 2),    1e-1);
	check("Hermite_3v-464 Differential trasp", 	XADDiffs.get(5, 0, 4),   Otab.get_slopes_d1().get(3, 0, 2) + Otab.get_slopes_d2().get(3, 0, 2) + Otab.get_slopes_d3().get(3, 0, 2),    1e-1);
	check("Hermite_3v-465 Differential trasp", 	XADDiffs.get(0, 2, 4),   Otab.get_slopes_d1().get(0, 1, 2) + Otab.get_slopes_d2().get(0, 1, 2) + Otab.get_slopes_d3().get(0, 1, 2),    1e-1);
	check("Hermite_3v-466 Differential trasp", 	XADDiffs.get(2, 2, 4),   Otab.get_slopes_d1().get(1, 1, 2) + Otab.get_slopes_d2().get(1, 1, 2) + Otab.get_slopes_d3().get(1, 1, 2),    1e-1);
	check("Hermite_3v-467 Differential trasp", 	XADDiffs.get(4, 2, 4),   Otab.get_slopes_d1().get(2, 1, 2) + Otab.get_slopes_d2().get(2, 1, 2) + Otab.get_slopes_d3().get(2, 1, 2),    1e-1);
	check("Hermite_3v-468 Differential trasp", 	XADDiffs.get(5, 2, 4),   Otab.get_slopes_d1().get(3, 1, 2) + Otab.get_slopes_d2().get(3, 1, 2) + Otab.get_slopes_d3().get(3, 1, 2),    1e-1);
	check("Hermite_3v-469 Differential trasp", 	XADDiffs.get(0, 3, 4),   Otab.get_slopes_d1().get(0, 2, 2) + Otab.get_slopes_d2().get(0, 2, 2) + Otab.get_slopes_d3().get(0, 2, 2),    1e-1);
	check("Hermite_3v-470 Differential trasp", 	XADDiffs.get(2, 3, 4),   Otab.get_slopes_d1().get(1, 2, 2) + Otab.get_slopes_d2().get(1, 2, 2) + Otab.get_slopes_d3().get(1, 2, 2),    1e-1);
	check("Hermite_3v-471 Differential trasp", 	XADDiffs.get(4, 3, 4),   Otab.get_slopes_d1().get(2, 2, 2) + Otab.get_slopes_d2().get(2, 2, 2) + Otab.get_slopes_d3().get(2, 2, 2),    1e-1);
	check("Hermite_3v-472 Differential trasp", 	XADDiffs.get(5, 3, 4),   Otab.get_slopes_d1().get(3, 2, 2) + Otab.get_slopes_d2().get(3, 2, 2) + Otab.get_slopes_d3().get(3, 2, 2),    1e-1);
	check("Hermite_3v-473 Differential trasp", 	XADDiffs.get(0, 5, 4),   Otab.get_slopes_d1().get(0, 3, 2) + Otab.get_slopes_d2().get(0, 3, 2) + Otab.get_slopes_d3().get(0, 3, 2),    1e-1);
	check("Hermite_3v-474 Differential trasp", 	XADDiffs.get(2, 5, 4),   Otab.get_slopes_d1().get(1, 3, 2) + Otab.get_slopes_d2().get(1, 3, 2) + Otab.get_slopes_d3().get(1, 3, 2),    1e-1);
	check("Hermite_3v-475 Differential trasp", 	XADDiffs.get(4, 5, 4),   Otab.get_slopes_d1().get(2, 3, 2) + Otab.get_slopes_d2().get(2, 3, 2) + Otab.get_slopes_d3().get(2, 3, 2),	   1e-1);
	check("Hermite_3v-476 Differential trasp", 	XADDiffs.get(5, 5, 4),   Otab.get_slopes_d1().get(3, 3, 2) + Otab.get_slopes_d2().get(3, 3, 2) + Otab.get_slopes_d3().get(3, 3, 2),    1e-1);
	check("Hermite_3v-477 Differential trasp", 	XADDiffs.get(0, 6, 4),   Otab.get_slopes_d1().get(0, 4, 2) + Otab.get_slopes_d2().get(0, 4, 2) + Otab.get_slopes_d3().get(0, 4, 2),    1e-1);
	check("Hermite_3v-478 Differential trasp", 	XADDiffs.get(2, 6, 4),   Otab.get_slopes_d1().get(1, 4, 2) + Otab.get_slopes_d2().get(1, 4, 2) + Otab.get_slopes_d3().get(1, 4, 2),    1e-1);
	check("Hermite_3v-479 Differential trasp", 	XADDiffs.get(4, 6, 4),   Otab.get_slopes_d1().get(2, 4, 2) + Otab.get_slopes_d2().get(2, 4, 2) + Otab.get_slopes_d3().get(2, 4, 2),    1e-1);
	check("Hermite_3v-480 Differential trasp", 	XADDiffs.get(5, 6, 4),   Otab.get_slopes_d1().get(3, 4, 2) + Otab.get_slopes_d2().get(3, 4, 2) + Otab.get_slopes_d3().get(3, 4, 2),    1e-1);
	check("Hermite_3v-481 Differential trasp", 	XADDiffs.get(0, 8, 4),   Otab.get_slopes_d1().get(0, 5, 2) + Otab.get_slopes_d2().get(0, 5, 2) + Otab.get_slopes_d3().get(0, 5, 2),    1e-1);
	check("Hermite_3v-482 Differential trasp", 	XADDiffs.get(2, 8, 4),   Otab.get_slopes_d1().get(1, 5, 2) + Otab.get_slopes_d2().get(1, 5, 2) + Otab.get_slopes_d3().get(1, 5, 2),    1e-1);
	check("Hermite_3v-483 Differential trasp", 	XADDiffs.get(4, 8, 4),   Otab.get_slopes_d1().get(2, 5, 2) + Otab.get_slopes_d2().get(2, 5, 2) + Otab.get_slopes_d3().get(2, 5, 2),    1e-1);
	check("Hermite_3v-484 Differential trasp", 	XADDiffs.get(5, 8, 4),   Otab.get_slopes_d1().get(3, 5, 2) + Otab.get_slopes_d2().get(3, 5, 2) + Otab.get_slopes_d3().get(3, 5, 2),    1e-1);
	
	// Compare 3D table differentials with its simmetrical both at the nodes and at other values
	for (int i = 0; i != XpointsE.size1(); ++i) {
		for (int j = 0; j != XpointsF.size1(); ++j) {
			for (int k = 0; k != XpointsG.size1(); ++k) {
				check("Hermite_3v-500 Differential symm ", 	XDDiffs.get(k, j, i),   XADDiffs.get(k, j, i),     1e-4);		
			}
		}
	}
	
} // closes test_hermite_d

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////









