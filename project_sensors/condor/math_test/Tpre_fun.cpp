 #include "Tpre_fun.h"

math::test::Tpre_fun::Tpre_fun(jail::counter& Ocounter)
: ::jail::unit_test(Ocounter) {
}
/* constructor based on counter */

void math::test::Tpre_fun::run() {
	::jail::unit_test::run();

    test0v_f_null();
    test0v_f_constant();
    test1v_f_lineal();
    test1v_f_parabolic();
    test1v_f_cubic();
    test1v_f_steps();
    test1v_f_table1V();
    test1v_f_table1Veq();
    test1v_f_table1V_spl();
    test2v_f_lineal_double();
    test2v_f_table2V();
    test2v_f_table2Veq();
    test2v_f_tabular2V();
    test3v_f_lineal_triple();
    test3v_f_table3V();
    test3v_f_table3Veq();
    test3v_f_tabular3V();
    test4v_f_table4Veq();
    test4v_f_table4V();

	finished();
}
/* execute tests and write results on console */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test0v_f_null() {
	math::f_null P1;
	double val = P1.value();
	check("0v f_null value ", val, 0., 1e-12);
	double der = P1.d_dt();
	check("0v f_null d_dt  ", der, 0., 1e-12);
} // closes test0v_f_null

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test0v_f_constant() {
	double f0 = 7.32;
	f_constant P1(f0);
    double val = P1.value();
    check("0v f_constant value ", val, 7.32, 1e-12);
	double der = P1.d_dt();
    check("0v f_constant d_dt  ", der, 0., 1e-12);
} // closes test0v_f_constant

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_lineal() {
    double x1 = 30.0;
    double dx1_dt = 5.0;
	double f0 = 3.27;
	double f1 = -8.41;
	f_lineal P1(f0, f1);
	double val = P1.value(x1);
	double val1 = f0 + f1 * x1;
    check("1v f_lineal value ", val, val1, 1e-12);
    double der = P1.d_dt(x1, dx1_dt);
	check("1v f_lineal d_dt  ", der, f1 * dx1_dt, 1e-12);
} // closes test1v_f_lineal

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_parabolic() {
    double x1 = 120.0;
    double dx1_dt = 5.0;
	double f0 = 14.97;
	double f1 = -18.41;
	double f2 = 2.21;
	f_parabolic P1(f0, f1, f2);
	double val = P1.value(x1);
    double val1 = f0 + f1 * x1 + f2 * x1 * x1;
    check("1v f_parabolic value  ", val, val1, 1e-11);
    double der = P1.d_dt(x1, dx1_dt);
    double der1 = f1 * dx1_dt + 2 * f2 * x1 * dx1_dt;
	check("1v f_parabolic d_dt   ", der, der1, 1e-12);
} // closes test1v_f_parabolic

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_cubic() {
    double x1 = 120.0;
    double dx1_dt = 5.0;
	double f0 = 14.97;
	double f1 = -18.41;
	double f2 = 2.21;
	double f3 = 0.02;
	f_cubic P1(f0, f1, f2, f3);
    double val = P1.value(x1);
	double val1 = f0 + f1 * x1 + f2 * pow(x1, 2) + f3 * pow(x1, 3);
	check("1v f_cubic value     ", val, val1, 1e-11);
    double der = P1.d_dt(x1, dx1_dt);
	double der1 = f1 * dx1_dt + 2 * f2 * x1 * dx1_dt + 3 * f3 * x1 * x1 * dx1_dt;
	check("1v f_cubic d_dt      ", der, der1, 1e-10);
} // closes test1v_f_cubic

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_steps() {
    double dx1_dt = 5.0;
	double p0	= -2.6;
	double p1	= 7.8;
	double p2	= 14.7;
	double v0	= 3.4;
	double v1	= 23.4;
	vec1 points(3);
	points.set(0, p0);
	points.set(1, p1);
	points.set(2, p2);
	vec1 values(2);
	values.set(0, v0);
	values.set(1, v1);
	f_steps P1(points, values);

    double x1 = -2.6;
	double val1 = P1.value(x1);
    x1 = 0.;
    double val2 = P1.value(x1);
    x1 = 7.3;
    double val3 = P1.value(x1);
    x1 = 7.8;
    double val4 = P1.value(x1);
    x1 = 14.08;
    double val5 = P1.value(x1);
    x1 = -12.6;
    double val7 = P1.value(x1);
    x1 = 28.5;
    double val8 = P1.value(x1);
    x1 = 14.7;
    double val6 = P1.value(x1);

	check("1v f_steps value        ", val1, v0,	1e-12);
	check("1v f_steps value        ", val2, v0,	1e-12);
	check("1v f_steps value        ", val3, v0,	1e-12);
	check("1v f_steps value        ", val4, v1,	1e-11);
	check("1v f_steps value        ", val5, v1,	1e-11);
	check("1v f_steps value        ", val6, v1,	1e-11);
	check("1v f_steps value        ", val7, v0,	1e-12);
	check("1v f_steps value        ", val8, v1,	1e-11);

    double der = P1.d_dt(x1, dx1_dt);
	check("1v f_steps d_dt         ", der,	0,	1e-16);
	// ========================================================================
} // closes test1v_f_steps

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_table1V() {
    double dx1_dt = 0.5;
    double p0	= -2.;
	double p1	= 8.;
	double p2	= 14.;
	double p3	= 22.;
	vec1 points(4);
	points.set(0, p0);
	points.set(1, p1);
	points.set(2, p2);
	points.set(3, p3);
	double v0	= 340.;
	double v1	= 20.;
	double v2	= 220.;
	double v3   = 200.;
	vec1 values(4);
	values.set(0, v0);
	values.set(1, v1);
	values.set(2, v2);
	values.set(3, v3);

	f_table1V P1(points, values, math::logic::lagrange_first_precompute);
	f_table1V P2(points, values, math::logic::lagrange_first);
	f_table1V P3(points, values, math::logic::lagrange_second);
	f_table1V P4(points, values, math::logic::lagrange_third);
	f_table1V P5(points, values, math::logic::spline);
	f_table1V P6(points, values, math::logic::hermite_first);
	f_table1V P7(points, values, math::logic::hermite_second);

	// ===== ===== (p0) ===== =====
	double x1 = -2.;
	double val1 = P1.value(x1);
	double der1 = P1.d_dt(x1, dx1_dt);
	double res1 = 340.; //v0
	check("1v f_table1V precom value      ", val1, res1,	1e-12);
	check("1v f_table1V precom d_dt       ", der1, 0.5 * (-320.0 / 10.0),	1e-12);
    val1 = P2.value(x1);
	der1 = P2.d_dt(x1, dx1_dt);
	res1 = 340.;
	check("1v f_table1V linear value      ", val1, res1,	1e-12);
	check("1v f_table1V linear d_dt       ", der1,	0.5 * (-320.0 / 10.0), 1e-11);
    val1 = P3.value(x1);
	res1 = 340.;
	check("1v f_table1V quadra value      ", val1, res1,	1e-12);
    val1 = P4.value(x1);
	res1 = 340.;
	check("1v f_table1V cubic valueu      ", val1, res1,	1e-12);
    val1 = P5.value(x1);
	res1 = 340.;
	check("1v f_table1V spline value      ", val1, res1,	1e-12);
    val1 = P6.value(x1);
	res1 = 340.;
	check("1v f_table1V hermite value     ", val1, res1,	1e-12);

	// ===== ===== 20% between p0 and p1 ===== ======
    x1 = 0.;
    double val2 = P1.value(x1);
    double der2 = P1.d_dt(x1, dx1_dt);
	double res2 = 276.;
	check("1v f_table1V precom value      ", val2, res2,	1e-12);
	check("1v f_table1V precom d_dt       ", der2,	0.5 * (-320.0 / 10.0),			1e-12);
    val2 = P2.value(x1);
    der2 = P2.d_dt(x1, dx1_dt);
	res2 = 276.;
	check("1v f_table1V linear value      ", val2, res2,	1e-12);
	check("1v f_table1V linear d_dt       ", der2,	0.5 * (-320.0 / 10.0),			1e-11);


	// ===== ===== 60% between p0 and p1 ===== =====
	x1 = 4.;
    double val3 = P1.value(x1);
	double der3 = P1.d_dt(x1, dx1_dt);
	double res3 = 148.;
	check("1v f_table1V precom value      ", val3, res3,	1e-12);
	check("1v f_table1V precom d_dt       ", der3,	0.5 * (-320.0 / 10.0),			1e-12);
    val3 = P2.value(x1);
	der3 = P2.d_dt(x1, dx1_dt);
	res3 = 148.;
	check("1v f_table1V value     ", val3, res3,	1e-12);
	check("1v f_table1V d_dt          ", der3,	0.5 * (-320.0 / 10.0),			1e-11);

	 // ===== ===== p1 ===== =====
    x1 = 8.;
    double val4 = P1.value(x1);
	double der4 = P1.d_dt(x1, dx1_dt);
	double res4 = 20.;
	check("1v f_table1V precom value      ", val4, res4,	1e-12);
	check("1v f_table1V precom d_dt       ", der4,	0.5 * (200.0 / 6.0),			1e-12);
    val4 = P2.value(x1);
    der4 = P2.d_dt(x1, dx1_dt);
	res4 = 20.;
	check("1v f_table1V linear value      ", val4, res4,	1e-12);
	check("1v f_table1V linear d_dt       ", der4,	0.5 * (200.0 / 6.0),			1e-12);
    val4 = P3.value(x1);
	res4 = 20.;
	check("1v f_table1V quadra value      ", val4, res4,	1e-12);
    val4 = P4.value(x1);
	res4 = 20.;
	check("1v f_table1V cubic  value      ", val4, res4,	1e-12);
    val4 = P5.value(x1);
	res4 = 20.;
	check("1v f_table1V spline value      ", val4, res4,	1e-12);
    val4 = P6.value(x1);
	res4 = 20.;
	check("1v f_table1V hermite value     ", val4, res4,	1e-12);
    val4 = P7.value(x1);
	res4 = 20.;
	check("1v f_table1V hermite value     ", val4, res4,	1e-12);

	// ===== ===== 50% between p1 and p2 ===== =====
    x1 = 11.;
    double val5 = P1.value(x1);
    double der5 = P1.d_dt(x1, dx1_dt);
	double res5 = 120.;
	check("1v f_table1V precom value      ", val5, res5,	1e-12);
	check("1v f_table1V precom d_dt       ", der5,	0.5 * (200.0 / 6.0),			1e-12);
    val5 = P2.value(x1);
	der5 = P2.d_dt(x1, dx1_dt);
	res5 = 120.;
	check("1v f_table1V linear value      ", val5, res5,	1e-12);
	check("1v f_table1V linear d_dt       ", der5,	0.5 * (200.0 / 6.0),			1e-11);

	// ===== ===== p3 ===== =====
    x1 = 22.;
    double val6 = P1.value(x1);
	double der6 = P1.d_dt(x1, dx1_dt);
	double res6 = 200.;
	check("1v f_table1V precom value      ", val6, res6,	1e-12);
	check("1v f_table1V precom d_dt       ", der6,	0.5 * (-20.0 / 8.0),			1e-12);
    val6 = P2.value(x1);
    der6 = P2.d_dt(x1, dx1_dt);
	res6 = 200.;
	check("1v f_table1V linear value      ", val6, res6,	1e-12);
	check("1v f_table1V linear d_dt       ", der6,	0.5 * (-20.0 / 8.0),			1e-11);
    val6 = P3.value(x1);
	res6 = 200.;
	check("1v f_table1V quadra value      ", val6, res6,	1e-12);
    val6 = P4.value(x1);
	res6 = 200.;
	check("1v f_table1V cubic  value      ", val6, res6,	1e-12);
    val6 = P5.value(x1);
	res6 = 200.;
	check("1v f_table1V spline value      ", val6, res6,	1e-12);
    val6 = P6.value(x1);
	res6 = 200.;
	check("1v f_table1V hermite value     ", val6, res6,	1e-12);

	// ===== ===== less than p0 ===== =====
    x1 = -3.0;
    double val7 = P1.value(x1);
	double der7 = P1.d_dt(x1, dx1_dt);
	double res7 = 372;
	check("1v f_table1V precom value      ", val7, res7,	1e-12);
	check("1v f_table1V precom d_dt       ", der7,	0.5 * (-320.0 / 10.0),			1e-12);
    val7 = P2.value(x1);
    der7 = P2.d_dt(x1, dx1_dt);
	res7 = 372;
	check("1v f_table1V linear value      ", val7, res7,	1e-12);
	check("1v f_table1V linear d_dt       ", der7,	0.5 * (-320.0 / 10.0),			1e-10);

	// ===== ===== more than p3 ===== =====
    x1 = 24.0;
    double val8 = P1.value(x1);
    double der8 = P1.d_dt(x1, dx1_dt);
	double res8 = 200. + 0.25 * (-20.0);
	check("1v f_table1V precom value      ", val8, res8,	1e-12);
	check("1v f_table1V precom d_dt       ", der8,	0.5 * (-20.0 / 8.0),			1e-12);
    val8 = P2.value(x1);
    der8 = P2.d_dt(x1, dx1_dt);
	res8 = 200. + 0.25 * (-20.0);
	check("1v f_table1V linear value      ", val8, res8,	1e-12);
	check("1v f_table1V linear d_dt       ", der8,	0.5 * (-20.0 / 8.0),			1e-10);
} // closes test1v_f_table1V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_table1Veq() {
    double dx1_dt = 5.0;

	double p0	= 200.0;
	double p1	= 300.0;
	double p2	= 400.0;
	double v0	= 10.;
	double v1	= 20.;
	double v2	= 40.;
	vec1 points(3);
	points.set(0, p0);
	points.set(1, p1);
	points.set(2, p2);
	vec1 values(3);
	values.set(0, v0);
	values.set(1, v1);
	values.set(2, v2);
	f_table1V P1(points, values, math::logic::lagrange_first_precompute);
	check("1v-gcx f_table1Veq is equispaced        ", P1.get_equi1(), true);
	f_table1V P2(points, values, math::logic::lagrange_first);
	check("1v-gcx f_table1Veq is equispaced        ", P2.get_equi1(), true);
	f_table1V P6(points, values, math::logic::hermite_first);
	check("1v-gcx f_table1Veq is equispaced        ", P6.get_equi1(), true);
	f_table1V P7(points, values, math::logic::hermite_second);
	check("1v-gcx f_table1Veq is equispaced        ", P7.get_equi1(), true);

	double xder1 =  10.0 / 100.;
	double xder2 =  20.0 / 100.;

	// ===== Obtain value
	double x1 = 200.0; // (p0)
	double val1 = P1.value(x1);
    double der1 = P1.d_dt(x1, dx1_dt);
	double res1 = 10.;
	check("1v f_table1Veq value    ", val1, res1,	1e-12);
	check("1v f_table1Veq d_dt     ", der1,	5 * xder1,	1e-12);
    val1 = P2.value(x1);
    der1 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val1, res1,	1e-12);
    check("1v f_table1Veq d_dt     ", der1,	5 * xder1,	1e-10);

	x1 = 225.0; // 25% between p0 and p1
    double val2 = P1.value(x1);
    double der2 = P1.d_dt(x1, dx1_dt);
	double res2 = 12.5;
	check("1v f_table1Veq value    ", val2, res2,	1e-12);
	check("1v f_table1Veq d_dt     ", der2,	5 * xder1,	1e-12);
    val2 = P2.value(x1);
    der2 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val2, res2,	1e-12);
	check("1v f_table1Veq d_dt     ", der2,	5 * xder1,	1e-10);

	x1 = 275.0; // 75% between p0 and p1
    double val3 = P1.value(x1);
    double der3 = P1.d_dt(x1, dx1_dt);
	double res3 = 17.5;
	check("1v f_table1Veq value    ", val3, res3,	1e-12);
	check("1v f_table1Veq d_dt     ", der3,	5 * xder1,	1e-12);
    val3 = P2.value(x1);
    der3 = P2.d_dt(x1, dx1_dt);
	check("1v f_table1Veq value    ", val3, res3,	1e-12);
    check("1v f_table1Veq d_dt     ", der3,	5 * xder1,	1e-10);

	x1 = 300.0; // p1
    double val4 = P1.value(x1);
    double der4 = P1.d_dt(x1, dx1_dt);
	double res4 = 20.;
	check("1v f_table1Veq value    ", val4, res4,	1e-12);
	check("1v f_table1Veq d_dt     ", der4,	5 * xder2,	1e-12);
    val4 = P2.value(x1);
    der4 = P2.d_dt(x1, dx1_dt);
	check("1v f_table1Veq value    ", val4, res4,	1e-12);
    check("1v f_table1Veq d_dt     ", der4,	5 * xder2,	1e-10);
    val4 = P6.value(x1);
	check("1v f_table1Veq value    ", val4, res4,	1e-12);
    val4 = P7.value(x1);
	check("1v f_table1Veq value    ", val4, res4,	1e-12);

	x1 = 350.0; // 50% between p1 and p2
    double val5 = P1.value(x1);
    double der5 = P1.d_dt(x1, dx1_dt);
    double res5 = 30;
	check("1v f_table1Veq value    ", val5, res5,	1e-12);
	check("1v f_table1Veq d_dt     ", der5,	5 * xder2,	1e-12);
    val5 = P2.value(x1);
    der5 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val5, res5,	1e-12);
    check("1v f_table1Veq d_dt     ", der5,	5 * xder2,	1e-10);

	x1 = 175.0;
    double val7 = P1.value(x1);
    double der7 = P1.d_dt(x1, dx1_dt);
    double res7 = 7.5;
	//double res7 = -150. - 1.25 * 170.0 + 360.;
	check("1v f_table1Veq value    ", val7, res7,	1e-12);
	check("1v f_table1Veq d_dt     ", der7,	5 * xder1,	1e-12);
    val7 = P2.value(x1);
    der7 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val7, res7,	1e-12);
    check("1v f_table1Veq d_dt     ", der7,	5 * xder1,	1e-10);

	x1 = 450.0; // more than (p2)
    double val8 = P1.value(x1);
    double der8 = P1.d_dt(x1, dx1_dt);
	double res8 = 40. + 0.5 * 20.;
	check("1v f_table1Veq value    ", val8, res8,	1e-12);
	check("1v f_table1Veq d_dt     ", der8,	5 * xder2,	1e-12);
    val8 = P2.value(x1);
    der8 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val8, res8,	1e-12);
    check("1v f_table1Veq d_dt     ", der8,	5 * xder2,	1e-10);

	x1 = 400.0; // (p2)
    double val6 = P1.value(x1);
	double der6 = P1.d_dt(x1, dx1_dt);
	double res6 = 40.;
	check("1v f_table1Veq value    ", val6, res6,	1e-12);
	check("1v f_table1Veq d_dt     ", der6,	5 * xder2,	1e-12);
    val6 = P2.value(x1);
    der6 = P2.d_dt(x1, dx1_dt);
    check("1v f_table1Veq value    ", val6, res6,	1e-12);
    check("1v f_table1Veq d_dt     ", der6,	5 * xder2,	1e-10);
} // closes test1v_f_table1Veq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test1v_f_table1V_spl() {
	// f_table1V_num 340. at -2., 20. in 8., 220. at 14., 250. at 20., 56. at 24.
	// ==========================================================================
    double dx1_dt = 0.5;

	// Create f_table1V object
	double p0	= -2.;
    double p1	= 8.;
    double p2	= 14.;
    double p3 = 20.;
    double p4 = 24.;
	double v0	= 340.;
    double v1	= 20.;
    double v2	= 220.;
    double v3 = 250.;
    double v4 = 56.;
	std::vector<double>* points = new std::vector<double>(5);
	(*points)[0] = p0;
    (*points)[1] = p1;
    (*points)[2] = p2;
    (*points)[3] = p3;
    (*points)[4] = p4;
	std::vector<double>* values = new std::vector<double>(5);
	(*values)[0] = v0;
    (*values)[1] = v1;
    (*values)[2] = v2;
    (*values)[3] = v3;
    (*values)[4] = v4;

	vec1 Points(5);
	Points.set(0, -2.);
	Points.set(1,  8.);
	Points.set(2, 14.);
	Points.set(3, 20.);
	Points.set(4, 24.);
	vec1 Values(5);
	Values.set(0, 340.);
	Values.set(1,  20.);
	Values.set(2, 220.);
	Values.set(3, 250.);
	Values.set(4,  56.);

	f_table1V     P1(Points, Values, math::logic::spline);
	f_table1V_spl P2(points, values);

	// Differential with time holders
    double dd_dt1(0.), dd_dt2(0.);

	// Obtain value
	double x_tas = -2.; // (p0)
	double val11 = P1.value(x_tas);
	double val21 = P2.value(x_tas);
	double der11 = P1.d_dt(x_tas, dx1_dt);
    double der21 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 0.; // 20% between p0 and p1
    double val12 = P1.value(x_tas);
    double val22 = P2.value(x_tas);
    double der12 = P1.d_dt(x_tas, dx1_dt);
    double der22 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 4.; // 60% between p0 and p1
    double val13 = P1.value(x_tas);
    double val23 = P2.value(x_tas);
    double der13 = P1.d_dt(x_tas, dx1_dt);
    double der23 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 8.; // p1
    double val14 = P1.value(x_tas);
    double val24 = P2.value(x_tas);
    double der14 = P1.d_dt(x_tas, dx1_dt);
    double der24 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 11.; // 50% between p1 and p2
    double val15 = P1.value(x_tas);
    double val25 = P2.value(x_tas);
    double der15 = P1.d_dt(x_tas, dx1_dt);
    double der25 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 14.; // (p2)
    double val16 = P1.value(x_tas);
    double val26 = P2.value(x_tas);
    double der16 = P1.d_dt(x_tas, dx1_dt);
    double der26 = P2.d_dt(x_tas, dx1_dt);

	x_tas = 26.; // higher than (p4)
    double val17 = P1.value(x_tas);
    double val27 = P2.value(x_tas);
    double der17 = P1.d_dt(x_tas, dx1_dt);
    double der27 = P2.d_dt(x_tas, dx1_dt);

	x_tas = -6.; // lower than (p0)
    double val18 = P1.value(x_tas);
    double val28 = P2.value(x_tas);
    double der18 = P1.d_dt(x_tas, dx1_dt);
    double der28 = P2.d_dt(x_tas, dx1_dt);

	// The accuracy gets manually degraded as many digits as has the result
	check("1v-220 f_table1V_spl value      ", val11, val21,	1e-13);
	check("1v-221 f_table1V_spl value      ", val12, val22,	1e-13);
	check("1v-222 f_table1V_spl value      ", val13, val23,	1e-13);
	check("1v-223 f_table1V_spl value      ", val14, val24,	1e-13);
	check("1v-224 f_table1V_spl value      ", val15, val25,	1e-13);
	check("1v-225 f_table1V_spl value      ", val16, val26,	1e-13);
	check("1v-226 f_table1V_spl value      ", val17, val27,	1e-13);
	check("1v-227 f_table1V_spl value      ", val18, val28,	1e-13);

	check("1v-230 f_table1V_spl d_dt           ", der11, der21,	1e-10);
	check("1v-231 f_table1V_spl d_dt           ", der12, der22,	1e-10);
	check("1v-232 f_table1V_spl d_dt           ", der13, der23,	1e-10);
	check("1v-233 f_table1V_spl d_dt           ", der14, der24,	1e-10);
	check("1v-234 f_table1V_spl d_dt           ", der15, der25,	1e-10);
	check("1v-235 f_table1V_spl d_dt           ", der16, der26,	1e-10);
	check("1v-236 f_table1V_spl d_dt           ", der17, der27,	1e-10);
	check("1v-237 f_table1V_spl d_dt           ", der18, der28,	1e-10);

} // closes test1v_f_table1V_spl

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test2v_f_lineal_double() {
    double x1 = 30.0;
    double x2 = 2500.0;
    double dx1_dt = 5.0;
    double dx2_dt = 25.0;

	double f0	= 3.3;
	double f1x	= 2.7;
	double f1y	= 0.4;
	double f1xy	= 0.001;

	f_lineal_double P1(f0, f1x, f1y, f1xy);
    double val1 = P1.value(x1, x2);
	double res1	= f0 + f1x * x1 + f1y * x2 + f1xy * x1 * x2;
	check("2v f_lineal_double value ", val1, res1, 1e-11);

    double der1 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double xder1 =  (f1x + f1xy * x2) * dx1_dt + (f1y + f1xy * x1) * dx2_dt;
	check("2v f_lineal_double d_dt  ", der1, xder1,	1e-14);
} // closes test2v_f_lineal_double

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test2v_f_table2V() {

	// f_table2V of Hp[ft] as function of m [lb] and phi [deg]
	// =======================================================

    double dx1_dt = 25.0;
    double dx2_dt = 3.0;

	double m0	= -23.6;
	double m1	= -13.6;
	double m2	= -1.6;
	vec1 points1(3);
	points1.set(0, m0);
	points1.set(1, m1);
	points1.set(2, m2);

	double p0	= 12.1;
	double p1	= 52.1;
	double p2	= 72.1;
	double p3	= 82.1;
	vec1 points2(4);
	points2.set(0, p0);
	points2.set(1, p1);
	points2.set(2, p2);
	points2.set(3, p3);

	double a00	= 10;	double a01	= 12;	double a02	= 15;	double a03	= 19;
	double a10	= 20;	double a11	= 32;	double a12	= 45;	double a13	= 69;
	double a20	= 40;	double a21	= 22;	double a22	= 40;	double a23	= 139;
	vec2 Values(4, 3);
	Values.set(0, 0, a00);
	Values.set(1, 0, a01);
	Values.set(2, 0, a02);
	Values.set(3, 0, a03);
	Values.set(0, 1, a10);
	Values.set(1, 1, a11);
	Values.set(2, 1, a12);
	Values.set(3, 1, a13);
	Values.set(0, 2, a20);
	Values.set(1, 2, a21);
	Values.set(2, 2, a22);
	Values.set(3, 2, a23);

	f_table2V P1(points1, points2, Values, math::logic::lagrange_first_precompute);
	f_table2V P2(points1, points2, Values, math::logic::lagrange_first);
	f_table2V P6(points1, points2, Values, math::logic::hermite_first);
	f_table2V P7(points1, points2, Values, math::logic::hermite_second);

	// Evaluate function at different points
	double x1 = -28.6;		// 50% less than (x0)
	double x2 = 2.1;		// 25% less than (y0)
	double val1 = P1.value(x1, x2);
    double der1 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double temp1		= a00 - 0.25 * (a01 - a00);
	double temp2		= a10 - 0.25 * (a11 - a10);
	double res1	= temp1 - 0.50 * (temp2 - temp1);
	double diff_d1_1a	= (a10 - a00) / (m1 - m0);
	double diff_d1_1b	= (a11 - a01) / (m1 - m0);
	double diff_d1_1	= diff_d1_1a - 0.25 * (diff_d1_1b - diff_d1_1a);
	double diff_d2_1a	= (a01 - a00) / (p1 - p0);
	double diff_d2_1b	= (a11 - a10) / (p1 - p0);
	double diff_d2_1	= diff_d2_1a -0.5 * (diff_d2_1b - diff_d2_1a);
	double xder1	= diff_d1_1 * dx1_dt + diff_d2_1 * dx2_dt;
	check("2v f_table2V value      ", val1, res1,	1e-12);
	check("2v f_table2V d_dt       ", der1,	xder1,	1e-11);
    val1 = P2.value(x1, x2);
    der1 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val1, res1,	1e-12);
    check("2v f_table2V d_dt       ", der1,	xder1,	1e-10);

	x1 = -23.6;		// (x0)
	x2 = 2.1;		// 25% less than (y0)
    double val2 = P1.value(x1, x2);
    double der2 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res2 = a00 - 0.25 * (a01 - a00);
	double diff_d1_2a	= (a10 - a00) / (m1 - m0);
	double diff_d1_2b	= (a11 - a01) / (m1 - m0);
	double diff_d1_2	= diff_d1_2a - 0.25 * (diff_d1_2b - diff_d1_2a);
	double diff_d2_2	= (a01 - a00) / (p1 - p0);
	double xder2	= diff_d1_2 * dx1_dt + diff_d2_2 * dx2_dt;
	check("2v f_table2V value      ", val2, res2,	1e-12);
	check("2v f_table2V d_dt       ", der2,	xder2,	1e-11);
    val2 = P2.value(x1, x2);
    der2 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val2, res2,	1e-12);
    check("2v f_table2V d_dt       ", der2,	xder2,	1e-10);

	x1 = -28.6;		// 50% less than (x0)
	x2 = 12.1;		// (y0)
    double val3 = P1.value(x1, x2);
    double der3 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res3	= a00 - 0.5 * (a10 - a00);
		double diff_d1_3	= (a10 - a00) / (m1 - m0);
	double diff_d2_3a	= (a01 - a00) / (p1 - p0);
	double diff_d2_3b	= (a11 - a10) / (p1 - p0);
	double diff_d2_3	= diff_d2_3a -0.5 * (diff_d2_3b - diff_d2_3a);
	double xder3	= diff_d1_3 * dx1_dt + diff_d2_3 * dx2_dt;
	check("2v f_table2V value      ", val3, res3,	1e-12);
	check("2v f_table2V d_dt       ", der3,	xder3,	1e-11);
    val3 = P2.value(x1, x2);
    der3 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val3, res3,	1e-12);
    check("2v f_table2V d_dt       ", der3,	xder3,	1e-10);

	x1 = -23.6;		// (x0)
	x2 = 12.1;		// (y0)
    double val4 = P1.value(x1, x2);
    double der4 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res4 = a00;
	double diff_d1_4	= (a10 - a00) / (m1 - m0);
	double diff_d2_4	= (a01 - a00) / (p1 - p0);
	double AAAA = (a01 - a00) / (p1 - p0) * dx2_dt;
	double xder4	= diff_d1_4 * dx1_dt + diff_d2_4 * dx2_dt;
	check("2v f_table2V value      ", val4, res4,	1e-12);
	check("2v f_table2V d_dt       ", der4,	xder4,	1e-11);
    val4 = P2.value(x1, x2);
    der4 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val4, res4,	1e-12);
    check("2v f_table2V d_dt       ", der4,	xder4,	1e-10);

	x1 = -7.6;		// 50% between (x1) and (x2)
	x2 = 52.1;	// (y1)
    double val5 = P1.value(x1, x2);
    double der5 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res5 = a11 + 0.5 * (a21 - a11);
	double diff_d1_5	= (a21 - a11) / (m2 - m1);
	double diff_d2_5a	= (a12 - a11) / (p2 - p1);
	double diff_d2_5b	= (a22 - a21) / (p2 - p1);
	double diff_d2_5	=  diff_d2_5a + 0.5 * (diff_d2_5b - diff_d2_5a);
	double xder5	= diff_d1_5 * dx1_dt + diff_d2_5 * dx2_dt;
	check("2v f_table2V value      ", val5, res5,	1e-12);
	check("2v f_table2V d_dt       ", der5,	xder5,	1e-11);
    val5 = P2.value(x1, x2);
    der5 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val5, res5,	1e-12);
    check("2v f_table2V d_dt       ", der5,	xder5,	1e-10);

	x1 = -13.6;		// (x1)
	x2 = 57.1;		// 25% between (y1) and (y2)
    double val6 = P1.value(x1, x2);
    double der6 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res6 = a11 + 0.25 * (a12 - a11);
	double diff_d1_6a	= (a21 - a11) / (m2 - m1);
	double diff_d1_6b	= (a22 - a12) / (m2 - m1);
	double diff_d1_6	= diff_d1_6a + 0.25 * (diff_d1_6b - diff_d1_6a);
	double diff_d2_6	= (a12 - a11) / (p2 - p1);
	double xder6	= diff_d1_6 * dx1_dt + diff_d2_6 * dx2_dt;
	check("2v f_table2V value      ", val6, res6,	1e-12);
	check("2v f_table2V d_dt       ", der6,	xder6,	1e-11);
    val6 = P2.value(x1, x2);
    der6 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val6, res6,	1e-12);
	check("2v f_table2V d_dt       ", der6,	xder6,	1e-10);

	x1 = -7.6;		// 50% between (x1) and (x2)
	x2 = 57.1;	// 25% between (y1) and (y2)
    double val7 = P1.value(x1, x2);
    double der7 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	temp1		= a11 + 0.25 * (a12 - a11);
	temp2		= a21 + 0.25 * (a22 - a21);
	double res7	= temp1 + 0.5 * (temp2 - temp1);
	double diff_d1_7a	= (a21 - a11) / (m2 - m1);
	double diff_d1_7b	= (a22 - a12) / (m2 - m1);
	double diff_d1_7	= diff_d1_7a + 0.25 * (diff_d1_7b - diff_d1_7a);
	double diff_d2_7a	= (a12 - a11) / (p2 - p1);
	double diff_d2_7b	= (a22 - a21) / (p2 - p1);
	double diff_d2_7	=  diff_d2_7a + 0.5 * (diff_d2_7b - diff_d2_7a);
	double xder7	= diff_d1_7 * dx1_dt + diff_d2_7 * dx2_dt;
	check("2v f_table2V value      ", val7, res7,	1e-12);
	check("2v f_table2V d_dt       ", der7,	xder7,	1e-10);
    val7 = P2.value(x1, x2);
    der7 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val7, res7,	1e-12);
    check("2v f_table2V d_dt       ", der7,	xder7,	1e-10);

	x1 = -7.6;		// 50% between (x1) and (x2)
	x2 = 82.1;	// (y3)
    double val8 = P1.value(x1, x2);
    double der8 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res8 = a13 + 0.5 * (a23 - a13);
	double diff_d1_8	= (a23 - a13) / (m2 - m1);
	double diff_d2_8a	= (a13 - a12) / (p3 - p2);
	double diff_d2_8b	= (a23 - a22) / (p3 - p2);
	double diff_d2_8	=  diff_d2_8a + 0.5 * (diff_d2_8b - diff_d2_8a);
	double xder8	= diff_d1_8 * dx1_dt + diff_d2_8 * dx2_dt;
	check("2v f_table2V value      ", val8, res8,	1e-12);
	check("2v f_table2V d_dt       ", der8,	xder8,	1e-10);
    val8 = P2.value(x1, x2);
    der8 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val8, res8,	1e-12);
    check("2v f_table2V d_dt       ", der8,	xder8,	1e-9);

	x1 = -1.6;		// (x2)
	x2 = 57.1;	// 25% between (y1) and (y2)
    double val9 = P1.value(x1, x2);
    double der9 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res9 = a21 + 0.25 * (a22 - a21);
	double diff_d1_9a	= (a21 - a11) / (m2 - m1);
	double diff_d1_9b	= (a22 - a12) / (m2 - m1);
	double diff_d1_9	= diff_d1_9a + 0.25 * (diff_d1_9b - diff_d1_9a);
	double diff_d2_9	= (a22 - a21) / (p2 - p1);
	double xder9	= diff_d1_9 * dx1_dt + diff_d2_9 * dx2_dt;
	check("2v f_table2V value      ", val9, res9,	1e-12);
	check("2v f_table2V d_dt       ", der9,	xder9,	1e-10);
    val9 = P2.value(x1, x2);
    der9 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val9, res9,	1e-12);
    check("2v f_table2V d_dt       ", der9,	xder9,	1e-10);

	x1 = -1.6;		// (x2)
	x2 = 82.1;	// (y3)
    double val10 = P1.value(x1, x2);
    double der10 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res10 = a23;
	double diff_d1_10	= (a23 - a13) / (m2 - m1);
	double diff_d2_10	= (a23 - a22) / (p3 - p2);
	double xder10	= diff_d1_10 * dx1_dt +	diff_d2_10 * dx2_dt;
	check("2v f_table2V value      ", val10, res10,	1e-12);
	check("2v f_table2V d_dt       ", der10, xder10,	1e-10);
    val10 = P2.value(x1, x2);
    der10 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val10, res10,	1e-12);
	check("2v f_table2V d_dt       ", der10, xder10,	1e-9);

    val10 = P6.value(x1, x2);
	check("2v f_table2V value      ", val10, res10,	1e-12);

    val10 = P7.value(x1, x2);
	check("2v f_table2V value      ", val10, res10,	1e-12);

	x1 = 1.4;		// 25% higher than (x2)
	x2 = 82.1;	// (y3)
    double val11 = P1.value(x1, x2);
    double der11 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res11 = a13 + 1.25 * (a23 - a13);
	double diff_d1_11	= (a23 - a13) / (m2 - m1);
	double diff_d2_11a	= (a13 - a12) / (p3 - p2);
	double diff_d2_11b	= (a23 - a22) / (p3 - p2);
	double diff_d2_11	=  diff_d2_11a + 1.25 * (diff_d2_11b - diff_d2_11a);
	double xder11	= diff_d1_11 * dx1_dt +	diff_d2_11 * dx2_dt;
	check("2v f_table2V value      ", val11, res11,	1e-12);
	check("2v f_table2V d_dt       ", der11, xder11,	1e-10);
    val11 = P2.value(x1, x2);
    der11 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val11, res11,	1e-12);
    check("2v f_table2V d_dt       ", der11, xder11,	1e-8);

	x1 = -1.6;		// (x2)
	x2 = 89.6;	// 75% higher than (y3)
    double val12 = P1.value(x1, x2);
    double der12 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res12 = a22 + 1.75 * (a23 - a22);
	double diff_d1_12a	= (a22 - a12) / (m2 - m1);
	double diff_d1_12b	= (a23 - a13) / (m2 - m1);
	double diff_d1_12	= diff_d1_12a + 1.75 * (diff_d1_12b - diff_d1_12a);
	double diff_d2_12	= (a23 - a22) / (p3 - p2);
	double xder12	= diff_d1_12 * dx1_dt + diff_d2_12 * dx2_dt;
	check("2v f_table2V value      ", val12, res12,	1e-12);
	check("2v f_table2V d_dt       ", der12, xder12,	1e-11);
    val12 = P2.value(x1, x2);
    der12 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val12, res12,	1e-12);
    check("2v f_table2V d_dt       ", der12, xder12,	1e-9);

	x1 = -28.6;	// 50% less than (x0)
	x2 = 82.1;	// (y3)
    double val13 = P1.value(x1, x2);
    double der13 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res13 = a03 - 0.5 * (a13 - a03);
	double diff_d1_13	= (a13 - a03) / (m1 - m0);
	double diff_d2_13a	= (a03 - a02) / (p3 - p2);
	double diff_d2_13b	= (a13 - a12) / (p3 - p2);
	double diff_d2_13	=  diff_d2_13a - 0.50 * (diff_d2_13b - diff_d2_13a);
	double xder13	= diff_d1_13 * dx1_dt +	diff_d2_13 * dx2_dt;
	check("2v f_table2V value      ", val13, res13,	1e-12);
	check("2v f_table2V d_dt       ", der13, xder13,	1e-11);
    val13 = P2.value(x1, x2);
    der13 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v f_table2V value      ", val13, res13,	1e-12);
    check("2v f_table2V d_dt       ", der13, xder13,	1e-10);

	x1 = -1.6;		// (x2)
	x2 = 2.1;	// 25% less than (y0)
    double val14 = P1.value(x1, x2);
    double der14 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res14 = a20 - 0.25 * (a21 - a20);
	double diff_d1_14a	= (a20 - a10) / (m2 - m1);
	double diff_d1_14b	= (a21 - a11) / (m2 - m1);
	double diff_d1_14	= diff_d1_14a - 0.25 * (diff_d1_14b - diff_d1_14a);
	double diff_d2_14	= (a21 - a20) / (p1 - p0);
	double xder14	= diff_d1_14 * dx1_dt + diff_d2_14 * dx2_dt;
	check("2v-151 f_table2V value      ", val14, res14,	1e-12);
	check("2v-153 f_table2V d_dt       ", der14, xder14,	1e-10);
    val14 = P2.value(x1, x2);
    der14 = P2.d_dt(x1, x2, dx1_dt, dx2_dt);
    check("2v-151 f_table2V value      ", val14, res14,	1e-12);
    check("2v-155 f_table2V d_dt       ", der14, xder14,	1e-9);
} // closes test2v_f_table2V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test2v_f_table2Veq() {
	// f_table2Veq of Hp[ft] as function of m [lb] and phi [deg]
	// =========================================================

    double dx1_dt = 25.0;
    double dx2_dt = 3.0;

	double m0	= -23.6;
	double m1	= -13.6;
	double m2	= -3.6;
	vec1 points1(3);
	points1.set(0, m0);
	points1.set(1, m1);
	points1.set(2, m2);

	double p0	= 12.1;
	double p1	= 32.1;
	double p2	= 52.1;
	double p3	= 72.1;
	vec1 points2(4);
	points2.set(0, p0);
	points2.set(1, p1);
	points2.set(2, p2);
	points2.set(3, p3);

	double a00	= 10;	double a01	= 12;	double a02	= 15;	double a03	= 19;
	double a10	= 20;	double a11	= 32;	double a12	= 45;	double a13	= 69;
	double a20	= 40;	double a21	= 22;	double a22	= 40;	double a23	= 139;
	vec2 Values(4, 3);
	Values.set(0, 0, a00);
	Values.set(1, 0, a01);
	Values.set(2, 0, a02);
	Values.set(3, 0, a03);
	Values.set(0, 1, a10);
	Values.set(1, 1, a11);
	Values.set(2, 1, a12);
	Values.set(3, 1, a13);
	Values.set(0, 2, a20);
	Values.set(1, 2, a21);
	Values.set(2, 2, a22);
	Values.set(3, 2, a23);

	// Create fun object
	f_table2V P1(points1, points2, Values, math::logic::lagrange_first_precompute);
	check("2v-000 f_table2V is equispaced      ", P1.get_equi1(), true);
	check("2v-000 f_table2V is equispaced      ", P1.get_equi2(), true);
	f_table2V P6(points1, points2, Values, math::logic::hermite_first);
	check("2v-000 f_table2V is equispaced      ", P6.get_equi1(), true);
	check("2v-000 f_table2V is equispaced      ", P6.get_equi2(), true);
	f_table2V P7(points1, points2, Values, math::logic::hermite_second);
	check("2v-000 f_table2V is equispaced      ", P7.get_equi1(), true);
	check("2v-000 f_table2V is equispaced      ", P7.get_equi2(), true);

	double x1 = -23.6;		// (x0)
	double x2 = 12.1;		// (y0)
	double val1 = P1.value(x1, x2);
    double der1 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res1 = a00;
	double diff_d1_1	= (a10 - a00) / (m1 - m0);
	double diff_d2_1	= (a01 - a00) / (p1 - p0);
	double xder1 = diff_d1_1 * dx1_dt +	diff_d2_1 * dx2_dt;
	check("2v f_table2Veq value    ", val1, res1,	1e-12);
	check("2v f_table2Veq d_dt     ", der1,	xder1,	1e-11);

	x1 = -11.1;		// 25% between (x1) and (x2)
	x2 = 32.1;		// (y1)
    double val2 = P1.value(x1, x2);
    double der2 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res2 = a11 + 0.25 * (a21 - a11);
	double diff_d1_2	= (a21 - a11) / (m2 - m1);
	double diff_d2_2a	= (a12 - a11) / (p2 - p1);
	double diff_d2_2b	= (a22 - a21) / (p2 - p1);
	double diff_d2_2	=  diff_d2_2a + 0.25 * (diff_d2_2b - diff_d2_2a);
	double xder2 = diff_d1_2 * dx1_dt +	diff_d2_2 * dx2_dt;
	check("2v f_table2Veq value    ", val2, res2,	1e-12);
	check("2v f_table2Veq d_dt     ", der2,	xder2,	1e-10);

	x1 = -13.6;		// (x1)
	x2 = 34.1;		// 10% between (y1) and (y2)
    double val3 = P1.value(x1, x2);
    double der3 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res3 = a11 + 0.1 * (a12 - a11);
	double diff_d1_3a	= (a21 - a11) / (m2 - m1);
	double diff_d1_3b	= (a22 - a12) / (m2 - m1);
	double diff_d1_3	= diff_d1_3a + 0.1 * (diff_d1_3b - diff_d1_3a);
	double diff_d2_3	= (a12 - a11) / (p2 - p1);
	double xder3	= diff_d1_3 * dx1_dt + diff_d2_3 * dx2_dt;
	check("2v f_table2Veq value    ", val3, res3,	1e-12);
	check("2v f_table2Veq d_dt     ", der3,	xder3,	1e-11);

	x1 = -11.1;		// 25% between (x1) and (x2)
	x2 = 34.1;		// 10% between (y1) and (y2)
    double val4 = P1.value(x1, x2);
    double der4 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double temp1		= a11 + 0.1 * (a12 - a11);
	double temp2		= a21 + 0.1 * (a22 - a21);
	double res4	= temp1 + 0.25 * (temp2 - temp1);
	double diff_d1_4a	= (a21 - a11) / (m2 - m1);
	double diff_d1_4b	= (a22 - a12) / (m2 - m1);
	double diff_d1_4	= diff_d1_4a + 0.1 * (diff_d1_4b - diff_d1_4a);
	double diff_d2_4a	= (a12 - a11) / (p2 - p1);
	double diff_d2_4b	= (a22 - a21) / (p2 - p1);
	double diff_d2_4	=  diff_d2_4a + 0.25 * (diff_d2_4b - diff_d2_4a);
	double xder4	= diff_d1_4 * dx1_dt + diff_d2_4 * dx2_dt;
	check("2v f_table2Veq value    ", val4, res4,	1e-12);
	check("2v f_table2Veq d_dt     ", der4,	xder4,	1e-10);

	x1 = -11.1;		// 25% between (x1) and (x2)
	x2 = 72.1;		// (y3)
    double val5 = P1.value(x1, x2);
    double der5 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res5 = a13 + 0.25 * (a23 - a13);
	double diff_d1_5	= (a23 - a13) / (m2 - m1);
	double diff_d2_5a	= (a13 - a12) / (p3 - p2);
	double diff_d2_5b	= (a23 - a22) / (p3 - p2);
	double diff_d2_5	=  diff_d2_5a + 0.25 * (diff_d2_5b - diff_d2_5a);
	double xder5	= diff_d1_5 * dx1_dt + diff_d2_5 * dx2_dt;
	check("2v f_table2Veq value    ", val5, res5,	1e-12);
	check("2v f_table2Veq d_dt     ", der5,	xder5,	1e-11);

	x1 = -3.6;		// (x2)
	x2 = 34.1;		// 10% between (y1) and (y2)
    double val6 = P1.value(x1, x2);
    double der6 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res6 = a21 + 0.1 * (a22 - a21);
	double diff_d1_6a	= (a21 - a11) / (m2 - m1);
	double diff_d1_6b	= (a22 - a12) / (m2 - m1);
	double diff_d1_6	= diff_d1_6a + 0.1 * (diff_d1_6b - diff_d1_6a);
	double diff_d2_6	= (a22 - a21) / (p2 - p1);
	double xder6	= diff_d1_6 * dx1_dt + diff_d2_6 * dx2_dt;
	check("2v f_table2Veq value    ", val6, res6,	1e-12);
	check("2v f_table2Veq d_dt     ", der6,	xder6,	1e-11);

	x1 = -3.6;		// (x2)
	x2 = 72.1;		// (y3)
    double val7 = P1.value(x1, x2);
    double der7 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res7 = a23;
	double diff_d1_7	= (a23 - a13) / (m2 - m1);
	double diff_d2_7	= (a23 - a22) / (p3 - p2);
	double xder7	= diff_d1_7 * dx1_dt + diff_d2_7 * dx2_dt;
	check("2v f_table2Veq value    ", val7, res7,	1e-12);
	check("2v f_table2Veq d_dt     ", der7,	xder7,	1e-10);

    val7 = P6.value(x1, x2);
	check("2v f_table2Veq value    ", val7, res7,	1e-12);

    val7 = P7.value(x1, x2);
	check("2v f_table2Veq value    ", val7, res7,	1e-12);

	x1 = -28.6;		// 50% less than (x0)
	x2 = 7.1;		// 25% less than (y0)
    double val8 = P1.value(x1, x2);
    double der8 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	temp1		= a00 - 0.25 * (a01 - a00);
	temp2		= a10 - 0.25 * (a11 - a10);
	double res8	= temp1 - 0.5 * (temp2 - temp1);
	double diff_d1_8a	= (a10 - a00) / (m1 - m0);
	double diff_d1_8b	= (a11 - a01) / (m1 - m0);
	double diff_d1_8	= diff_d1_8a - 0.25 * (diff_d1_8b - diff_d1_8a);
	double diff_d2_8a	= (a01 - a00) / (p1 - p0);
	double diff_d2_8b	= (a11 - a10) / (p1 - p0);
	double diff_d2_8	= diff_d2_8a -0.5 * (diff_d2_8b - diff_d2_8a);
	double xder8	= diff_d1_8 * dx1_dt + diff_d2_8 * dx2_dt;
	check("2v f_table2Veq value    ", val8, res8,	1e-12);
	check("2v f_table2Veq d_dt     ", der8,	xder8,	1e-11);

	x1 = -28.6;		// 50% less than (x0)
	x2 = 12.1;		// (y0)
    double val9 = P1.value(x1, x2);
    double der9 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res9 = a00 - 0.50 * (a10 - a00);
	double diff_d1_9	= (a10 - a00) / (m1 - m0);
	double diff_d2_9a	= (a01 - a00) / (p1 - p0);
	double diff_d2_9b	= (a11 - a10) / (p1 - p0);
	double diff_d2_9	= diff_d2_9a -0.5 * (diff_d2_9b - diff_d2_9a);
	double xder9	= diff_d1_9 * dx1_dt + diff_d2_9 * dx2_dt;
	check("2v f_table2Veq value    ", val9, res9,	1e-12);
	check("2v f_table2Veq d_dt     ", der9,	xder9,	1e-11);

	x1 = -23.6;		// (x0)
	x2 = 7.1;		// 25% less than (y0)
    double val10 = P1.value(x1, x2);
    double der10 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res10 = a00 - 0.25 * (a01 - a00);
	double diff_d1_10a	= (a10 - a00) / (m1 - m0);
	double diff_d1_10b	= (a11 - a01) / (m1 - m0);
	double diff_d1_10	= diff_d1_10a - 0.25 * (diff_d1_10b - diff_d1_10a);
	double diff_d2_10	= (a01 - a00) / (p1 - p0);
	double xder10	= diff_d1_10 * dx1_dt + diff_d2_10 * dx2_dt;
	check("2v f_table2Veq value    ", val10, res10,	1e-12);
	check("2v f_table2Veq d_dt     ", der10, xder10,	1e-11);

	x1 = -2.6;		// 10% more than (x2)
	x2 = 72.1;	// (y3)
    double val11 = P1.value(x1, x2);
    double der11 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res11 = a13 + 1.10 * (a23 - a13);
	double diff_d1_11	= (a23 - a13) / (m2 - m1);
	double diff_d2_11a	= (a13 - a12) / (p3 - p2);
	double diff_d2_11b	= (a23 - a22) / (p3 - p2);
	double diff_d2_11	= diff_d2_11a + 1.10 * (diff_d2_11b - diff_d2_11a);
	double xder11	= diff_d1_11 * dx1_dt + diff_d2_11 * dx2_dt;
	check("2v f_table2Veq value    ", val11, res11,	1e-12);
	check("2v f_table2Veq d_dt     ", der11, xder11,	1e-10);

	x1 = -3.6;		// (x2)
	x2 = 82.1;	// 50% more than (y3)
    double val12 = P1.value(x1, x2);
    double der12 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	double res12 = a22 + 1.50 * (a23 - a22);
	double diff_d1_12a	= (a22 - a12) / (m2 - m1);
	double diff_d1_12b	= (a23 - a13) / (m2 - m1);
	double diff_d1_12	= diff_d1_12a + 1.50 * (diff_d1_12b - diff_d1_12a);
	double diff_d2_12	= (a23 - a22) / (p3 - p2);
	double xder12	= diff_d1_12 * dx1_dt + diff_d2_12 * dx2_dt;
	check("2v f_table2Veq value    ", val12, res12,	1e-12);
	check("2v f_table2Veq d_dt     ", der12, xder12,	1e-11);

	x1 = -2.6;		// 10% more than (x2)
	x2 = 82.1;	// 50% more than (y3)
    double val13 = P1.value(x1, x2);
    double der13 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	temp1		= a12 + 1.50 * (a13 - a12);
	temp2		= a22 + 1.50 * (a23 - a22);
	double res13	= temp1 + 1.10 * (temp2 - temp1);
	double diff_d1_13a	= (a22 - a12) / (m2 - m1);
	double diff_d1_13b	= (a23 - a13) / (m2 - m1);
	double diff_d1_13	= diff_d1_13a + 1.50 * (diff_d1_13b - diff_d1_13a);
	double diff_d2_13a	= (a13 - a12) / (p3 - p2);
	double diff_d2_13b	= (a23 - a22) / (p3 - p2);
	double diff_d2_13	= diff_d2_13a + 1.10 * (diff_d2_13b - diff_d2_13a);
	double xder13	= diff_d1_13 * dx1_dt + diff_d2_13 * dx2_dt;
	check("2v f_table2Veq value    ", val13, res13,	1e-12);
	check("2v f_table2Veq d_dt     ", der13, xder13,	1e-10);

	x1 = -28.6;	// 50% less than (x0)
	x2 = 82.1;	// 50% more than (y3)
    double val14 = P1.value(x1, x2);
    double der14 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	temp1		= a02 + 1.50 * (a03 - a02);
	temp2		= a12 + 1.50 * (a13 - a12);
	double res14	= temp1 - 0.50 * (temp2 - temp1);
	double diff_d1_14a	= (a12 - a02) / (m1 - m0);
	double diff_d1_14b	= (a13 - a03) / (m1 - m0);
	double diff_d1_14	= diff_d1_14a + 1.50 * (diff_d1_14b - diff_d1_14a);
	double diff_d2_14a	= (a03 - a02) / (p3 - p2);
	double diff_d2_14b	= (a13 - a12) / (p3 - p2);
	double diff_d2_14	= diff_d2_14a - 0.50 * (diff_d2_14b - diff_d2_14a);
	double xder14	= diff_d1_14 * dx1_dt + diff_d2_14 * dx2_dt;
	check("2v f_table2Veq value    ", val14, res14,	1e-12);
	check("2v f_table2Veq d_dt     ", der14, xder14,	1e-11);

	x1 = -2.6;		// 10% more than (x2)
	x2 = 7.1;	// 25% less than (y0)
    double val15 =  P1.value(x1, x2);
    double der15 = P1.d_dt(x1, x2, dx1_dt, dx2_dt);
	temp1		= a10 - 0.25 * (a11 - a10);
	temp2		= a20 - 0.25 * (a21 - a20);
	double res15	= temp1 + 1.10 * (temp2 - temp1);
	double diff_d1_15a	= (a20 - a10) / (m2 - m1);
	double diff_d1_15b	= (a21 - a11) / (m2 - m1);
	double diff_d1_15	= diff_d1_15a - 0.25 * (diff_d1_15b - diff_d1_15a);
	double diff_d2_15a	= (a11 - a10) / (p1 - p0);
	double diff_d2_15b	= (a21 - a20) / (p1 - p0);
	double diff_d2_15	= diff_d2_15a + 1.10 * (diff_d2_15b - diff_d2_15a);
	double xder15	= diff_d1_15 * dx1_dt + diff_d2_15 * dx2_dt;
	check("2v f_table2Veq value    ", val15, res15,	1e-12);
	check("2v f_table2Veq d_dt     ", der15, xder15,	1e-11);
} // closes test2v_f_table2Veq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test2v_f_tabular2V() {
    double dx1_dt = 5.0;
    double dx2_dt = 3.0;

	double p0	= -33.6;
	double p1	= -13.6;
	double p2	= -3.6;
	vec1 points1(3);
	points1.set(0, p0);
	points1.set(1, p1);
	points1.set(2, p2);
	
	double d0	=  82.1;
	double d1	= 132.1;
	double d2	= 232.1;
	double d3	= 332.1;
	vec1 points2(4);
	points2.set(0, d0);
	points2.set(1, d1);
	points2.set(2, d2);
	points2.set(3, d3);

	double a00	= 10;	double a01	= 12;	double a02	= 15;	double a03 = 52;
	double a10	= 20;	double a11	= 25;	double a12	= 7;	double a13 = 27;
	double a20	= 40;	double a21	= 21;	double a22	= 35;	double a23 = 36;
	vec1 Values0(4), Values1(4), Values2(4);
	Values0.set(0, a00);
	Values0.set(1, a01);
	Values0.set(2, a02);
	Values0.set(3, a03);
	Values1.set(0, a10);
	Values1.set(1, a11);
	Values1.set(2, a12);
	Values1.set(3, a13);
	Values2.set(0, a20);
	Values2.set(1, a21);
    Values2.set(2, a22);
	Values2.set(3, a23);
	
	// Create fun object
	f_table1V h0(points2, Values0, math::logic::lagrange_first);
	f_table1V h1(points2, Values1, math::logic::lagrange_first);
	f_table1V h2(points2, Values2, math::logic::lagrange_first);
	std::vector<f_table1V*> tables(3);
	tables[0] = &h0;
	tables[1] = &h1;
	tables[2] = &h2;
	f_tabular2V P1(points1, tables, math::logic::lagrange_first);

	double x2 = -33.6;	// p0
	double x1 = 82.1;	// d0
	double val1 = P1.value(x2, x1);
    double der1 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res1 = a00;
	double diff_d1_1	= (a10 - a00) / (p1 - p0);
	double diff_d2_1	= (a01 - a00) / (d1 - d0);
	double xder1	= diff_d1_1 * dx2_dt + diff_d2_1 * dx1_dt;
	check("2v f_tabular2V value      ", val1, res1,	1e-12);
	check("2v f_tabular2V d_dt       ", der1, xder1,1e-10);

	x2 = -13.6;	// p1
	x1 = 82.1;	// d0
    double val2 = P1.value(x2, x1);
    double der2 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res2 = a10;
	double diff_d1_2	= (a20 - a10) / (p2 - p1);
	double diff_d2_2	= (a11 - a10) / (d1 - d0);
	double xder2	= diff_d1_2 * dx2_dt + diff_d2_2 * dx1_dt;
	check("2v f_tabular2V value      ", val2, res2,	1e-12);
	check("2v f_tabular2V d_dt       ", der2, xder2,	1e-9);
	
	x2 = -11.1;	// 25% between p1 and p2
	x1 = 82.1;	// d0
    double val3 = P1.value(x2, x1);
    double der3 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res3 = a10 + 0.25 * (a20 - a10);
	double diff_d1_3	= (a20 - a10) / (p2 - p1);
	double diff_d2_3a	= (a11 - a10) / (d1 - d0);
	double diff_d2_3b	= (a21 - a20) / (d1 - d0);
	double diff_d2_3	= diff_d2_3a + 0.25 * (diff_d2_3b - diff_d2_3a);
	double xder3	= diff_d1_3 * dx2_dt + diff_d2_3 * dx1_dt;
	check("2v f_tabular2V value      ", val3, res3,	1e-12);
	check("2v f_tabular2V d_dt       ", der3, xder3,	1e-9);

	x2 = -13.6;	// p1
	x1 = 92.1;	// 20% between d0 and d1
    double val4 = P1.value(x2, x1);
    double der4 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res4 = a10 + 0.2 * (a11 - a10);
	double diff_d1_4a	= (a20 - a10) / (p2 - p1);
	double diff_d1_4b	= (a21 - a11) / (p2 - p1);
	double diff_d1_4	= diff_d1_4a + 0.2 * (diff_d1_4b - diff_d1_4a);
	double diff_d2_4	= (a11 - a10) / (d1 - d0);
	double xder4	= diff_d1_4 * dx2_dt + diff_d2_4 * dx1_dt;
	check("2v f_tabular2V value      ", val4, res4,	1e-12);
	check("2v f_tabular2V d_dt       ", der4, xder4,	1e-10);

	x2 = -11.1;	// 25% between p1 and p2
	x1 = 92.1;	// 20% between d0 and d1
    double val5 = P1.value(x2, x1);
    double der5 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double temp1	= a10 + 0.2  * (a11 - a10);
	double temp2	= a20 + 0.2  * (a21 - a20);
	double res5		= temp1 + 0.25 * (temp2 - temp1);
	double diff_d1_5a	= (a20 - a10) / (p2 - p1);
	double diff_d1_5b	= (a21 - a11) / (p2 - p1);
	double diff_d1_5	= diff_d1_5a + 0.2 * (diff_d1_5b - diff_d1_5a);
	double diff_d2_5a	= (a11 - a10) / (d1 - d0);
	double diff_d2_5b	= (a21 - a20) / (d1 - d0);
	double diff_d2_5	= diff_d2_5a + 0.25 * (diff_d2_5b - diff_d2_5a);
	double xder5	= diff_d1_5 * dx2_dt +	diff_d2_5 * dx1_dt;
	check("2v f_tabular2V value      ", val5, res5,	1e-12);
	check("2v f_tabular2V d_dt       ", der5, xder5,	1e-10);

	x2 = -11.1;	// 25% between p1 and p2
	x1 = 332.1;  // d3
    double val7 = P1.value(x2, x1);
    double der7 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res7       	= a13 + 0.25 * (a23 - a13);
	double diff_d1_7	= (a23 - a13) / (p2 - p1);
	double diff_d2_7a	= (a13 - a12) / (d3 - d2);
	double diff_d2_7b	= (a23 - a22) / (d3 - d2);
	double diff_d2_7	= diff_d2_7a + 0.25 * (diff_d2_7b - diff_d2_7a);
	double xder7	= diff_d1_7 * dx2_dt +	diff_d2_7 * dx1_dt;
	check("2v f_tabular2V value      ", val7, res7,	1e-12);
	check("2v f_tabular2V d_dt       ", der7, xder7,	1e-9);

	x2 = -3.6;	// p2
	x1 = 92.1;	// 20% between d0 and d1
    double val8 = P1.value(x2, x1);
    double der8 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res8			= a20 + 0.2 * (a21 - a20);
	double diff_d1_8a	= (a20 - a10) / (p2 - p1);
	double diff_d1_8b	= (a21 - a11) / (p2 - p1);
	double diff_d1_8	= diff_d1_8a + 0.2 * (diff_d1_8b - diff_d1_8a);
	double diff_d2_8	= (a21 - a20) / (d1 - d0);
	double xder8	= diff_d1_8 * dx2_dt +	diff_d2_8 * dx1_dt;
	check("2v f_tabular2V value      ", val8, res8,	1e-12);
	check("2v f_tabular2V d_dt       ", der8, xder8,	1e-9);

	x2 = -3.6;	// p2
	x1 = 332.1;	// d3
    double val10 = P1.value(x2, x1);
    double der10 = P1.d_dt(x2, x1, dx2_dt, dx1_dt);
	double res10	= a23;
	double diff_d1_10	= (a23 - a13) / (p2 - p1);
	double diff_d2_10	= (a23 - a22) / (d3 - d2);
	double xder10	= diff_d1_10 * dx2_dt +	diff_d2_10 * dx1_dt;
	check("2v f_tabular2V value      ", val10, res10,	1e-12);
	check("2v f_tabular2V d_dt       ", der10, xder10,	1e-8);
} // closes test2v_f_tabular2V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test3v_f_lineal_triple() {

    double x3 = 120.0;
    double x2 = 2500.0;
    double x1 = constant::PI()/4;
    double dx3_dt = 5.0;
    double dx2_dt = 25.0;
    double dx1_dt = 3.0;

	double f0	 = 1.;
	double f1x	 = 2.;
	double f1y	 = 3.;
	double f1z   = 4.;
	double f1xy  = 5.;
	double f1xz	 = 6.;
	double f1yz  = 7.;
	double f1xyz = 8.;

	f_lineal_triple P1(f0, f1x, f1y, f1z,f1xy, f1xz, f1yz, f1xyz);

	// Obtain value
	double val1 = P1.value(x3, x2, x1);
	double res1	= f0 + f1x * x3 + f1y * x2 + f1z * x1 + f1xy * x3 * x2
					 + f1xz * x3 *  x1 + f1yz *   x2 *  x1 + f1xyz*   x3 *   x2	* x1;
	check("3v f_lineal_triple value_", val1, res1,	1e-4);
    double der1 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double xder1 = (f1x + f1xy * x2 + f1xz * x1 + f1xyz * x2 * x1) * dx3_dt +
                   (f1y + f1xy * x3 + f1yz * x1 + f1xyz * x3 * x1) * dx2_dt +
                   (f1z + f1xz * x3 + f1yz * x2 + f1xyz * x3 * x2) * dx1_dt;
	check("3v f_lineal_triple d_dt  ", der1,	xder1,	1e-4);
} // closes test3v_f_lineal_triple

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test3v_f_table3V() {
	double dx1_dt = 5.0;
    double dx2_dt = 3.0;
    double dx3_dt = 25.0;

	double m0	= 1200.0;
	double m1	= 1300.0;
	vec1 points1(2);
	points1.set(0, m0);
	points1.set(1, m1);

	double p0	= -33.6;
	double p1	= -13.6;
	double p2	= -3.6;
	vec1 points2(3);
	points2.set(0, p0);
	points2.set(1, p1);
	points2.set(2, p2);

	double d0	=  82.1;
	double d1	= 132.1;
	double d2	= 232.1;
	double d3	= 332.1;
	vec1 points3(4);
	points3.set(0, d0);
	points3.set(1, d1);
	points3.set(2, d2);
	points3.set(3, d3);

	double a000	= 10;	double a001	= 12;	double a002	= 15;	double a003 = 52;
	double a010	= 20;	double a011	= 25;	double a012	= 7;	double a013 = 27;
	double a020	= 40;	double a021	= 21;	double a022	= 35;	double a023 = 36;
	double a100	= -100;	double a101	= -120;	double a102	= -150;	double a103 = -520;
	double a110	= -200;	double a111	= -250;	double a112	= -70;	double a113 = -270;
	double a120	= -400;	double a121	= -210;	double a122	= -350;	double a123 = -360;
	vec3 VValues(4, 3, 2);
	VValues.set(0, 0, 0, a000);
	VValues.set(1, 0, 0, a001);
	VValues.set(2, 0, 0, a002);
	VValues.set(3, 0, 0, a003);
	VValues.set(0, 1, 0, a010);
	VValues.set(1, 1, 0, a011);
	VValues.set(2, 1, 0, a012);
	VValues.set(3, 1, 0, a013);
	VValues.set(0, 2, 0, a020);
	VValues.set(1, 2, 0, a021);
	VValues.set(2, 2, 0, a022);
	VValues.set(3, 2, 0, a023);
	VValues.set(0, 0, 1, a100);
    VValues.set(1, 0, 1, a101);
	VValues.set(2, 0, 1, a102);
	VValues.set(3, 0, 1, a103);
	VValues.set(0, 1, 1, a110);
	VValues.set(1, 1, 1, a111);
	VValues.set(2, 1, 1, a112);
	VValues.set(3, 1, 1, a113);
	VValues.set(0, 2, 1, a120);
	VValues.set(1, 2, 1, a121);
	VValues.set(2, 2, 1, a122);
	VValues.set(3, 2, 1, a123);

	f_table3V P1(points1, points2, points3, VValues, math::logic::lagrange_first_precompute);
	f_table3V P2(points1, points2, points3, VValues, math::logic::lagrange_first);
	f_table3V P6(points1, points2, points3, VValues, math::logic::hermite_first);

	double x3 = 1200.0;	// (x0)
	double x2 = -33.6;	// (y0)
	double x1 = 82.1;	// (z0)
    double val1 = P1.value(x3, x2, x1);
    double der1 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
    double res1 = a000;
	double diff_d1_1	= (a100 - a000) / (m1 - m0);
	double diff_d2_1	= (a010 - a000) / (p1 - p0);
	double diff_d3_1	= (a001 - a000) / (d1 - d0);
	double xder1	= diff_d1_1 * dx3_dt + diff_d2_1 * dx2_dt + diff_d3_1 * dx1_dt;
	check("3v f_table3V value      ", val1, res1,	1e-12);
	check("3v f_table3V d_dt       ", der1,	xder1,	1e-11);
    val1 = P2.value(x3, x2, x1);
    der1 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val1, res1,1e-12);
	check("3v f_table3V d_dt       ", der1,	xder1,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -13.6;	// (y1)
	x1 = 82.1;			// (z0)
    double val2 = P1.value(x3, x2, x1);
    double der2 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res2 = a010;
	double diff_d1_2	= (a110 - a010) / (m1 - m0);
	double diff_d2_2	= (a020 - a010) / (p2 - p1);
	double diff_d3_2	= (a011 - a010) / (d1 - d0);
	double xder2	= diff_d1_2 * dx3_dt + diff_d2_2 * dx2_dt + diff_d3_2 * dx1_dt;
	check("3v f_table3V value      ", val2, res2,	1e-12);
	check("3v f_table3V d_dt       ", der2,	xder2,	1e-11);
    val2 = P2.value(x3, x2, x1);
    der2 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val2, res2,1e-12);
	check("3v f_table3V d_dt       ", der2,	xder2,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 82.1;			// (z0)
    double val3 = P1.value(x3, x2, x1);
    double der3 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res3 = a010 + 0.25 * (a020 - a010);
	double diff_d1_3a	= (a110 - a010) / (m1 - m0);
	double diff_d1_3b	= (a120 - a020) / (m1 - m0);
	double diff_d1_3	= diff_d1_3a + 0.25 * (diff_d1_3b - diff_d1_3a);
	double diff_d2_3	= (a020 - a010) / (p2 - p1);
	double diff_d3_3a	= (a011 - a010) / (d1 - d0);
	double diff_d3_3b	= (a021 - a020) / (d1 - d0);
	double diff_d3_3	= diff_d3_3a + 0.25 * (diff_d3_3b - diff_d3_3a);
	double xder3	= diff_d1_3 * dx3_dt + diff_d2_3 * dx2_dt + diff_d3_3 * dx1_dt;
	check("3v f_table3V value      ", val3, res3,	1e-12);
	check("3v f_table3V d_dt       ", der3,	xder3,	1e-11);
    val3 = P2.value(x3, x2, x1);
    der3 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val3, res3,1e-12);
	check("3v f_table3V d_dt       ", der3,	xder3,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -13.6;	// (y1)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val4 = P1.value(x3, x2, x1);
    double der4 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res4 = a010 + 0.2 * (a011 - a010);
	double diff_d1_4a	= (a110 - a010) / (m1 - m0);
	double diff_d1_4b	= (a111 - a011) / (m1 - m0);
	double diff_d1_4	= diff_d1_4a + 0.2 * (diff_d1_4b - diff_d1_4a);
	double diff_d2_4a	= (a020 - a010) / (p2 - p1);
	double diff_d2_4b	= (a021 - a011) / (p2 - p1);
	double diff_d2_4	= diff_d2_4a + 0.2 * (diff_d2_4b - diff_d2_4a);
	double diff_d3_4	= (a011 - a010) / (d1 - d0);
	double xder4	= diff_d1_4 * dx3_dt + diff_d2_4 * dx2_dt + diff_d3_4 * dx1_dt;
	check("3v f_table3V value      ", val4, res4,	1e-12);
	check("3v f_table3V d_dt       ", der4,	xder4,	1e-11);
    val4 = P2.value(x3, x2, x1);
    der4 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val4, res4,1e-12);
	check("3v f_table3V d_dt       ", der4,	xder4,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val5 = P1.value(x3, x2, x1);
    double der5 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double temp1	= a010 + 0.2  * (a011 - a010);
	double temp2	= a020 + 0.2  * (a021 - a020);
	double res5		= temp1 + 0.25 * (temp2 - temp1);
	double diff_d1_5a	= (a110 - a010) / (m1 - m0);
	double diff_d1_5b	= (a111 - a011) / (m1 - m0);
	double diff_d1_5c	= (a120 - a020) / (m1 - m0);
	double diff_d1_5d	= (a121 - a021) / (m1 - m0);
	double diff_d1_5ab	= diff_d1_5a + 0.2 * (diff_d1_5b - diff_d1_5a);
	double diff_d1_5cd	= diff_d1_5c + 0.2 * (diff_d1_5d - diff_d1_5c);
	double diff_d1_5	= diff_d1_5ab + 0.25 * (diff_d1_5cd - diff_d1_5ab);
	double diff_d2_5a	= (a020 - a010) / (p2 - p1);
	double diff_d2_5b	= (a021 - a011) / (p2 - p1);
	double diff_d2_5	= diff_d2_5a + 0.2 * (diff_d2_5b - diff_d2_5a);
	double diff_d3_5a	= (a011 - a010) / (d1 - d0);
	double diff_d3_5b	= (a021 - a020) / (d1 - d0);
	double diff_d3_5	= diff_d3_5a + 0.25 * (diff_d3_5b - diff_d3_5a);
	double xder5	= diff_d1_5 * dx3_dt + diff_d2_5 * dx2_dt + diff_d3_5 * dx1_dt;
	check("3v f_table3V value      ", val5, res5,	1e-12);
	check("3v f_table3V d_dt       ", der5,	xder5,	1e-11);
    val5 = P2.value(x3, x2, x1);
    der5 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val5, res5,1e-12);
	check("3v f_table3V d_dt       ", der5,	xder5,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val6 = P1.value(x3, x2, x1);
    double der6 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x6temp1a		= a010 + 0.2 * (a011 - a010);
	double x6temp1b		= a020 + 0.2 * (a021 - a020);
	double x6temp1c		= a110 + 0.2 * (a111 - a110);
	double x6temp1d		= a120 + 0.2 * (a121 - a120);
	double x6temp1ab	= x6temp1a + 0.25 * (x6temp1b - x6temp1a);
	double x6temp1cd	= x6temp1c + 0.25 * (x6temp1d - x6temp1c);
	double res6			= x6temp1ab + 0.5 * (x6temp1cd - x6temp1ab);
	double diff_d1_6a	= (a110 - a010) / (m1 - m0);
	double diff_d1_6b	= (a111 - a011) / (m1 - m0);
	double diff_d1_6c	= (a120 - a020) / (m1 - m0);
	double diff_d1_6d	= (a121 - a021) / (m1 - m0);
	double diff_d1_6ab	= diff_d1_6a + 0.2 * (diff_d1_6b - diff_d1_6a);
	double diff_d1_6cd	= diff_d1_6c + 0.2 * (diff_d1_6d - diff_d1_6c);
	double diff_d1_6	= diff_d1_6ab + 0.25 * (diff_d1_6cd - diff_d1_6ab);
	double diff_d2_6a	= (a020 - a010) / (p2 - p1);
	double diff_d2_6b	= (a021 - a011) / (p2 - p1);
	double diff_d2_6c	= (a120 - a110) / (p2 - p1);
	double diff_d2_6d	= (a121 - a111) / (p2 - p1);
	double diff_d2_6ac	= diff_d2_6a + 0.5 * (diff_d2_6c - diff_d2_6a);
	double diff_d2_6bd	= diff_d2_6b + 0.5 * (diff_d2_6d - diff_d2_6b);
	double diff_d2_6	= diff_d2_6ac + 0.2 * (diff_d2_6bd - diff_d2_6ac);
	double diff_d3_6a	= (a011 - a010) / (d1 - d0);
	double diff_d3_6b	= (a021 - a020) / (d1 - d0);
	double diff_d3_6c	= (a111 - a110) / (d1 - d0);
	double diff_d3_6d	= (a121 - a120) / (d1 - d0);
	double diff_d3_6ac	= diff_d3_6a + 0.5 * (diff_d3_6c - diff_d3_6a);
	double diff_d3_6bd	= diff_d3_6b + 0.5 * (diff_d3_6d - diff_d3_6b);
	double diff_d3_6	= diff_d3_6ac + 0.25 * (diff_d3_6bd - diff_d3_6ac);
	double xder6	= diff_d1_6 * dx3_dt + diff_d2_6 * dx2_dt + diff_d3_6 * dx1_dt;
	check("3v f_table3V value      ", val6, res6,	1e-12);
	check("3v f_table3V d_dt       ", der6,	xder6,	1e-09);
    val6 = P2.value(x3, x2, x1);
    der6 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val6, res6,1e-12);
	check("3v f_table3V d_dt       ", der6,	xder6,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 332.1;		// (z3)
    double val7 = P1.value(x3, x2, x1);
    double der7 = P1.d_dt( x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x7temp1ab	= a013 + 0.25 * (a023 - a013);
	double x7temp1cd	= a113 + 0.25 * (a123 - a113);
	double res7			= x7temp1ab + 0.5 * (x7temp1cd - x7temp1ab);
	double diff_d1_7a	= (a113 - a013) / (m1 - m0);
	double diff_d1_7b	= (a123 - a023) / (m1 - m0);
	double diff_d1_7	= diff_d1_7a + 0.25 * (diff_d1_7b - diff_d1_7a);
	double diff_d2_7a	= (a023 - a013) / (p2 - p1);
	double diff_d2_7b	= (a123 - a113) / (p2 - p1);
	double diff_d2_7	= diff_d2_7a + 0.5 * (diff_d2_7b - diff_d2_7a);
	double diff_d3_7a	= (a013 - a012) / (d3 - d2);
	double diff_d3_7b	= (a023 - a022) / (d3 - d2);
	double diff_d3_7c	= (a113 - a112) / (d3 - d2);
	double diff_d3_7d	= (a123 - a122) / (d3 - d2);
	double diff_d3_7ac	= diff_d3_7a + 0.5 * (diff_d3_7c - diff_d3_7a);
	double diff_d3_7bd	= diff_d3_7b + 0.5 * (diff_d3_7d - diff_d3_7b);
	double diff_d3_7	= diff_d3_7ac + 0.25 * (diff_d3_7bd - diff_d3_7ac);
	double xder7	= diff_d1_7 * dx3_dt + diff_d2_7 * dx2_dt + diff_d3_7 * dx1_dt;
	check("3v f_table3V value      ", val7, res7,	1e-11);
	check("3v f_table3V d_dt       ", der7,	xder7,	1e-10);
    val7 = P2.value(x3, x2, x1);
    der7 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val7, res7,1e-12);
	check("3v f_table3V d_dt       ", der7,	xder7,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val8 = P1.value(x3, x2, x1);
    double der8 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x8temp1b		= a020 + 0.2 * (a021 - a020);
	double x8temp1d		= a120 + 0.2 * (a121 - a120);
	double res8			= x8temp1b + 0.5 * (x8temp1d - x8temp1b);
	double diff_d1_8a	= (a120 - a020) / (m1 - m0);
	double diff_d1_8b	= (a121 - a021) / (m1 - m0);
	double diff_d1_8	= diff_d1_8a + 0.2 * (diff_d1_8b - diff_d1_8a);
	double diff_d2_8a	= (a020 - a010) / (p2 - p1);
	double diff_d2_8b	= (a021 - a011) / (p2 - p1);
	double diff_d2_8c	= (a120 - a110) / (p2 - p1);
	double diff_d2_8d	= (a121 - a111) / (p2 - p1);
	double diff_d2_8ac	= diff_d2_8a + 0.5 * (diff_d2_8c - diff_d2_8a);
	double diff_d2_8bd	= diff_d2_8b + 0.5 * (diff_d2_8d - diff_d2_8b);
	double diff_d2_8	= diff_d2_8ac + 0.2 * (diff_d2_8bd - diff_d2_8ac);
	double diff_d3_8a	= (a021 - a020) / (d1 - d0);
	double diff_d3_8b	= (a121 - a120) / (d1 - d0);
	double diff_d3_8	= diff_d3_8a + 0.5 * (diff_d3_8b - diff_d3_8a);
	double xder8	= diff_d1_8 * dx3_dt + diff_d2_8 * dx2_dt + diff_d3_8 * dx1_dt;
	check("3v f_table3V value      ", val8, res8,	1e-12);
	check("3v f_table3V d_dt       ", der8,	xder8,	1e-10);
    val8 = P2.value(x3, x2, x1);
    der8 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val8, res8,1e-12);
	check("3v f_table3V d_dt       ", der8,	xder8,	1e-9);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val9 = P1.value(x3, x2, x1);
    double der9 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x9temp1	= a110 + 0.2  * (a111 - a110);
	double x9temp2	= a120 + 0.2  * (a121 - a120);
	double res9		= x9temp1 + 0.25 * (x9temp2 - x9temp1);
	double diff_d1_9a	= (a110 - a010) / (m1 - m0);
	double diff_d1_9b	= (a111 - a011) / (m1 - m0);
	double diff_d1_9c	= (a120 - a020) / (m1 - m0);
	double diff_d1_9d	= (a121 - a021) / (m1 - m0);
	double diff_d1_9ab	= diff_d1_9a + 0.2 * (diff_d1_9b - diff_d1_9a);
	double diff_d1_9cd	= diff_d1_9c + 0.2 * (diff_d1_9d - diff_d1_9c);
	double diff_d1_9	= diff_d1_9ab + 0.25 * (diff_d1_9cd - diff_d1_9ab);
	double diff_d2_9a	= (a120 - a110) / (p2 - p1);
	double diff_d2_9b	= (a121 - a111) / (p2 - p1);
	double diff_d2_9	= diff_d2_9a + 0.2 * (diff_d2_9b - diff_d2_9a);
	double diff_d3_9a	= (a111 - a110) / (d1 - d0);
	double diff_d3_9b	= (a121 - a120) / (d1 - d0);
	double diff_d3_9	= diff_d3_9a + 0.25 * (diff_d3_9b - diff_d3_9a);
	double xder9	= diff_d1_9 * dx3_dt + diff_d2_9 * dx2_dt + diff_d3_9 * dx1_dt;
	check("3v f_table3V value      ", val9, res9,	1e-12);
	check("3v f_table3V d_dt       ", der9,	xder9,	1e-09);
    val9 = P2.value(x3, x2, x1);
    der9 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val9, res9,1e-12);
	check("3v f_table3V d_dt       ", der9,	xder9,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 332.1;		// (z3)
    double val10 = P1.value(x3, x2, x1);
    double der10 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res10		= a023 + 0.5 * (a123 - a023);
	double diff_d1_10	= (a123 - a023) / (m1 - m0);
	double diff_d2_10a	= (a023 - a013) / (p2 - p1);
	double diff_d2_10b	= (a123 - a113) / (p2 - p1);
	double diff_d2_10	= diff_d2_10a + 0.5 * (diff_d2_10b - diff_d2_10a);
	double diff_d3_10a	= (a023 - a022) / (d3 - d2);
	double diff_d3_10b	= (a123 - a122) / (d3 - d2);
	double diff_d3_10	= diff_d3_10a + 0.5 * (diff_d3_10b - diff_d3_10a);
	double xder10	= diff_d1_10 * dx3_dt + diff_d2_10 * dx2_dt + diff_d3_10 * dx1_dt;
	check("3v f_table3V value      ", val10, res10,	1e-12);
	check("3v f_table3V d_dt       ", der10, xder10,	1e-10);
    val10 = P2.value(x3, x2, x1);
    der10 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val10, res10,1e-12);
	check("3v f_table3V d_dt       ", der10, xder10,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 332.1;		// (z3)
    double val11 = P1.value(x3, x2, x1);
    double der11 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res11		= a113 + 0.25 * (a123 - a113);
	double diff_d1_11a	= (a113 - a013) / (m1 - m0);
	double diff_d1_11b	= (a123 - a023) / (m1 - m0);
	double diff_d1_11	= diff_d1_11a + 0.25 * (diff_d1_11b - diff_d1_11a);
	double diff_d2_11	= (a123 - a113) / (p2 - p1);
	double diff_d3_11a	= (a113 - a112) / (d3 - d2);
	double diff_d3_11b	= (a123 - a122) / (d3 - d2);
	double diff_d3_11	= diff_d3_11a + 0.25 * (diff_d3_11b - diff_d3_11a);
	double xder11	= diff_d1_11 * dx3_dt + diff_d2_11 * dx2_dt + diff_d3_11 * dx1_dt;
	check("3v f_table3V value      ", val11, res11,	1e-12);
	check("3v f_table3V d_dt       ", der11, xder11,	1e-10);
    val11 = P2.value(x3, x2, x1);
    der11 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val11, res11,1e-12);
	check("3v f_table3V d_dt       ", der11, xder11,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val12 = P1.value(x3, x2, x1);
    double der12= P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res12		= a120 + 0.2 * (a121 - a120);
	double diff_d1_12a	= (a120 - a020) / (m1 - m0);
	double diff_d1_12b	= (a121 - a021) / (m1 - m0);
	double diff_d1_12	= diff_d1_12a + 0.2 * (diff_d1_12b - diff_d1_12a);
	double diff_d2_12a	= (a120 - a110) / (p2 - p1);
	double diff_d2_12b	= (a121 - a111) / (p2 - p1);
	double diff_d2_12	= diff_d2_12a + 0.2 * (diff_d2_12b - diff_d2_12a);
	double diff_d3_12   = (a121 - a120) / (d1 - d0);
	double xder12	= diff_d1_12 * dx3_dt + diff_d2_12 * dx2_dt + diff_d3_12 * dx1_dt;
	check("3v f_table3V value      ", val12, res12,	1e-12);
	check("3v f_table3V d_dt       ", der12, xder12,	1e-09);
    val12 = P2.value(x3, x2, x1);
    der12 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val12, res12,1e-12);
	check("3v f_table3V d_dt       ", der12, xder12,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 332.1;		// (z3)
    double val13 = P1.value(x3, x2, x1);
    double der13 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res13		= a123;
	double diff_d1_13	= (a123 - a023) / (m1 - m0);
	double diff_d2_13	= (a123 - a113) / (p2 - p1);
	double diff_d3_13   = (a123 - a122) / (d3 - d2);
	double xder13	= diff_d1_13 * dx3_dt + diff_d2_13 * dx2_dt + diff_d3_13 * dx1_dt;
	check("3v f_table3V value      ", val13, res13,	1e-12);
	check("3v f_table3V d_dt       ", der13, xder13,	1e-9);
    val13 = P2.value(x3, x2, x1);
    der13= P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val13, res13,1e-12);
	check("3v f_table3V d_dt       ", der13, xder13,	1e-9);

    val13 = P6.value(x3, x2, x1);
	check("3v f_table3V value      ", val13, res13,	1e-12);
	
	x3 = 1180.0;	// 20% less than (x0)
	x2 = -33.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val14 = P1.value(x3, x2, x1);
    double der14 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res14		= a000 - 0.2 * (a100 - a000);
	double diff_d1_14	= (a100 - a000) / (m1 - m0);
	double diff_d2_14a	= (a010 - a000) / (p1 - p0);
	double diff_d2_14b	= (a110 - a100) / (p1 - p0);
	double diff_d2_14	= diff_d2_14a - 0.2 * (diff_d2_14b - diff_d2_14a);
	double diff_d3_14a	= (a001 - a000) / (d1 - d0);
	double diff_d3_14b	= (a101 - a100) / (d1 - d0);
	double diff_d3_14	= diff_d3_14a - 0.2 * (diff_d3_14b - diff_d3_14a);
	double xder14	= diff_d1_14 * dx3_dt + diff_d2_14 * dx2_dt + diff_d3_14 * dx1_dt;
	check("3v f_table3V value      ", val14, res14,	1e-12);
	check("3v f_table3V d_dt       ", der14,	xder14,	1e-10);
    val14 = P2.value(x3, x2, x1);
    der14 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val14, res14,1e-12);
	check("3v f_table3V d_dt       ", der14, xder14,	1e-9);

	x3 = 1340.0;	// 40% more than (x1)
	x2 = -33.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val15 = P1.value(x3, x2, x1);
    double der15 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res15		= a000 + 1.4 * (a100 - a000);
	double diff_d1_15	= (a100 - a000) / (m1 - m0);
	double diff_d2_15a	= (a010 - a000) / (p1 - p0);
	double diff_d2_15b	= (a110 - a100) / (p1 - p0);
	double diff_d2_15	= diff_d2_15a + 1.4 * (diff_d2_15b - diff_d2_15a);
	double diff_d3_15a	= (a001 - a000) / (d1 - d0);
	double diff_d3_15b	= (a101 - a100) / (d1 - d0);
	double diff_d3_15	= diff_d3_15a + 1.4 * (diff_d3_15b - diff_d3_15a);
	double xder15	= diff_d1_15 * dx3_dt + diff_d2_15 * dx2_dt + diff_d3_15 * dx1_dt;
	check("3v f_table3V value      ", val15, res15,	1e-12);
	check("3v f_table3V d_dt       ", der15, xder15,	1e-10);
    val15 = P2.value(x3, x2, x1);
    der15 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
    check("3v f_table3V value      ", val15, res15,1e-12);
	check("3v f_table3V d_dt       ", der15, xder15,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -13.6;	// (y1)
	x1 = 332.1;	// (z3)
    double val16 = P1.value(x3, x2, x1);
    double der16 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res16		= a013 - 0.2 * (a113 - a013);
	double diff_d1_16	= (a113 - a013) / (m1 - m0);
	double diff_d2_16a	= (a023 - a013) / (p2 - p1);
	double diff_d2_16b	= (a123 - a113) / (p2 - p1);
	double diff_d2_16	= diff_d2_16a - 0.2 * (diff_d2_16b - diff_d2_16a);
	double diff_d3_16a	= (a013 - a012) / (d3 - d2);
	double diff_d3_16b	= (a113 - a112) / (d3 - d2);
	double diff_d3_16	= diff_d3_16a - 0.2 * (diff_d3_16b - diff_d3_16a);
	double xder16	= diff_d1_16 * dx3_dt + diff_d2_16 * dx2_dt + diff_d3_16 * dx1_dt;
	check("3v f_table3V value      ", val16, res16,	1e-12);
	check("3v f_table3V d_dt       ", der16, xder16,	1e-10);
    val16 = P2.value(x3, x2, x1);
    der16 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val16, res16,1e-12);
	check("3v f_table3V d_dt       ", der16, xder16,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 132.1;	// (z1)
    double val17 = P1.value(x3, x2, x1);
    double der17 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res17		= a021 - 0.2 * (a121 - a021);
	double diff_d1_17	= (a121 - a021) / (m1 - m0);
	double diff_d2_17a	= (a021 - a011) / (p2 - p1);
	double diff_d2_17b	= (a121 - a111) / (p2 - p1);
	double diff_d2_17	= diff_d2_17a - 0.2 * (diff_d2_17b - diff_d2_17a);
	double diff_d3_17a	= (a022 - a021) / (d2 - d1);
	double diff_d3_17b	= (a122 - a121) / (d2 - d1);
	double diff_d3_17	= diff_d3_17a - 0.2 * (diff_d3_17b - diff_d3_17a);
	double xder17	= diff_d1_17 * dx3_dt + diff_d2_17 * dx2_dt + diff_d3_17 * dx1_dt;
	check("3v f_table3V value      ", val17, res17,	1e-11);
	check("3v f_table3V d_dt       ", der17, xder17,	1e-10);
    val17 = P2.value(x3, x2, x1);
    der17 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val17, res17,1e-12);
	check("3v f_table3V d_dt       ", der17, xder17,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 332.1;	// (z3)
    double val18 = P1.value(x3, x2, x1);
    double der18 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res18		= a023 - 0.2 * (a123 - a023);
	double diff_d1_18	= (a123 - a023) / (m1 - m0);
	double diff_d2_18a	= (a023 - a013) / (p2 - p1);
	double diff_d2_18b	= (a123 - a113) / (p2 - p1);
	double diff_d2_18	= diff_d2_18a - 0.2 * (diff_d2_18b - diff_d2_18a);
	double diff_d3_18a	= (a023 - a022) / (d3 - d2);
	double diff_d3_18b	= (a123 - a122) / (d3 - d2);
	double diff_d3_18	= diff_d3_18a - 0.2 * (diff_d3_18b - diff_d3_18a);
	double xder18	= diff_d1_18 * dx3_dt + diff_d2_18 * dx2_dt + diff_d3_18 * dx1_dt;
	check("3v f_table3V value      ", val18, res18,	1e-11);
	check("3v f_table3V d_dt       ", der18, xder18,	1e-10);
    val18 = P2.value(x3, x2, x1);
    der18 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3V value      ", val18, res18,1e-12);
	check("3v f_table3V d_dt       ", der18, xder18,	1e-9);

} // closes test3v_f_table3V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	
void math::test::Tpre_fun::test3v_f_table3Veq() {
    double dx1_dt = 5.0;
    double dx2_dt = 3.0;
    double dx3_dt = 25.0;

	double m0	= 1200.0;
	double m1	= 1300.0;
	vec1 points1(2);
	points1.set(0, m0);
	points1.set(1, m1);

	double p0	= -23.6;
	double p1	= -13.6;
	double p2	= -3.6;
	vec1 points2(3);
	points2.set(0, p0);
	points2.set(1, p1);
	points2.set(2, p2);

	double d0	=  82.1;
	double d1	= 132.1;
	double d2	= 182.1;
	double d3	= 232.1;
	vec1 points3(4);
	points3.set(0, d0);
	points3.set(1, d1);
	points3.set(2, d2);
	points3.set(3, d3);

	double a000	= 10;	double a001	= 12;	double a002	= 15;	double a003 = 52;
	double a010	= 20;	double a011	= 25;	double a012	= 7;	double a013 = 27;
	double a020	= 40;	double a021	= 21;	double a022	= 35;	double a023 = 36;
	double a100	= -100;	double a101	= -120;	double a102	= -150;	double a103 = -520;
	double a110	= -200;	double a111	= -250;	double a112	= -70;	double a113 = -270;
	double a120	= -400;	double a121	= -210;	double a122	= -350;	double a123 = -360;
	vec3 VValues(4, 3, 2);
	VValues.set(0, 0, 0, a000);
	VValues.set(1, 0, 0, a001);
	VValues.set(2, 0, 0, a002);
	VValues.set(3, 0, 0, a003);
	VValues.set(0, 1, 0, a010);
	VValues.set(1, 1, 0, a011);
	VValues.set(2, 1, 0, a012);
	VValues.set(3, 1, 0, a013);
	VValues.set(0, 2, 0, a020);
	VValues.set(1, 2, 0, a021);
	VValues.set(2, 2, 0, a022);
	VValues.set(3, 2, 0, a023);
	VValues.set(0, 0, 1, a100);
	VValues.set(1, 0, 1, a101);
	VValues.set(2, 0, 1, a102);
	VValues.set(3, 0, 1, a103);
	VValues.set(0, 1, 1, a110);
	VValues.set(1, 1, 1, a111);
	VValues.set(2, 1, 1, a112);
	VValues.set(3, 1, 1, a113);
	VValues.set(0, 2, 1, a120);
	VValues.set(1, 2, 1, a121);
	VValues.set(2, 2, 1, a122);
	VValues.set(3, 2, 1, a123);

	// Create fun object
	f_table3V P1(points1, points2, points3, VValues, math::logic::lagrange_first_precompute);
	check("3v-000 f_table3V is equispaced      ", P1.get_equi1(), true);
	check("3v-000 f_table3V is equispaced      ", P1.get_equi2(), true);
	check("3v-000 f_table3V is equispaced      ", P1.get_equi3(), true);
	f_table3V P2(points1, points2, points3, VValues, math::logic::lagrange_first);
	check("3v-000 f_table3V is equispaced      ", P2.get_equi1(), true);
	check("3v-000 f_table3V is equispaced      ", P2.get_equi2(), true);
	check("3v-000 f_table3V is equispaced      ", P2.get_equi3(), true);
	f_table3V P6(points1, points2, points3, VValues, math::logic::hermite_first);
	check("3v-000 f_table3V is equispaced      ", P6.get_equi1(), true);
	check("3v-000 f_table3V is equispaced      ", P6.get_equi2(), true);
	check("3v-000 f_table3V is equispaced      ", P6.get_equi3(), true);


	double x3 = 1200.0;		// (x0)
	double x2 = -23.6;		// (y0)
	double x1 = 82.1;		// (z0)
    double val1 = P1.value(x3, x2, x1);
    double der1 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res1 = a000;
	double diff_d1_1	= (a100 - a000) / (m1 - m0);
	double diff_d2_1	= (a010 - a000) / (p1 - p0);
	double diff_d3_1	= (a001 - a000) / (d1 - d0);
	double xder1	= diff_d1_1 * dx3_dt + diff_d2_1 * dx2_dt + diff_d3_1 * dx1_dt;
	check("3v f_table3Veq value    ", val1, res1,	1e-12);
	check("3v f_table3Veq d_dt     ", der1,	xder1,	1e-11);
    val1 = P2.value(x3, x2, x1);
	der1 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val1, res1,1e-12);
	check("3v f_table3Veq d_dt     ", der1,	xder1,	1e-9);

	x3 = 1200.0;	// (x0)
	x2 = -13.6;	// (y1)
	x1 = 82.1;	// (z0)
    double val2 = P1.value(x3, x2, x1);
    double der2 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res2 = a010;
	double diff_d1_2	= (a110 - a010) / (m1 - m0);
	double diff_d2_2	= (a020 - a010) / (p2 - p1);
	double diff_d3_2	= (a011 - a010) / (d1 - d0);
	double xder2	= diff_d1_2 * dx3_dt + 	diff_d2_2 * dx2_dt + diff_d3_2 * dx1_dt;
	check("3v f_table3Veq value    ", val2, res2,	1e-12);
	check("3v f_table3Veq d_dt     ", der2,	xder2,	1e-11);
    val2 = P2.value(x3, x2, x1);
	der2 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val2, res2,1e-12);
	check("3v f_table3Veq d_dt     ", der2,	xder2,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 82.1;			// (z0)
    double val3 = P1.value(x3, x2, x1);
	double der3 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res3 = a010 + 0.25 * (a020 - a010);
	double diff_d1_3a	= (a110 - a010) / (m1 - m0);
	double diff_d1_3b	= (a120 - a020) / (m1 - m0);
	double diff_d1_3	= diff_d1_3a + 0.25 * (diff_d1_3b - diff_d1_3a);
	double diff_d2_3	= (a020 - a010) / (p2 - p1);
	double diff_d3_3a	= (a011 - a010) / (d1 - d0);
	double diff_d3_3b	= (a021 - a020) / (d1 - d0);
	double diff_d3_3	= diff_d3_3a + 0.25 * (diff_d3_3b - diff_d3_3a);
	double xder3	= diff_d1_3 * dx3_dt + diff_d2_3 * dx2_dt + diff_d3_3 * dx1_dt;
	check("3v f_table3Veq value    ", val3, res3,	1e-12);
	check("3v f_table3Veq d_dt     ", der3,	xder3,	1e-11);
    val3 = P2.value(x3, x2, x1);
	der3 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val3, res3,1e-12);
	check("3v f_table3Veq d_dt     ", der3,	xder3,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -13.6;	// (y1)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val4 = P1.value(x3, x2, x1);
	double der4 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res4 = a010 + 0.2 * (a011 - a010);
	double diff_d1_4a	= (a110 - a010) / (m1 - m0);
	double diff_d1_4b	= (a111 - a011) / (m1 - m0);
	double diff_d1_4	= diff_d1_4a + 0.2 * (diff_d1_4b - diff_d1_4a);
	double diff_d2_4a	= (a020 - a010) / (p2 - p1);
	double diff_d2_4b	= (a021 - a011) / (p2 - p1);
	double diff_d2_4	= diff_d2_4a + 0.2 * (diff_d2_4b - diff_d2_4a);
	double diff_d3_4	= (a011 - a010) / (d1 - d0);
	double xder4	= diff_d1_4 * dx3_dt + diff_d2_4 * dx2_dt + diff_d3_4 * dx1_dt;
	check("3v f_table3Veq value    ", val4, res4,	1e-12);
	check("3v f_table3Veq d_dt     ", der4,	xder4,	1e-11);
    val4 = P2.value(x3, x2, x1);
	der4 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val4, res4,1e-12);
	check("3v f_table3Veq d_dt     ", der4,	xder4,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val5 = P1.value(x3, x2, x1);
	double der5 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double temp1	= a010 + 0.2  * (a011 - a010);
	double temp2	= a020 + 0.2  * (a021 - a020);
	double res5		= temp1 + 0.25 * (temp2 - temp1);
	double diff_d1_5a	= (a110 - a010) / (m1 - m0);
	double diff_d1_5b	= (a111 - a011) / (m1 - m0);
	double diff_d1_5c	= (a120 - a020) / (m1 - m0);
	double diff_d1_5d	= (a121 - a021) / (m1 - m0);
	double diff_d1_5ab	= diff_d1_5a + 0.2 * (diff_d1_5b - diff_d1_5a);
	double diff_d1_5cd	= diff_d1_5c + 0.2 * (diff_d1_5d - diff_d1_5c);
	double diff_d1_5	= diff_d1_5ab + 0.25 * (diff_d1_5cd - diff_d1_5ab);
	double diff_d2_5a	= (a020 - a010) / (p2 - p1);
	double diff_d2_5b	= (a021 - a011) / (p2 - p1);
	double diff_d2_5	= diff_d2_5a + 0.2 * (diff_d2_5b - diff_d2_5a);
	double diff_d3_5a	= (a011 - a010) / (d1 - d0);
	double diff_d3_5b	= (a021 - a020) / (d1 - d0);
	double diff_d3_5	= diff_d3_5a + 0.25 * (diff_d3_5b - diff_d3_5a);
	double xder5	= diff_d1_5 * dx3_dt + diff_d2_5 * dx2_dt + diff_d3_5 * dx1_dt;
	check("3v f_table3Veq value    ", val5, res5,	1e-12);
	check("3v f_table3Veq d_dt     ", der5,	xder5,	1e-11);
    val5 = P2.value(x3, x2, x1);
	der5 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val5, res5,1e-12);
	check("3v f_table3Veq d_dt     ", der5,	xder5,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val6 = P1.value(x3, x2, x1);
	double der6 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x6temp1a		= a010 + 0.2 * (a011 - a010);
	double x6temp1b		= a020 + 0.2 * (a021 - a020);
	double x6temp1c		= a110 + 0.2 * (a111 - a110);
	double x6temp1d		= a120 + 0.2 * (a121 - a120);
	double x6temp1ab	= x6temp1a + 0.25 * (x6temp1b - x6temp1a);
	double x6temp1cd	= x6temp1c + 0.25 * (x6temp1d - x6temp1c);
	double res6			= x6temp1ab + 0.5 * (x6temp1cd - x6temp1ab);
	double diff_d1_6a	= (a110 - a010) / (m1 - m0);
	double diff_d1_6b	= (a111 - a011) / (m1 - m0);
	double diff_d1_6c	= (a120 - a020) / (m1 - m0);
	double diff_d1_6d	= (a121 - a021) / (m1 - m0);
	double diff_d1_6ab	= diff_d1_6a + 0.2 * (diff_d1_6b - diff_d1_6a);
	double diff_d1_6cd	= diff_d1_6c + 0.2 * (diff_d1_6d - diff_d1_6c);
	double diff_d1_6	= diff_d1_6ab + 0.25 * (diff_d1_6cd - diff_d1_6ab);
	double diff_d2_6a	= (a020 - a010) / (p2 - p1);
	double diff_d2_6b	= (a021 - a011) / (p2 - p1);
	double diff_d2_6c	= (a120 - a110) / (p2 - p1);
	double diff_d2_6d	= (a121 - a111) / (p2 - p1);
	double diff_d2_6ac	= diff_d2_6a + 0.5 * (diff_d2_6c - diff_d2_6a);
	double diff_d2_6bd	= diff_d2_6b + 0.5 * (diff_d2_6d - diff_d2_6b);
	double diff_d2_6	= diff_d2_6ac + 0.2 * (diff_d2_6bd - diff_d2_6ac);
	double diff_d3_6a	= (a011 - a010) / (d1 - d0);
	double diff_d3_6b	= (a021 - a020) / (d1 - d0);
	double diff_d3_6c	= (a111 - a110) / (d1 - d0);
	double diff_d3_6d	= (a121 - a120) / (d1 - d0);
	double diff_d3_6ac	= diff_d3_6a + 0.5 * (diff_d3_6c - diff_d3_6a);
	double diff_d3_6bd	= diff_d3_6b + 0.5 * (diff_d3_6d - diff_d3_6b);
	double diff_d3_6	= diff_d3_6ac + 0.25 * (diff_d3_6bd - diff_d3_6ac);
	double xder6	= diff_d1_6 * dx3_dt + diff_d2_6 * dx2_dt + diff_d3_6 * dx1_dt;
	check("3v f_table3Veq value    ", val6, res6,	1e-12);
	check("3v f_table3Veq d_dt     ", der6,	xder6,	1e-9);
	val6 = P2.value(x3, x2, x1);
	der6 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val6, res6,1e-12);
	check("3v f_table3Veq d_dt     ", der6,	xder6,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 232.1;		// (z3)
    double val7 = P1.value(x3, x2, x1);
	double der7 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x7temp1ab	= a013 + 0.25 * (a023 - a013);
	double x7temp1cd	= a113 + 0.25 * (a123 - a113);
	double res7			= x7temp1ab + 0.5 * (x7temp1cd - x7temp1ab);
	double diff_d1_7a	= (a113 - a013) / (m1 - m0);
	double diff_d1_7b	= (a123 - a023) / (m1 - m0);
	double diff_d1_7	= diff_d1_7a + 0.25 * (diff_d1_7b - diff_d1_7a);
	double diff_d2_7a	= (a023 - a013) / (p2 - p1);
	double diff_d2_7b	= (a123 - a113) / (p2 - p1);
	double diff_d2_7	= diff_d2_7a + 0.5 * (diff_d2_7b - diff_d2_7a);
	double diff_d3_7a	= (a013 - a012) / (d3 - d2);
	double diff_d3_7b	= (a023 - a022) / (d3 - d2);
	double diff_d3_7c	= (a113 - a112) / (d3 - d2);
	double diff_d3_7d	= (a123 - a122) / (d3 - d2);
	double diff_d3_7ac	= diff_d3_7a + 0.5 * (diff_d3_7c - diff_d3_7a);
	double diff_d3_7bd	= diff_d3_7b + 0.5 * (diff_d3_7d - diff_d3_7b);
	double diff_d3_7	= diff_d3_7ac + 0.25 * (diff_d3_7bd - diff_d3_7ac);
	double xder7	= diff_d1_7 * dx3_dt + diff_d2_7 * dx2_dt + diff_d3_7 * dx1_dt;
	check("3v f_table3Veq value    ", val7, res7,	1e-12);
	check("3v f_table3Veq d_dt     ", der7,	xder7,	1e-10);
    val7 = P2.value(x3, x2, x1);
	der7 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val7, res7,1e-12);
	check("3v f_table3Veq d_dt     ", der7,	xder7,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val8 = P1.value(x3, x2, x1);
	double der8 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x8temp1b		= a020 + 0.2 * (a021 - a020);
	double x8temp1d		= a120 + 0.2 * (a121 - a120);
	double res8			= x8temp1b + 0.5 * (x8temp1d - x8temp1b);
	double diff_d1_8a	= (a120 - a020) / (m1 - m0);
	double diff_d1_8b	= (a121 - a021) / (m1 - m0);
	double diff_d1_8	= diff_d1_8a + 0.2 * (diff_d1_8b - diff_d1_8a);
	double diff_d2_8a	= (a020 - a010) / (p2 - p1);
	double diff_d2_8b	= (a021 - a011) / (p2 - p1);
	double diff_d2_8c	= (a120 - a110) / (p2 - p1);
	double diff_d2_8d	= (a121 - a111) / (p2 - p1);
	double diff_d2_8ac	= diff_d2_8a + 0.5 * (diff_d2_8c - diff_d2_8a);
	double diff_d2_8bd	= diff_d2_8b + 0.5 * (diff_d2_8d - diff_d2_8b);
	double diff_d2_8	= diff_d2_8ac + 0.2 * (diff_d2_8bd - diff_d2_8ac);
	double diff_d3_8a	= (a021 - a020) / (d1 - d0);
	double diff_d3_8b	= (a121 - a120) / (d1 - d0);
	double diff_d3_8	= diff_d3_8a + 0.5 * (diff_d3_8b - diff_d3_8a);
	double xder8	= diff_d1_8 * dx3_dt + diff_d2_8 * dx2_dt + diff_d3_8 * dx1_dt;
	check("3v f_table3Veq value    ", val8, res8,	1e-12);
	check("3v f_table3Veq d_dt     ", der8,	xder8,	1e-10);
    val8 = P2.value(x3, x2, x1);
	der8 = P2.d_dt( x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val8, res8,1e-12);
	check("3v f_table3Veq d_dt     ", der8,	xder8,	1e-9);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val9 = P1.value(x3, x2, x1);
	double der9 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x9temp1	= a110 + 0.2  * (a111 - a110);
	double x9temp2	= a120 + 0.2  * (a121 - a120);
	double res9		= x9temp1 + 0.25 * (x9temp2 - x9temp1);
	double diff_d1_9a	= (a110 - a010) / (m1 - m0);
	double diff_d1_9b	= (a111 - a011) / (m1 - m0);
	double diff_d1_9c	= (a120 - a020) / (m1 - m0);
	double diff_d1_9d	= (a121 - a021) / (m1 - m0);
	double diff_d1_9ab	= diff_d1_9a + 0.2 * (diff_d1_9b - diff_d1_9a);
	double diff_d1_9cd	= diff_d1_9c + 0.2 * (diff_d1_9d - diff_d1_9c);
	double diff_d1_9	= diff_d1_9ab + 0.25 * (diff_d1_9cd - diff_d1_9ab);
	double diff_d2_9a	= (a120 - a110) / (p2 - p1);
	double diff_d2_9b	= (a121 - a111) / (p2 - p1);
	double diff_d2_9	= diff_d2_9a + 0.2 * (diff_d2_9b - diff_d2_9a);
	double diff_d3_9a	= (a111 - a110) / (d1 - d0);
	double diff_d3_9b	= (a121 - a120) / (d1 - d0);
	double diff_d3_9	= diff_d3_9a + 0.25 * (diff_d3_9b - diff_d3_9a);
	double xder9	= diff_d1_9 * dx3_dt + diff_d2_9 * dx2_dt + diff_d3_9 * dx1_dt;
	check("3v f_table3Veq value    ", val9, res9,	1e-12);
	check("3v f_table3Veq d_dt     ", der9,	xder9,	1e-9);
    val9 = P2.value(x3, x2, x1);
	der9 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val9, res9,1e-12);
	check("3v f_table3Veq d_dt     ", der9,	xder9,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 232.1;		// (z3)
	double val10 = P1.value(x3, x2, x1);
	double der10 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res10		= a023 + 0.5 * (a123 - a023);
	double diff_d1_10	= (a123 - a023) / (m1 - m0);
	double diff_d2_10a	= (a023 - a013) / (p2 - p1);
	double diff_d2_10b	= (a123 - a113) / (p2 - p1);
	double diff_d2_10	= diff_d2_10a + 0.5 * (diff_d2_10b - diff_d2_10a);
	double diff_d3_10a	= (a023 - a022) / (d3 - d2);
	double diff_d3_10b	= (a123 - a122) / (d3 - d2);
	double diff_d3_10	= diff_d3_10a + 0.5 * (diff_d3_10b - diff_d3_10a);
	double xder10	= diff_d1_10 * dx3_dt + diff_d2_10 * dx2_dt + diff_d3_10 * dx1_dt;
	check("3v f_table3Veq value    ", val10, res10,	1e-12);
	check("3v f_table3Veq d_dt     ", der10, xder10,	1e-10);
    val10 = P2.value(x3, x2, x1);
	der10 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val10, res10,1e-12);
	check("3v f_table3Veq d_dt     ", der10, xder10,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 232.1;		// (z3)
    double val11 = P1.value(x3, x2, x1);
	double der11 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res11		= a113 + 0.25 * (a123 - a113);
	double diff_d1_11a	= (a113 - a013) / (m1 - m0);
	double diff_d1_11b	= (a123 - a023) / (m1 - m0);
	double diff_d1_11	= diff_d1_11a + 0.25 * (diff_d1_11b - diff_d1_11a);
	double diff_d2_11	= (a123 - a113) / (p2 - p1);
	double diff_d3_11a	= (a113 - a112) / (d3 - d2);
	double diff_d3_11b	= (a123 - a122) / (d3 - d2);
	double diff_d3_11	= diff_d3_11a + 0.25 * (diff_d3_11b - diff_d3_11a);
	double xder11	= diff_d1_11 * dx3_dt + diff_d2_11 * dx2_dt + diff_d3_11 * dx1_dt;
	check("3v f_table3Veq value    ", val11, res11,	1e-12);
	check("3v f_table3Veq d_dt     ", der11, xder11,	1e-10);
    val11 = P2.value(x3, x2, x1);
	der11 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val11, res11,1e-12);
	check("3v f_table3Veq d_dt     ", der11, xder11,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val12 = P1.value(x3, x2, x1);
	double der12 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res12		= a120 + 0.2 * (a121 - a120);
	double diff_d1_12a	= (a120 - a020) / (m1 - m0);
	double diff_d1_12b	= (a121 - a021) / (m1 - m0);
	double diff_d1_12	= diff_d1_12a + 0.2 * (diff_d1_12b - diff_d1_12a);
	double diff_d2_12a	= (a120 - a110) / (p2 - p1);
	double diff_d2_12b	= (a121 - a111) / (p2 - p1);
	double diff_d2_12	= diff_d2_12a + 0.2 * (diff_d2_12b - diff_d2_12a);
	double diff_d3_12   = (a121 - a120) / (d1 - d0);
	double xder12	= diff_d1_12 * dx3_dt + diff_d2_12 * dx2_dt + diff_d3_12 * dx1_dt;
	check("3v f_table3Veq value    ", val12, res12,	1e-12);
	check("3v f_table3Veq d_dt     ", der12, xder12,	1e-9);
    val12 = P2.value(x3, x2, x1);
	der12 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val12, res12,1e-12);
	check("3v f_table3Veq d_dt     ", der12, xder12,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 232.1;		// (z3)
    double val13 = P1.value(x3, x2, x1);
	double der13 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res13		= a123;
	double diff_d1_13	= (a123 - a023) / (m1 - m0);
	double diff_d2_13	= (a123 - a113) / (p2 - p1);
	double diff_d3_13   = (a123 - a122) / (d3 - d2);
	double xder13	= diff_d1_13 * dx3_dt + diff_d2_13 * dx2_dt + diff_d3_13 * dx1_dt;
	check("3v f_table3Veq value    ", val13, res13,	1e-12);
	check("3v f_table3Veq d_dt     ", der13, xder13,	1e-10);
    val13 = P2.value(x3, x2, x1);
	der13 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val13, res13,1e-12);
	check("3v f_table3Veq d_dt     ", der13, xder13,	1e-9);

    val13 = P6.value(x3, x2, x1);
	check("3v f_table3Veq value    ", val13, res13,	1e-12);

	x3 = 1180.0;	// 20% less than (x0)
	x2 = -23.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val14 = P1.value(x3, x2, x1);
	double der14 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res14		= a000 - 0.2 * (a100 - a000);
	double diff_d1_14	= (a100 - a000) / (m1 - m0);
	double diff_d2_14a	= (a010 - a000) / (p1 - p0);
	double diff_d2_14b	= (a110 - a100) / (p1 - p0);
	double diff_d2_14	= diff_d2_14a - 0.2 * (diff_d2_14b - diff_d2_14a);
	double diff_d3_14a	= (a001 - a000) / (d1 - d0);
	double diff_d3_14b	= (a101 - a100) / (d1 - d0);
	double diff_d3_14	= diff_d3_14a - 0.2 * (diff_d3_14b - diff_d3_14a);
	double xder14	= diff_d1_14 * dx3_dt + diff_d2_14 * dx2_dt + diff_d3_14 * dx1_dt;
	check("3v f_table3Veq value    ", val14, res14,	1e-12);
	check("3v f_table3Veq d_dt     ", der14, xder14,	1e-10);
    val14 = P2.value(x3, x2, x1);
	der14 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val14, res14,1e-12);
	check("3v f_table3Veq d_dt     ", der14, xder14,	1e-9);

	x3 = 1340.0;	// 40% more than (x1)
	x2 = -23.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val15 = P1.value(x3, x2, x1);
	double der15 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res15		= a000 + 1.4 * (a100 - a000);
	double diff_d1_15	= (a100 - a000) / (m1 - m0);
	double diff_d2_15a	= (a010 - a000) / (p1 - p0);
	double diff_d2_15b	= (a110 - a100) / (p1 - p0);
	double diff_d2_15	= diff_d2_15a + 1.4 * (diff_d2_15b - diff_d2_15a);
	double diff_d3_15a	= (a001 - a000) / (d1 - d0);
	double diff_d3_15b	= (a101 - a100) / (d1 - d0);
	double diff_d3_15	= diff_d3_15a + 1.4 * (diff_d3_15b - diff_d3_15a);
	double xder15	= diff_d1_15 * dx3_dt + diff_d2_15 * dx2_dt + diff_d3_15 * dx1_dt;
	check("3v f_table3Veq value    ", val15, res15,	1e-12);
	check("3v f_table3Veq d_dt     ", der15, xder15,	1e-10);
    val15 = P2.value(x3, x2, x1);
	der15 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val15, res15,1e-12);
	check("3v f_table3Veq d_dt     ", der15, xder15,	1e-8);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -13.6;	// (y1)
	x1 = 232.1;	// (z3)
    double val16 = P1.value(x3, x2, x1);
	double der16 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res16	= a013 - 0.2 * (a113 - a013);
	double diff_d1_16	= (a113 - a013) / (m1 - m0);
	double diff_d2_16a	= (a023 - a013) / (p2 - p1);
	double diff_d2_16b	= (a123 - a113) / (p2 - p1);
	double diff_d2_16	= diff_d2_16a - 0.2 * (diff_d2_16b - diff_d2_16a);
	double diff_d3_16a	= (a013 - a012) / (d3 - d2);
	double diff_d3_16b	= (a113 - a112) / (d3 - d2);
	double diff_d3_16	= diff_d3_16a - 0.2 * (diff_d3_16b - diff_d3_16a);
	double xder16	= diff_d1_16 * dx3_dt + diff_d2_16 * dx2_dt + diff_d3_16 * dx1_dt;
	check("3v f_table3Veq value    ", val16, res16, 1e-11);
	check("3v f_table3Veq d_dt     ", der16, xder16,	1e-10);
    val16 = P2.value(x3, x2, x1);
	der16 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val16, res16,1e-12);
	check("3v f_table3Veq d_dt     ", der16, xder16,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 132.1;	// (z1)
    double val17 = P1.value(x3, x2, x1);
	double der17 = P1.d_dt( x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res17		= a021 - 0.2 * (a121 - a021);
	double diff_d1_17	= (a121 - a021) / (m1 - m0);
	double diff_d2_17a	= (a021 - a011) / (p2 - p1);
	double diff_d2_17b	= (a121 - a111) / (p2 - p1);
	double diff_d2_17	= diff_d2_17a - 0.2 * (diff_d2_17b - diff_d2_17a);
	double diff_d3_17a	= (a022 - a021) / (d2 - d1);
	double diff_d3_17b	= (a122 - a121) / (d2 - d1);
	double diff_d3_17	= diff_d3_17a - 0.2 * (diff_d3_17b - diff_d3_17a);
	double xder17	= diff_d1_17 * dx3_dt + diff_d2_17 * dx2_dt + diff_d3_17 * dx1_dt;
	check("3v f_table3Veq value    " , val17, res17,	1e-11);
	check("3v f_table3Veq d_dt     ", der17, xder17,	1e-10);
    val17 = P2.value(x3, x2, x1);
	der17 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val17, res17,1e-12);
	check("3v f_table3Veq d_dt     ", der17, xder17,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 232.1;	// (z3)
    double val18 = P1.value(x3, x2, x1);
	double der18 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res18		= a023 - 0.2 * (a123 - a023);
	double diff_d1_18	= (a123 - a023) / (m1 - m0);
	double diff_d2_18a	= (a023 - a013) / (p2 - p1);
	double diff_d2_18b	= (a123 - a113) / (p2 - p1);
	double diff_d2_18	= diff_d2_18a - 0.2 * (diff_d2_18b - diff_d2_18a);
	double diff_d3_18a	= (a023 - a022) / (d3 - d2);
	double diff_d3_18b	= (a123 - a122) / (d3 - d2);
	double diff_d3_18	= diff_d3_18a - 0.2 * (diff_d3_18b - diff_d3_18a);
	double xder18	= diff_d1_18 * dx3_dt +  diff_d2_18 * dx2_dt + diff_d3_18 * dx1_dt;
	check("3v f_table3Veq value    ", val18, res18,	1e-11);
	check("3v f_table3Veq d_dt     ", der18, xder18,	1e-10);
    val18 = P2.value(x3, x2, x1);
	der18 = P2.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	check("3v f_table3Veq value    ", val18, res18,1e-12);
	check("3v f_table3Veq d_dt     ", der18, xder18,	1e-9);
} // closes test3v_f_table3Veq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test3v_f_tabular3V() {
    double dx1_dt = 5.0;
    double dx2_dt = 3.0;
    double dx3_dt = 25.0;

	double m0	= 1200.0;
	double m1	= 1300.0;
	vec1 points1(2);
	points1.set(0, m0);
	points1.set(1, m1);

	double p0	= -33.6;
	double p1	= -13.6;
	double p2	= -3.6;
	vec1 points20(3), points21(3);
	points20.set(0, p0);
	points20.set(1, p1);
	points20.set(2, p2);
	points21.set(0, p0);
	points21.set(1, p1);
	points21.set(2, p2);

	double d0	=  82.1;
	double d1	= 132.1;
	double d2	= 232.1;
	double d3	= 332.1;
	vec1 points30(4), points31(4);
	points30.set(0, d0);
	points30.set(1, d1);
	points30.set(2, d2);
	points30.set(3, d3);
	points31.set(0, d0);
	points31.set(1, d1);
	points31.set(2, d2);
	points31.set(3, d3);

	double a000	= 10;	double a001	= 12;	double a002	= 15;	double a003 = 52;
	double a010	= 20;	double a011	= 25;	double a012	= 7;	double a013 = 27;
	double a020	= 40;	double a021	= 21;	double a022	= 35;	double a023 = 36;
	double a100	= -100;	double a101	= -120;	double a102	= -150;	double a103 = -520;
	double a110	= -200;	double a111	= -250;	double a112	= -70;	double a113 = -270;
	double a120	= -400;	double a121	= -210;	double a122	= -350;	double a123 = -360;
	vec2 Values0(4, 3), Values1(4, 3);
	Values0.set(0, 0, a000);
	Values0.set(1, 0, a001);
	Values0.set(2, 0, a002);
	Values0.set(3, 0, a003);
	Values0.set(0, 1, a010);
	Values0.set(1, 1, a011);
	Values0.set(2, 1, a012);
	Values0.set(3, 1, a013);
	Values0.set(0, 2, a020);
	Values0.set(1, 2, a021);
	Values0.set(2, 2, a022);
	Values0.set(3, 2, a023);
	Values1.set(0, 0, a100);
	Values1.set(1, 0, a101);
	Values1.set(2, 0, a102);
	Values1.set(3, 0, a103);
	Values1.set(0, 1, a110);
	Values1.set(1, 1, a111);
	Values1.set(2, 1, a112);
	Values1.set(3, 1, a113);
	Values1.set(0, 2, a120);
	Values1.set(1, 2, a121);
	Values1.set(2, 2, a122);
	Values1.set(3, 2, a123);

	f_table2V h0(points20, points30, Values0, math::logic::lagrange_first);
	f_table2V h1(points21, points31, Values1, math::logic::lagrange_first);
	std::vector<f_table2V*> tables(2);
	tables[0] = &h0;
	tables[1] = &h1;
	f_tabular3V P1(points1, tables, math::logic::lagrange_first);


	double x3 = 1200.0;	// (x0)
	double x2 = -33.6;	// (y0)
	double x1 = 82.1;	// (z0)
    double val1 = P1.value(x3, x2, x1);
	double der1 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res1 = a000;
	double diff_d1_1	= (a100 - a000) / (m1 - m0);
	double diff_d2_1	= (a010 - a000) / (p1 - p0);
	double diff_d3_1	= (a001 - a000) / (d1 - d0);
	double xder1	= diff_d1_1 * dx3_dt + 	diff_d2_1 * dx2_dt + diff_d3_1 * dx1_dt;
	check("3v f_tabular3V value      ", val1, res1,	1e-12);
	check("3v f_tabular3V d_dt       ", der1, xder1,	1e-10);

	x3 = 1200.0;	// (x0)
	x2 = -13.6;	// (y1)
	x1 = 82.1;	// (z0)
    double val2 = P1.value(x3, x2, x1);
	double der2 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res2 = a010;
	double diff_d1_2	= (a110 - a010) / (m1 - m0);
	double diff_d2_2	= (a020 - a010) / (p2 - p1);
	double diff_d3_2	= (a011 - a010) / (d1 - d0);
	double xder2	= diff_d1_2 * dx3_dt + diff_d2_2 * dx2_dt + diff_d3_2 * dx1_dt;
	check("3v f_tabular3V value      ", val2, res2,	1e-12);
	check("3v f_tabular3V d_dt       ", der2, xder2,	1e-9);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 82.1;			// (z0)
    double val3 = P1.value(x3, x2, x1);
	double der3 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res3 = a010 + 0.25 * (a020 - a010);
	double diff_d1_3a	= (a110 - a010) / (m1 - m0);
	double diff_d1_3b	= (a120 - a020) / (m1 - m0);
	double diff_d1_3	= diff_d1_3a + 0.25 * (diff_d1_3b - diff_d1_3a);
	double diff_d2_3	= (a020 - a010) / (p2 - p1);
	double diff_d3_3a	= (a011 - a010) / (d1 - d0);
	double diff_d3_3b	= (a021 - a020) / (d1 - d0);
	double diff_d3_3	= diff_d3_3a + 0.25 * (diff_d3_3b - diff_d3_3a);
	double xder3	= diff_d1_3 * dx3_dt + diff_d2_3 * dx2_dt + diff_d3_3 * dx1_dt;
	check("3v f_tabular3V value      ", val3, res3,	1e-12);
	check("3v f_tabular3V d_dt       ", der3, xder3,	1e-10);

	x3 = 1200.0;		// (x0)
	x2 = -13.6;	// (y1)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val4 = P1.value(x3, x2, x1);
	double der4 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res4 = a010 + 0.2 * (a011 - a010);
	double diff_d1_4a	= (a110 - a010) / (m1 - m0);
	double diff_d1_4b	= (a111 - a011) / (m1 - m0);
	double diff_d1_4	= diff_d1_4a + 0.2 * (diff_d1_4b - diff_d1_4a);
	double diff_d2_4a	= (a020 - a010) / (p2 - p1);
	double diff_d2_4b	= (a021 - a011) / (p2 - p1);
	double diff_d2_4	= diff_d2_4a + 0.2 * (diff_d2_4b - diff_d2_4a);
	double diff_d3_4	= (a011 - a010) / (d1 - d0);
	double xder4	= diff_d1_4 * dx3_dt + diff_d2_4 * dx2_dt + diff_d3_4 * dx1_dt;
	check("3v f_tabular3V value      ", val4, res4,	1e-12);
	check("3v f_tabular3V d_dt       ", der4, xder4,	1e-10);

	x3 = 1200.0;		// (x0)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val5 = P1.value(x3, x2, x1);
	double der5 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double temp1	= a010 + 0.2  * (a011 - a010);
	double temp2	= a020 + 0.2  * (a021 - a020);
	double res5		= temp1 + 0.25 * (temp2 - temp1);
	double diff_d1_5a	= (a110 - a010) / (m1 - m0);
	double diff_d1_5b	= (a111 - a011) / (m1 - m0);
	double diff_d1_5c	= (a120 - a020) / (m1 - m0);
	double diff_d1_5d	= (a121 - a021) / (m1 - m0);
	double diff_d1_5ab	= diff_d1_5a + 0.2 * (diff_d1_5b - diff_d1_5a);
	double diff_d1_5cd	= diff_d1_5c + 0.2 * (diff_d1_5d - diff_d1_5c);
	double diff_d1_5	= diff_d1_5ab + 0.25 * (diff_d1_5cd - diff_d1_5ab);
	double diff_d2_5a	= (a020 - a010) / (p2 - p1);
	double diff_d2_5b	= (a021 - a011) / (p2 - p1);
	double diff_d2_5	= diff_d2_5a + 0.2 * (diff_d2_5b - diff_d2_5a);
	double diff_d3_5a	= (a011 - a010) / (d1 - d0);
	double diff_d3_5b	= (a021 - a020) / (d1 - d0);
	double diff_d3_5	= diff_d3_5a + 0.25 * (diff_d3_5b - diff_d3_5a);
	double xder5	= diff_d1_5 * dx3_dt + diff_d2_5 * dx2_dt + diff_d3_5 * dx1_dt;
	check("3v f_tabular3V value      ", val5, res5,	1e-12);
	check("3v f_tabular3V d_dt       ", der5, xder5,	1e-10);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val6 = P1.value(x3, x2, x1);
	double der6 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x6temp1a		= a010 + 0.2 * (a011 - a010);
	double x6temp1b		= a020 + 0.2 * (a021 - a020);
	double x6temp1c		= a110 + 0.2 * (a111 - a110);
	double x6temp1d		= a120 + 0.2 * (a121 - a120);
	double x6temp1ab	= x6temp1a + 0.25 * (x6temp1b - x6temp1a);
	double x6temp1cd	= x6temp1c + 0.25 * (x6temp1d - x6temp1c);
	double res6			= x6temp1ab + 0.5 * (x6temp1cd - x6temp1ab);
	double diff_d1_6a	= (a110 - a010) / (m1 - m0);
	double diff_d1_6b	= (a111 - a011) / (m1 - m0);
	double diff_d1_6c	= (a120 - a020) / (m1 - m0);
	double diff_d1_6d	= (a121 - a021) / (m1 - m0);
	double diff_d1_6ab	= diff_d1_6a + 0.2 * (diff_d1_6b - diff_d1_6a);
	double diff_d1_6cd	= diff_d1_6c + 0.2 * (diff_d1_6d - diff_d1_6c);
	double diff_d1_6	= diff_d1_6ab + 0.25 * (diff_d1_6cd - diff_d1_6ab);
	double diff_d2_6a	= (a020 - a010) / (p2 - p1);
	double diff_d2_6b	= (a021 - a011) / (p2 - p1);
	double diff_d2_6c	= (a120 - a110) / (p2 - p1);
	double diff_d2_6d	= (a121 - a111) / (p2 - p1);
	double diff_d2_6ac	= diff_d2_6a + 0.5 * (diff_d2_6c - diff_d2_6a);
	double diff_d2_6bd	= diff_d2_6b + 0.5 * (diff_d2_6d - diff_d2_6b);
	double diff_d2_6	= diff_d2_6ac + 0.2 * (diff_d2_6bd - diff_d2_6ac);
	double diff_d3_6a	= (a011 - a010) / (d1 - d0);
	double diff_d3_6b	= (a021 - a020) / (d1 - d0);
	double diff_d3_6c	= (a111 - a110) / (d1 - d0);
	double diff_d3_6d	= (a121 - a120) / (d1 - d0);
	double diff_d3_6ac	= diff_d3_6a + 0.5 * (diff_d3_6c - diff_d3_6a);
	double diff_d3_6bd	= diff_d3_6b + 0.5 * (diff_d3_6d - diff_d3_6b);
	double diff_d3_6	= diff_d3_6ac + 0.25 * (diff_d3_6bd - diff_d3_6ac);
	double xder6	= diff_d1_6 * dx3_dt + diff_d2_6 * dx2_dt + diff_d3_6 * dx1_dt;
	check("3v f_tabular3V value      ", val6, res6,	1e-12);
	check("3v f_tabular3V d_dt       ", der6, xder6,	1e-09);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 332.1;		// (z3)
    double val7 = P1.value(x3, x2, x1);
	double der7 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x7temp1ab	= a013 + 0.25 * (a023 - a013);
	double x7temp1cd	= a113 + 0.25 * (a123 - a113);
	double res7			= x7temp1ab + 0.5 * (x7temp1cd - x7temp1ab);
	double diff_d1_7a	= (a113 - a013) / (m1 - m0);
	double diff_d1_7b	= (a123 - a023) / (m1 - m0);
	double diff_d1_7	= diff_d1_7a + 0.25 * (diff_d1_7b - diff_d1_7a);
	double diff_d2_7a	= (a023 - a013) / (p2 - p1);
	double diff_d2_7b	= (a123 - a113) / (p2 - p1);
	double diff_d2_7	= diff_d2_7a + 0.5 * (diff_d2_7b - diff_d2_7a);
	double diff_d3_7a	= (a013 - a012) / (d3 - d2);
	double diff_d3_7b	= (a023 - a022) / (d3 - d2);
	double diff_d3_7c	= (a113 - a112) / (d3 - d2);
	double diff_d3_7d	= (a123 - a122) / (d3 - d2);
	double diff_d3_7ac	= diff_d3_7a + 0.5 * (diff_d3_7c - diff_d3_7a);
	double diff_d3_7bd	= diff_d3_7b + 0.5 * (diff_d3_7d - diff_d3_7b);
	double diff_d3_7	= diff_d3_7ac + 0.25 * (diff_d3_7bd - diff_d3_7ac);
	double xder7	= diff_d1_7 * dx3_dt + diff_d2_7 * dx2_dt + diff_d3_7 * dx1_dt;
	check("3v f_tabular3V value      ", val7, res7,	1e-11);
	check("3v f_tabular3V d_dt       ", der7, xder7,	1e-9);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val8 = P1.value(x3, x2, x1);
	double der8 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x8temp1b		= a020 + 0.2 * (a021 - a020);
	double x8temp1d		= a120 + 0.2 * (a121 - a120);
	double res8			= x8temp1b + 0.5 * (x8temp1d - x8temp1b);
	double diff_d1_8a	= (a120 - a020) / (m1 - m0);
	double diff_d1_8b	= (a121 - a021) / (m1 - m0);
	double diff_d1_8	= diff_d1_8a + 0.2 * (diff_d1_8b - diff_d1_8a);
	double diff_d2_8a	= (a020 - a010) / (p2 - p1);
	double diff_d2_8b	= (a021 - a011) / (p2 - p1);
	double diff_d2_8c	= (a120 - a110) / (p2 - p1);
	double diff_d2_8d	= (a121 - a111) / (p2 - p1);
	double diff_d2_8ac	= diff_d2_8a + 0.5 * (diff_d2_8c - diff_d2_8a);
	double diff_d2_8bd	= diff_d2_8b + 0.5 * (diff_d2_8d - diff_d2_8b);
	double diff_d2_8	= diff_d2_8ac + 0.2 * (diff_d2_8bd - diff_d2_8ac);
	double diff_d3_8a	= (a021 - a020) / (d1 - d0);
	double diff_d3_8b	= (a121 - a120) / (d1 - d0);
	double diff_d3_8	= diff_d3_8a + 0.5 * (diff_d3_8b - diff_d3_8a);
	double xder8	= diff_d1_8 * dx3_dt + diff_d2_8 * dx2_dt + diff_d3_8 * dx1_dt;
	check("3v f_tabular3V value      ", val8, res8,	1e-12);
	check("3v f_tabular3V d_dt       ", der8, xder8,	1e-9);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val9 = P1.value(x3, x2, x1);
	double der9 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double x9temp1	= a110 + 0.2  * (a111 - a110);
	double x9temp2	= a120 + 0.2  * (a121 - a120);
	double res9		= x9temp1 + 0.25 * (x9temp2 - x9temp1);
	double diff_d1_9a	= (a110 - a010) / (m1 - m0);
	double diff_d1_9b	= (a111 - a011) / (m1 - m0);
	double diff_d1_9c	= (a120 - a020) / (m1 - m0);
	double diff_d1_9d	= (a121 - a021) / (m1 - m0);
	double diff_d1_9ab	= diff_d1_9a + 0.2 * (diff_d1_9b - diff_d1_9a);
	double diff_d1_9cd	= diff_d1_9c + 0.2 * (diff_d1_9d - diff_d1_9c);
	double diff_d1_9	= diff_d1_9ab + 0.25 * (diff_d1_9cd - diff_d1_9ab);
	double diff_d2_9a	= (a120 - a110) / (p2 - p1);
	double diff_d2_9b	= (a121 - a111) / (p2 - p1);
	double diff_d2_9	= diff_d2_9a + 0.2 * (diff_d2_9b - diff_d2_9a);
	double diff_d3_9a	= (a111 - a110) / (d1 - d0);
	double diff_d3_9b	= (a121 - a120) / (d1 - d0);
	double diff_d3_9	= diff_d3_9a + 0.25 * (diff_d3_9b - diff_d3_9a);
	double xder9	= diff_d1_9 * dx3_dt + diff_d2_9 * dx2_dt + diff_d3_9 * dx1_dt;
	check("3v f_tabular3V value      ", val9, res9,	1e-12);
	check("3v f_tabular3V d_dt       ", der9, xder9,	1e-09);

	x3 = 1250.0;		// 50% between (x0) and (x1)
	x2 = -3.6;		// (y2)
	x1 = 332.1;		// (z3)
    double val10 = P1.value(x3, x2, x1);
	double der10 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res10		= a023 + 0.5 * (a123 - a023);
	double diff_d1_10	= (a123 - a023) / (m1 - m0);
	double diff_d2_10a	= (a023 - a013) / (p2 - p1);
	double diff_d2_10b	= (a123 - a113) / (p2 - p1);
	double diff_d2_10	= diff_d2_10a + 0.5 * (diff_d2_10b - diff_d2_10a);
	double diff_d3_10a	= (a023 - a022) / (d3 - d2);
	double diff_d3_10b	= (a123 - a122) / (d3 - d2);
	double diff_d3_10	= diff_d3_10a + 0.5 * (diff_d3_10b - diff_d3_10a);
	double xder10	= diff_d1_10 * dx3_dt + 	diff_d2_10 * dx2_dt + diff_d3_10 * dx1_dt;
	check("3v f_tabular3V value      ", val10, res10,	1e-12);
	check("3v f_tabular3V d_dt       ", der10, xder10,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -11.1;	// 25% between (y1) and (y2)
	x1 = 332.1;		// (z3)
    double val11 = P1.value(x3, x2, x1);
	double der11 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res11		= a113 + 0.25 * (a123 - a113);
	double diff_d1_11a	= (a113 - a013) / (m1 - m0);
	double diff_d1_11b	= (a123 - a023) / (m1 - m0);
	double diff_d1_11	= diff_d1_11a + 0.25 * (diff_d1_11b - diff_d1_11a);
	double diff_d2_11	= (a123 - a113) / (p2 - p1);
	double diff_d3_11a	= (a113 - a112) / (d3 - d2);
	double diff_d3_11b	= (a123 - a122) / (d3 - d2);
	double diff_d3_11	= diff_d3_11a + 0.25 * (diff_d3_11b - diff_d3_11a);
	double xder11	= diff_d1_11 * dx3_dt + diff_d2_11 * dx2_dt + diff_d3_11 * dx1_dt;
	check("3v f_tabular3V value      ", val11, res11,	1e-12);
	check("3v f_tabular3V d_dt       ", der11, xder11,	1e-8);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 92.1;			// 20% between (z0) and (z1)
    double val12 = P1.value(x3, x2, x1);
	double der12 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res12		= a120 + 0.2 * (a121 - a120);
	double diff_d1_12a	= (a120 - a020) / (m1 - m0);
	double diff_d1_12b	= (a121 - a021) / (m1 - m0);
	double diff_d1_12	= diff_d1_12a + 0.2 * (diff_d1_12b - diff_d1_12a);
	double diff_d2_12a	= (a120 - a110) / (p2 - p1);
	double diff_d2_12b	= (a121 - a111) / (p2 - p1);
	double diff_d2_12	= diff_d2_12a + 0.2 * (diff_d2_12b - diff_d2_12a);
	double diff_d3_12   = (a121 - a120) / (d1 - d0);
	double xder12	= diff_d1_12 * dx3_dt + diff_d2_12 * dx2_dt + diff_d3_12 * dx1_dt;
	check("3v f_tabular3V value      ", val12, res12,	1e-12);
	check("3v f_tabular3V d_dt       ", der12, xder12,	1e-08);

	x3 = 1300.0;		// (x1)
	x2 = -3.6;		// (y2)
	x1 = 332.1;		// (z3)
    double val13 = P1.value(x3, x2, x1);
	double der13 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res13		= a123;
	double diff_d1_13	= (a123 - a023) / (m1 - m0);
	double diff_d2_13	= (a123 - a113) / (p2 - p1);
	double diff_d3_13   = (a123 - a122) / (d3 - d2);
	double xder13	= diff_d1_13 * dx3_dt + diff_d2_13 * dx2_dt + diff_d3_13 * dx1_dt;
	check("3v f_tabular3V value      ", val13, res13,	1e-12);
	check("3v f_tabular3V d_dt       ", der13, xder13,	1e-9);
	
	x3 = 1180.0;	// 20% less than (x0)
	x2 = -33.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val14 = P1.value(x3, x2, x1);
	double der14 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res14		= a000 - 0.2 * (a100 - a000);
	double diff_d1_14	= (a100 - a000) / (m1 - m0);
	double diff_d2_14a	= (a010 - a000) / (p1 - p0);
	double diff_d2_14b	= (a110 - a100) / (p1 - p0);
	double diff_d2_14	= diff_d2_14a - 0.2 * (diff_d2_14b - diff_d2_14a);
	double diff_d3_14a	= (a001 - a000) / (d1 - d0);
	double diff_d3_14b	= (a101 - a100) / (d1 - d0);
	double diff_d3_14	= diff_d3_14a - 0.2 * (diff_d3_14b - diff_d3_14a);
	double xder14	= diff_d1_14 * dx3_dt + diff_d2_14 * dx2_dt + diff_d3_14 * dx1_dt;
	check("3v f_tabular3V value      ", val14, res14,	1e-12);
	check("3v f_tabular3V d_dt       ", der14, xder14,	1e-9);

	x3 = 1340.0;	// 40% more than (x1)
	x2 = -33.6;	// (y0)
	x1 = 82.1;	// (z0)
    double val15 = P1.value(x3, x2, x1);
	double der15 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res15		= a000 + 1.4 * (a100 - a000);
	double diff_d1_15	= (a100 - a000) / (m1 - m0);
	double diff_d2_15a	= (a010 - a000) / (p1 - p0);
	double diff_d2_15b	= (a110 - a100) / (p1 - p0);
	double diff_d2_15	= diff_d2_15a + 1.4 * (diff_d2_15b - diff_d2_15a);
	double diff_d3_15a	= (a001 - a000) / (d1 - d0);
	double diff_d3_15b	= (a101 - a100) / (d1 - d0);
	double diff_d3_15	= diff_d3_15a + 1.4 * (diff_d3_15b - diff_d3_15a);
	double xder15	= diff_d1_15 * dx3_dt + diff_d2_15 * dx2_dt + diff_d3_15 * dx1_dt;
	check("3v f_tabular3V value      ", val15, res15,	1e-12);
	check("3v f_tabular3V d_dt       ", der15, xder15,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -13.6;	// (y1)
	x1 = 332.1;	// (z3)
    double val16 = P1.value(x3, x2, x1);
	double der16 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res16		= a013 - 0.2 * (a113 - a013);
	double diff_d1_16	= (a113 - a013) / (m1 - m0);
	double diff_d2_16a	= (a023 - a013) / (p2 - p1);
	double diff_d2_16b	= (a123 - a113) / (p2 - p1);
	double diff_d2_16	= diff_d2_16a - 0.2 * (diff_d2_16b - diff_d2_16a);
	double diff_d3_16a	= (a013 - a012) / (d3 - d2);
	double diff_d3_16b	= (a113 - a112) / (d3 - d2);
	double diff_d3_16	= diff_d3_16a - 0.2 * (diff_d3_16b - diff_d3_16a);
	double xder16	= diff_d1_16 * dx3_dt + diff_d2_16 * dx2_dt + diff_d3_16 * dx1_dt;
	check("3v f_tabular3V value      ", val16, res16,	1e-11);
	check("3v f_tabular3V d_dt       ", der16, xder16,	1e-10);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 132.1;	// (z1)
    double val17 = P1.value(x3, x2, x1);
	double der17 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res17		= a021 - 0.2 * (a121 - a021);
	double diff_d1_17	= (a121 - a021) / (m1 - m0);
	double diff_d2_17a	= (a021 - a011) / (p2 - p1);
	double diff_d2_17b	= (a121 - a111) / (p2 - p1);
	double diff_d2_17	= diff_d2_17a - 0.2 * (diff_d2_17b - diff_d2_17a);
	double diff_d3_17a	= (a022 - a021) / (d2 - d1);
	double diff_d3_17b	= (a122 - a121) / (d2 - d1);
	double diff_d3_17	= diff_d3_17a - 0.2 * (diff_d3_17b - diff_d3_17a);
	double xder17	= diff_d1_17 * dx3_dt + diff_d2_17 * dx2_dt + diff_d3_17 * dx1_dt;
	check("3v f_tabular3V value      ", val17, res17,	1e-11);
	check("3v f_tabular3V d_dt       ", der17, xder17,	1e-9);

	x3 = 1180.0;		// 20% less than (x0)
	x2 = -3.6;	// (y2)
	x1 = 332.1;	// (z3)
    double val18 = P1.value(x3, x2, x1);
	double der18 = P1.d_dt(x3, x2, x1, dx3_dt, dx2_dt, dx1_dt);
	double res18		= a023 - 0.2 * (a123 - a023);
	double diff_d1_18	= (a123 - a023) / (m1 - m0);
	double diff_d2_18a	= (a023 - a013) / (p2 - p1);
	double diff_d2_18b	= (a123 - a113) / (p2 - p1);
	double diff_d2_18	= diff_d2_18a - 0.2 * (diff_d2_18b - diff_d2_18a);
	double diff_d3_18a	= (a023 - a022) / (d3 - d2);
	double diff_d3_18b	= (a123 - a122) / (d3 - d2);
	double diff_d3_18	= diff_d3_18a - 0.2 * (diff_d3_18b - diff_d3_18a);
	double xder18	= diff_d1_18 * dx3_dt + diff_d2_18 * dx2_dt + diff_d3_18 * dx1_dt;
	check("3v f_tabular3V value      ", val18, res18,	1e-11);
	check("3v f_tabular3V d_dt       ", der18, xder18,	1e-10);
} // closes test3v_f_tabular3V

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test4v_f_table4Veq() {
    double dx4_dt = 25.0;
    double dx3_dt = 3.0;
    double dx2_dt = 5.0;
    double dx1_dt = 4.0;

    vec1 points1(2);
    double m0	= 1200.0;	points1.set(0, m0);
    double m1	= 1300.0;	points1.set(1, m1);

    vec1 points2(2);
    double r0	= 25.4;		points2.set(0, r0);
    double r1	= 35.4;		points2.set(1, r1);

    vec1 points3(3);
    double d0	=  82.1;	points3.set(0, d0);
    double d1	= 132.1;	points3.set(1, d1);
    double d2	= 182.1;	points3.set(2, d2);

    vec1 points4(2);
    double t0	=  3.0;		points4.set(0, t0);
    double t1	=  7.0;		points4.set(1, t1);

    double a0000 = 10, a0001 = 12, a0010 = 15, a0011 = 18, a0020 = 25, a0021 = 29;
    double a0100 = 60, a0101 = 56, a0110 = 65, a0111 = 58, a0120 = 85, a0121 = 94;
    double a1000 =210, a1001 =216, a1010 =215, a1011 =211, a1020 =230, a1021 =245;
    double a1100 =260, a1101 =215, a1110 =265, a1111 =287, a1120 =294, a1121 =384;
    vec4 VVValues(2, 3, 2, 2);
    VVValues.set(0, 0, 0, 0, a0000);	VVValues.set(1, 0, 0, 0, a0001);
    VVValues.set(0, 1, 0, 0, a0010);	VVValues.set(1, 1, 0, 0, a0011);
    VVValues.set(0, 2, 0, 0, a0020);	VVValues.set(1, 2, 0, 0, a0021);
    VVValues.set(0, 0, 1, 0, a0100);	VVValues.set(1, 0, 1, 0, a0101);
    VVValues.set(0, 1, 1, 0, a0110);	VVValues.set(1, 1, 1, 0, a0111);
    VVValues.set(0, 2, 1, 0, a0120);	VVValues.set(1, 2, 1, 0, a0121);
    VVValues.set(0, 0, 0, 1, a1000);	VVValues.set(1, 0, 0, 1, a1001);
    VVValues.set(0, 1, 0, 1, a1010);	VVValues.set(1, 1, 0, 1, a1011);
    VVValues.set(0, 2, 0, 1, a1020);	VVValues.set(1, 2, 0, 1, a1021);
    VVValues.set(0, 0, 1, 1, a1100);	VVValues.set(1, 0, 1, 1, a1101);
    VVValues.set(0, 1, 1, 1, a1110);	VVValues.set(1, 1, 1, 1, a1111);
    VVValues.set(0, 2, 1, 1, a1120);	VVValues.set(1, 2, 1, 1, a1121);

    f_table4V P1(points1, points2, points3, points4, VVValues, math::logic::lagrange_first_precompute);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi1(), true);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi2(), true);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi4(), true);
    f_table4V P2(points1, points2, points3, points4, VVValues, math::logic::lagrange_first);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi1(), true);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi2(), true);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi4(), true);
    f_table4V P6(points1, points2, points3, points4, VVValues, math::logic::hermite_first);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi1(), true);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi2(), true);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi4(), true);


    double x4 = 1200.0;		// (x0)
    double x3 =   25.4;	// (y0)
    double x2 =   82.1;		// (z0)
    double x1 =    3.0;	// (q0)
    double val1 = P1.value(x4, x3, x2, x1);
    double der1 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res1 = a0000;
    double diff_d1_1	= (a1000 - a0000) / (m1 - m0);
    double diff_d2_1	= (a0100 - a0000) / (r1 - r0);
    double diff_d3_1	= (a0010 - a0000) / (d1 - d0);
    double diff_d4_1	= (a0001 - a0000) / (t1 - t0);
    double xder1	= diff_d1_1 * dx4_dt + diff_d2_1 * dx3_dt +
                    diff_d3_1 * dx2_dt + diff_d4_1 * dx1_dt;
    check("4v f_table4Veq value    ", val1, res1,	1e-15);
    check("4v f_table4Veq d_dt     ", der1,	xder1,	1e-11);
    val1 = P2.value(x4, x3, x2, x1);
    der1 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val1, res1,1e-13);
    check("4v f_table4Veq d_dt     ", der1,	xder1,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =   82.1;		// (z0)
    x1 =    3.0;	// (q0)
    double val2 = P1.value(x4, x3, x2, x1);
    double der2 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res2 = a0100;
    double diff_d1_2	= (a1100 - a0100) / (m1 - m0);
    double diff_d2_2	= (a0100 - a0000) / (r1 - r0);
    double diff_d3_2	= (a0110 - a0100) / (d1 - d0);
    double diff_d4_2	= (a0101 - a0100) / (t1 - t0);
    double xder2	= diff_d1_2 * dx4_dt + diff_d2_2 * dx3_dt +
                    diff_d3_2 * dx2_dt + diff_d4_2 * dx1_dt;
    check("4v f_table4Veq value    ", val2, res2,	1e-14);
    check("4v f_table4Veq d_dt     ", der2,	xder2,	1e-11);
    val2 = P2.value(x4, x3, x2, x1);
    der2 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val2, res2,1e-13);
    check("4v f_table4Veq d_dt     ", der2,	xder2,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    3.0;	// (q0)
    double val3 = P1.value(x4, x3, x2, x1);
    double der3 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res3 = a0020;
    double diff_d1_3	= (a1020 - a0020) / (m1 - m0);
    double diff_d2_3	= (a0120 - a0020) / (r1 - r0);
    double diff_d3_3	= (a0020 - a0010) / (d2 - d1);
    double diff_d4_3	= (a0021 - a0020) / (t1 - t0);
    double xder3	= diff_d1_3 * dx4_dt + diff_d2_3 * dx3_dt +
                    diff_d3_3 * dx2_dt + diff_d4_3 * dx1_dt;
    check("4v f_table4Veq value    ", val3, res3,	1e-14);
    check("4v f_table4Veq d_dt     ", der3,	xder3,	1e-11);
    val3 = P2.value(x4, x3, x2, x1);
    der3 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val3, res3,1e-13);
    check("4v f_table4Veq d_dt     ", der3,	xder3,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  132.1;		// (z1)
    x1 =    7.0;	// (q1)
    double val4 = P1.value(x4, x3, x2, x1);
    double der4 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res4 = a0011;
    double diff_d1_4	= (a1011 - a0011) / (m1 - m0);
    double diff_d2_4	= (a0111 - a0011) / (r1 - r0);
    double diff_d3_4	= (a0021 - a0011) / (d2 - d1);
    double diff_d4_4	= (a0011 - a0010) / (t1 - t0);
    double xder4	= diff_d1_4 * dx4_dt + diff_d2_4 * dx3_dt +
                    diff_d3_4 * dx2_dt + diff_d4_4 * dx1_dt;
    check("4v f_table4Veq value    ", val4, res4,	1e-14);
    check("4v f_table4Veq d_dt     ", der4,	xder4,	1e-11);
    val4 = P2.value(x4, x3, x2, x1);
    der4 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val4, res4,1e-13);
    check("4v f_table4Veq d_dt     ", der4,	xder4,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    3.0;	// (q0)
    double val5 = P1.value(x4, x3, x2, x1);
    double der5 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res5 = a0120;
    double diff_d1_5	= (a1120 - a0120) / (m1 - m0);
    double diff_d2_5	= (a0120 - a0020) / (r1 - r0);
    double diff_d3_5	= (a0120 - a0110) / (d2 - d1);
    double diff_d4_5	= (a0121 - a0120) / (t1 - t0);
    double xder5	= diff_d1_5 * dx4_dt +   diff_d2_5 * dx3_dt +
                    diff_d3_5 * dx2_dt + diff_d4_5 * dx1_dt;
    check("4v f_table4Veq value    ", val5, res5,	1e-14);
    check("4v f_table4Veq d_dt     ", der5,	xder5,	1e-11);
    val5 = P2.value(x4, x3, x2, x1);
    der5 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val5, res5,1e-13);
    check("4v f_table4Veq d_dt     ", der5,	xder5,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  132.1;		// (z1)
    x1 =    7.0;	// (q1)
    double val6 = P1.value(x4, x3, x2, x1);
    double der6 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res6 = a0111;
    double diff_d1_6	= (a1111 - a0111) / (m1 - m0);
    double diff_d2_6	= (a0111 - a0011) / (r1 - r0);
    double diff_d3_6	= (a0121 - a0111) / (d2 - d1);
    double diff_d4_6	= (a0111 - a0110) / (t1 - t0);
    double xder6	= diff_d1_6 * dx4_dt + diff_d2_6 * dx3_dt +
                    diff_d3_6 * dx2_dt +  diff_d4_6 * dx1_dt;
    check("4v f_table4Veq value    ", val6, res6,	1e-14);
    check("4v f_table4Veq d_dt     ", der6,	xder6,	1e-11);
    val6 = P2.value(x4, x3, x2, x1);
    der6 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val6, res6,1e-13);
    check("4v f_table4Veq d_dt     ", der6,	xder6,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val7 = P1.value(x4, x3, x2, x1);
    double der7 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res7 = a0021;
    double diff_d1_7	= (a1021 - a0021) / (m1 - m0);
    double diff_d2_7	= (a0121 - a0021) / (r1 - r0);
    double diff_d3_7	= (a0021 - a0011) / (d2 - d1);
    double diff_d4_7	= (a0021 - a0020) / (t1 - t0);
    double xder7	= diff_d1_7 * dx4_dt +  diff_d2_7 * dx3_dt +
                    diff_d3_7 * dx2_dt +  diff_d4_7 * dx1_dt;
    check("4v f_table4Veq value    ", val7, res7,	1e-14);
    check("4v f_table4Veq d_dt     ", der7,	xder7,	1e-11);
    val7 = P2.value(x4, x3, x2, x1);
    der7 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val7, res7,1e-13);
    check("4v f_table4Veq d_dt     ", der7,	xder7,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val8 = P1.value(x4, x3, x2, x1);
    double der8 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res8 = a0121;
    double diff_d1_8	= (a1121 - a0121) / (m1 - m0);
    double diff_d2_8	= (a0121 - a0021) / (r1 - r0);
    double diff_d3_8	= (a0121 - a0111) / (d2 - d1);
    double diff_d4_8	= (a0121 - a0120) / (t1 - t0);
    double xder8	= diff_d1_8 * dx4_dt +  diff_d2_8 * dx3_dt +
                    diff_d3_8 * dx2_dt +  diff_d4_8 * dx1_dt;
    check("4v f_table4Veq value    ", val8, res8,	1e-14);
    check("4v f_table4Veq d_dt     ", der8,	xder8,	1e-10);
    val8 = P2.value(x4, x3, x2, x1);
    der8 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val8, res8,1e-13);
    check("4v f_table4Veq d_dt     ", der8,	xder8,	1e-9);

    val8 = P6.value(x4, x3, x2, x1);
    check("4v f_table4Veq value    ", val8, res8,	1e-12);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =   82.1;		// (z0)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val9 = P1.value(x4, x3, x2, x1);
    double der9 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res9 = a0000 + 0.25 * (a0001 - a0000);
    double diff_d1_9_a	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_9_b	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_9	= diff_d1_9_a + 0.25 * (diff_d1_9_b - diff_d1_9_a);
    double diff_d2_9_a	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_9_b	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_9	= diff_d2_9_a + 0.25 * (diff_d2_9_b - diff_d2_9_a);
    double diff_d3_9_a	= (a0010 - a0000) / (d1 - d0);
    double diff_d3_9_b	= (a0011 - a0001) / (d1 - d0);
    double diff_d3_9	= diff_d3_9_a + 0.25 * (diff_d3_9_b - diff_d3_9_a);
    double diff_d4_9	= (a0001 - a0000) / (t1 - t0);
    double xder9	= diff_d1_9 * dx4_dt +  diff_d2_9 * dx3_dt +
                    diff_d3_9 * dx2_dt + diff_d4_9 * dx1_dt;
    check("4v f_table4Veq value    ", val9, res9,	1e-14);
    check("4v f_table4Veq d_dt     ", der9,	xder9,	1e-11);
    val9 = P2.value(x4, x3, x2, x1);
    der9 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val9, res9,1e-13);
    check("4v f_table4Veq d_dt     ", der9,	xder9,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    3.0;	// (q0)
    double val10 = P1.value(x4, x3, x2, x1);
    double der10 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res10 = a0100 + 0.50 * (a0110 - a0100);
    double diff_d1_10_a	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_10_b	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_10	= diff_d1_10_a + 0.50 * (diff_d1_10_b - diff_d1_10_a);
    double diff_d2_10_a	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_10_b	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_10	= diff_d2_10_a + 0.50 * (diff_d2_10_b - diff_d2_10_a);
    double diff_d3_10	= (a0110 - a0100) / (d1 - d0);
    double diff_d4_10_a	= (a0101 - a0100) / (t1 - t0);
    double diff_d4_10_b	= (a0111 - a0110) / (t1 - t0);
    double diff_d4_10	= diff_d4_10_a + 0.50 * (diff_d4_10_b - diff_d4_10_a);
    double xder10	= diff_d1_10 * dx4_dt +  diff_d2_10 * dx3_dt +
                     diff_d3_10 * dx2_dt +  diff_d4_10 * dx1_dt;
    check("4v f_table4Veq value    ", val10, res10,	1e-14);
    check("4v f_table4Veq d_dt     ", der10, xder10,	1e-11);
    val10 = P2.value(x4, x3, x2, x1);
    der10 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val10, res10,1e-13);
    check("4v f_table4Veq d_dt     ", der10, xder10,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val11 = P1.value(x4, x3, x2, x1);
    double der11 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double temp1 = a0100 + 0.25 * (a0101 - a0100);
    double temp2 = a0110 + 0.25 * (a0111 - a0110);
    double res11 = temp1 + 0.50 * (temp2 - temp1);
    double diff_d1_11_aa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_11_ab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_11_a		= diff_d1_11_aa + 0.25 * (diff_d1_11_ab - diff_d1_11_aa);
    double diff_d1_11_ba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_11_bb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_11_b		= diff_d1_11_ba + 0.25 * (diff_d1_11_bb - diff_d1_11_ba);
    double diff_d1_11		= diff_d1_11_a  + 0.50 * (diff_d1_11_b  - diff_d1_11_a);
    double diff_d2_11_aa	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_11_ab	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_11_a		= diff_d2_11_aa + 0.25 * (diff_d2_11_ab - diff_d2_11_aa);
    double diff_d2_11_ba	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_11_bb	= (a0111 - a0011) / (r1 - r0);
    double diff_d2_11_b		= diff_d2_11_ba + 0.25 * (diff_d2_11_bb - diff_d2_11_ba);
    double diff_d2_11		= diff_d2_11_a  + 0.50 * (diff_d2_11_b  - diff_d2_11_a);
    double diff_d3_11_a		= (a0110 - a0100) / (d1 - d0);
    double diff_d3_11_b		= (a0111 - a0101) / (d1 - d0);
    double diff_d3_11		= diff_d3_11_a + 0.25 * (diff_d3_11_b - diff_d3_11_a);
    double diff_d4_11_a		= (a0101 - a0100) / (t1 - t0);
    double diff_d4_11_b		= (a0111 - a0110) / (t1 - t0);
    double diff_d4_11		= diff_d4_11_a + 0.50 * (diff_d4_11_b - diff_d4_11_a);
    double xder11	= diff_d1_11 * dx4_dt + diff_d2_11 * dx3_dt +
                     diff_d3_11 * dx2_dt + diff_d4_11 * dx1_dt;
    check("4v f_table4Veq value    ", val11, res11,	1e-14);
    check("4v f_table4Veq d_dt     ", der11, xder11,	1e-11);
    val11 = P2.value(x4, x3, x2, x1);
    der11 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val11, res11,1e-13);
    check("4v f_table4Veq d_dt     ", der11, xder11,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val12 = P1.value(x4, x3, x2, x1);
    double der12 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    double temp3 = a0100 + 0.25 * (a0101 - a0100);
    double temp4 = a0110 + 0.25 * (a0111 - a0110);
    double temp5 = temp1 + 0.50 * (temp2 - temp1);
    double temp6 = temp3 + 0.50 * (temp4 - temp3);
    double res12 = temp5 + 0.80 * (temp6 - temp5);
    double diff_d1_12_aaa	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_12_aab	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_12_aa	= diff_d1_12_aaa + 0.25 * (diff_d1_12_aab - diff_d1_12_aaa);
    double diff_d1_12_aba	= (a1010 - a0010) / (m1 - m0);
    double diff_d1_12_abb	= (a1011 - a0011) / (m1 - m0);
    double diff_d1_12_ab	= diff_d1_12_aba + 0.25 * (diff_d1_12_abb - diff_d1_12_aba);
    double diff_d1_12_baa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_12_bab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_12_ba	= diff_d1_12_baa + 0.25 * (diff_d1_12_bab - diff_d1_12_baa);
    double diff_d1_12_bba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_12_bbb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_12_bb	= diff_d1_12_bba + 0.25 * (diff_d1_12_bbb - diff_d1_12_bba);
    double diff_d1_12_a		= diff_d1_12_aa + 0.50 * (diff_d1_12_ab - diff_d1_12_aa);
    double diff_d1_12_b		= diff_d1_12_ba + 0.50 * (diff_d1_12_bb - diff_d1_12_ba);
    double diff_d1_12		= diff_d1_12_a  + 0.80 * (diff_d1_12_b - diff_d1_12_a);
    double diff_d2_12_aa	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_12_ab	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_12_a		= diff_d2_12_aa + 0.25 * (diff_d2_12_ab - diff_d2_12_aa);
    double diff_d2_12_ba	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_12_bb	= (a0111 - a0011) / (r1 - r0);
    double diff_d2_12_b		= diff_d2_12_ba + 0.25 * (diff_d2_12_bb - diff_d2_12_ba);
    double diff_d2_12		= diff_d2_12_a  + 0.50 * (diff_d2_12_b  - diff_d2_12_a);
    double diff_d3_12_aa	= (a0010 - a0000) / (d1 - d0);
    double diff_d3_12_ab	= (a0011 - a0001) / (d1 - d0);
    double diff_d3_12_a		= diff_d3_12_aa + 0.25 * (diff_d3_12_ab - diff_d3_12_aa);
    double diff_d3_12_ba	= (a0110 - a0100) / (d1 - d0);
    double diff_d3_12_bb	= (a0111 - a0101) / (d1 - d0);
    double diff_d3_12_b		= diff_d3_12_ba + 0.25 * (diff_d3_12_bb - diff_d3_12_ba);
    double diff_d3_12		= diff_d3_12_a  + 0.80 * (diff_d3_12_b  - diff_d3_12_a);
    double diff_d4_12_aa	= (a0001 - a0000) / (t1 - t0);
    double diff_d4_12_ab	= (a0011 - a0010) / (t1 - t0);
    double diff_d4_12_a		= diff_d4_12_aa + 0.50 * (diff_d4_12_ab - diff_d4_12_aa);
    double diff_d4_12_ba	= (a0101 - a0100) / (t1 - t0);
    double diff_d4_12_bb	= (a0111 - a0110) / (t1 - t0);
    double diff_d4_12_b		= diff_d4_12_ba + 0.50 * (diff_d4_12_bb - diff_d4_12_ba);
    double diff_d4_12		= diff_d4_12_a  + 0.80 * (diff_d4_12_b  - diff_d4_12_a);
    double xder12	= diff_d1_12 * dx4_dt +  diff_d2_12 * dx3_dt +
                     diff_d3_12 * dx2_dt +  diff_d4_12 * dx1_dt;
    check("4v f_table4Veq value    ", val12, res12,	1e-14);
    check("4v f_table4Veq d_dt     ", der12, xder12,	1e-11);
    val12 = P2.value(x4, x3, x2, x1);
    der12 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val12, res12,1e-13);
    check("4v f_table4Veq d_dt     ", der12, xder12,	1e-9);

    x4 = 1210.0;	// 10% between (x0) and (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val13 = P1.value(x4, x3, x2, x1);
    // differential too complicated to compute manually
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    temp3 = a0100 + 0.25 * (a0101 - a0100);
    temp4 = a0110 + 0.25 * (a0111 - a0110);
    temp5 = a1000 + 0.25 * (a1001 - a1000);
    temp6 = a1010 + 0.25 * (a1011 - a1010);
    double temp7 = a1100 + 0.25 * (a1101 - a1100);
    double temp8 = a1110 + 0.25 * (a1111 - a1110);
    double temp9  = temp1 + 0.50 * (temp2 - temp1);
    double temp10 = temp3 + 0.50 * (temp4 - temp3);
    double temp11 = temp5 + 0.50 * (temp6 - temp5);
    double temp12 = temp7 + 0.50 * (temp8 - temp7);
    double temp13 = temp9  + 0.80 * (temp10 - temp9);
    double temp14 = temp11 + 0.80 * (temp12 - temp11);
    double res13 = temp13 + 0.10 * (temp14 - temp13);
    check("4v f_table4Veq value    ", val13, res13,	1e-14);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val14 = P1.value(x4, x3, x2, x1);
    double der14 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a1000 + 0.25 * (a1001 - a1000);
    temp2 = a1010 + 0.25 * (a1011 - a1010);
    temp3 = a1100 + 0.25 * (a1101 - a1100);
    temp4 = a1110 + 0.25 * (a1111 - a1110);
    temp5 = temp1 + 0.50 * (temp2 - temp1);
    temp6 = temp3 + 0.50 * (temp4 - temp3);
    double res14 = temp5 + 0.80 * (temp6 - temp5);
    double diff_d1_14_aaa	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_14_aab	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_14_aa	= diff_d1_14_aaa + 0.25 * (diff_d1_14_aab - diff_d1_14_aaa);
    double diff_d1_14_aba	= (a1010 - a0010) / (m1 - m0);
    double diff_d1_14_abb	= (a1011 - a0011) / (m1 - m0);
    double diff_d1_14_ab	= diff_d1_14_aba + 0.25 * (diff_d1_14_abb - diff_d1_14_aba);
    double diff_d1_14_baa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_14_bab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_14_ba	= diff_d1_14_baa + 0.25 * (diff_d1_14_bab - diff_d1_14_baa);
    double diff_d1_14_bba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_14_bbb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_14_bb	= diff_d1_14_bba + 0.25 * (diff_d1_14_bbb - diff_d1_14_bba);
    double diff_d1_14_a		= diff_d1_14_aa + 0.50 * (diff_d1_14_ab - diff_d1_14_aa);
    double diff_d1_14_b		= diff_d1_14_ba + 0.50 * (diff_d1_14_bb - diff_d1_14_ba);
    double diff_d1_14		= diff_d1_14_a  + 0.80 * (diff_d1_14_b - diff_d1_14_a);
    double diff_d2_14_aa	= (a1100 - a1000) / (r1 - r0);
    double diff_d2_14_ab	= (a1101 - a1001) / (r1 - r0);
    double diff_d2_14_a		= diff_d2_14_aa + 0.25 * (diff_d2_14_ab - diff_d2_14_aa);
    double diff_d2_14_ba	= (a1110 - a1010) / (r1 - r0);
    double diff_d2_14_bb	= (a1111 - a1011) / (r1 - r0);
    double diff_d2_14_b		= diff_d2_14_ba + 0.25 * (diff_d2_14_bb - diff_d2_14_ba);
    double diff_d2_14		= diff_d2_14_a  + 0.50 * (diff_d2_14_b  - diff_d2_14_a);
    double diff_d3_14_aa	= (a1010 - a1000) / (d1 - d0);
    double diff_d3_14_ab	= (a1011 - a1001) / (d1 - d0);
    double diff_d3_14_a		= diff_d3_14_aa + 0.25 * (diff_d3_14_ab - diff_d3_14_aa);
    double diff_d3_14_ba	= (a1110 - a1100) / (d1 - d0);
    double diff_d3_14_bb	= (a1111 - a1101) / (d1 - d0);
    double diff_d3_14_b		= diff_d3_14_ba + 0.25 * (diff_d3_14_bb - diff_d3_14_ba);
    double diff_d3_14		= diff_d3_14_a  + 0.80 * (diff_d3_14_b  - diff_d3_14_a);
    double diff_d4_14_aa	= (a1001 - a1000) / (t1 - t0);
    double diff_d4_14_ab	= (a1011 - a1010) / (t1 - t0);
    double diff_d4_14_a		= diff_d4_14_aa + 0.50 * (diff_d4_14_ab - diff_d4_14_aa);
    double diff_d4_14_ba	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_14_bb	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_14_b		= diff_d4_14_ba + 0.50 * (diff_d4_14_bb - diff_d4_14_ba);
    double diff_d4_14		= diff_d4_14_a  + 0.80 * (diff_d4_14_b  - diff_d4_14_a);
    double xder14	= diff_d1_14 * dx4_dt +  diff_d2_14 * dx3_dt +
                     diff_d3_14 * dx2_dt +   diff_d4_14 * dx1_dt;
    check("4v f_table4Veq value    ", val14, res14,	1e-13);
    check("4v f_table4Veq d_dt     ", der14, xder14,	1e-10);
    val14 = P2.value(x4, x3, x2, x1);
    der14 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val14, res14,1e-13);
    check("4v f_table4Veq d_dt     ", der14, xder14,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    3.0;	// (q0)
    double val15 = P1.value(x4, x3, x2, x1);
    double der15 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res15 = a1100 + 0.50 * (a1110 - a1100);
    double diff_d1_15_a	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_15_b	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_15	= diff_d1_15_a + 0.50 * (diff_d1_15_b - diff_d1_15_a);
    double diff_d2_15_a	= (a1100 - a1000) / (r1 - r0);
    double diff_d2_15_b	= (a1110 - a1010) / (r1 - r0);
    double diff_d2_15	= diff_d2_15_a + 0.50 * (diff_d2_15_b - diff_d2_15_a);
    double diff_d3_15	= (a1110 - a1100) / (d1 - d0);
    double diff_d4_15_a	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_15_b	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_15	= diff_d4_15_a + 0.50 * (diff_d4_15_b - diff_d4_15_a);
    double xder15	= diff_d1_15 * dx4_dt + diff_d2_15 * dx3_dt +
                     diff_d3_15 * dx2_dt +  diff_d4_15 * dx1_dt;
    check("4v f_table4Veq value    ", val15, res15,	1e-13);
    check("4v f_table4Veq d_dt     ", der15, xder15,	1e-10);
    val15 = P2.value(x4, x3, x2, x1);
    der15 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val15, res15,1e-13);
    check("4v f_table4Veq d_dt     ", der15, xder15,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val16 = P1.value(x4, x3, x2, x1);
    double der16= P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res16 = a1020 + 0.25 * (a1021 - a1020);
    double diff_d1_16_a	= (a1020 - a0020) / (m1 - m0);
    double diff_d1_16_b	= (a1021 - a0021) / (m1 - m0);
    double diff_d1_16	= diff_d1_16_a + 0.25 * (diff_d1_16_b - diff_d1_16_a);
    double diff_d2_16_a	= (a1120 - a1020) / (r1 - r0);
    double diff_d2_16_b	= (a1121 - a1021) / (r1 - r0);
    double diff_d2_16	= diff_d2_16_a + 0.25 * (diff_d2_16_b - diff_d2_16_a);
    double diff_d3_16_a	= (a1020 - a1010) / (d1 - d0);
    double diff_d3_16_b	= (a1021 - a1011) / (d1 - d0);
    double diff_d3_16	= diff_d3_16_a + 0.25 * (diff_d3_16_b - diff_d3_16_a);
    double diff_d4_16	= (a1021 - a1020) / (t1 - t0);
    double xder16	= diff_d1_16 * dx4_dt +   diff_d2_16 * dx3_dt +
                     diff_d3_16 * dx2_dt +  diff_d4_16 * dx1_dt;
    check("4v f_table4Veq value    ", val16, res16,	1e-13);
    check("4v f_table4Veq d_dt     ", der16, xder16,	1e-10);
    val16 = P2.value(x4, x3, x2, x1);
    der16 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val16, res16,1e-13);
    check("4v f_table4Veq d_dt     ", der16, xder16,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =   82.1;		// (z0)
    x1 =    7.0;	// (q1)
    double val17 = P1.value(x4, x3, x2, x1);
    double der17 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res17 = a1001 + 0.80 * (a1101 - a1001);
    double diff_d1_17_a	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_17_b	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_17	= diff_d1_17_a + 0.80 * (diff_d1_17_b - diff_d1_17_a);
    double diff_d2_17	= (a1101 - a1001) / (r1 - r0);
    double diff_d3_17_a	= (a1011 - a1001) / (d1 - d0);
    double diff_d3_17_b	= (a1111 - a1101) / (d1 - d0);
    double diff_d3_17	= diff_d3_17_a + 0.80 * (diff_d3_17_b - diff_d3_17_a);
    double diff_d4_17_a	= (a1001 - a1000) / (t1 - t0);
    double diff_d4_17_b	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_17	= diff_d4_17_a + 0.80 * (diff_d4_17_b - diff_d4_17_a);
    double xder17	= diff_d1_17 * dx4_dt + diff_d2_17 * dx3_dt +
                     diff_d3_17 * dx2_dt + diff_d4_17 * dx1_dt;
    check("4v f_table4Veq value    ", val17, res17,	1e-13);
    check("4v f_table4Veq d_dt     ", der17, xder17,	1e-10);
    val17 = P2.value(x4, x3, x2, x1);
    der17 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val17, res17,1e-13);
    check("4v f_table4Veq d_dt     ", der17, xder17,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val18 = P1.value(x4, x3, x2, x1);
    double der18 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res18 = a1021 + 0.80 * (a1121 - a1021);
    double diff_d1_18_a	= (a1021 - a0021) / (m1 - m0);
    double diff_d1_18_b	= (a1121 - a0121) / (m1 - m0);
    double diff_d1_18	= diff_d1_18_a + 0.80 * (diff_d1_18_b - diff_d1_18_a);
    double diff_d2_18	= (a1121 - a1021) / (r1 - r0);
    double diff_d3_18_a	= (a1021 - a1011) / (d2 - d1);
    double diff_d3_18_b	= (a1121 - a1111) / (d2 - d1);
    double diff_d3_18	= diff_d3_18_a + 0.80 * (diff_d3_18_b - diff_d3_18_a);
    double diff_d4_18_a	= (a1021 - a1020) / (t1 - t0);
    double diff_d4_18_b	= (a1121 - a1120) / (t1 - t0);
    double diff_d4_18	= diff_d4_18_a + 0.80 * (diff_d4_18_b - diff_d4_18_a);
    double xder18	= diff_d1_18 * dx4_dt +  diff_d2_18 * dx3_dt +
                     diff_d3_18 * dx2_dt + diff_d4_18 * dx1_dt;
    check("4v f_table4Veq value    ", val18, res18,	1e-13);
    check("4v f_table4Veq d_dt     ", der18, xder18,	1e-10);
    val18 = P2.value(x4, x3, x2, x1);
    der18 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val18, res18,1e-13);
    check("4v f_table4Veq d_dt     ", der18, xder18,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    7.0;	// (q1)
    double val19 = P1.value(x4, x3, x2, x1);
    double der19 = P1.d_dt( x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res19 = a1101 + 0.50 * (a1111 - a1101);
    double diff_d1_19_a	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_19_b	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_19	= diff_d1_19_a + 0.50 * (diff_d1_19_b - diff_d1_19_a);
    double diff_d2_19_a	= (a1101 - a1001) / (r1 - r0);
    double diff_d2_19_b	= (a1111 - a1011) / (r1 - r0);
    double diff_d2_19	= diff_d2_19_a + 0.50 * (diff_d2_19_b - diff_d2_19_a);
    double diff_d3_19	= (a1111 - a1101) / (d1 - d0);
    double diff_d4_19_a	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_19_b	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_19	= diff_d4_19_a + 0.50 * (diff_d4_19_b - diff_d4_19_a);
    double xder19	= diff_d1_19 * dx4_dt +  diff_d2_19 * dx3_dt +
                     diff_d3_19 * dx2_dt + diff_d4_19 * dx1_dt;
    check("4v f_table4Veq value    ", val19, res19,	1e-13);
    check("4v f_table4Veq d_dt     ", der19, xder19,	1e-10);
    val19 = P2.value(x4, x3, x2, x1);
    der19 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val19, res19,1e-13);
    check("4v f_table4Veq d_dt     ", der19, xder19,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val20 = P1.value(x4, x3, x2, x1);
    double der20 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res20 = a1120 + 0.25 * (a1121 - a1120);
    double diff_d1_20_a	= (a1120 - a0120) / (m1 - m0);
    double diff_d1_20_b	= (a1121 - a0121) / (m1 - m0);
    double diff_d1_20	= diff_d1_20_a + 0.25 * (diff_d1_20_b - diff_d1_20_a);
    double diff_d2_20_a	= (a1120 - a1020) / (r1 - r0);
    double diff_d2_20_b	= (a1121 - a1021) / (r1 - r0);
    double diff_d2_20	= diff_d2_20_a + 0.25 * (diff_d2_20_b - diff_d2_20_a);
    double diff_d3_20_a	= (a1120 - a1110) / (d1 - d0);
    double diff_d3_20_b	= (a1121 - a1111) / (d1 - d0);
    double diff_d3_20	= diff_d3_20_a + 0.25 * (diff_d3_20_b - diff_d3_20_a);
    double diff_d4_20	= (a1121 - a1120) / (t1 - t0);
    double xder20	= diff_d1_20 * dx4_dt + diff_d2_20 * dx3_dt +
                     diff_d3_20 * dx2_dt +  diff_d4_20 * dx1_dt;
    check("4v f_table4Veq value    ", val20, res20,	1e-13);
    check("4v f_table4Veq d_dt     ", der20, xder20,	1e-10);
    val20 = P2.value(x4, x3, x2, x1);
    der20 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val20, res20,1e-13);
    check("4v f_table4Veq d_dt     ", der20, xder20,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val21 = P1.value(x4, x3, x2, x1);
    double der21 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res21 = a1121;
    double diff_d1_21	= (a1121 - a0121) / (m1 - m0);
    double diff_d2_21	= (a1121 - a1021) / (r1 - r0);
    double diff_d3_21	= (a1121 - a1111) / (d2 - d1);
    double diff_d4_21	= (a1121 - a1120) / (t1 - t0);
    double xder21	= diff_d1_21 * dx4_dt +   diff_d2_21 * dx3_dt +
                     diff_d3_21 * dx2_dt +  diff_d4_21 * dx1_dt;
    check("4v f_table4Veq value    ", val21, res21,	1e-13);
    check("4v f_table4Veq d_dt     ", der21, xder21,	1e-10);
    val21 = P2.value(x4, x3, x2, x1);
    der21 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val21, res21,1e-13);
    check("4v f_table4Veq d_dt     ", der21, xder21,	1e-9);

    x4 = 1175.0;	// 25% less than (x0)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  107.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val22 = P1.value(x4, x3, x2, x1);
    double der22 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    temp3 = a0100 + 0.25 * (a0101 - a0100);
    temp4 = a0110 + 0.25 * (a0111 - a0110);
    temp5 = a1000 + 0.25 * (a1001 - a1000);
    temp6 = a1010 + 0.25 * (a1011 - a1010);
    temp7 = a1100 + 0.25 * (a1101 - a1100);
    temp8 = a1110 + 0.25 * (a1111 - a1110);
    temp9  = temp1 + 0.50 * (temp2 - temp1);
    temp10 = temp3 + 0.50 * (temp4 - temp3);
    temp11 = temp5 + 0.50 * (temp6 - temp5);
    temp12 = temp7 + 0.50 * (temp8 - temp7);
    temp13 = temp9  + 0.80 * (temp10 - temp9);
    temp14 = temp11 + 0.80 * (temp12 - temp11);
    double res22 = temp13 - 0.25 * (temp14 - temp13);
    // differential too complicated to compute manually
    check("4v f_table4Veq value    ", val22, res22,	1e-13);

} // closes test3v_f_table4Veq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void math::test::Tpre_fun::test4v_f_table4V() {


    // f_table4V of Hp[m] as function of m [kg], windi [mps], d[m], and windii[mps]
    // =============================================================================
    double dx2_dt = 5.0;				// dd_dt	= 5 mps
    double dx4_dt = 25.0;				// dm_dt	= 25 kg/sec
    double dx3_dt = 3.0;			// dDeltaT_dt
    double dx1_dt = 4.0;			// dDeltap_dt
    // Create mag objects
    double Hp1, Hp2;

    // 1st independent magnitude
    vec1 points1(2);
    double m0	= 1200.0;	points1.set(0, m0);
    double m1	= 1300.0;	points1.set(1, m1);

    // 2nd independent magnitude
    vec1 points2(2);
    double r0	= 25.4;		points2.set(0, r0);
    double r1	= 35.4;		points2.set(1, r1);

    // 3rd independent magnitude
    vec1 points3(3);
    double d0	=  82.1;	points3.set(0, d0);
    double d1	= 122.1;	points3.set(1, d1);
    double d2	= 182.1;	points3.set(2, d2);

    // 4th independent magnitude
    vec1 points4(2);
    double t0	=  3.0;		points4.set(0, t0);
    double t1	=  7.0;		points4.set(1, t1);

    // Dependent magnitudes
    double a0000 = 10, a0001 = 12, a0010 = 15, a0011 = 18, a0020 = 25, a0021 = 29;
    double a0100 = 60, a0101 = 56, a0110 = 65, a0111 = 58, a0120 = 85, a0121 = 94;
    double a1000 =210, a1001 =216, a1010 =215, a1011 =211, a1020 =230, a1021 =245;
    double a1100 =260, a1101 =215, a1110 =265, a1111 =287, a1120 =294, a1121 =384;
    vec4 VVValues(2, 3, 2, 2);
    VVValues.set(0, 0, 0, 0, a0000);	VVValues.set(1, 0, 0, 0, a0001);
    VVValues.set(0, 1, 0, 0, a0010);	VVValues.set(1, 1, 0, 0, a0011);
    VVValues.set(0, 2, 0, 0, a0020);	VVValues.set(1, 2, 0, 0, a0021);
    VVValues.set(0, 0, 1, 0, a0100);	VVValues.set(1, 0, 1, 0, a0101);
    VVValues.set(0, 1, 1, 0, a0110);	VVValues.set(1, 1, 1, 0, a0111);
    VVValues.set(0, 2, 1, 0, a0120);	VVValues.set(1, 2, 1, 0, a0121);
    VVValues.set(0, 0, 0, 1, a1000);	VVValues.set(1, 0, 0, 1, a1001);
    VVValues.set(0, 1, 0, 1, a1010);	VVValues.set(1, 1, 0, 1, a1011);
    VVValues.set(0, 2, 0, 1, a1020);	VVValues.set(1, 2, 0, 1, a1021);
    VVValues.set(0, 0, 1, 1, a1100);	VVValues.set(1, 0, 1, 1, a1101);
    VVValues.set(0, 1, 1, 1, a1110);	VVValues.set(1, 1, 1, 1, a1111);
    VVValues.set(0, 2, 1, 1, a1120);	VVValues.set(1, 2, 1, 1, a1121);

    // Create fun object
    f_table4V P1(points1, points2, points3, points4, VVValues, math::logic::lagrange_first_precompute);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi4(), true);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi2(), false);
    check("4v-000 f_table4V is equispaced      ", P1.get_equi1(), true);
    f_table4V P2(points1, points2, points3, points4, VVValues, math::logic::lagrange_first);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi4(), true);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi2(), false);
    check("4v-000 f_table4V is equispaced      ", P2.get_equi1(), true);
    f_table4V P6(points1, points2, points3, points4, VVValues, math::logic::hermite_first);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi4(), true);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi3(), true);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi2(), false);
    check("4v-000 f_table4V is equispaced      ", P6.get_equi1(), true);

    double x4 = 1200.0;		// (x0)
    double x3 =   25.4;	// (y0)
    double x2 =   82.1;		// (z0)
    double x1 =    3.0;	// (q0)
    double val1 = P1.value(x4, x3, x2, x1);
    double der1 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res1 = a0000;
    double diff_d1_1	= (a1000 - a0000) / (m1 - m0);
    double diff_d2_1	= (a0100 - a0000) / (r1 - r0);
    double diff_d3_1	= (a0010 - a0000) / (d1 - d0);
    double diff_d4_1	= (a0001 - a0000) / (t1 - t0);
    double xder1	= diff_d1_1 * dx4_dt +  diff_d2_1 * dx3_dt +
                    diff_d3_1 * dx2_dt +  diff_d4_1 * dx1_dt;
    check("4v f_table4V value      ", val1, res1,	1e-15);
    check("4v f_table4V d_dt       ", der1,	xder1,	1e-11);
    val1 = P2.value(x4, x3, x2, x1);
    der1 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val1, res1,1e-13);
    check("4v f_table4Veq d_dt     ", der1,	xder1,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =   82.1;		// (z0)
    x1 =    3.0;	// (q0)
    double val2 = P1.value(x4, x3, x2, x1);
    double der2 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res2 = a0100;
    double diff_d1_2	= (a1100 - a0100) / (m1 - m0);
    double diff_d2_2	= (a0100 - a0000) / (r1 - r0);
    double diff_d3_2	= (a0110 - a0100) / (d1 - d0);
    double diff_d4_2	= (a0101 - a0100) / (t1 - t0);
    double xder2	= diff_d1_2 * dx4_dt + diff_d2_2 * dx3_dt +
                    diff_d3_2 * dx2_dt +   diff_d4_2 * dx1_dt;
    check("4v f_table4V value      ", val2, res2,	1e-14);
    check("4v f_table4V d_dt       ", der2,	xder2,	1e-11);
    val2 = P2.value(x4, x3, x2, x1);
    der2 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val2, res2,1e-13);
    check("4v f_table4Veq d_dt     ", der2,	xder2,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    3.0;	// (q0)
    double val3 = P1.value(x4, x3, x2, x1);
    double der3 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res3 = a0020;
    double diff_d1_3	= (a1020 - a0020) / (m1 - m0);
    double diff_d2_3	= (a0120 - a0020) / (r1 - r0);
    double diff_d3_3	= (a0020 - a0010) / (d2 - d1);
    double diff_d4_3	= (a0021 - a0020) / (t1 - t0);
    double xder3	= diff_d1_3 * dx4_dt +   diff_d2_3 * dx3_dt +
                    diff_d3_3 * dx2_dt +   diff_d4_3 * dx1_dt;
    check("4v f_table4V value      ", val3, res3,	1e-14);
    check("4v f_table4V d_dt       ", der3,	xder3,	1e-11);
    val3 = P2.value(x4, x3, x2, x1);
    der3 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val3, res3,1e-13);
    check("4v f_table4Veq d_dt     ", der3, xder3,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  122.1;		// (z1)
    x1 =    7.0;	// (q1)
    double val4 = P1.value(x4, x3, x2, x1);
    double der4 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res4 = a0011;
    double diff_d1_4	= (a1011 - a0011) / (m1 - m0);
    double diff_d2_4	= (a0111 - a0011) / (r1 - r0);
    double diff_d3_4	= (a0021 - a0011) / (d2 - d1);
    double diff_d4_4	= (a0011 - a0010) / (t1 - t0);
    double xder4	= diff_d1_4 * dx4_dt + diff_d2_4 * dx3_dt +
                    diff_d3_4 * dx2_dt +   diff_d4_4 * dx1_dt;
    check("4v f_table4V value      ", val4, res4,	1e-14);
    check("4v f_table4V d_dt       ", der4,	xder4,	1e-11);
    val4 = P2.value(x4, x3, x2, x1);
    der4 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val4, res4,1e-13);
    check("4v f_table4Veq d_dt     ", der4,	xder4,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    3.0;	// (q0)
    double val5 = P1.value(x4, x3, x2, x1);
    double der5 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res5 = a0120;
    double diff_d1_5	= (a1120 - a0120) / (m1 - m0);
    double diff_d2_5	= (a0120 - a0020) / (r1 - r0);
    double diff_d3_5	= (a0120 - a0110) / (d2 - d1);
    double diff_d4_5	= (a0121 - a0120) / (t1 - t0);
    double xder5	= diff_d1_5 * dx4_dt +    diff_d2_5 * dx3_dt +
                    diff_d3_5 * dx2_dt +  diff_d4_5 * dx1_dt;
    check("4v f_table4V value      ", val5, res5,	1e-14);
    check("4v f_table4V d_dt       ", der5,	xder5,	1e-11);
    val5 = P2.value(x4, x3, x2, x1);
    der5 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val5, res5,1e-13);
    check("4v f_table4Veq d_dt     ", der5,	xder5,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  122.1;		// (z1)
    x1 =    7.0;	// (q1)
    double val6 = P1.value(x4, x3, x2, x1);
    double der6 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res6 = a0111;
    double diff_d1_6	= (a1111 - a0111) / (m1 - m0);
    double diff_d2_6	= (a0111 - a0011) / (r1 - r0);
    double diff_d3_6	= (a0121 - a0111) / (d2 - d1);
    double diff_d4_6	= (a0111 - a0110) / (t1 - t0);
    double xder6	= diff_d1_6 * dx4_dt +  diff_d2_6 * dx3_dt +
                    diff_d3_6 * dx2_dt + diff_d4_6 * dx1_dt;
    check("4v f_table4V value      ", val6, res6,	1e-14);
    check("4v f_table4V d_dt       ", der6,	xder6,	1e-11);
    val6 = P2.value(x4, x3, x2, x1);
    der6 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val6, res6,1e-13);
    check("4v f_table4Veq d_dt     ", der6,	xder6,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val7 = P1.value(x4, x3, x2, x1);
    double der7 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res7 = a0021;
    double diff_d1_7	= (a1021 - a0021) / (m1 - m0);
    double diff_d2_7	= (a0121 - a0021) / (r1 - r0);
    double diff_d3_7	= (a0021 - a0011) / (d2 - d1);
    double diff_d4_7	= (a0021 - a0020) / (t1 - t0);
    double xder7	= diff_d1_7 * dx4_dt +  diff_d2_7 * dx3_dt +
                    diff_d3_7 * dx2_dt +  diff_d4_7 * dx1_dt;
    check("4v f_table4V value      ", val7, res7,	1e-14);
    check("4v f_table4V d_dt       ", der7,	xder7,	1e-11);
    val7 = P2.value(x4, x3, x2, x1);
    der7 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val7, res7,1e-13);
    check("4v f_table4Veq d_dt     ", der7,	xder7,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val8 = P1.value(x4, x3, x2, x1);
    double der8 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res8 = a0121;
    double diff_d1_8	= (a1121 - a0121) / (m1 - m0);
    double diff_d2_8	= (a0121 - a0021) / (r1 - r0);
    double diff_d3_8	= (a0121 - a0111) / (d2 - d1);
    double diff_d4_8	= (a0121 - a0120) / (t1 - t0);
    double xder8	= diff_d1_8 * dx4_dt +   diff_d2_8 * dx3_dt +
                    diff_d3_8 * dx2_dt +  diff_d4_8 * dx1_dt;
    check("4v f_table4V value      ", val8, res8,	1e-14);
    check("4v f_table4V d_dt       ", der8,	xder8,	1e-10);
    val8 = P2.value(x4, x3, x2, x1);
    der8 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val8, res8,1e-13);
    check("4v f_table4Veq d_dt     ", der8,	xder8,	1e-9);

    val8 = P6.value(x4, x3, x2, x1);
    check("4v f_table4V value      ", val8, res8,	1e-12);

    x4 = 1200.0;	// (x0)
    x3 =   25.4;	// (y0)
    x2 =   82.1;		// (z0)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val9 = P1.value(x4, x3, x2, x1);
    double der9 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res9 = a0000 + 0.25 * (a0001 - a0000);
    double diff_d1_9_a	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_9_b	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_9	= diff_d1_9_a + 0.25 * (diff_d1_9_b - diff_d1_9_a);
    double diff_d2_9_a	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_9_b	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_9	= diff_d2_9_a + 0.25 * (diff_d2_9_b - diff_d2_9_a);
    double diff_d3_9_a	= (a0010 - a0000) / (d1 - d0);
    double diff_d3_9_b	= (a0011 - a0001) / (d1 - d0);
    double diff_d3_9	= diff_d3_9_a + 0.25 * (diff_d3_9_b - diff_d3_9_a);
    double diff_d4_9	= (a0001 - a0000) / (t1 - t0);
    double xder9	= diff_d1_9 * dx4_dt +   diff_d2_9 * dx3_dt +
                    diff_d3_9 * dx2_dt + diff_d4_9 * dx1_dt;
    check("4v f_table4V value      ", val9, res9,	1e-14);
    check("4v f_table4V d_dt       ", der9,	xder9,	1e-11);
    val9 = P2.value(x4, x3, x2, x1);
    der9 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val9, res9,1e-13);
    check("4v f_table4Veq d_dt     ", der9,	xder9,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    3.0;	// (q0)
    double val10 = P1.value(x4, x3, x2, x1);
    double der10 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res10 = a0100 + 0.50 * (a0110 - a0100);
    double diff_d1_10_a	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_10_b	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_10	= diff_d1_10_a + 0.50 * (diff_d1_10_b - diff_d1_10_a);
    double diff_d2_10_a	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_10_b	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_10	= diff_d2_10_a + 0.50 * (diff_d2_10_b - diff_d2_10_a);
    double diff_d3_10	= (a0110 - a0100) / (d1 - d0);
    double diff_d4_10_a	= (a0101 - a0100) / (t1 - t0);
    double diff_d4_10_b	= (a0111 - a0110) / (t1 - t0);
    double diff_d4_10	= diff_d4_10_a + 0.50 * (diff_d4_10_b - diff_d4_10_a);
    double xder10	= diff_d1_10 * dx4_dt +  diff_d2_10 * dx3_dt +
                     diff_d3_10 * dx2_dt + diff_d4_10 * dx1_dt;
    check("4v f_table4V value      ", val10, res10,	1e-14);
    check("4v f_table4V d_dt       ", der10, xder10,	1e-11);
    val10 = P2.value(x4, x3, x2, x1);
    der10 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val10, res10,1e-13);
    check("4v f_table4Veq d_dt     ", der10, xder10,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   35.4;	// (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val11 = P1.value(x4, x3, x2, x1);
    double der11 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double temp1 = a0100 + 0.25 * (a0101 - a0100);
    double temp2 = a0110 + 0.25 * (a0111 - a0110);
    double res11 = temp1 + 0.50 * (temp2 - temp1);
    double diff_d1_11_aa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_11_ab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_11_a		= diff_d1_11_aa + 0.25 * (diff_d1_11_ab - diff_d1_11_aa);
    double diff_d1_11_ba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_11_bb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_11_b		= diff_d1_11_ba + 0.25 * (diff_d1_11_bb - diff_d1_11_ba);
    double diff_d1_11		= diff_d1_11_a  + 0.50 * (diff_d1_11_b  - diff_d1_11_a);
    double diff_d2_11_aa	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_11_ab	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_11_a		= diff_d2_11_aa + 0.25 * (diff_d2_11_ab - diff_d2_11_aa);
    double diff_d2_11_ba	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_11_bb	= (a0111 - a0011) / (r1 - r0);
    double diff_d2_11_b		= diff_d2_11_ba + 0.25 * (diff_d2_11_bb - diff_d2_11_ba);
    double diff_d2_11		= diff_d2_11_a  + 0.50 * (diff_d2_11_b  - diff_d2_11_a);
    double diff_d3_11_a		= (a0110 - a0100) / (d1 - d0);
    double diff_d3_11_b		= (a0111 - a0101) / (d1 - d0);
    double diff_d3_11		= diff_d3_11_a + 0.25 * (diff_d3_11_b - diff_d3_11_a);
    double diff_d4_11_a		= (a0101 - a0100) / (t1 - t0);
    double diff_d4_11_b		= (a0111 - a0110) / (t1 - t0);
    double diff_d4_11		= diff_d4_11_a + 0.50 * (diff_d4_11_b - diff_d4_11_a);
    double xder11	= diff_d1_11 * dx4_dt +  diff_d2_11 * dx3_dt +
                     diff_d3_11 * dx2_dt + diff_d4_11 * dx1_dt;
    check("4v f_table4V value      ", val11, res11,	1e-14);
    check("4v f_table4V d_dt       ", der11, xder11,	1e-11);
    val11 = P2.value(x4, x3, x2, x1);
    der11 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val11, res11,1e-13);
    check("4v f_table4Veq d_dt     ", der11, xder11,	1e-9);

    x4 = 1200.0;	// (x0)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val12 = P1.value(x4, x3, x2, x1);
    double der12 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    double temp3 = a0100 + 0.25 * (a0101 - a0100);
    double temp4 = a0110 + 0.25 * (a0111 - a0110);
    double temp5 = temp1 + 0.50 * (temp2 - temp1);
    double temp6 = temp3 + 0.50 * (temp4 - temp3);
    double res12 = temp5 + 0.80 * (temp6 - temp5);
    double diff_d1_12_aaa	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_12_aab	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_12_aa	= diff_d1_12_aaa + 0.25 * (diff_d1_12_aab - diff_d1_12_aaa);
    double diff_d1_12_aba	= (a1010 - a0010) / (m1 - m0);
    double diff_d1_12_abb	= (a1011 - a0011) / (m1 - m0);
    double diff_d1_12_ab	= diff_d1_12_aba + 0.25 * (diff_d1_12_abb - diff_d1_12_aba);
    double diff_d1_12_baa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_12_bab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_12_ba	= diff_d1_12_baa + 0.25 * (diff_d1_12_bab - diff_d1_12_baa);
    double diff_d1_12_bba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_12_bbb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_12_bb	= diff_d1_12_bba + 0.25 * (diff_d1_12_bbb - diff_d1_12_bba);
    double diff_d1_12_a		= diff_d1_12_aa + 0.50 * (diff_d1_12_ab - diff_d1_12_aa);
    double diff_d1_12_b		= diff_d1_12_ba + 0.50 * (diff_d1_12_bb - diff_d1_12_ba);
    double diff_d1_12		= diff_d1_12_a  + 0.80 * (diff_d1_12_b - diff_d1_12_a);
    double diff_d2_12_aa	= (a0100 - a0000) / (r1 - r0);
    double diff_d2_12_ab	= (a0101 - a0001) / (r1 - r0);
    double diff_d2_12_a		= diff_d2_12_aa + 0.25 * (diff_d2_12_ab - diff_d2_12_aa);
    double diff_d2_12_ba	= (a0110 - a0010) / (r1 - r0);
    double diff_d2_12_bb	= (a0111 - a0011) / (r1 - r0);
    double diff_d2_12_b		= diff_d2_12_ba + 0.25 * (diff_d2_12_bb - diff_d2_12_ba);
    double diff_d2_12		= diff_d2_12_a  + 0.50 * (diff_d2_12_b  - diff_d2_12_a);
    double diff_d3_12_aa	= (a0010 - a0000) / (d1 - d0);
    double diff_d3_12_ab	= (a0011 - a0001) / (d1 - d0);
    double diff_d3_12_a		= diff_d3_12_aa + 0.25 * (diff_d3_12_ab - diff_d3_12_aa);
    double diff_d3_12_ba	= (a0110 - a0100) / (d1 - d0);
    double diff_d3_12_bb	= (a0111 - a0101) / (d1 - d0);
    double diff_d3_12_b		= diff_d3_12_ba + 0.25 * (diff_d3_12_bb - diff_d3_12_ba);
    double diff_d3_12		= diff_d3_12_a  + 0.80 * (diff_d3_12_b  - diff_d3_12_a);
    double diff_d4_12_aa	= (a0001 - a0000) / (t1 - t0);
    double diff_d4_12_ab	= (a0011 - a0010) / (t1 - t0);
    double diff_d4_12_a		= diff_d4_12_aa + 0.50 * (diff_d4_12_ab - diff_d4_12_aa);
    double diff_d4_12_ba	= (a0101 - a0100) / (t1 - t0);
    double diff_d4_12_bb	= (a0111 - a0110) / (t1 - t0);
    double diff_d4_12_b		= diff_d4_12_ba + 0.50 * (diff_d4_12_bb - diff_d4_12_ba);
    double diff_d4_12		= diff_d4_12_a  + 0.80 * (diff_d4_12_b  - diff_d4_12_a);
    double xder12	= diff_d1_12 * dx4_dt +  diff_d2_12 * dx3_dt +
                     diff_d3_12 * dx2_dt +   diff_d4_12 * dx1_dt;
    check("4v f_table4V value      ", val12, res12,	1e-14);
    check("4v f_table4V d_dt       ", der12, xder12,	1e-11);
    val12 = P2.value(x4, x3, x2, x1);
    der12 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val12, res12,1e-13);
    check("4v f_table4Veq d_dt     ", der12, xder12,	1e-9);

    x4 = 1210.0;	// 10% between (x0) and (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val13 = P1.value(x4, x3, x2, x1);
    // differential too complicated to compute manually
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    temp3 = a0100 + 0.25 * (a0101 - a0100);
    temp4 = a0110 + 0.25 * (a0111 - a0110);
    temp5 = a1000 + 0.25 * (a1001 - a1000);
    temp6 = a1010 + 0.25 * (a1011 - a1010);
    double temp7 = a1100 + 0.25 * (a1101 - a1100);
    double temp8 = a1110 + 0.25 * (a1111 - a1110);
    double temp9  = temp1 + 0.50 * (temp2 - temp1);
    double temp10 = temp3 + 0.50 * (temp4 - temp3);
    double temp11 = temp5 + 0.50 * (temp6 - temp5);
    double temp12 = temp7 + 0.50 * (temp8 - temp7);
    double temp13 = temp9  + 0.80 * (temp10 - temp9);
    double temp14 = temp11 + 0.80 * (temp12 - temp11);
    double res13 = temp13 + 0.10 * (temp14 - temp13);
    check("4v f_table4V value      ", val13, res13,	1e-14);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val14 = P1.value(x4, x3, x2, x1);
    double der14 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a1000 + 0.25 * (a1001 - a1000);
    temp2 = a1010 + 0.25 * (a1011 - a1010);
    temp3 = a1100 + 0.25 * (a1101 - a1100);
    temp4 = a1110 + 0.25 * (a1111 - a1110);
    temp5 = temp1 + 0.50 * (temp2 - temp1);
    temp6 = temp3 + 0.50 * (temp4 - temp3);
    double res14 = temp5 + 0.80 * (temp6 - temp5);
    double diff_d1_14_aaa	= (a1000 - a0000) / (m1 - m0);
    double diff_d1_14_aab	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_14_aa	= diff_d1_14_aaa + 0.25 * (diff_d1_14_aab - diff_d1_14_aaa);
    double diff_d1_14_aba	= (a1010 - a0010) / (m1 - m0);
    double diff_d1_14_abb	= (a1011 - a0011) / (m1 - m0);
    double diff_d1_14_ab	= diff_d1_14_aba + 0.25 * (diff_d1_14_abb - diff_d1_14_aba);
    double diff_d1_14_baa	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_14_bab	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_14_ba	= diff_d1_14_baa + 0.25 * (diff_d1_14_bab - diff_d1_14_baa);
    double diff_d1_14_bba	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_14_bbb	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_14_bb	= diff_d1_14_bba + 0.25 * (diff_d1_14_bbb - diff_d1_14_bba);
    double diff_d1_14_a		= diff_d1_14_aa + 0.50 * (diff_d1_14_ab - diff_d1_14_aa);
    double diff_d1_14_b		= diff_d1_14_ba + 0.50 * (diff_d1_14_bb - diff_d1_14_ba);
    double diff_d1_14		= diff_d1_14_a  + 0.80 * (diff_d1_14_b - diff_d1_14_a);
    double diff_d2_14_aa	= (a1100 - a1000) / (r1 - r0);
    double diff_d2_14_ab	= (a1101 - a1001) / (r1 - r0);
    double diff_d2_14_a		= diff_d2_14_aa + 0.25 * (diff_d2_14_ab - diff_d2_14_aa);
    double diff_d2_14_ba	= (a1110 - a1010) / (r1 - r0);
    double diff_d2_14_bb	= (a1111 - a1011) / (r1 - r0);
    double diff_d2_14_b		= diff_d2_14_ba + 0.25 * (diff_d2_14_bb - diff_d2_14_ba);
    double diff_d2_14		= diff_d2_14_a  + 0.50 * (diff_d2_14_b  - diff_d2_14_a);
    double diff_d3_14_aa	= (a1010 - a1000) / (d1 - d0);
    double diff_d3_14_ab	= (a1011 - a1001) / (d1 - d0);
    double diff_d3_14_a		= diff_d3_14_aa + 0.25 * (diff_d3_14_ab - diff_d3_14_aa);
    double diff_d3_14_ba	= (a1110 - a1100) / (d1 - d0);
    double diff_d3_14_bb	= (a1111 - a1101) / (d1 - d0);
    double diff_d3_14_b		= diff_d3_14_ba + 0.25 * (diff_d3_14_bb - diff_d3_14_ba);
    double diff_d3_14		= diff_d3_14_a  + 0.80 * (diff_d3_14_b  - diff_d3_14_a);
    double diff_d4_14_aa	= (a1001 - a1000) / (t1 - t0);
    double diff_d4_14_ab	= (a1011 - a1010) / (t1 - t0);
    double diff_d4_14_a		= diff_d4_14_aa + 0.50 * (diff_d4_14_ab - diff_d4_14_aa);
    double diff_d4_14_ba	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_14_bb	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_14_b		= diff_d4_14_ba + 0.50 * (diff_d4_14_bb - diff_d4_14_ba);
    double diff_d4_14		= diff_d4_14_a  + 0.80 * (diff_d4_14_b  - diff_d4_14_a);
    double xder14	= diff_d1_14 * dx4_dt +   diff_d2_14 * dx3_dt +
                     diff_d3_14 * dx2_dt +diff_d4_14 * dx1_dt;
    check("4v f_table4V value      ", val14, res14,	1e-13);
    check("4v f_table4V d_dt       ", der14, xder14,	1e-10);
    val14 = P2.value(x4, x3, x2, x1);
    der14 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val14, res14,1e-13);
    check("4v f_table4Veq d_dt     ", der14, xder14,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    3.0;	// (q0)
    double val15 = P1.value(x4, x3, x2, x1);
    double der15 = P1.d_dt( x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res15 = a1100 + 0.50 * (a1110 - a1100);
    double diff_d1_15_a	= (a1100 - a0100) / (m1 - m0);
    double diff_d1_15_b	= (a1110 - a0110) / (m1 - m0);
    double diff_d1_15	= diff_d1_15_a + 0.50 * (diff_d1_15_b - diff_d1_15_a);
    double diff_d2_15_a	= (a1100 - a1000) / (r1 - r0);
    double diff_d2_15_b	= (a1110 - a1010) / (r1 - r0);
    double diff_d2_15	= diff_d2_15_a + 0.50 * (diff_d2_15_b - diff_d2_15_a);
    double diff_d3_15	= (a1110 - a1100) / (d1 - d0);
    double diff_d4_15_a	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_15_b	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_15	= diff_d4_15_a + 0.50 * (diff_d4_15_b - diff_d4_15_a);
    double xder15	= diff_d1_15 * dx4_dt + diff_d2_15 * dx3_dt +
                     diff_d3_15 * dx2_dt +diff_d4_15 * dx1_dt;
    check("4v f_table4V value      ", val15, res15,	1e-13);
    check("4v f_table4V d_dt       ", der15, xder15,	1e-10);
    val15 = P2.value(x4, x3, x2, x1);
    der15 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val15, res15,1e-13);
    check("4v f_table4Veq d_dt     ", der15, xder15,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   25.4;	// (y0)
    x2 =  182.1;		// (z2)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val16 = P1.value(x4, x3, x2, x1);
    double der16 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res16 = a1020 + 0.25 * (a1021 - a1020);
    double diff_d1_16_a	= (a1020 - a0020) / (m1 - m0);
    double diff_d1_16_b	= (a1021 - a0021) / (m1 - m0);
    double diff_d1_16	= diff_d1_16_a + 0.25 * (diff_d1_16_b - diff_d1_16_a);
    double diff_d2_16_a	= (a1120 - a1020) / (r1 - r0);
    double diff_d2_16_b	= (a1121 - a1021) / (r1 - r0);
    double diff_d2_16	= diff_d2_16_a + 0.25 * (diff_d2_16_b - diff_d2_16_a);
    double diff_d3_16_a	= (a1020 - a1010) / (d2 - d1);
    double diff_d3_16_b	= (a1021 - a1011) / (d2 - d1);
    double diff_d3_16	= diff_d3_16_a + 0.25 * (diff_d3_16_b - diff_d3_16_a);
    double diff_d4_16	= (a1021 - a1020) / (t1 - t0);
    double xder16	= diff_d1_16 * dx4_dt +  diff_d2_16 * dx3_dt +
                     diff_d3_16 * dx2_dt +  diff_d4_16 * dx1_dt;
    check("4v f_table4V value      ", val16, res16,	1e-13);
    check("4v f_table4V d_dt       ", der16, xder16,	1e-10);
    val16 = P2.value(x4, x3, x2, x1);
    der16 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val16, res16,1e-13);
    check("4v f_table4Veq d_dt     ", der16, xder16,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =   82.1;		// (z0)
    x1 =    7.0;	// (q1)
    double val17 = P1.value(x4, x3, x2, x1);
    double der17 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res17 = a1001 + 0.80 * (a1101 - a1001);
    double diff_d1_17_a	= (a1001 - a0001) / (m1 - m0);
    double diff_d1_17_b	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_17	= diff_d1_17_a + 0.80 * (diff_d1_17_b - diff_d1_17_a);
    double diff_d2_17	= (a1101 - a1001) / (r1 - r0);
    double diff_d3_17_a	= (a1011 - a1001) / (d1 - d0);
    double diff_d3_17_b	= (a1111 - a1101) / (d1 - d0);
    double diff_d3_17	= diff_d3_17_a + 0.80 * (diff_d3_17_b - diff_d3_17_a);
    double diff_d4_17_a	= (a1001 - a1000) / (t1 - t0);
    double diff_d4_17_b	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_17	= diff_d4_17_a + 0.80 * (diff_d4_17_b - diff_d4_17_a);
    double xder17	= diff_d1_17 * dx4_dt +  diff_d2_17 * dx3_dt +
                     diff_d3_17 * dx2_dt +   diff_d4_17 * dx1_dt;
    check("4v f_table4V value      ", val17, res17,	1e-13);
    check("4v f_table4V d_dt       ", der17, xder17,	1e-10);
    val17 = P2.value(x4, x3, x2, x1);
    der17 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val17, res17,1e-13);
    check("4v f_table4Veq d_dt     ", der17, xder17,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val18 = P1.value(x4, x3, x2, x1);
    double der18 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res18 = a1021 + 0.80 * (a1121 - a1021);
    double diff_d1_18_a	= (a1021 - a0021) / (m1 - m0);
    double diff_d1_18_b	= (a1121 - a0121) / (m1 - m0);
    double diff_d1_18	= diff_d1_18_a + 0.80 * (diff_d1_18_b - diff_d1_18_a);
    double diff_d2_18	= (a1121 - a1021) / (r1 - r0);
    double diff_d3_18_a	= (a1021 - a1011) / (d2 - d1);
    double diff_d3_18_b	= (a1121 - a1111) / (d2 - d1);
    double diff_d3_18	= diff_d3_18_a + 0.80 * (diff_d3_18_b - diff_d3_18_a);
    double diff_d4_18_a	= (a1021 - a1020) / (t1 - t0);
    double diff_d4_18_b	= (a1121 - a1120) / (t1 - t0);
    double diff_d4_18	= diff_d4_18_a + 0.80 * (diff_d4_18_b - diff_d4_18_a);
    double xder18	= diff_d1_18 * dx4_dt +   diff_d2_18 * dx3_dt +
                     diff_d3_18 * dx2_dt + diff_d4_18 * dx1_dt;
    check("4v f_table4V value      ", val18, res18,	1e-13);
    check("4v f_table4V d_dt       ", der18, xder18,	1e-10);
    val18 = P2.value(x4, x3, x2, x1);
    der18 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val18, res18,1e-13);
    check("4v f_table4Veq d_dt     ", der18, xder18,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    7.0;	// (q1)
    double val19 = P1.value(x4, x3, x2, x1);
    double der19 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res19 = a1101 + 0.50 * (a1111 - a1101);
    double diff_d1_19_a	= (a1101 - a0101) / (m1 - m0);
    double diff_d1_19_b	= (a1111 - a0111) / (m1 - m0);
    double diff_d1_19	= diff_d1_19_a + 0.50 * (diff_d1_19_b - diff_d1_19_a);
    double diff_d2_19_a	= (a1101 - a1001) / (r1 - r0);
    double diff_d2_19_b	= (a1111 - a1011) / (r1 - r0);
    double diff_d2_19	= diff_d2_19_a + 0.50 * (diff_d2_19_b - diff_d2_19_a);
    double diff_d3_19	= (a1111 - a1101) / (d1 - d0);
    double diff_d4_19_a	= (a1101 - a1100) / (t1 - t0);
    double diff_d4_19_b	= (a1111 - a1110) / (t1 - t0);
    double diff_d4_19	= diff_d4_19_a + 0.50 * (diff_d4_19_b - diff_d4_19_a);
    double xder19	= diff_d1_19 * dx4_dt + diff_d2_19 * dx3_dt +
                     diff_d3_19 * dx2_dt +  diff_d4_19 * dx1_dt;
    check("4v f_table4V value      ", val19, res19,	1e-13);
    check("4v f_table4V d_dt       ", der19, xder19,	1e-10);
    val19 = P2.value(x4, x3, x2, x1);
    der19 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val19, res19,1e-13);
    check("4v f_table4Veq d_dt     ", der19, xder19,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val20 = P1.value(x4, x3, x2, x1);
    double der20 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res20 = a1120 + 0.25 * (a1121 - a1120);
    double diff_d1_20_a	= (a1120 - a0120) / (m1 - m0);
    double diff_d1_20_b	= (a1121 - a0121) / (m1 - m0);
    double diff_d1_20	= diff_d1_20_a + 0.25 * (diff_d1_20_b - diff_d1_20_a);
    double diff_d2_20_a	= (a1120 - a1020) / (r1 - r0);
    double diff_d2_20_b	= (a1121 - a1021) / (r1 - r0);
    double diff_d2_20	= diff_d2_20_a + 0.25 * (diff_d2_20_b - diff_d2_20_a);
    double diff_d3_20_a	= (a1120 - a1110) / (d2 - d1);
    double diff_d3_20_b	= (a1121 - a1111) / (d2 - d1);
    double diff_d3_20	= diff_d3_20_a + 0.25 * (diff_d3_20_b - diff_d3_20_a);
    double diff_d4_20	= (a1121 - a1120) / (t1 - t0);
    double xder20	= diff_d1_20 * dx4_dt +  diff_d2_20 * dx3_dt +
                     diff_d3_20 * dx2_dt +    diff_d4_20 * dx1_dt;
    check("4v f_table4V value      ", val20, res20,	1e-13);
    check("4v f_table4V d_dt       ", der20, xder20,	1e-10);
    val20 = P2.value(x4, x3, x2, x1);
    der20 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val20, res20,1e-13);
    check("4v f_table4Veq d_dt     ", der20, xder20,	1e-9);

    x4 = 1300.0;	// (x1)
    x3 =   35.4;	// (y1)
    x2 =  182.1;		// (z2)
    x1 =    7.0;	// (q1)
    double val21 = P1.value(x4, x3, x2, x1);
    double der21 = P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    double res21 = a1121;
    double diff_d1_21	= (a1121 - a0121) / (m1 - m0);
    double diff_d2_21	= (a1121 - a1021) / (r1 - r0);
    double diff_d3_21	= (a1121 - a1111) / (d2 - d1);
    double diff_d4_21	= (a1121 - a1120) / (t1 - t0);
    double xder21	= diff_d1_21 * dx4_dt + diff_d2_21 * dx3_dt +
                     diff_d3_21 * dx2_dt +  diff_d4_21 * dx1_dt;
    check("4v f_table4V value      ", val21, res21,	1e-13);
    check("4v f_table4V d_dt       ", der21, xder21,	1e-10);
    val21 = P2.value(x4, x3, x2, x1);
    der21 = P2.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    check("4v f_table4Veq value    ", val21, res21,1e-13);
    check("4v f_table4Veq d_dt     ", der21, xder21,	1e-9);

    x4 = 1175.0;	// 25% less than (x0)
    x3 =   33.4;	// 80% between (y0) and (y1)
    x2 =  102.1;		// 50% between (z0) and (z1)
    x1 =    4.0;	// 25% between (q0) and (q1)
    double val22 = P1.value(x4, x3, x2, x1);
    double der22= P1.d_dt(x4, x3, x2, x1, dx4_dt, dx3_dt, dx2_dt, dx1_dt);
    temp1 = a0000 + 0.25 * (a0001 - a0000);
    temp2 = a0010 + 0.25 * (a0011 - a0010);
    temp3 = a0100 + 0.25 * (a0101 - a0100);
    temp4 = a0110 + 0.25 * (a0111 - a0110);
    temp5 = a1000 + 0.25 * (a1001 - a1000);
    temp6 = a1010 + 0.25 * (a1011 - a1010);
    temp7 = a1100 + 0.25 * (a1101 - a1100);
    temp8 = a1110 + 0.25 * (a1111 - a1110);
    temp9  = temp1 + 0.50 * (temp2 - temp1);
    temp10 = temp3 + 0.50 * (temp4 - temp3);
    temp11 = temp5 + 0.50 * (temp6 - temp5);
    temp12 = temp7 + 0.50 * (temp8 - temp7);
    temp13 = temp9  + 0.80 * (temp10 - temp9);
    temp14 = temp11 + 0.80 * (temp12 - temp11);
    double res22 = temp13 - 0.25 * (temp14 - temp13);
    // differential too complicated to compute manually
    check("4v f_table4V value      ", val22, res22,	1e-13);


} // closes test3v_f_table4Veq

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

















