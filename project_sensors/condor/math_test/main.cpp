
#include "Tpre_fun.h"
#include "Tinterp.h"
#include "Tlow_pass.h"
#include "Tmetrics.h"
#include "Tdistributions.h"
#include "Tdecompose.h"

int main(int argc, char **argv) {

    jail::counter Ocounter;

    {math::test::Tinterp       	o(Ocounter); o.run(); }	// interpolation
    {math::test::Tpre_fun     	o(Ocounter); o.run(); }	// functions and predicates
    {math::test::Tlow_pass     	o(Ocounter); o.run(); }	// low pass filter
    {math::test::Tmetrics      	o(Ocounter); o.run(); }	// metrics
    {math::test::Tdistributions	o(Ocounter); o.run(); }	// statistical distributions
    {math::test::Tdecompose 	o(Ocounter); o.run(); }	// linear solving, least squares

    Ocounter.write_results();

	return 0;
}



