
#include "Tso3.h"
#include "Tse3.h"
#include "Thowtouse.h"
#include "Ttools.h"
#include "Tintegration.h"
#include "Tscrew.h"
#include "Tso3_jacobian.h"
#include "Tse3_jacobian.h"
#include "Tintegr_compare.h"


#include <iostream>
#include <iomanip>


int main(int argc, char **argv) {

    jail::counter Ocounter;

    {ang::test::Tso3        	    o(Ocounter); o.run();}
    {ang::test::Tse3           	    o(Ocounter); o.run();}
    {ang::test::Tintegration   	    o(Ocounter); o.run();}
    {ang::test::Tscrew              o(Ocounter); o.run();}
    {ang::test::Thowtouse      	    o(Ocounter); o.run();}
    {ang::test::Ttools      	    o(Ocounter); o.run();}
    {ang::test::Tso3_jacobian 	    o(Ocounter); o.run();}
    {ang::test::Tse3_jacobian 	    o(Ocounter); o.run();}
    {ang::test::Tintegr_compare   	o(Ocounter); o.run();}

	Ocounter.write_results();

	return 0;
}



