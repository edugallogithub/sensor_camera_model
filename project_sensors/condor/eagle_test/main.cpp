#include "Teagle.h"
#include "Tplanet_eagle.h"
#include "Tcamera.h"
#include "Teagle_zones.h"
#include "Tmode.h"

#include <iostream>

int main(int argc, char **argv) {

    jail::counter Ocounter;

    //{eagle::test::Tcamera           o(Ocounter); o.run(); }
    //{eagle::test::Teagle            o(Ocounter); o.run(&argc, argv); }
    {eagle::test::Teagle_zones      o(Ocounter); o.run(&argc, argv); }
    //{eagle::test::Tplanet_eagle     o(Ocounter); o.run(&argc, argv); }
    //{eagle::test::Tmode            o(Ocounter); o.run(&argc, argv); }

    Ocounter.write_results();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

