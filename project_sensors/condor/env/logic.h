#ifndef ENV_LOGIC
#define ENV_LOGIC

#include "env.h"

namespace env {

namespace logic {

    /**< Enumeration that contains the different offsets models */
    enum OFFSETS_ID {
        offsets_id00   = 0,
        offsets_id01   = 1,
        offsets_guid01 = 2,
        offsets_guid03 = 3,
        offsets_guid04 = 4,
        offsets_size   = 5
    };

    /**< Enumeration that contains the different wind models */
    enum WIND_ID {
        wind_id00   = 0,
        wind_id01   = 1,
        wind_id02   = 2,
        wind_guid01 = 3,
        wind_guid03 = 4,
        wind_guid04 = 5,
        wind_size   = 6
    };

    /**< enumeration that contains the realism modes for gravitation and magnetism */
    enum REALISM_ID {
        realism_grav_yes_magn_yes = 0,
        realism_grav_yes_magn_no  = 1,
        realism_grav_no_magn_yes  = 2,
        realism_grav_no_magn_no   = 3,
        realism_size              = 4
    };

    /**< enumeration that contains the different flight zones */
    enum ZONE_ID {
        zone_default        = 0,
        zone_desert         = 1,
        zone_urban          = 2,
        zone_everglades     = 3,
        zone_forest         = 4,
        zone_farm           = 5,
        zone_mix            = 6,
        zone_mexico         = 7,
        zone_colombia       = 8,
        zone_argentina      = 9,
        zone_canada         = 10,
        zone_rozas          = 11,
        zone_wisconsin      = 12,
        zone_pntt           = 13,
        zone_river          = 14,
        zone_size           = 15
    };
} // closes namespace logic

} // closes namespace env

#endif








