#ifndef MATH_SHARE
#define MATH_SHARE

#include "../math.h"
#include "logic.h"
#include <string>

namespace math {

// CLASS SHARE
// ===========
// ===========

class MATH_API share {
public:
    /**< location of CONDOR_INPUT global variable */
    static std::string condor_input;
    /**< location of CONDOR_OUTPUT_HARD_DISK global variable */
    static std::string condor_output_hard_disk;
    /**< location of CONDOR_OUTPUT_DISK2 global variable */
    static std::string condor_output_disk2;
    /**< location of CONDOR_OUTPUT_DISK4 global variable */
    static std::string condor_output_disk4;
    /**< location of CONDOR_OUTPUT_DISK5 global variable */
    static std::string condor_output_disk5;
    /**< location of CONDOR_OUTPUT_DISK51 global variable */
    static std::string condor_output_disk51;
    /**< location of CONDOR_OUTPUT_DISK52 global variable */
    static std::string condor_output_disk52;
    /**< location of CONDOR_OUTPUT_DISK53 global variable */
    static std::string condor_output_disk53;
    /**< location of CONDOR_OUTPUT_DISK54 global variable */
    static std::string condor_output_disk54;
    /**< location of CONDOR_OUTPUT_DISK6 global variable */
    static std::string condor_output_disk6;
    /**< location of CONDOR_OUTPUT_DISK7 global variable */
    static std::string condor_output_disk7;

    /**< return string with folder location based on enumeration */
    static std::string get_location(math::logic::FOLDER folder);
}; // closes class share

} // closes namespace math

#endif
