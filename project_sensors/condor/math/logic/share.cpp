#include "share.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CLASS SHARE
// ===========
// ===========

#ifdef _WIN32
#	define PATHSEPARATOR '\\'
#else
#	define PATHSEPARATOR '/'
#endif

std::string get_env_path(const std::string& key) {
	char *val;
	std::string path;

	val = std::getenv(key.c_str()); // get environment variable
	if (val != 0) {
		path = val;
		if (path.empty() == false) {
			if (path[path.size()-1] != PATHSEPARATOR) {
				path += PATHSEPARATOR;
			}
		}
	}
	if (path.empty() == true) {
		std::cerr << "CONFIGURATION_PREFIX environment variable not set." << std::endl;
		throw "";
	}
	return path;
}
/* returns the full path of global variable of name "key" */
//} // closes generic namespace

std::string math::share::condor_input = get_env_path("CONDOR_INPUT");
/* location of CONDOR_INPUT global variable */

std::string math::share::condor_output_hard_disk = get_env_path("CONDOR_OUTPUT_HARD_DISK");
/* location of CONDOR_OUTPUT_HARD_DISK global variable */

std::string math::share::condor_output_disk2 = get_env_path("CONDOR_OUTPUT_DISK2");
/* location of CONDOR_OUTPUT_DISK2 global variable */

std::string math::share::condor_output_disk4 = get_env_path("CONDOR_OUTPUT_DISK4");
/* location of CONDOR_OUTPUT_DISK4 global variable */

std::string math::share::condor_output_disk5 = get_env_path("CONDOR_OUTPUT_DISK5");
/* location of CONDOR_OUTPUT_DISK5 global variable */

std::string math::share::condor_output_disk51 = get_env_path("CONDOR_OUTPUT_DISK51");
/* location of CONDOR_OUTPUT_DISK51 global variable */

std::string math::share::condor_output_disk52 = get_env_path("CONDOR_OUTPUT_DISK52");
/* location of CONDOR_OUTPUT_DISK52 global variable */

std::string math::share::condor_output_disk53 = get_env_path("CONDOR_OUTPUT_DISK53");
/* location of CONDOR_OUTPUT_DISK53 global variable */

std::string math::share::condor_output_disk54 = get_env_path("CONDOR_OUTPUT_DISK54");
/* location of CONDOR_OUTPUT_DISK54 global variable */

std::string math::share::condor_output_disk6 = get_env_path("CONDOR_OUTPUT_DISK6");
/* location of CONDOR_OUTPUT_DISK6 global variable */

std::string math::share::condor_output_disk7 = get_env_path("CONDOR_OUTPUT_DISK7");
/* location of CONDOR_OUTPUT_DISK7 global variable */

std::string math::share::get_location(math::logic::FOLDER folder) {
    switch (folder) {
        case math::logic::folder_hard_disk:
            return math::share::condor_output_hard_disk;
        case math::logic::folder_disk2:
            return math::share::condor_output_disk2;
        case math::logic::folder_disk4:
            return math::share::condor_output_disk4;
        case math::logic::folder_disk5:
            return math::share::condor_output_disk5;
        case math::logic::folder_disk51:
            return math::share::condor_output_disk51;
        case math::logic::folder_disk52:
            return math::share::condor_output_disk52;
        case math::logic::folder_disk53:
            return math::share::condor_output_disk53;
        case math::logic::folder_disk54:
            return math::share::condor_output_disk54;
        case math::logic::folder_disk6:
            return math::share::condor_output_disk6;
        case math::logic::folder_disk7:
            return math::share::condor_output_disk7;
        default:
            throw std::runtime_error("Incorrect output folder choice.");
    }
}
/* return string with folder location based on enumeration */

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////