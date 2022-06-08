#ifndef MATH_LOGIC_H
#define MATH_LOGIC_H

#include "../math.h"

namespace math {
namespace logic {

// ENUM PRED_NAME
// ==============
// ==============

enum PRED_NAME {
    f_null				= 0,	// 0 variables: null
    f_constant			= 1,	// 0 variables: constant
    f_lineal			= 2,	// 1 variable:  lineal
    f_parabolic			= 3,	// 1 variable:  parabolic
    f_cubic			    = 4,	// 1 variable:  cubic
    f_steps				= 5,	// 1 variable:  constant steps
    f_table1V			= 6,	// 1 variable:  uni dimensional table
    f_cat				= 7,	// 1 variable:  concatenation
    f_cat_partial		= 8,	// 1 variable:  concatenation #2
    f_table1V_num		= 9,	// 1 variable:: numeric uni dimensional table with linear interpolation (avoid use if possible)
    f_table1V_spl		= 10,	// 1 variable:: numeric uni dimensional table with spline interpolation (avoid use if possible)
    f_lineal_double		= 11,	// 2 variables: lineal
    f_table2V			= 12,	// 2 variables: bi dimensional table
    f_lineal_triple		= 13,	// 3 variables: lineal
    f_table3V			= 14,	// 3 variables: tri dimensional table
    f_tabular3V			= 15,	// 3 variables: vector of bidimensional tables
    f_table4V			= 16,	// 4 variables: quatri dimensional table
    f_tabular2V         = 17,
    pred_size			= 18	// contains size of enumeration
};
/**< enumeration that contains the names of the different predicates accepted by the fun objects. */

// ENUM INTERP_MODE
// ================
// ================

enum INTERP_MODE {
	lagrange_first_precompute	= 0, // Lagrange 1st order interpolation (based on 2 magnitudes) with the differentials precomputed
	lagrange_first				= 1, // Lagrange 1st order interpolation (based on 2 magnitudes)
	lagrange_second				= 2, // Lagrange 2nd order interpolation (based on 3 magnitudes)
	lagrange_third				= 3, // Lagrange 3rd order interpolation (based on 4 magnitudes)
	biparabolic					= 4, // biparabolic interpolation (based on 4 magnitudes)
	hermite_first				= 5, // Hermite cubic interpolation (differentials based on averages of Lagrange 1st order) (based on 2 magnitudes)
	hermite_second				= 6, // Hermite cubic interpolation (differentials based on Lagrange 2nd order) (based on 2 magnitudes)
	spline						= 7  // splines (based on 4 magnitudes) (differentials at border based on Lagrange 3rd order)
	//interpolation_size	= 8 // BETTER WITHOUT IT
};
/**< enumeration settings that indicates the interpolation method to employ for tables */

// ENUM FOLDER
// ===========
// ===========

enum FOLDER {
    folder_hard_disk = 0, // output files saved in CONDOR_OUTPUT_HARD_DISK global variable
    folder_disk2     = 1, // output files saved in CONDOR_OUTPUT_DISK2 global variable
    folder_disk4     = 2, // output files saved in CONDOR_OUTPUT_DISK4 global variable
    folder_disk5     = 3, // output files saved in CONDOR_OUTPUT_DISK5 global variable
    folder_disk51    = 4, // output files saved in CONDOR_OUTPUT_DISK51 global variable
    folder_disk52    = 5, // output files saved in CONDOR_OUTPUT_DISK52 global variable
    folder_disk53    = 6, // output files saved in CONDOR_OUTPUT_DISK53 global variable
    folder_disk54    = 7, // output files saved in CONDOR_OUTPUT_DISK54 global variable
    folder_disk6     = 8, // output files saved in CONDOR_OUTPUT_DISK6 global variable
    folder_disk7     = 9  // output files saved in CONDOR_OUTPUT_DISK7 global variable
};
/**< enumeration to specify the different places where output files may be saved */

} // closes namespace logic
} // closes namespace math

#endif













