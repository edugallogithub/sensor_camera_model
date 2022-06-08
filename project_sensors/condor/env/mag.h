#ifndef ENV_MAG
#define ENV_MAG

#include "env.h"
#include "logic.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <random>

namespace env {
    class geodetic_coord;
}
namespace ang {
    class euler;
}

namespace env {

/* The data is taken from https://www.ngdc.noaa.gov/geomag-web/#igrfwmm
 *
 * ===== ===== zone_desert ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -111.998815  +32.157903  +3000   01-Jan-2020	    23997.2     4200.7     39574.8      init    24362.1     46472.4     +09.9289    +58.3837
 * -112.0       +33.0       +3000   01-Jan-2020     23708.8     4205.8     40325.4      ul      24079.0     46967.4
 * -111.0       +33.0       +3000   01-Jan-2020     23672.4  	4064.6     40535.5      ur      24018.8     47117.2
 * -112.0       +32.0       +3000   01-Jan-2020     24050.6     4199.7     39432.9      ll      24414.5     46379.1
 * -111.0       +32.0       +3000   01-Jan-2020     24017.7     4060.8     39642.7      lr      24358.6     46528.3
 *
 * ===== ===== zone_urban ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -118.200269  +33.924426  +3000   01-Jan-2020     23615.6     4938.8      39741.3     init    24126.5     46491.5     +11.8122    +58.7385
 * -119.0       +34.0       +3000   01-Jan-2020     23618.6     5015.2      39617.6     ul      24145.2     46395.6
 * -118.0       +34.0       +3000   01-Jan-2020     23583.6     4920.5      39854.7     ur      24091.4     46570.3
 * -119.0       +33.0       +3000   01-Jan-2020     23940.3     4996.0      38740.2     ll      24456.0     45813.8
 * -118.0       +33.0       +3000   01-Jan-2020     23910.0     4903.5      38976.2     lr      24407.6     45987.8
 *
 * ===== ===== zone_everglades ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -080.911166  +25.855172  +3000   01-Jan-2020     25111.7     -2814.4     35592.9     init    25268.9		43650.6     -06.3947	+54.6275
 * -081.0       +26.0       +3000   01-Jan-2020     25077.3     -2795.3 	35743.1     ul      24921.0     43573.2
 * -080.0       +26.0       +3000   01-Jan-2020     25069.0     -3065.0 	35586.0     ur      24880.9     43421.5
 * -081.0       +25.0       +3000   01-Jan-2020     25314.9  	-2759.2 	34794.9     ll      25164.1     42940.9
 * -080.0       +25.0       +3000   01-Jan-2020     25302.5 	-3030.2 	34639.6     lr      25120.4     42789.4
 *
 * ===== ===== zone_forest ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -072.509195  +43.354486  +3000   01-Jan-2020     19085.9  	-4778.7 	48594.2     init    19675.1		52426.2     -14.0568	+67.9577
 * -073.0       +44.0       +3000   01-Jan-2020     18731.8 	-4671.3 	49143.0     ul      18140.0     52384.1
 * -072.0       +44.0       +3000   01-Jan-2020     18812.9 	-4846.3 	48922.8     ur      18178.0     52190.8
 * -073.0       +43.0       +3000   01-Jan-2020    	19218.1 	-4701.3 	48456.7     ll      18634.2     51916.1
 * -072.0       +43.0       +3000   01-Jan-2020     19297.3 	-4878.9 	48232.0     lr      18670.3     51719.5
 *
 * ===== ===== zone_farm ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -087.877629  +38.865625  +3000   01-Jan-2020     20565.2     -1213.0     47629.0     init    20601.0     51893.3     -03.3754	+66.6101
 * -088.0       +39.0       +3000   01-Jan-2020     20504.5 	-1182.5 	47737.6     ul      20470.4     51941.4
 * -087.0       +39.0       +3000   01-Jan-2020     20506.9 	-1449.3 	47679.8     ur      20455.6     51882.5
 * -088.0       +38.0       +3000   01-Jan-2020     20950.1 	-1165.7 	46967.4     ll      20917.6     51414.8
 * -087.0       +38.0       +3000   01-Jan-2020     20951.2 	-1434.3 	46909.8     lr      20902.0     51355.9
 *
 * ===== ===== zone_mix ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -089.015462  +34.720636  +3000   01-Jan-2020     22322.9 	- 827.9 	44347.8     init  	22338.2 	49656.0     -02.1239	+63.2654
 * -089.0       +35.0       +3000   01-Jan-2020     22211.8	    - 837.9 	44581.8     ul      22196.0     49801.6
 * -088.0       +35.0       +3000   01-Jan-2020     22204.8 	-1110.6 	44539.6     ur      22177.0     49755.3
 * -089.0       +34.0       +3000   01-Jan-2020     22603.9 	- 816.6 	43736.2     ll      22589.1     49225.2
 * -088.0       +34.0       +3000   01-Jan-2020     22595.3 	-1090.1 	43694.4     lr      22569.0     49178.9
 *
 * ===== ===== zone_mexico ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -100.0       +20.0       +3000   01-Jan-2020     27148.9     2274.2      29855.4     init    27244.0		40417.6     +04.7884	+47.6185
 * -100.5       +20.5       +3000   01-Jan-2020   	27056.8 	2357.6  	30307.8     ul      27159.4	 	40696.3
 * -099.5       +20.5       +3000   01-Jan-2020    	27017.5 	2161.0  	30432.4     ur	    27103.7     40752.2
 * -100.5       +19.5       +3000   01-Jan-2020    	27277.5 	2384.1  	29273.4     ll	    27381.5	    40083.4
 * -099.5       +19.5       +3000   01-Jan-2020   	27236.4 	2190.9      29399.1     lr	    27324.4     40136.4
 *
 * ===== ===== zone_colombia ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -070.0       +05.0       +3000   01-Jan-2020     26564.7	    -4850.2 	13468.4     init    27003.8		30176.2     -10.3472	+26.5082
 * -070.5       +05.5       +3000   01-Jan-2020     26615.5     -4728.6     14064.6     ul      27032.3		30472.3
 * -069.5       +05.5       +3000   01-Jan-2020     26571.5 	-5023.6 	13799.1     ur      27042.3		30359.5
 * -070.5       +04.5       +3000   01-Jan-2020     26558.1 	-4673.3 	13133.1     ll      26966.1		29994.1
 * -069.5       +04.5       +3000   01-Jan-2020     26509.5 	-4971.9 	12869.5     lr      26971.7 	29884.7
 *
 * ===== ===== zone_argentina ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -060.0       -25.0       +3000   01-Jan-2020     18914.5 	-4239.0 	-10858.6    init    19383.7		22217.9     -12.6320	-29.2572
 * -060.5       -24.5       +3000   01-Jan-2020     19120.4 	-4205.9 	-10532.6    ul  	19577.5     22230.9
 * -059.5       -24.5       +3000   01-Jan-2020     18926.5 	-4447.3 	-10769.9    ur  	19442.0 	22225.7
 * -060.5       -25.5       +3000   01-Jan-2020 	18905.6 	-4026.8 	-10953.1    ll  	19329.7 	22217.3
 * -059.5       -25.5       +3000   01-Jan-2020 	18711.0 	-4268.6 	-11178.7    lr      19191.7	   	22210.1
 *
 * ===== ===== zone_canada ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -105.0       +52.0       +3000   01-Jan-2020     14291.0	    2296.8  	55034.0     init    14474.4 	56905.6	    +09.1303    +75.2644
 * -105.5       +52.5       +3000   01-Jan-2020   	14038.7 	2354.1  	55209.6     ul      14234.7     57015.2
 * -104.5       +52.5       +3000   01-Jan-2020     13961.6 	2164.1  	55301.2     ur      14128.3 	57077.4
 * -105.5       +51.5       +3000   01-Jan-2020 	14619.9 	2428.6  	54753.5     ll      14820.3 	56723.8
 * -104.5	    +51.5       +3000   01-Jan-2020     14543.3 	2236.5  	54851.0     lr      14714.2 	56790.3
 *
 *
 *
 * ===== ===== zone_default ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)
 * -87.0        +45.0       +0      01-Jan-2018     17478.0     -1531.5     52134.2
 *
 * ===== ===== zone_wisconsin ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square
 * -90.0        +45.0       +3000   01-Jan-2018     17455.7      -764.6     52192.7     ul
 * -89.0        +45.0       +3000   01-Jan-2018     17449.5     -1019.6     52158.5     ur
 * -90.0        +44.0       +3000   01-Jan-2018     17983.9      -739.1     51571.5     ll
 * -89.0        +44.0       +3000   01-Jan-2018     17977.2      -997.0     51538.4     lr
 *
 * ===== ===== zone_rozas ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)
 * -7.5652        43.1142   1300    10-Jun-2019     24175.8     -749.4      38813.1
 *
 * ===== ===== zone_pntt ===== =====
 * Long (deg)   Lat (deg)   Alt (m) Date            Bnorth(nT)  Beast(nT)   Bdown(nT)   Square  Bhor(nT)    Btotal(nT)  Deast(de)   Idown(deg)
 * -91.1162     +38.7319    ¡0      01-Jan-2020     20677.1 	-346.8  	47696.0	            20680.0 	51986.3     -00.9610	+66.5595
 *
 */

/* ===== ===== ===== MAGNETISM ===== ===== =====
 * It is important to distinguish between the true magnetic field and that provided by a model, as the
 * former is the one measured by the magnetometers while the later represents our best knowledge of
 * it based on position and time.
 * The magnetic model should be computed according to the WMM, but instead of that I just evaluate
 * it (based on external implementation) at either a single point or the four corners of an area
 * encompassing the flight area, and then keep it constant or perform bidimensional linear interpolation
 * in between. The reason not to implement WMM myself is programming complexity and computing resources.
 * If that were the case, the results would be closer to reality, but that is not the point.
 * The important fact is to recognize that we have a model, not reality. We can assume with no loss of
 * generality that the differences between the model and reality are those of the complete WMM model,
 * instead of our interpolated implementation. Complying with this is just a matter of implementing WMM.
 * I base my modeling of the differences between the real and modelled magnetism on the numbers
 * provided by "The US/UK World Magnetic Model for 2015-2020", which provides standard deviations
 * of (138, 89, 165) [nT] for the magnetic field viewed in NED. I apply independent normal distributions
 * to the four corners of the grid (a total of 12 executions) and then proceed exactly as in the
 * case of the model, interpolating in two dimensions in between.
 */

// CLASS MAG
// =========
// =========

class ENV_API mag {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    mag() = default;
    /**< copy constructor */
    mag(const mag&) = delete;
    /**< move constructor */
    mag(mag&&) = delete;
    /**< destructor */
    virtual ~mag() = default;
    /**< copy assignment */
    mag& operator=(const mag&) = delete;
    /**< move assignment */
    mag& operator=(mag&&) = delete;

    /**< computes the MODEL Earth magnetic field based on time and geodetic coordinates */
    virtual Eigen::Vector3d compute_B_n_nT_model(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const = 0;
    /**< computes the REAL Earth magnetic field based on time and geodetic coordinates */
    virtual Eigen::Vector3d compute_B_n_nT_truth(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const = 0;
    /**< returns the magnetic Euler angles (declination, minus inclination, 0) based on the magnetic field */
    static ang::euler compute_euler_mag(const Eigen::Vector3d& mag_n_nT);

    /**< return pointer to magnetic object based on location, realism id, normal distribution, and seed generator. */
    static env::mag* create_mag(env::logic::ZONE_ID, env::logic::REALISM_ID, std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen);
    /**< describe magnetic model and realism in stream */
    virtual void create_text(std::ostream& Ostream, const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const = 0;
}; // closes class mag

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS MAG_CONSTANT
// ==================
// ==================

class ENV_API mag_constant : public mag {
private:
    /**< constant magnetic field in NED (representing model) */
    Eigen::Vector3d _B_n_nT_model;
    /**< constant magnetic field in NED (representing reality) */
    Eigen::Vector3d _B_n_nT_truth;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    mag_constant() = delete;
    /**< constructor based on NED magnetic field. Model and real magnetic fields both coincide with input. */
    explicit mag_constant(const Eigen::Vector3d& B_n_nT_model);
    /**< constructor based on NED magnetic field, normal distribution, and seed generator, and three NED standard deviations. */
    mag_constant(const Eigen::Vector3d& B_n_nT_model, std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen, const double& sigma_B_nT_north, const double& sigma_B_nT_east, const double& sigma_B_nT_down);
    /**< copy constructor */
    mag_constant(const mag_constant&) = delete;
    /**< move constructor */
    mag_constant(mag_constant&&) = delete;
    /**< destructor */
    virtual ~mag_constant() = default;
    /**< copy assignment */
    mag_constant& operator=(const mag_constant&) = delete;
    /**< move assignment */
    mag_constant& operator=(mag_constant&&) = delete;

    /**< computes the MODEL Earth magnetic field based on time and geodetic coordinates */
    Eigen::Vector3d compute_B_n_nT_model(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
    /**< computes the REAL Earth magnetic field based on time and geodetic coordinates */
    Eigen::Vector3d compute_B_n_nT_truth(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
    /**< describe magnetic model and realism in stream */
    void create_text(std::ostream& Ostream, const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
}; // closes class mag_constant

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// CLASS MAG_LINEAR
// ================
// ================

class ENV_API mag_linear : public mag {
private:
    /**< magnetic field in NED at lower left corner (representing model) */
    Eigen::Vector3d _B_n_nT_model_ll;
    /**< magnetic field in NED at lower right corner (representing model) */
    Eigen::Vector3d _B_n_nT_model_lr;
    /**< NED magnetic field difference between upper left and lower left corners (representing model) */
    Eigen::Vector3d _B_n_nT_model_ul_ll;
    /**< NED magnetic field difference between upper right and lower right corners (representing model) */
    Eigen::Vector3d _B_n_nT_model_ur_lr;

    /**< latitude of lower square corners */
    double _phi_rad_low;
    /**< longitude of left square corners */
    double _lambda_rad_left;
    /**< longitude difference between right and left square corner */
    double _Delta_lambda_rad;
    /**< latitude difference between upper and lower square corners */
    double _Delta_phi_rad;

    /**< magnetic field in NED at lower left corner (representing reality) */
    Eigen::Vector3d _B_n_nT_truth_ll;
    /**< magnetic field in NED at lower right corner (representing reality) */
    Eigen::Vector3d _B_n_nT_truth_lr;
    /**< NED magnetic field difference between upper left and lower left corners (representing reality) */
    Eigen::Vector3d _B_n_nT_truth_ul_ll;
    /**< NED magnetic field difference between upper right and lower right corners (representing reality) */
    Eigen::Vector3d _B_n_nT_truth_ur_lr;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**< default constructor */
    mag_linear() = delete;
    /**< constructor based on NED magnetic field and geodetic coordinates at four points forming square. Model and real magnetic
     * fields both coincide with each other and input. Altitude neglected. Linear interpolation in between and also outside square */
    mag_linear(const Eigen::Vector3d& B_n_nT_model_ul, const env::geodetic_coord& x_gdt_rad_m_ul,
               const Eigen::Vector3d& B_n_nT_model_ur, const env::geodetic_coord& x_gdt_rad_m_ur,
               const Eigen::Vector3d& B_n_nT_model_ll, const env::geodetic_coord& x_gdt_rad_m_ll,
               const Eigen::Vector3d& B_n_nT_model_lr, const env::geodetic_coord& x_gdt_rad_m_lr);
    /**< constructor based on NED magnetic field and geodetic coordinates at four points forming square, normal distribution,
     * seed generator, and three NED standard deviations. Altitude neglected. Linear interpolation in between and also outside square */
    mag_linear(const Eigen::Vector3d& B_n_nT_model_ul, const env::geodetic_coord& x_gdt_rad_m_ul,
               const Eigen::Vector3d& B_n_nT_model_ur, const env::geodetic_coord& x_gdt_rad_m_ur,
               const Eigen::Vector3d& B_n_nT_model_ll, const env::geodetic_coord& x_gdt_rad_m_ll,
               const Eigen::Vector3d& B_n_nT_model_lr, const env::geodetic_coord& x_gdt_rad_m_lr,
               std::normal_distribution<double>& dist_normal, std::ranlux24_base& gen,
               const double& sigma_B_nT_north, const double& sigma_B_nT_east, const double& sigma_B_nT_down);
    /**< copy constructor */
    mag_linear(const mag_linear&) = delete;
    /**< move constructor */
    mag_linear(mag_linear&&) = delete;
    /**< destructor */
    virtual ~mag_linear() = default;
    /**< copy assignment */
    mag_linear& operator=(const mag_linear&) = delete;
    /**< move assignment */
    mag_linear& operator=(mag_linear&&) = delete;

    /**< computes the MODEL Earth magnetic field based on time and geodetic coordinates */
    Eigen::Vector3d compute_B_n_nT_model(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
    /**< computes the REAL Earth magnetic field based on time and geodetic coordinates */
    Eigen::Vector3d compute_B_n_nT_truth(const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
    /**< describe magnetic model and realism in stream */
    void create_text(std::ostream& Ostream, const double& t_sec, const env::geodetic_coord& x_gdt_rad_m) const override;
}; // closes class mag_linear

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

} // closes namespace env

#endif
