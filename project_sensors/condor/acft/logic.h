#ifndef ACFT_LOGIC
#define ACFT_LOGIC

#include "acft.h"

namespace control {
namespace logic {
    /**< enumeration that contains the different controls */
    enum CNTR_ID {
        cntr_THR     = 0,
        cntr_ELV     = 1,
        cntr_AIL     = 2,
        cntr_RUD     = 3,
        cntr_id_size = 4
    };
    /**< Enumeration that contains the different throttle guidance modes */
    enum THR_ID {
        thr_vtas_mps = 0,
        thr_deltaT   = 1,
        thr_id_size  = 2
    };
    /**< Enumeration that contains the different elevator guidance modes */
    enum ELV_ID {
        elv_theta_deg    = 0,
        elv_Hp_m         = 1,
        elv_h_m          = 2,
        elv_gammaTAS_deg = 3,
        elv_id_size      = 4
    };
    /**< Enumeration that contains the different ailerons guidance modes */
    enum AIL_ID {
        ail_xi_deg     = 0,
        ail_chi_deg    = 1,
        ail_psi_deg    = 2,
        ail_muTAS_deg  = 3,
        ail_chiTAS_deg = 4,
        ail_id_size    = 5
    };
    /**< Enumeration that contains the different rudder guidance modes */
    enum RUD_ID {
        rud_beta_deg = 0,
        rud_id_size  = 1
    };
    /**< Enumeration that contains the different trigger guidance modes
     * (use trgg instead of trg so name does not coincide with trigger classes) */
    enum TRG_ID {
        trgg_t_sec      = 0,
        trgg_Deltat_sec = 1,
        trgg_h_m        = 2,
        trgg_Hp_m       = 3,
        trgg_gamma_deg  = 4,
        trgg_chi_deg    = 5,
        trgg_psi_deg    = 6,
        trgg_id_size    = 7
    };
    /**< enumeration that contains the different primary PID control loops */
    enum PID_PRIM_ID {
        pid_THR_vtas_mps  = 0,
        pid_ELV_theta_deg = 1,
        pid_AIL_xi_deg    = 2,
        pid_RUD_beta_deg  = 3,
        pid_prim_id_size  = 4
    };
    /**< enumeration that contains the different secondary PID control loops */
    enum PID_SECND_ID {
        pid_ELV_Hp_m          = 0,
        pid_ELV_h_m           = 1,
        pid_ELV_gammaTAS_deg  = 2,
        pid_AIL_chi_deg       = 3,
        pid_AIL_psi_deg       = 4,
        pid_AIL_muTAS_deg     = 5,
        pid_AIL_chiTAS_deg    = 6,
        pid_secnd_id_size     = 7
    };
    /**< Enumeration that contains the different motion types (implies whether navigation runs
     * or not and whether visual odometry takes input from objects or text files) */
    enum MOTION_TYPE {
        motion_default = 0,
        motion_dummy   = 1,
        motion_size    = 2
    };
    /**< Enumeration that contains the different visual and visual inertial algorithms */
    enum VISION_TYPE {
        vision_00       = 0,
        vision_01       = 1,
        vision_02       = 2,
        vision_03       = 3,
        vision_04       = 4,
        vision_05       = 5,
        vision_06       = 6,
        vision_CL       = 7,
        vision_99       = 8,
        vision_size     = 9
    };

} // closes namespace logic
} // closes namespace control

namespace sens {
namespace logic {
    /**< Enumeration that specifies the different terms considered in the sensors */
    enum SENS_COMPLETE_ID {
        sens_complete_noise             = 0,
        sens_complete_bias              = 1,
        sens_complete_bias_scale        = 2,
        sens_complete_bias_scale_par    = 3,
        sens_complete_size              = 4
    };
    /**< Enumeration that contains the random walk oscillation bands for inertial sensors */
    enum BAND_ID {
        band_id_base    = 0,
        band_id_bigger  = 1,
        band_id_biggest = 2,
        band_id_size    = 3
    };
    /**< Enumeration that contains the different accelerometer models */
    enum ACC_ID {
        acc_id_zero       = 0,
        acc_id_base       = 1,
        acc_id_better     = 2,
        acc_id_worse      = 3,
        acc_id_worst      = 4,
        acc_id53 = 5,
        acc_id_noise_only = 6,
        acc_id_base_old   = 7,
        acc_size          = 8
    };
    /**< Enumeration that contains the different gyroscope models */
    enum GYR_ID {
        gyr_id_zero       = 0,
        gyr_id_base       = 1,
        gyr_id_better     = 2,
        gyr_id_worse      = 3,
        gyr_id_worst      = 4,
        gyr_id53 = 5,
        gyr_id54 = 6,
        gyr_id_noise_only = 7,
        gyr_id_base_old   = 8,
        gyr_size          = 9
    };
    /**< Enumeration that contains the different magnetometer models */
    enum MAG_ID {
        mag_id_zero     = 0,
        mag_id_base     = 1,
        mag_id_better   = 2,
        mag_id_worse    = 3,
        mag_id_base_old = 4,
        mag_size        = 5
    };
    /**< Enumeration that contains the different outside static pressure sensor models */
    enum OSP_ID {
        osp_id_zero  = 0,
        osp_id_base  = 1,
        osp_id_worse = 2,
        osp_id_worst = 3,
        osp_size     = 4
    };
    /**< Enumeration that contains the different outside air temperature sensor models */
    enum OAT_ID {
        oat_id_zero  = 0,
        oat_id_base  = 1,
        oat_id_worse = 2,
        oat_id_worst = 3,
        oat_size     = 4
    };
    /**< Enumeration that contains the different true airspeed sensor models */
    enum TAS_ID {
        tas_id_zero   = 0,
        tas_id_base   = 1,
        tas_id_worse  = 2,
        tas_id_better = 3,
        tas_size      = 4
    };
    /**< Enumeration that contains the different angle of attack sensor models */
    enum AOA_ID {
        aoa_id_zero   = 0,
        aoa_id_base   = 1,
        aoa_id_worse  = 2,
        aoa_id_better = 3,
        aoa_size      = 4
    };
    /**< Enumeration that contains the different angle of sideslip sensor models */
    enum AOS_ID {
        aos_id_zero   = 0,
        aos_id_base   = 1,
        aos_id_worse  = 2,
        aos_id_better = 3,
        aos_size      = 4
    };
    /**< Enumeration that contains the different GPS receiver models */
    enum GPS_ID {
        gps_id_zero = 0,
        gps_id_base = 1,
        gps_size    = 2
    };

} // closes namespace logic
} // closes namespace sens

namespace st {
namespace logic {
    /**< Enumeration that contains the different initial Euler angle error generators */
    enum INITEUL_ID {
        initeul_id_zero   = 0,
        initeul_id_base   = 1,
        initeul_id_better = 2,
        initeul_id_worse  = 3,
        initeul_id_worst  = 4,
        initeul_size      = 5
    };

    /**< Enumeration that contains the different initial accelerometer error generators */
    enum INITACC_ID {
        initacc_id_zero  = 0,
        initacc_id_base  = 1,
        initacc_id_worse = 2,
        initacc_id_worst = 3,
        initacc_size     = 4
    };

    /**< Enumeration that contains the different initial gyroscope error generators */
    enum INITGYR_ID {
        initgyr_id_zero  = 0,
        initgyr_id_base  = 1,
        initgyr_id_worse = 2,
        initgyr_id_worst = 3,
        initgyr_size     = 4
    };

    /**< Enumeration that contains the different initial magnetometer error generators */
    enum INITMAG_ID {
        initmag_id_zero  = 0,
        initmag_id_base  = 1,
        initmag_id_worse = 2,
        initmag_id_worst = 3,
        initmag_size     = 4
    };

    /**< Enumeration that contains the different initial NED magnetic field error generators */
    enum INITMGN_ID {
        initmgn_id_zero  = 0,
        initmgn_id_base  = 1,
        initmgn_id_worse = 2,
        initmgn_id_worst = 3,
        initmgn_size     = 4
    };
    /**< enumeration that contains the different initial height error generators */
    enum INITHGH_ID {
        inithgh_id_zero  = 0,
        inithgh_id_base  = 1,
        inithgh_id_worse = 2,
        inithgh_size     = 3
    };
    /**< enumeration that contains the different starting points */
    enum STI_ID {
        sti_h3000_tas25_psi00 = 0,  // h [m] 3000.0, tas [mps] 25.0, psi [deg]   0.0
        sti_h3000_tas30_psi00 = 1,  // h [m] 3000.0, tas [mps] 30.0, psi [deg]   0.0
        sti_h3000_tas30_psi30 = 2,  // h [m] 3000.0, tas [mps] 30.0, psi [deg]  30.0
        sti_Hp3000_tas30_psi00 = 3,  // Hp [m] 3000.0, tas [mps] 30.0, psi [deg]  0.0
        sti_Hp3000_tas30_psi90 = 4,  // Hp [m] 3000.0, tas [mps] 30.0, psi [deg] 90.0
        sti_Hp3000_tas30_chi00 = 5,  // Hp [m] 3000.0, tas [mps] 30.0, chi [deg]  0.0
        sti_Hp3000_tas30_chi90 = 6,  // Hp [m] 3000.0, tas [mps] 30.0, chi [deg] 90.0
        sti_Hp3000_tas30_chiXX = 7,  // Hp [m] 3000.0, tas [mps] 30.0, chi [deg] from guidance
        sti_Hp3000_tas30_psiXX = 8,  // Hp [m] 3000.0, tas [mps] 30.0, psi [deg] from guidance
        sti_HpXX_tasXX_chiXX = 9,  // Hp [m], tas [mps], and chi [deg] from guidance
        sti_HpXX_tasXX_psiXX = 10, // Hp [m], tas [mps], and psi [deg] from guidance
        sti_rozas = 11,
        sti_size = 12
    };
} // closes namespace logic
} // closes namespace st

namespace vis {
    namespace logic {

    /**< Enumeration that contains the different frame handling stages. */
    enum HANDLING_STAGE {
        handling_stage_first_frame     = 0,
        handling_stage_second_frame    = 1,
        handling_stage_normal_frame    = 2,
        handling_stage_relocalizing    = 3,
        handling_stage_size            = 4
    };

    } // closes namespace logic
} // closes namespace vis

#endif
