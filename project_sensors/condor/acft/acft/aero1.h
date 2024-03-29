#ifndef ACFT_AERO1
#define ACFT_AERO1

#include "../acft.h"
#include "aero.h"
#include "math/vec/vec1.h"
#include <string>

namespace math {
    class f_table2V;
}
namespace ang {
    class so3_tangent;
}
namespace acft {

// CLASS AERO1
// ===========
// ===========

class ACFT_API aero1 : public aero {
private:
    /**< name required internally */
    static const std::string _filename;
    /**< vector of angle of attack positions */
    math::vec1 _vec_alpha_rad;
    /**< vector of elevator positions */
    math::vec1 _vec_deltaE_deg;

    /**< amount to add to force coefficient along xBFS */
    double _cf_bfsi_extra;

    /**< force coefficient */
    math::f_table2V* _Pfun_cf_bfsi;
    math::f_table2V* _Pfun_cf_bfsii;
    math::f_table2V* _Pfun_cf_bfsiii;

    /**< force coefficient relationship with sideslip */
    math::f_table2V* _Pfun_dcf_bfsi_dbeta_rad;
    math::f_table2V* _Pfun_dcf_bfsii_dbeta_rad;
    math::f_table2V* _Pfun_dcf_bfsiii_dbeta_rad;

    /**< force coefficient relationship with wi_bfsned_bfs */
    math::f_table2V* _Pfun_dcf_bfsi_dwi_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsii_dwi_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsiii_dwi_bfsnedbfs_rps;

    /**< force coefficient relationship with wii_bfsned_bfs */
    math::f_table2V* _Pfun_dcf_bfsi_dwii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsii_dwii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsiii_dwii_bfsnedbfs_rps;

    /**< force coefficient relationship with wiii_bfsned_bfs */
    math::f_table2V* _Pfun_dcf_bfsi_dwiii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsii_dwiii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcf_bfsiii_dwiii_bfsnedbfs_rps;

    /**< force coefficient relationship with surfi_deg */
    math::f_table2V* _Pfun_dcf_bfsi_dsurfi_deg;
    math::f_table2V* _Pfun_dcf_bfsii_dsurfi_deg;
    math::f_table2V* _Pfun_dcf_bfsiii_dsurfi_deg;

    /**< force coefficient relationship with surfii_deg */
    math::f_table2V* _Pfun_dcf_bfsi_dsurfii_deg;
    math::f_table2V* _Pfun_dcf_bfsii_dsurfii_deg;
    math::f_table2V* _Pfun_dcf_bfsiii_dsurfii_deg;

    /**< force coefficient relationship with surfiii_deg */
    math::f_table2V* _Pfun_dcf_bfsi_dsurfiii_deg;
    math::f_table2V* _Pfun_dcf_bfsii_dsurfiii_deg;
    math::f_table2V* _Pfun_dcf_bfsiii_dsurfiii_deg;

    /**< moment coefficient */
    math::f_table2V* _Pfun_cm_bfsi;
    math::f_table2V* _Pfun_cm_bfsii;
    math::f_table2V* _Pfun_cm_bfsiii;

    /**< moment coefficient relationship with sideslip */
    math::f_table2V* _Pfun_dcm_bfsi_dbeta_rad;
    math::f_table2V* _Pfun_dcm_bfsii_dbeta_rad;
    math::f_table2V* _Pfun_dcm_bfsiii_dbeta_rad;

    /**< moment coefficient relationship with wi_bfsned_bfs */
    math::f_table2V* _Pfun_dcm_bfsi_dwi_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsii_dwi_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsiii_dwi_bfsnedbfs_rps;

    /**< moment coefficient relationship with wii_bfsned_bfs */
    math::f_table2V* _Pfun_dcm_bfsi_dwii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsii_dwii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsiii_dwii_bfsnedbfs_rps;

    /**< moment coefficient relationship with wiii_bfsned_bfs */
    math::f_table2V* _Pfun_dcm_bfsi_dwiii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsii_dwiii_bfsnedbfs_rps;
    math::f_table2V* _Pfun_dcm_bfsiii_dwiii_bfsnedbfs_rps;

    /**< moment coefficient relationship with surfi_deg */
    math::f_table2V* _Pfun_dcm_bfsi_dsurfi_deg;
    math::f_table2V* _Pfun_dcm_bfsii_dsurfi_deg;
    math::f_table2V* _Pfun_dcm_bfsiii_dsurfi_deg;

    /**< moment coefficient relationship with surfii_deg */
    math::f_table2V* _Pfun_dcm_bfsi_dsurfii_deg;
    math::f_table2V* _Pfun_dcm_bfsii_dsurfii_deg;
    math::f_table2V* _Pfun_dcm_bfsiii_dsurfii_deg;

    /**< moment coefficient relationship with surfiii_deg */
    math::f_table2V* _Pfun_dcm_bfsi_dsurfiii_deg;
    math::f_table2V* _Pfun_dcm_bfsii_dsurfiii_deg;
    math::f_table2V* _Pfun_dcm_bfsiii_dsurfiii_deg;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**< constructor based on interpolation method */
    aero1(math::logic::INTERP_MODE interp_mode = math::logic::biparabolic);
    /**< copy constructor */
    aero1(const aero1&) = delete;
    /**< move constructor */
    aero1(aero1&&) = delete;
    /**< destructor */
    ~aero1() override;
    /**< copy assignment */
    aero1& operator=(const aero1&) = delete;
    /**< move assignment */
    aero1& operator=(aero1&&) = delete;

    /**< computes the aerodynamic forces coefficient based on the true airspeed, the Euler angles
    from WFS to BFS (angles of attack and sideslip), the control parameters, and the angular
    velocity of BFS over NED expressed in BFS. */
    Eigen::Vector3d obtain_cf_aer(const double& vtas_mps, const ang::euler& euler_wfsbfs_rad, const Eigen::Array4d& delta_control, const ang::so3_tangent& w_nedbfsbfs_rps) const override;
    /**< computes the aerodynamics moments coefficient based on the true airspeed, the Euler angles
    from WFS to BFS (angles of attack and sideslip), the control parameters, and the angular
    velocity of BFS over NED expressed in BFS. */
    Eigen::Vector3d obtain_cm_aer(const double& vtas_mps, const ang::euler& euler_wfsbfs_rad, const Eigen::Array4d& delta_control, const ang::so3_tangent& w_nedbfsbfs_rps) const override;
}; // closes class aero1

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

} // closes namespace acft

#endif

