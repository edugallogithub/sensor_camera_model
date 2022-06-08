#ifndef EAGLE_LOGIC
#define EAGLE_LOGIC

#include "eagle.h"

namespace eagle {

namespace logic {
    /**< Enumeration that contains the eagle image creation methods */
    enum EAGLE_TYPE {
        eagle_final    = 0,
        eagle_att      = 1,
        eagle_osg      = 2,
        eagle_blind    = 3,
        eagle_size     = 4
    };

    /**< Enumeration that contains the different converter options */
    enum CONVERTER_ID {
        converter_osg_vis_opencv_no         = 0,
        converter_osg_sav_opencv_no         = 1,
        converter_osg_sav_opencv_vis        = 2,
        converter_osg_sav_opencv_sav        = 3,
        converter_osg_sav_opencv_vis_stream = 4,
        converter_osg_sav_opencv_sav_stream = 5,
        converter_osg_str_opencv_no         = 6,
        converter_osg_str_opencv_vis        = 7,
        converter_osg_str_opencv_sav        = 8,
        converter_osg_savstr_opencv_vis     = 9,
        converter_osg_savstr_opencv_sav     = 10,
        converter_size                      = 11
    };
    /**< Enumeration that contains the different image file formats */

    enum IMAGE_ID {
        image_jpg  = 0,
        image_png  = 1,
        image_bmp  = 2,
        image_size = 3
    };

} // closes namespace logic

} // closes namespace vis

#endif








