
Six motion representations (speu_dcm, speu_rodrigues, homogeneous, trfv, screw, dual)
- All include constructors, move constructors, assigments, and move assignments:
- Transformations                            --> * and /, inverse, power, sclerp, plus
- Specific transformations                   --> negative for rodrigues, factor * and factor / for rotv
- Operations                                 --> * and /, ^ and &, minus, log map (no trfv nor screw), exp map (only trfv and screw)
- Adjoint                                    --> | and % (no screw) --> transforms twists or motion velocities
- Twist or Motion Velocity - Time Derivative --> space and body, both directions (no trfv, no screw)
- Linear algebra                             --> to obtain underlying Eigen classes (homogeneous and trfv only)
- Hat and Wedge                              --> back and forth from vector form (no screwl)
- Jacobians                                  --> no screw, no dual

The se3_tangent, se3_tangent_homo, and se3_tangent_dual represent the twist or motion velocity.





