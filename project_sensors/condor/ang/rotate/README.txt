
Four rotation representations (dcm, rodrigues, rotv, euler)
- All include constructors, move constructors, assigments, and move assignments:
- Transformations                    --> * and /, inverse, power, slerp, plus (euler none)
- Specific transformations           --> negative for rodrigues, factor * and factor / for rotv
- Operations                         --> * and /, minus (no euler), log map (no rotv, euler), exp map (only rotv)
- Adjoint                            --> | and % (euler none) --> transforms angular velocities
- Angular Velocity - Time Derivative --> space and body, both directions
- Linear algebra                     --> to obtain underlying Eigen classes (euler none)
- Hat and Wedge                      --> back and forth from vector form (euler none)
- Obtain Individual Euler angles     --> euler also set, rotv none
- Jacobians                          --> euler none
- Specific euler only methods

The so3_tangent and so3_tangent_skew represent the angular velocity.



