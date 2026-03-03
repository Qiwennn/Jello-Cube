# 3D Soft-Body Physics Simulator

# Demo Video



# What’s Implemented

# 1) External force field evaluation
- Where: `physics.cpp` (`getForceField()`).
- What:
  - Maps 3D world coordinates to 3D grid indices for the external force field array.
  - Implements trilinear interpolation (computing weights `A000` to `A111`) to sample the exact external force vector at any particle's continuous position.
  - Includes robust boundary checks to clamp grid indices within the valid range of the force field.

# 2) Structural mass-spring network
- Where: `physics.cpp` (`computeAcceleration()`).
- What:
  - Iterates through the 8x8x8 jello control point grid to accumulate forces.
  - Computes structural spring forces using Hooke's Law (`kElastic`) and damping (`dElastic`).
  - Connects each particle to its 6 immediate orthogonal neighbors (up, down, left, right, front, back) with a rest length of `1/7`.

# 3) Shear and bend springs
- Where: `physics.cpp` (`computeAcceleration()`).
- What:
  - Implements shear springs across face diagonals (12 neighbors, rest length `sqrt(2)/7`) and body diagonals (4 neighbors, rest length `sqrt(3)/7`) to maintain jello volume.
  - Implements bend springs connecting every second particle (6 neighbors, rest length `2/7`) to resist folding.
  - Accumulates both elastic and damping forces for all spring types into a total internal force per particle.


# 4) Bounding box collision detection
- Where: `physics.cpp` (`computeAcceleration()`).
- What:
  - Performs boundary checks against a defined cubic bounding box (from `-2.0` to `2.0` on the X, Y, and Z axes).
  - Applies a penalty method upon collision: computes a restoring elastic force (`kCollision`) based on penetration depth and a damping force (`dCollision`) based on velocity to push particles back into bounds.

# 5) Acceleration computation
- Where: `physics.cpp` (`computeAcceleration()`).
- What:
  - Sums up all structural, shear, bend, bounding box collision, and external field forces for each control point.
  - Converts the accumulated total force into acceleration by dividing by the particle's mass (`1.0 / jello->mass`).

# 6) Explicit Euler integration
- Where: `physics.cpp` (`Euler()`).
- What:
  - Implements a first-order explicit numerical integration step.
  - Calls `computeAcceleration()` to get the current state's derivatives.
  - Updates positions and velocities directly by multiplying current velocities and accelerations by `dt`.

# 7) 4th-Order Runge-Kutta (RK4) integration
- Where: `physics.cpp` (`RK4()`).
- What:
  - Implements a highly stable fourth-order integration method to handle the stiff spring dynamics of the soft-body simulation.
  - Uses a buffer to evaluate derivatives at 4 intermediate states (`F1`, `F2`, `F3`, `F4`) across the timestep.
  - Combines the 4 states using the standard RK4 weighted sum `(1/6 * (F1 + 2*F2 + 2*F3 + F4))` to properly update the final positions and velocities, drastically reducing explosive instabilities.
