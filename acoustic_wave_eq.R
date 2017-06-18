# Basic application of the FEM to solve evolutions of the 1D acoustic wave
# equation u_tt - c^2 u_xx = 0 i.e. the space and time second derivatives are
# proportional. u(x,t) is the displacement of the wave at position x and time t.
# c^2 is the "wave speed". We can derive this relation as a limit of point
# masses connected by springs obeying Hooke's law. That is, the force from a
# string is f = kx, etc. on Wikipedia.

### SOURCES:
# 1. Introduction: https://www.geophysik.uni-muenchen.de/~igel/Lectures/NMG/07_finite_elements_weq.pdf
# 2. Dynamic mesh: http://epubs.siam.org/doi/abs/10.1137/0903002
# 3. Basic introduction: https://www.geophysik.uni-muenchen.de/~igel/downloads/nmgiifembasic.pdf
