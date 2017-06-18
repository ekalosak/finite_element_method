# Basic application of the FEM to solve evolutions of the 1D acoustic wave
# equation u_tt - c^2 u_xx = 0 i.e. the space and time second derivatives are
# proportional. u(x,t) is the displacement of the wave at position x and time t.
# c^2 is the "wave speed". We can derive this relation as a limit of point
# masses connected by springs obeying Hooke's law. That is, the force from a
# string is f = kx, etc. on Wikipedia.

# As a beginner problem, I'll solve for u(x) when d2/dx2 u = f where
# f(x) = x/6 over [0,1]. The analytical solution is simply
# u = int int f dx dx = int x^2/2 + c dx = let c and cc = 0 so
# u = x^3 and hence it's easy to validate my solution.

### SOURCES:
# 1. Introduction: https://www.geophysik.uni-muenchen.de/~igel/Lectures/NMG/07_finite_elements_weq.pdf
# 2. Dynamic mesh: http://epubs.siam.org/doi/abs/10.1137/0903002
# 3. Basic introduction: https://www.geophysik.uni-muenchen.de/~igel/downloads/nmgiifembasic.pdf
