### POISSON FORMULATION
## - u_xx = f

# As a beginner problem, I'll solve for u(x) when d2/dx2 u = f where
# f(x) = x/6 over [0,1]. The analytical solution is simply
# u = int int f dx dx = int x^2/2 + c dx = let c and cc = 0 so
# u = x^3 and hence it's easy to validate my solution.

## Parameters
d = c(0,1) # target function domain
nx = 25 # number of collocating nodes
hs = rep(abs(d[1]-d[2])/(nx-1), nx-1) # distance between nodes, equally spaced

## Collocate nodes
xs = c(d[1], Reduce(sum, hs, accumulate=T))

## Define negative second derivative of the displacement function
f = function(x){
    return(x/6)
}

u_analytic = function(x){
    return(x^3)
}

## Using a simple set of basis functions i.e. triangles with peaks at target
## nodes with phi_i(x_i) = 1 and bases at neighboring nodes with phi_i(x_j) = 0
## for all j \ne i, the stiffness matrix is rather simple:
## Also, ignoring boundary conditions at first
A = diag(nx) # set diagonals to 1
diag_ix = c(1,  Reduce(sum, rep(dim(A)[1]+1,dim(A)[1]-1), accumulate=T) + 1)
A[diag_ix[2:length(diag_ix)-1]] = 2 # set non-corner diagonals to 2
sub_diag_ix = diag_ix[1:length(diag_ix)-1] + 1
A[sub_diag_ix] = -1 # the subdiagonal entries must be -1
sup_diag_ix = diag_ix[2:length(diag_ix)] - 1
A[sup_diag_ix] = -1 # above teh diagonal must be -1 as well
