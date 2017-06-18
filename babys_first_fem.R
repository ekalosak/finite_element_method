rm(list=ls())
set.seed(33)

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
K = 100 # resolution of plotting xs

## Collocate nodes
xs = c(d[1], Reduce(sum, hs, accumulate=T))
xs[nx] = d[2]
pltxs = seq(d[1], d[2], length.out=K)

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
h = hs[1]
A = (1/h)*A # this step relies on the fact that the nodes are equally spaced

## Generate basis functions
# Here, these are the simple triangles described above
phi1 = function(x){
    if(x < xs[1]){
        return(0)
    }else if(x >= xs[2]){
        return(0)
    }else{
        return(-(x-xs[2])/h)
    }
}

phin = function(x){
    if(x <= xs[nx-1]){
        return(0)
    }else if(x > xs[nx]){
        return(0)
    }else{
        return((x-xs[nx-1])/h)
    }
}

phi = function(x, i){
    if(x < xs[i-1]){
        return(0)
    }else if(x > xs[i+1]){
        return(0)
    }else if(x < xs[i]){
        return((x-xs[i-1])/h)
    }else{
        return(1 - (x-xs[i])/h)
    }
}
phiv = Vectorize(phi, "x")

## Solve b*A = g
fs = f(xs)
b = solve(A) %*% fs # these are the coefficients of the basis functions

## Regress the FEM solution
ys = rep(0, K)
for(a in 1:K){
    # for each plotting x value
    for(i in 1:nx){
        # for each basis function
        # calculate the value of the function at that point and add its weighted
        # result to the ys vector
        if(i == 1){
            ys[a] = ys[a] + b[1]*phi1(pltxs[a])
        }else if(i == nx){
            ys[a] = ys[a] + b[nx]*phin(pltxs[a])
        }else{
            ys[a] = ys[a] + b[i]*phi(pltxs[a], i)
        }
    }
}

## Plot a couple basis functions
# plot(pltxs, apply(matrix(pltxs), 1, phi1))
plt_phi_edges = ggplot() +
    geom_line(data=data.frame(x=pltxs, y=apply(matrix(pltxs), 1, phin)),
              aes(x=x, y=y), color="steelblue") +
    geom_line(data=data.frame(x=pltxs, y=apply(matrix(pltxs), 1, phi1)),
              aes(x=x, y=y), color="darkviolet")

