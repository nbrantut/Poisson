include("exportfig.jl")
include("dem.jl")
using PyPlot

#take an array of initial Poisson's ratios
nu0 = [0.15, 0.2, 0.25, 0.3, 0.35]
#an array of aspect ratios
a0 = [0.001, 0.1, 1.0-1e-12]
#an array of compressibility ratios
z0 = [0.001, 0.01, 0.1]

#Solve DEM equations
#initial conditions
m0 = [1.0;1.0]
#porosity range
phispan = (0.0,0.3)
#tolerance for integration
tol = 1e-8

#initialise solution arrays
NU  = Array{Array{Float64,1}}(length(nu0),length(a0),length(z0))
PHI = Array{Array{Float64,1}}(length(nu0),length(a0),length(z0))

for i in eachindex(nu0), j in eachindex(a0), k in eachindex(z0)
    
    #parameters:p=[r0,alpha,zeta]
    p = [GoverK(nu0[i]),a0[j],z0[k]]
    y,phi = ode_rk4(DEMKGunrelaxed,m0,phispan,tol,p)

    NU[i,j,k]  = Nu(getitem(y,1),getitem(y,2),nu0[i])
    PHI[i,j,k] = phi
    
end
