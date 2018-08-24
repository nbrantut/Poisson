include("exportfig.jl")
include("dem.jl")
using PyPlot

#take an array of initial Poisson's ratios
nu0 = [0.15, 0.2, 0.25, 0.3, 0.35]
Nn = length(nu0)
#an array of aspect ratios
a0 = [0.00101, 0.1, 1.001] #note that we avoid alpha=zeta or alpha=1 for numerical reasons. this has no impact (negligible differences in the plots)
Na = length(a0)
#an array of compressibility ratios
z0 = [0.001, 0.01, 0.1]
Nz = length(z0)

#Solve DEM equations
#initial conditions
m0 = [1.0;1.0]
#porosity range
phispan = (0.0,0.95)
#tolerance for integration
tol = 1e-9

#initialise solution arrays
NU  = Array{Array{Float64,1}}(undef,Nn,Na,Nz)
PHI = Array{Array{Float64,1}}(undef,Nn,Na,Nz)

for i in eachindex(nu0), j in eachindex(a0), k in eachindex(z0)
    
    #parameters:p=[r0,alpha,zeta]
    p = [GoverK(nu0[i]),a0[j],z0[k]]
    y,phi = ode_rk4(DEMKGunrelaxed,m0,phispan,tol,p)

    NU[i,j,k]  = Nu(getitem(y,1),getitem(y,2),nu0[i])
    PHI[i,j,k] = phi
    
end

fig = figure()
cla()

for k=1:Nz, j=1:Na

    ipl = Na*(k-1) + j
    
    ax = subplot(Nz,Na,ipl)

    for i=1:Nn
        plot(PHI[i,j,k]*100, NU[i,j,k], "k", linewidth=1)
        ylim(0,0.5)
    end

    #if we have cracks just show 1% range in porosity, else show all.
    if a0[j]<=0.01
        xlim(0,1)
    else
        xlim(phispan[1]*100,phispan[2]*100)
    end
    
    if rem(ipl,Na)!=1
        ax[:get_yaxis]()[:set_ticklabels]("")
    else
        ylabel("Poisson's ratio, \$\\nu\$", usetex="true")
    end
    

    if ipl<=Na*(Nz-1)
        ax[:get_xaxis]()[:set_ticklabels]("")
    else
        xlabel("porosity, \$\\phi\$ (\\%)", usetex="true")
    end
    
    
end

#manual annotations...
subplot(3,3,1)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-3}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
         usetex="true")
subplot(3,3,2)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-3}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
         usetex="true")
subplot(3,3,3)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-3}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
             usetex="true")
subplot(3,3,4)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-2}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
         usetex="true")
subplot(3,3,5)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-2}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
             usetex="true")
subplot(3,3,6)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-2}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
         usetex="true")
subplot(3,3,7)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-1}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
             usetex="true")
subplot(3,3,8)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-1}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
         usetex="true")
subplot(3,3,9)
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-1}\$",
             xy=[xl[2],yl[1]],
             ha="right",
             va="bottom",
             usetex="true")

exportfig(fig, "all_unrelaxed",xsize=15,ysize=12)

