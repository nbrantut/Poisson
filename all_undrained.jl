include("exportfig.jl")
include("dem.jl")
using PyPlot

#take an array of initial Poisson's ratios
nu0 = [0.15, 0.2, 0.25, 0.3, 0.35]
Nn = length(nu0)
#an array of aspect ratios
a0 = [0.001, 0.1, 1.001] #note that we avoid alpha=zeta or alpha=1 for numerical reasons. this has no impact (negligible differences in the plots)
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
tol = 1e-6

#initialise solution arrays
NU  = Array{Array{Float64,1}}(Nn,Na,Nz)
#NUDRY = Array{Array{Float64,1}}(Nn,Na)
PHI = Array{Array{Float64,1}}(Nn,Na,Nz)
AX = Array{PyCall.PyObject,1}(Nn*Na*Nz)

dndp_crack(n0,a,z) = (20-34n0)./(45pi*a) + (1-n0)/3*(1-1./z) 
dndp_sphere(n0,z)  = 1.5*(1-5n0-n0^2+5n0^3)/(7-5n0) + 0.75*(1-n0-n0^2+n0^3)/(1-2n0)*z
dndp_needle(n0,z)  = (5-23n0-12n0^2+16n0^3)/15 + (25-15n0-24n0^2+16n0^3)/(27*(1-2n0))*z


for i in eachindex(nu0), j in eachindex(a0), k in eachindex(z0)
    
    #parameters:p=[r0,alpha,zeta]
    p = [GoverK(nu0[i]),a0[j],z0[k]]
    y,phi = ode_rk4(DEMKGdry,m0,phispan,tol,p)

    kg = KGassmann(phi, z0[k], getitem(y,1))

    #NUDRY[i,j] = Nu(getitem(y,1),getitem(y,2),nu0[i])
    NU[i,j,k]  = Nu(kg,getitem(y,2),nu0[i])
    PHI[i,j,k] = phi
    
end

fig = figure()
cla()

for k=1:Nz, j=1:Na

    ipl = Na*(k-1) + j
    
    AX[ipl] = subplot(Nz,Na,ipl)
    ax = AX[ipl]

    for i=1:Nn
        #plot(PHI[i,j,k]*100, NUDRY[i,j], color=0.8*[1,1,1],linewidth=1)
        plot(PHI[i,j,k]*100, NU[i,j,k], "k", linewidth=1)
        ylim(0,0.5)
    end

    #if we have cracks just show 1% range in porosity, and plot asymptote. else show all.
    if a0[j]<=0.01
        xlim(0,1)
        if a0[j]<z0[k]
            phirange= linspace(0,.25e-2)
            for i=1:Nn
                plot(phirange*100,
                     nu0[i] + phirange*dndp_crack(nu0[i],a0[j],z0[k]),
                     "k--",
                     linewidth=1)
            end
        end
    else
        xlim(phispan[1]*100,phispan[2]*100)
    end

    if a0[j]>=1
        phirange= linspace(0,50e-2)
        for i=1:Nn
            plot(phirange*100,
                 nu0[i] + phirange*dndp_sphere(nu0[i],z0[k]),
                 "k--",
                 linewidth=1)
        end
    end
    
    if rem(ipl,Na)!=1
        ax[:get_yaxis]()[:set_ticklabels]("")
    else
        ylabel("Poisson's ratio, \$\\nu\$", usetex="true")
    end

    if rem(ipl,Na)==0
        ax2 = ax[:twinx]()
        
        ax2[:get_yaxis]()[:set_ticks_position]("right")
        ax2[:get_yaxis]()[:set_ticks](poisson([1.5, 1.6, 1.7, 1.8, 2.0, 2.4, 3.0]))
        ax2[:get_yaxis]()[:set_ticklabels](("1.5","1.6", "1.7", "1.8", "2.0", "2.4", "3.0"))
        ax2[:get_yaxis]()[:set_label_position]("right")
        ax2[:set_ylabel]("\$V_\\mathrm{P}/V_\\mathrm{S}\$ ratio", usetex="true")
        ax2[:set_ylim](ax[:get_ylim]())
        ax2[:set_visible]("true")

    end
    

    if ipl<=Na*(Nz-1)
        ax[:get_xaxis]()[:set_ticklabels]("")
    else
        ax[:set_xlabel]("porosity, \$\\phi\$ (\\%)", usetex="true")
        if a0[j]<0.01
            ax[:get_xaxis]()[:set_ticklabels](("0.00","0.25", "0.50", "0.75", "1.0"))
        end
        
    end
    
    
end

#manual annotations...
sca(AX[1])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-3}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(a)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[2])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-3}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(b)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[3])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-3}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(c)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[4])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-2}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(d)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[5])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-2}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(e)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[6])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-2}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(f)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[7])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-3}\$\n\$\\zeta=10^{-1}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(g)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[8])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=10^{-1}\$\n\$\\zeta=10^{-1}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(h)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")
sca(AX[9])
xl = xlim()
yl = ylim()
annotate("\$\\alpha=1\$\n\$\\zeta=10^{-1}\$",
         xy=[xl[2],yl[1]],
         ha="right",
         va="bottom",
         usetex="true")
annotate("(i)",
         xy=[xl[1],yl[2]],
         ha="left",
         va="top")



exportfig(fig, "all_undrained",xsize=16.8,ysize=12.6)

