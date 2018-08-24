include("exportfig.jl")
include("dem.jl")
using PyPlot

z0=[0.001 0.01 0.1]
Nz0=length(z0)

Nnu=250
tabnu = linspace(0,0.5,Nnu)
Nalpha=250
tabalpha = logspace(-3,2,Nalpha)
dnu=zeros(Nnu,Nalpha,Nz0)
for j=1:Nnu, i=1:Nalpha, k=1:Nz0
    dnu[i,j,k] = (1-z0[k])*Pu(tabnu[j],tabalpha[i],z0[k]) - Qu(tabnu[j],tabalpha[i],z0[k])
end

nuc_crack(a,z) = 8 .*z./(27pi.*a)
nuc_sphere(z) = 0.2*(1+4z)
nuc_needle(z) = (7-sqrt(29))/8 + (551+91*sqrt(29)).*z/1392

fig=figure()

for k=1:Nz0
    semilogx(logspace(-3,-0.5,50), nuc_crack(logspace(-3,-0.5,50),z0[k]), "k--", linewidth=1)
    plot(1, nuc_sphere(z0[k]), "k.", linewidth=1)
    plot([5;100], nuc_needle(z0[k])*[1;1], "k--", linewidth=1)
    contour(tabalpha,tabnu,permutedims(dnu[:,:,k]),[0], colors="k",linewidths=1)
    plot([1;1],[0.47,0.5],"k-",linewidth=0.5)
end

ylim(0,0.5)
xlabel("aspect ratio, \$\\alpha\$",usetex="true")
ylabel("critical Poisson's ratio, \$\\nu_{0,\\mathrm{crit}}\$",usetex="true")

annotate("unrelaxed limit",usetex="true",
         xy=[100, 0],
         ha="right",
         va="bottom")

annotate("cracks",
         usetex="true",
         xy=[0.03,0.48],
         ha="center",
         va="center")


annotate("needles",
         usetex="true",
         xy=[10,0.48],
         ha="center",
         va="center")

exportfig(fig, "nucrit_unrelaxed")
