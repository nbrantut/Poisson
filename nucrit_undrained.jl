include("exportfig.jl")
include("dem.jl")
using PyPlot

z0=[0 0.001 0.01 0.05 0.1]
Nz0=length(z0)

Nnu=250
tabnu = linspace(0,0.495,Nnu)
Nalpha=250
tabalpha = logspace(-3,2,Nalpha)
dnu=zeros(Nnu,Nalpha,Nz0)
for j=1:Nnu, i=1:Nalpha, k=1:Nz0
    dnu[i,j,k] = (1-z0[k])*Pu(tabnu[j],tabalpha[i],z0[k]) - Q(tabnu[j],tabalpha[i])
end

nuc_crack(a,z) = 40.*z./(81pi.*a)
nuc_sphere(z) = 0.2*(1+4z)
nuc_needle(z) = (7-sqrt(29))/8 + (203+36*sqrt(29)).*z/522

fig = figure()
ax = gca()


contour(tabalpha,tabnu,dnu[:,:,1]',[0], colors=[0.6 0.6 0.6],linewidths=1)

for k=2:Nz0
    semilogx(logspace(-3,-1,50), nuc_crack(logspace(-3,-1,50),z0[k]), "k--", linewidth=1)
    #plot(10.^[-.5;.5], nuc_sphere(z0[k])*[1;1], "k--", linewidth=1)
    plot(1, nuc_sphere(z0[k]), "k.", linewidth=1)
    plot([5;100], nuc_needle(z0[k])*[1;1], "k--", linewidth=1)
    contour(tabalpha,tabnu,dnu[:,:,k]',[0], colors="k",linewidths=1)
end

plot([1;1],[0.47,0.5],"k-",linewidth=0.5)

ylim(0,0.5)
xlabel("aspect ratio, \$\\alpha\$",usetex="true")
ylabel("critical Poisson's ratio, \$\\nu_{0,\\mathrm{crit}}\$",usetex="true")

ax2 = ax[:twinx]()
ax2[:get_yaxis]()[:set_ticks_position]("right")
ax2[:get_yaxis]()[:set_ticks](poisson([1.5, 1.6, 1.7, 1.8, 2.0, 2.4, 3.0]))
ax2[:get_yaxis]()[:set_ticklabels](("1.5","1.6", "1.7", "1.8", "2.0", "2.4", "3.0"))
ax2[:get_yaxis]()[:set_label_position]("right")
ax2[:set_ylabel]("critical \$V_\\mathrm{P}/V_\\mathrm{S}\$ ratio", usetex="true")
ax2[:set_ylim](ax[:get_ylim]())
ax2[:set_visible]("true")

#annotate("undrained limit",usetex="true",
#         xy=[100, 0],
#         ha="right",
#         va="bottom")

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

annotate(L"\zeta=0.1",
         xy=[0.1,0.3],
         ha="left",
         va="bottom",
         usetex="true")
annotate(L"0.05",
         xy=[0.04,0.25],
         ha="left",
         va="bottom",
         usetex="true")
annotate(L"0.01",
         xy=[0.01, 0.2],
         ha="left",
         va="bottom",
         usetex="true")
annotate(L"10^{-3}",
         xy=[0.0017, 0.1],
         ha="left",
         va="bottom",
         usetex="true")
annotate("dry (\$\\zeta=0\$)",
         xy=[0.0011, 0.002],
         usetex="true",
         ha="left",
         va="bottom",
         color=(0.6,0.6,0.6))

exportfig(fig, "nucrit_undrained", xsize=12, ysize=9)
