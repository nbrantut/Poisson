include("exportfig.jl")
include("dem.jl")
using PyPlot

#z0=[0.097 0.0175 0.073 0.049 0.18 0.27]
z0=[0.025 0.05 0.075 0.1 0.2]
Nz0=length(z0)

Nnu=250
tabnu = linspace(0.2,0.4,Nnu)
Nalpha=250
tabalpha = logspace(-2,0,Nalpha)
dnu=zeros(Nnu,Nalpha,Nz0)
for j=1:Nnu, i=1:Nalpha, k=1:Nz0
    dnu[i,j,k] = (1-z0[k])*Pu(tabnu[j],tabalpha[i],z0[k]) - Q(tabnu[j],tabalpha[i])
end

#nuc_crack(a,z) = 40.*z./(81pi.*a)
nuc_sphere(z) = 0.2*(1+4z)
#nuc_needle(z) = (7-sqrt(29))/8 + (203+36*sqrt(29)).*z/522

fig = figure()
ax = gca()

for k=1:Nz0
    contour(tabalpha,tabnu,dnu[:,:,k]',[0], colors="k",linewidths=1)
end

for k=1:Nz0-1
    annotate(string(z0[k]),
             xy=[1,nuc_sphere(z0[k])],
             ha="right",
             va="bottom")
end

annotate("\$\\zeta="*string(z0[Nz0])*"\$",
         usetex="true",
         xy=[1,nuc_sphere(z0[Nz0])],
         ha="right",
         va="bottom")

rec = Array{Any,1}(4)
#gyp dhy
rec[1] = matplotlib[:patches][:Rectangle](
    xy=[0.04,0.33],
    width=0.06,
    height=0.02,
    alpha=0.6,
    facecolor=(0.6,0.6,0.6),
    edgecolor=(0.3,0.3,0.3)
)
annotate("gyp",
         xy=[0.04,0.33],
         ha="left",
         va="bottom")
#liz dhy (1GPa)
rec[2] = matplotlib[:patches][:Rectangle](
    xy=[0.04,0.31],
    width=0.06,
    height=0.02,
    alpha=0.6,
    facecolor=(0.6,0.6,0.6),
    edgecolor=(0.3,0.3,0.3)
)
annotate("liz",
         xy=[0.04,0.31],
         ha="left",
         va="bottom")
#atg dhy (1GPa)
rec[3] = matplotlib[:patches][:Rectangle](
    xy=[0.04,0.27],
    width=0.06,
    height=0.02,
    alpha=0.6,
    facecolor=(0.6,0.6,0.6),
    edgecolor=(0.3,0.3,0.3)
)
annotate("atg",
         xy=[0.04,0.27],
         ha="left",
         va="bottom")
#melt (2GPa)
rec[4] = matplotlib[:patches][:Rectangle](
    xy=[0.1,0.25],
    width=0.4,
    height=0.05,
    alpha=0.6,
    facecolor=(0.6,0.6,0.6),
    edgecolor=(0.3,0.3,0.3)
)
annotate("melt",
         xy=[0.4,0.30],
         ha="right",
         va="top")

for k=1:4
    ax[:add_patch](rec[k])
end

ax[:set_xscale]("log")

ylim(0.2,0.4)
xlabel("aspect ratio, \$\\alpha\$",usetex="true")
ylabel("critical Poisson's ratio, \$\\nu_{0,\\mathrm{crit}}\$",usetex="true")

ax2 = ax[:twinx]()
ax2[:get_yaxis]()[:set_ticks_position]("right")
ax2[:get_yaxis]()[:set_ticks](poisson([1.6, 1.7, 1.8, 2.0, 2.4, 3.0]))
ax2[:get_yaxis]()[:set_ticklabels](("1.6", "1.7", "1.8", "2.0", "2.4", "3.0"))
ax2[:get_yaxis]()[:set_label_position]("right")
ax2[:set_ylabel]("critical \$V_\\mathrm{P}/V_\\mathrm{S}\$ ratio", usetex="true")
ax2[:set_ylim](ax[:get_ylim]())
ax2[:set_visible]("true")



exportfig(fig, "nucrit_examples", xsize=8.4, ysize=6.5)
