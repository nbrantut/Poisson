include("exportfig.jl")
include("dem.jl")
using PyPlot

Nnu=250
tabnu = linspace(0,0.25,Nnu)
Nalpha=250
tabalpha = logspace(-3,2,Nalpha)
dnu=zeros(Nnu,Nalpha)
for j=1:Nnu, i=1:Nalpha
    dnu[i,j] = P(tabnu[j],tabalpha[i]) - Q(tabnu[j],tabalpha[i])
end

nu_crack(a) = (4/(3pi)+5pi/36)*a + (-254/81+80/(27pi^2)+29pi^2/864)*a.^2 + (1228800 + 1660160pi^2+165504pi^4+315pi^6)/(186624pi^3).*a.^3
nu_sphere(a) = 0.2 - 16/875*(1-a).^2 - (1-a).^3 *(3017088*(5751377+23283*sqrt(59385)))/(42875*(135+sqrt(59385))^4)
nu_needle = (7-sqrt(29))/8

fig=figure()

semilogx(logspace(-3,-0.5,50), nu_crack(logspace(-3,-0.5,50)), "k--", linewidth=1)
plot(logspace(-0.7,0.3), nu_sphere(logspace(-0.7,0.3)), "k--", linewidth=1)
plot([5;100], nu_needle*[1;1], "k--", linewidth=1)
contour(tabalpha,tabnu,dnu',[0], colors="k",linewidths=1)
plot([1;1],[0.23,0.25],"k-",linewidth=0.5)

xlabel("aspect ratio, \$\\alpha\$",usetex="true")
ylabel("critical Poisson's ratio, \$\\nu_\\mathrm{fixed}\$",usetex="true")

annotate("dry limit\n\$\\zeta\\rightarrow 0\$",usetex="true",
         xy=[100, 0],
         ha="right",
         va="bottom")

annotate("cracks",
         usetex="true",
         xy=[0.03,0.24],
         ha="center",
         va="center")


annotate("needles",
         usetex="true",
         xy=[10,0.24],
         ha="center",
         va="center")

exportfig(fig, "dry")
