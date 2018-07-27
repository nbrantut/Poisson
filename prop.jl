function V2moduli(vp,vs,rho)
    G = rho*vs^2
    K = rho*(vp^2 - 4*vs^2/3)
    nu = (3K-2G)/(6K+2G)
    return K,G,nu
end



function meltcompress(rho, rho0, KT0, KTp)
    P = 1.5*KT0*((rho/rho0).^(7/3) - (rho/rho0).^(5/3)).*(1 - 0.75(4-KTp)*((rho/rho0).^(2/3)-1))
    K = 1.5*KT0*(
        ((7/3)*(rho/rho0).^(4/3) - (5/3)*(rho/rho0).^(2/3)).*(1 - 0.75(4-KTp)*((rho/rho0).^(2/3)-1)) +
        ((rho/rho0).^(7/3) - (rho/rho0).^(5/3)).*(0.5*KTp*(rho/rho0).^(-1/3))
    )
    return P,K
end

#list of rock VpVs from Christensen 1996

println("Rock properties at 1 GPa")
println("")
println("AND | (K,G,nu)="*sprint(showcompact,V2moduli(5.940,3.177,2.627)))
println("BAS | (K,G,nu)="*sprint(showcompact,V2moduli(6.118,3.291,2.882)))
println("DIA | (K,G,nu)="*sprint(showcompact,V2moduli(6.814,3.766,2.936)))
println("GRA | (K,G,nu)="*sprint(showcompact,V2moduli(6.372,3.726,2.652)))
println("DIO | (K,G,nu)="*sprint(showcompact,V2moduli(6.675,3.756,2.810)))
println("GAB | (K,G,nu)="*sprint(showcompact,V2moduli(7.299,3.929,2.968)))
println("MGW | (K,G,nu)="*sprint(showcompact,V2moduli(6.139,3.512,2.682)))
println("SLT | (K,G,nu)="*sprint(showcompact,V2moduli(6.379,3.432,2.807)))
println("PHY | (K,G,nu)="*sprint(showcompact,V2moduli(6.398,3.608,2.738)))
println("BZE | (K,G,nu)="*sprint(showcompact,V2moduli(6.530,3.493,2.915)))
println("BPP | (K,G,nu)="*sprint(showcompact,V2moduli(6.571,3.623,2.835)))
println("BGR | (K,G,nu)="*sprint(showcompact,V2moduli(6.983,3.955,2.978)))
println("GGN | (K,G,nu)="*sprint(showcompact,V2moduli(6.271,3.627,2.643)))
println("BGN | (K,G,nu)="*sprint(showcompact,V2moduli(6.366,3.636,2.742)))
println("QSC | (K,G,nu)="*sprint(showcompact,V2moduli(6.523,3.654,2.824)))
println("AMP | (K,G,nu)="*sprint(showcompact,V2moduli(7.046,3.987,2.996)))
println("FGR | (K,G,nu)="*sprint(showcompact,V2moduli(6.571,3.667,2.758)))
println("PGR | (K,G,nu)="*sprint(showcompact,V2moduli(6.497,3.658,2.761)))
println("AGR | (K,G,nu)="*sprint(showcompact,V2moduli(7.114,3.810,2.763)))
println("MGR | (K,G,nu)="*sprint(showcompact,V2moduli(7,3.849,2.971)))
println("GGR | (K,G,nu)="*sprint(showcompact,V2moduli(7.324,4.052,3.111)))
println("ECL | (K,G,nu)="*sprint(showcompact,V2moduli(8.198,4.594,3.485)))
println("SER | (K,G,nu)="*sprint(showcompact,V2moduli(5.607,2.646,2.566)))
println("QTZ | (K,G,nu)="*sprint(showcompact,V2moduli(6.091,4.054,2.652)))
println("MBL | (K,G,nu)="*sprint(showcompact,V2moduli(6.985,3.794,2.721)))
println("ANO | (K,G,nu)="*sprint(showcompact,V2moduli(7.124,3.717,2.730)))
println("HBL | (K,G,nu)="*sprint(showcompact,V2moduli(7.317,4.160,3.248)))
println("PYX | (K,G,nu)="*sprint(showcompact,V2moduli(7.935,4.519,3.267)))
println("DUN | (K,G,nu)="*sprint(showcompact,V2moduli(8.399,4.783,3.31)))

