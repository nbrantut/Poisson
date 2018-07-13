function spheroidgeometry(alpha)
    if (alpha < 1)
        g = alpha * (acos(alpha) - alpha * sqrt(0.1e1 - alpha ^ 2)) * (0.1e1 - alpha ^ 2) ^ (-0.3e1 / 0.2e1)
    end
    if (alpha > 1)
        g = alpha * (alpha * sqrt((alpha ^ 2 - 1)) - acosh(alpha)) * ((alpha ^ 2 - 1) ^ (-0.3e1 / 0.2e1))
    end 
    if (alpha == 1.0)
        g=2/3
    end
    return g
end


function Pu(nu,alpha,Zeta)
    g = spheroidgeometry(alpha)
    return ((12 * nu * alpha ^ 2 - 24 * alpha ^ 2 - 12 * nu - 3) * g - 4 * nu * alpha ^ 2 + 14 * alpha ^ 2 + 4 * nu + 4) * (1 - nu) / ((12 * nu ^ 2 - 12 * nu ^ 2 * alpha ^ 2 - 6 * nu * alpha ^ 2 + 12 * Zeta * nu ^ 2 * alpha ^ 2 - 12 * Zeta * nu ^ 2 + 6 * Zeta * nu * alpha ^ 2 - 6 * Zeta * alpha ^ 2 - 6 + 6 * alpha ^ 2 + 6 * nu + 6 * Zeta - 6 * Zeta * nu) * g ^ 2 + (-9 * Zeta - 12 * Zeta * nu * alpha ^ 2 - 12 * nu + 6 + 48 * nu * alpha ^ 2 + 12 * Zeta * nu ^ 2 - 24 * alpha ^ 2 - 12 * Zeta * nu ^ 2 * alpha ^ 2 + 3 * Zeta * nu) * g + 4 * Zeta * nu ^ 2 * alpha ^ 2 - 4 * Zeta * nu ^ 2 - 24 * nu * alpha ^ 2 + 4 * Zeta + 2 * Zeta * alpha ^ 2 + 12 * alpha ^ 2 + 6 * Zeta * nu * alpha ^ 2)
end

function Qu(nu,alpha,Zeta)
    g = spheroidgeometry(alpha)
    return  0.4e1 / 0.5e1 * (-1 + nu) * (-32 * Zeta * nu ^ 2 + 544 * alpha ^ 4 * nu - 256 * nu ^ 2 * alpha ^ 4 + 88 * Zeta * nu ^ 2 * alpha ^ 2 + 128 * alpha ^ 6 * nu ^ 2 - 352 * alpha ^ 6 * nu + 16 * alpha ^ 6 * Zeta + 40 * Zeta * alpha ^ 4 + 96 * Zeta * nu ^ 3 * alpha ^ 4 - 136 * Zeta * alpha ^ 4 * nu - 32 * Zeta * nu ^ 3 * alpha ^ 6 + 96 * Zeta * nu * alpha ^ 2 + 24 * alpha ^ 6 * Zeta * nu ^ 2 - 88 * Zeta * alpha ^ 2 + 32 * Zeta + 144 * alpha ^ 6 - 80 * Zeta * nu ^ 2 * alpha ^ 4 + 72 * alpha ^ 6 * Zeta * nu + (80 + 52 * Zeta * nu ^ 2 - 1800 * alpha ^ 4 * nu + 1056 * nu ^ 2 * alpha ^ 4 - 210 * Zeta * nu ^ 2 * alpha ^ 2 - 424 * alpha ^ 6 * nu ^ 2 + 1012 * alpha ^ 6 * nu - 2 * alpha ^ 6 * Zeta - 168 * Zeta * alpha ^ 4 - 312 * Zeta * nu ^ 3 * alpha ^ 4 + 408 * Zeta * alpha ^ 4 * nu + 104 * Zeta * nu ^ 3 * alpha ^ 6 - 96 * nu ^ 3 * alpha ^ 4 - 300 * Zeta * nu * alpha ^ 2 - 106 * alpha ^ 6 * Zeta * nu ^ 2 + 96 * nu ^ 3 * alpha ^ 2 + 222 * Zeta * alpha ^ 2 - 52 * Zeta - 32 * nu ^ 3 + 208 * nu ^ 2 + 32 * nu ^ 3 * alpha ^ 6 - 404 * alpha ^ 6 + 264 * Zeta * nu ^ 2 * alpha ^ 4 - 212 * alpha ^ 6 * Zeta * nu + 648 * alpha ^ 4 - 256 * nu - 324 * alpha ^ 2 + 104 * Zeta * nu - 840 * nu ^ 2 * alpha ^ 2 - 104 * Zeta * nu ^ 3 + 1044 * nu * alpha ^ 2 + 312 * Zeta * nu ^ 3 * alpha ^ 2) * g - 208 * alpha ^ 4 + 64 * alpha ^ 2 - 32 * Zeta * nu + (-6 - 36 * Zeta * nu ^ 2 - 162 * alpha ^ 4 * nu - 108 * nu ^ 2 * alpha ^ 4 + 72 * alpha ^ 6 * nu ^ 2 + 72 * alpha ^ 6 * nu + 48 * alpha ^ 6 * Zeta - 90 * Zeta * alpha ^ 4 - 144 * Zeta * nu ^ 3 * alpha ^ 4 + 162 * Zeta * alpha ^ 4 * nu + 48 * Zeta * nu ^ 3 * alpha ^ 6 + 144 * nu ^ 3 * alpha ^ 4 - 108 * Zeta * nu * alpha ^ 2 - 72 * alpha ^ 6 * Zeta * nu ^ 2 - 144 * nu ^ 3 * alpha ^ 2 + 36 * Zeta * alpha ^ 2 + 6 * Zeta + 48 * nu ^ 3 + 36 * nu ^ 2 - 48 * nu ^ 3 * alpha ^ 6 - 48 * alpha ^ 6 + 108 * Zeta * nu ^ 2 * alpha ^ 4 - 72 * alpha ^ 6 * Zeta * nu + 90 * alpha ^ 4 - 18 * nu - 36 * alpha ^ 2 + 18 * Zeta * nu - 48 * Zeta * nu ^ 3 + 108 * nu * alpha ^ 2 + 144 * Zeta * nu ^ 3 * alpha ^ 2) * g ^ 3 + 128 * nu ^ 2 * alpha ^ 2 + 32 * Zeta * nu ^ 3 - 192 * nu * alpha ^ 2 + (-70 - 9 * Zeta * nu ^ 2 + 1026 * alpha ^ 4 * nu - 252 * nu ^ 2 * alpha ^ 4 + 198 * Zeta * nu ^ 2 * alpha ^ 2 + 72 * alpha ^ 6 * nu ^ 2 - 588 * alpha ^ 6 * nu - 88 * alpha ^ 6 * Zeta + 291 * Zeta * alpha ^ 4 + 408 * Zeta * nu ^ 3 * alpha ^ 4 - 486 * Zeta * alpha ^ 4 * nu - 136 * Zeta * nu ^ 3 * alpha ^ 6 - 192 * nu ^ 3 * alpha ^ 4 + 378 * Zeta * nu * alpha ^ 2 + 180 * alpha ^ 6 * Zeta * nu ^ 2 + 192 * nu ^ 3 * alpha ^ 2 - 228 * Zeta * alpha ^ 2 + 25 * Zeta - 64 * nu ^ 3 - 108 * nu ^ 2 + 64 * nu ^ 3 * alpha ^ 6 + 268 * alpha ^ 6 - 369 * Zeta * nu ^ 2 * alpha ^ 4 + 228 * alpha ^ 6 * Zeta * nu - 426 * alpha ^ 4 + 210 * nu + 228 * alpha ^ 2 - 120 * Zeta * nu + 288 * nu ^ 2 * alpha ^ 2 + 136 * Zeta * nu ^ 3 - 648 * nu * alpha ^ 2 - 408 * Zeta * nu ^ 3 * alpha ^ 2) * g ^ 2 - 96 * Zeta * nu ^ 3 * alpha ^ 2) / ((alpha ^ 2 - nu + nu * alpha ^ 2 + 2) * g - 2 * alpha ^ 2) / ((7 - 4 * alpha ^ 2 - 8 * nu + 8 * nu * alpha ^ 2) * g - 8 + 6 * alpha ^ 2 + 8 * nu - 8 * nu * alpha ^ 2) / ((12 * nu ^ 2 - 12 * nu ^ 2 * alpha ^ 2 - 6 * nu * alpha ^ 2 + 12 * Zeta * nu ^ 2 * alpha ^ 2 - 12 * Zeta * nu ^ 2 + 6 * Zeta * nu * alpha ^ 2 - 6 * Zeta * alpha ^ 2 - 6 + 6 * alpha ^ 2 + 6 * nu + 6 * Zeta - 6 * Zeta * nu) * g ^ 2 + (-9 * Zeta - 12 * Zeta * nu * alpha ^ 2 - 12 * nu + 6 + 48 * nu * alpha ^ 2 + 12 * Zeta * nu ^ 2 - 24 * alpha ^ 2 - 12 * Zeta * nu ^ 2 * alpha ^ 2 + 3 * Zeta * nu) * g + 4 * Zeta * nu ^ 2 * alpha ^ 2 - 4 * Zeta * nu ^ 2 - 24 * nu * alpha ^ 2 + 4 * Zeta + 2 * Zeta * alpha ^ 2 + 12 * alpha ^ 2 + 6 * Zeta * nu * alpha ^ 2)
end

function P(nu,alpha)
    return Pu(nu,alpha,0)
end

function Q(nu,alpha)
    return Qu(nu,alpha,0)
end

function vpvs(nu)
    return sqrt((2nu-1)/(2nu-1))
end

function GoverK(nu)
    return 3*(1-2nu)/(2*(1+nu))
end

function Nu(k,g,nu0)
    return (3.0.*k - 2.0.*GoverK(nu0).*g)./(6.0.*k + 2.0*GoverK(nu0).*g)
end

function KGassmann(phi,z,Kd)
    return Kd*(phi*(1-1/z) + 1-1/Kd)/(phi*(1-1/z)+Kd-1)
end


function DEMKGdry(dm,m,p,phi)
    #k=m[1]
    #g=m[2]
    #r0=G0overK0=p[1]
    #alpha=p[2]

    nu = (3.0*m[1]-2.0*p[1]*m[2])/(6.0*m[1]+2.0*p[1]*m[2])
    dm[1] = -P(nu,p[2])*m[1]/(1-phi)
    dm[2] = -Q(nu,p[2])*m[2]/(1-phi)
end

function DEMnudry(nu,a,phi)
    dnudphi = (1.0/(1.0-phi))*(1.0+nu)*(1.0-2.0*nu)/3.0 * (Q(nu,a)-P(nu,a))
    return dnudphi
end

function DEMKGunrelaxed(dm,m,p,phi)
    #k=m[1]
    #g=m[2]
    #r0=G0overK0=p[1]
    #alpha=p[2]
    #zeta=p[3]

    nu = (3.0*m[1]-2.0*p[1]*m[2])/(6.0*m[1]+2.0*p[1]*m[2])
    dm[1] = -Pu(nu,p[2],p[3]/m[1])*m[1]*(1.0-p[3]/m[1])/(1-phi)
    dm[2] = -Qu(nu,p[2],p[3]/m[1])*m[2]/(1-phi)
end

function getitem(A::Array{Array{Float64,1},1},i::Int)
    a = Array{typeof(A[1][i])}(length(A))
    for (ind,el) in enumerate(A)
        a[ind] = el[i]
    end
    return a
end
        

function ode_rk4(func,u0,tspan,tol,p)
    t = Float64[]  
    y = Array{typeof(u0)}(0)

    #initialise
    push!(t,tspan[1])
    push!(y,u0)

    #initialise bunch of shit
    k1=zeros(u0)
    k2=zeros(u0)
    k3=zeros(u0)
    k4=zeros(u0)

    dt = (tspan[2]-tspan[1])/100
    n=length(t)
    #loop
    while t[n]<tspan[2]

        n = length(t)

        #get RK steps using dt
        y0 = y[n] + rk4step(k1, k2, k3, k4, func, p, y[n], t[n], dt)

        #get RK steps using dt/2
        dthalf = 0.5*dt
        yhalf = y[n] + rk4step(k1, k2, k3, k4, func, p, y[n], t[n], dthalf)
        y1 = yhalf + rk4step(k1, k2, k3, k4, func, p, yhalf, t[n]+dthalf, dthalf)
        
        #check if in tolerance
        err = norm(y0-y1)
        if (err>tol) | (err<tol/10)
            dt = dt*exp(0.2*log(tol/err))
            y0 = y[n] + rk4step(k1, k2, k3, k4, func, p, y[n], t[n], dt)
        end
        push!(t,t[n]+dt)
        push!(y,y0)
    end
    
    return y,t
end

function rk4step(k1, k2, k3, k4, fun, p, y, t, dt)
    fun(k1,y,p,t)
    fun(k2,y+dt*(k1/2.0),p,t+(dt/2.0))
    fun(k3,y+dt*(k2/2.0),p,t+(dt/2.0))
    fun(k4,y+dt*k3,p,t+dt)

    return dt*((k1/6)+(k2/3)+(k3/3)+(k4/6))
end
