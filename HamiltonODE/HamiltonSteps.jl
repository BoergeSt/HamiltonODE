function euler_step(q,p,dt,HS)
    Q = q+dt*HS.Minv*p
    P = p+dt*HS.F(q)
    return (dt,Q,P)
end

function i_euler_step(q,p,dt,HS;tol = 1e-6,max_iter=100,no_warning=false)
    Q=q;P=p
    for i in 1:max_iter
        Qn = q+dt*HS.Minv*P
        Pn = p+dt*HS.F(Q)
        if norm(Qn-Q)<=tol && norm(Pn-P)<=tol
            return (dt,Qn,Pn)
        end
        Q=Qn
        P=Pn
    end
    if !no_warning
        println("WARNING: Implicity euler step did not reach tol.")
    end
    return (dt,Q,P)
end

function theta_step(q,p,dt,HS;tol = 1e-6,max_iter=Int(1e6),theta=0.5)
    _,Qe,Pe = euler_step(q,p,dt,HS)
    _,Qi,Pi = i_euler_step(q,p,dt,HS,tol = tol,max_iter = max_iter)
    Q,P = theta.*(Qe,Pe).+(1-theta).*(Qi,Pi)
    return (dt,Q,P)
end

function s_euler_step(q,p,dt,HS)
    P = p+dt*HS.F(q)
    Q = q+dt*HS.Minv*P
    return (dt,Q,P)
end

function verlet_step(q,p,dt,HS)
    Q = q+dt*HS.Minv*p+dt^2/2*HS.Minv*HS.F(q)
    P = p+dt/2*(HS.F(q)+HS.F(Q))
    return (dt,Q,P)
end

function rk4_step(q,p,dt,HS)
    k1q = dt*HS.Minv*p
    k1p = dt*HS.F(q)

    k2q = dt*HS.Minv*(p+k1p/2)
    k2p = dt*HS.F(q+k1q/2)

    k3q = dt*HS.Minv*(p+k2p/2)
    k3p = dt*HS.F(q+k2q/2)

    k4q = dt*HS.Minv*(p+k3p)
    k4p = dt*HS.F(q+k3q)

    Q = q + 1/6*(k1q + 2*k2q + 2*k3q + k4q)
    P = p + 1/6*(k1p + 2*k2p + 2*k3p + k4p)
    return (dt,Q,P)
end



function projected_step(q,p,dt,HS,tol=1e-5,sigma=0.5,max_iter=100,step=rk4_step)
    dt,Q,P = step(q,p,dt,HS)
    E0 = getEnergy(HS,HS.q0,HS.p0)
    E=getEnergy(HS,Q,P)-E0
    i=0;
    while abs(E)>tol
        if i > max_iter
            println("WARNING: projected rk4 did not reach tol")
            return (dt,Q,P)
        end
        dq = HS.F(Q)*E
        dp = -HS.Minv*P*E
        s=1/sigma
        diff = 0
        E_new = E
        j = 0
        while diff>=0
            if j>=max_iter
                println("WARNING: could not find sigma small enough!")
                return (dt,Q,P)
            end
            s*=sigma
            E_new = getEnergy(HS,(Q+s*dq),(P+s*dp))-E0
            diff = E_new^2-E^2
            j+=1
        end
        E=E_new
        Q=Q+s*dq
        P=P+s*dp
        i+=1
    end
    return (dt,Q,P)
end

println("Finished loading HamiltonSteps")
