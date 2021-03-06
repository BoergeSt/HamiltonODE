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
        println("WARNING: Implicit euler step did not reach tol.")
    end
    return (dt,Q,P)
end

function midpoint_step(q,p,dt,HS)
    Q = q+dt*HS.Minv*p
    P = p+dt*HS.F(q)
    Q = q+dt*HS.Minv*((p+P)/2)
    P = p+dt*HS.F((q+Q)/2)
    return (dt,Q,P)
end

function i_midpoint_step(q,p,dt,HS;tol = 1e-6,max_iter=100,no_warning=false)
    Q=q;P=p
    for i in 1:max_iter
        Qn = q+dt*HS.Minv*((p+P)/2)
        Pn = p+dt*HS.F((q+Q)/2)
        if norm(Qn-Q)<=tol && norm(Pn-P)<=tol
            return (dt,Qn,Pn)
        end
        Q=Qn
        P=Pn
    end
    if !no_warning
        println("WARNING: Implicit midpoint step did not reach tol.")
    end
    return (dt,Q,P)
end


function theta_step(q,p,dt,HS;theta=0.5, kwargs...)
    _,Qe,Pe = euler_step(q,p,dt,HS)
    _,Qi,Pi = i_euler_step(q,p,dt,HS;kwargs...)
    Q,P = theta.*(Qe,Pe).+(1-theta).*(Qi,Pi)
    return (dt,Q,P)
end

function s_euler_step(q,p,dt,HS)
    P = p+dt*HS.F(q)
    Q = q+dt*HS.Minv*P
    return (dt,Q,P)
end

function sa_euler_step(q,p,dt,HS)
    Q = q+dt*HS.Minv*p
    P = p+dt*HS.F(Q)
    return (dt,Q,P)
end

function verlet_step(q,p,dt,HS)
    Q = q+dt*HS.Minv*p+dt^2/2*HS.Minv*HS.F(q)
    P = p+dt/2*(HS.F(q)+HS.F(Q))
    return (dt,Q,P)
end

function q_split_verlet_step(q,p,dt,HS)
    Q = q + dt/2*HS.Minv*p
    P = p + dt*(HS.F(Q))
    Q = Q + dt/2*HS.Minv*P
    return (dt,Q,P)
end

function p_split_verlet_step(q,p,dt,HS)
    P = p + dt/2*(HS.F(q))
    Q = q + dt*HS.Minv*P
    P = P + dt/2*(HS.F(Q))
    return (dt,Q,P)
end


function puls_projected_step(q,p,dt,HS;step=verlet_step,kwargs...)
    (dt,Q,P) = verlet_step(q,p,dt,HS;kwargs...);
    E = getEnergy(HS,HS.q0,HS.p0)
    if P'*HS.Minv*P<0.001 || U(Q)>E
        gamma = 1
    else
        gamma = sqrt(2*(E-HS.U(Q))/(P'*HS.Minv*P))
    end
    return (dt,Q,gamma*P)
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



function projected_step(q,p,dt,HS;tol=1e-5,sigma=0.5,max_iter=100,step=rk4_step,no_warning=false,kwargs...)
    dt,Q,P = step(q,p,dt,HS,kwargs...)
    E0 = getEnergy(HS,HS.q0,HS.p0)
    E=getEnergy(HS,Q,P)-E0
    i=0;
    while abs(E)>tol
        if i > max_iter
            if !no_warning
                println("WARNING: projected step did not reach tol")
            end
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
                if !no_warning
                    println("WARNING: projected step could not find sigma small enough!")
                end
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
