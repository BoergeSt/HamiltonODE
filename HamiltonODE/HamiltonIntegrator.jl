using Printf

struct Integrator
    step
    dt
    T
    const_dt
    kwargs
end
Integrator(step::Function,dt,T,const_dt=false;kwargs...) =
    Integrator(step,dt,T,const_dt,kwargs)


function integrate(HS,step::Function,dt,T;kwargs...)
    t=[.0]
    Q=[HS.q0]
    P=[HS.p0]
    N = Int(ceil(T/dt/100))
    progress = Progress(N, desc="Integrating using $(typeof(step).name.mt.name) dt = $dt...")
    lastupdate=0
    while(t[end]<T)
        DT,q,p=step(Q[end],P[end],dt,HS;kwargs...)
        push!(Q,q)
        push!(P,p)
        push!(t,t[end]+DT)
        if t[end]-lastupdate>T/N
            lastupdate=t[end]
            ProgressMeter.next!(progress)
        end
    end
    ProgressMeter.finish!(progress)
    return (t,hcat(Q...),hcat(P...))
end

function integrate_constant_dt(HS,step::Function,dt,T;kwargs...)
    N = Int(ceil(T/dt))
    Q=zeros(length(HS.q0),N+1)
    P=zeros(length(HS.q0),N+1)
    t=[dt*i for i in 0:N]
    Q[:,1]=HS.q0
    P[:,1]=HS.p0
    progress = Progress(N, desc="Integrating using $(typeof(step).name.mt.name) dt = $dt...")
    lastupdate=0
    for i in 1:N
        DT,q,p=step(Q[:,i],P[:,i],dt,HS;kwargs...)
        Q[:,i+1]=q
        P[:,i+1]=p
        ProgressMeter.next!(progress)
    end
    ProgressMeter.finish!(progress)
    return (t,Q,P)
end

integrate(HS,I::Integrator) = I.const_dt ?
    integrate_constant_dt(HS,I.step,I.dt,I.T;I.kwargs...) :
    integrate(HS,I.step,I.dt,I.T;I.kwargs...)


function integrate(HS,Is)
    ts = []
    Qs = []
    Ps = []
    for I in Is
        t,Q,P = integrate(HS,I)
        push!(ts,t)
        push!(Qs,Q)
        push!(Ps,P)
    end
    return (ts,Qs,Ps)
end

function getNames(Is)
     return ["$(typeof(I.step).name.mt.name) with dt=$(@sprintf("%.2e", I.dt))$( length(I.kwargs)>0 ? "\n$(join(["$k: $(I.kwargs[k])" for k in keys(I.kwargs)],", "))" : "")" for I in Is]
end

println("Finished loading HamiltonIntegrator")
