using LinearAlgebra
using ProgressMeter


struct HamiltonSystem
    d
    m
    F
    q0
    p0
    Minv
    sol
    U
    HamiltonSystem(d,m,F,q0,p0,Minv,sol,U) =
        length(q0)==d*length(m) && length(p0)==d*length(m) ?
        new(d,m,F,q0,p0,Minv,sol,U) : error("dimension mismatch")
end

function getMinv(m,d)
    m_long=zeros(d*length(m))
    for (i,n) in enumerate(m)
        m_long[(i-1)*d+1:(i)*d] = length(n)==d ? collect(n) : [n for i in 1:d]
    end
    return Diagonal(1 ./ m_long)
end

HamiltonSystem(d,m,q0,p0,F=missing,U=missing,sol=missing) = HamiltonSystem(d,m,F,q0,p0,getMinv(m,d),sol,U)


function getEnergy(HS,q,p)
    return 1/2*p'*HS.Minv*p + HS.U(q)
end

function getEnergyNorm(HS,Q,P)
    QA = [Q[:,i] for i in 1:size(Q,2)]
    PA = [P[:,i] for i in 1:size(P,2)]
    E = map(getEnergy,[HS for q in QA],QA,PA)
    E = E.-E[1]
    return map(abs,E)
end

getEnergyErrors(HS,Qs,Ps) = [getEnergyNorm(HS,Q,P) for (Q,P) in zip(Qs,Ps)]
#getMaximumErrors(HS,ts,Qs)= [[norm((Q-hcat(map(HS.sol,t)...))[:,i]) for i in 1:length(t)] for (t,Q) in zip(ts,Qs)]
function getMaximumErrors(HS,ts,Qs)
    ME = []
    for (t,Q) in zip(ts,Qs)
        QA = QA = [Q[:,i] for i in 1:size(Q,2)]
        E=map(norm,QA-map(HS.sol,t))
        push!(ME,E)
    end
    return ME
end

println("Finished loading HamiltonSystem")
