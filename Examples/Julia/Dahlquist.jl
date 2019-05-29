include("../../HamiltonODE.jl")

using Plots
pyplot()

d = 2
m = (1)
q0 = [1.,0.]
p0 = [0.,1.]

F(q)= -q
U(q)= q⋅q
sol(q)=[cos(q),sin(q)]

HS = HamiltonSystem(d,m,q0,p0,F,U,sol)

T=2*pi
dt = 1e-1

t,Q,P=integrate(HS,Integrator(verlet_step,dt,T));


p1 = plot(t,Q')
p2 = plot(Q[1,:],Q[2,:])
plot(p1,p2)
