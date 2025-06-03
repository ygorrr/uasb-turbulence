using Gridap
n = 30

Ly=2
domain = (0,1,0,Ly)
partition = (n,Ly*n)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

# For convenience, we create two new boundary tags,  namely `"diri1"` and `"diri0"`, one for the top side of the square (where the velocity is non-zero), and another for the rest of the boundary (where the velocity is zero).

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[6,])
add_tag_from_tags!(labels,"botton",[5,])

# ## FE spaces
#
# For the velocities, we need to create a conventional vector-valued continuous Lagrangian FE space. In this example, we select a second order interpolation.

order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
# V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["diri00","botton"])#,"diri1"])
V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["botton"])#,"diri1"])

# The interpolation space for the pressure is built as follows

reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2,dirichlet_tags=["top"])#constraint=:zeromean)

reffec = ReferenceFE(lagrangian,Float64,order)
R = TestFESpace(model,reffec,conformity=:H1,labels=labels,dirichlet_tags=["botton"])#constraint=:zeromean)

reffek = ReferenceFE(lagrangian,Float64,order)
S = TestFESpace(model,reffec,conformity=:H1,labels=labels,dirichlet_tags=["botton"])#constraint=:zeromean)

reffeepsilon = ReferenceFE(lagrangian,Float64,order)
W = TestFESpace(model,reffec,conformity=:H1,labels=labels,dirichlet_tags=["botton"])#constraint=:zeromean)

# With the options `:Lagrangian`, `space=:P`, `valuetype=Float64`, and `order=order-1`, we select the local polynomial space $P_{k-1}(T)$ on the cells $T\in\mathcal{T}$. With the symbol `space=:P` we specifically chose a local Lagrangian interpolation of type "P". Without using `space=:P`, would lead to a local Lagrangian of type "Q" since this is the default for quadrilateral or hexahedral elements. On the other hand, `constraint=:zeromean` leads to a FE space, whose functions are constrained to have mean value equal to zero, which is just what we need for the pressure space. With these objects, we build the trial multi-field FE spaces

uD0 = VectorValue(0,0)
uD1 = VectorValue(0,1)
cD1 = 1.0

# U = TrialFESpace(V,[uD0,uD1])
g0(x,t::Real) = uD0 #[1] #1-2*x[1] 1-2*x[1] #
g0(t::Real) = x -> g0(x,t)
# g1(x,t::Real) = VectorValue(0,(1.0-x[2]/Ly)*(0.5*(1.0-(cos(5*2*pi*x[1])))^4+0.1*rand())+x[2]/Ly )
gtv(a,b) = (a>b) ? 1 : 0
# g1(x,t::Real) = VectorValue(0,10.0 * (1.0-x[2]/Ly)*(((0.5*(1.0-(cos(3*2*pi*x[1])))+0.025*rand()) > 0.99 ? 1.0 : 0.0 ))+x[2]/Ly  )#[1] #1-2*x[1] 1-2*x[1] #
c0(x,t::Real) = 0.0 #[1] #1-2*x[1] 1-2*x[1] #
c0(t::Real) = x -> c0(x,t)


g1(x,t::Real) = VectorValue(0,10.0 * 
(1.0-sqrt(x[2]/Ly))*
(((0.5*(1.0-(cos(1*2*pi*x[1])))+0.0*rand()) > 0.99 ? 1.0 : 0.0 ))+
sqrt(x[2]/Ly)  )#[1] #1-2*x[1] 1-2*x[1] #
g1(t::Real) = x -> g1(x,t)
p1(t::Real) = x -> 0
c1(x,t::Real) = 1.0
c1(t::Real) = x -> c1(x,t)

k1(t::Real) = x -> 1.0
epsilon1(t::Real) = x -> 1.0

U = TransientTrialFESpace(V ,[g1])
P = TransientTrialFESpace(Q)
C = TransientTrialFESpace(R,[c1] )
K = TransientTrialFESpace(S,[k1] )
E = TransientTrialFESpace(W,[epsilon1] )


Y = MultiFieldFESpace([V, Q, R, S, W])
X = TransientMultiFieldFESpace([U, P, C, K, E])

# ## Triangulation and integration quadrature
#
# From the discrete model we can define the triangulation and integration measure

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

Γₕ = BoundaryTriangulation(model,tags="botton")
dΓₕ = Measure(Γₕ,degree)
n_Γₕ = get_normal_vector(Γₕ)


# ## Nonlinear weak form
#
# The different terms of the nonlinear weak form for this example are defined following an approach similar to the one discussed for the $p$-Laplacian equation, but this time using the notation for multi-field problems.
Sc=1.0
Re = 300.0
beta=0.1
wc=VectorValue(0.0,1.1)

conv(u,∇u) = (∇u')⋅u
# convc(u,∇c) = ∇c[1]* u[1] +  ∇c[2]* (u[2]. - 2)
buss(v,c) = beta * v[2]*c

wwc(u)= u - wc

convc(u,∇c) = (∇c)⋅u
# convc(u,∇c) =  (∇c')  ⊙  u 

# dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)
nuEddy(u3, u4) = 0*Cμ * u3 * u3 / max(u4, minVal)
kProduction(∇u1, u3, u4) = nuEddy(u3, u4) * (∇u1 + ∇u1') ⊙ ∇u1

# The bilinear form reads
a((u,p,c,k,epsilon),(v,q,r,s,w)) = 
∫( 1/Re*(∇(v)⊙∇(u)) - (∇⋅v)*p )dΩ + 
∫( q*(∇⋅u) )dΩ + 
∫(v ⋅ (conv∘(u,∇(u))))dΩ +
∫(buss∘(v,c))dΩ +
∫( 1/(Re*Sc)*(∇(r)⊙∇(c)) )dΩ +
∫(r ⋅ inner(wwc(u),∇(c)))dΩ +
∫( (1/Re)*(∇(s)⊙∇(k)) )dΩ +
∫(s ⋅ inner(u,∇(k)))dΩ +
0*∫( s * (kProduction∘(∇(u),c,epsilon)) )dΩ
∫( 1/Re*(∇(w)⊙∇(epsilon)) )dΩ +
∫(w ⋅ inner(u,∇(epsilon)))dΩ


at(t,(u,p,c,k,epsilon),(v,q,r,s,w)) = ∫( ∂t(u)⋅v)dΩ + ∫( ∂t(c)⋅r)dΩ + ∫( ∂t(k)⋅s)dΩ + ∫( ∂t(epsilon)⋅w)dΩ


# Finally, the Navier-Stokes weak form residual and Jacobian can be defined as
res(t,(u,p,c,k,epsilon),(v,q,r,s,w)) = at(t,(u,p,c,k,epsilon),(v,q,r,s,w)) + a((u,p,c,k,epsilon),(v,q,r,s,w)) # + c((u,p,c),(v,q,r)) ##+ cc(u,c,r)
# jac((u,p,c),(du,dp,dc),(v,q,r)) = a((du,dp,dc),(v,q,r)) #+ dc(u,du,v)
# jac_t(t,(u,p,c),(dut,dpt,dct),(v,q,r)) = a((dut,dpt,dct),(v,q,r)) #+ dc(u,du,v)

# res(t,u,v) = ∫( ∂t(u)*v + κ(t)*(∇(u)⋅∇(v)) - f(t)*v )dΩ
# jac(t,u,du,v) = ∫( κ(t)*(∇(du)⋅∇(v)) )dΩ
# jac_t(t,u,duₜ,v) = ∫( duₜ*v )dΩ
# op = TransientFEOperator(res,jac,jac_t,X,Y)



# With the functions `res`, and `jac` representing the weak residual and the Jacobian, we build the nonlinear FE problem:
# op = FEOperator(res,jac,X,Y)
op = TransientFEOperator(res,X,Y)

# ## Nonlinear solver phase
#
# To finally solve the problem, we consider the same nonlinear solver as previously considered for the  $p$-Laplacian equation.

using LineSearches: BackTracking
nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=20)

# Then, we define the ODE solver. That is, the scheme that will be used for the time integration. In this tutorial we use the `ThetaMethod` with $\theta = 0.5$, resulting in a 2nd order scheme. The `ThetaMethod` function receives the linear solver, the time step size $\Delta t$ (constant) and the value of $\theta $.
CFL=1.0/10
Δt = CFL*1/n
θ = 1

# ode_solver = ThetaMethod(linear_solver,Δt,θ)
ode_solver = ThetaMethod(nls,Δt,θ)


# Finally, we define the solution using the `solve` function, giving the ODE solver, the FE operator, an initial solution, an initial time and a final time. To construct the initial condition we interpolate the initial value (in that case a constant value of 0.0) into the FE space $U(t)$ at $t=0.0$.

U₀ = interpolate_everywhere([g1(0),p1(0),c0(0),c0(0),c0(0)],X(0.0))
t₀ = 0.0
T = 4.0
uₕₜ = solve(ode_solver,op,t₀,T,U₀)
it=0
uh, ph, ch, kh, epsilonh = U₀
#  pvd[t] = createvtk(Ω,"burgers2D_$t"*".vtu",cellfields=["u"=>uₕ])
writevtk(Ωₕ,(@__DIR__)*"/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph, "ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh]) #"omega"=>curl(uh)])

it=1
for (t,uₕ) in uₕₜ
  global it, vec1
  local uh,ph,ch,kh,epsilonh
  uh, ph, ch, kh, epsilonh = uₕ
  #  pvd[t] = createvtk(Ω,"burgers2D_$t"*".vtu",cellfields=["u"=>uₕ])
  # writevtk(Ωₕ,"results/tns-resultsC$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"omega"=>curl(uh)])
  if(mod(it,10)==0)
    writevtk(Ωₕ,(@__DIR__)*"/results/uasbcp$(div(it,10)).vtu",cellfields=["uh"=>uh,"ph"=>ph,"ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh]) #"omega"=>curl(uh)])
  end

  it=it+1
  display(it)
end