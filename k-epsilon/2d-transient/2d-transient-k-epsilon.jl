using Gridap

n = 30
Lx = 1.0
Ly = 2.0
domain = (0,Lx,0,Ly)
partition = (n,Ly*n)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[6,])
add_tag_from_tags!(labels,"bottom",[5,])

#=
Incógnitas:
u1 -> vetor velocidade média
u2 -> pressão média
U3 -> concentração média de partículas
u4 -> k
u5 -> epsilon

Espaços de funções de teste: V1, V2, V3, V4, V5
Espaços de de funções de ensaio: S1, S2, S3, S4, S5
=#

order = 2
reffe_u1 = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V1 = TestFESpace(model, reffe_u1, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

reffe_u2 = ReferenceFE(lagrangian, Float64, order-1; space=:P)
V2 = TestFESpace(model, reffe_u2, conformity=:L2, dirichlet_tags=["top"])

reffe_u3 = ReferenceFE(lagrangian, Float64, order)
V3 = TestFESpace(model, reffe_u3, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

reffe_u4 = ReferenceFE(lagrangian, Float64, order)
V4 = TestFESpace(model, reffe_u4, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

reffe_u5 = ReferenceFE(lagrangian, Float64, order)
V5 = TestFESpace(model, reffe_u5, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

# -----------------------
# Condições de contorno
# -----------------------

function u1BC(x, t::Real)
  global Lx, Ly
  ω = 2.0
  treshold = 0.89
  velocityAmplitude = 1.0
  velocityY = 0.0
  
  if -cos(2π*ω*x[1]/Lx) >= treshold
    velocityY = velocityAmplitude * (1 - sqrt(x[2]/Ly))
  end

  return VectorValue(0.0, velocityY)
end

u1BC(t::Real) = x -> u1BC(x,t)
# u1IC(x, t::Real) = u1BC(x,t)

# g1(x,t::Real) = VectorValue(0, 10.0 * (1.0-sqrt(x[2]/Ly)) * (
#     (0.5 * (1.0 - (cos(1*2*pi*x[1]))) + 0.0*rand() ) > 0.99 ? 1.0 : 0.0
#   ) + sqrt(x[2]/Ly))

u2BC(x, t::Real) = 0.0
u2BC(t::Real) = x -> u2BC(x,t)
u2IC(x, t::Real) = 0.0
u2IC(t::Real) = x -> u2IC(x,t)

u3BC(x, t::Real) = 1.0
u3BC(t::Real) = x -> u3BC(x,t)
u3IC(x, t::Real) = 0.0
u3IC(t::Real) = x -> u3IC(x,t)

# Condições de Dirichlet para k e epsilon
Uref = 1.0  # Velocidade de referência
I = 0.05      # Intensidade de turbulência
ell = 0.05    # Comprimento de mistura (tamanho do jato)
u4BC(x, t::Real) = 500 #1.5 * (I * Uref)^2
u4BC(t::Real) = x -> u4BC(x,t)
u4IC(x, t::Real) = 0.0
u4IC(t::Real) = x -> u4IC(x,t)

u5BC(x, t::Real) = 1200 #(Cμ)^(3/4) * k_inlet^(3/2) / ell
u5BC(t::Real) = x -> u5BC(x,t)
u5IC(x, t::Real) = 0.0
u5IC(t::Real) = x -> u5IC(x,t)

S1 = TransientTrialFESpace(V1, [u1BC])
S2 = TransientTrialFESpace(V2, [u2BC])
S3 = TransientTrialFESpace(V3, [u3BC])
S4 = TransientTrialFESpace(V4, [u4BC])
S5 = TransientTrialFESpace(V5, [u5BC])

Y = MultiFieldFESpace([V1, V2, V3, V4, V5])
X = TransientMultiFieldFESpace([S1, S2, S3, S4, S5])

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

wc = VectorValue(0.0, -1.1)
Us = VectorValue(0.0, -1.1) # Velocidade de decantação de partículas

β = 0.1
g = VectorValue(0.0, -1.0)  # Vetor aceleração unitário
nu = 1.0e-3  # Viscosidade cinemática
ρ = 1.0  # Densidade do fluido
Sc = 1.0
Re = 30.0
# Ri = β*g*(Cref/Lref)/(Uref/Lref)^2
Ri = 1e-1

Cμ = 0.09  # Constante de turbulência
Cϵ1 = 1.44  # Constante de produção de epsilon
Cϵ2 = 1.92  # Constante de destruição de epsilon
σk = 1.0  # Constante de difusão de k
σϵ = 1.3  # Constante de difusão de epsilon
σ = 1.0

minVal = 1e-8 # Constante de proteção contra divisão por zero
nuT(u4, u5) = Cμ * u4 * u4 / max(u5, minVal)
kProduction(∇u1, u4, u5) = nuEddy(u4, u5) * ((∇u1 + ∇u1') ⊙ ∇u1)
epsilonProduction(∇u1, u4, u5) = Cϵ1 * kProduction(∇u1, u4, u5) * u5 / max(u4, minVal)
# epsilonProduction(∇u1, u4, u5) = Cϵ1 * (nuEddy(u4, u5) * ((∇u1 + ∇u1') ⊙ ∇u1)) * u5 / max(u4, minVal)
epsilonDestruction(u4, u5) = Cϵ2 * u5 * u5 / max(u4, minVal)

# -----------------------
# Resíduos
# -----------------------

# Equação de Navier-Stokes 
resNS(t, u1, u2, u3, u4, u5, v1) = 
  ∫( v1 ⋅ ∂t(u1) )dΩ +
  ∫( v1 ⋅ (∇(u1)' ⋅ u1) )dΩ +
  ∫( (1/Re + 0*(nuT∘(u4,u5))) * ∇(u1)⊙∇(v1) )dΩ - 
  ∫( (∇ ⋅ v1) * u2 )dΩ -
  ∫( Ri * u3 * (v1 ⋅ g) )dΩ

# Equação da continuidade
resCont(u1, v2) =
  ∫( v2 * (∇ ⋅ u1) )dΩ

# Equação de transporte de partículas
resC(t, u1, u3, u4, u5, v3) =
  ∫( v3 * ∂t(u3) )dΩ +
  ∫( v3 ⋅ inner(∇(u3),(u1 + 0*Us)) )dΩ +
  ∫( (1/Re * 1/Sc + 0*(nuT∘(u4,u5))/σ) * ∇(u3)⊙∇(v3) )dΩ

# Equação de k
resk(t, u1, u3, u4, u5, v4) = 
  ∫( v4 * ∂t(u4) )dΩ +
  ∫( v4 ⋅ (∇(u4)' ⋅ u1) )dΩ +
  ∫( (nuT∘(u4,u5)) / σk * ∇(u4)⊙∇(v4) )dΩ -
  ∫( v4 * (kProduction∘(∇(u1),u4,u5)) )dΩ +
  ∫( v4 * u4 )dΩ +
  ∫( v4 * Ri * (nuT∘(u4,u5))/σ * g ⋅ ∇(u3) )dΩ

# Equação de epsilon
resEpsilon(t, u1, u4, u5, v5) =
  ∫( v5 * ∂t(u5) )dΩ +
  ∫( v5 ⋅ (∇(u5)' ⋅ u1) )dΩ +
  ∫( (nuT∘(u4,u5)) / σϵ * ∇(u5)⊙∇(v5) )dΩ -
  ∫( v5 * (epsilonProduction∘(∇(u1),u4,u5)) )dΩ +
  ∫( v5 * (epsilonDestruction∘(u4,u5)) )dΩ

res(t, (u1, u2, u3, u4, u5), (v1, v2, v3, v4, v5)) =
  resNS(t, u1, u2, u3, u4, u5, v1) +
  resCont(u1, v2) +
  0*resk(t, u1, u3, u4, u5, v4) +
  0*resEpsilon(t, u1, u4, u5, v5)

# conv(u,∇u) = (∇u')⋅u
# buss(v,c) = beta * v[2]*c

# wwc(u)= u - wc

# convc(u,∇c) = (∇c)⋅u

# a((u,p,c,k,epsilon),(v,q,r,s,w)) = 
# ∫( 1/Re*(∇(v)⊙∇(u)) - (∇⋅v)*p )dΩ + 
# ∫( q*(∇⋅u) )dΩ + 
# ∫(v ⋅ (conv∘(u,∇(u))))dΩ +
# ∫(buss∘(v,c))dΩ +
# ∫( 1/(Re*Sc)*(∇(r)⊙∇(c)) )dΩ +
# ∫(r ⋅ inner(wwc(u),∇(c)))dΩ +
# ∫( (1/Re)*(∇(s)⊙∇(k)) )dΩ +
# ∫(s ⋅ inner(u,∇(k)))dΩ +
# # 0*∫( s * (kProduction∘(∇(u),c,epsilon)) )dΩ +
# ∫( 1/Re*(∇(w)⊙∇(epsilon)) )dΩ +
# ∫(w ⋅ inner(u,∇(epsilon)))dΩ


# at(t,(u,p,c,k,epsilon),(v,q,r,s,w)) = ∫( ∂t(u)⋅v)dΩ + ∫( ∂t(c)⋅r)dΩ + ∫( ∂t(k)⋅s)dΩ + ∫( ∂t(epsilon)⋅w)dΩ

# res(t,(u,p,c,k,epsilon),(v,q,r,s,w)) = at(t,(u,p,c,k,epsilon),(v,q,r,s,w)) + a((u,p,c,k,epsilon),(v,q,r,s,w)) 

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

U₀ = interpolate_everywhere([u1BC(0),u2IC(0),u3IC(0),u4IC(0),u5IC(0)],X(0.0))
t₀ = 0.0
T = 4.0
uₕₜ = solve(ode_solver,op,t₀,T,U₀)
it=0
uh, ph, ch, kh, epsilonh = U₀

if !isdir((@__DIR__)*"/results")
  mkdir((@__DIR__)*"/results")
end

#  pvd[t] = createvtk(Ω,"burgers2D_$t"*".vtu",cellfields=["u"=>uₕ])
writevtk(Ω,(@__DIR__)*"/results/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph, "ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh]) #"omega"=>curl(uh)])

it=1
for (t,uₕ) in uₕₜ
  global it, vec1
  local uh,ph,ch,kh,epsilonh
  uh, ph, ch, kh, epsilonh = uₕ
  #  pvd[t] = createvtk(Ω,"burgers2D_$t"*".vtu",cellfields=["u"=>uₕ])
  # writevtk(Ωₕ,"results/tns-resultsC$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"omega"=>curl(uh)])
  writevtk(Ω,(@__DIR__)*"/results/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh]) #"omega"=>curl(uh)])
  # if(mod(it,1)==0)
  # end

  it=it+1
  display(it)
end