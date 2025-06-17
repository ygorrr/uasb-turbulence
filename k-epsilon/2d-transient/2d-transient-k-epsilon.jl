using Gridap
using LineSearches: BackTracking
using Plots

uasbDiam = 1.8 # m
uasbHeight = 2.0 # m
uasbArea = π * (uasbDiam/2)^2 # m²
HDT = 6*60*60 # 6 horas em segundos
volFlowRate = uasbArea * uasbHeight / hdt # m³/s
meanVelocity = volFlowRate / uasbArea # m/s
molecularViscosity = 1.0e-3 # viscosidade cinemática do fluido (água) em m²/s

# Número de Reynolds do UASB
Re_uasb = meanVelocity * uasbDiam / molecularViscosity 

jetDiameter = 0.05 # m
jetArea = π * (jetDiameter/2)^2 # m²
meanJetVelocity = volFlowRate / jetArea # m/s

# Número de Reynolds do jato
Re_jet = meanJetVelocity * jetDiameter / molecularViscosity 

n = 30
Lx = 10.0
Ly = 2.0 * Lx
domain = (-Lx/2, Lx/2, 0, Ly)
partition = (n, 2*n)
model = CartesianDiscreteModel(domain, partition; isperiodic=(true,false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "top", [6,])
add_tag_from_tags!(labels, "bottom", [5,])

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

#=-----------------------
  Condições iniciais e de contorno
  -----------------------=#

function singleJet(x)
  global jetDiameter
  local a = 75.0
  local xTreshold = log(1999)/(2*a)
  local amplitude = 1.0
  Δx = -jetDiameter/2 + xTreshold
  return amplitude * (tanh(a*(x-Δx)) * -tanh(a*(x+Δx)) + 1) / 2
end

function jetProfile(x, t::Real)
  global xJets
  local velocityY = 0.0
  for xJet in xJets
    velocityY += singleJet(x - xJet)
  end
  return velocityY
end

# Visualização do perfil de velocidade do jato
begin
  x = -1:0.001:1.0
  xJets = [-1, -0.5, 0.0, 0.5, 1, 0.25]
  jetContour = xJets[3] .+ [- jetDiameter / 2, jetDiameter / 2]
  
  plot(x, jetProfile.(x,0), xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile", ylims=(-0.01,1.0))
  
  scatter!(jetContour, [0.0, 0.0], label="jetContour")
end

u1BC(x, t::Real) = jetProfile(x,t)
u1BC(t::Real) = x -> u1BC(x,t)
u1IC(x, t::Real) = exp(-5.0 * x[2] / Ly) * jetProfile(x,t)

u2BC(x, t::Real) = 0.0
u2BC(t::Real) = x -> u2BC(x,t)
u2IC(x, t::Real) = 0.0
u2IC(t::Real) = x -> u2IC(x,t)

u3BC(x, t::Real) = jetProfile(x,t)
u3BC(t::Real) = x -> u3BC(x,t)
u3IC(x, t::Real) = 0.0
u3IC(t::Real) = x -> u3IC(x,t)

#=
  Velocidade de referência, vel. média do jato
=#
# nJets = length(xJets)
# Ujet = volFlowRate / nJets / jetArea
# Uref = Ujet / Uchar
Uref = 1.0

#=
  Intensidade de turbulência. Valores típicos:
  Jato livre: 0.05 ~ 0.1
  Escoamento altamente turbulento: até 0.2
=#
I = 0.1 

#= 
  Estimativa de k na saída do jato
=#
kEstimate = (3/2) * (I * Uref)^2

u4BC(x, t::Real) = kEstimate * jetProfile(x,t)
u4BC(t::Real) = x -> u4BC(x,t)
u4IC(x, t::Real) = 0.01 * kEstimate
u4IC(t::Real) = x -> u4IC(x,t)

#= 
  Estimativa dimensional de epsilon na saída do jato.
  Supõe-se que u' ~ √k e que νₜ ~ u' * ℓ. Usando a relação entre k e ϵ:
    ϵ = Cμ^(3/4) * k^(3/2) / ℓ
  onde Cμ é uma constante de turbulência e ℓ é uma escala de comprimento. O expoente de Cμ é empírico e não vem diretamente da teoria.
  ℓ pode ser estimada por uma fração do diâmetro do jato:
    ℓ = α * φ
  onde α é uma constante, tipicamente entre 0.07 e 0.1.
=#
α = 0.07
ϵEstimate = Cμ^(3/4) * kEstimate^(3/2) / (α * jetDiameter)

#= 
  Estimativa baseada em argumentos de escala:
    νₜ ~ Ujet * φ
  Essa estimativa vem de argumentos de escala, supondo que as flutuações da velocidade e comprimento de mistura sejam proporcionais à velocidade média do jato e ao diâmetro do jato, ou seja, νₜ ~ u'ℓ ~ Ujet * φ.
=#
ϵEstimate = Cμ * kEstimate / Ujet / jetDiameter

u5BC(x, t::Real) = ϵEstimate * jetProfile(x,t)
u5BC(t::Real) = x -> u5BC(x,t)
u5IC(x, t::Real) = ϵEstimate * exp(-5.0 * x[2] / Ly) * jetProfile(x,t)
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