using Gridap
using LineSearches: BackTracking
using Plots

#=------------------------
  Propriedades físicas
------------------------=#
β = 0.1 # Coeficiente de dilatação mássica
ρ = 1.0  # Densidade do fluido
nu = 1.0e-3 # Viscosidade cinemática do fluido (água) em m²/s
g = VectorValue(0.0, -1.0)  # Vetor aceleração unitário

#=------------------------
  Constantes do modelo k-epsilon
------------------------=#
Cμ = 0.09  # Constante de turbulência
Cϵ1 = 1.44  # Constante de produção de epsilon
Cϵ2 = 1.92  # Constante de destruição de epsilon
σk = 1.0  # Constante de difusão de k
σϵ = 1.3  # Constante de difusão de epsilon
σ = 1.0  # Constante de difusão de concentração

#=------------------------
  Parâmetros geométricos e operacionais do UASB
------------------------=#
uasbDiameter = 1.8 # m
uasbHeight = 3 * uasbDiameter # m
uasbArea = π * (uasbDiameter/2)^2 # m²
HDT = 1*60*60 # s
volFlowRate = uasbArea * uasbHeight / HDT # m³/s
meanVelocity = uasbHeight / HDT # m/s

jetDiameter = 0.2 # m
jetArea = π * (jetDiameter/2)^2 # m²
xJets = [-2, 2]
nJets = length(xJets) # Número de jatos
meanJetVelocity = volFlowRate / (nJets * jetArea) # m/s

#=------------------------
  Grandezas características do modelo
------------------------=#
Uc = meanJetVelocity
Lc = jetDiameter
Cc = 1.0

#=------------------------
  Números adimensionais do modelo
------------------------=#
# Número de Reynolds do UASB
Re_uasb = meanVelocity * uasbDiameter / nu 

# Número de Reynolds do jato
Re_jet = meanJetVelocity * jetDiameter / nu 

Re = Uc * Lc / nu
# Re = 300
Sc = 1.0
# Ri = β*g*(Cref/Lref)/(Uref/Lref)^2
Ri = 1e-1

n = 50
Lx = uasbDiameter / Lc
Ly = uasbHeight / Lc
domain = (-Lx/2, Lx/2, 0, Ly)
partition = (n, 2*n)
model = CartesianDiscreteModel(domain, partition; isperiodic=(true,false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "top", [6,])
add_tag_from_tags!(labels, "bottom", [5,])

#=------------------------
Incógnitas:
u1 -> vetor velocidade média
u2 -> pressão média
U3 -> concentração média de partículas
u4 -> k
u5 -> epsilon

Espaços de funções de teste: V1, V2, V3, V4, V5
Espaços de de funções de ensaio: S1, S2, S3, S4, S5
------------------------=#

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

#=------------------------
  Condições iniciais e de contorno
------------------------=#

function unitJet(x)
  global jetDiameter
  local a = 75.0/5
  local xThreshold = log(1999)/(2*a)
  Δx = -jetDiameter / (2*Lc) + xThreshold
  return (tanh(a*(x[1]-Δx)) * -tanh(a*(x[1]+Δx)) + 1) / 2
end

function jetProfile(x, t::Real)
  global Uc
  global xJets
  local amplitude = meanJetVelocity / Uc
  # amplitude = 1.0 # Forçar amplitude para 1.0
  local velocityY = 0.0
  for xJet in xJets
    velocityY += amplitude * unitJet(x[1] - xJet)
  end
  return VectorValue(0.0, velocityY)
end

function scalarProfile(x, t::Real)
  global xJets
  local scalarVal = 0.0
  for xJet in xJets
    scalarVal += unitJet(x[1] - xJet)
  end
  return scalarVal
end

#=------------------------
  Visualização do perfil de velocidade do jato.
  Melhor comentar esse bloco begin/end se for rodar no cluster!
------------------------=#
begin
  jetContour = xJets[1] .+ [- jetDiameter / (2*Lc), jetDiameter / (2*Lc)]
  xRange = -Lx/2:0.001:Lx/2
  yComponents = []
  for x in xRange
    append!(yComponents, jetProfile(x, 0)[2])
  end

  plot(xRange, yComponents, xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile")

  scatter!(jetContour, [0.0, 0.0], label="jetContour")
end

u1BC(x, t::Real) = jetProfile(x,t)
u1BC(t::Real) = x -> u1BC(x,t)
u1IC(x, t::Real) = exp(-5.0 * x[2] / Ly) * jetProfile(x,t)
u1IC(t::Real) = x -> u1IC(x,t)

u2BC(x, t::Real) = 0.0
u2BC(t::Real) = x -> u2BC(x,t)
u2IC(x, t::Real) = 0.0
u2IC(t::Real) = x -> u2IC(x,t)

u3BC(x, t::Real) = scalarProfile(x,t)
u3BC(t::Real) = x -> u3BC(x,t)
u3IC(x, t::Real) = 0.0
u3IC(t::Real) = x -> u3IC(x,t)

#=------------------------
  Estimativa de k na saída do jato
------------------------=#

#=------------------------
Intensidade de turbulência. Valores típicos:
Jato livre: 0.05 ~ 0.1
Escoamento altamente turbulento: até 0.2
I é definida como a razão entre a velocidade de flutuação e a velocidade média do jato.
I = u' / Uref
------------------------=#

Uref = meanJetVelocity / Uc
I = 0.1 
kEstimate = (3/2) * (I * Uref)^2
kEstimate = 0.015 # Forçar k para 0.015

u4BC(x, t::Real) = kEstimate * scalarProfile(x,t)
u4BC(t::Real) = x -> u4BC(x,t)
u4IC(x, t::Real) = kEstimate * exp(-5.0 * x[2] / Ly) * scalarProfile(x,t)
u4IC(t::Real) = x -> u4IC(x,t)

#=------------------------
  Adote apenas uma das estratégias abaixo para a estimativa de epsilon!
------------------------=#

#=------------------------
  Estimativa dimensional de epsilon na saída do jato.
  Supõe-se que u' ~ √k e que νₜ ~ u' * ℓ. Usando a relação entre k e ϵ:
    ϵ = Cμ^(3/4) * k^(3/2) / ℓ
  onde Cμ é uma constante de turbulência e ℓ é uma escala de comprimento. O expoente de Cμ é empírico e não vem diretamente da teoria.
  ℓ pode ser estimada por uma fração do diâmetro do jato:
    ℓ = α * φ
  onde α é uma constante, tipicamente entre 0.07 e 0.1.
------------------------=#
α = 0.07
ϵEstimate = Cμ^(3/4) * kEstimate^(3/2) / (α * jetDiameter)

#=------------------------
  Estimativa baseada em argumentos de escala:
    νₜ ~ Ujet * φ
  Essa estimativa vem de argumentos de escala, supondo que as flutuações da velocidade e comprimento de mistura sejam proporcionais à velocidade média do jato e ao diâmetro do jato, ou seja, νₜ ~ u'ℓ ~ Ujet * φ.
------------------------=#
ϵEstimate = Cμ * kEstimate^2 / (Uref * jetDiameter/Lc)

ϵEstimate = 0.002025 # Forçar epsilon para 0.002025

u5BC(x, t::Real) = ϵEstimate * scalarProfile(x,t)
u5BC(t::Real) = x -> u5BC(x,t)
u5IC(x, t::Real) = ϵEstimate * exp(-5.0 * x[2] / Ly) * scalarProfile(x,t)
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

Us = VectorValue(0.0, -0.1) # Velocidade de decantação de partículas
minVal = 1e-2 # Constante de proteção contra divisão por zero

w(u) = u + Us
nuT(u4, u5) = Cμ * u4 * u4 / max(u5, minVal)
kProduction(∇u1, u4, u5) = nuT(u4, u5) * ((∇u1 + ∇u1') ⊙ ∇u1)
epsilonProduction(∇u1, u4, u5) = Cϵ1 * kProduction(∇u1, u4, u5) * max(u5, minVal) / max(u4, minVal)
epsilonDestruction(u4, u5) = Cϵ2 * u5 * u5 / max(u4, minVal)

#=------------------------
  Resíduos
------------------------=#

# Equação de Navier-Stokes 
resNS(t, u1, u2, u3, u4, u5, v1) = 
  ∫( v1 ⋅ ∂t(u1) )dΩ +
  ∫( v1 ⋅ (∇(u1)' ⋅ u1) )dΩ +
  ∫( (1/Re + (nuT∘(u4,u5))) * ∇(v1)⊙∇(u1) )dΩ - 
  ∫( (∇ ⋅ v1) * u2 )dΩ -
  ∫( Ri * (v1 ⋅ g) * u3 )dΩ

# Equação da continuidade
resCont(t, u1, v2) =
  ∫( v2 * (∇ ⋅ u1) )dΩ

# Equação de transporte de partículas
resC(t, u1, u3, u4, u5, v3) =
  ∫( v3 * ∂t(u3) )dΩ +
  ∫( v3 ⋅ inner(w(u1), ∇(u3)) )dΩ +
  ∫( (1/Re * 1/Sc + (nuT∘(u4,u5))/σ) * ∇(u3)⊙∇(v3) )dΩ

# Equação de k
resk(t, u1, u3, u4, u5, v4) = 
  ∫( v4 * ∂t(u4) )dΩ +
  ∫( v4 ⋅ (∇(u4)' ⋅ u1) )dΩ +
  ∫( (nuT∘(u4,u5)) / σk * ∇(v4)⊙∇(u4) )dΩ -
  ∫( v4 * (kProduction∘(∇(u1),u4,u5)) )dΩ +
  ∫( v4 * u5 * (tanh∘(10.0*u4/kEstimate)) )dΩ +
  ∫( v4 * Ri * (nuT∘(u4,u5)) / σ * (g ⋅ ∇(u3)) * (tanh∘(10.0*u4/kEstimate)) )dΩ

# Equação de epsilon
resEpsilon(t, u1, u4, u5, v5) =
  ∫( v5 * ∂t(u5) )dΩ +
  ∫( v5 ⋅ (∇(u5)' ⋅ u1) )dΩ +
  ∫( (nuT∘(u4,u5)) / σϵ * ∇(u5)⊙∇(v5) )dΩ -
  ∫( v5 * (epsilonProduction∘(∇(u1),u4,u5)) )dΩ +
  ∫( v5 * (epsilonDestruction∘(u4,u5)) * (tanh∘(10.0*u5/ϵEstimate)) )dΩ

res(t, (u1, u2, u3, u4, u5), (v1, v2, v3, v4, v5)) =
  resNS(t, u1, u2, u3, u4, u5, v1) +
  resCont(t, u1, v2) +
  resC(t, u1, u3, u4, u5, v3) +
  resk(t, u1, u3, u4, u5, v4) +
  resEpsilon(t, u1, u4, u5, v5)

op = TransientFEOperator(res,X,Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)

CFL = 1.0/2.0
Δt = CFL*Lx/n
θ = 1

ode_solver = ThetaMethod(nls,Δt,θ)

U₀ = interpolate_everywhere([u1IC(0),u2IC(0),u3IC(0),u4IC(0),u5IC(0)],X(0.0))
t₀ = 0.0
T = 100.0
uₕₜ = solve(ode_solver,op,t₀,T,U₀)
it = 0
uh, ph, ch, kh, epsilonh = U₀

if !isdir((@__DIR__)*"/results")
  mkdir((@__DIR__)*"/results")
end

writevtk(Ω,(@__DIR__)*"/results/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph, "ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])

it = 1
totalIts = T/Δt
for (t,uₕ) in uₕₜ
  global it
  local uh,ph,ch,kh,epsilonh
  uh, ph, ch, kh, epsilonh = uₕ
  
  println("Iteration $it/$totalIts")
  
  if(mod(it,1)==0)
    writevtk(Ω,(@__DIR__)*"/results/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])
  end

  it = it + 1
end