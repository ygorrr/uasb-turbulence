#=
Simulação de 1/9 de um tanque.
Cada contêiner tem 1 uasbc
Cada uasb tem 9 furos
=#

using Gridap
using LineSearches: BackTracking
using Plots

#=------------------------
  Propriedades físicas
------------------------=#
β = 0.1 # Coeficiente de dilatação mássica
ρ = 1000.0  # Densidade do fluido
μ = 1.0e-3 # Viscosidade dinâmica do fluido (água) em Pa.s
g = 9.81 # Aceleração da gravidade em m/s²
ĝ = VectorValue(0.0, -1.0)  # Vetor aceleração unitário

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
φ_uasb = 1.8 # m
H_uasb = 2 # m
A_uasb = π*(φ_uasb/2)^2 # m²
Q = 0.500/1000 # m³/s
HDT = A_uasb*H_uasb / Q # s
HDT_horas = HDT/3600
U_uasb = Q/A_uasb # m/s

φ_jet = 0.2 # m
A_jet = π*(φ_jet/2)^2 # m²
xJets = [0.0] # Posições dos jatos ao longo de x
nJets = length(xJets) # Número de jatos
U_jet = Q/(9*A_jet) # m/s

#=------------------------
  Grandezas características do modelo
------------------------=#
Uc = U_jet
Lc = φ_jet
Cc = 1.0

#=------------------------
  Números adimensionais do modelo
------------------------=#
# Número de Reynolds do UASB
Re_uasb = ρ*U_uasb*φ_uasb/μ

# Número de Reynolds do jato
Re_jet = ρ*U_jet*φ_jet/μ 

# Número de Reynolds do modelo
Re = ρ*Uc*Lc/μ
Sc = 1.0
# Ri = β*g*(Cc/Lc)/(Uc/Lc)^2
Ri = 1e-1

n = 50
Lx = 0.6 / Lc
Ly = 2 * Lx
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
  global φ_jet
  local a = 75.0/5
  local xThreshold = log(1999)/(2*a)
  Δx = - φ_jet / (2*Lc) + xThreshold
  return (tanh(a*(x[1]-Δx)) * -tanh(a*(x[1]+Δx)) + 1) / 2
end

function jetProfile(x, t::Real)
  global Uc
  global xJets
  local amplitude = U_jet / Uc
  local velocityY = 0.0
  for xJet in xJets
    velocityY += amplitude * unitJet(x[1]/2 - xJet / Lc)
  end
  return velocityY
end

function scalarProfile(x, t::Real)
  global xJets
  local scalarVal = 0.0
  for xJet in xJets
    scalarVal += unitJet(x[1]/2 - xJet / Lc)
  end
  return scalarVal
end

#=------------------------
  Visualização do perfil de velocidade do jato.
  Melhor comentar esse bloco begin/end se for rodar no cluster.
------------------------=#
begin
  jetContour = ([-1, 1] * φ_jet / 2 .+ xJets[1]) / Lc
  xRange = -Lx/2:0.001:Lx/2
  yComponents = []
  for x in xRange
    append!(yComponents, jetProfile(x, 0))
  end

  plot(xRange, yComponents, xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile")

  scatter!(jetContour, [0.0, 0.0], label="jetContour")
end

# --- Velocidade média na entrada ---
function calcMeanValue(x_min, x_max, f)
  # x_min = -Lx/2
  # x_max = Lx/2
  N = 1000
  x_vals = range(x_min, x_max; length=N)
  dx = (x_max - x_min)/(N-1)
  
  # Método dos trapézios
  integral = f((x_vals[1], 0.0),0)
  for x in x_vals[2:end-1]
    integral += 2 * f((x, 0.0),0)
  end
  integral += (f((x_vals[end], 0.0),0))
  integral = integral * dx/2

  return integral / (x_max - x_min)
end

x_min = -Lx/2
x_max = Lx/2
# Primitiva do WolframAlpha
a = 75.0/5
xThreshold = log(1999)/(2*a)
b = a*(- φ_jet / (2*Lc) + xThreshold)
F(x) = (-1/2/a) * coth(2*b)*(log(cosh(b-a*x)) - log(cosh(a*x+b)))
(F(x_max) - F(x_min))/(x_max - x_min)

# Primitiva do ChatGPT
F(x) = (coth(2*b)/a) * atanh(tanh(b) * tanh(a*x))
(F(x_max) - F(x_min))/(x_max - x_min)

α = 0.15
f(t) = 1-exp(-t/α)
xRange = range(0,5;length=100)
plot(xRange, f)

u1BC(x, t::Real) = VectorValue(0, jetProfile(x,t)) * f(t)
u1BC(t::Real) = x -> u1BC(x,t)
u1IC(x, t::Real) = VectorValue(0, 0)
u1IC(t::Real) = x -> u1IC(x,t)

u2BC(x, t::Real) = 0.0
u2BC(t::Real) = x -> u2BC(x,t)
u2IC(x, t::Real) = 0.0
u2IC(t::Real) = x -> u2IC(x,t)

u3BC(x, t::Real) = scalarProfile(x,t) * f(t)
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

Uref = U_jet / Uc
I = 0.1 
kEstimate = (3/2) * (I * Uref)^2
# kEstimate = 0.015 # Forçar k para 0.015

u4BC(x, t::Real) = kEstimate * scalarProfile(x,t) * f(t)
u4BC(t::Real) = x -> u4BC(x,t)
u4IC(x, t::Real) = kEstimate * scalarProfile(x,t) * f(t)
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
γ = 0.07
ϵEstimate = Cμ^(3/4) * kEstimate^(3/2) / (γ * φ_jet/Lc)

#=------------------------
  Estimativa baseada em argumentos de escala:
    νₜ ~ Ujet * φ
  Essa estimativa vem de argumentos de escala, supondo que as flutuações da velocidade e comprimento de mistura sejam proporcionais à velocidade média do jato e ao diâmetro do jato, ou seja, νₜ ~ u'ℓ ~ Ujet * φ.
------------------------=#
# ϵEstimate = Cμ * kEstimate^2 / (Uref * φ_jet/Lc)

# ϵEstimate = 0.002025 # Forçar epsilon para 0.002025

u5BC(x, t::Real) = ϵEstimate * scalarProfile(x,t) * f(t)
u5BC(t::Real) = x -> u5BC(x,t)
u5IC(x, t::Real) = ϵEstimate * scalarProfile(x,t) * f(t)
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

Us = VectorValue(0.0, -0.4) # Velocidade de decantação de partículas
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
# Experimentando a função de parte simétrica de tensor ε(u1) = ∇(u1) + (∇(u1)') / 2
resNS(t, u1, u2, u3, u4, u5, v1) = 
  ∫( v1 ⋅ ∂t(u1) )dΩ +
  ∫( v1 ⋅ (∇(u1)' ⋅ u1) )dΩ +
  ∫( (1/Re + (nuT∘(u4,u5))) * (∇(v1)⊙ε(u1))*2 )dΩ - 
  ∫( (∇ ⋅ v1) * u2 )dΩ -
  ∫( Ri * (v1 ⋅ ĝ) * u3 )dΩ

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
  ∫( v4 * Ri * (nuT∘(u4,u5)) / σ * (ĝ ⋅ ∇(u3)) * (tanh∘(10.0*u4/kEstimate)) )dΩ

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

#=------------------------
  Preparo do diretório de resultados
------------------------=#

caseComment = """
----------------------------
      Case description
----------------------------
Spatial dimensions: 2
Transient:          true
Reynolds number:    $ReStr
Richardson number:  $Ri
Number of jets:     $nJets
Settling velocity:  $Us
"""

dirPath = joinpath(@__DIR__, "output")
if !isdir(dirPath)
  mkdir(dirPath)
end

ReStr = Int(round(Re, RoundDown))
dirPath = joinpath(dirPath,"Re$ReStr-nJets$nJets-Us$Us")
if isdir(dirPath)
  rm(dirPath, recursive=true)
end
mkdir(dirPath)

filePath = joinpath(dirPath, "case-description.txt")
touch(filePath)
open(filePath, "w") do file
    write(file, caseComment)
end

#=------------------------
  Solução e escrita dos resultados
------------------------=#

writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets-Us$Us/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph, "ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])

it = 1
totalIts = T/Δt
for (t,uₕ) in uₕₜ
  global it
  local uh,ph,ch,kh,epsilonh
  uh, ph, ch, kh, epsilonh = uₕ
  
  println("Iteration $it/$totalIts")
  
  if(mod(it,1)==0)
    writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets-Us$Us/uasbcp$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])
  end

  it = it + 1
end