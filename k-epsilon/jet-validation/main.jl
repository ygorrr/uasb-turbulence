using Gridap
using LineSearches: BackTracking
using Plots

#=------------------------
  Propriedades físicas
------------------------=#
ρ = 1.0  # Densidade do fluido
nu = 1.0e-3 # Viscosidade cinemática do fluido (água) em m²/s
g = 9.81 # Aceleração da gravidade em m/s²
gHat = VectorValue(0.0, -1.0)  # Vetor aceleração unitário

#=------------------------
  Constantes do modelo k-epsilon
------------------------=#
Cμ = 0.09  # Constante de turbulência
Cϵ1 = 1.44  # Constante de produção de epsilon
Cϵ2 = 1.92  # Constante de destruição de epsilon
σk = 1.0  # Constante de difusão de k
σϵ = 1.3  # Constante de difusão de epsilon
σ = 1.0  # Constante de difusão de concentração

slitWidth = 0.1 # m

Re = 200

Uj = Re * nu / slitWidth  # Velocidade do jato

n = 100
Lx = 20 * slitWidth
Ly = 200 * slitWidth
domain = (-Lx/2, Lx/2, 0, Ly)
partition = (2*n, n)
model = CartesianDiscreteModel(domain, partition)
writevtk(model, @__DIR__)

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "top", [6,])
add_tag_from_tags!(labels, "bottom", [5,])
add_tag_from_tags!(labels, "left", [1, 3, 7])
add_tag_from_tags!(labels, "right", [2, 4, 8])

neumannTags = ["left", "right", "top"]

#=------------------------
Incógnitas:

Espaços de funções de teste: V1, V2, V3, V4, V5
Espaços de de funções de ensaio: S1, S2, S3, S4, S5
------------------------=#

order = 2

reffe_U = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
VU = TestFESpace(model, reffe_U, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

reffe_P = ReferenceFE(lagrangian, Float64, order-1; space=:P)
VP = TestFESpace(model, reffe_P, conformity=:L2, dirichlet_tags=["top", "left", "right"])

reffe_k = ReferenceFE(lagrangian, Float64, order)
Vk = TestFESpace(model, reffe_k, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

reffe_ε = ReferenceFE(lagrangian, Float64, order)
Vε = TestFESpace(model, reffe_ε, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

#=------------------------
  Condições iniciais e de contorno
------------------------=#

function unitJet(x)
  global slitWidth
  local a = 75.0*3
  local xThreshold = log(1999)/(2*a)
  Δx = -slitWidth / 2 + xThreshold
  return (tanh(a*(x[1]-Δx)) * -tanh(a*(x[1]+Δx)) + 1) / 2
end

function jetProfile(x)
  global Uj
  # amplitude = 1.0 # Forçar amplitude para 1.0
  local velocityY = 0.0
  velocityY += Uj * unitJet(x[1])
  return VectorValue(0.0, velocityY)
end

function scalarProfile(x)
  local scalarVal = 0.0
  scalarVal += unitJet(x[1])
  return scalarVal
end

#=------------------------
  Visualização do perfil de velocidade do jato.
  Melhor comentar esse bloco begin/end se for rodar no cluster!
------------------------=#
begin
  jetContour = [- slitWidth / 2, slitWidth / 2]
  xRange = -Lx/2:0.001:Lx/2
  yComponents = []
  for x in xRange
    append!(yComponents, jetProfile(x)[2])
  end

  plot(xRange, yComponents, xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile")

  scatter!(jetContour, [0.0, 0.0], label="jetContour")
end

UBC(x) = jetProfile(x)

PBC(x) = 0.0

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

I = 0.1 
kEstimate = (3/2) * (I * Uj)^2
# kEstimate = 0.015 # Forçar k para 0.015

kBC(x) = kEstimate * scalarProfile(x)

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
ϵEstimate = Cμ^(3/4) * kEstimate^(3/2) / (α * slitWidth)

#=------------------------
  Estimativa baseada em argumentos de escala:
    νₜ ~ Ujet * φ
  Essa estimativa vem de argumentos de escala, supondo que as flutuações da velocidade e comprimento de mistura sejam proporcionais à velocidade média do jato e ao diâmetro do jato, ou seja, νₜ ~ u'ℓ ~ Ujet * φ.
------------------------=#
# ϵEstimate = Cμ * kEstimate^2 / (Uref * jetDiameter/Lc)

# ϵEstimate = 0.002025 # Forçar epsilon para 0.002025

εBC(x) = ϵEstimate * scalarProfile(x)

SU = TrialFESpace(VU, [UBC])
SP = TrialFESpace(VP, [PBC])
Sk = TrialFESpace(Vk, [kBC])
Sε = TrialFESpace(Vε, [εBC])

Y = MultiFieldFESpace([VU, VP, Vk, Vε])
X = MultiFieldFESpace([SU, SP, Sk, Sε])

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

Γ = BoundaryTriangulation(Ω, neumannTags)
dΓ = Measure(Γ, degree)

minVal = 1e-2 # Constante de proteção contra divisão por zero

nuT(k, ε) = Cμ * k * k / max(ε, minVal)
kProduction(∇U, k, ε) = nuT(k, ε) * ((∇U + ∇U') ⊙ ∇U)
epsilonProduction(∇U, k, ε) = Cϵ1 * kProduction(∇U, k, ε) * max(ε, minVal) / max(k, minVal)
epsilonDestruction(k, ε) = Cϵ2 * ε * ε / max(k, minVal)

#=------------------------
  Resíduos
------------------------=#
# Equação de Navier-Stokes
# Experimentando a função de parte simétrica de tensor ε(u1) = ∇(u1) + (∇(u1)') / 2

resNS(U, P, k, ε, vU) = 
  ∫( vU ⋅ (∇(U)' ⋅ U) )dΩ +
  ∫( (1/Re + (nuT∘(k,ε))) * (∇(vU)⊙(∇(U) + (∇(U))')) )dΩ - 
  ∫( (∇ ⋅ vU) * P )dΩ 

# Equação da continuidade
resCont(U, vP) =
  ∫( vP * (∇ ⋅ U) )dΩ

# Equação de k
resk(U, k, ε, vk) = 
  ∫( vk ⋅ (∇(k)' ⋅ U) )dΩ +
  ∫( (nuT∘(k,ε)) / σk * ∇(vk)⊙∇(k) )dΩ -
  ∫( vk * (kProduction∘(∇(U),k,ε)) )dΩ +
  ∫( vk * ε * (tanh∘(10.0*k/kEstimate)) )dΩ

# Equação de epsilon
resEpsilon(U, k, ε, vε) =
  ∫( vε ⋅ (∇(ε)' ⋅ U) )dΩ +
  ∫( (nuT∘(k,ε)) / σϵ * ∇(vε)⊙∇(ε) )dΩ -
  ∫( vε * (epsilonProduction∘(∇(U),k,ε)) )dΩ +
  ∫( vε * (epsilonDestruction∘(k,ε)) * (tanh∘(10.0*ε/ϵEstimate)) )dΩ

res((U, P, k, ε), (vU, vP, vk, vε)) =
  resNS(U, P, k, ε, vU) +
  resCont(U, vP) +
  resk(U, k, ε, vk) +
  resEpsilon(U, k, ε, vε)

op = FEOperator(res,X,Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)
solver = FESolver(nls)
Uh, Ph, kh, εh = solve(solver,op)

#=------------------------
  Preparo do diretório de resultados
------------------------=#

dirPath = joinpath(@__DIR__, "output")
if !isdir(dirPath)
  mkdir(dirPath)
end

ReStr = Int(round(Re, RoundDown))
dirPath = joinpath(dirPath,"Re$ReStr")
if !isdir(dirPath)
  mkdir(dirPath)
end

caseComment = """
----------------------------
      Case description
----------------------------
Spatial dimensions: 2
Transient:          False
Reynolds number:    $ReStr
"""

filePath = joinpath(dirPath, "case-description.txt")
touch(filePath)
open(filePath, "w") do file
    write(file, caseComment)
end

#=------------------------
  Solução e escrita dos resultados
------------------------=#

writevtk(Ω,(@__DIR__)*"/output/Re$ReStr/plane-jet.vtu",cellfields=["Uh"=>Uh,"Ph"=>Ph, "kh"=>kh, "epsilonh"=>εh])