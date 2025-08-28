using Gridap
using LineSearches: BackTracking
using Plots

# --- Parâmetros físicos ---
ρ = 1.0                     # Densidade do fluido [kg/m³]
nu = 1.0e-3                 # Viscosidade cinemática do fluido [m²/s]
g = 9.81                    # Aceleração da gravidade [m/s²]
ĝ = VectorValue(0.0, -1.0)  # Vetor aceleração gravitacional unitário

# --- Parâmetros geométricos ---
Lj = 0.1      # Largura de entrada do jato [m]
xJets = [0.0] # Coordenadas dos centros dos jatos em x

# --- Constantes do modelo k-ϵ ---
Cμ = 0.09     # Constante de turbulência
Cϵ1 = 1.44    # Constante de produção de ϵ
Cϵ2 = 1.92    # Constante de destruição de ϵ
σk = 1.0      # Número de Prandtl turbulento de k
σϵ = 1.3      # Número de Prandtl turbulento de ϵ
σ = 1.0       # Número de Prandtl turbulento de C

Re = 200

Uj = Re * nu / Lj  # Velocidade do jato

# --- Função de perfil do jato ---
# Perfil da componente vertical um jato unitário
function unitJet(x)
  global Lj
  local a = 75.0*2
  local λ = 0.999
  local x0 = atanh(λ) / a
    
  return (tanh(a*(x[1]-x0+Lj/2)) * - tanh(a*(x[1]+x0-Lj/2)) +1)/2
end

# Perfil de velocidade dos jatos
function jetProfile(x)
  global xJets
  local Uy = 0.0
  
  for xJet in xJets
    Uy += unitJet(x[1] - xJet)
  end

  Uy = VectorValue(0.0, Uy)
  
  return Uy
end

# Perfil de velocidade dos jatos
function scalarProfile(x)
  global xJets
  local θ = 0.0
  
  for xJet in xJets
    θ += unitJet(x[1] - xJet)
  end
  
  return θ
end

# Função para visualização do perfil do jato
function plotJetProfile()
  local jetContour = [- Lj/2, Lj/2] .+ xJets[1]
  local xRange = -Lx/2:0.001:Lx/2
  local Uy = []

  for x in xRange
    append!(Uy, jetProfile(x)[2])
  end

  plt = plot(xRange, Uy, xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile")

  scatter!(plt, jetContour, [0.0, 0.0], label="jetContour")
end

# Estimativa de k no jato
I = 0.1               # Intensidade de turbulência
kj = (3/2)*(I*Uj)^2   # Estimativa de k na saída do jato [m²/s²]

# Estimativa de ϵ no jato
α = 0.07
ϵj = Cμ^(3/4) * kj^(3/2) / (α*Lj)   # Estimativa de ϵ na saída do jato [m²/s³]

# --- Malha e modelo geométrico ---
n = 100
Lx = 10 * Lj
Ly = 10 * Lj
domain = (-Lx/2, Lx/2, 0, Ly)
partition = (n, n)
model = CartesianDiscreteModel(domain, partition)
writevtk(model, @__DIR__)

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "top", [6,])
add_tag_from_tags!(labels, "bottom", [5,])
add_tag_from_tags!(labels, "left", [1, 3, 7])
add_tag_from_tags!(labels, "right", [2, 4, 8])

neumannTags = ["left", "right", "top"]

# --- Método dos Elementos Finitos ---
order = 2

reffe_U = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
reffe_P = ReferenceFE(lagrangian, Float64, order-1; space=:P)
reffe_k = ReferenceFE(lagrangian, Float64, order)
reffe_ϵ = ReferenceFE(lagrangian, Float64, order)

VU = TestFESpace(model, reffe_U, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])
VP = TestFESpace(model, reffe_P, conformity=:L2, dirichlet_tags=["top", "left", "right"])
Vk = TestFESpace(model, reffe_k, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])
Vϵ = TestFESpace(model, reffe_ϵ, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

UBC(x) = Uj * jetProfile(x)
PBC(x) = 0.0
kBC(x) = kj * scalarProfile(x)
ϵBC(x) = ϵj * scalarProfile(x)

SU = TrialFESpace(VU, [UBC])
SP = TrialFESpace(VP, [PBC, PBC, PBC])
Sk = TrialFESpace(Vk, [kBC])
Sϵ = TrialFESpace(Vϵ, [ϵBC])

Y = MultiFieldFESpace([VU, VP, Vk, Vϵ])
X = MultiFieldFESpace([SU, SP, Sk, Sϵ])

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

Γ = BoundaryTriangulation(Ω, neumannTags)
dΓ = Measure(Γ, degree)

minVal = 1e-2 # Constante de regularização

nuT(k, ϵ) = Cμ * k * k / max(ϵ, minVal)
kProduction(∇U, k, ϵ) = nuT(k, ϵ) * ((∇U + ∇U') ⊙ ∇U)
epsilonProduction(∇U, k, ϵ) = Cϵ1 * kProduction(∇U, k, ϵ) * max(ϵ, minVal) / max(k, minVal)
epsilonDestruction(k, ϵ) = Cϵ2 * ϵ * ϵ / max(k, minVal)

#=------------------------
  Resíduos
------------------------=#
# Equação de Navier-Stokes
# Experimentando a função de parte simétrica de tensor ε(u1) = ∇(u1) + (∇(u1)') / 2

resNS(U, P, k, ϵ, vU) = 
  ∫( vU ⋅ (∇(U)' ⋅ U) )dΩ +
  ∫( (1/Re + (nuT∘(k,ϵ))) * (∇(vU)⊙(∇(U) + (∇(U))')) )dΩ - 
  ∫( (∇ ⋅ vU) * P )dΩ 

# Equação da continuidade
resCont(U, vP) =
  ∫( vP * (∇ ⋅ U) )dΩ

# Equação de k
resk(U, k, ϵ, vk) = 
  ∫( vk ⋅ (∇(k)' ⋅ U) )dΩ +
  ∫( (nuT∘(k,ϵ)) / σk * ∇(vk)⊙∇(k) )dΩ -
  ∫( vk * (kProduction∘(∇(U),k,ϵ)) )dΩ +
  ∫( vk * ϵ * (tanh∘(10.0*k/kj)) )dΩ

# Equação de epsilon
resEpsilon(U, k, ϵ, vϵ) =
  ∫( vϵ ⋅ (∇(ϵ)' ⋅ U) )dΩ +
  ∫( (nuT∘(k,ϵ)) / σϵ * ∇(vϵ)⊙∇(ϵ) )dΩ -
  ∫( vϵ * (epsilonProduction∘(∇(U),k,ϵ)) )dΩ +
  ∫( vϵ * (epsilonDestruction∘(k,ϵ)) * (tanh∘(10.0*ϵ/ϵj)) )dΩ

res((U, P, k, ϵ), (vU, vP, vk, vϵ)) =
  resNS(U, P, k, ϵ, vU) +
  resCont(U, vP) +
  resk(U, k, ϵ, vk) +
  resEpsilon(U, k, ϵ, vϵ)

op = FEOperator(res,X,Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)
solver = FESolver(nls)

# Chutes iniciais para P, U, k e ϵ
P0 = interpolate_everywhere(PBC, VP)
P0_array = get_free_dof_values(P0)

U_gess(x) = UBC(x) * exp(-5*x[2]/Ly)
U0 = interpolate_everywhere(U_gess, VU)
U0_array = get_free_dof_values(U0)

k_gess(x) = kBC(x) * exp(-5*x[2]/Ly)
k0 = interpolate_everywhere(k_gess, Vk)
k0_array = get_free_dof_values(k0)

ϵ_gess(x) = ϵBC(x) * exp(-5*x[2]/Ly)
ϵ0 = interpolate_everywhere(ϵ_gess, Vϵ)
ϵ0_array = get_free_dof_values(ϵ0)

x = vcat(U0_array, P0_array, k0_array, ϵ0_array)
uh0 = FEFunction(X, x)
usol = solve!(uh0, solver, op)

Uh, Ph, kh, ϵh = usol[1]

#=------------------------
  Preparo do diretório de resultados
------------------------=#
begin
  dirPath = joinpath(@__DIR__, "output")
  if !isdir(dirPath)
    mkdir(dirPath)
  end
  
  dirPath = joinpath(dirPath,"Re$(round(Re))")
  if !isdir(dirPath)
    mkdir(dirPath)
  end
  
  caseComment = """
  ----------------------------
  Case description
  ----------------------------
  Spatial dimensions: 2
  Transient:          False
  Reynolds number:    $(round(Re))
  """
  
  filePath = joinpath(dirPath, "case-description.txt")
  touch(filePath)
  open(filePath, "w") do file
    write(file, caseComment)
  end
end

#=------------------------
  Solução e escrita dos resultados
------------------------=#

writevtk(Ω,(@__DIR__)*"/output/Re$ReStr/plane-jet.vtu",cellfields=["Uh"=>Uh,"Ph"=>Ph, "kh"=>kh, "ϵh"=>ϵh])