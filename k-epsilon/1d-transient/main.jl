using Gridap
using Plots
using LineSearches: BackTracking

# -------------------------------------------------------------------
#                      Parâmetros e funções
# -------------------------------------------------------------------

# Parâmetros geométricos
const reactorDiam = 1.8 # m
const reactorHeight = 3 * reactorDiam # m
const reactorCrossArea = π*reactorDiam^2/4
const hdt = 1*60*60 # s
const reactorVelocity = reactorDiam/hdt # m/s

const jetDiam = 0.2 # m
const nJets = 1
const totalJetArea = (π*jetDiam^2/4)*nJets
const totalJetPerimeter = (π*jetDiam)*nJets
const jetVelocity = reactorVelocity * (reactorCrossArea/totalJetArea) # m/s

# Parâmetros físicos
const gravityAccel = 9.81 # m/s²
const fluidSpecificMass = 997.0 # kg/m³
const fluidDynamicViscosity = 0.00089 # Pa*s
const fluidKinematicViscosity = fluidDynamicViscosity/fluidSpecificMass
const molecularDiffC = 1.0 # m²/s
const beta = 0.1

# Constantes do modelo k-epsilon
const Cμ = 0.09
const Cϵ1 = 1.44
const Cϵ2 = 1.92
const σ = 1.0
const σk = 1.0
const σϵ = 1.3

# Grandezas características
const URef = jetVelocity
const LRef = jetDiam
const CRef = 1

# Grupos adimensionais
# const Re = cVelocity*cLength/fluidKinematicViscosity
const Re = 300
# const Sc = fluidKinematicViscosity/molecularDiffC
const Sc = 1.0
# const Fr = cVelocity/sqrt(9.81*cLength)
# const Ri = (1/Fr^2)*cConcentration*beta
const Ri = 0.1

U = VectorValue(reactorVelocity / URef)
Us = VectorValue(-0.1)
W = U + Us
ĝ = VectorValue(-1.0)

eddyVisc(k, ϵ) = Cμ*k*k/max(ϵ, 0.01)

res(t, (C, k, ϵ), (v1, v2, v3)) =
    ∫(v1*∂t(C) + v1*W⋅∇(C) + (1/Re/Sc + (eddyVisc∘(k,ϵ))/σ)*∇(v1)⋅∇(C))dΩ +
    ∫(v2*∂t(k) + v2*U⋅∇(k) + (eddyVisc∘(k,ϵ))/σk*∇(v2)⋅∇(k) + v2*ϵ*(tanh∘(10.0*k/dirichlet_k)) + v2*(Ri*(eddyVisc∘(k,ϵ))/σ)*(tanh∘(10.0*k/dirichlet_k)))dΩ +
    ∫(v3*∂t(ϵ) + v3*U⋅∇(ϵ) + (eddyVisc∘(k, ϵ))/σϵ*∇(v3)⋅∇(ϵ) + v3*Cϵ2*ϵ*ϵ/k*(tanh∘(10.0*ϵ/dirichlet_epsilon)))dΩ

# -------------------------------------------------------------------
#                      Formulação de MEF
# -------------------------------------------------------------------

n = 150
model_Ω = CartesianDiscreteModel((0, reactorHeight/LRef), (n))
Ω = Interior(model_Ω)
dΩ = Measure(Ω, 2)

# Tipo de elemento, tipo de função de forma, grau da função de forma
refFE = ReferenceFE(lagrangian, Float64, 1)

# Espaços de funções de teste
testSpace_C = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="boundary")
testSpace_k = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="boundary")
testSpace_epsilon = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="boundary")

# Condições de contorno
bc_C(x, t::Real) = (1 - LRef/reactorHeight*x[1])
bc_C(t::Real) = x -> bc_C(x, t)

I = 0.1
dirichlet_k = (3/2)*(I*jetVelocity/URef)^2 #(1/2)*(reactorCrossArea/totalJetArea)^2
bc_k(x, t::Real) = dirichlet_k*(1 - LRef/reactorHeight*x[1])
bc_k(t::Real) = x -> bc_k(x, t)

dirichlet_epsilon = Cμ*dirichlet_k/(0.01*jetVelocity*jetDiam/URef/LRef) #2*(reactorDiam/jetDiam)*(reactorCrossArea/totalJetArea)^3
bc_epsilon(x, t::Real) = dirichlet_epsilon*(1 - LRef/reactorHeight*x[1])
bc_epsilon(t::Real) = x -> bc_epsilon(x, t)

# Espaços de funções candidatas
trialSpace_C = TransientTrialFESpace(testSpace_C, bc_C)
trialSpace_k = TransientTrialFESpace(testSpace_k, bc_k)
trialSpace_epsilon = TransientTrialFESpace(testSpace_epsilon, bc_epsilon)

Y = MultiFieldFESpace([testSpace_C, testSpace_k, testSpace_epsilon])
X = TransientMultiFieldFESpace([trialSpace_C, trialSpace_k, trialSpace_epsilon])

op = TransientFEOperator(res, X, Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)

CFL = 0.25
dx = reactorHeight/LRef/n
Δt = CFL*dx
θ = 1
ode_solver = ThetaMethod(nls,Δt,θ)

u_θ1 = interpolate_everywhere(bc_C(0), trialSpace_C(0.0))
u_θ2 = interpolate_everywhere(bc_k(0), trialSpace_k(0.0))
u_θ3 = interpolate_everywhere(bc_epsilon(0), trialSpace_epsilon(0.0))
u_θ = interpolate_everywhere([u_θ1, u_θ2, u_θ3], X(0.0))
t_θ = 0.0
T = 50.0

u_ht = solve(ode_solver, op, t_θ, T, u_θ)

#=------------------------
  Preparo do diretório de resultados
------------------------=#

caseDesc = """
----------------------------
Case description
----------------------------
Spatial dimensions: 1
Transient:          true
Reynolds number:    $ReStr
Richardson number:  $Ri
Number of jets:     $nJets
Settling velocity:  $Us
"""

path = joinpath(@__DIR__, "case-description.txt")
if isfile(path)
  rm(path)
end
touch(path)
print(path, caseDesc)

path = joinpath(@__DIR__,"case.log")
if isfile(path)
  rm(path)
end
touch(path)

path = joinpath(@__DIR__, "output")
if !isdir(dirPath)
  mkdir(dirPath)
end

ReStr = Int(round(Re, RoundDown))
path = joinpath(path,"Re$ReStr-nJets$nJets")
if isdir(dirPath)
  rm(dirPath, recursive=true)
end
mkdir(dirPath)

#=------------------------
  Solução e escrita dos resultados
------------------------=#

ch, kh, epsilonh = u_θ
it = 0

writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets/uasbcp$it.vtu",cellfields=["ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh, "nuT"=>eddyVisc∘(kh, epsilonh)])

totalIts = Int(round(T/Δt, RoundUp))
for (t,u_h) in u_ht
  global it
  local ch, kh, epsilonh
  it = it + 1
  println("Iterations completed: $it/$totalIts")
  ch, kh, epsilonh = u_h 
  if (mod(it,1)==0)
    writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets/uasbcp$it.vtu",cellfields=["ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh, "nuT"=>eddyVisc∘(kh, epsilonh)])
  end
end

println("Simulation completed without errors.")