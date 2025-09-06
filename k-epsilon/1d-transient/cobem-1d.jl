using Gridap
using Plots
using LineSearches: BackTracking

# -------------------------------------------------------------------
#                      Parâmetros e funções
# -------------------------------------------------------------------

# Parâmetros geométricos
φ_uasb = 1.8 # m
H_uasb = 2 # m
A_uasb = π*(φ_uasb/2)^2
Q = 0.5 / 1000 # m³/s
HDT = H_uasb * A_uasb / Q # s
HDT_horas = HDT / 3600
U_uasb = H_uasb/HDT # m/s

φ_jet = 0.2 # m
A_jet = π*(φ_jet/2)^2 # m²
U_jet = Q/(9*A_jet) # m/s (9 orifícios por UASB)
nJets = 1

# Parâmetros físicos
g = 9.81 # m/s²
ρ = 1000.0 # kg/m³
μ = 1.0e-3 # Pa*s
β = 0.1

# Constantes do modelo k-epsilon
Cμ = 0.09
Cϵ1 = 1.44
Cϵ2 = 1.92
σ = 1.0
σk = 1.0
σϵ = 1.3

# Grandezas características
Uc = U_jet
Lc = φ_jet
Cc = 1

# Grupos adimensionais
Re = ρ*U_jet*φ_jet/μ
# Re = 300
# const Sc = fluidKinematicViscosity/molecularDiffC
Sc = 1.0
# const Fr = cVelocity/sqrt(9.81*cLength)
# const Ri = (1/Fr^2)*cConcentration*beta
Ri = 0.1


U = VectorValue(0.33) # Valor adimensional calculado do perfil do jato no código 2d
Us = VectorValue(-0.36)
W = U + Us
ĝ = VectorValue(-1.0)

nu_T(k, ϵ) = Cμ*k*k/max(ϵ, 0.01)

resC(t, C, k, ϵ, v1) =
  ∫( v1 * ∂t(C) )dΩ +
  ∫( v1 * W⋅∇(C) )dΩ +
  ∫( (1/Re/Sc + (nu_T∘(k,ϵ)) / σ) * ∇(v1)⋅∇(C) )dΩ

resk(t, k, ϵ, v2) =
  ∫( v2 * ∂t(k) )dΩ +
  ∫( v2 * U⋅∇(k) )dΩ +
  ∫( (nu_T∘(k,ϵ)) / σk * ∇(v2)⋅∇(k) )dΩ +
  ∫( v2 * ϵ * (tanh∘(10.0*k/k_ΓD)) )dΩ +
  ∫( v2 * (Ri * (nu_T∘(k,ϵ))/σ) * (tanh∘(10.0*k/k_ΓD)) )dΩ

resEpsilon(t, k, ϵ, v3) =
  ∫( v3 * ∂t(ϵ) )dΩ +
  ∫( v3 * U⋅∇(ϵ) )dΩ + 
  ∫( (nu_T∘(k, ϵ))/σϵ * ∇(v3)⋅∇(ϵ) )dΩ +
  ∫( v3 * Cϵ2 * ϵ * ϵ/k * (tanh∘(10.0*ϵ/ϵ_ΓD)) )dΩ

res(t, (C, k, ϵ), (v1, v2, v3)) =  
  resC(t, C, k, ϵ, v1) +
  resk(t, k, ϵ, v2) +
  resEpsilon(t, k, ϵ, v3)

# -------------------------------------------------------------------
#                      Formulação de MEF
# -------------------------------------------------------------------

n = 150
model_Ω = CartesianDiscreteModel((0, H_uasb/Lc), (n))
Ω = Interior(model_Ω)
dΩ = Measure(Ω, 2)

# Tipo de elemento, tipo de função de forma, grau da função de forma
refFE = ReferenceFE(lagrangian, Float64, 1)

# Espaços de funções de teste
testSpace_C = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="tag_1")
testSpace_k = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="tag_1")
testSpace_epsilon = TestFESpace(model_Ω, refFE, conformity=:H1, dirichlet_tags="tag_1")


# --- Condições de contorno e condições iniciais ---
#=
06/09/2025
- Condições de contorno na saída agora são Neumann homogêneas
- Condições de contorno agora possuem inicialização suave no tempo
- Condições iniciais são pequenas, mas não nulas. Código divergiu com nulas
=#

α = 0.15
f(t) = 1-exp(-t/α)

bc_C(x, t::Real) = 1.0 * f(t)
bc_C(t::Real) = x -> bc_C(x, t)
ic_C(x, t::Real) = 0.01
ic_C(t::Real) = x -> ic_C(x, t)

I = 0.1
k_ΓD = (3/2)*(I*U_jet/Uc)^2 / 3
bc_k(x, t::Real) = k_ΓD * f(t)
bc_k(t::Real) = x -> bc_k(x, t)
ic_k(x, t::Real) = 0.01 * k_ΓD
ic_k(t::Real) = x -> ic_k(x, t)

γ = 0.07
ϵ_ΓD = Cμ^(3/4) * k_ΓD^(3/2) / (γ * φ_jet/Lc) / 3
bc_epsilon(x, t::Real) = ϵ_ΓD * f(t)
bc_epsilon(t::Real) = x -> bc_epsilon(x, t)
ic_epsilon(x, t::Real) = 0.01 * ϵ_ΓD
ic_epsilon(t::Real) = x -> ic_epsilon(x, t)

# Espaços de funções candidatas
trialSpace_C = TransientTrialFESpace(testSpace_C, bc_C)
trialSpace_k = TransientTrialFESpace(testSpace_k, bc_k)
trialSpace_epsilon = TransientTrialFESpace(testSpace_epsilon, bc_epsilon)

Y = MultiFieldFESpace([testSpace_C, testSpace_k, testSpace_epsilon])
X = TransientMultiFieldFESpace([trialSpace_C, trialSpace_k, trialSpace_epsilon])

op = TransientFEOperator(res, X, Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)

CFL = 0.75
dx = H_uasb/Lc/n
Δt = CFL*dx
θ = 1
ode_solver = ThetaMethod(nls,Δt,θ)

u_θ1 = interpolate_everywhere(ic_C(0), trialSpace_C(0.0))
u_θ2 = interpolate_everywhere(ic_k(0), trialSpace_k(0.0))
u_θ3 = interpolate_everywhere(ic_epsilon(0), trialSpace_epsilon(0.0))
u_θ = interpolate_everywhere([u_θ1, u_θ2, u_θ3], X(0.0))
t_θ = 0.0
T = 500.0

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
Reynolds number:    $Re
Richardson number:  $Ri
Number of jets:     $nJets
Settling velocity:  $Us
"""

path = joinpath(@__DIR__, "output")
if !isdir(path)
  mkdir(path)
end

ReStr = Int(round(Re, RoundDown))
path = joinpath(path,"Re$ReStr-nJets$nJets-Us$Us")
if isdir(path)
  rm(path, recursive=true)
end
mkdir(path)

path = joinpath(path, "case-description.txt")
if isfile(path)
  rm(path)
end
touch(path)
open(path, "w") do file
    write(file, caseDesc)
end

#=------------------------
Solução e escrita dos resultados
------------------------=#

ch, kh, epsilonh = u_θ
it = 0

writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets/uasbcp$it.vtu",cellfields=["ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh, "nu_T"=>nu_T∘(kh, epsilonh)])

totalIts = Int(round(T/Δt, RoundUp))
for (t,u_h) in u_ht
  global it
  local ch, kh, epsilonh
  it = it + 1
  println("Iterations completed: $it/$totalIts")
  ch, kh, epsilonh = u_h 
  if (mod(it,1)==0)
    writevtk(Ω,(@__DIR__)*"/output/Re$ReStr-nJets$nJets/uasbcp$it.vtu",cellfields=["ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh, "nu_T"=>nu_T∘(kh, epsilonh)])
  end
end