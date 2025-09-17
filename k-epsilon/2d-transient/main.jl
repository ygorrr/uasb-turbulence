using Gridap
using GridapGmsh
using LineSearches: BackTracking
using Plots
using Dates

#--- Propriedades físicas ---
β = 0.1 # Coeficiente de dilatação volumétrica
ρ = 1000.0  # Densidade do fluido
μ = 1.0e-3 # Viscosidade dinâmica da água em Pa.s
g = 9.81 # Aceleração da gravidade em m/s²
ĝ = VectorValue(0.0, -1.0)  # Vetor aceleração unitário

#--- Constantes do modelo k-epsilon ---
Cμ = 0.09
Cϵ1 = 1.44
Cϵ2 = 1.92
σk = 1.0
σϵ = 1.3
σ = 1.0

#--- Parâmetros geométricos e operacionais do UASB ---
φ_uasb = 1.8              # m
H_uasb = 2         # m
A_uasb = π*(φ_uasb/2)^2   # m²
Q = 0.5/1000    # m³/s
HDT = H_uasb * A_uasb / Q
HDT_horas = HDT / 3600
U_uasb = Q / A_uasb       # m/s

φ_jet = 0.2               # m
A_jet = π*(φ_jet/2)^2     # m²
xJets = [0.0]     # Posições dos jatos ao longo de x
nJets = length(xJets)     # Número de jatos
U_jet = Q/(9*A_jet)   # m/s

Us = VectorValue(0.0, -0.1) # Velocidade de decantação de partículas
minVal = 1e-2 # Constante de proteção contra divisão por zero

#--- Grandezas características do modelo ---
Uc = U_jet
Lc = φ_jet
Cc = 1.0

#--- Números adimensionais do modelo ---
# Número de Reynolds do UASB
Re_uasb = ρ * U_uasb * φ_uasb / μ

# Número de Reynolds do jato
Re_jet = ρ * U_jet * φ_jet / μ

# Número de Reynolds do modelo
Re = ρ * Uc * Lc / μ
Sc = 1.0
# Ri = β*g*(Cc/Lc)/(Uc/Lc)^2
Ri = 1e-1

#=
Intensidade de turbulência. Valores típicos:
Jato livre: 0.05 ~ 0.1
Escoamento altamente turbulento: até 0.2
I é definida como a razão entre a velocidade de flutuação e a velocidade média do jato.
I = u' / Uref
=#

Uref = U_jet / Uc
I = 0.1 
k_jet = (3/2) * (I * Uref)^2
k_jet = 0.015 # Forçar k para 0.015

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
ε_jet = Cμ^(3/4) * k_jet^(3/2) / (α * φ_jet/Lc)

#=
Estimativa baseada em argumentos de escala:
  νₜ ~ U_jet * φ_jet
Essa estimativa vem de argumentos de escala, supondo que as flutuações da velocidade e comprimento de mistura sejam proporcionais à velocidade média do jato e ao diâmetro do jato, ou seja, νₜ ~ u'ℓ ~ U_jet * φ_jet.
=#
# ε_jet = Cμ * k_jet^2 / (Uref * φ_jet/Lc)

#--- Condições iniciais e de contorno ---
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
    velocityY += amplitude * unitJet(x[1] - xJet / Lc)
  end
  return velocityY
end

function scalarProfile(x, t::Real)
  global xJets
  local scalarVal = 0.0
  for xJet in xJets
    scalarVal += unitJet(x[1] - xJet / Lc)
  end
  return scalarVal
end

#=
  Visualização do perfil de velocidade do jato.
  Melhor comentar esse bloco begin/end se for rodar no cluster.
=#
begin
  jetContour = ([-1, 1] * φ_jet / 2 .+ xJets[1]) / Lc
  xRange = -φ_uasb/Lc/2:0.001:φ_uasb/Lc/2
  yComponents = []
  for x in xRange
    append!(yComponents, jetProfile(x, 0))
  end

  plot(xRange, yComponents, xlabel="x", ylabel="Velocidade em y", title="Perfil de velocidade de jato", label="jetProfile")

  scatter!(jetContour, [0.0, 0.0], label="jetContour")
end

# Função utilitária para calcular valor médio de uma função f no intervalo [x_min, x_max]
function calcMeanValue(x_min, x_max, f)
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

function W(U)
  global Us
  return U + Us
end

function nu_T(k, ε)
  global Cμ, minVal
  return Cμ * k * k / max(ε, minVal)
end

function kProd(∇U, k,ε)
  return nu_T(k, ε) * ((∇U + ∇U') ⊙ ∇U)
end

function εProd(∇U, k, ε)
  global Cϵ1, minVal
  return Cϵ1 * kProd(∇U, k, ε) * max(ε, minVal) / max(k, minVal)  
end

function εDest(k, ε)
  global Cϵ2, minVal
  return Cϵ2 * ε * ε / max(k, minVal)  
end

#--- Gridap ---
n = 50
Lx = φ_uasb / Lc
Ly = H_uasb / Lc
domain = (-Lx/2, Lx/2, 0, Ly)
partition = (n, 2*n)
# model = CartesianDiscreteModel(domain, partition; isperiodic=(true,false))
model = GmshDiscreteModel(joinpath(@__DIR__,"uasb.msh"); isperiodic=(true,false))

# labels = get_face_labeling(model)
# add_tag_from_tags!(labels, "top", [6,])
# add_tag_from_tags!(labels, "bottom", [5,])

order = 2

reffe_U = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
reffe_P = ReferenceFE(lagrangian, Float64, order-1; space=:P)
reffe_C = ReferenceFE(lagrangian, Float64, order)
reffe_k = ReferenceFE(lagrangian, Float64, order)
reffe_ε = ReferenceFE(lagrangian, Float64, order)

V_U = TestFESpace(model, reffe_U, conformity=:H1, dirichlet_tags=["bottom"])
V_P = TestFESpace(model, reffe_P, conformity=:L2, dirichlet_tags=["top"])
V_C = TestFESpace(model, reffe_C, conformity=:H1, dirichlet_tags=["bottom"])
V_k = TestFESpace(model, reffe_k, conformity=:H1, dirichlet_tags=["bottom"])
V_ε = TestFESpace(model, reffe_ε, conformity=:H1, dirichlet_tags=["bottom"])

BC_U(x, t::Real) = VectorValue(0, jetProfile(x,t))
BC_U(t::Real) = x -> BC_U(x,t)
IC_U(x, t::Real) = exp(-5.0 * x[2] / Ly) * VectorValue(0, jetProfile(x,t))
IC_U(t::Real) = x -> IC_U(x,t)

BC_P(x, t::Real) = 0.0
BC_P(t::Real) = x -> BC_P(x,t)
IC_P(x, t::Real) = 0.0
IC_P(t::Real) = x -> IC_P(x,t)

BC_C(x, t::Real) = scalarProfile(x,t)
BC_C(t::Real) = x -> BC_C(x,t)
IC_C(x, t::Real) = 0.0
IC_C(t::Real) = x -> IC_C(x,t)

BC_k(x, t::Real) = k_jet * scalarProfile(x,t)
BC_k(t::Real) = x -> BC_k(x,t)
IC_k(x, t::Real) = k_jet * exp(-5.0 * x[2] / Ly) * scalarProfile(x,t)
IC_k(t::Real) = x -> IC_k(x,t)

BC_ε(x, t::Real) = ε_jet * scalarProfile(x,t)
BC_ε(t::Real) = x -> BC_ε(x,t)
IC_ε(x, t::Real) = ε_jet * exp(-5.0 * x[2] / Ly) * scalarProfile(x,t)
IC_ε(t::Real) = x -> IC_ε(x,t)

S_U = TransientTrialFESpace(V_U, [BC_U])
S_P = TransientTrialFESpace(V_P, [BC_P])
S_C = TransientTrialFESpace(V_C, [BC_C])
S_k = TransientTrialFESpace(V_k, [BC_k])
S_ε = TransientTrialFESpace(V_ε, [BC_ε])

Y = MultiFieldFESpace([V_U, V_P, V_C, V_k, V_ε])
X = TransientMultiFieldFESpace([S_U, S_P, S_C, S_k, S_ε])

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

#--- Resíduos ---
# Equação de Navier-Stokes
# Experimentando a função de parte simétrica de tensor ε(U) = ∇(U) + (∇(U)') / 2
resNS(t, U, P, C, k, ε, v1) = 
  ∫( v1 ⋅ ∂t(U) )dΩ +
  ∫( v1 ⋅ (∇(U)' ⋅ U) )dΩ +
  ∫( (1/Re + (nu_T∘(k,ε))) * (∇(v1) ⊙ (∇(U) + (∇(U))')) )dΩ - 
  ∫( (∇ ⋅ v1) * P )dΩ -
  ∫( Ri * (v1 ⋅ ĝ) * C )dΩ

# Equação da continuidade
resCont(t, U, v2) =
  ∫( v2 * (∇ ⋅ U) )dΩ

# Equação de transporte de partículas
resC(t, U, C, k, ε, v3) =
  ∫( v3 * ∂t(C) )dΩ +
  ∫( v3 ⋅ inner(W(U), ∇(C)) )dΩ +
  ∫( (1/Re * 1/Sc + (nu_T∘(k,ε))/σ) * ∇(C)⊙∇(v3) )dΩ

# Equação de k
resk(t, U, C, k, ε, v4) = 
  ∫( v4 * ∂t(k) )dΩ +
  ∫( v4 ⋅ (∇(k)' ⋅ U) )dΩ +
  ∫( (nu_T∘(k,ε)) / σk * ∇(v4)⊙∇(k) )dΩ -
  ∫( v4 * (kProd∘(∇(U),k,ε)) )dΩ +
  ∫( v4 * ε * (tanh∘(10.0*k/k_jet)) )dΩ +
  ∫( v4 * Ri * (nu_T∘(k,ε)) / σ * (ĝ ⋅ ∇(C)) * (tanh∘(10.0*k/k_jet)) )dΩ

# Equação de epsilon
resEpsilon(t, U, k, ε, v5) =
  ∫( v5 * ∂t(ε) )dΩ +
  ∫( v5 ⋅ (∇(ε)' ⋅ U) )dΩ +
  ∫( (nu_T∘(k,ε)) / σϵ * ∇(ε)⊙∇(v5) )dΩ -
  ∫( v5 * (εProd∘(∇(U),k,ε)) )dΩ +
  ∫( v5 * (εDest∘(k,ε)) * (tanh∘(10.0*ε/ε_jet)) )dΩ

res(t, (U, P, C, k, ε), (v1, v2, v3, v4, v5)) =
  resNS(t, U, P, C, k, ε, v1) +
  resCont(t, U, v2) +
  resC(t, U, C, k, ε, v3) +
  resk(t, U, C, k, ε, v4) +
  resEpsilon(t, U, k, ε, v5)

op = TransientFEOperator(res,X,Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=10)

CFL = 1.0/2.0
Δt = CFL*Lx/n
θ = 1

ode_solver = ThetaMethod(nls,Δt,θ)

U₀ = interpolate_everywhere([IC_U(0),IC_P(0),IC_C(0),IC_k(0),IC_ε(0)],X(0.0))
t₀ = 0.0
T = 100.0
uₕₜ = solve(ode_solver,op,t₀,T,U₀)
it = 0
uh, ph, ch, kh, epsilonh = U₀

#--- Diretório de saída ---
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

info = """
Spatial dimensions: 2
Transient:          true
Reynolds number:    $(round(Re, RoundDown))
Richardson number:  $Ri
Number of jets:     $nJets
Settling velocity:  $Us
Simulation time:    $T
Time step:          $Δt
"""

path = joinpath(@__DIR__, "output")
isdir(path) ? nothing : mkdir(path)

path = joinpath(path, timestamp)
isdir(path) ? nothing : mkdir(path)

touch(joinpath(path, "info.txt"))
open(joinpath(path, "info.txt"), "w") do file
    write(file, info)
end

#--- Solução e escrita dos resultados ---
writevtk(Ω, path*"/0.vtu",cellfields=["uh"=>uh,"ph"=>ph, "ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])

it = 1
totalIts = T/Δt
for (t,uₕ) in uₕₜ
  global it
  local uh,ph,ch,kh,epsilonh
  uh, ph, ch, kh, epsilonh = uₕ
  
  println("Iteration $it/$totalIts")
  
  if(mod(it,1)==0)
    writevtk(Ω, path*"/$it.vtu",cellfields=["uh"=>uh,"ph"=>ph,"ch"=>ch, "kh"=>kh, "epsilonh"=>epsilonh])
  end

  it = it + 1
end