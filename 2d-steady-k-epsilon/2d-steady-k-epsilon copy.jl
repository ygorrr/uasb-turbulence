using Gridap, GridapGmsh
using LineSearches: BackTracking

g = VectorValue(0.0, -9.81)  # Aceleração da gravidade
nu = 1.0e-1  # Viscosidade cinemática
ρ = 1.0  # Densidade do fluido

#  --- Formas compactas ---
a(u,v,D) = ∫( D*∇(u)⊙∇(v) )dΩ
c(u,v,a) = ∫( v ⋅ ( (∇(u)') ⋅ a) )dΩ
m(u,v) = ∫( v ⋅ u )dΩ

#  --- Resíduos ---
#=
u1 -> vetor velocidade média
u2 -> pressão média
u3 -> k
u4 -> epsilon
=#

minVal = 1e-8 # Constante de proteção contra divisão por zero

Cμ = 0.09  # Constante de turbulência
Cϵ1 = 1.44  # Constante de produção de epsilon
Cϵ2 = 1.92  # Constante de destruição de epsilon
σk = 1.0  # Constante de difusão de k
σϵ = 1.3  # Constante de difusão de epsilon
nuEddy(u3, u4) = 0*Cμ * u3 * u3 / max(u4, minVal)
kProduction(∇u1, u3, u4) = - (2/3) * u3 + nuEddy(u3, u4) * (∇u1 + ∇u1') ⊙ ∇u1
epsilonProduction(∇u1, u3, u4) = Cϵ1 * kProduction(∇u1, u3, u4) * u4 / max(u3, minVal)
epsilonDestruction(u3, u4) = Cϵ2 * u4 * u4 / max(u3, minVal)

# Equação de Navier-Stokes 
resNS(u1, u2, u3, u4, v1) = 
  c(u1, v1, u1) +
  ∫( (nu + nuEddy∘(u3,u4)) * ∇(u1)⊙∇(v1) )dΩ + 
  ∫( v1 ⋅ (∇(u2) / ρ - g) )dΩ

# Equação da continuidade
resCont(u1, v2) =
  ∫( v2 * (∇ ⋅ u1) )dΩ

# Equação de k
resk(u1, u3, u4, v3) = 0
  # c(u3, v3, u1) +
  # ∫( (nu + nuEddy∘(u3,u4)) / σk * ∇(u3)⊙∇(v3) )dΩ -
  # ∫( v3 * (kProduction∘(∇(u1),u3,u4) - u4) )dΩ

# Equação de epsilon
resEpsilon(u1, u3, u4, v4) = 0
  # c(u4, v4, u1) +
  # ∫( (nu + nuEddy∘(u3,u4)) / σϵ * ∇(u4)⊙∇(v4) )dΩ -
  # ∫( v4 * (epsilonProduction∘(∇(u1),u3,u4) - epsilonDestruction∘(u3,u4)) )dΩ

res((u1, u2, u3, u4), (v1, v2, v3, v4)) =
  resNS(u1, u2, u3, u4, v1) +
  resCont(u1, v2) +
  resk(u1, u3, u4, v3) +
  resEpsilon(u1, u3, u4, v4)

# --- Malha ---
n = 50
Ly = 3
domain = (0,1,0,Ly)
partition = (n,Ly*n)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[6,])
add_tag_from_tags!(labels,"bottom",[5,])

# msh_file = (@__DIR__)*"/mesh/jato.msh"
# model = GmshDiscreteModel(msh_file)

if !isdir((@__DIR__)*"/results")
  mkdir((@__DIR__)*"/results")
end

writevtk(model, (@__DIR__)*"/results/model")

#  --- Condições de contorno ---
# Velocidade do jato
Uref = 1.0
# inlet_velocity(x) = VectorValue(0.0, Uref)  # Jato entrando na vertical

inlet_velocity(x) = VectorValue(0, 1.0*(1.0-sqrt(x[2]/Ly))*(((0.5*(1.0-(cos(1*2*pi*x[1])))+0.0*rand()) > 0.89 ? 1.0 : 0.0 )) + sqrt(x[2]/Ly))


# using Plots
# x=collect(range(0,1,100))
# g2(x) =  VectorValue(0.0, 1.0*(((0.5*(1.0-(cos(1*2*pi*x[1])))+0.0*rand()) > 0.89 ? 1.0 : 0.0 )))

# plot(x, g2.(x), label="Velocidade do jato", xlabel="x", ylabel="Velocidade (m/s)", title="Perfil de Velocidade do Jato")

# Função constante p = 0.0
# inlet_pressure(x) = 0.0
outlet_pressure(x) = 0.0

# Condições de Dirichlet para k e epsilon
I = 0.05      # Intensidade de turbulência
ell = 0.05    # Comprimento de mistura (tamanho do jato)
k_inlet = 1.5 * (I * Uref)^2
ϵ_inlet = (Cμ)^(3/4) * k_inlet^(3/2) / ell

inlet_k(x) = k_inlet
inlet_ϵ(x) = ϵ_inlet

#  --- Espaços de funções ---
order = 2

# Espaço de funções de teste para a velocidade
reffe_u1 = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V_u1 = TestFESpace(model, reffe_u1, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

# Espaço de funções de teste para a pressão
reffe_u2 = ReferenceFE(lagrangian, Float64, order-1)
V_u2 = TestFESpace(model, reffe_u2, conformity=:L2, dirichlet_tags=["top"])

# Espaço de funções de teste para k
reffe_u3 = ReferenceFE(lagrangian, Float64, order)
V_u3 = TestFESpace(model, reffe_u3, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

# Espaço de funções de teste para epsilon
reffe_u4 = ReferenceFE(lagrangian, Float64, order)
V_u4 = TestFESpace(model, reffe_u4, conformity=:H1, labels=labels, dirichlet_tags=["bottom"])

U_u1 = TrialFESpace(V_u1, inlet_velocity)
U_u2 = TrialFESpace(V_u2, outlet_pressure)
U_u3 = TrialFESpace(V_u3, inlet_k)
U_u4 = TrialFESpace(V_u4, inlet_ϵ)

Y = MultiFieldFESpace([V_u1, V_u2, V_u3, V_u4])
X = MultiFieldFESpace([U_u1, U_u2, U_u3, U_u4])

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

op = FEOperator(res,X,Y)

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
)

solver = FESolver(nls)

uh, ph, kh, epsilonh = solve(solver,op)

writevtk(Ω,(@__DIR__)*"/results",cellfields=["uh"=>uh,"ph"=>ph,"kh"=>kh,"epsilonh"=>epsilonh])