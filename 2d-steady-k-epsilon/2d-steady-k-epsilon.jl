using Gridap

g = VectorValue(0.0, -9.81)  # Aceleração da gravidade
nu = 1.0e-3  # Viscosidade cinemática

#  --- Formas compactas ---
a(u,v,D) = ∫( D*∇(u)⊙∇(v) )dΩ
c(u,v,a) = ∫( v ⋅ (∇(u) * a) )dΩ
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
nuEddy(u3, u4) = Cμ * u3^2 / max(u4, minVal)

# Equação de Navier-Stokes 
resNS(u1, u2, u3, u4, v1) = 
  c(u1, v1, u1) +
  a(u1, v1, nu + nuEddy(u3, u4)) + 
  ∫( v1 ⋅ (∇(u2) / ρ - g) )dΩ

# Equação da continuidade
resCont(u1, v2) =
  ∫( v2 * (∇ ⋅ u1) )dΩ

# Equação de k
resk(u1, u3, u4, v3) = 
  c(u3, v3, u1) +
  a(u3, v3, nuEddy(u3, u4) / σk) -
  ∫( v3 * (source_u3 - u4) )dΩ

# Equação de epsilon
resEpsilon(u1, u3, u4, v4) =
  c(u4, v4, u1) +
  a(u4, v4, nuEddy(u3, u4) / σε) -
  ∫( v4 * (Cϵ1 * source_u3 * u4 / max(u3, minVal)) - Cϵ2 * u4^2 / max(u3, minVal) )dΩ

# --- Malha ---
L = 1.0
partition = (50, 50)
domain = (0.0, L, 0.0, L)

model_base = CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)

if !isdir((@__DIR__)*"/results")
  mkdir((@__DIR__)*"/results")
end

writevtk(model, (@__DIR__)*"/results/model")

tags = Dict(
  :inlet  => x -> isapprox(x[2], 0.0; atol=1e-8) && 0.4 ≤ x[1] ≤ 0.6,
  :walls  => x -> isapprox(x[2], 0.0; atol=1e-8) && !(0.4 ≤ x[1] ≤ 0.6),
  :sides  => x -> isapprox(x[1], 0.0; atol=1e-8) || isapprox(x[1], L; atol=1e-8),
  :top    => x -> isapprox(x[2], L; atol=1e-8),
  :anchor => x -> isapprox(x[1], 0.5; atol=1e-8) && isapprox(x[2], 0.0; atol=1e-8)
)

# Cria o modelo diretamente com os rótulos de contorno aplicados
model = CartesianDiscreteModel(domain, partition; boundary_tags=tags)

#  --- Condições de contorno ---
# Velocidade do jato
Uref = 1.0
inlet_velocity(x) = VectorValue(0.0, Uref)  # Jato entrando na vertical

u1_dirichlet = DirichletBoundary(inlet_velocity, tags = ["inlet"])

# Ponto de ancoragem da pressão
anchor_tag = Dict(:anchor => x -> isapprox(x[1], 0.5; atol=1e-8) && isapprox(x[2], 0.0; atol=1e-8))
anchor_tags = generate_boundary_tags(model, anchor_tag)

# Função constante p = 0.0
pressure_ref(x) = 0.0
p_dirichlet = DirichletBoundary(pressure_ref, tags = ["anchor"])

# Condições de Dirichlet para k e epsilon
I = 0.05      # Intensidade de turbulência
ell = 0.05    # Comprimento de mistura (tamanho do jato)

k_inlet = 1.5 * (I * Uref)^2
ϵ_inlet = (Cμ)^(3/4) * k_inlet^(3/2) / ell

inlet_k(x) = k_inlet
inlet_ϵ(x) = ϵ_inlet

k_dirichlet = DirichletBoundary(inlet_k, tags = ["inlet"])
ϵ_dirichlet = DirichletBoundary(inlet_ϵ, tags = ["inlet"])

#  --- Espaços de funções ---
order = 1

# Espaço de funções de teste para a velocidade
reffe_u1 = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V_u1 = TestFESpace(model, reffe_u1, conformity=:H1, labels=labels, dirichlet_tags=["inlet"])

# Espaço de funções de teste para a pressão
reffe_u2 = ReferenceFE(lagrangian, Float64, order-1; space=:P)
V_u2 = TestFESpace(model, reffe_u2, conformity=:L2, dirichlet_tags=["inlet"])

# Espaço de funções de teste para k
reffe_u3 = ReferenceFE(lagrangian, Float64, order)
V_u3 = TestFESpace(model, reffe_u3, conformity=:H1, labels=labels, dirichlet_tags=["inlet"])

# Espaço de funções de teste para epsilon
reffe_u4 = ReferenceFE(lagrangian, Float64, order)
V_u4 = TestFESpace(model, reffe_u4, conformity=:H1, labels=labels, dirichlet_tags=["inlet"])

U_u1 = TrialFESpace(V_u1, u1_dirichlet)
U_u2 = TrialFESpace(V_u2, p_dirichlet)
U_u3 = TrialFESpace(V_u3, k_dirichlet)
U_u4 = TrialFESpace(V_u4, ϵ_dirichlet)

Y = MultiFieldFESpace([V_u1, V_u2, V_u3, V_u4])
X = MultiFieldFESpace([U_u1, U_u2, U_u3, U_u4])

degree = 2*order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

conv(u,∇u) = Re*(∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

m(t, dtu, v) = ∫( Re*v⋅dtu )dΩ

a(t, (u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ

c(t, u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(t, u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

res(t, (u,p), (v,q)) = m(t, u, v) + a(t,(u,p),(v,q)) + c(t,u,v)
jac(t,(u,p),(du,dp),(v,q)) = a(t,(du,dp),(v,q)) + dc(t,u,du,v)

# op = TransientFEOperator(res,jac,X,Y)
op = TransientFEOperator(res,X,Y)

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
)

Δt = 0.05
θ = 0.5
solver = ThetaMethod(nls, Δt, θ)
# solver = FESolver(nls)

t0, tF = 0.0, 1.5
# uh0 = interpolate_everywhere(u0(t0), U(t0))
# ph0 = interpolate_everywhere(p0(t0), P(t0))
X0 = interpolate_everywhere([u0(0), p0(0)], X(t0))
u_ht = solve(solver, op, t0, tF, X0)
# uh, ph = solve(solver, op, t0, tF, X0)
# uh, ph = solve(solver,op)

# it = 0
# writevtk(Ωₕ,(@__DIR__)*"/ins-results$it.vtu",cellfields=["uh"=>uh,"ph"=>ph])

# for (u_h, t) in u_ht
#   global it
#   it += 1
#   uh, ph = u_h
#   writevtk(Ωₕ,(@__DIR__)*"/ins-results$it.vtu",cellfields=["uh"=>uh,"ph"=>ph])
# end

if !isdir((@__DIR__)*"/tmp")
  mkdir((@__DIR__)*"/tmp")
end

createpvd((@__DIR__)*"/tmp/results") do pvd
  uh0, ph0 = X0
  pvd[0] = createvtk(Ωₕ, (@__DIR__)*"/tmp/results_0" * ".vtu", cellfields=["u" => uh0, "p" => ph0])
  for (tn, u_hn) in u_ht
    uh, ph = u_hn
    pvd[tn] = createvtk(Ωₕ, (@__DIR__)*"/tmp/results_$tn" * ".vtu", cellfields=["u" => uh, "p" => ph])
  end
end