#=
Código para escoamento laminar em canal.
Adaptado do tutorial de Navier-Stokes incompressível com escoamento em cavidade.

https://gridap.github.io/Tutorials/dev/pages/t008_inc_navier_stokes/#inc_navier_stokes.jl-1
=#

using Gridap
using LineSearches: BackTracking

#  --- Propriedades físicas e parâmetros adimensionais ---
const Re = 10.0

# --- Parâmetros de domínio e malha ---
n = 50
domainLength = 10.0
domainHeight = 1.0
domain = (0, domainLength, 0, domainHeight)
partition = (5*n, n)

#  --- Condições de contorno ---
inletVelocity(t) = x -> VectorValue(10.0, 0)
wallVelocity(t) = x -> VectorValue(0, 0)
outletPressure(t) = x -> 0

#  --- Condição inicial ---
u0(t) = x -> VectorValue(0.0, 0.0)
p0(t) = x -> 0.0

model = CartesianDiscreteModel(domain,partition)
# writevtk(model, (@__DIR__)*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"inlet",[7,])
add_tag_from_tags!(labels,"outlet",[8,])
add_tag_from_tags!(labels,"walls",[1,2,3,4,5,6])

order = 2
reffeᵤ = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, labels=labels, dirichlet_tags=["inlet", "walls"])

reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, dirichlet_tags=["outlet"])

U = TransientTrialFESpace(V, [inletVelocity, wallVelocity])
P = TransientTrialFESpace(Q, [outletPressure])

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = 2*order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

conv(u,∇u) = Re*(∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

m(t, dtu, v) = ∫( v⋅dtu )dΩ

a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

res(t, (u,p), (v,q)) = m(t, u, v) + a((u,p),(v,q)) + c(u,v)
jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)

# op = TransientFEOperator(res,jac,X,Y)
op = TransientFEOperator(res,X,Y)

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
)

Δt = 0.05
θ = 0.5
solver = ThetaMethod(nls, Δt, θ)
# solver = FESolver(nls)

t0, tF = 0.0, 1.0
# uh0 = interpolate_everywhere(u0(t0), U(t0))
# ph0 = interpolate_everywhere(p0(t0), P(t0))
X0 = interpolate_everywhere([u0(0), p0(0)], X(t0))
u_ht = solve(solver, op, t0, tF, X0)
# uh, ph = solve(solver, op, t0, tF, X0)
# uh, ph = solve(solver,op)

it = 0
writevtk(Ωₕ,(@__DIR__)*"/ins-results$it.vtu",cellfields=["uh"=>uh,"ph"=>ph])

# for (u_h, t) in u_ht
#   global it
#   it += 1
#   uh, ph = u_h
#   writevtk(Ωₕ,(@__DIR__)*"/ins-results$it.vtu",cellfields=["uh"=>uh,"ph"=>ph])
# end

if !isdir((@__DIR__)*"/tmp")
  mkdir((@__DIR__)*"/tmp")
end

createpvd("results") do pvd
  uh0, ph0 = X0
  pvd[0] = createvtk(Ωₕ, (@__DIR__)*"/tmp/results_0" * ".vtu", cellfields=["u" => uh0, "p" => ph0])
  for (tn, u_hn) in uh
    uh, ph = u_hn
    pvd[tn] = createvtk(Ωₕ, (@__DIR__)*"/tmp/results_$tn" * ".vtu", cellfields=["u" => uhn, "p" => phn])
  end
end