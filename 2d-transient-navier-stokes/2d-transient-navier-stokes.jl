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

a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)
jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)

op = TransientFEOperator(res,jac,X,Y)

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
)

Δt = 0.05
θ = 0.5
solver = ThetaMethod(nls, Δt, θ)
# solver = FESolver(nls)

t0, tF = 0.0, 10.0
uh0 = interpolate_everywhere(g(t0), Ug(t0))
uh, ph = solve(solver, op, t0, tF, uh0)
# uh, ph = solve(solver,op)

writevtk(Ωₕ,(@__DIR__)*"/ins-results",cellfields=["uh"=>uh,"ph"=>ph])