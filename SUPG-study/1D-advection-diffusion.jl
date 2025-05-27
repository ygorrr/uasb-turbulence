#= Código para estudar a estabilização do transporte de grandeza escalar u usando a técnica de SUPG. =#

using Gridap

n = 5
domain = (0, 1)
partition = (n)

#  --- Condições de contorno ---
u_0 = 0.0
u_1 = 0.0

ν = 1.0
advecVel = 1.0
source = 1.0

model = CartesianDiscreteModel(domain,partition)

if !isdir((@__DIR__)*"/results")
  mkdir((@__DIR__)*"/results")
end

writevtk(model, (@__DIR__)*"/results/model")

order = 1
reffeᵤ = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["tag_1","tag_2"])

U = TrialFESpace(V, [u_0, u_1])

degree = 2*order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

c(u,v) = ∫( v⋅(∇(u))'⋅advecVel )dΩ
# dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

a(u,v) = ∫( ν*∇(v)⊙∇(u) )dΩ

# c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
# dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

res(u,v) = a(u,v) + c(u,v) - ∫( v*source )dΩ
# jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)

op = AffineFEOperator(res,U,V)

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
)
  
solver = FESolver(nls)

uh, ph = solve(solver,op)

writevtk(Ωₕ,(@__DIR__)*"/ins-results",cellfields=["uh"=>uh,"ph"=>ph])