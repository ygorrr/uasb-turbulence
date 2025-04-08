using Gridap
using LineSearches: BackTracking

n = 100
domain = (0,20,0,1)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)
# writevtk(model, (@__DIR__)*"/model")

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"inlet",[7,])
add_tag_from_tags!(labels,"outlet",[8,])
add_tag_from_tags!(labels,"walls",[1,2,3,4,5,7])

order = 2
reffeᵤ = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, labels=labels, dirichlet_tags=["inlet", "walls"])

reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean, dirichlet_tags=["outlet"])

# With the options `:Lagrangian`, `space=:P`, `valuetype=Float64`, and `order=order-1`, we select the local polynomial space $P_{k-1}(T)$ on the cells $T\in\mathcal{T}$. With the symbol `space=:P` we specifically chose a local Lagrangian interpolation of type "P". Without using `space=:P`, would lead to a local Lagrangian of type "Q" since this is the default for quadrilateral or hexahedral elements. On the other hand, `constraint=:zeromean` leads to a FE space, whose functions are constrained to have mean value equal to zero, which is just what we need for the pressure space. With these objects, we build the trial multi-field FE spaces

inletVelocity = x -> VectorValue(1.0, 0.0)
wallVelocity = x -> VectorValue(0.0, 0.0)
outletPressure = x -> 0.0
U = TransientTrialFESpace(V, [inletVelocity, wallVelocity])
P = TransientTrialFESpace(Q, [outletPressure])

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

# ## Triangulation and integration quadrature
#
# From the discrete model we can define the triangulation and integration measure

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

# ## Nonlinear weak form
#
# The different terms of the nonlinear weak form for this example are defined following an approach similar to the one discussed for the $p$-Laplacian equation, but this time using the notation for multi-field problems.

const Re = 10.0
conv(u,∇u) = Re*(∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

# The bilinear form reads
a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ

# The nonlinear term and its Jacobian are given by
c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

# Finally, the Navier-Stokes weak form residual and Jacobian can be defined as
res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)
jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)

# With the functions `res`, and `jac` representing the weak residual and the Jacobian, we build the nonlinear FE problem:
op = FEOperator(res,jac,X,Y)

# ## Nonlinear solver phase
#
# To finally solve the problem, we consider the same nonlinear solver as previously considered for the  $p$-Laplacian equation.

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

# In this example, we solve the problem without providing an initial guess (a default one equal to zero will be generated internally)

uh, ph = solve(solver,op)

# Finally, we write the results for visualization (see next figure).

writevtk(Ωₕ,"ins-results",cellfields=["uh"=>uh,"ph"=>ph])

# ![](../assets/inc_navier_stokes/ins_solution.png)
#
#  ## References
#
#  [1] H. C. Elman, D. J. Silvester, and A. J. Wathen. *Finite elements and fast iterative solvers: with applications in incompressible fluid dynamics*. Oxford University Press, 2005.