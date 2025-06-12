# filepath: advection-diffusion-project/advection-diffusion-project/src/main.jl
using Gridap

# Define the domain and partition
n = 20
domain = (0.0, 1.0)
partition = (n,)

# Boundary conditions
u_0 = 0.0  # Dirichlet boundary condition at x=0
u_1 = 0.0  # Dirichlet boundary condition at x=1

# Parameters
D = 1.0  # Diffusion coefficient
advecVel = VectorValue(100.0)  # Advection velocity
source = 1.0  # Source term

# Create the discrete model
model = CartesianDiscreteModel(domain, partition)

# Create results directory if it doesn't exist
if !isdir((@__DIR__)*"/results")
    mkdir((@__DIR__)*"/results")
end

# Write the model to a VTK file
writevtk(model, (@__DIR__)*"/results/model")

# Define finite element spaces
order = 1
reffeᵤ = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["tag_1", "tag_2"])
U = TrialFESpace(V, [u_0, u_1])

# Define measures
degree = 2 * order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ, degree)

# Define stabilization parameter (SUPG parameter τ)
function compute_tau(h, advecVel, D)
    return h / (2 * norm(advecVel)) * (coth(Pe) - 1 / Pe)
end

h = 1.0 / n  # Element size (assume uniform mesh)

Pe = norm(advecVel) * h / (2 * D)  # Péclet number
@info "Péclet number: " Pe

τ = compute_tau(h, advecVel, D)

# Define variational forms with SUPG stabilization
a(u, v) = ∫(v * (advecVel ⋅ ∇(u)) + D * ∇(v) ⋅ ∇(u))dΩ +
∫(τ * (advecVel ⋅ ∇(v)) * (advecVel ⋅ ∇(u)))dΩ

b(v) = ∫(v * source)dΩ +
∫(τ * (advecVel ⋅ ∇(v)) * source)dΩ

# Define the operator
op = AffineFEOperator(a, b, U, V)

# Solve the problem
uh = solve(op)

# Write the solution to a VTK file
writevtk(Ωₕ, (@__DIR__)*"/results/solution", cellfields=["uh" => uh])