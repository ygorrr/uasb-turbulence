using Gridap, GridapGmsh

function g(t, x)
  if x[1] == 0
    return 1
  else
    return 0
  end
end

g(t) = x -> g(t, x)

domain = (0, 10, 0, 1)
partition = (20, 20)
model = CartesianDiscreteModel(domain, partition)
# model = GmshDiscreteModel("rectangular-domain.msh")
writevtk(model, (@__DIR__)*"/model")
# writevtk(Ω, (@__DIR__)*"/model")

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)

V0 = TestFESpace(model, reffe, dirichlet_tags=["tag_1", "tag_3", "tag_7"])
Ug = TransientTrialFESpace(V0, g)

# ## Triangulation and quadrature

# As usual, we equip the model with an integration mesh and a measure

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

neumanntags = ["tag_2", "tag_4", "tag_5", "tag_6", "tag_8"]
Γ = BoundaryTriangulation(model, tags=neumanntags)
dΓ = Measure(Γ, degree)

# ## Weak form
# We define the thermal diffusivity $\alpha$ and the rate of external temperature generation $f$.

# α(t) = x -> 1 + sin(t) * (x[1]^2 + x[2]^2) / 4
α(t) = 1
# f(t) = x -> sin(t) * sinpi(x[1]) * sinpi(x[2])
f(t) = 0

# We are going to construct a transient linear FEOperator by providing the bilinear forms associated to $\partial_{t} u$ and $u$, as well as the forcing term. Note that they now receive time as an additional argument, and the time derivative operator is `∂t`.

m(t, dtu, v) = ∫(v * dtu)dΩ
a(t, u, v) = ∫(α(t) * ∇(v) ⋅ ∇(u))dΩ
l(t, v) = ∫(v * f(t))dΩ
op = TransientLinearFEOperator((a, m), l, Ug, V0)

# In our case, the mass term ($m(t, \cdot, \cdot)$) is constant in time. We can take advantage of that to save some computational effort, and indicate it to Gridap as follows
# op_opt = TransientLinearFEOperator((a, m), l, Ug, V0, constant_forms=(true, false))

# If the stiffness term ($a(t, \cdot, \cdot)$) had been constant in time, we could have set `constant_forms=(true, true)`.

# ## Transient solver

# Once we have defined the FE operator, we proceed with the definition of the ODE solver, i.e. the scheme that will be used for the integration in time. In this tutorial, we use the `ThetaMethod` with $\theta = 1/2$, resulting in a second-order scheme. The `ThetaMethod` function receives a solver for systems of equations, the time step size $\Delta t$ (constant) and the value of $\theta \in [0, 1]$. Since the ODE is linear the systems of equation that will arise in the time-marching scheme will be linear so we can provide `ThetaMethod` with a linear solver.

ls = LUSolver()
Δt = 0.05
θ = 0.5
solver = ThetaMethod(ls, Δt, θ)

# Gridap also implements explicit and diagonally-implicit Runge-Kutta schemes. One can access the full list of available Butcher tableaus through the exported constant `available_tableaus`. There are also constructors for explicit 2- and 3-stage schemes: `EXRK22(α)` and `EXRK33(α, β)`, `EXRK33_1(α)`, `EXRK33_2(α)` respectively, and diagonally-implicit 1- and 2-stage schemes: `SDIRK11(α)`, `SDIRK12()`, `SDIRK22(α, β, γ)`, `SDIRK23(λ)`. See the documentation of [Runge-Kutta schemes in Gridap](https://gridap.github.io/Gridap.jl/dev/ODEs/#Runge-Kutta) for a description of the corresponding tableaus. For example, one could have chosen a two-stage singly-diagonally-implicit scheme (of order 2) as follows.
# tableau = :SDIRK_2_2
# solver_rk = RungeKutta(ls, ls, Δt, tableau)

# Let $t_{F} > t_{0}$ be a final time, until when we want to evolve the problem. We define the solution using the `solve` function, giving the ODE solver, the transient FE operator, the initial and final times, and the initial solution. To construct the initial condition we interpolate the initial function $u_{0}$ onto the FE space $U_{g}$ at the initial time. In our case, $u_{0}$ is simply $g(t_{0})$.

t0, tF = 0.0, 10.0
uh0 = interpolate_everywhere(g(t0), Ug(t0))
uh = solve(solver, op, t0, tF, uh0)

# ## Postprocessing

# We highlight that `uh` is an iterable function and the result at each time step is only computed lazily when iterating over it. We can post-process the results and generate the corresponding `vtk` files using the `createpvd` and `createvtk` functions. The former will create a `.pvd` file with the collection of `.vtu` files saved at each time step by `createvtk`. The computation of the problem solutions will be triggered in the following loop:

tmpStr = (@__DIR__)*"/tmp"

if !isdir(tmpStr)
  mkdir(tmpStr)
end

createpvd(tmpStr * "/results") do pvd
  pvd[0] = createvtk(Ω, tmpStr * "/results_0" * ".vtu", cellfields=["u" => uh0])
  for (tn, uhn) in uh
    pvd[tn] = createvtk(Ω, tmpStr * "/results_$tn" * ".vtu", cellfields=["u" => uhn])
  end
end

# ![](../assets/transient_linear/result.gif)