# Advection-Diffusion Transport Problem

This project implements a numerical solution to the advection-diffusion transport problem using the Gridap package in Julia. The problem involves the transport of a scalar quantity \( u \) in a one-dimensional domain, subject to advection and diffusion processes, with a source term included. The solution is computed under homogeneous Dirichlet boundary conditions.

## Overview

The advection-diffusion equation can be expressed as:

\[
\frac{\partial u}{\partial t} + \mathbf{v} \cdot \nabla u = \nu \nabla^2 u + f
\]

where:
- \( u \) is the scalar quantity being transported,
- \( \mathbf{v} \) is the advection velocity,
- \( \nu \) is the diffusion coefficient,
- \( f \) is the source term.

## Project Structure

- `src/main.jl`: Contains the main implementation of the advection-diffusion transport problem.
- `Project.toml`: Project configuration file specifying dependencies.
- `Manifest.toml`: Locks the dependencies to specific versions for reproducibility.
- `README.md`: Documentation for the project.

## Setup Instructions

1. **Install Julia**: Ensure you have Julia installed on your machine. You can download it from [the official Julia website](https://julialang.org/downloads/).

2. **Clone the Repository**: Clone this repository to your local machine.

3. **Navigate to the Project Directory**:
   ```bash
   cd advection-diffusion-project
   ```

4. **Activate the Project**:
   ```julia
   using Pkg
   Pkg.activate(".")
   ```

5. **Install Dependencies**:
   ```julia
   Pkg.instantiate()
   ```

## Running the Project

To run the advection-diffusion simulation, execute the following command in the Julia REPL:

```julia
include("src/main.jl")
```

This will set up the problem, solve it, and output the results.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.