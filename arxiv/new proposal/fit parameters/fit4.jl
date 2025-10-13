# fit_herbivory.jl

############################################################
# 0) (First time only) — Activate project & ensure packages
############################################################
# import Pkg
# Pkg.activate(@__DIR__)    # activate a project in this folder
# Pkg.instantiate()         # install packages from Manifest.toml if present
# Pkg.add([
#   "CSV", "DataFrames",
#   "DifferentialEquations", "SciMLSensitivity",
#   "Turing", "Distributions", "MCMCChains"
# ])

############################################################
# 1) Load all needed libraries
############################################################
using CSV, DataFrames
using DifferentialEquations, SciMLSensitivity
using Turing, Distributions
using MCMCChains

############################################################
# 2) Read in your data CSV
############################################################
df = CSV.read("herbivore_dynamics_memory2.csv", DataFrame)

# identify the distinct runs (groups) and initial‐condition IDs
runs = sort(unique(df.run))      # e.g. [1, 2]
ics  = sort(unique(df.ic_id))    # e.g. [1, 2]

############################################################
# 3) Define the integer‐order herbivory ODE
############################################################
function herbivory!(du, u, p, t)
    r, K, A, b = p
    x = u[1]
    du[1] = r * x * (1 - x/K) - b * x / (x + A)
end

############################################################
# 4) Build the Turing model for joint inference
############################################################
Turing.setadbackend(:forwarddiff)   # use ForwardDiff for ODE AD

@model function fit_herbivory(df)
    # ===== Priors =====
    r  ~ truncated(Normal(1.0, 0.5), 0, Inf)    # intrinsic growth rate
    K  ~ truncated(Normal(50.0, 20.0), 0, Inf)  # carrying capacity
    A  ~ truncated(Normal(10.0, 5.0),  0, Inf)  # half‐saturation constant
    b1 ~ truncated(Normal(1.0, 0.5),  0, Inf)   # herbivory rate for run=1
    b2 ~ truncated(Normal(1.0, 0.5),  0, Inf)   # herbivory rate for run=2
    σ  ~ InverseGamma(2, 3)                    # observation noise

    # ===== Loop over each (run, ic) trajectory =====
    for run in runs
      for ic in ics
        idx   = findall((df.run .== run) .& (df.ic_id .== ic))
        if isempty(idx)  # skip if no data for this pair
          continue
        end

        times = df.t[idx]
        obs   = df.x[idx]
        u0    = obs[1]                   # treat first obs as known IC
        bval  = run == 1 ? b1 : b2       # select b1 or b2 by run

        # solve ODE at the observation times
        prob = ODEProblem(herbivory!, [u0], (times[1], times[end]),
                          (r, K, A, bval))
        # sol  = solve(prob, Tsit5(); saveat=times)
        sol  = solve(
        prob,
        Tsit5();
            saveat=times,
        reltol=1e-6,
        abstol=1e-6,
        maxiters=1_000_000,
        # sensealg=ForwardDiffSensitivity(),
        )
        # likelihood: skip j=1 (no noise on IC)
        for j in 2:length(times)
          obs[j] ~ Normal(sol[j][1], σ)
        end
      end
    end
end

############################################################
# 5) Run the sampler
############################################################
model = fit_herbivory(df)
# Use NUTS; 4 chains × 2000 samples
chains = sample(model, NUTS(0.65), 2000; chains=4)

############################################################
# 6) Print a summary of the posterior samples
############################################################
println(describe(chains))
