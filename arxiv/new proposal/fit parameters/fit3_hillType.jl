# 1. Load herbivore dynamics data from CSV
using CSV, DataFrames
using DifferentialEquations
using Turing, Distributions
using MCMCChains, StatsPlots
using LinearAlgebra, StaticArrays

# Read the data file (assumes it's in the current directory or provide full path)
data = CSV.read("herbivore_dynamics_memory3.csv", DataFrame)

# 2. Reshape data: group each pair of trajectories with the same `b` value
# The data has columns: run (group ID), b (grazing rate per group), ic_id (trajectory ID), t (time), x (population)
# We pivot the data to have one row per time point per run, with separate columns for each trajectory (x1, x2)
df_wide = unstack(data, [:run, :t], :ic_id, :x)
rename!(df_wide, Symbol("1") => :x1, Symbol("2") => :x2)  # rename pivoted columns to x1, x2
df_wide = sort(df_wide, [:run, :t])  # ensure data is sorted by run and time

# Extract unique runs and time points
unique_runs = collect(sort(unique(df_wide.run)))
n_groups = length(unique_runs)
times = collect(sort(unique(df_wide.t)))    # observation times (assumed common to all groups)
n_times = length(times)

# Organize observed trajectories into a 3D array: obs[group, replicate, time_index]
obs = Array{Float64}(undef, n_groups, 2, n_times)
for (i, grp) in enumerate(groupby(df_wide, :run))
    # Each group has two trajectories: x1 and x2 at each time
    obs[i, 1, :] = grp.x1
    obs[i, 2, :] = grp.x2
end

# 3. Define the memoryless herbivory model ODE
# This is a logistic growth with grazing term model (no memory term).
# We combine the two trajectories for each group into a 2D system for efficiency.
###############################################################################
# MEMORY‑LESS HILL + θ‑LOGISTIC HERBIVORY MODEL
###############################################################################
# Return x^p if x>0, else 0, without branching that breaks AD
@inline pow_pos(x, p) = ifelse(x > 0, x^p, zero(x))

function herbivory_Hθ!(du, u, p, t)
    r, K, A, h, θ, b = p
    K_inv  = 1 / K
    Ah     = A^h           # pre‑compute, cheaper & AD‑friendly

    @inbounds for j in 1:2
        x = u[j]

        # SAFE powers — no DomainError when x≤0
        growth  = r * x               * (1 - pow_pos(x * K_inv, θ))
        grazing = b * pow_pos(x, h+1) / (pow_pos(x, h) + Ah)

        du[j] = growth - grazing
    end
end

# Prepare ODE problem for each group with its two initial conditions
# (We use a dummy parameter vector; actual parameters will be supplied in the solver.)
tspan     = (minimum(times), maximum(times))
par_guess = [1, 1, 1, 1, 1, 1]                 # 6 dummy params

problems  = Vector{ODEProblem}(undef, n_groups)
for i in 1:n_groups
    u0            = copy(obs[i, :, 1])         # initial conditions
    problems[i]   = ODEProblem(
                        herbivory_Hθ!, u0, tspan, par_guess
                    )
end


# 4. Build a hierarchical Turing model with global parameters r, K, A; group-specific grazing rates b[i]; and shared noise σ
@model function fit_Hθ(obs, times, problems)
    # GLOBAL
    r   ~ LogNormal(log(0.8), 0.15)
    K   ~ LogNormal(log(6.0), 0.25)
    A   ~ LogNormal(log(0.1), 0.7)

    h   ~ LogNormal(log(1.5), 0.4)      # Hill exponent  (≈1–3 typical)
    θ   ~ LogNormal(log(1.0), 0.35)     # θ‑logistic curvature (≈0.5–2)

    # GROUP‑SPECIFIC
    b   ~ filldist(LogNormal(log(1.0), 0.3), size(obs,1))

    # OBSERVATION NOISE
    σ   ~ Truncated(Cauchy(0,1), 0, Inf)

    # 5. Use DifferentialEquations.jl to solve the 2D ODE for each group, combining each pair of trajectories
    for i in 1:size(obs, 1)  # loop over groups
        # Solve the ODE for group i with current parameters [r, K, A, b[i]]
        # Using Tsit5() solver, at observation times for efficiency (saveat = times)
        # sol = solve(problems[i], Tsit5(); p=[r, K, A, b[i]], saveat=times)
        sol = solve(problems[i], Tsit5();
            p = (r, K, A, h, θ, b[i]),
            saveat = times,
            isoutofdomain = (u, p, t) -> any(u .< 0),  # shrink step instead of going negative
            abstol = 1e-8, reltol = 1e-6, maxiters = 1e7)

        for j in 2:length(times)
        # get the two‐element prediction at t = times[j]
        pred = sol(times[j])      # this is a Vector{<:Real}, length 2

        # collect your two observations into a vector
        y_obs = @SVector [obs[i,1,j], obs[i,2,j]]  # or simply [obs1[j], obs2[j]]

        # joint likelihood
        y_obs ~ MvNormal(pred, σ^2 * I)
        end
    end
end

# Instantiate the Turing model with our data
model = fit_Hθ(obs, times, problems)

# 6. Sample from the posterior using NUTS (No-U-Turn Sampler) with multiple chains
# (We run, for example, 4 chains of 2000 iterations each, with NUTS adaptation on each chain)
nsamples = 2000
nchains = 4
# Set Turing to use multiple threads for parallel chains if available (optional)
# Turing.setadbackend(:forwarddiff)  # default AD backend; adjust if needed
chain = sample(model, NUTS(0.65), nsamples; nchains=nchains, progress=true)

# Summarize posterior results (optional printout)
println(chain)

# 7. Save a corner plot of the posterior distributions for parameters
corner(chain; size=(900, 800))
savefig("posterior_corner.png")  # saved corner plot image of parameter posteriors

# 8. Posterior predictive plots: compare fitted trajectories vs observed data for each group
# We will draw several samples from the posterior and simulate the ODE to see model-predicted dynamics.
# Prepare to overlay many posterior sample trajectories ("posterior predictive") with the observed data.
# Extract parameter samples for r, K, A, and all b[i] from the chain:
param_names = [:r, :K, :A]  # base global parameters
for i in 1:n_groups
    push!(param_names, Symbol("b[$i]"))
end
# Take a subset of posterior samples for plotting (e.g., 100 random draws)
posterior_samples = sample(chain[param_names], 100; replace=false)
posterior_array = Array(posterior_samples)  # convert to matrix (rows = samples, columns = parameters)

# Retrieve actual b values per group from the data (for labeling plots)
true_b = [ first(data[data.run .== run, :].b) for run in unique_runs ]

# Loop over each group to generate its predictive plot

# assume
# posterior_array :: Matrix{Float64} with size (n_samples, 3 + n_groups)
#    columns 1:3 = r,K,A and columns 4:4+n_groups-1 = b_1 ... b_n_groups
# problems        :: Vector{ODEProblem} of length n_groups
# unique_runs     :: Vector of group IDs (1:n_groups)
# true_b          :: Vector of length n_groups (the “true” b for titles)
# obs             :: Array{Float64}(n_groups, n_reps, n_times)
# times           :: Vector{Float64} of length n_times

n_groups = length(unique_runs)
_, n_reps, _ = size(obs)

for (g, run_id) in enumerate(unique_runs)
    # new figure for group g
    p = plot(title = "Group $run_id (true b = $(true_b[g]))",
             xlabel = "Time", ylabel = "Population density",
             legend = :topright)

    # background: sample predictive trajectories
    for s in 1:size(posterior_array, 1)
        # unpack global params
        r_s, K_s, A_s = posterior_array[s, 1:3]
        # grab the b for this group
        b_s = posterior_array[s, 3 + g]

        # solve the ODEProblem for group g with these params
        sol = solve(problems[g],
                    Tsit5();
                    p = (r_s, K_s, A_s, b_s),
                    saveat = times)

        # plot each replicate's trajectory in gray
        for rep in 1:n_reps
            y_pred = getindex.(sol.u, rep)    # extract component `rep` from each solution state
            plot!(p, sol.t, y_pred;
                  color = :gray, alpha = 0.3, label = false)
        end
    end

    # overlay the actual data points
    for rep in 1:n_reps
        scatter!(p, times, obs[g, rep, :];
                 label = "Replicate $rep (obs)")
    end

    savefig(p, "posterior_predictive_group$(run_id).png")
end
