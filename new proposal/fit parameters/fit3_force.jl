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
function herbivory_model!(du, u, p, t)
    # Parameters: p = (r, K, A, b) where 
    # r = intrinsic growth rate, K = carrying capacity, A = half-saturation constant for grazing, b = grazing rate for this group.
    @inbounds begin
        let r = p[1], K = p[2], A = p[3], b = p[4]
            # u[1] and u[2] are the two replicate population densities for this group.
            du[1] = r * u[1] * (1 - u[1] / K) - b * u[1] / (A + u[1])
            du[2] = r * u[2] * (1 - u[2] / K) - b * u[2] / (A + u[2])
        end
    end
end

# Prepare ODE problem for each group with its two initial conditions
# (We use a dummy parameter vector; actual parameters will be supplied in the solver.)
tspan = (minimum(times), maximum(times))
problems = Vector{ODEProblem}(undef, n_groups)
for i in 1:n_groups
    # Initial conditions for the two trajectories in group i (at t=0)
    # We take the observed initial values (assuming no observation noise at t=0)
    u0 = copy(obs[i, :, 1])  # initial state vector of length 2 (first time point of each replicate)
    # Create an ODEProblem for this group (with dummy parameters, to be overridden in the solver)
    problems[i] = ODEProblem(herbivory_model!, u0, tspan, [1.0, 1.0, 1.0, 1.0])
end

# 4. Build a hierarchical Turing model with global parameters r, K, A; group-specific grazing rates b[i]; and shared noise σ
@model function fit_herbivory(obs, times, problems)
    # Priors for global parameters (constrained to positive values using log-normal or similar distributions)
    ### global parameters
    r ~ truncated(Normal(0.8, 0.05),      0,  1.2)   # 95 % of mass in [0.7, 0.9]
    K ~ truncated(Normal(6.0, 0.3),       0,  10)    # keeps K quite close to 6
    A ~ truncated(Normal(0.1, 0.03),      0,  0.25)  # almost all mass in [0.04, 0.16]

    ### group‑specific grazing rates
    b ~ filldist(truncated(Normal(0.8, 0.2), 0, 1.2), size(obs,1))

    ### shared observation noise
    σ ~ truncated(Cauchy(0,1), 0, Inf)

    # 5. Use DifferentialEquations.jl to solve the 2D ODE for each group, combining each pair of trajectories
    for i in 1:size(obs, 1)  # loop over groups
        # Solve the ODE for group i with current parameters [r, K, A, b[i]]
        # Using Tsit5() solver, at observation times for efficiency (saveat = times)
        # sol = solve(problems[i], Tsit5(); p=[r, K, A, b[i]], saveat=times)
        sol = solve(problems[i], Tsit5(); p=[r,K,A,b[i]],
            saveat=times, abstol=1e-8, reltol=1e-6,
                maxiters=1e7        # raise the default cap (from ~1e5) to 10 million
                )
        # sol = solve(
        # problems[i], Rodas5();
        # p      = [r, K, A, b[i]],
        # saveat = times,
        # abstol  = 1e-8,
        # reltol  = 1e-6)
        # The solution `sol` yields a vector of states for each time point (sol[j] is a 2-element vector at time times[j]).
        # 6. Likelihood: compare predicted vs observed values with noise σ
        # We exclude t=0 from the likelihood since initial conditions are known and fixed.
        # for j in 2:length(times)   # start from index 2 (first index is t=0 initial condition)
        #     pred = sol(times[j])        # interpolate at the exact time
        #     obs[i, 1, j] ~ Normal(pred[1], σ)   # observed trajectory 1 at time j
        #     obs[i, 2, j] ~ Normal(pred[2], σ)   # observed trajectory 2 at time j
        # end
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
model = fit_herbivory(obs, times, problems)

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
