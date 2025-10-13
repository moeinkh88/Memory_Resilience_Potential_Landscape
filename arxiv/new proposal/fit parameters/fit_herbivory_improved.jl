# fit_herbivory_improved.jl

using CSV, DataFrames
using DifferentialEquations, StaticArrays
using Turing, Distributions, LinearAlgebra
using Random
using StatsPlots   # for cornerplot
using Plots        # for posterior‑predictive plot

# ----------------------------------------
# 1 | Load & Organize Data
# ----------------------------------------
df = CSV.read("herbivore_dynamics_memory2.csv", DataFrame)

# group each (run, ic_id) pair into one trajectory
gdf     = groupby(df, [:run, :ic_id])
t_all   = [Float64.(g.t)          for g in gdf]   # ensure times are Float64
x_all   = [g.x[:]                 for g in gdf]
x0_all  = [first(g.ic_val)        for g in gdf]
run_idx = [first(g.run)           for g in gdf]
ngroup  = length(unique(run_idx))

# center priors for b on the empirical means per run
b_prior = [mean(df.b[df.run .== r]) for r in 1:ngroup]


# ----------------------------------------
# 2 | Define ODE (out‑of‑place, StaticArrays)
# ----------------------------------------
function herbivory(u, p, t)
    r, K, A, b = p
    return @SVector [ r * u[1] * (1 - u[1]/K) - b * u[1] / (u[1] + A) ]
end

# dummy problem for reuse via `remake`
tspan_base = (minimum(vcat(t_all...)), maximum(vcat(t_all...)))
u0_base    = SVector(x0_all[1])               # placeholder u0
p0_base    = (0.8, 6.0, 0.1, b_prior[1])      # placeholder p
base_prob  = ODEProblem(herbivory, u0_base, tspan_base, p0_base)


# ----------------------------------------
# 3 | Helper: solve one trajectory
# ----------------------------------------
function predict_trajectory(r, K, A, b, x0, times)
    t_f = Float64.(times)
    prob = remake(base_prob;
                  p      = (r, K, A, b),
                  u0     = SVector(x0),
                  tspan  = (minimum(t_f), maximum(t_f)))
    sol = solve(prob, Tsit5();
                saveat = t_f,
                abstol = 1e-6,
                reltol = 1e-6)

    # check return code against symbol :Success
    if sol.retcode != :Success
        error("ODE solver failed with retcode $(sol.retcode)")
    end

    return Array(sol)[:]   # extract the 1D solution as Vector{Float64}
end


# ----------------------------------------
# 4 | Turing Model
# ----------------------------------------
@model function fit_model(x_obs, x0_vec, tvecs, run_idx, b_prior)
    # global priors
    r ~ LogNormal(log(0.8), 0.1)
    K ~ LogNormal(log(6.0), 0.2)
    A ~ LogNormal(log(0.1), 0.7)
    σ ~ truncated(Cauchy(0, 1), 0, Inf)

    # hierarchical grazing rates, one b[j] per run
    m = maximum(run_idx)
    b = Vector{Float64}(undef, m)
    for j in 1:m
        b[j] ~ LogNormal(log(b_prior[j] + 1e-6), 0.3)
    end

    # likelihood: each trajectory is one multivariate Normal
    for i in 1:length(x_obs)
        x_pred = predict_trajectory(
            r, K, A, b[run_idx[i]],
            x0_vec[i], tvecs[i]
        )
        x_obs[i] ~ MvNormal(x_pred, σ^2 * I)
    end
end


# ----------------------------------------
# 5 | Run NUTS Sampling
# ----------------------------------------
Random.seed!(2025)
model  = fit_model(x_all, x0_all, t_all, run_idx, b_prior)
chains = sample(model, NUTS(0.65), MCMCSerial(), 500, 2; progress=true)


# ----------------------------------------
# 6 | Diagnostics & Summaries
# ----------------------------------------
describe(chains)

# corner plot for top‑level parameters
cornerplot(chains; vars=[:r, :K, :A, :σ]) |> savefig("posterior_corner.png")
println("Saved → posterior_corner.png")

# posterior‑predictive check for the first trajectory
posterior = Array(chains)

# extract posterior means
r̂ = mean(posterior[:, chainvarinfo(chains).ind[:r]])
K̂ = mean(posterior[:, chainvarinfo(chains).ind[:K]])
Â = mean(posterior[:, chainvarinfo(chains).ind[:A]])
b̂ = mean(posterior[:, chainvarinfo(chains).ind[Symbol("b[1]")]])

sol̂ = predict_trajectory(r̂, K̂, Â, b̂, x0_all[1], t_all[1])

scatter(t_all[1], x_all[1], label="data", ms=3)
plot!(t_all[1], sol̂, lw=2, label="posterior mean")
savefig("traj_fit_run1.png")
println("Saved → traj_fit_run1.png")
