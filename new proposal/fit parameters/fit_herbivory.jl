# fit_herbivory.jl
using CSV, DataFrames
using DifferentialEquations
using Turing, Distributions, LinearAlgebra
using StatsPlots
using Random

# ------------------------------------------------------------
# 1 | Load and reshape the data
# ------------------------------------------------------------
df = CSV.read("herbivore_dynamics_memory2.csv", DataFrame)

# Treat each (run, ic_id) pair as a separate trajectory
gdf = groupby(df, [:run, :ic_id])

# Make vectors of vectors for easy indexing later
t_all  = [g.t               for g in gdf]         # Vector{Vector{Float64}}
x_all  = [g.x               for g in gdf]
x0_all = [first(g.ic_val)   for g in gdf]         # initial condition supplied
run_id  = [first(g.run)    for g in gdf]     # e.g. 1,1,2,2
nruns   = length(t_all)                      # 4 trajectories
ngroups = maximum(run_id)                    # 2 unique runs (= 2 b's)

# ------------------------- helper: prior guesses per run ----
b_guess = [mean(df.b[df.run .== k]) for k in 1:ngroups]   # [0.2, 0.8]

# ------------------------------------------------------------
# 2 | Define the ODE
# ------------------------------------------------------------
function herbivory!(dx, x, p, t)
    r, K, A, b = p               # p is a 4‑tuple
    dx[1] = r * x[1]*(1 - x[1]/K) - b * x[1] / (x[1] + A)
end

pvec=[.8, 6, .1, .2] # [r K A b[vectors but not repeated]]
prob = ODEProblem(herbivory!, [x0_all[1]], (minimum(t_all[1]), maximum(t_all[1])), pvec)
# ------------------------------------------------------------
# 3 | Turing model
# ------------------------------------------------------------
@model function fit_model(data, prob,x0_vec, b_known, run_idx, t_all)    
    # ---- global priors ----
    r ~ LogNormal(log(0.8), 0.1)          # r > 0
    K ~ LogNormal(log(6.0), 0.2)          # K > 0
    A ~ LogNormal(log(0.1), 0.7)         # small positive
    σ ~ Truncated(Cauchy(0, 1), 0.0, Inf)

    # ---- hierarchical: one grazing parameter per run ----
    nx0 = length(x0_vec)
    
    m = maximum(run_idx)                    # number of groups (runs)
    b = Vector{Float64}(undef, m)

    for i in 1:m
        # informative prior centred near the experimental value
        b[i] ~ LogNormal(log(b_known[i]+1e-6), 0.3)
    end

    # # ---- likelihood: solve ODE for every trajectory ----
    # for i in 1:nx0
    #     t_i  = t_all[i]
    #     bid   = b[run_idx[i]]                   # shared b for this trajectory
    #     x0   = x0_vec[i]
    #     pvec = [r, K, A, bid]

    #     prob = remake(prob; p = pvec, u0 = [x0])
    #     sol  = solve(prob, Tsit5(); saveat=t_i, abstol=1e-8, reltol=1e-6)

    #     x_pred = Array(sol)[:]              # length(t_i) = 201
    # end

    # for i in 1:length(x_pred)
    #     # likelihood: data ~ predicted
    #     data[:,i] ~ MvNormal(x_pred[i], σ^2 * I)  # I is identity matrix
    # end
    # ── forward model for every trajectory ────────────────────────
    nx0      = length(x0_vec)
    ntime    = length(t_all[1])            # 201 time points (all equal)

    sol0 = solve(remake(prob; u0=[x0_vec[1]], p=[r,K,A,b[run_idx[1]]]),
             Tsit5(); saveat=t_all[1], abstol=1e-8, reltol=1e-6)
    Tdual = eltype(sol0.u[1])                   # usually ForwardDiff.Dual{…}
    preds_mt = Matrix{Tdual}(undef, nx0, ntime)
    for j in 1:nx0
        pvec = [r, K, A, b[run_idx[j]]]
        probj = remake(prob; p = pvec, u0 = [x0_vec[j]])
        sol   = solve(probj, Tsit5(); saveat = t_all[j],
                      abstol = 1e-8, reltol = 1e-6)
        XX=hcat(sol.u...)
        preds_mt[j, :] = vec(XX)     # store the whole row
    end

    # ── reshape to “vector‑of‑vectors” (one vector per time point) ─
    x_pred = [Vector(preds_mt[:, t]) for t in 1:ntime]  # length = 201

    # ── likelihood: data and predictions now share the same shape ─
        data_mt = hcat(data...)'                # size (nx0, ntime)
    for t in 1:ntime
        data_mt[:, t] ~ MvNormal(x_pred[t], σ^2 * I)   # I chooses dim automatically
    end
    return nothing
end

# ------------------------------------------------------------
# 4 | Run NUTS
# ------------------------------------------------------------
Random.seed!(2025)
model  = fit_model(x_all, prob, x0_all, b_guess, run_id, t_all)
nchain = 2
samples = sample(model, NUTS(0.65), MCMCSerial(), 500, nchain; progress=true)

# ------------------------------------------------------------
# 5 | Summaries & quick diagnostics
# ------------------------------------------------------------
describe(samples)

# corner plot of the top‑level parameters
p = cornerplot(samples; vars=[:r, :K, :A, :σ])
savefig(p, "posterior_corner.png");  println("Saved ↪ posterior_corner.png")

# ------------------------------------------------------------
# 6 | Posterior predictive check for the first trajectory
# ------------------------------------------------------------
posterior = Matrix(samples)
using Statistics
r̂ = mean(posterior[:, sampler_varindex(samples, :r)])
K̂ = mean(posterior[:, sampler_varindex(samples, :K)])
Â = mean(posterior[:, sampler_varindex(samples, :A)])
b̂ = mean(posterior[:, sampler_varindex(samples, Symbol("b[1]"))])

prob̂ = ODEProblem(herbivory!, [x0_all[1]], (minimum(t_all[1]), maximum(t_all[1])),
                   (r̂, K̂, Â, b̂))
sol̂  = solve(prob̂, Tsit5(); saveat=t_all[1])

plot(t_all[1], x_all[1], seriestype=:scatter, label="data", markersize=3)
plot!(t_all[1], sol̂[:], label="posterior mean", lw=2)
savefig("traj_fit_run1.png"); println("Saved ↪ traj_fit_run1.png")
