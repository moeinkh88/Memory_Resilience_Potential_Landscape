###############################################################
# 0.  Packages
###############################################################
using CSV, DataFrames
using DifferentialEquations
using Turing, Distributions
using MCMCChains, StatsPlots
using LinearAlgebra, StaticArrays
Turing.setprogress!(true)          # progress bars on/off

###############################################################
# 1.  Read & reshape the data  (unchanged)
###############################################################
data     = CSV.read("herbivore_dynamics_memory3.csv", DataFrame)
df_wide  = unstack(data, [:run, :t], :ic_id, :x)
rename!(df_wide, Symbol("1") => :x1, Symbol("2") => :x2)
df_wide  = sort(df_wide, [:run, :t])

unique_runs = collect(sort(unique(df_wide.run)))
n_groups    = length(unique_runs)
times       = collect(sort(unique(df_wide.t)))
n_times     = length(times)

obs = Array{Float64}(undef, n_groups, 2, n_times)   # (group, replicate, time)
for (i, grp) in enumerate(groupby(df_wide, :run))
    obs[i,1,:] = grp.x1
    obs[i,2,:] = grp.x2
end

###############################################################
# 2.  ODE with perceived biomass Y
###############################################################
"""
    herbivory_perceived!(du,u,p,t)

States per group:  u = [X₁, X₂, Y₁, Y₂]
Parameters:        p = (r,K,A,τ,b)
"""
function herbivory_perceived!(du,u,p,t)
    r, K, A, τ, b = p

    X1, X2, Y1, Y2 = u
    # safe denominators
    denom1 = max(A + Y1, 1e-6)
    denom2 = max(A + Y2, 1e-6)

    # dX/dt
    du[1] = r*X1*(1 - X1/K) - b*X1/denom1
    du[2] = r*X2*(1 - X2/K) - b*X2/denom2
    # dY/dt (first‑order low‑pass)
    du[3] = (X1 - Y1)/τ
    du[4] = (X2 - Y2)/τ
end

###############################################################
# 3.  One ODEProblem per group
###############################################################
tspan    = (minimum(times), maximum(times))
problems = Vector{ODEProblem}(undef, n_groups)

for i in 1:n_groups
    x10, x20 = obs[i,1,1], obs[i,2,1]     # initial true biomass
    u0       = [x10, x20, x10, x20]       # assume Y(0)=X(0)
    problems[i] = ODEProblem(herbivory_perceived!,
                             u0, tspan,
                             [1.0,1.0,1.0,1.0,1.0])  # dummy p
end

###############################################################
# 4.  Hierarchical Turing model
###############################################################
@model function fit_herbivory(obs, times, problems)
    # ---------- global priors ----------
    r   ~ LogNormal(log(0.8), 0.10)
    K   ~ LogNormal(log(6.0), 0.20)
    A   ~ LogNormal(log(0.1), 0.70)
    τ   ~ LogNormal(log(5.0), 0.40)      # time‑constant (days)

    # group‑specific grazing rates
    b   ~ filldist(LogNormal(log(1.0), 0.30), size(obs,1))

    # observation noise
    σ   ~ Truncated(Cauchy(0,1), 0, Inf)

    # ---------- likelihood ----------
    for i in 1:size(obs,1)          # loop over groups
        sol = solve(problems[i], Tsit5();
                    p = (r,K,A,τ,b[i]),
                    saveat = times,
                    abstol = 1e-8, reltol = 1e-6,
                    maxiters = 1e7)

        for j in 2:length(times)    # skip t₀ (fixed IC)
            pred      = sol(times[j])          # 4‑vector
            y_pred    = @SVector [pred[1], pred[2]]
            y_obs     = @SVector [obs[i,1,j], obs[i,2,j]]
            y_obs ~ MvNormal(y_pred, σ^2 * I)  # 2‑variate
        end
    end
end

model = fit_herbivory(obs, times, problems)

###############################################################
# 5.  Posterior sampling
###############################################################
nsamples, nchains = 2000, 4
chain = sample(model, NUTS(0.65), nsamples; nchains=nchains, progress=true)
println(chain)

###############################################################
# 6.  Posterior‑corner plot
###############################################################
corner(chain; size=(900,800))
savefig("posterior_corner.png")

###############################################################
# 7.  Posterior predictive checks
###############################################################
# collect parameters we need for prediction
param_names = [:r,:K,:A,:τ]
for i in 1:n_groups; push!(param_names, Symbol("b[$i]")); end

posterior_samples = sample(chain[param_names], 100; replace=false)
posterior_array   = Array(posterior_samples)

true_b = [first(data[data.run .== run, :].b) for run in unique_runs]

n_groups = length(unique_runs)
_, n_reps, _ = size(obs)

for (g, run_id) in enumerate(unique_runs)
    pplt = plot(title = "Group $run_id (true b = $(true_b[g]))",
                xlabel="Time", ylabel="Population",
                legend=:topright)

    for s in 1:size(posterior_array,1)
        r_s, K_s, A_s, τ_s = posterior_array[s, 1:4]
        b_s                = posterior_array[s, 4 + g]

        sol = solve(problems[g], Tsit5();
                    p=(r_s,K_s,A_s,τ_s,b_s), saveat=times)

        for rep in 1:n_reps
            y_pred = getindex.(sol.u, rep)   # rep=1→X₁, rep=2→X₂
            plot!(pplt, sol.t, y_pred; color=:gray, alpha=0.3, label=false)
        end
    end

    # overlay data
    for rep in 1:n_reps
        scatter!(pplt, times, obs[g,rep,:]; label="Rep $rep (obs)")
    end
    savefig(pplt, "posterior_predictive_group$(run_id).png")
end
