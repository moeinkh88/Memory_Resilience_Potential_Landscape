########################################################################
# herbivory_fit_optim.jl
#
# Fit the memory‑less Hill + θ‑logistic herbivory model
# by non‑linear least squares using Optim.jl.
#
# GLOBAL parameters:         r, K, A, h, θ
# GROUP‑specific parameters: b[1 : n_groups]
#
# All parameters are optimised in log‑space to enforce positivity.
#
# ── Required packages ────────────────────────────────────────────────
#   import Pkg
#   Pkg.add([
#       "CSV", "DataFrames", "DifferentialEquations",
#       "Optim", "ForwardDiff", "StaticArrays", "LinearAlgebra"
#   ])
########################################################################

using CSV, DataFrames
using DifferentialEquations
using Optim
using ForwardDiff
using StaticArrays, LinearAlgebra

# -----------------------------------------------------------
# 1. READ & RESHAPE DATA
# -----------------------------------------------------------
data = CSV.read("herbivore_dynamics_memory3.csv", DataFrame)

df_wide = unstack(data, [:run, :t], :ic_id, :x)
rename!(df_wide, Symbol("1") => :x1, Symbol("2") => :x2)
sort!(df_wide, [:run, :t])

runs      = sort(unique(df_wide.run))
n_groups  = length(runs)
times     = sort(unique(df_wide.t))
n_times   = length(times)

# 3‑D array: obs[group, replicate, time_index]
obs = Array{Float64}(undef, n_groups, 2, n_times)
for (i, gdf) in enumerate(groupby(df_wide, :run))
    obs[i,1,:] = gdf.x1
    obs[i,2,:] = gdf.x2
end

init_conds = [obs[i,:,1] for i in 1:n_groups]      # initial state for each group

# -----------------------------------------------------------
# 2. MODEL DEFINITION
# -----------------------------------------------------------
@inline pow_pos(x,p) = ifelse(x>0, x^p, zero(x))

function herbivory_Hθ!(du,u,p,t)
    r, K, A, h, θ, b = p
    Kinv = 1/K
    Ah   = A^h
    @inbounds for j in 1:2
        x        = u[j]
        growth   = r * x * (1 - pow_pos(x*Kinv, θ))
        grazing  = b * pow_pos(x, h+1) / (pow_pos(x,h) + Ah)
        du[j]    = growth - grazing
    end
end

function simulate_group(global_pars::NTuple{5,Float64},
                        b::Float64,
                        u0::AbstractVector,
                        ts::AbstractVector;
                        abstol=1e-8, reltol=1e-6)

    r,K,A,h,θ = global_pars
    prob = ODEProblem(herbivory_Hθ!, u0, (ts[1], ts[end]),
                      (r,K,A,h,θ,b))
    sol  = solve(prob, Tsit5();
                 saveat        = ts,
                 abstol        = abstol,
                 reltol        = reltol,
                 isoutofdomain = (u,p,t)->any(u.<0))
    return hcat(sol.u...)               # 2 × n_times matrix
end

# -----------------------------------------------------------
# 3. LEAST‑SQUARES OBJECTIVE (in log‑space)
# -----------------------------------------------------------
function build_objective(obs, times, init_conds)
    n_groups, _, _ = size(obs)

    function obj(p_log::AbstractVector)          # log‑parameters
        pars = exp.(p_log)                       # back‑transform
        gp   = pars[1:5]                         # global: r,K,A,h,θ
        bs   = view(pars, 6:length(pars))        # group‑specific b’s

        err  = 0.0
        @inbounds for i in 1:n_groups
            pred = simulate_group(tuple(gp...), bs[i], init_conds[i], times)
            # skip first point (model matches IC exactly)
            err += sum(abs2, pred[:,2:end] .- obs[i,:,2:end])
        end
        return err
    end
    return obj
end

obj = build_objective(obs, times, init_conds)

# -----------------------------------------------------------
# 4. OPTIMISATION SET‑UP
# -----------------------------------------------------------
p0 = vcat(log.([0.8, 6.0, 0.1, 1.5, 1.0]),   # r,K,A,h,θ  (guesses)
          fill(log(1.0), n_groups))          # b[1:n_groups]

lower = fill(-10.0, length(p0))              # exp(−10) ≈ 4.5 × 10⁻⁵
upper = fill( 10.0, length(p0))              # exp(10)  ≈ 2.2 × 10⁴

opt_state = optimize(
    obj,
    lower, upper,
    p0,
    Fminbox(LBFGS(
        m = 10,     # memory size for L-BFGS
        ftol = 1e-6, # function tolerance
        gtol = 1e-6, # gradient tolerance
    )),
    Optim.Options(iterations      = 10_000,
                  show_trace      = true,
                  show_every      = 50,
                  allow_f_increases = false);
    autodiff = :forward
)

println("\n-------- Optimisation summary --------")
println(opt_state)

# -----------------------------------------------------------
# 5. RESULTS
# -----------------------------------------------------------
p_hat = exp.(Optim.minimizer(opt_state))
r̂, K̂, Â, ĥ, θ̂ = p_hat[1:5]
b̂               = p_hat[6:end]

println("\nFitted global parameters:")
@printf("   r   = %.5g\n   K   = %.5g\n   A   = %.5g\n   h   = %.5g\n   θ   = %.5g\n",
        r̂, K̂, Â, ĥ, θ̂)

println("\nFitted grazing coefficients (b) per group:")
for (i, bval) in enumerate(b̂)
    @printf("   Group %s : %.5g\n", runs[i], bval)
end

# -----------------------------------------------------------
# 6. (Optional) DIAGNOSTIC PLOTS  ---------------------------
# Uncomment if Plots.jl is available
# using Plots
# for i in 1:n_groups
#     pred = simulate_group((r̂, K̂, Â, ĥ, θ̂), b̂[i], init_conds[i], times)
#     p = plot(title = "Group $(runs[i])",
#              xlabel="Time", ylabel="Population density")
#     scatter!(times, vec(obs[i,1,:]); label="Obs rep 1")
#     scatter!(times, vec(obs[i,2,:]); label="Obs rep 2")
#     plot!(times, vec(pred[1,:]);       label="Fit rep 1")
#     plot!(times, vec(pred[2,:]);       label="Fit rep 2")
#     display(p)
# end
########################################################################
