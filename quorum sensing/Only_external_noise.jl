using FdeSolver              # pkg> add FdeSolver
using Plots                  # pkg> add Plots
using Random                 # for reproducibility
using DataFrames, CSV        # pkg> add DataFrames CSV
using SpecialFunctions       # pkg> add SpecialFunctions

# --------------------------
# global model parameters
# --------------------------
V   = 3.0
K   = 1.0
A0p = 0.05

# decay d(rho)
decay(ρ) = 0.1 + (1 - ρ) / (ρ * (2 - ρ))      # singular at ρ -> 0

# time grid and measurement settings
t_end         = 200.0
tSpan         = (0.0, t_end)
sample_times  = collect(0.0:1.0:t_end)         # a few snapshots up to 15
h             = 0.01                          # solver step
nc_corr       = 3
replicates    = 1
noise_sigma   = .5                          # white noise std at measurement times

# fractional derivative orders

# single scenario with constant rho, no perturbation
rho_base1 = 0.32 
orders = 1.0
const scenario = (label="S", rho_base=rho_base1, Ainit=2.0, pert=nothing)

# rho schedule (constant)
make_rho_fun(rho_base; pert=nothing) = t -> rho_base

# right hand side
function QS_RHS(t, y, par)
    V, K, A0, ρ_fun, alpha, h, noise_sigma, rng = par
    A = y
    ρ = ρ_fun(t)
    det = V .* A.^2 ./ (K .+ A.^2) .+ A0 .- decay(ρ) .* A
    if noise_sigma > 0
        noise = noise_sigma .* gamma(alpha .+ 1) ./ (h .^ (alpha .- 0.5)) .* randn(rng)
        return det .+ noise
    else
        return det
    end
end

# simple linear interpolation for measurements at requested times
function interp_at(ts_in, ys_in, tq_in)
    ts = collect(ts_in)              # ensure Vector
    ys = vec(ys_in)                  # ensure Vector
    tq = collect(tq_in)              # ensure Vector
    n = length(ts)
    @assert n == length(ys)
    out = similar(tq)
    for (j, t) in enumerate(tq)
        if t ≤ ts[1]
            out[j] = ys[1]
        elseif t ≥ ts[end]
            out[j] = ys[end]
        else
            i = searchsortedlast(ts, t)
            i = max(1, min(i, n - 1))
            t1, t2 = ts[i], ts[i+1]
            y1, y2 = ys[i], ys[i+1]
            w = (t - t1) / (t2 - t1)
            out[j] = (1 - w) * y1 + w * y2
        end
    end
    return out
end


# run one experiment for given alpha and replicates
function run_experiment(exp_id::Int, α::Float64, sc; seed=42)
    rng = MersenneTwister(seed)
    ρ_fun = make_rho_fun(sc.rho_base; pert=sc.pert)

    # deterministic run for A_true
    par_det = (V, K, A0p, ρ_fun, α, h, 0.0, rng)
    t, Aapp = FDEsolver(QS_RHS, [tSpan[1], tSpan[2]], sc.Ainit, α, par_det; h=h, nc=nc_corr)

    # noiseless truth at sample times
    A_true = interp_at(t, Aapp, sample_times)

    # build dataframe rows for replicates with measurement noise
    rows = DataFrame(exp_id=Int[], scenario=String[], alpha=Float64[], replicate=Int[],
                     t=Float64[], rho=Float64[], A_true=Float64[], A_obs=Float64[],
                     rho_base=Float64[], pert_on=Float64[], pert_off=Float64[], rho_pert=Float64[])

    # collect rho value at each sample time for reference
    rho_at_t = [ρ_fun(tt) for tt in sample_times]

    for r in 1:replicates
        # new rng per replicate
        rng_r = MersenneTwister(seed + r)
        par_stoch = (V, K, A0p, ρ_fun, α, h, noise_sigma, rng_r)
        t_r, Aapp_r = FDEsolver(QS_RHS, [tSpan[1], tSpan[2]], sc.Ainit, α, par_stoch; h=h, nc=nc_corr)
        A_stoch = interp_at(t_r, Aapp_r, sample_times)
        A_obs = clamp.(A_stoch, 0.0, Inf)   # force nonnegative
        for (k, tt) in enumerate(sample_times)
            ρp = sc.pert === nothing ? NaN : sc.pert[3]
            t_on = sc.pert === nothing ? NaN : sc.pert[1]
            t_off = sc.pert === nothing ? NaN : sc.pert[2]
            push!(rows, (exp_id, sc.label, α, r,
                         tt, rho_at_t[k], A_true[k], A_obs[k],
                         sc.rho_base, t_on, t_off, ρp))
        end
    end
    return rows
end

# assemble experiments for each alpha
all_rows = DataFrame()
exp_counter = 0
for α in orders
    exp_counter += 1
    df = run_experiment(exp_counter, α, scenario; seed=1000 + 17*exp_counter)
    append!(all_rows, df)
end

using Plots
mkpath("plots")

# Compute potential landscape (independent of alpha)
rho = copy(rho_base1)
d = decay(rho)
A_range = collect(0:0.01:4)
U = - ( V .* (A_range .- sqrt(K) .* atan.(A_range ./ sqrt(K))) .+ A0p .* A_range .- (d / 2) .* A_range.^2 )

exp_ids = sort(unique(all_rows.exp_id))

for ex_id in exp_ids
    df1 = all_rows[all_rows.exp_id .== ex_id, :]

    scen  = df1.scenario[1]
    αval  = df1.alpha[1]
    reps  = sort(unique(df1.replicate))

    # Dynamics plot
    plt1 = plot(title = "Dynamics, alpha=$(round(αval,digits=2))",
               xlabel = "time", ylabel = "A", legend = :topright)

    # replicates
    for r in reps
        dfr = df1[df1.replicate .== r, :]
        plot!(plt1, dfr.t, dfr.A_obs, markersize=4, label=false)
    end

    # Potential plot
    plt2 = plot(A_range, U, title = "Potential Landscape", xlabel = "A", ylabel = "U(A)", lw=3, legend=false)

    # Combine into subplots
    plt = plot(plt1, plt2, layout = (1, 2), size=(800, 400))

    # save + show
    scen_safe = replace(string(scen), r"[^\w]+" => "_")
    α_safe    = replace(string(round(αval,digits=2)), "." => "-")
    
    display(plt)
end