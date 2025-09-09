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
t_end         = 20.0
tSpan         = (0.0, t_end)
sample_times  = collect(0.0:1.0:t_end)         # a few snapshots up to 15
h             = 0.01                          # solver step
nc_corr       = 3
replicates    = 50
noise_sigma   = 0.05                          # white noise std at measurement times

# fractional derivative orders
orders = [0.7, 0.8, 0.9, 0.97, 1.0]

# scenarios
# 1 to 4: rho=0.35, A(0)=2.0, pulse rho=0.25 from t in [3,5]
# 5 to 8: rho=0.40, A(0)=0.05, pulse rho=0.55 from t in [3,5]
# 9 to 12: rho=0.45, A(0)=2.0, no pulse
const scenarios = [
    (label="S1", rho_base=0.35, Ainit=2.0,  pert=(3.0, 5.0, 0.2)),
    (label="S2", rho_base=0.40, Ainit=0.05, pert=(3.0, 6.0, 0.85)),
    (label="S3", rho_base=0.45, Ainit=2.0,  pert=nothing)
]

# rho schedule
make_rho_fun(rho_base; pert=nothing) = t -> begin
    if pert === nothing
        return rho_base
    else
        t_on, t_off, rho_pert = pert
        return (t ≥ t_on && t ≤ t_off) ? rho_pert : rho_base
    end
end

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


# run one experiment for given alpha, scenario, and replicates
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

# assemble all 12 experiments
all_rows = DataFrame()
exp_counter = 0
for (gi, sc) in enumerate(scenarios)
    for α in orders
        exp_counter += 1
        df = run_experiment(exp_counter, α, sc; seed=1000 + 17*exp_counter)
        append!(all_rows, df)
    end
end

# save CSV
CSV.write("qs_bistable_dataset1.csv", all_rows)
println("Saved qs_bistable_dataset1.csv with ", nrow(all_rows), " rows.")

using Plots
mkpath("plots")

exp_ids = sort(unique(all_rows.exp_id))

for ex_id in exp_ids
    df1 = all_rows[all_rows.exp_id .== ex_id, :]

    scen  = df1.scenario[1]
    αval  = df1.alpha[1]
    reps  = sort(unique(df1.replicate))

    plt = plot(title = "Experiment $(ex_id)  scenario=$(scen)  alpha=$(round(αval,digits=2))",
               xlabel = "time", ylabel = "A", legend = :topright)

    # true (noiseless) series once
    true_series = combine(groupby(df1, :t), :A_true => first => :A_true)
    # plot!(plt, true_series.t, true_series.A_true, lw=3, label="true")
    plot!(plt, true_series.t, true_series.A_true, lw=3, label=false)

    # replicates
    for r in reps
        dfr = df1[df1.replicate .== r, :]
        # scatter!(plt, dfr.t, dfr.A_obs, markersize=4, label="rep $r")
        # scatter!(plt, dfr.t, dfr.A_obs, markersize=4, label=false)
        plot!(plt, dfr.t, dfr.A_obs, markersize=4, label=false)
    end

    # mark perturbation window if present
    has_pert = any(.!isnan.(df1.pert_on)) && any(.!isnan.(df1.pert_off))
    if has_pert
        t_on  = first(df1.pert_on[.!isnan.(df1.pert_on)])
        t_off = first(df1.pert_off[.!isnan.(df1.pert_off)])
        vline!(plt, [t_on, t_off], lc=:gray, ls=:dash, label=false)
    end

    # save + show
    scen_safe = replace(string(scen), r"[^\w]+" => "_")
    α_safe    = replace(string(round(αval,digits=2)), "." => "-")
    fn = "plots/exp_$(lpad(ex_id,2,'0'))_$(scen_safe)_alpha$(α_safe).png"
    savefig(plt, fn)
    display(plt)
    @info "saved $fn"
end