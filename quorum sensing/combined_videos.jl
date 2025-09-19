using FdeSolver              # pkg> add FdeSolver
using Plots                  # pkg> add Plots
using Random                 # for reproducibility
using DataFrames, CSV        # pkg> add DataFrames CSV
using SpecialFunctions       # pkg> add SpecialFunctions
using Polynomials            # pkg> add Polynomials

# --------------------------
# global model parameters
# --------------------------
V   = 3.0
K   = 1.0
A0p = 0.05

# decay d(rho)
decay(ρ) = 0.1 + (1 - ρ) / (ρ * (2 - ρ))      # singular at ρ -> 0

# time grid and measurement settings
t_end         = 650.0
tSpan         = (0.0, t_end)
sample_times  = collect(0.0:1.0:t_end)         # a few snapshots up to 15
h             = 0.01                          # solver step
nc_corr       = 3
replicates    = 1
noise_sigma   = 1                          # white noise std at measurement times

# fractional derivative orders
orders = [0.8, 1.0]

# single scenario with constant rho, no perturbation
const scenario = (label="S", rho_base=0.35, Ainit=2.0, pert=nothing)

# rho schedule (stepwise changing)
function rho_fun(t)
    if t < 40
        0.4
    elseif t < 180
        0.35
    elseif t < 280
        0.3
    elseif t < 450
        0.3
    elseif t < 550
        0.3
    else
        0.35
    end
end

make_rho_fun(rho_base; pert=nothing) = rho_fun

# right hand side
function QS_RHS(t, y, par)
    V, K, A0, ρ_fun, alpha, h, noise_sigma, rng = par
    A = y
    ρ = ρ_fun(t)
    det = V .* A.^2 ./ (K .+ A.^2) .+ A0 .- decay(ρ) .* A
    if noise_sigma > 0
        # noise = noise_sigma .* gamma(alpha .+ 1) ./ (h .^ (alpha .- 0.5)) .* randn(rng)
                noise = noise_sigma .* randn(rng)
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
gr()
mkpath("plots")

# Pre-compute bifurcation diagram
rhos = collect(range(0.05, 0.6, length=400))  # avoid rho near zero singularity
ρ_stable   = Float64[]
A_stable   = Float64[]
ρ_unstable = Float64[]
A_unstable = Float64[]
tol = 1e-9

for ρ in rhos
    d = decay(ρ)
    c0, c1, c2, c3 = (A0p*K, -d*K, (V + A0p), -d)  # cubic_coeffs(d)
    p = Polynomial([c0, c1, c2, c3])
    rts = roots(p)
    for z in rts
        if abs(imag(z)) < tol
            A2 = real(z)
            if A2 >= 0
                fp = 2*V*K*A2 / (K + A2^2)^2 - d
                if fp < 0
                    push!(ρ_stable, ρ);   push!(A_stable, A2)
                else
                    push!(ρ_unstable, ρ); push!(A_unstable, A2)
                end
            end
        end
    end
end

A_range = collect(0:0.01:4)

# Extract data for both alphas
df08 = all_rows[all_rows.alpha .== 0.8, :]
dfr08 = df08[df08.replicate .== 1, :]
t_vals08 = dfr08.t
A_obs08 = dfr08.A_obs
rho_at_t08 = dfr08.rho

df10 = all_rows[all_rows.alpha .== 1.0, :]
dfr10 = df10[df10.replicate .== 1, :]
t_vals10 = dfr10.t
A_obs10 = dfr10.A_obs
rho_at_t10 = dfr10.rho

# Assume t_vals are the same
t_vals = t_vals08

anim = @animate for i = 1:15:length(t_vals)  # step to speed up animation
    # Current rho (same for both)
    current_rho = rho_at_t08[i]
    d = decay(current_rho)

    # Compute current potential (same for both since independent of alpha)
    U = - ( V .* (A_range .- sqrt(K) .* atan.(A_range ./ sqrt(K))) .+ A0p .* A_range .- (d / 2) .* A_range.^2 )
    U = U .- minimum(U) # shift to min 0

    # Left: alpha=0.8
    # Bifurcation left
    p1_left = plot(rhos, ones(length(rhos)), label="", alpha=0, xlabel="ρ", ylabel="Equilibria A",
                   title="Bifurcation (alpha=0.8)", ylim=(0, 4), legend=false)
    scatter!(p1_left, ρ_stable, A_stable, label="", markersize=1, color=:green)
    plot!(p1_left, ρ_unstable, A_unstable, label="", lw=2, color=:red, linestyle=:dash)
    vline!(p1_left, [current_rho], label="", color=:gray, lw=2)
    scatter!(p1_left, [current_rho], [A_obs08[i]], markersize=7, color=:black, label="")

    # Potential left
    p2_left = plot(A_range, U, title="Potential (rho=$(round(current_rho, digits=3)))", xlabel="A", ylabel="U(A)", lw=2, color=:black, legend=false,
                   ylims=(-0.05, maximum(U)+0.05))
    plot!(p2_left, A_range, U, fillrange=minimum(U)-0.05, fillalpha=0.3, c=:purple)
    scatter!(p2_left, [A_obs08[i]], [U[argmin(abs.(A_range .- A_obs08[i]))]], markersize=10, color=:black)

    # Dynamics left
    p3_left = plot(A_obs08[1:i], t_vals[1:i], title="Dynamics", xlabel="A", ylabel="Time", lw=2, color=:blue, legend=false,
                   yflip=true, xlims=(0, 4), ylims=(0, t_end))

    # Right: alpha=1.0
    # Bifurcation right
    p1_right = plot(rhos, ones(length(rhos)), label="", alpha=0, xlabel="ρ", ylabel="Equilibria A",
                    title="Bifurcation (alpha=1.0)", ylim=(0, 4), legend=false)
    scatter!(p1_right, ρ_stable, A_stable, label="", markersize=1, color=:green)
    plot!(p1_right, ρ_unstable, A_unstable, label="", lw=2, color=:red, linestyle=:dash)
    vline!(p1_right, [current_rho], label="", color=:gray, lw=2)
    scatter!(p1_right, [current_rho], [A_obs10[i]], markersize=7, color=:black, label="")

    # Potential right
    p2_right = plot(A_range, U, title="Potential (rho=$(round(current_rho, digits=3)))", xlabel="A", ylabel="U(A)", lw=2, color=:black, legend=false,
                    ylims=(-0.05, maximum(U)+0.05))
    plot!(p2_right, A_range, U, fillrange=minimum(U)-0.05, fillalpha=0.3, c=:purple)
    scatter!(p2_right, [A_obs10[i]], [U[argmin(abs.(A_range .- A_obs10[i]))]], markersize=10, color=:black)

    # Dynamics right
    p3_right = plot(A_obs10[1:i], t_vals[1:i], title="Dynamics", xlabel="A", ylabel="Time", lw=2, color=:blue, legend=false,
                    yflip=true, xlims=(0, 4), ylims=(0, t_end))

    # Combine: left and right columns
    plot(p1_left, p1_right, p2_left, p2_right, p3_left, p3_right, layout=(3,2), size=(1600, 900))
end

# Save as GIF
fn = "plots/combined_dynamics.gif"
gif(anim, fn, fps=15)
@info "saved $fn"