using FdeSolver, Random, Statistics, Plots
using Printf

# Fix the random seed for reproducibility
Random.seed!(123)

# Model parameters
α_mem   = 0.7    # fractional order (memory case, 0<α<1)
α_nomem = 1.0    # standard case (no memory, α=1)
y0      = 1.0    # initial state (start near one stable equilibrium)

# Simulation horizon and time step
T       = 1000.0             # total simulation time
h       = 0.01               # time step for integration (solver step size)
tSpan   = [0.0, T]           # time interval

# Noise parameters for additive and multiplicative cases
σ_add       = 0.1            # noise intensity for additive noise (small, within-basin variability)
σ_add_large = 0.3            # larger noise for triggering switches between basins
σ_mult      = 0.2            # base noise intensity for multiplicative noise

# Pre-generate piecewise-constant noise values for each time step interval
Nsteps          = Int(floor(T/h))               # number of intervals
noise_vals_add  = [ σ_add/sqrt(h) * randn() for k in 1:Nsteps ];       # small additive noise
noise_vals_add_large = [ σ_add_large/sqrt(h) * randn() for k in 1:Nsteps ];  # larger additive noise

# Define noise functions that return the pre-sampled noise value at time t
noise_func_add  = t -> begin 
    # Determine index of interval for time t
    let idx = min(floor(Int, t/h) + 1, Nsteps)
        return noise_vals_add[idx]
    end
end

noise_func_add_large = t -> begin 
    let idx = min(floor(Int, t/h) + 1, Nsteps)
        return noise_vals_add_large[idx]
    end
end

using Random

"Return a callable noise function ξ(t) for [t0,t1] with step h and sd σ (white-noise approx)."
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=123)
    N = Int(ceil((t1 - t0)/h))                # cover full span
    rng = MersenneTwister(seed)
    vals = (σ/√h) .* randn(rng, N)            # piecewise-constant increments
    t0_, h_ = float(t0), float(h)
    return (t::Real)-> begin
        idx = clamp(Int(floor((t - t0_)/h_)) + 1, 1, N)
        return vals[idx]
    end
end

# Define a state-dependent noise amplitude function for multiplicative noise.
# Here we choose noise amplitude that is **higher near the unstable state x≈0** and lower near stable states x≈±1.
noise_amp = x -> begin 
    # Example: amplitude = 1 - 0.8*x^2 (clamped to non-negative)
    local amp = 1.0 - 0.8*x^2
    return amp < 0 ? 0.0 : amp
end

# Pre-generate noise for multiplicative case as well (to have a deterministic sample path for fair comparison)
noise_vals_mult = [ σ_mult/sqrt(h) * randn() for k in 1:Nsteps ]
noise_func_mult = t -> begin 
    let idx = min(floor(Int, t/h) + 1, Nsteps)
        return noise_vals_mult[idx]
    end
end


# Define the system dynamics (double-well) with additive noise
# dx/dt = x - x^3 + noise(t)
noise_add = make_step_noise(tSpan[1], tSpan[2], h, σ_add; seed=1)
# Solve for the memory-free case (α = 1.0) and memory case (α = 0.5)
function F_additive(t, x)
    (x .- x.^3) .+ noise_add(t)   # noise is scalar; broadcast handles array x if needed
end

t_vals, X_nomem = FDEsolver(F_additive, tSpan, y0, α_nomem; h=h)
_,      X_mem   = FDEsolver(F_additive, tSpan, y0, α_mem;   h=h)


# Compute sample statistics for the two cases (exclude initial transient if needed)
steady_index = Int(floor(100.0/h))  # skip initial 100 time units as transient (if any)
X_nomem_steady = X_nomem[steady_index:end]
X_mem_steady   = X_mem[steady_index:end]

# Calculate sample variance in each case
var_nomem = var(X_nomem_steady)
var_mem   = var(X_mem_steady)
println("Sample variance (no memory)   = ", var_nomem)
println("Sample variance (with memory) = ", var_mem)

# Compute autocorrelation functions for each case (to a limited lag)
maxlag = 1000  # compute autocorrelation up to 1000 steps (lag = 10 time units)
function autocorr(data, maxlag)
    n = length(data)
    m = mean(data)
    v = var(data)
    acf = Float64[]
    for τ in 0:maxlag
        # covariance at lag τ
        if τ == 0
            push!(acf, 1.0)  # autocorr(0) = 1
        else
            # covariance for lag τ
            cov_τ = 0.0
            for i in 1:(n-τ)
                cov_τ += (data[i] - m) * (data[i+τ] - m)
            end
            cov_τ /= (n-τ)
            push!(acf, cov_τ / v)
        end
    end
    return acf
end

acf_nomem = autocorr(X_nomem_steady, maxlag)
acf_mem   = autocorr(X_mem_steady, maxlag)

# Print an example of autocorrelation at some lags
println("Autocorrelation at lag 1 (no memory)   = ", acf_nomem[2])
println("Autocorrelation at lag 1 (with memory) = ", acf_mem[2])
println("Autocorrelation at lag 50 (no memory)  = ", acf_nomem[51])
println("Autocorrelation at lag 50 (with memory)= ", acf_mem[51])

maxlag_steps = round(Int, 20000)   # 20000*0.01 = 200 time units
acf_nomem = autocorr(X_nomem_steady, maxlag_steps)
acf_mem   = autocorr(X_mem_steady,   maxlag_steps)

# log–log view to reveal power-law tail
lags_time = (0:maxlag_steps) .* h
plot(lags_time[2:end], acf_nomem[2:end], xscale=:log10, yscale=:log10, label="no memory")
plot!(lags_time[2:end], acf_mem[2:end],   xscale=:log10, yscale=:log10, label="with memory")
xlabel!("lag (time)"); ylabel!("ACF")


# Plot a short segment of the trajectories for visual comparison
t_plot_start = 200.0   # start time for plotting segment
t_plot_end   = 250.0   # end time for plotting segment
i_start = Int(floor(t_plot_start/h))
i_end   = Int(floor(t_plot_end/h))
plot(t_vals[i_start:i_end], X_mem[i_start:i_end], label="With Memory (α=$α_mem)", color=:red)
plot!(t_vals[i_start:i_end], X_nomem[i_start:i_end], label="No Memory (α=1)", color=:blue)
xlabel!("Time"), ylabel!("State x(t)")

function tau_int(acf, h, W_time)
    K = min(length(acf)-1, round(Int, W_time/h))
    return h * (1 + 2*sum(acf[2:K+1]))   # time units
end

println("τ_int(20) no memory   = ", tau_int(acf_nomem, h, 20.0))
println("τ_int(20) with memory = ", tau_int(acf_mem,   h, 20.0))
println("τ_int(50) no memory   = ", tau_int(acf_nomem, h, 50.0))
println("τ_int(50) with memory = ", tau_int(acf_mem,   h, 50.0))

# Solve again with larger noise to observe switching between basins
T_switch     = 5000.0
tSpan_switch = [0.0, T_switch]

σ_add_large = 0.3
# same seed across α to ensure identical noise path for fair comparison
noise_add_large = make_step_noise(tSpan_switch[1], tSpan_switch[2], h, σ_add_large; seed=2)

function F_additive_large(t, x)
    (x .- x.^3) .+ noise_add_large(t)   # NOTE: no dot on the function call
end

t_switch, X_nomem_switch = FDEsolver(F_additive_large, tSpan_switch, y0, α_nomem; h=h)
_,        X_mem_switch   = FDEsolver(F_additive_large, tSpan_switch, y0, α_mem;   h=h)

# Determine dwell times in each stable basin for each case
function compute_dwell_times(t_series, x_series)
    dwell_times = Float64[]
    current_sign = sign(x_series[1])
    start_time = t_series[1]
    for i in 2:length(x_series)
        s = sign(x_series[i])
        if s != current_sign && s != 0
            # state sign has flipped (left one basin for the other)
            push!(dwell_times, t_series[i] - start_time)
            # reset for next basin
            current_sign = s
            start_time = t_series[i]
        end
    end
    return dwell_times
end

dwell_nomem = compute_dwell_times(t_switch, X_nomem_switch)
dwell_mem   = compute_dwell_times(t_switch, X_mem_switch)

println("Number of switches (no memory)   = ", length(dwell_nomem))
println("Number of switches (with memory) = ", length(dwell_mem))
println("Mean dwell time (no memory)      = ", mean(dwell_nomem))
println("Mean dwell time (with memory)    = ", mean(dwell_mem))
println("Dwell time distribution (no memory)   – quartiles: ", quantile(dwell_nomem, [0.25, 0.5, 0.75]))
println("Dwell time distribution (with memory) – quartiles: ", quantile(dwell_mem, [0.25, 0.5, 0.75]))

# Plot histogram of dwell times for memory vs no-memory
histogram(dwell_nomem, bins=20, normalize=:pdf, label="No Memory (α=1)")
histogram!(dwell_mem, bins=20, normalize=:pdf, alpha=0.6, label="With Memory (α=$α_mem)")
xlabel!("Dwell time duration"), ylabel!("Probability Density")

using Plots, Statistics

# --- helper: robust basin state with small deadband δ around 0
basin(x; δ=0.02) = x >  δ ? 1 : (x < -δ ? -1 : 0)

# --- find switch times/indices with deadband
function find_switches(t, x; δ=0.02)
    n = length(x)
    state = basin(x[1]; δ=δ)
    ts, idxs = Float64[], Int[]
    for i in 2:n
        st = basin(x[i]; δ=δ)
        if st != 0 && state != 0 && st != state
            push!(ts, t[i]); push!(idxs, i)
        end
        if st != 0
            state = st
        end
    end
    return ts, idxs
end

# --- split a series into colored segments by basin (NaN-masking)
function split_by_basin(x; δ=0.02)
    y_pos = copy(x);  y_neg = copy(x)
    for i in eachindex(x)
        if basin(x[i]; δ=δ) == 1
            y_neg[i] = NaN
        elseif basin(x[i]; δ=δ) == -1
            y_pos[i] = NaN
        else
            y_pos[i] = NaN; y_neg[i] = NaN
        end
    end
    return y_pos, y_neg
end

# --- plot one panel
function plot_basin_switches(t, x; αlabel="α=?",
                             δ=0.02, ylims=(-1.6, 1.6),
                             c_pos=:blue, c_neg=:red)
    p = plot(legend=:topright, xlabel="time", ylabel="x(t)",
             title=αlabel, ylim=ylims)
    y_pos, y_neg = split_by_basin(x; δ=δ)
    plot!(p, t, y_pos, lw=1.8, color=c_pos, label="right basin (x>0)")
    plot!(p, t, y_neg, lw=1.8, color=c_neg, label="left basin (x<0)")
    hline!(p, [0.0], c=:black, ls=:dash, lw=1.2, label="unstable x=0")

    # ts, idxs = find_switches(t, x; δ=δ)
    # if !isempty(ts)
    #     vline!(p, ts, c=:gray, ls=:dash, alpha=0.5, label="switch")
    #     scatter!(p, t[idxs], x[idxs], m=:diamond, ms=4.5, c=:black,
    #              label="crossing")
    # end
    return p
end

# --- build the two-panel figure (use your existing t_switch, X_*_switch)
p1 = plot_basin_switches(t_switch, X_nomem_switch;
                         αlabel="No memory (α=1)")
p2 = plot_basin_switches(t_switch, X_mem_switch;
                         αlabel="With memory (α=$α_mem)")
plot(p1, p2, layout=(2,1), size=(950,650))

# ---------- parameters for robust detection ----------
# x_core   = 0.7          # core threshold (inside the valley)
# τ_hold   = 0.20         # minimum hold time to count entry (time units)
# hstep(t) = t[2]-t[1]    # helper

# # Label: +1 in right core, -1 in left core, 0 otherwise
# label_core(x; xc=x_core) = x ≥ xc ?  1 : (x ≤ -xc ? -1 : 0)

# # Find committed entries into basin cores with a hold-time requirement
# function core_entries(t, x; xc=x_core, τhold=τ_hold)
#     n = length(x)
#     Δ = hstep(t)
#     hold = max(1, round(Int, τhold/Δ))
#     lab = [label_core(x[i]; xc=xc) for i in 1:n]

#     idx, lab_keep = Int[], Int[]
#     i = 1
#     while i ≤ n-hold+1
#         ℓ = lab[i]
#         if ℓ == 0
#             i += 1
#             continue
#         end
#         # must remain in the same core for 'hold' steps
#         ok = all(lab[j] == ℓ for j in i:i+hold-1)
#         if ok
#             push!(idx, i); push!(lab_keep, ℓ)
#             # skip forward until we leave this core (prevents dense duplicates)
#             k = i + hold
#             while k ≤ n && lab[k] == ℓ
#                 k += 1
#             end
#             i = k
#         else
#             i += 1
#         end
#     end
#     return idx, lab_keep
# end

# # Dwell times between committed entries into opposite cores
# function dwell_times_committed(t, x; xc=x_core, τhold=τ_hold)
#     idx, labs = core_entries(t, x; xc=xc, τhold=τhold)
#     dwell = Float64[]
#     for k in 2:length(idx)
#         if labs[k] != labs[k-1]                 # only opposite-core entries
#             push!(dwell, t[idx[k]] - t[idx[k-1]])
#         end
#     end
#     return dwell, idx, labs
# end

# # Compute committed-switch dwell times for both cases
# dwell_nomem_c, idx_nomem, labs_nomem = dwell_times_committed(t_switch, X_nomem_switch)
# dwell_mem_c,   idx_mem,   labs_mem   = dwell_times_committed(t_switch, X_mem_switch)

# println("Committed switches (no mem)   = ", length(dwell_nomem_c))
# println("Committed switches (with mem) = ", length(dwell_mem_c))
# println("Mean committed dwell (no mem) = ", mean(dwell_nomem_c))
# println("Mean committed dwell (mem)    = ", mean(dwell_mem_c))

# # -------- plots --------

# # 1) Two-panel time series with cores and committed switch markers
# function plot_switch_panels(t, x; title_str="", xc=x_core, τhold=τ_hold)
#     p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
#     plot!(p, t, x, lw=0.8, color=:gray70, label=false)
#     hline!(p, [0.0], ls=:dash, c=:black, label="unstable x=0")
#     hline!(p, [ xc], ls=:dot,  c=:blue,  label="core +x")
#     hline!(p, [-xc], ls=:dot,  c=:red,   label="core -x")
#     idx, labs = core_entries(t, x; xc=xc, τhold=τhold)
#     if !isempty(idx)
#         scatter!(p, t[idx], x[idx], m=:diamond, ms=1, c=:purple, label="committed entry")
#     end
#     p
# end

# p1 = plot_switch_panels(t_switch, X_nomem_switch; title_str="No memory (α=1.0)")
# p2 = plot_switch_panels(t_switch, X_mem_switch;   title_str="With memory (α=$α_mem)")
# plot(p1, p2, layout=(2,1), size=(950,650))


# ---------- parameters ----------
x_core   = 0.7          # core threshold (inside the valley)
τ_hold   = 0.50         # minimum hold time to count entry (time units)
hstep(t) = t[2]-t[1]    # helper

# Label: +1 in right core, -1 in left core, 0 otherwise
label_core(x; xc=x_core) = x ≥ xc ?  1 : (x ≤ -xc ? -1 : 0)

# EARLY ANCHOR: first committed entry (optional helper)
function first_committed_entry(t, x; xc=x_core, τhold=τ_hold)
    Δ = t[2]-t[1]; hold = max(1, round(Int, τhold/Δ))
    lab = [label_core(x[i]; xc=xc) for i in eachindex(x)]
    for i in 1:length(x)-hold+1
        ℓ = lab[i]
        if ℓ != 0 && all(lab[j]==ℓ for j in i:i+hold-1)
            return i, ℓ    # index and which core (+1 or -1)
        end
    end
    return nothing, 0
end

# Committed entries into basin cores (hold-time required)
function core_entries(t, x; xc=x_core, τhold=τ_hold)
    n = length(x)
    Δ = hstep(t)
    hold = max(1, round(Int, τhold/Δ))
    lab = [label_core(x[i]; xc=xc) for i in 1:n]

    idx, lab_keep = Int[], Int[]
    i = 1
    while i ≤ n - hold + 1
        ℓ = lab[i]
        if ℓ == 0
            i += 1
            continue
        end
        ok = all(lab[j] == ℓ for j in i:i+hold-1)
        if ok
            push!(idx, i); push!(lab_keep, ℓ)
            k = i + hold
            while k ≤ n && lab[k] == ℓ
                k += 1
            end
            i = k
        else
            i += 1
        end
    end
    return idx, lab_keep
end

# Dwell times between committed entries into opposite cores,
# with optional anchoring at the very first committed entry.
function committed_durations_censored(t, x; xc=x_core, τhold=τ_hold, anchor_at_first::Bool=true)
    if anchor_at_first
        i0, _ = first_committed_entry(t, x; xc=xc, τhold=τhold)
        if i0 !== nothing
            t = t[i0:end]; x = x[i0:end]   # clip series to start at first committed entry
        end
    end

    idx_all, labs_all = core_entries(t, x; xc=xc, τhold=τhold)
    times, cens = Float64[], Int[]   # 0 = observed, 1 = right-censored

    # completed dwells only (opposite-core entries)
    for k in 2:length(idx_all)
        if labs_all[k] != labs_all[k-1]
            push!(times, t[idx_all[k]] - t[idx_all[k-1]])
            push!(cens, 0)
        end
    end
    # add the final ongoing dwell as right-censored
    if !isempty(idx_all)
        push!(times, t[end] - t[idx_all[end]])
        push!(cens, 1)
    end
    return times, cens
end


# Convenience: get only the switch indices for plotting
committed_switch_indices(t, x; xc=x_core, τhold=τ_hold) =
    last(dwell_times_committed(t, x; xc=xc, τhold=τhold))


    # List all committed core-entry times and labels
idx_all, labs_all = core_entries(t_switch, X_nomem_switch; xc=x_core, τhold=τ_hold)
E = t_switch[idx_all]
println("Committed entries (no-mem):")
for k in 1:length(E)
    println(@sprintf("  k=%2d  t=%.3f  label=%+d", k, E[k], labs_all[k]))
end

println("\nCommitted dwells (no-mem):")
for k in 2:length(E)
    if labs_all[k] != labs_all[k-1]
        Δt = E[k] - E[k-1]
        println(@sprintf("  %2d: %.3f time units  (~%d steps)", k-1, Δt, round(Int, Δt/0.01)))
    end
end

function print_dwells(t, x; xc=x_core, τhold=τ_hold, h=0.01)
    idx_all, labs_all = core_entries(t, x; xc=xc, τhold=τhold)
    E = t[idx_all]
    println("Committed dwells:")
    for k in 2:length(E)
        if labs_all[k] != labs_all[k-1]
            Δt = E[k] - E[k-1]
            println(@sprintf("  %.3f time units  (~%d steps)", Δt, round(Int, Δt/h)))
        end
    end
end

# Example:
print_dwells(t_switch, X_nomem_switch; h=0.01)
print_dwells(t_switch, X_mem_switch;   h=0.01)


# -------------- plotting (only true committed switches) --------------

function plot_switch_panels(t, x; title_str="", xc=x_core, τhold=τ_hold,
                            show_vlines=true, ms=5)
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x, lw=0.8, color=:gray70, label=false)
    hline!(p, [0.0], ls=:dash, c=:black, label="unstable x=0")
    hline!(p, [ xc], ls=:dot,  c=:blue,  label="core +x")
    hline!(p, [-xc], ls=:dot,  c=:red,   label="core -x")

    sw_idx = committed_switch_indices(t, x; xc=xc, τhold=τhold)
    if !isempty(sw_idx)
        if show_vlines
            vline!(p, t[sw_idx], c=:gray, ls=:dash, alpha=0.35, label=false)
        end
        scatter!(p, t[sw_idx], x[sw_idx], m=:diamond, ms=ms, c=:purple,
                 label="committed switch")
    end
    return p
end

# -------- usage ----------
dwell_nomem_c, swidx_nomem = dwell_times_committed(t_switch, X_nomem_switch)
dwell_mem_c,   swidx_mem   = dwell_times_committed(t_switch, X_mem_switch)

println("Committed switches (no mem)   = ", length(swidx_nomem))
println("Committed switches (with mem) = ", length(swidx_mem))
println("Mean committed dwell (no mem) = ", mean(dwell_nomem_c))
println("Mean committed dwell (mem)    = ", mean(dwell_mem_c))

p1 = plot_switch_panels(t_switch, X_nomem_switch; title_str="No memory (α=1.0)")
p2 = plot_switch_panels(t_switch, X_mem_switch;   title_str="With memory (α=$(α_mem))")
plot(p1, p2, layout=(2,1), size=(950,650), legendposition=:topleft)

# 2) Dwell-time histograms with logarithmic bins
function logbins(x; nb=15)
    xs = sort(x); lo = minimum(xs[xs .> 0]); hi = maximum(xs)
    exp.(range(log(lo), log(hi), length=nb))
end
bins_nomem = logbins(dwell_nomem_c)
bins_mem   = logbins(dwell_mem_c)

histogram(dwell_nomem_c, bins=bins_nomem, normalize=:pdf, label="No memory (committed)", alpha=0.6)
histogram!(dwell_mem_c,   bins=bins_mem,   normalize=:pdf, label="With memory (committed)", alpha=0.6)
xlabel!("Dwell time"); ylabel!("PDF"); title!("Committed dwell-time distributions")

##lets correct dwell time
# Durations + censoring from committed entries
# returns: times (durations), cens (0=event observed, 1=right-censored)
function committed_durations_censored(t, x; xc=x_core, τhold=τ_hold)
    idx_all, labs_all = core_entries(t, x; xc=xc, τhold=τhold)
    times, cens = Float64[], Int[]
    for k in 2:length(idx_all)
        if labs_all[k] != labs_all[k-1]
            push!(times, t[idx_all[k]] - t[idx_all[k-1]]); push!(cens, 0)
        end
    end
    if !isempty(idx_all)
        # add the final ongoing dwell as right-censored
        push!(times, t[end] - t[idx_all[end]]); push!(cens, 1)
    end
    return times, cens
end

# Kaplan–Meier estimator of survival S(t)
function kaplan_meier(times::Vector{<:Real}, cens::Vector{<:Integer})
    p = sortperm(times); t = collect(times[p]); c = collect(cens[p])
    n = length(t); at_risk = n
    T, S = Float64[], Float64[]
    s_prev = 1.0; i = 1
    while i ≤ n
        ti = t[i]; d = 0; cc = 0
        while i ≤ n && t[i] == ti
            (c[i] == 0) ? (d += 1) : (cc += 1); i += 1
        end
        if d > 0
            s_prev *= (at_risk - d)/at_risk
            push!(T, ti); push!(S, s_prev)
        end
        at_risk -= (d + cc)
    end
    return T, S
end

# Restricted mean (area under survival to W)
function restricted_mean(T::Vector{<:Real}, S::Vector{<:Real}, W::Real)
    # piecewise-constant survival between event times
    tgrid = [0.0; T; W]
    sgrid = [1.0; S; (isempty(S) ? 1.0 : S[end])]
    rmst = 0.0
    for j in 1:length(tgrid)-1
        Δ = tgrid[j+1] - tgrid[j]
        rmst += sgrid[j]*Δ
    end
    rmst
end

times_no, cens_no = committed_durations_censored(t_switch, X_nomem_switch)
times_me, cens_me = committed_durations_censored(t_switch, X_mem_switch)

Tno, Sno = kaplan_meier(times_no, cens_no)
Tme, Sme = kaplan_meier(times_me, cens_me)

println("KM median (no mem): ", Tno[findfirst(<(0.5), Sno)] )
println("KM median (mem):    ", Tme[findfirst(<(0.5), Sme)] )

println("no-mem dwell (steps): ",
        round.(Int, dwell_nomem_c ./ h))
println("mem dwell (steps): ",
        round.(Int, dwell_mem_c   ./ h))
println("no-mem min/med/max (time): ",
        minimum(dwell_nomem_c), " / ",
        median(dwell_nomem_c), " / ",
        maximum(dwell_nomem_c))


println("Restricted mean to W=200 (no mem): ", restricted_mean(Tno, Sno, 200.0))
println("Restricted mean to W=200 (mem):    ", restricted_mean(Tme, Sme, 200.0))

rate_no = length(times_no[cens_no .== 0]) / (t_switch[end] - t_switch[1])
rate_me = length(times_me[cens_me .== 0]) / (t_switch[end] - t_switch[1])
println("Committed switch rate  no mem = ", rate_no)
println("Committed switch rate  mem    = ", rate_me)

# Occupancy fraction based on instantaneous core label (no hold)
lab = x -> label_core(x; xc=x_core)
frac_in_core(x) = mean([lab(xi)!=0 for xi in x])
println("Occupancy in cores no mem = ", frac_in_core(X_nomem_switch))
println("Occupancy in cores mem    = ", frac_in_core(X_mem_switch))

# 3) Survival curves (1 - ECDF) on log–log axes to show tails
function survival_curve(x; npts=20)
    xs = sort(x)
    n = length(xs)
    # stepwise survival at event times: S(t_k) = 1 - k/n
    tk = xs
    Sk = 1 .- (1:n)./n
    (tk, Sk)
end
tk1, S1 = survival_curve(dwell_nomem_c)
tk2, S2 = survival_curve(dwell_mem_c)
# S1 = collect(S1)        # convert StepRangeLen → Vector{Float64}
# S1[end] = S1[end] + 1e-10   # or .+ (both work, but + is cleaner)
# S2 = collect(S2)        # convert StepRangeLen → Vector{Float64}
# S2[end] = S2[end] + 1e-10   # or .+ (both work, but + is cleaner)

# plot(tk1, S1, xscale=:log10, yscale=:log10, lw=2, label="No memory")
# plot!(tk2, S2, xscale=:log10, yscale=:log10, lw=2, label="With memory")
plot(tk1, S1, lw=2, label="No memory")
plot!(tk2, S2, lw=2, label="With memory")
xlabel!("t"); ylabel!("Survival 1-ECDF"); title!("Committed dwell-time survival")



###test: Linearized check to match theory even more clearly: 
# λ = 2.0
# function F_add_lin(t, y)   # D^α y = -λ y + ξ(t)
#     return -λ.*y .+ noise_func_add.(t)
# end
# t_lin, Y_nomem = FDEsolver(F_add_lin, [0.0, 1000.0], 0.0, 1.0;  h=h)
# _,     Y_mem   = FDEsolver(F_add_lin, [0.0, 1000.0], 0.0, 0.7;  h=h)

# acf_nomem_lin = autocorr(Y_nomem[steady_index:end], 20000)
# acf_mem_lin   = autocorr(Y_mem[steady_index:end],   20000)
# # Plot as above; compute τ_int as above
# # log–log view to reveal power-law tail
# lags_time = (0:maxlag_steps) .* h
# mask = (lags_time .> 0) .& (acf_nomem_lin .> 0) .& (acf_mem_lin .> 0)

# plot(lags_time[mask], acf_nomem_lin[mask], xscale=:log10,yscale=:log10, label="no memory")
# plot!(lags_time[mask], acf_mem_lin[mask],  label="with memory")
# xlabel!("lag (time)"); ylabel!("ACF")

# plot(lags_time[2:end], acf_nomem_lin[2:end], xscale=:log10,yscale=:log10, label="no memory")
# plot!(lags_time[2:end], acf_mem_lin[2:end],  label="with memory")
# xlabel!("lag (time)"); ylabel!("ACF")


##
# block-average by factor m (low-pass + decimate)
function block_average(x, m)
    n = length(x) ÷ m
    [mean(view(x, (i-1)*m+1:i*m)) for i in 1:n]
end

m_list = [1, 5, 10, 20, 50, 100, 150, 200]      # 1 = raw; larger = coarser
W = 50.0                         # window (time units) for tau_int
function tau_int_from_series(x, h, W)
    maxlag = round(Int, W/h)
    # (use your autocorr() but limit maxlag ≪ length)
    acf = autocorr(x, maxlag)
    h * (1 + 2*sum(acf[2:end]))  # same formula you already use
end

for m in m_list
    Xn = block_average(X_nomem_steady, m)
    Xm = block_average(X_mem_steady,   m)
    h_eff = h*m
    println("m=$m  τ_int no-mem = ", tau_int_from_series(Xn, h_eff, W),
            "   τ_int mem = ",      tau_int_from_series(Xm, h_eff, W))
end


using FFTW
function acf_fft(x; maxlag::Int)
    xz = x .- mean(x)
    n  = length(xz)
    n2 = 1
    while n2 < 2n; n2 <<= 1; end
    fx = rfft([xz; zeros(n2-n)])
    s  = real(irfft(abs.(fx).^2, n2))
    c  = s[1:n] ./ (n .- (0:n-1))         # unbiased
    ac = c ./ c[1]
    ac[1:maxlag+1]
end

# assumes you already have: X_nomem_steady, X_mem_steady, h, autocorr(), tau_int()

function block_average(x, m)
    n = length(x) ÷ m
    [mean(view(x, (i-1)*m+1:i*m)) for i in 1:n]
end

m_list = [1,5,10,20,50,100,150,200]
W = 50.0  # time window for tau_int

tau_no = Float64[]
tau_me = Float64[]
for m in m_list
    Xn = block_average(X_nomem_steady, m)
    Xm = block_average(X_mem_steady,   m)
    h_eff = h*m
    acf_n = autocorr(Xn, round(Int, W/h_eff))
    acf_m = autocorr(Xm, round(Int, W/h_eff))
    push!(tau_no, h_eff * (1 + 2*sum(acf_n[2:end])))
    push!(tau_me, h_eff * (1 + 2*sum(acf_m[2:end])))
end

# tau_int vs coarse graining factor
plot(m_list, tau_no, lw=2, label="no memory")
plot!(m_list, tau_me, lw=2, label="with memory")
xlabel!("block size m"); ylabel!("tau_int over window W"); title!("Effective persistence vs coarse graining")

# also show the ratio to highlight the cross over
plot(m_list, tau_me ./ tau_no, lw=2, label="(with memory)/(no memory)")
hline!([1.0], ls=:dash, label="= 1")
xlabel!("block size m"); ylabel!("ratio of tau_int"); title!("Cross over of effective persistence")


var_no = Float64[]
var_me = Float64[]
for m in m_list
    Xn = block_average(X_nomem_steady, m)
    Xm = block_average(X_mem_steady,   m)
    push!(var_no, var(Xn))
    push!(var_me, var(Xm))
end
plot(m_list, var_no, lw=2, label="no memory")
plot!(m_list, var_me, lw=2, label="with memory")
xlabel!("block size m"); ylabel!("variance"); title!("Variance after coarse graining")


