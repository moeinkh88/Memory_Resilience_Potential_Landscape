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


# ---------- switching run (long span, larger noise) ----------
T_switch     = 5000.0
tSpan_switch = [0.0, T_switch]
σ_add_large  = 0.35

# identical noise path for both α
noise_add_large = make_step_noise(tSpan_switch[1], tSpan_switch[2], h, σ_add_large; seed=2)
F_additive_large(t, x) = (x .- x.^3) .+ noise_add_large(t)

t_switch, X_nomem_switch = FDEsolver(F_additive_large, tSpan_switch, y0, α_nomem; h=h)
_,        X_mem_switch   = FDEsolver(F_additive_large, tSpan_switch, y0, α_mem;   h=h)

# ---------- committed switching in *steps* (no time units) ----------
x_core      = 0.9      # basin core threshold
hold_steps  = 20       # must stay in core this many *steps* to count entry (20 steps = 0.2 time if h=0.01)

label_core(x; xc=x_core) = x ≥ xc ?  1 : (x ≤ -xc ? -1 : 0)

# committed entries into cores with a "hold_steps" requirement (returns indices and labels)
function core_entries_steps(x; xc=x_core, hold_steps::Int=hold_steps)
    n = length(x)
    lab = [label_core(x[i]; xc=xc) for i in 1:n]
    idx, labs = Int[], Int[]
    i = 1
    while i ≤ n - hold_steps + 1
        ℓ = lab[i]
        if ℓ == 0
            i += 1
            continue
        end
        ok = true
        @inbounds for j in i:i+hold_steps-1
            if lab[j] != ℓ; ok = false; break; end
        end
        if ok
            push!(idx, i); push!(labs, ℓ)
            k = i + hold_steps
            while k ≤ n && lab[k] == ℓ
                k += 1
            end
            i = k
        else
            i += 1
        end
    end
    return idx, labs
end

# dwell lengths in *steps* between committed entries into opposite cores
function dwell_steps_committed(x; xc=x_core, hold_steps::Int=hold_steps)
    idx_all, labs_all = core_entries_steps(x; xc=xc, hold_steps=hold_steps)
    dw_steps, switch_idx = Int[], Int[]
    for k in 2:length(idx_all)
        if labs_all[k] != labs_all[k-1]     # opposite core => a real switch
            push!(dw_steps, idx_all[k] - idx_all[k-1])
            push!(switch_idx, idx_all[k])   # mark the entry that completed the switch
        end
    end
    return dw_steps, switch_idx
end

# compute dwell (steps) and committed switch indices
dwell_nomem_steps, swidx_nomem = dwell_steps_committed(X_nomem_switch)
dwell_mem_steps,   swidx_mem   = dwell_steps_committed(X_mem_switch)

println("Committed switches (no mem)   = ", length(swidx_nomem))
println("Committed switches (with mem) = ", length(swidx_mem))
println("Mean committed dwell (steps)  no mem = ", mean(dwell_nomem_steps))
println("Mean committed dwell (steps)  mem    = ", mean(dwell_mem_steps))

# ---------- plots ----------

# two-panel time series with cores and committed switches (diamonds); x-axis in time for readability
function plot_switch_panels_steps(t, x, swidx; title_str="")
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x, lw=0.8, color=:gray70, label=false)
    hline!(p, [0.0], ls=:dash, c=:black, label="unstable x=0")
    hline!(p, [ x_core], ls=:dot,  c=:blue,  label="core +x")
    hline!(p, [-x_core], ls=:dot,  c=:red,   label="core -x")
    if !isempty(swidx)
        vline!(p, t[swidx], c=:gray, ls=:dash, alpha=0.35, label=false)
        scatter!(p, t[swidx], x[swidx], m=:diamond, ms=5, c=:purple, label="committed switch")
    end
    p
end

p1 = plot_switch_panels_steps(t_switch, X_nomem_switch, swidx_nomem; title_str="No memory (α=1.0)")
p2 = plot_switch_panels_steps(t_switch, X_mem_switch,   swidx_mem;   title_str="With memory (α=$(α_mem))")
plot(p1, p2, layout=(2,1), size=(950,650), legendposition=:topleft)

# histogram of dwell lengths in steps (log bins help if wide range)
using StatsBase
function logbins_steps(x::Vector{<:Integer}; nb=10)
    xs = sort(x); lo = minimum(xs[xs .> 0]); hi = maximum(xs)
    # ensure integer bin edges
    round.(Int, exp.(range(log(lo), log(hi), length=nb)))
end
bins1 = logbins_steps(dwell_nomem_steps)
bins2 = logbins_steps(dwell_mem_steps)

histogram(dwell_nomem_steps, bins=bins1, normalize=:pdf, alpha=0.6, label="No memory")
histogram!(dwell_mem_steps,   bins=bins2, normalize=:pdf, alpha=0.6, label="With memory")
xlabel!("Dwell length (steps)"); ylabel!("PDF"); title!("Committed dwell-time distributions (steps)")

# survival curve (1-ECDF) in steps
function survival_steps(d::Vector{<:Integer})
    xs = sort(d); n = length(xs)
    tk = Float64.(xs)
    Sk = 1 .- (1:n)./n
    (tk, Sk)
end
tk1, S1 = survival_steps(dwell_nomem_steps)
tk2, S2 = survival_steps(dwell_mem_steps)
plot(tk1, S1, lw=2, label="No memory")
plot!(tk2, S2, lw=2, label="With memory")
xlabel!("dwell length (steps)"); ylabel!("Survival 1-ECDF"); title!("Committed dwell survival (steps)")




##optional

using StatsBase, Printf, Plots

# --- helpers for nicer, comparable histograms ---

"Freedman–Diaconis bin width (in steps) on pooled data; rounded to an int and floored at min_bw."
function fd_binwidth_steps(data::AbstractVector{<:Real}; min_bw::Int=20)
    bw = 2iqr(data) / length(data)^(1/3)
    max(min_bw, max(1, round(Int, bw)))
end

"Common linear bin edges for two samples (in steps), using FD rule on pooled data."
function common_linear_edges_steps(a::Vector{<:Integer}, b::Vector{<:Integer}; min_bw::Int=10)
    pooled = vcat(a, b)
    bw = fd_binwidth_steps(pooled; min_bw=min_bw)
    lo = floor(Int, minimum(pooled))
    hi = ceil(Int,  maximum(pooled))
    collect(lo:bw:hi)
end

"Common log-spaced edges (in steps) for two samples; good for heavy tails."
function common_log_edges_steps(a::Vector{<:Integer}, b::Vector{<:Integer}; nbins::Int=18)
    pooled = vcat(a, b)
    lo = minimum(pooled[pooled .> 0])              # avoid log(0)
    hi = maximum(pooled)
    edges = round.(Int, exp.(range(log(lo), log(hi), length=nbins)))
    unique(sort(edges))
end

"Plot two overlaid histograms with shared edges, step outlines, transparent fill, and median lines."
function plot_dwell_hists_steps(a::Vector{<:Integer}, b::Vector{<:Integer};
                                edges::Vector{Int}, title_str::String)
    med_a, med_b = median(a), median(b)

    p = histogram(a; bins=edges, normalize=:pdf,
                  fillalpha=0.25, linealpha=0.9, linewidth=2,
                  linecolor=:blue, fillcolor=:blue, label="No memory")
    histogram!(p, b; bins=edges, normalize=:pdf,
               fillalpha=0.25, linealpha=0.9, linewidth=2,
               linecolor=:orange, fillcolor=:orange, label="With memory")

    vline!(p, [med_a], lc=:blue,   ls=:dash, lw=2, label="median (no mem)")
    vline!(p, [med_b], lc=:orange, ls=:dash, lw=2, label="median (mem)")

    xlabel!(p, "Dwell length (steps)")
    ylabel!(p, "PDF")
    title!(p, title_str)
    p
end

# --- make the nicer plots (pick one or show both) ---

# 1) Linear bins via FD rule (good for central mass)
edges_lin = common_linear_edges_steps(dwell_nomem_steps, dwell_mem_steps; min_bw=25)
p_lin = plot_dwell_hists_steps(dwell_nomem_steps, dwell_mem_steps;
                               edges=edges_lin,
                               title_str="Committed dwell-time distributions \n (steps, linear bins)")

# 2) Log-spaced bins (good for tail visualisation)
edges_log = common_log_edges_steps(dwell_nomem_steps, dwell_mem_steps; nbins=16)
p_log = plot_dwell_hists_steps(dwell_nomem_steps, dwell_mem_steps;
                               edges=edges_log,
                               title_str="Committed dwell-time distributions \n (steps, log-spaced bins)")

plot(p_lin, p_log, layout=(1,2), size=(1100,420), legend=:topright)


#
using KernelDensity, StatsBase, Plots

# shared linear bins (equal width) for a light background histogram
lo = min(minimum(dwell_nomem_steps), minimum(dwell_mem_steps))
hi = max(maximum(dwell_nomem_steps), maximum(dwell_mem_steps))
edges = collect(floor(Int, lo):100:ceil(Int, hi))   # 100-step bins; tweak as you like

histogram(dwell_nomem_steps; bins=edges, normalize=:pdf,
          fillalpha=0.15, linealpha=0.3, label="No memory (hist)")
histogram!(dwell_mem_steps;   bins=edges, normalize=:pdf,
           fillalpha=0.15, linealpha=0.3, label="With memory (hist)")

# KDE curves (on the raw steps)
kd_no = kde(Float64.(dwell_nomem_steps))
kd_me = kde(Float64.(dwell_mem_steps))
plot!(kd_no.x, kd_no.density; lw=3, label="No memory (KDE)")
plot!(kd_me.x, kd_me.density; lw=3, label="With memory (KDE)")

xlabel!("Dwell length (steps)"); ylabel!("PDF")
title!("Committed dwell-time distributions (steps)")

vline!([median(dwell_nomem_steps)], ls=:dash, lc=:blue,  label="median no-mem")
vline!([median(dwell_mem_steps)],   ls=:dash, lc=:orange, label="median mem")

##
using KernelDensity, Plots

# KDE in log-domain, then map back: f_X(x) = f_Y(log x) / x
kd_no_log = kde(log.(Float64.(dwell_nomem_steps)))
kd_me_log = kde(log.(Float64.(dwell_mem_steps)))

x_no = exp.(kd_no_log.x);  pdf_no = kd_no_log.density ./ x_no
x_me = exp.(kd_me_log.x);  pdf_me = kd_me_log.density ./ x_me

# (optional) faint equal-width histogram in background
lo = min(minimum(dwell_nomem_steps), minimum(dwell_mem_steps))
hi = max(maximum(dwell_nomem_steps), maximum(dwell_mem_steps))
edges = collect(floor(Int, lo):100:ceil(Int, hi))
histogram(dwell_nomem_steps; bins=edges, normalize=:pdf,
          fillalpha=0.1, linealpha=0.0, label="No memory (hist)")
histogram!(dwell_mem_steps;   bins=edges, normalize=:pdf,
           fillalpha=0.1, linealpha=0.0, label="With memory (hist)")

# smooth curves overlaid
plot!(x_no, pdf_no; lw=3, label="No memory (log-KDE)")
plot!(x_me, pdf_me; lw=3, label="With memory (log-KDE)")

xlabel!("Dwell length (steps)"); ylabel!("PDF")
title!("Committed dwell-time distributions (steps, log-KDE)")

vline!([median(dwell_nomem_steps)], ls=:dash, lc=:blue,  label="median no-mem")
vline!([median(dwell_mem_steps)],   ls=:dash, lc=:orange, label="median mem")
