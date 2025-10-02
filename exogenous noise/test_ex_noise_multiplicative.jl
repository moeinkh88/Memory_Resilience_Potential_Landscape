using FdeSolver, Random, Statistics, StatsBase, KernelDensity, Plots
using Printf

# -------------------------- setup --------------------------

Random.seed!(123)

# Model parameters
α_mem   = 0.8     # fractional order (memory case, 0<α<1)
α_nomem = 1.0     # standard case (no memory, α=1)
y0      = 1.0     # initial state (near one stable equilibrium)

# Time grid
T     = 1000.0
h     = 0.01
tSpan = [0.0, T]

# Multiplicative noise base intensity (for the short run analyses)
σ_mult = 0.2

# State-dependent noise amplitude: higher near unstable x≈0, lower near ±1
noise_amp = x -> begin
    amp = 1.0 - 0.8*x^2
    amp < 0 ? 0.0 : amp
end

"Return a callable step-noise ξ(t) for [t0,t1] with step h and sd σ (white-noise approx)."
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=123)
    N   = Int(ceil((t1 - t0)/h))              # cover full span
    rng = MersenneTwister(seed)
    vals = (σ/√h) .* randn(rng, N)            # piecewise-constant increments
    t0_, h_ = float(t0), float(h)
    return (t::Real)-> begin
        idx = clamp(Int(floor((t - t0_)/h_)) + 1, 1, N)
        return vals[idx]
    end
end

# One fixed noise path used for *both* α to keep driving noise comparable
ξ_mult = make_step_noise(tSpan[1], tSpan[2], h, σ_mult; seed=1)

# System with multiplicative noise: dx/dt = x - x^3 + noise_amp(x)*ξ(t)
F_multiplicative(t, x) = (x .- x.^3) .+ noise_amp.(x) .* ξ_mult(t)

# -------------------------- simulate (short run) --------------------------

t_vals, X_nomem = FDEsolver(F_multiplicative, tSpan, y0, α_nomem; h=h)
_,      X_mem   = FDEsolver(F_multiplicative, tSpan, y0, α_mem;   h=h)

# Remove transient
steady_index   = Int(floor(100.0/h))
X_nomem_steady = X_nomem[steady_index:end]
X_mem_steady   = X_mem[steady_index:end]

# Sample variance
var_nomem = var(X_nomem_steady)
var_mem   = var(X_mem_steady)
@printf "Sample variance (no memory)   = %.6f\n" var_nomem
@printf "Sample variance (with memory) = %.6f\n" var_mem

# Autocorrelation
function autocorr_vec(data, maxlag::Int)
    n = length(data)
    m = mean(data)
    v = var(data)
    acf = Vector{Float64}(undef, maxlag+1)
    acf[1] = 1.0
    @inbounds for τ in 1:maxlag
        s = 0.0
        for i in 1:(n-τ)
            s += (data[i] - m) * (data[i+τ] - m)
        end
        acf[τ+1] = (s / (n-τ)) / v
    end
    return acf
end

maxlag_steps = 20000              # 20000*0.01 = 200 time units
acf_nomem = autocorr_vec(X_nomem_steady, maxlag_steps)
acf_mem   = autocorr_vec(X_mem_steady,   maxlag_steps)

# log–log ACF plot
lags_time = (0:maxlag_steps) .* h
plot(lags_time[2:end], acf_nomem[2:end], xscale=:log10, yscale=:log10, label="no memory")
plot!(lags_time[2:end], acf_mem[2:end],   xscale=:log10, yscale=:log10, label="with memory")
xlabel!("lag (time)"); ylabel!("ACF"); title!("ACF (log–log) under multiplicative noise")

# Short segment plot
t_plot_start = 200.0
t_plot_end   = 250.0
i_start = Int(floor(t_plot_start/h))
i_end   = Int(floor(t_plot_end/h))
plot(t_vals[i_start:i_end], X_mem[i_start:i_end],   label="With Memory (α=$(α_mem))",  color=:red)
plot!(t_vals[i_start:i_end], X_nomem[i_start:i_end], label="No Memory (α=1)",          color=:blue)
xlabel!("Time"); ylabel!("x(t)"); title!("Trajectories (multiplicative noise)")

# Integrated autocorrelation time over a finite window W_time
function tau_int(acf::AbstractVector{<:Real}, h::Real, W_time::Real)
    K = min(length(acf)-1, round(Int, W_time/h))
    return h * (1 + 2*sum(acf[2:K+1]))
end

@printf "τ_int(20) no memory   = %.6f\n" tau_int(acf_nomem, h, 20.0)
@printf "τ_int(20) with memory = %.6f\n" tau_int(acf_mem,   h, 20.0)
@printf "τ_int(50) no memory   = %.6f\n" tau_int(acf_nomem, h, 50.0)
@printf "τ_int(50) with memory = %.6f\n" tau_int(acf_mem,   h, 50.0)

# -------------------------- switching experiment --------------------------

T_switch     = 5000.0
tSpan_switch = [0.0, T_switch]
σ_mult_large = 0.5   # increase if you get too few switches (e.g., try 0.5)

ξ_mult_large = make_step_noise(tSpan_switch[1], tSpan_switch[2], h, σ_mult_large; seed=2)
F_multiplicative_large(t, x) = (x .- x.^3) .+ noise_amp.(x) .* ξ_mult_large(t)

t_switch, X_nomem_switch = FDEsolver(F_multiplicative_large, tSpan_switch, y0, α_nomem; h=h)
_,        X_mem_switch   = FDEsolver(F_multiplicative_large, tSpan_switch, y0, α_mem;   h=h)

# Core-based committed switching in *steps*
x_core     = 0.9
hold_steps = 20

label_core(x; xc=x_core) = x ≥ xc ?  1 : (x ≤ -xc ? -1 : 0)

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

function dwell_steps_committed(x; xc=x_core, hold_steps::Int=hold_steps)
    idx_all, labs_all = core_entries_steps(x; xc=xc, hold_steps=hold_steps)
    dw_steps, switch_idx = Int[], Int[]
    for k in 2:length(idx_all)
        if labs_all[k] != labs_all[k-1]
            push!(dw_steps, idx_all[k] - idx_all[k-1])
            push!(switch_idx, idx_all[k])
        end
    end
    return dw_steps, switch_idx
end

dwell_nomem_steps, swidx_nomem = dwell_steps_committed(X_nomem_switch)
dwell_mem_steps,   swidx_mem   = dwell_steps_committed(X_mem_switch)

println("Committed switches (no mem)   = ", length(swidx_nomem))
println("Committed switches (with mem) = ", length(swidx_mem))
println("Mean committed dwell (steps)  no mem = ", mean(dwell_nomem_steps))
println("Mean committed dwell (steps)  mem    = ", mean(dwell_mem_steps))

# -------------------------- plots: switching panels --------------------------

function plot_switch_panels_steps(t, x, swidx; title_str="")
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x, lw=0.8, color=:gray70, label=false)
    hline!(p, [0.0],  ls=:dash, c=:black, label="unstable x=0")
    hline!(p, [ x_core], ls=:dot,  c=:blue,  label="core +x")
    hline!(p, [-x_core], ls=:dot,  c=:red,   label="core -x")
    if !isempty(swidx)
        vline!(p, t[swidx], c=:gray, ls=:dash, alpha=0.35, label=false)
        scatter!(p, t[swidx], x[swidx], m=:diamond, ms=5, c=:purple, label="committed switch")
    end
    p
end

p1 = plot_switch_panels_steps(t_switch, X_nomem_switch, swidx_nomem; title_str="No memory (α=1.0), multiplicative noise")
p2 = plot_switch_panels_steps(t_switch, X_mem_switch,   swidx_mem;   title_str="With memory (α=$(α_mem)), multiplicative noise")
plot(p1, p2, layout=(2,1), size=(950,650), legendposition=:topleft)

# -------------------------- dwell-time distributions --------------------------

# Simple log-spaced bins (integer edges)
function logbins_steps(x::Vector{<:Integer}; nb=10)
    xs = sort(x); lo = minimum(xs[xs .> 0]); hi = maximum(xs)
    unique(round.(Int, exp.(range(log(lo), log(hi), length=nb))))
end

bins1 = logbins_steps(dwell_nomem_steps)
bins2 = logbins_steps(dwell_mem_steps)

histogram(dwell_nomem_steps, bins=bins1, normalize=:pdf, alpha=0.6, label="No memory")
histogram!(dwell_mem_steps,   bins=bins2, normalize=:pdf, alpha=0.6, label="With memory")
xlabel!("Dwell length (steps)"); ylabel!("PDF"); title!("Committed dwell-time distributions (steps, multiplicative)")

# Survival curve (1 - ECDF)
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
xlabel!("dwell length (steps)"); ylabel!("Survival 1-ECDF"); title!("Committed dwell survival (steps, multiplicative)")

# -------------------------- optional: nicer shared-bin histograms + KDE --------------------------

"Freedman–Diaconis bin width (in steps) on pooled data; rounded to an int and floored at min_bw."
function fd_binwidth_steps(data::AbstractVector{<:Real}; min_bw::Int=20)
    bw = 2 * iqr(data) / length(data)^(1/3)
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
    lo = minimum(pooled[pooled .> 0])   # avoid log(0)
    hi = maximum(pooled)
    edges = unique(round.(Int, exp.(range(log(lo), log(hi), length=nbins))))
    sort(edges)
end

"Plot two overlaid histograms with shared edges, outlines, transparent fills, and median lines."
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

# 1) Linear shared bins
edges_lin = common_linear_edges_steps(dwell_nomem_steps, dwell_mem_steps; min_bw=25)
p_lin = plot_dwell_hists_steps(dwell_nomem_steps, dwell_mem_steps;
                               edges=edges_lin,
                               title_str="Committed dwell-time distributions (multiplicative)\n(steps, linear bins)")

# 2) Log-spaced shared bins
edges_log = common_log_edges_steps(dwell_nomem_steps, dwell_mem_steps; nbins=16)
p_log = plot_dwell_hists_steps(dwell_nomem_steps, dwell_mem_steps;
                               edges=edges_log,
                               title_str="Committed dwell-time distributions (multiplicative)\n(steps, log-spaced bins)")

plot(p_lin, p_log, layout=(1,2), size=(1100,420), legend=:topright)

# 3) KDE overlays on equal-width background hist
lo = min(minimum(dwell_nomem_steps), minimum(dwell_mem_steps))
hi = max(maximum(dwell_nomem_steps), maximum(dwell_mem_steps))
edges = collect(floor(Int, lo):100:ceil(Int, hi))   # equal-width 100-step bins

histogram(dwell_nomem_steps; bins=edges, normalize=:pdf,
          fillalpha=0.15, linealpha=0.3, label="No memory (hist)")
histogram!(dwell_mem_steps;   bins=edges, normalize=:pdf,
           fillalpha=0.15, linealpha=0.3, label="With memory (hist)")

kd_no = kde(Float64.(dwell_nomem_steps))
kd_me = kde(Float64.(dwell_mem_steps))
plot!(kd_no.x, kd_no.density; lw=3, label="No memory (KDE)")
plot!(kd_me.x, kd_me.density; lw=3, label="With memory (KDE)")
xlabel!("Dwell length (steps)"); ylabel!("PDF")
title!("Committed dwell-time distributions (multiplicative)")

vline!([median(dwell_nomem_steps)], ls=:dash, lc=:blue,   label="median no-mem")
vline!([median(dwell_mem_steps)],   ls=:dash, lc=:orange, label="median mem")

# 4) (Optional) Log-domain KDE for heavy tails: f_X(x) = f_Y(log x) / x
kd_no_log = kde(log.(Float64.(dwell_nomem_steps)))
kd_me_log = kde(log.(Float64.(dwell_mem_steps)))
x_no = exp.(kd_no_log.x);  pdf_no = kd_no_log.density ./ x_no
x_me = exp.(kd_me_log.x);  pdf_me = kd_me_log.density ./ x_me

histogram(dwell_nomem_steps; bins=edges, normalize=:pdf,
          fillalpha=0.1, linealpha=0.0, label="No memory (hist)")
histogram!(dwell_mem_steps;   bins=edges, normalize=:pdf,
           fillalpha=0.1, linealpha=0.0, label="With memory (hist)")
plot!(x_no, pdf_no; lw=3, label="No memory (log-KDE)")
plot!(x_me, pdf_me; lw=3, label="With memory (log-KDE)")
xlabel!("Dwell length (steps)"); ylabel!("PDF")
title!("Committed dwell-time distributions (multiplicative, log-KDE)")

vline!([median(dwell_nomem_steps)], ls=:dash, lc=:blue,   label="median no-mem")
vline!([median(dwell_mem_steps)],   ls=:dash, lc=:orange, label="median mem")



####

using FFTW, Statistics

"Return periodogram frequencies and power on zero mean data."
function periodogram_simple(x::AbstractVector{<:Real}, dt::Real)
    xz = x .- mean(x)
    n  = length(xz)
    Xf = abs2.(fft(xz)) ./ n
    half = fld(n, 2)
    freqs = (1:half) ./ (n*dt)   # one sided, drop DC
    power = Xf[2:half+1]
    return freqs, power
end

"High frequency power fraction above f0."
function fast_power_fraction(x::AbstractVector{<:Real}, dt::Real, f0::Real)
    f, P = periodogram_simple(x, dt)
    num = sum(P[f .> f0]); den = sum(P)
    return den == 0 ? 0.0 : num / den
end

"Law of total variance with bins in |x|. Returns (total, E[var|bin], var[E|bin], p, mu_k, v_k, edges)."
function variance_decomposition_abs(x::AbstractVector{<:Real}; nbins::Int=8)
    edges = range(0.0, stop=maximum(abs.(x)), length=nbins+1)
    binidx(v) = clamp(searchsortedfirst(edges, abs(v)) - 1, 1, nbins)
    b  = [binidx(xi) for xi in x]
    mu = mean(x)
    p    = [mean(b .== k) for k in 1:nbins]
    mu_k = [mean(x[b .== k]) for k in 1:nbins]
    v_k  = [var(x[b .== k])  for k in 1:nbins]
    Ev   = sum(p .* v_k)
    Vmu  = sum(p .* (mu_k .- mu).^2)
    return var(x), Ev, Vmu, p, mu_k, v_k, collect(edges)
end

# Stratonovich like midpoint evaluation for multiplicative term g(x)*xi(t)
# define g(x) and xi(t) in your main script
g(x) = max(0.0, 1 - 0.8x^2)

function F_mult_strat(t, x, xi, dt)
    x_pred = x .+ (x .- x.^3) .* dt .+ g.(x) .* xi(t) .* dt
    g_mid  = 0.5 .* (g.(x) .+ g.(x_pred))
    return (x .- x.^3) .+ g_mid .* xi(t)
end


# then compute steady segments and the same diagnostics on X_nomem_s and X_mem_s
# base intensity for multiplicative noise (short run)
σ_mult = 0.2

# piecewise-constant step noise for [tSpan[1], tSpan[2]] with step h
xi_mult = make_step_noise(tSpan[1], tSpan[2], h, σ_mult; seed=1)

# (optional) if you already have a Greek ξ_mult from elsewhere, you can just alias:
# xi_mult = ξ_mult

# Stratonovich-like RHS closure (midpoint for g)
rhs_strat(t, x) = F_mult_strat(t, x, xi_mult, h)

# now these will work:
t_vals_s, X_nomem_s = FDEsolver(rhs_strat, tSpan, y0, α_nomem; h=h)
_,         X_mem_s  = FDEsolver(rhs_strat, tSpan, y0, α_mem;   h=h)

# build steady segments the same way as before
steady_index = Int(floor(100.0/h))
X_nomem_s_steady = X_nomem_s[steady_index:end]
X_mem_s_steady   = X_mem_s[steady_index:end]

σ_mult_large  = 0.35          # bump if you get too few switches
T_switch      = 5000.0
tSpan_switch  = (0.0, T_switch)

xi_mult_large = make_step_noise(tSpan_switch[1], tSpan_switch[2], h, σ_mult_large; seed=2)
rhs_strat_big(t, x) = F_mult_strat(t, x, xi_mult_large, h)

t_sw, X_nomem_sw_s = FDEsolver(rhs_strat_big, [tSpan_switch...], y0, α_nomem; h=h)
_,    X_mem_sw_s   = FDEsolver(rhs_strat_big, [tSpan_switch...], y0, α_mem;   h=h)


# fast-content again (Stratonovich-like)
f0 = 0.2
phi_nomem_s = fast_power_fraction(X_nomem_s_steady, h, f0)
phi_mem_s   = fast_power_fraction(X_mem_s_steady,   h, f0)
println("fast power fraction (Strat-like) no mem = ", phi_nomem_s)
println("fast power fraction (Strat-like) mem    = ", phi_mem_s)

# variance decomposition (Strat-like)
V_nom_s, Ev_nom_s, Vmu_nom_s, _, _, _, _ = variance_decomposition_abs(X_nomem_s_steady; nbins=8)
V_mem_s, Ev_mem_s, Vmu_mem_s, _, _, _, _ = variance_decomposition_abs(X_mem_s_steady; nbins=8)
println("total var no mem (S) = ", V_nom_s, " | within = ", Ev_nom_s, " | between = ", Vmu_nom_s)
println("total var with mem (S) = ", V_mem_s, " | within = ", Ev_mem_s, " | between = ", Vmu_mem_s)
