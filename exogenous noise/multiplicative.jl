using FdeSolver, Random, Statistics, StatsBase, KernelDensity, Plots
using Printf

# ========================== setup ==========================

Random.seed!(123)

# model and time
const α_nomem = 1.0        # no memory
const α_mem   = 0.8        # memory case
const y0      = 1.0
const T       = 1000.0
const h       = 0.01
const tSpan   = (0.0, T)

# multiplicative noise base intensity for short run
const σ_mult  = 0.2

# state dependent noise amplitude, lower near the wells
noise_amp(x) = max(0.0, 1.0 - 0.8*x^2)

"""
step noise ξ(t) on [t0, t1] with step h and standard deviation σ
piecewise constant white noise approximation
"""
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=123)
    N   = Int(ceil((t1 - t0) / h))
    rng = MersenneTwister(seed)
    vals = (σ / √h) .* randn(rng, N)
    t0f, hf = float(t0), float(h)
    return (t::Real) -> begin
        idx = clamp(Int(floor((t - t0f) / hf)) + 1, 1, N)
        vals[idx]
    end
end

# one fixed noise path used for both α so the driving noise is comparable
const ξ_mult = make_step_noise(tSpan[1], tSpan[2], h, σ_mult; seed=1)

# system with multiplicative noise: dx/dt = x - x^3 + noise_amp(x)*ξ(t)
F_multiplicative(t, x) = (x .- x.^3) .+ noise_amp.(x) .* ξ_mult(t)

# ===================== helper functions =====================

"""
biased autocorrelation up to maxlag
returns a vector acf with acf[1] = 1
"""
function autocorr_biased(x::AbstractVector{<:Real}, maxlag::Int)
    n = length(x)
    μ = mean(x)
    σ2 = var(x)
    acf = Vector{Float64}(undef, maxlag + 1)
    acf[1] = 1.0
    @inbounds for τ in 1:maxlag
        s = 0.0
        for i in 1:(n - τ)
            s += (x[i] - μ) * (x[i + τ] - μ)
        end
        acf[τ + 1] = s / ((n - τ) * σ2)
    end
    return acf
end

"""
integrated autocorrelation time over a window W_time
"""
function tau_int(acf::AbstractVector{<:Real}, h::Real, W_time::Real)
    K = min(length(acf) - 1, Int(round(W_time / h)))
    h * (1 + 2 * sum(@view acf[2:K+1]))
end

"""
block average by factor m
drops any leftover samples that do not fill a full block
"""
function block_average(x::AbstractVector{<:Real}, m::Int)
    nblocks = length(x) ÷ m
    @inbounds [mean(@view x[(i - 1) * m + 1 : i * m]) for i in 1:nblocks]
end

# ========================= simulate =========================

t_vals, X_nomem = FDEsolver(F_multiplicative, collect(tSpan), y0, α_nomem; h=h)
_,      X_mem   = FDEsolver(F_multiplicative, collect(tSpan), y0, α_mem;   h=h)

# trim initial transient
t_trim = 100.0
i0 = max(1, Int(floor(t_trim / h)) + 1)
@views Xn = X_nomem[i0:end]
@views Xm = X_mem[i0:end]
@views tn = t_vals[i0:end]

# ======================= quick statistics =======================

var_nomem = var(Xn)
var_mem   = var(Xm)
@printf "Sample variance no memory   = %.6f\n" var_nomem
@printf "Sample variance with memory = %.6f\n" var_mem

# ===================== autocorrelation analysis =====================

# long acf for semilog plot
maxlag_steps_long = 20000                  # 200 time units
acf_nomem_long = autocorr_biased(Xn, maxlag_steps_long)
acf_mem_long   = autocorr_biased(Xm, maxlag_steps_long)
lags_time_long = (0:maxlag_steps_long) .* h

p_acf = plot(lags_time_long[2:end], acf_nomem_long[2:end],
             xscale=:log10, lw=2, label="no memory")
plot!(p_acf, lags_time_long[2:end], acf_mem_long[2:end],
      xscale=:log10, lw=2, label="with memory")
xlabel!(p_acf, "lag time")
ylabel!(p_acf, "ACF")
title!(p_acf, "ACF semilog under multiplicative noise")
display(p_acf)

# integrated acf over finite windows
for W in (20.0, 50.0)
    @printf "τ_int(%.0f) no memory   = %.6f\n" W  tau_int(acf_nomem_long, h, W)
    @printf "τ_int(%.0f) with memory = %.6f\n" W  tau_int(acf_mem_long,   h, W)
end

# ====================== short segment plot ======================

t_plot_start = 200.0
t_plot_end   = 250.0
i_start = max(1, Int(floor(t_plot_start / h)))
i_end   = min(length(t_vals), Int(floor(t_plot_end / h)))
p_seg = plot(t_vals[i_start:i_end], X_mem[i_start:i_end],
             lw=2, label="with memory α=$(α_mem)", color=:indianred3)
plot!(p_seg, t_vals[i_start:i_end], X_nomem[i_start:i_end],
      lw=2, label="no memory α=1", color=:dodgerblue1)
xlabel!(p_seg, "time")
ylabel!(p_seg, "x(t)")
title!(p_seg, "Trajectories under multiplicative noise")
display(p_seg)

# =================== coarse graining study ===================

m_list = [1, 5, 10, 20, 50, 100, 150, 200]
W_cg   = 50.0

tau_no = Float64[]; tau_me = Float64[]
var_no = Float64[]; var_me = Float64[]

for m in m_list
    Xn_b = block_average(Xn, m)
    Xm_b = block_average(Xm, m)
    h_eff = h * m
    maxlag_eff = Int(round(W_cg / h_eff))
    acf_n_b = autocorr_biased(Xn_b, maxlag_eff)
    acf_m_b = autocorr_biased(Xm_b, maxlag_eff)
    push!(tau_no, tau_int(acf_n_b, h_eff, W_cg))
    push!(tau_me, tau_int(acf_m_b, h_eff, W_cg))
    push!(var_no, var(Xn_b))
    push!(var_me, var(Xm_b))
end

p_tau = plot(m_list, tau_no, lw=2, marker=:o, label="no memory")
plot!(p_tau, m_list, tau_me, lw=2, marker=:o, label="with memory")
xlabel!(p_tau, "block size m")
ylabel!(p_tau, "τ_int over W=$(W_cg)")
title!(p_tau, "Persistence vs coarse graining")
display(p_tau)

p_var = plot(m_list, var_no, lw=2, marker=:o, label="no memory")
plot!(p_var, m_list, var_me, lw=2, marker=:o, label="with memory")
xlabel!(p_var, "block size m")
ylabel!(p_var, "variance")
title!(p_var, "Variance after coarse graining")
display(p_var)

# -------------------------- switching experiment --------------------------
###############################
# Multiplicative noise — transitions + residency (time, log-binned)
###############################

T_switch     = 5000.0
σ_mult_large = 0.5      # multiplicative noise intensity
SEED_NOISE   = 2

# transition detection
x_core     = 0.9                   # "in well" threshold
hold_time  = 0.20
hold_steps = max(1, Int(ceil(hold_time / h)))

# ---------- 1) Helpers ----------
drift(x) = x - x^3

# state-dependent noise amplitude (example; tweak if desired)
noise_amp(x) = max(0.0, 1.0 - 0.8*x^2)   # strong near 0, weaker near ±1

"Piecewise-constant white-noise ξ(t) over [t0,t1] with step h and sd σ."
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=1)
    N   = Int(ceil((t1 - t0)/h))
    rng = MersenneTwister(seed)
    vals = (σ/√h) .* randn(rng, N)        # so that ∫ ξ dt ≈ σ dW
    t0f, hf = float(t0), float(h)
    return (t::Real)-> begin
        idx = clamp(Int(floor((t - t0f)/hf)) + 1, 1, N)
        return vals[idx]
    end
end

"Multiplicative SDE RHS using a shared noise path."
function make_F_multiplicative(ξ)
    return (t, x) -> (x .- x.^3) .+ noise_amp.(x) .* ξ(t)
end

"Label for core membership at a sample."
label_core_scalar(x; xc::Real=x_core) = x ≥ xc ? 1 : (x ≤ -xc ? -1 : 0)

"All committed entries: first index of any run of length ≥ hold_steps inside ±core."
function committed_entries(x::AbstractVector{<:Real}; xc::Real=x_core, hold_steps::Int=hold_steps)
    n = length(x)
    labs = similar(x, Int)
    @inbounds for i in 1:n
        labs[i] = label_core_scalar(x[i]; xc=xc)
    end
    idx, lab = Int[], Int[]
    i = 1
    while i ≤ n - hold_steps + 1
        ℓ = labs[i]
        if ℓ == 0
            i += 1
            continue
        end
        ok = true
        @inbounds for j in i:i+hold_steps-1
            if labs[j] != ℓ; ok = false; break; end
        end
        if ok
            push!(idx, i); push!(lab, ℓ)
            # skip forward while staying in this core run
            k = i + hold_steps
            @inbounds while k ≤ n && labs[k] == ℓ
                k += 1
            end
            i = k
        else
            i += 1
        end
    end
    return idx, lab
end

"Keep only alternating entries (−→+ or +→−); drop same-side re-entries."
function alternating_transitions(idx::Vector{Int}, lab::Vector{Int})
    keep_idx = Int[]; keep_lab = Int[]
    for k in eachindex(idx)
        if isempty(keep_lab) || lab[k] != keep_lab[end]
            push!(keep_idx, idx[k]); push!(keep_lab, lab[k])
        end
    end
    return keep_idx, keep_lab
end

"Residency times (steps) = time between consecutive alternating transition starts."
residency_from_transitions(idx_alt::Vector{Int}) = length(idx_alt) ≥ 2 ? diff(idx_alt) : Int[]

"Simple survival (1-ECDF)."
function survival_curve(x::AbstractVector{<:Real})
    xs = sort(x); n = length(xs)
    τ  = Float64.(xs)
    S  = 1 .- (1:n)./n
    return τ, S
end

# ---------- 2) Simulate (shared noise path) ----------
ξ_mult1 = make_step_noise(0.0, T_switch, h, σ_mult_large; seed=SEED_NOISE)
F_mul  = make_F_multiplicative(ξ_mult1)

t_switch, X_nom = FDEsolver(F_mul, [0.0, T_switch], y0, α_nomem; h=h)
_,        X_mem = FDEsolver(F_mul, [0.0, T_switch], y0, α_mem;   h=h)

x_nom = vec(X_nom)
x_mem = vec(X_mem)

# ---------- 3) Transitions & residency ----------
# (a) committed entries
idx_nom_all, lab_nom_all = committed_entries(x_nom; xc=x_core, hold_steps=hold_steps)
idx_mem_all, lab_mem_all = committed_entries(x_mem; xc=x_core, hold_steps=hold_steps)

# (b) only true shifts (alternating)
idx_nom, lab_nom = alternating_transitions(idx_nom_all, lab_nom_all)
idx_mem, lab_mem = alternating_transitions(idx_mem_all, lab_mem_all)

# (c) residency (between transition starts)
res_nom_steps = residency_from_transitions(idx_nom)
res_mem_steps = residency_from_transitions(idx_mem)
res_nom_time  = h .* res_nom_steps
res_mem_time  = h .* res_mem_steps

# ---------- 3a) Text summary ----------
@printf "\n-- thresholds --\nx_core = %.3f, hold_steps = %d (≈ %.3f time)\n" x_core hold_steps hold_steps*h

println("\n-- transitions (alternating committed entries) --")
println("transitions: no-mem = $(length(idx_nom)), mem = $(length(idx_mem))")

println("\n-- residency (between transitions) --")
println("count: no-mem = $(length(res_nom_time)), mem = $(length(res_mem_time))   (should be transitions-1)")
println("median time:  no-mem = $(isempty(res_nom_time) ? NaN : median(res_nom_time)),  mem = $(isempty(res_mem_time) ? NaN : median(res_mem_time))")
@printf "mean time:     no-mem = %.3f, mem = %.3f\n" (isempty(res_nom_time) ? NaN : mean(res_nom_time)) (isempty(res_mem_time) ? NaN : mean(res_mem_time))

# ---------- 4) Plots ----------
# (a) time series with transitions marked
function plot_panel(t, x, idx_tr; title_str="")
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x; lw=0.8, color=:gray70, label="trajectory")
    hline!(p, [0.0];  ls=:dash, c=:black, label="x=0")
    hline!(p, [ x_core]; ls=:dot,  c=:blue,  label="+core")
    hline!(p, [-x_core]; ls=:dot,  c=:red,   label="-core")
    if !isempty(idx_tr)
        scatter!(p, t[idx_tr], x[idx_tr]; m=:diamond, ms=5, c=:purple, label="transition start")
    end
    return p
end

p1 = plot_panel(t_switch, X_nom, idx_nom; title_str="No memory (α=1.0)")
p2 = plot_panel(t_switch, X_mem, idx_mem; title_str="With memory (α=$(α_mem))")
display(plot(p1, p2, layout=(2,1), size=(960,650), legendposition=:topleft))

# (b) Residency histogram + KDE with LOG bins (time)
if isempty(res_nom_time) || isempty(res_mem_time)
    @warn "No residency intervals to plot."
else
    nb = 30
    vals = vcat(res_nom_time, res_mem_time)
    vals = vals[vals .> 0]                    # log bins need positive values
    lo, hi = minimum(vals), maximum(vals)
    edges = collect(exp.(range(log(lo), log(hi), length=nb)))

    p_res_kde = histogram(res_nom_time; bins=edges, normalize=:pdf,
                          fillalpha=0.35, linealpha=0.5, label="No memory (hist)")
    histogram!(p_res_kde, res_mem_time; bins=edges, normalize=:pdf,
               fillalpha=0.35, linealpha=0.5, label="With memory (hist)")

    # KDE overlays (truncate to x >= first edge)
    kd_no = kde(res_nom_time)
    kd_me = kde(res_mem_time)
    mask_no = kd_no.x .>= edges[1]
    mask_me = kd_me.x .>= edges[1]
    plot!(p_res_kde, kd_no.x[mask_no], kd_no.density[mask_no]; lw=3, label="No memory (KDE)")
    plot!(p_res_kde, kd_me.x[mask_me], kd_me.density[mask_me]; lw=3, label="With memory (KDE)")

    vline!(p_res_kde, [median(res_nom_time)]; ls=:dash, label="median no-mem")
    vline!(p_res_kde, [median(res_mem_time)]; ls=:dash, label="median mem")

    xlabel!(p_res_kde, "Residency between transitions (time)")
    ylabel!(p_res_kde, "PDF")
    title!(p_res_kde, "Residency distributions (log bins, time)")
    display(p_res_kde)
end

# (c) Residency survival (time)
τ1, S1 = survival_curve(res_nom_time)
τ2, S2 = survival_curve(res_mem_time)
p_surv = plot(τ1, S1; lw=2, label="No memory")
plot!(p_surv, τ2, S2; lw=2, label="With memory")
xlabel!(p_surv, "Residency between transitions (time)")
ylabel!(p_surv, "Survival 1-ECDF")
title!(p_surv, "Residency survival")
display(p_surv)
