###############################
# Switching analysis (clean)
###############################
using FdeSolver, Random, Statistics, StatsBase, KernelDensity, Plots

# ---------- 0) Settings ----------
Random.seed!(123)                 # reproducibility for anything else

# dynamics
α_nomem = 1.0                     # no-memory
α_mem   = 0.7                     # fractional order (memory)
y0      = 1.0                     # start near +well
h       = 0.01                    # solver step (time units)

# long run for switching
T_switch     = 6000.0             # total time for switching run
σ_add_large  = 0.35               # stronger additive noise to induce switches
SEED_NOISE   = 2                  # seed for the shared noise path

# switching detection (choose in time units; converted to steps)
x_core    = 0.9                   # core threshold for ± wells
hold_time = 0.20                  # must remain in core ≥ this time to confirm entry
hold_steps = max(1, Int(ceil(hold_time / h)))

# ---------- 1) Helpers ----------
"Piecewise-constant white-noise approximation ξ(t) over [t0,t1] with step h and sd σ."
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=1)
    N   = Int(ceil((t1 - t0)/h))
    rng = MersenneTwister(seed)
    vals = (σ/√h) .* randn(rng, N)             # so that ∫ ξ dt ≈ σ dW
    t0f, hf = float(t0), float(h)
    return (t::Real)-> begin
        idx = clamp(Int(floor((t - t0f)/hf)) + 1, 1, N)
        return vals[idx]
    end
end

drift(x) = x - x^3                            # double-well drift
rhs_with_noise(ξ) = (t, x) -> drift.(x) .+ ξ(t)

label_core_scalar(x; xc=x_core) = x ≥ xc ? 1 : (x ≤ -xc ? -1 : 0)

"Find committed entries into ±cores with a 'hold_steps' requirement."
function core_entries_steps(x::AbstractVector{<:Real}; xc::Real=x_core, hold_steps::Int=hold_steps)
    n = length(x)
    labs = Vector{Int}(undef, n)
    @inbounds @simd for i in 1:n
        labs[i] = label_core_scalar(x[i]; xc=xc)
    end
    idx, lab_out = Int[], Int[]
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
            push!(idx, i); push!(lab_out, ℓ)
            # skip forward while staying in the same core
            k = i + hold_steps
            @inbounds while k ≤ n && labs[k] == ℓ
                k += 1
            end
            i = k
        else
            i += 1
        end
    end
    return idx, lab_out
end

"Dwell lengths in *steps* between committed entries into opposite cores."
# accept vectors OR matrices; internally convert to a column vector
function dwell_steps_committed(x::AbstractArray{<:Real}; xc::Real=x_core, hold_steps::Int=hold_steps)
    xv = x isa AbstractVector ? x : vec(x)
    idx_all, labs_all = core_entries_steps(xv; xc=xc, hold_steps=hold_steps)
    dw_steps, switch_idx = Int[], Int[]
    for k in 2:length(idx_all)
        if labs_all[k] != labs_all[k-1]
            push!(dw_steps, idx_all[k] - idx_all[k-1])
            push!(switch_idx, idx_all[k])
        end
    end
    return dw_steps, switch_idx
end

"Make (shared) log bins in integer *steps* for two dwell vectors."
function logbins_steps_two(a::Vector{<:Integer}, b::Vector{<:Integer}; nb::Int=12)
    xs = vcat(a, b)
    xs = xs[xs .> 0]
    lo, hi = minimum(xs), maximum(xs)
    edges = round.(Int, exp.(range(log(lo), log(hi), length=nb)))
    return sort(unique(edges))
end

# ---------- 2) Simulate (shared noise path) ----------
ξL = make_step_noise(0.0, T_switch, h, σ_add_large; seed=SEED_NOISE)
F  = rhs_with_noise(ξL)

t_switch, X_nomem_switch = FDEsolver(F, [0.0, T_switch], y0, α_nomem; h=h)
_,        X_mem_switch   = FDEsolver(F, [0.0, T_switch], y0, α_mem;   h=h)

# ---------- 3) Switching stats ----------
x_nom = vec(X_nomem_switch)   # works for N×1, 1×N, Adjoint, etc.
x_mem = vec(X_mem_switch)

dwell_nomem_steps, swidx_nomem = dwell_steps_committed(x_nom; xc=x_core, hold_steps=hold_steps)
dwell_mem_steps,   swidx_mem   = dwell_steps_committed(x_mem; xc=x_core, hold_steps=hold_steps)

dwell_nomem_time = h .* dwell_nomem_steps
dwell_mem_time   = h .* dwell_mem_steps

println("-- committed switching summary --")
println("hold_steps = $hold_steps  (≈ $(hold_steps*h) time units),  x_core = $x_core")
println("switches (no mem)   = ", length(swidx_nomem), "   mean dwell (time) = ", mean(dwell_nomem_steps))
println("switches (with mem) = ", length(swidx_mem),   "   mean dwell (time) = ", mean(dwell_mem_steps))

# ---------- 4) Plots ----------
# (a) time series with cores and committed switches
function plot_switch_panel(t, x, swidx; title_str="")
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x, lw=0.8, color=:gray70, label=false)
    hline!(p, [0.0], ls=:dash, c=:black, label="unstable x=0")
    hline!(p, [ x_core], ls=:dot,  c=:blue,  label="core +x")
    hline!(p, [-x_core], ls=:dot,  c=:red,   label="core -x")
    if !isempty(swidx)
        vline!(p, t[swidx], c=:gray, ls=:dash, alpha=0.35, label=false)
        scatter!(p, t[swidx], x[swidx], m=:diamond, ms=5, c=:purple, label="committed switch")
    end
    return p
end

p_ts1 = plot_switch_panel(t_switch, X_nomem_switch, swidx_nomem; title_str="No memory (α=1.0)")
p_ts2 = plot_switch_panel(t_switch, X_mem_switch,   swidx_mem;   title_str="With memory (α=$(α_mem))")
display(plot(p_ts1, p_ts2, layout=(2,1), size=(960,650), legendposition=:topleft))

# (b) dwell-time distributions in *steps* (shared log bins)
edges_steps = logbins_steps_two(dwell_nomem_steps, dwell_mem_steps; nb=30)
p_hist = histogram(dwell_nomem_steps; bins=edges_steps, normalize=:pdf, alpha=0.6, label="No memory")
histogram!(p_hist, dwell_mem_steps;   bins=edges_steps, normalize=:pdf, alpha=0.6, label="With memory")
xlabel!(p_hist, "Dwell length (steps)"); ylabel!(p_hist, "PDF")
title!(p_hist, "Committed dwell-time distributions (steps)")
display(p_hist)

# (c) survival curves (1-ECDF) in *time units*
function survival_curve(x::AbstractVector{<:Real})
    xs = sort(x); n = length(xs)
    τ  = Float64.(xs)
    S  = 1 .- (1:n)./n
    return τ, S
end
τ1, S1 = survival_curve(dwell_nomem_time)
τ2, S2 = survival_curve(dwell_mem_time)
p_surv = plot(τ1, S1; lw=2, label="No memory")
plot!(p_surv, τ2, S2; lw=2, label="With memory")
xlabel!(p_surv, "Dwell length (time)"); ylabel!(p_surv, "Survival 1-ECDF")
title!(p_surv, "Committed dwell survival (time)")
display(p_surv)

# (d) histogram background + KDE (steps, truncated at x >= 0)
lo = min(minimum(dwell_nomem_steps), minimum(dwell_mem_steps))
hi = max(maximum(dwell_nomem_steps), maximum(dwell_mem_steps))
edges = collect(floor(Int, lo):100:ceil(Int, hi))   # equal-width 100-step bins

histogram(dwell_nomem_steps; bins=edges, normalize=:pdf,
          fillalpha=0.35, linealpha=0.5, label="No memory (hist)")
histogram!(dwell_mem_steps;   bins=edges, normalize=:pdf,
           fillalpha=0.35, linealpha=0.5, label="With memory (hist)")

kd_no = kde(Float64.(dwell_nomem_steps))
kd_me = kde(Float64.(dwell_mem_steps))

mask_no = kd_no.x .>= 0
mask_me = kd_me.x .>= 0

plot!(kd_no.x[mask_no], kd_no.density[mask_no]; lw=3, label="No memory (KDE)", color=:royalblue)
plot!(kd_me.x[mask_me], kd_me.density[mask_me]; lw=3, label="With memory (KDE)", color=:crimson)

xlabel!("Dwell length (steps)"); ylabel!("PDF")
title!("Committed dwell-time distributions")

vline!([median(dwell_nomem_steps)], ls=:dash, lc=:blue,   label="median no-mem")
vline!([median(dwell_mem_steps)],   ls=:dash, lc=:red, label="median mem")
