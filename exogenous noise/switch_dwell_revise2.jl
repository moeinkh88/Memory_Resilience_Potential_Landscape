###############################
# Switching analysis (simple) — additive noise, transitions + residency
###############################
using FdeSolver, Random, Statistics, KernelDensity, StatsBase, Plots, Printf

# ---------- 0) Settings ----------
Random.seed!(123)

# dynamics
α_nomem = 1.0          # no memory
α_mem   = 0.7          # with memory
y0      = 0.0
h       = 0.01

# long run
T_switch     = 5000.0
σ_add_large  = 0.34
SEED_NOISE   = 2

# transition detection
x_core     = 0.9                  # threshold for being "in a well"
hold_time  = 0.20                 # must stay this long to confirm
hold_steps = max(1, Int(ceil(hold_time / h)))

# ---------- 1) Helpers ----------
drift(x) = x - x^3

"Piecewise-constant white-noise ξ(t) over [t0,t1] with step h and sd σ."
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

rhs_with_noise(ξ) = (t, x) -> drift.(x) .+ ξ(t)

"Label for core membership at a sample."
label_core_scalar(x; xc::Real=x_core) = x ≥ xc ? 1 : (x ≤ -xc ? -1 : 0)

"All committed entries: first index of any run of length ≥ hold_steps inside ±core."
function committed_entries(x::AbstractVector{<:Real}; xc::Real=x_core, hold_steps::Int=hold_steps)
    n = length(x)
    labs = similar(x, Int); @inbounds for i in 1:n labs[i] = label_core_scalar(x[i]; xc=xc) end
    idx, lab = Int[], Int[]
    i = 1
    while i ≤ n - hold_steps + 1
        ℓ = labs[i]
        if ℓ == 0
            i += 1; continue
        end
        ok = true
        @inbounds for j in i:i+hold_steps-1
            if labs[j] != ℓ; ok = false; break; end
        end
        if ok
            push!(idx, i); push!(lab, ℓ)
            # skip forward while staying in this core run
            k = i + hold_steps
            @inbounds while k ≤ n && labs[k] == ℓ; k += 1; end
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

"Shared integer log-bins for two step vectors."
function logbins_steps_two(a::Vector{<:Integer}, b::Vector{<:Integer}; nb::Int=30)
    xs = vcat(a, b); xs = xs[xs .> 0]
    edges = round.(Int, exp.(range(log(minimum(xs)), log(maximum(xs)), length=nb)))
    sort(unique(edges))
end

"Simple survival (1-ECDF)."
function survival_curve(x::AbstractVector{<:Real})
    xs = sort(x); n = length(xs)
    τ  = Float64.(xs)
    S  = 1 .- (1:n)./n
    return τ, S
end

# ---------- 2) Simulate (shared noise path) ----------
ξ = make_step_noise(0.0, T_switch, h, σ_add_large; seed=SEED_NOISE)
F = rhs_with_noise(ξ)

t_switch, X_nom = FDEsolver(F, [0.0, T_switch], y0, α_nomem; h=h)
_,        X_mem = FDEsolver(F, [0.0, T_switch], y0, α_mem;   h=h)

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


# (b) Residency histogram + KDE (steps, truncated at x >= 0)
if isempty(res_nom_time) || isempty(res_mem_time)
    @warn "No residency intervals to plot."
else
    lo = min(minimum(res_nom_time), minimum(res_mem_time))
    hi = max(maximum(res_nom_time), maximum(res_mem_time))

    # Equal-width bins; tweak bin_width if you want more/less granularity
    bin_width = max(1, round(Int, (hi - lo) / 40))   # ~40 bins
    edges = collect(floor(Int, lo):bin_width:ceil(Int, hi))

    p_res_kde = histogram(res_nom_time; bins=edges, normalize=:pdf,
                          fillalpha=0.35, linealpha=0.5, label="No memory (hist)")
    histogram!(p_res_kde, res_mem_time; bins=edges, normalize=:pdf,
               fillalpha=0.35, linealpha=0.5, label="With memory (hist)")

    # KDE overlays (truncate to x >= 0)
    kd_no = kde(Float64.(res_nom_time))
    kd_me = kde(Float64.(res_mem_time))
    mask_no = kd_no.x .>= 0
    mask_me = kd_me.x .>= 0
    plot!(p_res_kde, kd_no.x[mask_no], kd_no.density[mask_no]; lw=3, label="No memory (KDE)", color=:royalblue)
    plot!(p_res_kde, kd_me.x[mask_me], kd_me.density[mask_me]; lw=3, label="With memory (KDE)", color=:crimson)

    # Medians
    # vline!(p_res_kde, [median(res_nom_time)]; ls=:dash, lc=:blue,   label="median no-mem")
    # vline!(p_res_kde, [median(res_mem_time)]; ls=:dash, lc=:red,    label="median mem")
    vline!(p_res_kde, [mean(res_nom_time)]; ls=:dash, lc=:blue,   label="mean no-mem")
    vline!(p_res_kde, [mean(res_mem_time)]; ls=:dash, lc=:red,    label="mean mem")

    xlabel!(p_res_kde, "Residency time between transitions")
    ylabel!(p_res_kde, "PDF")
    title!(p_res_kde, "Residency distributions (hist + KDE)")
    display(p_res_kde)
end




# (c) residency survival (time)
τ1, S1 = survival_curve(res_nom_time)
τ2, S2 = survival_curve(res_mem_time)
p_surv = plot(τ1, S1; lw=2, label="No memory")
plot!(p_surv, τ2, S2; lw=2, label="With memory")
xlabel!(p_surv, "Residency between transitions (time)")
ylabel!(p_surv, "Survival 1-ECDF")
title!(p_surv, "Residency survival")
display(p_surv)
