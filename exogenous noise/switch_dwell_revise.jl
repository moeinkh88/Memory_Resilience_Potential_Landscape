###############################
# Switching analysis (clean) — inter-switch + residence
###############################
using FdeSolver, Random, Statistics, StatsBase, KernelDensity, Plots, Printf

# ---------- 0) Settings ----------
Random.seed!(123)

# dynamics
α_nomem = 1.0
α_mem   = 0.7
y0      = 0.0
h       = 0.01

# long run for switching
T_switch     = 6000.0
σ_add_large  = 0.35
SEED_NOISE   = 2

# switching detection (choose in time units; converted to steps)
x_core     = 0.9                        # committed-entry threshold
hold_time  = 0.20
hold_steps = max(1, Int(ceil(hold_time / h)))
x_gate     = max(0.0, x_core - 0.001)    # <— hysteresis gate for *exit* (residence)

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

drift(x) = x - x^3
rhs_with_noise(ξ) = (t, x) -> drift.(x) .+ ξ(t)

label_core_scalar(x; xc=x_core) = x ≥ xc ? 1 : (x ≤ -xc ? -1 : 0)

"Committed entries into ±cores with a 'hold_steps' requirement."
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

# ---------- New: two metrics ----------
"Inter-switch intervals (steps): between committed entries into opposite cores."
function interswitch_steps(x::AbstractArray{<:Real}; xc::Real=x_core, hold_steps::Int=hold_steps)
    xv = x isa AbstractVector ? x : vec(x)
    idx, labs = core_entries_steps(xv; xc=xc, hold_steps=hold_steps)
    Δ, swidx = Int[], Int[]
    for k in 2:length(idx)
        if labs[k] != labs[k-1]
            push!(Δ, idx[k] - idx[k-1])
            push!(swidx, idx[k])
        end
    end
    return Δ, swidx
end

"Residence (dwell) with committed exit: entry (committed) → first exit past x_gate held for exit_hold_steps."
function residence_steps_committed(x::AbstractArray{<:Real};
        xc::Real=x_core, hold_steps::Int=hold_steps,
        xgate::Real=x_gate, exit_hold_steps::Int=hold_steps)

    xv = x isa AbstractVector ? x : vec(x)
    n = length(xv)
    entries, labs = core_entries_steps(xv; xc=xc, hold_steps=hold_steps)

    is_outside = (val, ℓ)-> (ℓ== 1 && val <  xgate) || (ℓ==-1 && val > -xgate)

    res = Int[]; ent_out = Int[]; lab_out = Int[]
    for (eidx, ℓ) in zip(entries, labs)
        j = eidx + 1
        consec = 0
        while j ≤ n
            if is_outside(xv[j], ℓ)
                consec += 1
                if consec ≥ exit_hold_steps
                    len_steps = (j - exit_hold_steps) - eidx + 1    # last inside index = j-exit_hold_steps
                    push!(res, max(len_steps, 1))
                    push!(ent_out, eidx); push!(lab_out, ℓ)
                    break
                end
            else
                consec = 0
            end
            j += 1
        end
        if j > n
            # right-censored at the end
            push!(res, n - eidx + 1)
            push!(ent_out, eidx); push!(lab_out, ℓ)
        end
    end
    return res, ent_out, lab_out
end



"Shared integer log-bins for two step vectors."
function logbins_steps_two(a::Vector{<:Integer}, b::Vector{<:Integer}; nb::Int=12)
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
ξL = make_step_noise(0.0, T_switch, h, σ_add_large; seed=SEED_NOISE)
F  = rhs_with_noise(ξL)

t_switch, X_nomem_switch = FDEsolver(F, [0.0, T_switch], y0, α_nomem; h=h)
_,        X_mem_switch   = FDEsolver(F, [0.0, T_switch], y0, α_mem;   h=h)

x_nom = vec(X_nomem_switch)
x_mem = vec(X_mem_switch)

# ---------- 3) Metrics ----------
# (A) Inter-switch intervals (unchanged)
is_nom_steps, swidx_nom = interswitch_steps(x_nom; xc=x_core, hold_steps=hold_steps)
is_mem_steps, swidx_mem = interswitch_steps(x_mem; xc=x_core, hold_steps=hold_steps)
is_nom_time = h .* is_nom_steps
is_mem_time = h .* is_mem_steps

# (B) Residence inside a core (committed exit; include initial if it qualifies)
res_nom_steps, ent_nom, labs_nom = residence_steps_committed(x_nom; xc=x_core, hold_steps=hold_steps,
                                                             xgate=x_gate, exit_hold_steps=hold_steps)
res_mem_steps, ent_mem, labs_mem = residence_steps_committed(x_mem; xc=x_core, hold_steps=hold_steps,
                                                             xgate=x_gate, exit_hold_steps=hold_steps)
res_nom_time = h .* res_nom_steps
res_mem_time = h .* res_mem_steps

# sanity: same arrays used below
println("\n-- thresholds --")
println("x_core = $x_core, hold_steps = $hold_steps (≈ $(hold_steps*h) time), x_gate = $x_gate")
println("residence exit-hold = $hold_steps steps")

println("\n-- inter-switch intervals (committed entries to opposite core) --")
println("count: no-mem = $(length(is_nom_steps)), mem = $(length(is_mem_steps))")
println("median (steps): no-mem = $(median(is_nom_steps)), mem = $(median(is_mem_steps))")
@printf "10%% trimmed mean (time): no-mem = %.3f, mem = %.3f\n" mean(is_nom_time) mean(is_mem_time)

println("\n-- residence inside core (committed entry → committed exit past gate) --")
println("count: no-mem = $(length(res_nom_steps)), mem = $(length(res_mem_steps))")
println("median (steps): no-mem = $(median(res_nom_steps)), mem = $(median(res_mem_steps))")
@printf "mean (time):  no-mem = %.3f, mem = %.3f\n" mean(res_nom_time) mean(res_mem_time)


# ---------- 4) Plots ----------
# (a) time series with cores and committed switches (diamonds)
function plot_switch_panel(t, x, swidx; title_str="")
    p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
    plot!(p, t, x, lw=0.8, color=:gray70, label=false)
    hline!(p, [0.0],  ls=:dash, c=:black, label="unstable x=0")
    hline!(p, [ x_core], ls=:dot,  c=:blue,  label="core +x")
    hline!(p, [-x_core], ls=:dot,  c=:red,   label="core -x")
    if !isempty(swidx)
        vline!(p, t[swidx], c=:gray, ls=:dash, alpha=0.35, label=false)
        scatter!(p, t[swidx], x[swidx], m=:diamond, ms=5, c=:purple, label="committed switch")
    end
    return p
end
p_ts1 = plot_switch_panel(t_switch, X_nomem_switch, swidx_nom; title_str="No memory (α=1.0)")
p_ts2 = plot_switch_panel(t_switch, X_mem_switch,   swidx_mem; title_str="With memory (α=$(α_mem))")
display(plot(p_ts1, p_ts2, layout=(2,1), size=(960,650), legendposition=:topleft))


# function residence_spans_committed(x; xc=x_core, hold_steps=hold_steps, xgate=x_gate, exit_hold_steps=hold_steps)
#     res, entries, labs = residence_steps_committed(x; xc=xc, hold_steps=hold_steps, xgate=xgate, exit_hold_steps=exit_hold_steps)
#     spans = Tuple{Int,Int,Int}[]
#     for (eidx, ℓ, L) in zip(entries, labs, res)
#         push!(spans, (eidx, eidx + L - 1, ℓ))
#     end
#     return spans
# end


# --- helper: compute residence spans as index ranges with labels (-1/+1)
# function residence_spans(x::AbstractVector{<:Real};
#                          xc::Real=x_core, hold_steps::Int=hold_steps, xgate::Real=x_gate)
#     xv = vec(x); n = length(xv)
#     entries, labs = core_entries_steps(xv; xc=xc, hold_steps=hold_steps)
#     out_of_gate = (val, ℓ)-> (ℓ== 1 && val <  xgate) || (ℓ==-1 && val > -xgate)

#     spans = Tuple{Int,Int,Int}[]  # (i_start, i_end, label)
#     for (eidx, ℓ) in zip(entries, labs)
#         j = eidx + 1
#         @inbounds while j ≤ n && !out_of_gate(xv[j], ℓ)
#             j += 1
#         end
#         push!(spans, (eidx, j-1, ℓ))  # last index still inside gate
#     end
#     return spans
# end

# # --- new colored panel
# function plot_switch_panel_colored(t, x, swidx; title_str="", palette::Symbol=:nomem)
#     # colors per palette
#     if palette == :nomem
#         c_neg = RGBA(0.55, 0.75, 1.00, 0.9)   # light blue  (− well)
#         c_pos = RGB(0.00, 0.20, 0.60)         # dark  blue  (+ well)
#     else
#         c_neg = RGBA(1.00, 0.65, 0.65, 0.9)   # light red   (− well)
#         c_pos = RGB(0.70, 0.05, 0.10)         # dark  red   (+ well)
#     end

#     p = plot(xlabel="time", ylabel="x(t)", title=title_str, legend=:topright)
#     # base trajectory in gray (covers inter-switch intervals too)
#     plot!(p, t, x; lw=0.8, color=:gray70, label="trajectory")

#     # core lines
#     hline!(p, [0.0];  ls=:dash, c=:black, label="unstable x=0")
#     hline!(p, [ x_core]; ls=:dot,  c=:blue,  label="core +x")
#     hline!(p, [-x_core]; ls=:dot,  c=:red,   label="core -x")

#     # overlay residence segments in color
#     for (i1, i2, ℓ) in residence_spans(x; xc=x_core, hold_steps=hold_steps, xgate=x_gate)
#         col = (ℓ == 1) ? c_pos : c_neg
#         plot!(p, t[i1:i2], x[i1:i2]; lw=2.2, color=col, label=false)
#     end

#     # add legend swatches (without drawing extra data)
#     plot!(p, [NaN], [NaN]; lw=3, color=c_pos, label="+ well (residence)")
#     plot!(p, [NaN], [NaN]; lw=3, color=c_neg, label="− well (residence)")

#     # committed switches
#     if !isempty(swidx)
#         vline!(p, t[swidx]; c=:gray, ls=:dash, alpha=0.35, label=false)
#         scatter!(p, t[swidx], x[swidx]; m=:diamond, ms=5, c=:purple, label="committed switch")
#     end
#     return p
# end

# # --- use it
# p_ts1 = plot_switch_panel_colored(t_switch, vec(X_nomem_switch), swidx_nom;
#                                   title_str="No memory (α=1.0)", palette=:nomem)
# p_ts2 = plot_switch_panel_colored(t_switch, vec(X_mem_switch),   swidx_mem;
#                                   title_str="With memory (α=$(α_mem))", palette=:mem)
# display(plot(p_ts1, p_ts2, layout=(2,1), size=(960,650), legendposition=:outertopright))


# Inter-switch PDF
edges_is = logbins_steps_two(is_nom_steps, is_mem_steps; nb=30)
p_is = histogram(is_nom_steps; bins=edges_is, normalize=:pdf, alpha=0.6, label="No memory")
histogram!(p_is, is_mem_steps; bins=edges_is, normalize=:pdf, alpha=0.6, label="With memory")
xlabel!(p_is, "Inter-switch interval (steps)"); ylabel!(p_is, "PDF")
title!(p_is, "Committed inter-switch distributions (steps)")
display(p_is)

# Inter-switch survival
τ1, S1 = survival_curve(is_nom_time)
τ2, S2 = survival_curve(is_mem_time)
p_surv_is = plot(τ1, S1; lw=2, label="No memory")
plot!(p_surv_is, τ2, S2; lw=2, label="With memory")
xlabel!(p_surv_is, "Inter-switch interval (time)"); ylabel!(p_surv_is, "Survival 1-ECDF")
title!(p_surv_is, "Inter-switch survival (time)")
display(p_surv_is)

# Residence PDF (now using committed arrays & shared log bins)
edges_res = logbins_steps_two(res_nom_steps, res_mem_steps; nb=30)
p_res = histogram(res_nom_steps; bins=edges_res, normalize=:pdf, alpha=0.6, label="No memory")
histogram!(p_res, res_mem_steps; bins=edges_res, normalize=:pdf, alpha=0.6, label="With memory")
xlabel!(p_res, "Residence time in core (steps)"); ylabel!(p_res, "PDF")
title!(p_res, "Residence-time distributions (gate = $(x_gate))")
display(p_res)

# Residence survival
ρ1, R1 = survival_curve(res_nom_time)
ρ2, R2 = survival_curve(res_mem_time)
p_surv_res = plot(ρ1, R1; lw=2, label="No memory")
plot!(p_surv_res, ρ2, R2; lw=2, label="With memory")
xlabel!(p_surv_res, "Residence time (time units)"); ylabel!(p_surv_res, "Survival 1-ECDF")
title!(p_surv_res, "Residence survival (gate = $(x_gate))")
display(p_surv_res)


# (f) Optional: KDE overlays for inter-switch, truncated to x >= 0
kd_no = kde(Float64.(is_nom_steps)); mask_no = kd_no.x .>= 0
kd_me = kde(Float64.(is_mem_steps)); mask_me = kd_me.x .>= 0
p_kde = histogram(is_nom_steps; bins=edges_is, normalize=:pdf,
                  fillalpha=0.35, linealpha=0.5, label="No memory (hist)")
histogram!(p_kde, is_mem_steps; bins=edges_is, normalize=:pdf,
           fillalpha=0.35, linealpha=0.5, label="With memory (hist)")
plot!(p_kde, kd_no.x[mask_no], kd_no.density[mask_no]; lw=3, label="No memory (KDE)")
plot!(p_kde, kd_me.x[mask_me], kd_me.density[mask_me]; lw=3, label="With memory (KDE)")
vline!([median(is_nom_steps)]; ls=:dash, label="median no-mem")
vline!([median(is_mem_steps)]; ls=:dash, label="median mem")
xlabel!(p_kde, "Inter-switch interval (steps)"); ylabel!(p_kde, "PDF")
title!(p_kde, "Committed inter-switch distributions")
display(p_kde)
