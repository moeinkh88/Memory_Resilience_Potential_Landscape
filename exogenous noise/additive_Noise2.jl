using FdeSolver, Random, Statistics, Plots, FFTW

# reproducibility
Random.seed!(123)

# model and simulation params
α_mem   = 0.7         # memory case
α_nomem = 1.0         # no memory case
y0      = 1.0
T       = 1000.0
h       = 0.01
tSpan   = (0.0, T)

# noise strength for additive forcing  dx/dt = x - x^3 + ξ(t)
σ_add   = 0.1

"Piecewise-constant white-noise approximation ξ(t) over [t0,t1] with step h and sd σ."
function make_step_noise(t0::Real, t1::Real, h::Real, σ::Real; seed::Integer=123)
    N   = Int(ceil((t1 - t0)/h))
    rng = MersenneTwister(seed)
    vals = (σ/√h) .* randn(rng, N)    # scale so that ∫ ξ dt ~ σ dW
    t0f, hf = float(t0), float(h)
    return (t::Real)-> begin
        idx = clamp(Int(floor((t - t0f)/hf)) + 1, 1, N)
        return vals[idx]
    end
end

# one shared noise path for both α to make a fair comparison
ξ = make_step_noise(tSpan[1], tSpan[2], h, σ_add; seed=1)

# double well drift
f(x) = x - x^3

# RHS wrapper
F(t, x) = f.(x) .+ ξ(t)  # works for scalar or vector x

# integrate
t_vals, X_nomem = FDEsolver(F, collect(tSpan), y0, α_nomem; h=h)
_,      X_mem   = FDEsolver(F, collect(tSpan), y0, α_mem;   h=h)

# drop an initial transient of 100 time units
t_trim = 100.0
i0 = max(1, Int(floor(t_trim/h)) + 1)
Xn = @view X_nomem[i0:end]
Xm = @view X_mem[i0:end]
tn = @view t_vals[i0:end]

# sample variance
var_nomem = var(Xn)
var_mem   = var(Xm)
println("variance no memory   = ", var_nomem)
println("variance with memory = ", var_mem)

# simple autocorrelation up to a chosen lag (biased estimator)
function autocorr_biased(x::AbstractVector{<:Real}, maxlag::Int)
    n = length(x)
    μ = mean(x)
    σ2 = var(x)
    acf = Vector{Float64}(undef, maxlag + 1)
    acf[1] = 1.0
    for τ in 1:maxlag
        s = 0.0
        @inbounds for i in 1:(n-τ)
            s += (x[i] - μ)*(x[i+τ] - μ)
        end
        acf[τ+1] = s / ((n-τ) * σ2)
    end
    return acf
end

# ACF window
W_time = 200.0                    # 200 time units
maxlag = Int(round(W_time / h))   # convert to steps

acf_n = autocorr_biased(Xn, maxlag)
acf_m = autocorr_biased(Xm, maxlag)

println("ACF lag 1   no memory   = ", acf_n[2])
println("ACF lag 1   with memory = ", acf_m[2])
println("ACF lag 50  no memory   = ", acf_n[51])
println("ACF lag 50  with memory = ", acf_m[51])

# integrated autocorrelation time over a window W_time
function tau_int(acf::AbstractVector{<:Real}, h::Real, W_time::Real)
    K = min(length(acf)-1, Int(round(W_time/h)))
    return h * (1 + 2*sum(@view acf[2:K+1]))
end

println("τ_int over ", W_time, "  no memory   = ", tau_int(acf_n, h, W_time))
println("τ_int over ", W_time, "  with memory = ", tau_int(acf_m, h, W_time))

# plotting

# 1) time segment
t0_plot, t1_plot = 400.0, 450.0
i_start = clamp(Int(floor((t0_plot - t_vals[1])/h))+1, 1, length(t_vals))
i_end   = clamp(Int(floor((t1_plot - t_vals[1])/h))+1, 1, length(t_vals))
p1 = plot(t_vals[i_start:i_end], X_mem[i_start:i_end], label="with memory α=$(α_mem)",  color=:indianred3)
plot!(p1, t_vals[i_start:i_end], X_nomem[i_start:i_end], label="no memory α=1",   color=:dodgerblue1)
xlabel!(p1, "time"); ylabel!(p1, "x(t)"); title!(p1, "trajectory segment")

# 2) ACF log log
lags_time = (0:maxlag) .* h
p2 = plot(lags_time[2:end], acf_n[2:end], xscale=:log10, lw=2, label="no memory")
plot!(p2, lags_time[2:end], acf_m[2:end], lw=2, label="with memory")
xlabel!(p2, "lag time"); ylabel!(p2, "ACF"); title!(p2, "autocorrelation")

# 3) coarse graining: variance and τ_int vs block size
function block_average(x::AbstractVector{<:Real}, m::Int)
    n = length(x) ÷ m
    return [mean(@view x[(i-1)*m+1 : i*m]) for i in 1:n]
end

m_list = [1, 5, 10, 20, 50, 100, 150, 200]
W_cg   = 50.0

tau_no = Float64[]; tau_me = Float64[]
var_no = Float64[]; var_me = Float64[]
for m in m_list
    Xn_b = block_average(Xn, m)
    Xm_b = block_average(Xm, m)
    h_eff = h * m
    acf_n_b = autocorr_biased(Xn_b, Int(round(W_cg / h_eff)))
    acf_m_b = autocorr_biased(Xm_b, Int(round(W_cg / h_eff)))
    push!(tau_no, tau_int(acf_n_b, h_eff, W_cg))
    push!(tau_me, tau_int(acf_m_b, h_eff, W_cg))
    push!(var_no, var(Xn_b))
    push!(var_me, var(Xm_b))
end

p3 = plot(m_list, tau_no, lw=2, marker=:o, label="no memory")
plot!(p3, m_list, tau_me, lw=2, marker=:o, label="with memory")
xlabel!(p3, "block size m"); ylabel!(p3, "τ_int over W="*string(W_cg)); title!(p3, "persistence vs coarse graining")

p4 = plot(m_list, var_no, lw=2, marker=:o, label="no memory")
plot!(p4, m_list, var_me, lw=2, marker=:o, label="with memory")
xlabel!(p4, "block size m"); ylabel!(p4, "variance"); title!(p4, "variance after coarse graining")

display(plot(p1, p2, p3, p4, layout=(2,2)))

savefig(p1, "plots//additive/trajectories_var.svg")
savefig(p2, "plots/additive/ACF.svg")
savefig(p3, "plots/additive/persist.svg")
savefig(p4, "plots/additive/var_coarse_grain.svg")

#optional

# -------- 3) Power spectral density on log log axes with reference slopes --------

"""
Welch PSD with Hann window.
Returns single sided frequency f (Hz) and PSD Pxx (power per Hz).
"""
function welch_psd(x::AbstractVector{<:Real}, h; seglen::Int=16384, overlap=0.5)
    fs = 1/h
    N  = length(x)
    seglen = seglen > N ? max(256, 2^(floor(Int, log2(N)))) : seglen
    step   = max(1, Int(round(seglen * (1 - overlap))))

    # Hann window and its power normalization
    w = 0.5 .* (1 .- cos.(2π .* (0:seglen-1) ./ (seglen-1)))
    U = sum(w.^2) / seglen

    acc = nothing
    count = 0
    for s in 1:step:(N - seglen + 1)
        seg = @view x[s:s+seglen-1]
        y   = (seg .- mean(seg)) .* w
        Y   = rfft(y)
        # single sided PSD, scaled to power per Hz
        Pseg = (abs.(Y).^2) ./ (fs * seglen * U)
        acc = acc === nothing ? Pseg : acc .+ Pseg
        count += 1
    end
    P = acc ./ count
    f = range(0, fs/2, length=length(P))
    return f, P
end

# --- 3) PSD with better visibility ------------------------------------------
# go longer segments => finer frequency grid (min f = fs/seglen)
SEG = min(65536, 2^(floor(Int, log2(length(Xn)))))   # ~0.0015 Hz with h=0.01
f_n, Pn = welch_psd(Xn, h; seglen=16384, overlap=0.5) 
f_m, Pm = welch_psd(Xm, h; seglen=16384, overlap=0.5)

# skip DC for log axes
i1 = 2:length(f_n)

# draw references from low f, in muted greys, then plot PSDs on top
function add_ref_slope!(plt, f, P; slope::Float64, label::String,
                        f_anchor::Union{Nothing,Real}=nothing, ls=:dash, kwargs...)
    # choose anchor point (by freq), default = 60% index
    j = isnothing(f_anchor) ? round(Int, 0.60*length(f)) : argmin(abs.(f .- f_anchor))
    f0, C = f[j], P[j]
    yref = C .* (f ./ f0) .^ (-slope)
    plot!(plt, f, yref; ls=ls, lw=2, label=label, kwargs...)
end

# build plot, skip DC at f = 0 for log axes 
i1 = 2:length(f_n) 
p3 = plot(f_n[i1], Pn[i1]; xscale=:log10, yscale=:log10, lw=2, label="no memory α=1") 
plot!(p3, f_m[i1], Pm[i1]; lw=2, label="with memory α=$(α_mem)", legendposition=:bottomleft) 
xlabel!(p3, "frequency") 
ylabel!(p3, "Power spectral density") 
# title!(p3, "power spectral density") 

# reference slopes: −2 for α=1 and −2α for memory 

# slope guides first (so they don’t hide PSD), anchored around 2 Hz
add_ref_slope!(p3, f_n[i1], Pn[i1]; slope=2.0,     label="slope −2 ref",
               f_anchor=2.0, color=:blue)
add_ref_slope!(p3, f_m[i1], Pm[i1]; slope=2*α_mem, label="slope −2α ref",
               f_anchor=2.0, ls=:dashdot, color=:red)


xlabel!(p3, "frequency")
ylabel!(p3, "Power spectral density")
# optional: tighten ranges so Nyquist bend doesn’t dominate
# ylims!(p3, 1e-10, 1e-4)

display(p3)
savefig(p3, "plots/additive/PSD.svg")