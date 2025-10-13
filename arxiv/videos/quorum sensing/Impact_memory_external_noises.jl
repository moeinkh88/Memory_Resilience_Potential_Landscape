###############################
# Fractional noise experiments
# One dimensional bistable system
# D^α x = F(x) + σ ξ(t)
###############################

using FdeSolver
using Random, Statistics
using FFTW, StatsBase
using Plots

# ---------------------------------------
# 1. Problem setup
# ---------------------------------------
SEED   = 42
Random.seed!(SEED)

α_list = [1.0, 0.8, 0.6]     # memory orders
σ      = 0.25                # external noise level
T      = 400.0               # total time
h      = 0.01                # step (Nyquist cutoff ~ π/h)
tSpan  = [0.0, T]
tgrid  = 0:h:T
N      = length(tgrid)

# Bistable drift: F(x) = x - x^3  (double well at ±1, barrier at 0)
F(x) = x .- x .^ 3

# Single bandlimited white sequence η used for all α
# Piecewise constant on [t_k, t_{k+1}), scaled as 1/√h
η = randn(N) ./ sqrt(h)

# helper to sample the piecewise constant forcing at time t
@inline function force_at(t::Float64, η::Vector{Float64}, h::Float64)
    k = Int(floor(t / h)) + 1
    k = clamp(k, 1, length(η))
    return η[k]
end

# ---------------------------------------
# 2. Additive noise experiment
# ---------------------------------------
x0 = 1.0  # start near the right well

sols = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()

for α in α_list
    rhs = (t, x) -> F(x) .+ σ .* force_at(t, η, h)   # additive forcing
    t, x = FDEsolver(rhs, tSpan, x0, α; h=h)       # direct call requested
    sols[α] = (t, vec(x))
end

# ---------------------------------------
# 3. Metrics: sample variance after burn in, ACF, PSD
# ---------------------------------------
burn   = 100.0
iburn  = searchsortedfirst(tgrid, burn):length(tgrid)

sample_variance(x) = var(view(x, iburn))

# simple autocorrelation up to maxlag seconds
function acf_simple(x::Vector{Float64}, h::Float64; maxlag_sec::Float64=10.0)
    xz = x .- mean(x)
    maxlag = Int(round(maxlag_sec / h))
    lagidx = collect(0:maxlag)                # <-- vector, not a single Int
    ac = autocor(xz, lagidx)                  # StatsBase expects a vector here
    lags = lagidx .* h
    return lags, ac
end

# one sided PSD via FFT
function psd_onesided(x::Vector{Float64}, h::Float64)
    n    = length(x)
    xz   = x .- mean(x)
    X    = rfft(xz)
    # Simple density style scaling so that area roughly reflects variance
    S    = (abs.(X).^2) .* (2h / n)
    freqs = range(0, stop=1/(2h), length=length(S))
    return freqs, S
end

println("\nWithin well sample variance after t ≥ $(burn):")
for α in reverse(α_list)
    t, x = sols[α]
    v = sample_variance(x)
    println("  α = $(α):  var ≈ $(round(v, sigdigits=4))")
end

# ---------------------------------------
# 4. Plots
# ---------------------------------------
default(marker=false, legend=:topright, linewidth=2)

# time series
plt_ts = plot(title="Sample paths, same noise, different memory", xlabel="time", ylabel="x(t)")
for α in reverse(α_list)
    t, x = sols[α]
    plot!(plt_ts, t, x, label="α = $(α)")
end
display(plt_ts)

# autocorrelation
plt_ac = plot(title="Autocorrelation (t ≥ $(burn))", xlabel="lag", ylabel="ACF")
for α in reverse(α_list)
    _, x = sols[α]
    lags, ac = acf_simple(x[iburn], h; maxlag_sec=50.0)
    plot!(plt_ac, lags, ac, label="α = $(α)")
end
display(plt_ac)

# replace acf_simple with this version that also computes IACT
function acf_and_iact(x::Vector{Float64}, h::Float64; maxlag_sec::Float64=20.0)
    n = length(x)
    maxlag = min(Int(round(maxlag_sec / h)), n-1)
    xz = x .- mean(x)
    c0 = sum(abs2, xz) / n
    lagidx = 0:maxlag
    ac = [ sum(@view(xz[1:n-k]) .* @view(xz[1+k:n])) / (n - k) / c0 for k in lagidx ]  
    lags = collect(lagidx) .* h
    # Integrated autocorrelation time (rectangle rule)
    iact = h * (sum(ac) - 0.5*(ac[1] + ac[end]))  # simple trapezoid
    return lags, ac, iact
end

# plotting (semilog y helps you see the “long tiny tail”)
# semilog y ACF (positive part only)
plt_ac1 = plot(title="Autocorrelation (log x)",
               xlabel="lag", ylabel="ACF", xscale=:log10)
for α in reverse(α_list)
    _, x = sols[α]
    lags, ac, iact = acf_and_iact(x[iburn], h; maxlag_sec=20.0)
    plot!(plt_ac1, lags[2:end], ac[2:end],
          label="α=$(α), IACT≈$(round(iact,sigdigits=3))")
end
display(plt_ac1)



function within_well_mask(t, x; thresh=0.4)
    return findall(xx -> xx > thresh, x)
end

println("\nWithin-well (x>0.4) sample variance after t≥$(burn):")
for α in α_list
    t, x = sols[α]
    idx = within_well_mask(t, x; thresh=0.4)
    idx = filter(i -> i ≥ first(iburn), idx)  # also apply burn-in
    if length(idx) > 200
        v = var(x[idx])
        println("  α=$(α): var≈", round(v, sigdigits=5), "  (n=", length(idx), ")")
    else
        println("  α=$(α): not enough within-well samples (n=", length(idx), ")")
    end
end


# power spectra on log log with reference slopes
plt_psd = plot(xscale=:log10, yscale=:log10,
               title="PSD, same noise, different memory",
               xlabel="frequency", ylabel="S(ω)")

ωref = 0.1
ref_curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()  # α => (f, ref)

for α in reverse(α_list)
    _, x = sols[α]
    f, S = psd_onesided(x[iburn], h)

    # PSD
    plot!(plt_psd, f[2:end], S[2:end], label="α = $(α)")

    # stash the reference slope for later
    idx = argmin(abs.(f .- ωref))
    if idx ≥ 2
        c   = S[idx] * (f[idx]^(2α))      # S ≈ c * ω^{-2α} through (ωref, Sref)
        ref = c .* (f .^ (-2α))
        ref_curves[α] = (f[2:end], ref[2:end])
    end
end

# draw all slope refs on top
for (α, (f2, ref2)) in sort(collect(ref_curves); by=first)
    plot!(plt_psd, f2, ref2, label="slope −2α ref (α=$(α))", ls=:dash, lw=1.8)
end

display(plt_psd)


# ---------------------------------------
# 5. Optional: multiplicative noise experiment
#    D^α x = F(x) + g(x) ξ(t)
# ---------------------------------------
do_multiplicative = true
if do_multiplicative
    g(x) = 0.25 .* (1.0 .+ 0.5 .* tanh.(x))   # smooth state dependent amplitude

    sols_mult = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for α in α_list
        rhs = (t, x) -> F(x) .+ g(x) .* force_at(t, η, h)
        t, x = FDEsolver(rhs, tSpan, x0, α; h=h)
        sols_mult[α] = (t, vec(x))
    end

    println("\nMultiplicative noise: effective level near x≈1 is g(1) = $(round(g(1.0), sigdigits=4))")
    println("Within well sample variance after t ≥ $(burn):")
    for α in α_list
        _, x = sols_mult[α]
        v = sample_variance(x)
        println("  α = $(α):  var ≈ $(round(v, sigdigits=4))")
    end

    plt_ts2 = plot(title="Multiplicative noise: sample paths", xlabel="time", ylabel="x(t)")
    for α in reverse(α_list)
        t, x = sols_mult[α]
        plot!(plt_ts2, t, x, label="α = $(α)")
    end
    display(plt_ts2)
end

# dwell times in the right basin (x > 0), simple sign-based wells
# classify sample with a guard band δ around 0
@inline basin_class(x; δ=0.1) = x >  δ ? :right :
                                x < -δ ? :left  : :middle

# dwell times in a chosen basin with hysteresis; do NOT push the last partial dwell
function dwell_times(x::Vector{Float64}, h::Float64; side::Symbol=:right, δ::Float64=0.1)
    want = side                      # :right or :left
    times = Float64[]
    inside = (basin_class(x[1]; δ=δ) == want)
    t_acc = 0.0
    for k in 2:length(x)
        bc = basin_class(x[k]; δ=δ)

        if inside
            t_acc += h
            # exit only when we fully cross into the opposite basin (skip the middle band)
            if (want == :right && bc == :left) || (want == :left && bc == :right)
                push!(times, t_acc)          # complete dwell
                t_acc = 0.0
                inside = false
            end
        else
            # start a dwell only when we fully enter the target basin
            if bc == want
                inside = true
                t_acc = 0.0
            end
        end
    end
    # do NOT push t_acc here: the last dwell is right-censored if still inside
    return times
end

# survival curve S(t) = P(T ≥ t), simple step estimator
function survival(times::Vector{Float64})
    isempty(times) && return (Float64[], Float64[])
    st = sort(times)
    n  = length(st)
    t  = unique(st)
    S  = [sum(st .>= τ)/n for τ in t]
    return t, S
end


plt_surv = plot(title="Dwell-time survival (right basin)",
                xlabel="time in basin", ylabel="S(t)", yscale=:log10)
for α in reverse(α_list)
    _, x = sols[α]
    times = dwell_times(x[iburn], h; side=:right)
    tS, S = survival(times)
    plot!(plt_surv, tS, S, label="α=$(α)")
end
display(plt_surv)


