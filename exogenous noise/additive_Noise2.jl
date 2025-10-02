using FdeSolver, Random, Statistics, Plots

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
p1 = plot(t_vals[i_start:i_end], X_mem[i_start:i_end], label="with memory α=$(α_mem)")
plot!(p1, t_vals[i_start:i_end], X_nomem[i_start:i_end], label="no memory α=1")
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


##optional
# τ_int as a function of the integration window W (in time units)
Wmax = 1000.0               # try 400 or 800 to see the flip clearly
L = Int(round(Wmax / h))   # max lag in steps

acf_n = autocorr_biased(Xn, L)
acf_m = autocorr_biased(Xm, L)

# cumulative trapezoid equivalent via cumulative sum of ACF (discrete form you used)
Ws  = (1:L) .* h
τnW = h .* (1 .+ 2 .* cumsum(acf_n[2:end]))   # length L
τmW = h .* (1 .+ 2 .* cumsum(acf_m[2:end]))

pτ = plot(Ws, τnW, lw=2, label="no memory")
plot!(pτ, Ws, τmW, lw=2, label="with memory")
xlabel!(pτ, "window W (time units)"); ylabel!(pτ, "τ_int(W)")
title!(pτ, "Integrated autocorrelation vs window")
display(pτ)


##
using Polynomials.PolyCompat
# pick a tail region (e.g., lags corresponding to 20–200 time units)
i1 = Int(round(20/h)); i2 = Int(round(200/h))
lags = (i1:i2) .* h
acf_tail = abs.(acf_m[i1+1:i2+1])            # avoid log of negative due to noise
coef = polyfit(log10.(lags), log10.(acf_tail), 1)
println("memory tail slope ≈ ", coef.coeffs[1], "  (ACF ~ t^{slope})")
