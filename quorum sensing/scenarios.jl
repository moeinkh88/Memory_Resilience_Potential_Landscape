#=
Three critical-scenario demos with memory (α<1) vs no-memory (α=1)

Model:   Dᵗ^α x(t) = a(t)*x - b*x^3 + u(t) + σ*η(t)
Double well when a(t)>0, b>0; unstable at x=0, stable near ±√(a/b)

Scenarios:
  1) Delayed outcomes after stress removal (temporary a(t) drop)
  2) Rollback after a near crossing (short u(t) pulse)
  3) Broadened hysteresis under slow ramps (compare α=1 vs α<1)
=#

using FdeSolver
using Random, Statistics
using Plots

# --- Utilities ---------------------------------------------------------------

"""
piecewise_constant_noise(t; t0, h, ξ) -> η(t)

Return η(t) that equals ξ[k] on grid cells of width h starting at t0.
"""
function piecewise_constant_noise(t; t0, h, ξ::AbstractVector)
    idx = clamp(Int(fld((t - t0)/h, 1)) + 1, 1, length(ξ))
    return ξ[idx]
end

"""
Build RHS f(t, x) for the cubic FDE with time-dependent a(t), u(t) and noise.
"""
function rhs_factory(a_fun, b::Real, u_fun, σ::Real, noise_fun)
    return (t, x) -> a_fun(t).*x .- b.*x.^3 .+ u_fun(t) .+ σ.*noise_fun(t)
end

# Helper to run a single trajectory
function run_traj(; α=0.6, tspan=[0.0, 150.0], h=0.02, x0=-1.0,
                  a_fun=t->1.0, b=1.0, u_fun=t->0.0, σ=0.02, rng_seed=42)

    # precompute noise on solver grid for reproducibility
    Random.seed!(rng_seed)
    nsteps = Int(cld(tspan[2]-tspan[1], h)) + 1
    ξ = randn(nsteps)                       # i.i.d. standard normal
    η(t) = piecewise_constant_noise(t; t0=tspan[1], h=h, ξ=ξ)

    F = rhs_factory(a_fun, b, u_fun, σ, η)

    tvec = [tspan[1], tspan[2]]
    t, x = FDEsolver(F, tvec, x0, α; h=h, nc=2)   # predictor-corrector
    return t, x
end

# Plot helpers
function plot_landscape!(t, a_fun; b=1.0, label_prefix="")
    # show instantaneous well locations ±√(a/b) when a>0
    a_of_t = a_fun.(t)
    m = @. sqrt(clamp(a_of_t/b, 0, Inf))
    plot!(t,  m, lw=1, ls=:dash, label="$(label_prefix)+√(a/b)")
    plot!(t, -m, lw=1, ls=:dash, label="$(label_prefix)-√(a/b)")
end

# --- Scenario 1: Delayed outcome after stress removal ------------------------

function scenario1()
    α = 0.6
    h = 0.02
    T = 150.0
    t1, t2 = 30.0, 70.0                # temporary stress window
    a0, alow = 1.0, 0.25               # lower barrier during stress
    a_fun = t -> (t ≥ t1 && t ≤ t2) ? alow : a0
    u_fun = t -> 0.0
    x0 = -sqrt(a0) + 0.05              # start in left valley
    σ = 0.03

    t, x = run_traj(α=α, tspan=[0.0, T], h=h, x0=x0,
                    a_fun=a_fun, b=1.0, u_fun=u_fun, σ=σ, rng_seed=11)

    p = plot(t, x, lw=2, label="x(t), α=$(α)",
             title="Scenario 1: Delayed outcome after stress removal",
             xlabel="time", ylabel="x")
    vline!([t1, t2], ls=:dot, lw=1, label=["stress on" "stress off"])
    plot_landscape!(t, a_fun; b=1.0)
    display(p)
    return t, x
end

# --- Scenario 2: Rollback after a near crossing ------------------------------

function scenario2()
    α = 0.6
    h = 0.02
    T = 120.0
    # Short pulse that temporarily tilts to the right, then vanishes
    tp, Δ = 20.0, 6.0
    a_fun = t -> 1.0
    u_fun = t -> (t ≥ tp && t ≤ tp+Δ) ? 0.9 : 0.0    # rightward push
    x0 = -1.0
    σ = 0.01

    t, x = run_traj(α=α, tspan=(0.0, T), h=h, x0=x0,
                    a_fun=a_fun, b=1.0, u_fun=u_fun, σ=σ, rng_seed=22)

    p = plot(t, x, lw=2, label="x(t), α=$(α)",
             title="Scenario 2: Rollback after a near crossing",
             xlabel="time", ylabel="x")
    vspan!([tp, tp+Δ], color=:gray, alpha=0.1, label="tilt pulse")
    plot_landscape!(t, a_fun; b=1.0)
    display(p)
    return t, x
end

# --- Scenario 3: Broadened hysteresis under slow ramps -----------------------

function scenario3()
    h = 0.02
    T = 240.0
    # Ramp a(t) up then down slowly
    a_min, a_max = 0.2, 1.2
    upT, holdT, downT = 90.0, 30.0, 90.0
    function a_fun(t)
        if t ≤ upT
            return a_min + (a_max - a_min)*t/upT
        elseif t ≤ upT + holdT
            return a_max
        elseif t ≤ upT + holdT + downT
            τ = t - (upT + holdT)
            return a_max - (a_max - a_min)*τ/downT
        else
            return a_min
        end
    end
    u_fun = t -> 0.0
    x0 = -sqrt(a_min) + 0.05
    σ = 0.02

    # Compare α=1 (no memory) vs α=0.6 (memory)
    t1, x1 = run_traj(α=1.0, tspan=(0.0, T), h=h, x0=x0,
                      a_fun=a_fun, b=1.0, u_fun=u_fun, σ=σ, rng_seed=33)
    t2, x2 = run_traj(α=0.6, tspan=(0.0, T), h=h, x0=x0,
                      a_fun=a_fun, b=1.0, u_fun=u_fun, σ=σ, rng_seed=33)

    p = plot(t1, x1, lw=2, label="α=1.0")
    plot!(t2, x2, lw=2, label="α=0.6",
          title="Scenario 3: Broadened hysteresis under slow ramps",
          xlabel="time", ylabel="x")
    plot_landscape!(t1, a_fun; b=1.0, label_prefix="")  # well positions
    display(p)

    # Optional: show x vs parameter to visualize the loop
    a_vals = a_fun.(t1)
    q = plot(a_vals, x1, seriestype=:scatter, ms=2, label="α=1.0",
             xlabel="a(t)", ylabel="x",
             title="Apparent hysteresis: x vs a(t)")
    plot!(a_fun.(t2), x2, seriestype=:scatter, ms=2, label="α=0.6")
    display(q)

    return (t1, x1, t2, x2)
end

# --- Run all -----------------------------------------------------------------

t1, x1 = scenario1()
t2, x2 = scenario2()
t3a, x3a, t3b, x3b = scenario3()
