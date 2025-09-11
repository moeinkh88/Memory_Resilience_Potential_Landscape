# =========================
# Quorum sensing toy ODE with FdeSolver
# dA/dt = V*A^2/(K + A^2) + A0 - d(rho)*A
# d(rho) = 1 - (1 - rho) / (rho * (2 - rho))
# =========================

using FdeSolver      # pkg> add FdeSolver
using Plots          # pkg> add Plots

# --- model parameters you can tweak
V   = 3.0
K   = 1.0            # try 1.0 or a small epsilon to mimic K≈0
A0p = 0.05           # basal production (named A0p to avoid clash with A(t))
ρs  = [0.1, 0.2, 0.3, 0.4, 0.5]   # densities to compare
tSpan = [0.0, 10.0]  # time window
Ainit = 2         # initial A
α = 1.0              # derivative order (ODE case)

# decay function exactly as in your message
decay(ρ) = 0.1 + (1 - ρ) / (ρ * (2 - ρ))     # note: singular at ρ→0

# If you actually intended kA = 1 with the original paper’s d(ρ) = 1 + (1-ρ)/(ρ(2-ρ)),
# then replace the line above with:
# decay(ρ) = 1 + (1 - ρ) / (ρ * (2 - ρ))

# RHS for FdeSolver: F(t, y, par) returns dy/dt
# par = (V, K, A0p, ρ)
function QS_RHS(t, y, par)
    V, K, A0p, ρ = par
    A = y
    return V .* A.^2 / (K .+ A.^2) .+ A0p .- decay(ρ) * A
end

# solve and plot for multiple rho
plt1 = plot(title = "A(t) for different densities ρ", xlabel = "time", ylabel = "A(t)")
for ρ in ρs
    par = (V, K, A0p, ρ)
    t, Aapp = FDEsolver(QS_RHS, tSpan, Ainit, α, par; h = 0.01, nc = 3)
    plot!(plt1, t, Aapp, label = "ρ = $(round(ρ, digits=2))")
end
display(plt1)

# Bifurcation diagram for A vs rho using cubic equilibria and stability test

using Polynomials          # pkg> add Polynomials

# cubic coefficients: -d*A^3 + (V + A0)*A^2 - d*K*A + A0*K = 0
# Polynomials.jl uses coeffs from constant term upward
function cubic_coeffs(d; V=V, K=K, A0=A0p)
    return (A0*K, -d*K, (V + A0), -d)
end

# stability test: f'(A) = 2 V K A /(K + A^2)^2 - d
fprime(A, d; V=V, K=K) = 2V*K*A / (K + A^2)^2 - d

# sweep rho and collect equilibria
rhos = collect(range(0.05, 0.6, length=400))  # avoid rho near zero singularity
ρ_stable   = Float64[]
A_stable   = Float64[]
ρ_unstable = Float64[]
A_unstable = Float64[]

tol = 1e-9
for ρ in rhos
    d = decay(ρ)
    c0, c1, c2, c3 = cubic_coeffs(d)
    p = Polynomial([c0, c1, c2, c3])
    rts = roots(p)

    for z in rts
        if abs(imag(z)) < 1e-9
            A2 = real(z)
            if A2 >= 0
                fp = fprime(A2, d)
                if fp < 0
                    push!(ρ_stable, ρ);   push!(A_stable, A2)
                else
                    push!(ρ_unstable, ρ); push!(A_unstable, A2)
                end
            end
        end
    end
end

# plot bifurcation diagram
plt = plot(legend=:topright, xlabel="ρ", ylabel="Equilibria A",
           title="Bifurcation diagram A versus ρ", legendposition=:topleft)
scatter!(plt, ρ_stable,   A_stable;   label="stable",   markersize=1, markerstrokecolor=:auto)
scatter!(plt, ρ_unstable, A_unstable; label="unstable", markersize=1, markerstrokecolor=:auto)
display(plt)


# X = A^2 at the folds from f=0 and f'=0
disc = V^2 - 8*A0p*V
X1 = (K*(V - 2A0p) + sqrt(disc)) / (2*(V + A0p))
X2 = (K*(V - 2A0p) - sqrt(disc)) / (2*(V + A0p))
A1, A2 = sqrt(X1), sqrt(X2)

# d at folds from f'(A)=0
d1 = 2V*K*A1 / (K + X1)^2
d2 = 2V*K*A2 / (K + X2)^2

# invert d(ρ) = 0.1 + (1-ρ)/(ρ(2-ρ))  ->  ρ(d)
ρ(d) = begin
    y = d - 0.1
    ((2y + 1) - sqrt(4y^2 + 1)) / (2y)
end

println("fold rhos ≈ ", (ρ(d1), ρ(d2)))