# === (1) Install once (or do in pkg mode) ===
# using Pkg
# Pkg.add(["FdeSolver", "Plots"])

using FdeSolver
using Plots

# === (2) Model RHS: Caputo D^α x = f(t,x,par) ===
# par = (V, Km, Ks, beta_loss)
substrate_inhibition!(t, x, par) = begin
    V, Km, Ks, βloss = par
    V .* x ./ (Km .+ x .+ (x.^2) ./ Ks) .- βloss .* x
end

# === (3) Parameters, IC, and fractional order ===
V      = 1.0          # max uptake rate
Km     = 3.0          # half-saturation
Ks     = 5.0          # inhibition constant
βloss  = 1.50         # linear loss/decay
par    = (V, Km, Ks, βloss)

α      = 0.85         # fractional order (0 < α ≤ 1); set α=1 for the ODE case
x0     = 1.05         # initial substrate
tspan  = [0.0, 150.0] # time span
h      = 0.02         # step size (reduce if you need more accuracy)

# === (4) Solve ===
t, x = FDEsolver(substrate_inhibition!, tspan, x0, α, par; h=h)

# === (5) Plot ===
plot(t, x, lw=3, xlabel="time", ylabel="x(t)",
     title = "Fractional substrate inhibition (α = $(α))",
     label = "x(t)")

# === (6) (Optional) Compare different α values ===
alphas = [0.6, 0.85, 1.0]
plt = plot(xlabel="time", ylabel="x(t)",
           title="Effect of fractional order α on dynamics")
for a in alphas
    tA, xA = FDEsolver(substrate_inhibition!, tspan, x0, a, par; h=h)
    plot!(plt, tA, xA, lw=2, label="α = $(a)")
end
display(plt)
