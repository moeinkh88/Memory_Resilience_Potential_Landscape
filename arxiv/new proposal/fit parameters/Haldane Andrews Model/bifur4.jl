
using Plots

# parameters (you can adjust!)
Km = 1.0
Ks = 0.1
β  = 0.5

# numerator polynomial p(x)
function p(x, V, Km, Ks, β)
    return -(β/Ks)*x^3 - β*x^2 + (V - β*Km)*x
end

# derivative for stability test
function dpdx(x, V, Km, Ks, β)
    return -(3β/Ks)*x^2 - 2β*x + (V - β*Km)
end


using Polynomials

function equilibria(V; Km=1.0, Ks=1.0, β=1.0)
    # coefficients of cubic: a3 x^3 + a2 x^2 + a1 x + a0
    a3 = -β/Ks
    a2 = -β
    a1 = V - β*Km
    a0 = 0.0
    poly = Polynomial([a0,a1,a2,a3])
    roots_real = real.(roots(poly))
    roots_pos = filter(x -> isreal(x) && x ≥ 0, roots_real)
    return roots_pos
end


Vvals = range(0.0, 5.0, length=200)   # scan V
x_stable = Float64[]
V_stable = Float64[]
x_unstable = Float64[]
V_unstable = Float64[]

for V in Vvals
    eqs = equilibria(V; Km=Km, Ks=Ks, β=β)
    for x in eqs
        slope = dpdx(x, V, Km, Ks, β)
        if slope < 0
            push!(x_stable, x)
            push!(V_stable, V)
        else
            push!(x_unstable, x)
            push!(V_unstable, V)
        end
    end
end


plot(V_stable, x_stable, seriestype=:scatter, color=:blue, label="Stable equilibria")
plot!(V_unstable, x_unstable, seriestype=:scatter, color=:red, marker=:x, label="Unstable equilibria")
xlabel!("V (max uptake)")
ylabel!("Equilibrium substrate/biomass (x*)")
title!("Bifurcation diagram of Haldane–Andrews model")
