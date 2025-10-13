using Plots

# fixed parameters except V
Km   = 0.3
Ks   = 2.0
beta = 1

Vcrit = beta * Km
Vmin, Vmax = 0.02, 10*Vcrit
Vgrid = range(Vmin, Vmax, length=700)

# storage
x_zero_stable   = fill(NaN, length(Vgrid))
x_zero_unstable = fill(NaN, length(Vgrid))
x_pos_stable    = fill(NaN, length(Vgrid))

# df/dx at (x, V)
fprime(x, V) = begin
    g = Km + x + x^2 / Ks
    V*(Km - x^2 / Ks)/g^2 - beta
end

for (i, V) in enumerate(Vgrid)
    # x = 0 branch
    if V < Vcrit
        x_zero_stable[i] = 0.0
    else
        x_zero_unstable[i] = 0.0
    end

    # positive equilibrium from (1/Ks)x^2 + x + (Km - V/beta) = 0
    c = Km - V/beta
    disc = 1.0 - 4.0*(c / Ks)          # = 1 - 4*(1/Ks)*(Km - V/beta)
    if disc >= -1e-12
        disc = max(disc, 0.0)
        xplus = (Ks/2) * (-1 + sqrt(disc))
        if xplus > 0
            if fprime(xplus, V) < 0     # stability in one dimension
                x_pos_stable[i] = xplus
            end
        end
    end
end

# plot
plt = plot(Vgrid, x_zero_stable, lw=2, label="x* = 0 stable",
           xlabel="V", ylabel="equilibrium x*", title="Bifurcation diagram vs V")
plot!(Vgrid, x_zero_unstable, lw=2, ls=:dash, label="x* = 0 unstable")
plot!(Vgrid, x_pos_stable, lw=2, label="x* > 0 stable")
vline!([Vcrit], ls=:dot, label="Vcrit = beta*Km")

display(plt)
println("Vcrit = beta*Km = ", Vcrit)
