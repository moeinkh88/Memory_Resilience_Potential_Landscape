using Plots

# --- fixed parameters except Km ---
V    = .6
Ks   = 2.0
beta = 0.10

Km_crit = V / beta                         # threshold where x*=0 changes stability
Km_min, Km_max = 0.02, 1.5*Km_crit
Km_grid = range(Km_min, Km_max, length=700)

# storage
x_zero_stable   = fill(NaN, length(Km_grid))
x_zero_unstable = fill(NaN, length(Km_grid))
x_pos_stable    = fill(NaN, length(Km_grid))   # positive physical root only

# derivative df/dx at (x, Km)
fprime(x, Km) = begin
    g = Km + x + x^2 / Ks
    V*(Km - x^2 / Ks)/g^2 - beta
end

for (i, Km) in enumerate(Km_grid)
    # x = 0 branch
    if Km > Km_crit
        x_zero_stable[i] = 0.0
    else
        x_zero_unstable[i] = 0.0
    end

    # non zero equilibria come from:
    # (1/Ks) x^2 + x + (Km - V/beta) = 0
    c = Km - V/beta
    disc = 1.0 - 4.0*(c / Ks)             # discriminant
    if disc >= -1e-12                      # allow tiny negatives from roundoff
        disc = max(disc, 0.0)
        xplus = (Ks/2) * (-1 + sqrt(disc))   # candidate positive root
        if xplus > 0
            if fprime(xplus, Km) < 0        # stability in one dimension
                x_pos_stable[i] = xplus
            end
        end
    end
end

# --- plot ---
plt = plot(Km_grid, x_zero_stable, lw=2, label="x* = 0 stable",
           xlabel="Km", ylabel="equilibrium x*", title="Bifurcation diagram vs Km")
plot!(Km_grid, x_zero_unstable, lw=2, ls=:dash, label="x* = 0 unstable")
plot!(Km_grid, x_pos_stable, lw=2, label="x* > 0 stable")
vline!([Km_crit], ls=:dot, label="Km = V/beta")

display(plt)
println("Km_crit = V/beta = ", Km_crit)
