using Plots

# --- fixed parameters ---
V  = 0.6
Km = 0.2
Ks = 2.0

βc = V / Km                    # threshold where x*=0 changes stability
βmin, βmax = 0.02, 1.5*βc
βgrid = range(βmin, βmax, length=600)

# storage
x_zero_stable   = fill(NaN, length(βgrid))
x_zero_unstable = fill(NaN, length(βgrid))
x_pos_stable    = fill(NaN, length(βgrid))   # positive physical root only

# helper: derivative df/dx at (x, β)
fprime(x, β) = begin
    g = Km + x + x^2 / Ks
    V*(Km - x^2 / Ks)/g^2 - β
end

for (i, β) in enumerate(βgrid)
    # x = 0 branch
    if β > βc
        x_zero_stable[i] = 0.0
    else
        x_zero_unstable[i] = 0.0
    end

    # nonzero equilibria from analytic quadratic:
    # (1/Ks) x^2 + x + (Km - V/β) = 0
    c = Km - V/β
    disc = Ks^2 - 4*Ks*c
    if disc >= 0
        xplus = (-Ks + sqrt(disc)) / 2      # the candidate positive root
        if xplus > 0
            # stability check (for 1D, sign of derivative decides)
            if fprime(xplus, β) < 0
                x_pos_stable[i] = xplus
            end
        end
    end
end

# --- plot bifurcation diagram ---
plt = plot(βgrid, x_zero_stable, lw=2, label="x* = 0 stable",
           xlabel="β (decay)", ylabel="equilibrium x*", title="Bifurcation diagram vs β")
plot!(βgrid, x_zero_unstable, lw=2, ls=:dash, label="x* = 0 unstable")
plot!(βgrid, x_pos_stable, lw=2, label="x* > 0 stable")
vline!([βc], ls=:dot, label="βc = V/Km")

display(plt)

# print the threshold for reference
println("βc = V/Km = ", βc)
