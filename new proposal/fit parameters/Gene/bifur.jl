using Plots, NLsolve

# ---------------- parameters struct ----------------
Base.@kwdef struct Params
    K::Float64
    c::Float64
    n::Float64
end

# ---------------- model ----------------
f(x, a, par::Params) = (a * x^par.n) / (par.K + x^par.n) - par.c * x

dfdx(x, a, par::Params) = begin
    if x == 0
        -par.c
    else
        (a * par.n * par.K * x^(par.n-1)) / (par.K + x^par.n)^2 - par.c
    end
end

# ---------------- equilibria solver ----------------
function equilibria_for_a(a; par=Params(K=1.0, c=0.5, n=2.0), xmax=50.0)
    roots = Float64[0.0]
    if isapprox(par.n, 2.0; atol=1e-12)
        # exact quadratic for n=2
        A, B, C = par.c, -a, par.c * par.K
        Δ = B^2 - 4A*C
        if Δ > 0
            for r in ((-B - sqrt(Δ))/(2A), (-B + sqrt(Δ))/(2A))
                if r > 0; push!(roots, r); end
            end
        end
        sort!(roots); return roots
    end

    guesses = range(1e-6, xmax; length=80)
    found = copy(roots)
    for x0 in guesses
        try
            sol = nlsolve(xx -> [f(xx[1], a, par)],
                          xx -> [dfdx(xx[1], a, par)],
                          [x0]; autodiff=:forward)
            x★ = sol.zero[1]
            if x★ ≥ 0 && all(abs(x★-r)>1e-4 for r in found)
                push!(found, x★)
            end
        catch; end
    end
    sort!(found); return found
end

is_stable(x, a, par) = dfdx(x, a, par) < 0

# ---------------- main sweep ----------------
par = Params(K=1.0, c=0.6, n=2.0)
athresh = 2 * par.c * sqrt(par.K)   # threshold for n=2
arange = range(0.1, 5.0; length=250)

a_stab = Float64[]; x_stab = Float64[]
a_unst = Float64[]; x_unst = Float64[]

for aval in arange   # <<< renamed loop variable, NOT "p"
    for xeq in equilibria_for_a(aval; par=par)
        if is_stable(xeq, aval, par)
            push!(a_stab, aval); push!(x_stab, xeq)
        else
            push!(a_unst, aval); push!(x_unst, xeq)
        end
    end
end

# ---------------- plot ----------------
plt = plot(legend=false, xlabel="a", ylabel="Equilibrium x*",
           title="Bifurcation diagram (K=$(par.K), c=$(par.c), n=$(par.n))")
scatter!(a_stab, x_stab, ms=3, color=:blue)
scatter!(a_unst, x_unst, ms=3, marker=:xcross, color=:red)
vline!([athresh], linestyle=:dash, color=:black)
display(plt)
