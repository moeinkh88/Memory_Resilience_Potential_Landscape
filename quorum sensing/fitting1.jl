using FdeSolver
using Optim, StatsBase
using Plots
using DataFrames, CSV

# Note: Assume the dataset CSV has been generated as per your original code.
# Load the dataset
all_rows = CSV.read("qs_bistable_dataset.csv", DataFrame)

# Filter for scenario S1, alpha=0.94, replicate=1, and use A_true as data
df = all_rows[(all_rows.scenario .== "S1") .& (all_rows.alpha .== 0.94), :]
dfr = df[df.replicate .== 1, :]
t_data = dfr.t
A_data = dfr.A_true
Ainit = A_data[1]  # Initial condition from data

# Known perturbation times
pert_times = (3.0, 5.0)

# Reuse definitions from your original code
decay(ρ) = 0.1 + (1 - ρ) / (ρ * (2 - ρ))

make_rho_fun(rho_base; pert=nothing) = t -> begin
    if pert === nothing
        return rho_base
    else
        t_on, t_off, rho_pert = pert
        return (t ≥ t_on && t ≤ t_off) ? rho_pert : rho_base
    end
end

function QS_RHS(t, y, par)
    V, K, A0, ρ_fun = par
    A = y
    ρ = ρ_fun(t)
    return V .* A.^2 ./ (K .+ A.^2) .+ A0 .- decay(ρ) .* A
end

# Reuse your interpolation function
function interp_at(ts_in, ys_in, tq_in)
    ts = collect(ts_in)
    ys = vec(ys_in)
    tq = collect(tq_in)
    n = length(ts)
    out = similar(tq)
    for (j, t) in enumerate(tq)
        if t ≤ ts[1]
            out[j] = ys[1]
        elseif t ≥ ts[end]
            out[j] = ys[end]
        else
            i = searchsortedlast(ts, t)
            i = max(1, min(i, n - 1))
            t1, t2 = ts[i], ts[i+1]
            y1, y2 = ys[i], ys[i+1]
            w = (t - t1) / (t2 - t1)
            out[j] = (1 - w) * y1 + w * y2
        end
    end
    return out
end

# Function to solve the model for given parameters
function solve_model(theta)
    V, K, A0, alpha, rho_base, rho_pert = theta
    rho_fun = make_rho_fun(rho_base; pert=(pert_times..., rho_pert))
    par = (V, K, A0, rho_fun)
    tSpan = [0.0, 20.0]
    h = 0.01
    nc = 3
    t, Aapp = FDEsolver(QS_RHS, tSpan, Ainit, alpha, par; h=h, nc=nc)
    A_model = interp_at(t, Aapp, t_data)
    return A_model
end

# Loss function (sum of squared errors)
function loss(theta)
    A_model = solve_model(theta)
    return     rmsd(A_model, A_data; normalize=:true)
end

# Initial guess for parameters [V, K, A0, alpha, rho_base, rho_pert]
theta0 = [1.0, 1.0, 0.1, 0.9, 0.3, 0.3]
p_lo_f_1=[.01, .01, 0.001, 0.2, 0.01, 0.08]; # lower bound 
p_up_f_1=[3.0, 3.0, 0.5, .999, 0.6, 0.3]; # upper bound 

# Optimize using BFGS (similar to Example 5 in FdeSolver docs)
res = optimize(loss,p_lo_f_1,p_up_f_1,theta0,Fminbox(LBFGS()), # LBFGS is suitable for large scale problems
# Result=optimize(loss,p_lo,p_up,pvec,SAMIN(rt=.99),
			Optim.Options(outer_iterations = 10,
						  iterations=20,
						  show_trace=true, # turn it true to see the optimization
						  show_every=1));

# Fitted parameters
fitted_params = Optim.minimizer(res)
println("Fitted parameters: ")
println("V = ", fitted_params[1])
println("K = ", fitted_params[2])
println("A0 = ", fitted_params[3])
println("alpha = ", fitted_params[4])
println("rho_base = ", fitted_params[5])
println("rho_pert = ", fitted_params[6])

# Optional: Plot fitted vs data for validation
A_model_fitted = solve_model(fitted_params)
plt = plot(t_data, A_data, label="Data (A_true)", marker=:circle)
plot!(plt, t_data, A_model_fitted, label="Fitted Model", lw=2)
display(plt)