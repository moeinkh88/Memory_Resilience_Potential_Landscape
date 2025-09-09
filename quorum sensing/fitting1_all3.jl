using FdeSolver
using Optim, StatsBase
using Plots
using DataFrames, CSV

# Load the dataset
all_rows = CSV.read("qs_bistable_dataset.csv", DataFrame)

# Define scenarios and their perturbation times
scenarios_fit = ["S1", "S2", "S3"]
pert_times = Dict("S1" => (3.0, 5.0), "S2" => (3.0, 6.0), "S3" => nothing)

# Extract data for each scenario with alpha=0.94, replicate=1, using A_true
data_dict = Dict()
t_data = nothing  # Will be the same for all
for sc in scenarios_fit
    df = all_rows[(all_rows.scenario .== sc) .& (all_rows.alpha .== 0.94), :]
    dfr = df[df.replicate .== 1, :]
    if t_data === nothing
        t_data = dfr.t
    end
    data_dict[sc] = (A_data = dfr.A_true, Ainit = dfr.A_true[1])
end

# Reuse definitions
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

# Interpolation function
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

# Solve model for a given scenario and parameters
function solve_model(sc, theta)
    V, K, A0, alpha, rho_b1, rho_p1, rho_b2, rho_p2, rho_b3 = theta
    rho_base = (sc == "S1" ? rho_b1 : (sc == "S2" ? rho_b2 : rho_b3))
    pert = pert_times[sc]
    rho_fun = if pert === nothing
        make_rho_fun(rho_base; pert=nothing)
    else
        rho_pert = (sc == "S1" ? rho_p1 : rho_p2)
        make_rho_fun(rho_base; pert=(pert..., rho_pert))
    end
    par = (V, K, A0, rho_fun)
    tSpan = [0.0, 20.0]
    h = 0.01
    nc = 3
    Ainit = data_dict[sc].Ainit
    t, Aapp = FDEsolver(QS_RHS, tSpan, Ainit, alpha, par; h=h, nc=nc)
    A_model = interp_at(t, Aapp, t_data)
    return A_model
end

# Loss function: sum of normalized RMSD over all scenarios
function loss(theta)
    total_loss = 0.0
    for sc in scenarios_fit
        A_model = solve_model(sc, theta)
        A_data = data_dict[sc].A_data
        total_loss += rmsd(A_model, A_data; normalize=true)
    end
    return total_loss
end

# Initial guess: [V, K, A0, alpha, rho_b1, rho_p1, rho_b2, rho_p2, rho_b3]
theta0 = [1.0, 1.0, 0.1, 0.9, 0.3, 0.3, 0.4, 0.5, 0.4]

# Bounds
lower = [0.01, 0.01, 0.001, 0.2, 0.01, 0.01, 0.01, 0.01, 0.01]
upper = [5.0, 5.0, 0.5, 0.999, 0.99, 0.99, 0.99, 0.99, 0.99]

# Optimize
res = optimize(loss, lower, upper, theta0, Fminbox(LBFGS()),
               Optim.Options(outer_iterations = 20,
                             iterations=20,
                             show_trace=true,
                             show_every=1));

# Fitted parameters
fitted_params = Optim.minimizer(res)
println("Fitted parameters: ")
println("V = ", fitted_params[1])
println("K = ", fitted_params[2])
println("A0 = ", fitted_params[3])
println("alpha = ", fitted_params[4])
println("rho_base_S1 = ", fitted_params[5])
println("rho_pert_S1 = ", fitted_params[6])
println("rho_base_S2 = ", fitted_params[7])
println("rho_pert_S2 = ", fitted_params[8])
println("rho_base_S3 = ", fitted_params[9])

# Optional: Plot fitted vs data for each scenario
for sc in scenarios_fit
    A_model_fitted = solve_model(sc, fitted_params)
    A_data = data_dict[sc].A_data
    plt = plot(t_data, A_data, label="Data (A_true)", marker=:circle, title="Scenario $sc")
    plot!(plt, t_data, A_model_fitted, label="Fitted Model", lw=2)
    display(plt)
end