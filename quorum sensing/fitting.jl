using FdeSolver              # For solving the ODE system
using Plots                  # For visualization (optional)
using Random                 # For reproducibility
using DataFrames             # For data handling
using CSV                    # For saving data
using Optim                  # For optimization


# --------------------------
# read the data
# --------------------------
dataset = CSV.read("qs_bistable_dataset.csv", DataFrame) # all data
# --------------------------
# Global model parameters
# --------------------------
V   = 3.0
K   = 1.0
A0p = 0.05
t_end = 20.0
tSpan = (0.0, t_end)
sample_times = collect(0.0:1.0:t_end)
h = 0.01
nc_corr = 3
noise_sigma = 0.05
replicates = 3

# Fractional derivative orders
orders = [0.94, 0.95, 0.96, 0.97]

# --------------------------
# Loss function for fitting
# --------------------------
function loss_function(params, exp_id, df, sample_times)
    # Extract parameters to be fitted
    K, V, rho_base, rho_pert, A0p, α = params
    
    # Define rho schedule function
    perturbation = (exp_id == 1) ? (3.0, 5.0, rho_pert) : nothing
    ρ_fun = make_rho_fun(rho_base; pert=perturbation)
    
    # Set up parameters for FdeSolver
    par = (V, K, A0p, ρ_fun)
    
    # Simulate the model using FdeSolver
    t, Aapp = FDEsolver(QS_RHS, [tSpan[1], tSpan[2]], 2.0, α, par; h=h, nc=nc_corr)
    
    # Interpolate Aapp at sample_times
    A_true = interp_at(t, Aapp, sample_times)
    
    # Calculate the sum of squared errors (SSE) between A_true and A_obs
    dfr = df[df.exp_id .== exp_id, :]
    observed_A = dfr.A_obs
    sse = sum((A_true .- observed_A).^2)
    
    return sse
end

# --------------------------
# Optimization to fit the model parameters
# --------------------------
function fit_single_scenario(exp_id, df)
    # Initial guess for the parameters (K, V, rho, perturbation rho, A0p, α)
    initial_params = [1.0, 1.0, 0.35, 0.2, 0.05, 0.95]  # Example initial values
    
    # Define the loss function
    objective_function = (params) -> loss_function(params, exp_id, df, sample_times)
    
    # Define optimization options
    opt = Optim.Options(maxiter=50, xtol=1e-6)
    
    # Perform optimization using Nelder-Mead method
    result = optimize(objective_function, initial_params, NelderMead(), options=opt)
    
    # Extract the optimized parameters
    optimized_params = result.minimizer
    return optimized_params
end
# --------------------------
# Example: Fit for the first scenario
# --------------------------
df = all_rows[all_rows.scenario .== "S1", :]  # Example: scenario S1 (rho=0.35, A(0)=2.0, perturbation rho=0.2)
exp_id = 1  # Experiment ID for "S1"
optimized_params = fit_single_scenario(exp_id, df)

println("Optimized Parameters for Experiment $exp_id: ")
println("K: ", optimized_params[1])
println("V: ", optimized_params[2])
println("rho_base: ", optimized_params[3])
println("rho_pert: ", optimized_params[4])
println("A0p: ", optimized_params[5])
println("α: ", optimized_params[6])
