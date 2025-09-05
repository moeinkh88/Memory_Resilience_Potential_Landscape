using CSV, DataFrames, Statistics
using Optim, FdeSolver, StatsBase


dataset = CSV.read("sheet_Y.csv", DataFrame) # all data of confirmed

# ============== load parameters ==============
# Prefer two csv files:
#   interaction.csv   -> 11x11 matrix M
#   growth.csv        -> one column named "mu" with 11 rows
# read the file
dfM = CSV.read("interaction.csv", DataFrame)

# if the first column holds taxa names, drop it
if eltype(dfM[!, 1]) <: AbstractString
    select!(dfM, Not(1))
end

# make sure every remaining column is Float64
for c in names(dfM)
    if eltype(dfM[!, c]) <: AbstractString
        # trim spaces and convert text to numbers
        dfM[!, c] = parse.(Float64, strip.(dfM[!, c]))
    else
        dfM[!, c] = Float64.(dfM[!, c])
    end
end

M = Matrix{Float64}(dfM)
μd = CSV.read("growth.csv", DataFrame)
μ  = Vector{Float64}(μd[:, 1])   # or μd.mu if column is named "mu"
Susc = CSV.read("sc.csv", DataFrame)


@assert size(M) == (11, 11)
@assert length(μ) == 11


# ============== gLV model for FdeSolver ==============
# dx_i/dt = μ_i * x_i + x_i * sum_j M_ij * x_j
# Vectorized form: dx = x .* (μ .+ M*x)
function gLV1(t, x, par)
    M = par.M
    μ = par.μ
    return x .* (μ .+ M * x)
end

function gLV2(t, x, par)
    M = par.M
    μ = par.μ
    return x .* (μ .+ M * x)
end

function gLV3(t, x, par)
    M = par.M
    μ = par.μ
    return x .* (μ .+ M * x)
end

# pack parameters in a NamedTuple
par = (M = M, μ = μ)

name = "Clostridium_difficile"
idx = findfirst(==(name), interaction.x1)

# ============== initial conditions and orders ==============
L = 11
x0 = fill(1.0, L)         # replace with your initial abundances
tspan = [0.0, 50.0]       # replace with your time window
α = fill(0.95, L)         # orders of derivatives, you can fit these later

# optional step size
h = 0.1