using CSV, DataFrames, Statistics
using Optim, FdeSolver, StatsBase


dataset = CSV.read("sheet_Y.csv", DataFrame) # all data of confirmed

# helper function to get the second pertubation value in population 3
function get_value(df, taxon::String; pop::Int, time::Int, rep::Int)
    # row of the taxon
    row_idx = findfirst(==(taxon), df[!, 1])

    # candidate columns starting with the right population number
    cols = filter(c -> startswith(String(c), string(pop)), names(df)[2:end])

    # select the one matching replicate and time
    for c in cols
        if df[1, c] == rep && df[3, c] == time
            return df[row_idx, c]
        end
    end
    error("No match for taxon=$taxon, pop=$pop, time=$time, rep=$rep")
end

# Example:
val = get_value(dataset, "Clostridium_difficile"; pop=3, time=2, rep=1)


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
Susc = CSV.read("suscep.csv", DataFrame)
Sus = Vector{Float64}(Susc[:, 1]) 

# ============== gLV model for FdeSolver ==============
# dx_i/dt = μ_i * x_i + x_i * sum_j M_ij * x_j
# Vectorized form: dx = x .* (μ .+ M*x)
function gLV1(t, x, par)
    M = par.M
    μ = par.μ
    return x .* (μ .+ M * x)
end

function gLV2(t, x, par)
    M, μ, Sus = par.M, par.μ, par.Sus
    μ1 = (t > 0.1 && t <= 1.0) ? (μ .- Sus) : μ
    
    return x .* (μ1 .+ M * x)
end

function gLV3(t, x, par)
    M, μ, Sus, value = par.M, par.μ, par.Sus, par.value
    μ1 = (t > 0.1 && t <= 1.0) ? (μ .- Sus) : μ
    if t > 2.1 && t <2.2
         x[idx] = value
    end
    return x .* (μ1 .+ M * x)
end

# pack parameters in a NamedTuple
par = (M = M, μ = μ, Sus = Sus)

# ============== initial conditions and orders ==============
# x0 = fill(1.0, L)         # replace with your initial abundances
# tspan = [0.0, 50.0]       # replace with your time window
# α = fill(0.95, L)         # orders of derivatives, you can fit these later

# optional step size
h = 0.1

## helper function to get the initial conditions for all taxa
# pick the column for (pop, time, rep)
select_col(df; pop::Int, time, rep) = begin
    cols = filter(c -> startswith(String(c), string(pop)), names(df)[2:end-2])
    for c in cols
        t = df[3, c]; r = df[1, c]
        if !ismissing(t) && !ismissing(r) && r == rep && t == time
            return c
        end
    end
    error("No column for pop=$pop, time=$time, rep=$rep")
end


function get_init_alltaxa(df; pop::Int, time, rep::Int)
    col = select_col(df; pop=pop, time=time, rep=rep)
    skip = Set(["Replicate","ID","time (in d)"])
    rows = [i for i in 1:nrow(df)-2 if
        df[i,1] isa AbstractString &&
        !(df[i,1] in skip) &&
        (df[i, col] isa Real || !ismissing(df[i, col]))
    ]
    names = String.(df[rows, 1])
    vals  = Float64.(df[rows, col])
    return names, vals
end

# Example:
taxa_names, x0 = get_init_alltaxa(dataset; pop=3, time=1, rep=1);


##try simulation