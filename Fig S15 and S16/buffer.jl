# pulse perturbation for Herbivory model to show how memory lead to shifting
using FdeSolver
using Plots
using LaTeXStrings
using LinearAlgebra
using Interpolations
using Random

## inputs
tSpan = [0, 100]   # time span
h = 0.01           # time step
# order of derivatives
β1 = 1.0 # order of the derivative
β2 = 0.8 # order of the derivative

X0 = 1.9568  # initial condition
# .8435, 1.9568

## System definition

# parametrisation
par = [.8,
       3.0,
       .2,
      .6]

# Ornstein–Uhlenbeck process parameters
Random.seed!(1234) 
theta = 0.05   # mean reversion speed
mu    = .05  # long-term mean (similar to the pulse offset)
sigma = 0.05  # noise intensity

# Simulate the OU process over the entire time span.
dt_noise = h
t_noise = 0:dt_noise:100
noise = zeros(length(t_noise))
noise[1] = 0.0  # initial noise value

# Discretized OU simulation using the Euler–Maruyama method
for i in 2:length(t_noise)
    noise[i] = noise[i-1] + theta*(mu - noise[i-1])*dt_noise + sigma*sqrt(dt_noise)*randn()
end

# Create an interpolation function for the OU process
ou_interp = LinearInterpolation(t_noise, noise, extrapolation_bc=Line())

function F(t, x, par)
    r = par[1] # 
    K = par[2] # 
    A = par[3] # 
    B = par[4] # 

    # Instead of a fixed pulse, add OU noise in the interval [50, 100]
    b=copy(B)
    if t > 0 && t < 100
        b = B + ou_interp(t)
    end

    # ODE
    Fun = r .* x .* (1 .- x ./ K) .- b .* x ./ (x .+ A);

    return Fun

end

t, x1 = FDEsolver(F, tSpan, X0, β1, par, h = h)
_, x2 = FDEsolver(F, tSpan, X0, β2, par, h = h)

## Plotting

myColor1 = [:antiquewhite1 :burlywood1]
myColor2 = [:royalblue2 :firebrick1]

P=plot(t[1:10:end], x1[1:10:end], color = myColor2[1], linewidth = 2, yaxis = "Population dynamics",
      labels = "X,      0", legendtitle = "Population, Memory", legendposition = :bottomleft, legendtitlefont = font(10))
plot!(t[1:10:end], x2[1:10:end], color = myColor2[2], xaxis = "Time",
      labels = "X,      0.2", linewidth = 2)

savefig(P, "BuferUni.svg")