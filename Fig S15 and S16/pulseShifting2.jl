# pulse perturbation for Herbivory model to show how memory lead to shifting
using FdeSolver
using Plots
using LaTeXStrings

## inputs
tSpan = [0, 100]   # time span
h = 0.01           # time step
# order of derivatives
β1 = 1.0 # order of the derivative
β2 = 0.9 # order of the derivative

X0 = .8433  # initial condition
# .8435, 1.9568

## System definition

# parametrisation
par = [.8,
       3.0,
       .2,
      .6]

function F(t, x, par)
    r = par[1] # 
    K = par[2] # 
    A = par[3] # 
    B = par[4] # 

    # Define perturbation parameters
    start_time = 45    # when perturbations start
    period = 100         # time between start of each pulse
    duration = 30       # length of each pulse
    pulse_amp = 0.1    # amplitude of perturbation

    # Only apply perturbation after start_time
    if t < start_time
        b = B
    else
        # Calculate time since perturbations started
        t_shifted = t - start_time
        t_mod = t_shifted % period
        
        if t_mod <= duration
            b = B + pulse_amp
        else
            b = B
        end
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

p1 = vspan([45, 75], color = myColor1[1], labels = :false)
plot!(t[1:10:end], x1[1:10:end], color = myColor2[1], linewidth = 2, yaxis = "Population dynamics",
      labels = "X,      0", legendtitle = "Population, Memory", legendposition = :bottomleft, legendtitlefont = font(10))
plot!(t[1:10:end], x2[1:10:end], color = myColor2[2], xaxis = "Time",
      labels = "X,      0.1", linewidth = 2)

savefig(p1, "LeadShiftUni.svg")