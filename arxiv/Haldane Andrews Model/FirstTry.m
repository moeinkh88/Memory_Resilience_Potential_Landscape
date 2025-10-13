% Fractional substrate inhibition model in MATLAB
% Caputo D^alpha x(t) = V*x/(Km + x + x.^2/Ks) - beta*x
% Requires: FDE_PI2_IM.m on your path (Garrappa)

clear; clc;

% --- parameters ---
V     = 1.2;      % max uptake rate
Km    = 0.3;      % half saturation
Ks    = 2.0;      % inhibition constant
beta  = 0.10;     % linear loss
par   = [V Km Ks beta];

% --- fractional order and numerics ---
al      = 0.85;     % 0 < alpha <= 1
X0      = 0.05;     % initial substrate
Tfinal  = 150;      % final time
h       = 0.02;     % step size

% --- RHS and Jacobian (scalar state) ---
F  = @(t,x,par)  par(1).*x ./ (par(2) + x + (x.^2)./par(3)) - par(4).*x;

% df/dx for f(x) = V*x/g - beta*x, g = Km + x + x^2/Ks
% Here df/dx = V*(Km - x.^2/Ks)/g^2 - beta
JF = @(t,x,par)  par(1).*( par(2) - (x.^2)./par(3) ) ./ ...
                ( par(2) + x + (x.^2)./par(3) ).^2 - par(4);

% --- solve ---
[t, X] = FDE_PI2_IM(al, F, JF, 0, Tfinal, X0, h, par);

% --- plot ---
figure; 
plot(t, X, 'LineWidth', 2);
xlabel('time'); ylabel('x(t)');
title(sprintf('Fractional substrate inhibition, \\alpha = %.2f', al));
grid on;

% --- optional: compare several alpha values ---
alphas = [0.6 0.85 1.0];
figure; hold on;
for a = alphas
    [tA, xA] = FDE_PI2_IM(a, F, JF, 0, Tfinal, X0, h, par);
    plot(tA, xA, 'LineWidth', 1.8, 'DisplayName', sprintf('\\alpha = %.2f', a));
end
xlabel('time'); ylabel('x(t)');
title('Effect of fractional order \alpha');
legend('show'); grid on;
