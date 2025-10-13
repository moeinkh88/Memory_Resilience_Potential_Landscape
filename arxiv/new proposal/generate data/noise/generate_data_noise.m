%% Herbivore fractional-order bistable model with additive noise (fixed indexing)
% ---------------------------------------------------------------
% * Fixes array‑index error when n = 1 in fractional derivative loop.
% * The loop now uses kmax = min(n‑1, M‑1) so that x(n‑k) ≥ x(1).
% * Everything else identical to the previous version.
% ---------------------------------------------------------------

%% House‑keeping
clc; clear; close all;
rng(123)                           % reproducible randomness

%% Global coefficients (visible to fun/Jfun)
global r K A
r = 0.8;
K = 6;
A = 0.2;

%% Numerical settings
al    = 0.8;      % fractional order (0<al<=1)
sigma  = 1e-2;    % noise strength (try 1e-3 … 1e-2)

t0 = 0;           % start time
T  = 2;         % final time
h  = 0.01;        % step size

%% --- Bifurcation diagram to locate multistable window --------
syms c xx
f_eqn = r*xx*(1-xx/K) - c*xx/(xx + A);
sc    = solve(f_eqn, xx);           % three equilibria (symbolic)

% Convert second and third roots to numeric fns
sc2 = matlabFunction(sc(2));
sc3 = matlabFunction(sc(3));

% Find saddle‑node intersection of sc2 and sc3
c_vec = 0:1e-6:5;
[~,idx] = min(abs(sc2(c_vec) - sc3(c_vec)));
c_sn    = c_vec(idx);              %#ok<NASGU> % not used later, but kept for completeness

% First c where sc2>0 (lower fold)
idx2 = find(sc2(c_vec)>0, 1, 'first');
c0   = c_vec(idx2);

%% Parameter sweep inside multistable region
bb = 0.2:0.9:1.2;                  % herbivory rates inside bistability
X0 = sc2(bb);                      % unstable mid‑branch as IC anchor

nRuns = 50;                       % ensemble size for statistics

sol = struct('b',[],'IC',[],'t',[],'x',[]);
sol = repmat(sol, numel(bb)*nRuns*2, 1);   % ×2 for ± perturbations
rec = 0;                           % record counter

for run = 1:nRuns
    rng(run)                       % new random seed each run
    for i = 1:numel(bb)
        b  = bb(i);
        xU = X0(i);                % mid equilibrium at this b
        for sgn = [-1 1]           % two close ICs bracketing xU
            rec = rec + 1;
            x0  = xU + sgn*1e-3;
            [t,x] = simulateFractional(al,@fun,t0,T,x0,h,b,sigma);

            sol(rec).b  = b;
            sol(rec).IC = x0;
            sol(rec).t  = t(:);
            sol(rec).x  = x(:);
        end
    end
end

%% Plots (first run only to avoid clutter)
figure; hold on
for i = 1:numel(bb)
    idx = ( (i-1)*2 + 1 );  % first two records correspond to run=1
    plot(sol(idx).t, sol(idx).x, '--');
    plot(sol(idx+1).t, sol(idx+1).x, '-');
end
xlabel('Time'); ylabel('x(t)'); box on

%% Save ensemble
save('herbivore_dynamics_memory.mat','sol')

% Flatten to a table for Python / R analysis
rows = sum(cellfun(@numel,{sol.t}));
tbl = table('Size',[rows 6], ...
    'VariableTypes',{'double','double','double','double','double','double'}, ...
    'VariableNames',{'run','b','ic','t','x','sigma'});
row = 0;
for k = 1:numel(sol)
    Nt = numel(sol(k).t);
    rng_id = ceil(k/(2*numel(bb)));    % back‑compute run index
    idx = row + (1:Nt);
    tbl.run(idx)  = rng_id;
    tbl.b(idx)    = sol(k).b;
    tbl.ic(idx)   = sol(k).IC;
    tbl.t(idx)    = sol(k).t;
    tbl.x(idx)    = sol(k).x;
    tbl.sigma(idx)= sigma;
    row = row + Nt;
end
writetable(tbl,'herbivore_dynamics_memory.csv');
fprintf('Ensemble saved to .mat and .csv (%d rows)\n', rows);

%% ---------------------------------------------------------------
%% Local functions (MATLAB R2016b+ supports local fns in scripts)
function [t,x] = simulateFractional(al,fun,t0,T,x0,h,b,sigma)
% Integrate a 1‑D Caputo fractional ODE with additive white noise
%   D^al x = f(x)    (Grünwald–Letnikov approximation)
%   x_{n+1} = x_n + h * fracDeriv + sigma*sqrt(h)*N(0,1)

    Nt = round((T - t0)/h) + 1;
    t  = linspace(t0,T,Nt)';
    x  = zeros(Nt,1); x(1) = x0;

    % GL weights w_k = (-1)^k * gamma(al+1)/(gamma(k+1)*gamma(al-k+1))
    M = Nt-1;                        % memory length (full past)
    w = ones(M,1);
    for k = 1:M-1
        w(k+1) = w(k)*( (k-al-1)/k );
    end
    w = w.*((-1).^(0:M-1)');        % include (-1)^k factor

    for n = 1:Nt-1
        kmax = min(n-1,M-1);        % **FIXED** ensure index >=1
        fracDeriv = 0;
        for k = 0:kmax
            fracDeriv = fracDeriv + w(k+1)*fun(x(n-k),b);
        end
        fracDeriv = fracDeriv / h^al;           % GL prefactor

        % Euler–Maruyama update
        x(n+1) = x(n) + h*fracDeriv + sigma*sqrt(h)*randn;
    end
end

function dx = fun(x,b)
% Herbivore model drift term (bistable)
    global r K A
    dx = r*x*(1 - x/K) - b*x/(x + A);
end
