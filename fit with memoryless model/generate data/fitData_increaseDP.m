%% fit_herbivore_integer.m
% Parameter estimation for the integer‑order herbivore model
%   dx/dt = r*x*(1 - x/K) - b*x/(x + A)
%
% The data file must be in the format produced by herbivore_dynamics_memory4.csv
% (columns: run, b, ic_id, ic_val, t, x).

% when memory increase the potential landscape

clc; clear; close all

dataFile = 'herbivore_dynamics_memory7.csv';
tbl      = readtable(dataFile);

%% ---------------- USER OPTIONS ------------------------------------------
fixed_rKA = false;      % true  → fit only the b‑values
                       % false → fit [r K A] + all the b‑values

% initial guesses (reasonable defaults for your data)
r0 = 80;  K0 = 3;  A0 = 0.2;
b0 = 50;                 % initial guess for every b_i

% solver / optimiser tolerances
odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-12);
lsqOpts = optimoptions('lsqnonlin',...
                       'Display','iter','MaxFunctionEvaluations',4e3);
%% -----------------------------------------------------------------------

runIDs = unique(tbl.run);
nruns  = numel(runIDs);

% Build the vector of unknowns:
%   p = [r K A b1 … b_nruns]  (when fixed_rKA is false)
%   p = [       b1 … b_nruns]  (when fixed_rKA is true)
if fixed_rKA
    p0 = repmat(b0,1,nruns);
    lb = zeros(1,nruns);                % positivity constraints
    ub = 10*ones(1,nruns);
else
    p0 = [r0 K0 A0 repmat(b0,1,nruns)];
    lb = [6  0  0   1*ones(1,nruns)];
    ub = [200 40  4  190*ones(1,nruns)];
end

% ---- main least‑squares fit -------------------------------------------
obj = @(p) allRunsResidual(p,tbl,runIDs,fixed_rKA,odeOpts);

[pEst,resnorm,~,exitflag] = lsqnonlin(obj,p0,lb,ub,lsqOpts);

% ---- unpack fitted parameters -----------------------------------------
if fixed_rKA
    rEst = r0; KEst = K0; AEst = A0;
    bEst = pEst(:);
else
    rEst = pEst(1);    KEst = pEst(2);    AEst = pEst(3);
    bEst = pEst(4:end).';
end

%% ---------------- results ----------------------------------------------
fprintf('\n=========  FIT COMPLETED  =========\n');
fprintf('exitflag: %d    resnorm: %.5g\n\n',exitflag,resnorm);
fprintf('r = %.6g\nK = %.6g\nA = %.6g\n\n',rEst,KEst,AEst);
for k = 1:nruns
    fprintf('b(%2d) for run %d :  %.6g\n',k,runIDs(k),bEst(k));
end

save('fit_integer_results_increaseDP.mat','pEst','rEst','KEst','AEst','bEst', ...
     'resnorm','exitflag');
fprintf('\nAll results saved in fit_integer_results.mat\n');

%% =========== helper functions ==========================================
function err = allRunsResidual(p,tbl,runIDs,fixed_rKA,odeOpts)
% Stack residuals from every trajectory into a single column vector
    if fixed_rKA
        r = 80; K = 3; A = 0.2;      % <-- EDIT here if your fixed values differ
        bAll = p(:);
    else
        r    = p(1);  K = p(2);  A = p(3);
        bAll = p(4:end);
    end

    err = [];

    for ir = 1:numel(runIDs)
        rows = tbl.run == runIDs(ir);
        tRun  = tbl.t(rows);
        xRun  = tbl.x(rows);
        icRun = tbl.ic_id(rows);
        icVal = tbl.ic_val(rows);

        for ic = [1 2]                         % each initial condition
            sel  = icRun == ic;
            t    = tRun(sel);
            xObs = xRun(sel);
            x0   = icVal(find(sel,1));         % scalar

            b = bAll(ir);                      % same b for both ICs of this run
            ode = @(t,x) r*x.*(1 - x./K) - b*x./(x + A);

            xSim = ode45_atTimes(ode,t,x0,odeOpts);

            err  = [err; xSim - xObs];         %#ok<AGROW>
        end
    end
end

function x = ode45_atTimes(odefun,tspan,x0,odeOpts)
% Integrate once over [t0 tf] and evaluate the solution at the
% exact measurement times using DEVAL — efficient and accurate.
    sol = ode45(odefun,[tspan(1) tspan(end)],x0,odeOpts);
    x   = deval(sol,tspan).';
end


%% ---------------- visual check (data vs. fitted model) ------------------
% We already have everything we need in memory:
%   • tbl         – the original data table
%   • rEst,K,ESt  – fitted global parameters
%   • bEst        – fitted b‑values (one per run)
%   • runIDs      – list of run indices
%   • nruns       – number of runs
%   • odeOpts     – solver tolerances

%%
figure; hold on
clr  = lines(nruns);               % colour per run
mkr  = {'o','s'};                  % IC‑1 → circle, IC‑2 → square
lsty = {'--','-'};                 % IC‑1 dashed, IC‑2 solid

for ir = 1:nruns
    rows = tbl.run == runIDs(ir);
    tRun  = tbl.t(rows);
    xRun  = tbl.x(rows);
    icRun = tbl.ic_id(rows);
    icVal = tbl.ic_val(rows);
    b     = bEst(ir);

    % ODE with fitted parameters for this run
    ode = @(t,x) rEst*x.*(1 - x./KEst) - b*x./(x + AEst);

    for ic = [1 2]
        sel  = icRun == ic;
        t    = tRun(sel);
        xObs = xRun(sel);
        x0   = icVal(find(sel,1));

        % integrate once for this IC
        sol  = ode45(ode,[t(1) t(end)],x0,odeOpts);
        xSim = deval(sol,t);

        % --- plot ---
        plot(t,xObs, mkr{ic}, ...
             'Color',clr(ir,:), 'MarkerFaceColor',clr(ir,:), ...
             'DisplayName',sprintf('data  run %d  IC %d',runIDs(ir),ic));
        plot(t,xSim, lsty{ic}, ...
             'Color',clr(ir,:), 'LineWidth',1.5, ...
             'DisplayName',sprintf('fit   run %d  IC %d',runIDs(ir),ic));
    end
end
xlabel('time');  ylabel('x(t)');
title('Herbivore model – fitted integer‑order trajectories vs. data');
legend('Location','bestoutside');
grid on;  box on
