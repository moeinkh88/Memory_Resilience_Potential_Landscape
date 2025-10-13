%% fit_internal_memory_integer.m
% Parameter estimation for the 2‑compartment "internal memory" herbivore model
%   dx/dt = r*x*(1 - x/K) - b*x/(z + A)
%   dz/dt = (x - z)/tau
%
% Data file format: herbivore_dynamics_memory6.csv
% columns: run, b, ic_id, ic_val, t, x

clc; clear; close all

dataFile = 'herbivore_dynamics_memory4.csv';
tbl      = readtable(dataFile);

%% ---------------- USER OPTIONS ------------------------------------------
fixed_rKA = false;      % true  → fit only the b‑values (and keep tau fixed)
                        % false → fit [r K A tau] + all b‑values

% initial guesses (tune as needed)
r0 = 0.8;  K0 = 6;  A0 = 0.1;
tau0 = 0.1;                    % time constant of the low‑pass compartment
b0 = 0.5;                      % initial guess for every b_i

% solver / optimiser tolerances
odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-12,...
                 'NonNegative',[1 2]);  % keep x,z ≥ 0
lsqOpts = optimoptions('lsqnonlin',...
                       'Display','iter','MaxFunctionEvaluations',4e3);
%% -----------------------------------------------------------------------

runIDs = unique(tbl.run);
nruns  = numel(runIDs);

% Build parameter vector:
%   p = [r K A tau  b1 … b_nruns]  (when fixed_rKA is false)
%   p = [               b1 … b_nruns] (when fixed_rKA is true)
if fixed_rKA
    p0 = repmat(b0,1,nruns);
    lb = zeros(1,nruns);
    ub = 10*ones(1,nruns);
else
    p0 = [r0 K0 A0 tau0 repmat(b0,1,nruns)];
    lb = [0  0  0  1e-10 zeros(1,nruns)];
    ub = [10 20 10  50   10*ones(1,nruns)];
end

% ---- main least‑squares fit -------------------------------------------
obj = @(p) allRunsResidual(p,tbl,runIDs,fixed_rKA,odeOpts);

[pEst,resnorm,~,exitflag] = lsqnonlin(obj,p0,lb,ub,lsqOpts);

% ---- unpack fitted parameters -----------------------------------------
if fixed_rKA
    rEst = r0; KEst = K0; AEst = A0; tauEst = tau0;
    bEst = pEst(:);
else
    rEst = pEst(1); KEst = pEst(2); AEst = pEst(3); tauEst = pEst(4);
    bEst = pEst(5:end).';
end

%% ---------------- results ----------------------------------------------
fprintf('\n=========  FIT COMPLETED  =========\n');
fprintf('exitflag: %d    resnorm: %.5g\n\n',exitflag,resnorm);
fprintf('r = %.6g\nK = %.6g\nA = %.6g\ntau = %.6g\n\n',rEst,KEst,AEst,tauEst);
for k = 1:nruns
    fprintf('b(%2d) for run %d :  %.6g\n',k,runIDs(k),bEst(k));
end

save('fit_internal_memory_results.mat','pEst','rEst','KEst','AEst','tauEst','bEst', ...
     'resnorm','exitflag');
fprintf('\nAll results saved in fit_internal_memory_results.mat\n');

%% =========== helper functions ==========================================
function err = allRunsResidual(p,tbl,runIDs,fixed_rKA,odeOpts)
% Stack residuals from every trajectory into a single column vector
    if fixed_rKA
        r = 0.8; K = 6; A = 0.1; tau = 1.0;   % <-- EDIT if your fixed values differ
        bAll = p(:);
    else
        r    = p(1);  K = p(2);  A = p(3);  tau = p(4);
        bAll = p(5:end);
    end

    err = [];

    for ir = 1:numel(runIDs)
        rows = tbl.run == runIDs(ir);
        tRun  = tbl.t(rows);
        xRun  = tbl.x(rows);
        icRun = tbl.ic_id(rows);
        icVal = tbl.ic_val(rows);

        for ic = [1 2]                 % each initial condition
            sel  = icRun == ic;
            t    = tRun(sel);
            xObs = xRun(sel);
            x0   = icVal(find(sel,1));         % scalar
            z0   = x0;                          % key choice: start filtered state at x0

            b = bAll(ir);                      % same b for both ICs of this run

            % 2‑state ODE
            ode = @(t,y) [ r*y(1).*(1 - y(1)./K) - b*y(1)./(y(2) + A);   % dx/dt
                           (y(1) - y(2))/tau ];                          % dz/dt

            xSim = ode45_x_atTimes(ode,t,x0,z0,odeOpts);

            if any(~isfinite(xSim))
                err = [err; 1e6*ones(size(xObs))]; %#ok<AGROW>  % robustify
            else
                err = [err; xSim - xObs];         %#ok<AGROW>
            end
        end
    end
end

function x = ode45_x_atTimes(odefun,tspan,x0,z0,odeOpts)
% Integrate the 2‑state system and return the x‑component at measurement times.
    y0  = [x0; z0];
    sol = ode45(odefun,[tspan(1) tspan(end)],y0,odeOpts);
    Y   = deval(sol,tspan).';
    x   = Y(:,1);
end

%% ---------------- visual check (data vs. fitted model) ------------------
figure; hold on
clr  = lines(nruns);
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
    ode = @(t,y) [ rEst*y(1).*(1 - y(1)./KEst) - b*y(1)./(y(2) + AEst);
                   (y(1) - y(2))/tauEst ];

    for ic = [1 2]
        sel  = icRun == ic;
        t    = tRun(sel);
        xObs = xRun(sel);
        x0   = icVal(find(sel,1));
        z0   = x0;

        % integrate once for this IC
        y0   = [x0; z0];
        sol  = ode45(ode,[t(1) t(end)],y0,odeOpts);
        Y    = deval(sol,t);
        xSim = Y(1,:).';

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
title('Two‑compartment herbivore model – fitted integer‑order trajectories vs. data');
legend('Location','bestoutside');
grid on; box on
