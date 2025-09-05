%% fit_theta_hill_herbivore.m
% Fits the θ‑logistic + Hill‑type grazing model
%      dX/dt = r*X*(1 - (X/K)^θ) - B*X^(q+1) / (A^q + X^q)
%
% The data must be in herbivore_dynamics_memory4.csv
% (columns: run, b, ic_id, ic_val, t, x).

clc; clear; close all

dataFile = 'herbivore_dynamics_memory2.csv';
tbl      = readtable(dataFile);

%% --------------- USER‑SET OPTIONS ---------------------------------------
% Initial guesses for the GLOBAL parameters
r0     = 0.8;
K0     = 6;
A0     = 0.1;
theta0 = 1;        % θ = 1  ⇒ standard logistic
q0     = 1;        % q = 1  ⇒ Holling‑II grazing

% Initial guess for every run‑specific B‑value
B0     = 0.5;

% Bounds (non‑negative, generous upper limits)
lbGlob = [0   0    0    0  0];   % r  K   A   θ   q
ubGlob = [100  200   50    50    50  ];

% Solver tolerances
odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);
lsqOpts = optimoptions('lsqnonlin',...
                       'Display','iter','MaxFunctionEvaluations',6e3);
%% -----------------------------------------------------------------------

runIDs = unique(tbl.run);
nruns  = numel(runIDs);

% Unknown vector   p = [ r K A θ q  B1 … B_nruns ]
p0 = [r0  K0  A0  theta0 q0 repmat(B0,1,nruns)];
lb = [lbGlob              zeros(1,nruns)];
ub = [ubGlob              10*ones(1,nruns)];

% -------------- least‑squares fit ---------------------------------------
obj = @(p) residualThetaHill(p,tbl,runIDs,odeOpts);

[pEst,resnorm,~,exitflag] = lsqnonlin(obj,p0,lb,ub,lsqOpts);

% -------------- unpack parameters ---------------------------------------
rEst     = pEst(1);
KEst     = pEst(2);
AEst     = pEst(3);
thetaEst = pEst(4);
qEst     = pEst(5);
BEst     = pEst(6:end).';

%% -------------- report --------------------------------------------------
fprintf('\n=========  θ‑HILL FIT COMPLETED  =========\n');
fprintf('exitflag: %d   resnorm: %.5g\n\n',exitflag,resnorm);
fprintf('r      = %.6g\n',rEst);
fprintf('K      = %.6g\n',KEst);
fprintf('A      = %.6g\n',AEst);
fprintf('theta  = %.6g\n',thetaEst);
fprintf('q      = %.6g\n\n',qEst);
for k = 1:nruns
    fprintf('B(%2d) for run %d :  %.6g\n',k,runIDs(k),BEst(k));
end

save('fit_theta_hill_results.mat','pEst','rEst','KEst','AEst', ...
     'thetaEst','qEst','BEst','resnorm','exitflag');
fprintf('\nAll results saved in fit_theta_hill_results.mat\n');

%% ---------------- VISUAL CHECK: data vs. fitted model -------------------
figure; hold on
clr  = lines(nruns);
mkr  = {'o','s'};     % IC‑1 circle, IC‑2 square
lsty = {'--','-'};    % IC‑1 dashed, IC‑2 solid

for ir = 1:nruns
    rows = tbl.run == runIDs(ir);
    tRun  = tbl.t(rows);
    xRun  = tbl.x(rows);
    icRun = tbl.ic_id(rows);
    icVal = tbl.ic_val(rows);
    B     = BEst(ir);

    ode = @(t,x) rEst*x.*(1 - (x./KEst).^thetaEst) ...
               - B*x.^(qEst+1) ./ (AEst.^qEst + x.^qEst);

    for ic = [1 2]
        sel  = icRun == ic;
        t    = tRun(sel);
        xObs = xRun(sel);
        x0   = icVal(find(sel,1));

        sol  = ode45(ode,[t(1) t(end)],x0,odeOpts);
        xSim = deval(sol,t);

        plot(t,xObs, mkr{ic}, ...
            'Color',clr(ir,:), 'MarkerFaceColor',clr(ir,:), ...
            'DisplayName',sprintf('data  run %d  IC %d',runIDs(ir),ic));
        plot(t,xSim, lsty{ic}, ...
            'Color',clr(ir,:), 'LineWidth',1.5, ...
            'DisplayName',sprintf('fit   run %d  IC %d',runIDs(ir),ic));
    end
end
xlabel('time');  ylabel('x(t)');
title('θ‑logistic + generalised grazing – fitted vs. data');
legend('Location','bestoutside');
grid on;  box on

%% ==================== helper function ==================================
function err = residualThetaHill(p,tbl,runIDs,odeOpts)
% Builds a column vector of residuals for every observed point
    r     = p(1);    K = p(2);    A = p(3);
    theta = p(4);    q = p(5);
    BAll  = p(6:end);

    err = [];

    for ir = 1:numel(runIDs)
        rows = tbl.run == runIDs(ir);
        tRun  = tbl.t(rows);
        xRun  = tbl.x(rows);
        icRun = tbl.ic_id(rows);
        icVal = tbl.ic_val(rows);

        B = BAll(ir);       % this run’s grazing rate

        ode = @(t,x) r*x.*(1 - (x./K).^theta) ...
                   - B*x.^(q+1) ./ (A.^q + x.^q);

        for ic = [1 2]
            sel  = icRun == ic;
            t    = tRun(sel);
            xObs = xRun(sel);
            x0   = icVal(find(sel,1));

            sol  = ode45(ode,[t(1) t(end)],x0,odeOpts);
            xSim = deval(sol,t);

            err  = [err; xSim - xObs];      %#ok<AGROW>
        end
    end
end
