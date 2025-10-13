function fitPulsesToData
% fitPulsesToData  Estimate the two pulse strengths Δb1 and Δb2
%                  for the integer-order herbivory model so that it
%                  matches a target trajectory (with known pulse timing).

%% ------------------------------------------------------------------------
%  USER INPUTS
% -------------------------------------------------------------------------
dataFile   = 'herbivoryTraj.csv';   % 2-column CSV produced earlier
pulseFile  = 'pulseWindows.csv';    % gives tStart, tEnd for each pulse
                                     
% Fixed integer-model parameters (from your earlier fit)
par.r  = 0.57261;
par.K  = 5.82416;
par.A  = 0.111232;
par.b0 = 0.237859;

%% ------------------------------------------------------------------------
%  LOAD DATA
% -------------------------------------------------------------------------
tbl      = readtable(dataFile);
tData    = tbl.time;
xData    = tbl.state;
x0       = xData(1);        % use first data point as initial condition

pulseTbl = readtable(pulseFile);
pulseTimes = [pulseTbl.tStart, pulseTbl.tEnd];   % N×2 matrix
if height(pulseTbl)~=2
    error('Expected exactly two pulses in pulseWindows.csv');
end

% ------------------------------------------------------------------
%  OPTIMISATION SETUP   (replace the old block with this one)
% ------------------------------------------------------------------
p0   = [1 1.5];          % initial guess  [Δb1 Δb2]
lb   = [0   0];            % lower bounds   (≥0)
ub   = [5   5];            % upper bounds   (arbitrary, change if needed)

% Define a handle that “freezes” the known inputs
modelResidual = @(p) residual(p,par,pulseTimes,tData,xData);  % <-- key line


% ---------- solver options ----------
lsqOpts = optimoptions('lsqnonlin', ...
        'Algorithm','trust-region-reflective', ...
        'Display','iter-detailed', ...
        'FiniteDifferenceType','central', ...
        'DiffMinChange',1e-18, ...     % smaller step → better derivative accuracy
        'StepTolerance',1e-5, ...
        'FunctionTolerance',1e-14, ...
        'OptimalityTolerance',1e-14, ...
        'MaxIterations',500, ...
        'ScaleProblem','jacobian');
% ---------- ODE solver options (passed via global) ----------
odeOpts = odeset('RelTol',1e-14,'AbsTol',1e-14,'NonNegative',1);
setappdata(0,'odeOpts',odeOpts);       % stash for residual() to pick up

% % ---------- run several starts to dodge local minima ----------
% bestCost = inf;   bestPulse = p0;
% for k = 1:5
%     if k>1, p0 = lb + (ub-lb).*rand(1,2); end   % random restart
%     [pEst,~,~,exitflag,output,residualVec] = lsqnonlin(modelResidual,p0,lb,ub,lsqOpts);
%     thisCost = norm(residualVec)^2;
%     if thisCost < bestCost
%         bestCost  = thisCost;
%         bestPulse = pEst;
%     end
% end
% pulseEst = bestPulse;

% fprintf('\nEstimated pulses (high-accuracy fit):  Δb1 = %.8f   Δb2 = %.8f\n', ...
%         pulseEst(1), pulseEst(2));

haveLSQ = license('test','optimization_toolbox') && exist('lsqnonlin','file');

if haveLSQ
    fprintf('Using lsqnonlin (bounded least squares)…\n');
    opts      = optimoptions('lsqnonlin','Display','iter','TolFun',1e-12);
    pulseEst  = lsqnonlin(modelResidual,p0,lb,ub,opts);
else
    fprintf('Using fminsearch (no bounds)…\n');
    pulseEst  = fminsearch(@(p) sum(modelResidual(p).^-2),p0);
end


%% ------------------------------------------------------------------------
%  DISPLAY FIT
% -------------------------------------------------------------------------
xSim = simulateHerbivory(par,pulseTimes,pulseEst,tData,x0);

figure; hold on
plot(tData,xData,'color',[0.8500    0.3250    0.0980],'LineWidth',2,'DisplayName','Data');
plot(tData,xSim ,'--','color',[0    0.4470    0.7410],'LineWidth',2,'DisplayName','Integer fit');
xlabel('Time'); ylabel('State'); box on; legend
title('Integer-order model fitted to memory-model trajectory');
% -----------------------------------------------
% after the optimisation block, add:
fprintf('\nEstimated pulses:  Δb1 = %.6f   Δb2 = %.6f\n', ...
        pulseEst(1), pulseEst(2));
% -----------------------------------------------

end
% ========================================================================

% ----------  residuals for least-squares  -------------------------------
function res = residual(pulses,par,pulseTimes,tData,xData,x0)
if nargin<6, x0 = xData(1); end
xSim = simulateHerbivory(par,pulseTimes,pulses,tData,x0);
res  = xSim - xData;        % vector of errors
end

% ----------  simulate integer-order herbivory model  ---------------------
function xOut = simulateHerbivory(par,pulseTimes,pulses,tEval,x0)
% piece-wise RHS with pulses
    function dx = rhs(t,x)
        b = par.b0;
        if t>=pulseTimes(1,1) && t<=pulseTimes(1,2)
            b = b + pulses(1);        % first pulse
        elseif t>=pulseTimes(2,1) && t<=pulseTimes(2,2)
            b = b + pulses(2);        % second pulse
        end
        dx = par.r*x*(1 - x/par.K) - b*x/(x + par.A);
    end

% Integrate (ode45 is fine for integer ODE)
odeOpts = getappdata(0,'odeOpts');      % retrieve high-precision opts
[~,X]   = ode15s(@rhs, tEval, x0, odeOpts);

xOut  = X;     % same length & ordering as tEval
end
