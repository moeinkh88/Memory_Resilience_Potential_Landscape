function runHerbivorySimulation
% runHerbivorySimulation  Integrate the fractional herbivory model and
% visualise the effect of two pulse perturbations.
%
%   Requires the solver FDE_PI2_IM by Garrappa (MATLAB File-Exchange #64845).

%% MODEL & NUMERICS --------------------------------------------------------
p               = struct();          % model constants
p.r             = 0.8;
p.K             = 6.0;
p.A             = 0.2;
p.B0            = .4;               % baseline herbivory rate

% pulse perturbations ------------------------------------------------------
p.pulses = [ ...
    struct("delta",1,"tStart",10,"tEnd",20)   % first pulse (b -> b+0.12)
    struct("delta",1.5,"tStart",90,"tEnd",130)   % second pulse (b -> b+0.16)
    ];

% solver options -----------------------------------------------------------
alpha   = 0.8;           % fractional order
tSpan   = [0 200];       % integration window
h       = 0.01;          % time step
x0      = 5.460549503;        % initial state (stable equilibrium)

%% SOLVE -------------------------------------------------------------------
F  = @(t,x)  rhsHerbivory(t,x,p);
JF = @(t,x) jacHerbivory(t,x,p);

% FDE_PI2_IM expects the parameter vector *P*, so we pass an empty placeholder
[tt,xx] = FDE_PI2_IM(alpha,F,JF,tSpan(1),tSpan(2),x0,h,[]); %#ok<ASGLU>

%% PLOT --------------------------------------------------------------------
figure; hold on
ylims = [0 1.05*max(xx(:))];

% shade each pulse interval
for k = 1:numel(p.pulses)
    patch([p.pulses(k).tStart p.pulses(k).tEnd p.pulses(k).tEnd p.pulses(k).tStart], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [0.9290 0.6940 0.1250], 'FaceAlpha',0.15, 'EdgeColor','none');
end

plot(tt,xx,'LineWidth',1.5,'Color',[0.8500    0.3250    0.0980]);
xlabel('Time'); ylabel('State x(t)');
axis([tSpan ylims]); box on
title('Fractional herbivory model with two pulse perturbations');

%% SAVE RESULTS ------------------------------------------------------------
outDir = 'results';                      % folder keeps things tidy
if ~exist(outDir,'dir');  mkdir(outDir); end

% 1) Everything together (best for re-loading in MATLAB)
save(fullfile(outDir,'herbivoryTraj.mat'), ...
     'tt','xx','p');                     % includes the p.pulses info too

% 2) Trajectory to CSV for quick inspection elsewhere
trajTbl = table(tt(:),xx(:), ...
                'VariableNames',{'time','state'});
writetable(trajTbl, fullfile(outDir,'herbivoryTraj.csv'));

% 3) Perturbation windows to CSV (start & end times plus Δb)
pulseTbl = struct2table(p.pulses);
writetable(pulseTbl, fullfile(outDir,'pulseWindows.csv'));
fprintf('Results saved to folder "%s"\n',outDir);

end
%% ------------------------------------------------------------------------

% ----------  model RHS (d/dt)^α x = f(t,x)  ------------------------------
function f = rhsHerbivory(t,x,p)
b = currentB(t,p);
f = p.r*x.*(1 - x/p.K) - b*x./(x + p.A);
end

% ----------  Jacobian ∂f/∂x  ---------------------------------------------
function J = jacHerbivory(t,x,p)
b = currentB(t,p);
J = p.r*(1 - 2*x/p.K) - b./(x + p.A) + b*x./(x + p.A).^2;
end

% ----------  helper: piecewise herbivory rate  ---------------------------
function b = currentB(t,p)
b = p.B0;                       % baseline value
for k = 1:numel(p.pulses)
    if t >= p.pulses(k).tStart && t <= p.pulses(k).tEnd
        b = p.B0 + p.pulses(k).delta;
        break
    end
end
end
