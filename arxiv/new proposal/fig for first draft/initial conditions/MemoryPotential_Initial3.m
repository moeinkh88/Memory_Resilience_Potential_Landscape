% dynamics_potential_memory_customIC.m
% --------------------------------------------------------------
% Plots system dynamics and potential energy for the *memory*
% case (fractional order alpha = 0.8).
% Besides the three canonical ICs (MLX0, MX0, RX0) we now add
% TWO user‑specified initial points **with their own start times**.
% The resulting trajectories launch at those absolute times, and
% the potential landscape shows markers at the same states.
% --------------------------------------------------------------

clc; clear;

%% Fractional derivative order (memory case only)
al  = 0.8;             % memory

%% Simulation grid
xt = 0:1e-5:9;         % fine grid for potential evaluation
F  = @funPoly_den;      % herbivory function
JF = @JfunPoly_den;     % Jacobian of F

%% Model parameters for *memory* dynamics
r = 0.8;   K = 6;        % intrinsic growth & carrying capacity
A = 0.2;   B = 0.4;      % herbivory parameters

%% Misc. numerical settings
h      = 0.01;          % time‑step for the PI2 solver
Tfinal = 50;           % integration horizon (absolute maximum)

%% --- Coefficients of the cubic polynomial ---
a1  = A*r - B;
a2  = r*(1 - A/K);
a3  = -r/K;
par = [a3 a2 a1 A];

%% --- Potential landscape of F(x) ---
U   = -polyint(par(1:3));      % -∫F(x) dx
U1  = polyval(U,xt);
TF  = islocalmin(U1);          % stable equilibria
TF1 = islocalmax(U1);          % unstable threshold
xmin = xt(TF);                 % valley bottom(s)
xmax = xt(TF1);                % ridge crest

%% --- Canonical initial states (unchanged) ---
X0   = 0;                % left valley (not plotted)
MLX0 = xmax - 0.002;     % just left of ridge
MX0  = xmax + 0.002;     % just right of ridge

% Approximate right‑hand valley as in the original script
idxMin = find(TF,1,'first');
ind    = 0; epsi2 = 1e-5;
while ~ind
    Indx2 = find( abs(U1(idxMin:end) - U1(TF1)) < epsi2, 1 );
    if isempty(Indx2)
        epsi2 = epsi2*2;
    else
        ind = 1;
    end
end
RX0 = xt(Indx2 + idxMin - 1)+1;

%% -----------------------------------------------------------
%% USER‑SPECIFIED EXTRA INITIAL STATES **WITH** START TIMES
%% -----------------------------------------------------------
extraICs   = [7 , 3];             % states
extraT0    = [0.9060023  , 19.22974022]; % corresponding launch times
assert(numel(extraICs)==numel(extraT0),"extraICs and extraT0 must match");

%% -----------------------------------------------------------
%% Integrate trajectories
%% -----------------------------------------------------------
ICs      = [MLX0, MX0, RX0, extraICs];   % the left‑valley X0 excluded to reproduce original plot
nTraj    = numel(ICs);
Toffsets = [0, 0, 0, extraT0];           % first three start at t=0

Traj  = cell(nTraj,1);
Tgrid = cell(nTraj,1);

for k = 1:nTraj
    t0   = Toffsets(k);
    tEnd = Tfinal - t0;  % integrate only up to Tfinal overall
    [t,x] = FDE_PI2_IM(al,F,JF,0,tEnd,ICs(k),h,par);
    Traj{k}  = x;
    Tgrid{k} = t + t0;   % shift to absolute time
end

%% --------------------------------------------------------------
%% Potential reconstruction (unchanged mathematics)
%% --------------------------------------------------------------
% We still use the four canonical ICs (including X0) for potential curve
[t_X0,X]   = FDE_PI2_IM(al,F,JF,0,Tfinal,X0,h,par);
[t_MLX,MLX] = FDE_PI2_IM(al,F,JF,0,Tfinal,MLX0,h,par);
[t_MX,MX]   = FDE_PI2_IM(al,F,JF,0,Tfinal,MX0,h,par);
[t_RX,RX]   = FDE_PI2_IM(al,F,JF,0,Tfinal,RX0,h,par);

dx   = diff(X')./h;
dxML = diff(MLX')./h;
dxM  = diff(MX')./h;
dxR  = diff(RX')./h;

J = 3:length(MX)-2;

dx(J-1)   = 1/(12*h).*(X(J-2)'   - 8.*X(J-1)'   + 8.*X(J+1)'   - X(J+2)');
dxML(J-1) = 1/(12*h).*(MLX(J-2)' - 8.*MLX(J-1)' + 8.*MLX(J+1)' - MLX(J+2)');
dxM(J-1)  = 1/(12*h).*(MX(J-2)'  - 8.*MX(J-1)'  + 8.*MX(J+1)'  - MX(J+2)');
dxR(J-1)  = 1/(12*h).*(RX(J-2)'  - 8.*RX(J-1)'  + 8.*RX(J+1)'  - RX(J+2)');

DX1   = cat(1,dx,flip(dxML),dxM,flip(dxR));
Xall1 = cat(1,X(1:end-1)',flip(MLX(1:end-1)'),MX(1:end-1)',flip(RX(1:end-1)'));

NX   = round((RX0-0)/(h/2));
[C,ia] = unique(Xall1);
Xall  = linspace(0,RX0,NX);
DX    = interp1(C,DX1(ia),Xall,'pchip');

Q = -cumtrapz(Xall,DX);     % potential values

% Potential values exactly at extraICs
Qextra = interp1(Xall,Q,extraICs,'pchip');

%% --------------------------------------------------------------
%% Plotting
%% --------------------------------------------------------------
figure('Color','w');

% (1) Dynamics --------------------------------------------------
subplot(1,2,1); hold on; box on;
col  = lines(nTraj);
for k = 1:nTraj
    plot(Tgrid{k},Traj{k},'LineWidth',1.6,'Color',col(k,:));
end
% Mark the starting dots for the two new ICs
scatter(extraT0,extraICs,60,col(4:5,:),"filled");

xlabel('Time'); ylabel('States');
legendStrings = ["IC near X_{u}^-","IC near X_{u}^+","IC right valley", ...
                 "IC = 5.9 (t_0=3.79)","IC = 2 (t_0=17.49)"];
legend(legendStrings,'Location','best');
axis tight;
title(['System dynamics (\alpha = ',num2str(al),')']);

% (2) Potential landscape --------------------------------------
subplot(1,2,2); hold on; box on;
plot(Xall,Q,'LineWidth',2,'Color',[0 0.4470 0.7410]);
% scatter(extraICs,Qextra,60,col(4:5,:),"filled");

xlabel('States'); ylabel('Potential');
axis tight;
title(['Potential landscape (',char(945),' = ',num2str(al),')']);


%% --------------------------------------------------------------
%% Local functions (unchanged) ----------------------------------
%% --------------------------------------------------------------
function dx = funPoly_den(~,x,par)
    a1=par(3); a2=par(2); a3=par(1); A=par(4);
    dx = (a3*x^3 + a2*x^2 + a1*x) / (x + A);
end

function dx = JfunPoly_den(~,x,par)
    a1=par(3); a2=par(2); a3=par(1); A=par(4);
    num = 2*a3*x.^3 + (3*a3*A + a2).*x.^2 + 2*a2*A*x + a1*A;
    den = (A + x).^2;
    dx  = num ./ den;
end
