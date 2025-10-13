% dynamics_potential_memory_multiIC.m
% --------------------------------------------------------------
% Plots system dynamics and potential energy for the *memory*
% case (fractional order alpha = 0.8) for MULTIPLE initial
% conditions.  Start by editing the "extraICs" vector and the
% script will do the rest (integration, plotting, legends).
% --------------------------------------------------------------

clc; clear;

%% Fractional derivative order (memory case only)
al  = 0.8;             % memory

%% Simulation grid
tx = 0:1e-5:9;         % fine grid for potential evaluation
F  = @funPoly_den;      % herbivory function
JF = @JfunPoly_den;     % Jacobian of F

%% Model parameters for *memory* dynamics
r = 0.8;   K = 6;        % intrinsic growth & carrying capacity
A = 0.2;   B = 0.4;      % herbivory parameters

%% Misc. numerical settings
h      = 0.01;          % time‑step for the PI2 solver
Tfinal = 300;           % integration horizon

%% Coefficients of the cubic polynomial (taken outside loop)
a1  = A*r - B;
a2  = r*(1 - A/K);
a3  = -r/K;
par = [a3 a2 a1 A];

%% Potential landscape of F(x) (unchanged)
U   = -polyint(par(1:3));      % -∫F(x) dx
U1  = polyval(U,tx);
TF  = islocalmin(U1);          % stable equilibria
TF1 = islocalmax(U1);          % unstable threshold
xmin = tx(TF);                 % valley bottom(s)
xmax = tx(TF1);                % ridge crest

%% --- Representative initial states from original script ---
X0   = 0;                % left valley
MLX0 = xmax - 0.002;     % just left of ridge
MX0  = xmax + 0.002;     % just right of ridge

% Approximate right‑hand valley as before
idxMin = find(TF,1,'first');
indx=find(TF==1);
ind=0; epsi1=0.0001;epsi2=0.00001;
while ind==0

Indx2=find((abs(U1(indx:end)-U1(TF1))<epsi2));

if isempty(Indx2)==1
    epsi2=epsi2*2;
else
    ind=1;
end
end
RX0    = tx(Indx2);

%% -----------------------------------------------------------
%% ADD **EXTRA** INITIAL CONDITIONS BELOW (edit only this line)
%% -----------------------------------------------------------
extraICs = [1.0 2.5 4.0 5.5];   % <--- your additional ICs
%% -----------------------------------------------------------

ICs = [X0, MLX0, MX0, RX0, extraICs];   % concatenate all ICs
nIC = numel(ICs);

%% Containers
Traj  = cell(nIC,1);
Tgrid = cell(nIC,1);

%% Integrate the FODE from every IC
for k = 1:nIC
    [t,x] = FDE_PI2_IM(al,F,JF,0,Tfinal,ICs(k),h,par);
    Traj{k}  = x;
    Tgrid{k} = t;
end

%% --------------------------------------------------------------
%% Potential reconstruction (based on original four ICs)
%% --------------------------------------------------------------
X   = Traj{1};  % X0
MLX = Traj{2};  % MLX0
MX  = Traj{3};  % MX0
RX  = Traj{4};  % RX0

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
Xall = linspace(0,RX0,NX);
DX   = interp1(C,DX1(ia),Xall,'pchip');

Q = -cumtrapz(Xall,DX);     % potential values

%% --------------------------------------------------------------
%% Plotting
%% --------------------------------------------------------------
figure('Color','w');

% (1) Dynamics
subplot(1,2,1); hold on;
pal = lines(nIC);  % distinct colours
for k = 1:nIC
    plot(Tgrid{k},Traj{k},'LineWidth',1.6,'Color',pal(k,:));
end
xlabel('Time'); ylabel('States');
legend( compose('IC = %.3g',ICs) ,'Location','best');
title(['System dynamics (\alpha = ',num2str(al),')']);
box on; axis tight;

% (2) Potential landscape
subplot(1,2,2); hold on;
plot(Xall,Q,'LineWidth',2,'Color',[0 0.4470 0.7410]);
xlabel('States'); ylabel('Potential');
box on; axis tight;

title(['Potential landscape (',char(945),' = ',num2str(al),')']);

%% --------------------------------------------------------------
%% Local functions
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
