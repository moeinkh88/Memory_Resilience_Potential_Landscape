% This code provides resistance
clear
clc

% Define the path to the data folder
dataPath = fullfile('..','..', 'data');

% Load all the .mat files from the data folder
load(fullfile(dataPath, 'PAR.mat')); % parameters for polynomial functions
load(fullfile(dataPath, 'CurvR.mat'));
load(fullfile(dataPath, 'Rdp.mat'));
load(fullfile(dataPath, 'xmax.mat')); % hilltop in potential landscapes (unstable points) of the models
load(fullfile(dataPath, 'xmin.mat')); % bottom of the valley in potential landscapes (positive stable points) of the models
load(fullfile(dataPath, 'Flat.mat'));
PAR(:,4)=[];

%% functions and conditions

F=@funPoly; % polynomial functions
JF=@JfunPoly; % Jacobian of functions

Order=1:-.2:.8; % order of derivatives
t0=0; % initial time
T=50; % final time
h=0.01; % step size for computation

J=1;
%%

for i=1:length(Order)
    for j=1:length(PAR)
tic


    p1=0.01; p2=10; %initial guess to perturbations later we update them and find the minimum perturbation that can shift the system
    pt0=10;pt1=20; % the time when perturbation start and ends

%  Find minimum pulse size that causes the shift
% ----------------------------------------------------------------
tol      = 5e-4;      % required accuracy  (0.0005)
maxIter  = 30;        % hard stop for safety
pLow     = p1;        % known NO-shift amplitude
pHigh    = p2;        % known shift amplitude
pCrit    = NaN;       % will hold the critical value

% helper that returns TRUE if the trajectory ends in the zero basin
isShifted = @(xEnd,xUnstable)  xEnd < xUnstable - 1e-8;   % small safety margin

% Make sure the initial bracket really encloses the threshold
% (grow pHigh until it does, shrink pLow if necessary)
while true
    % --- test lower end -----------------------------------------
    [~,xLow] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),pLow,pt0,pt1]);
    noShift  = ~isShifted(xLow(end),xmax(j));
    % --- test upper end -----------------------------------------
    [~,xHigh] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),pHigh,pt0,pt1]);
    madeShift =  isShifted(xHigh(end),xmax(j));

    if  noShift && madeShift
        break                                   % good bracket
    elseif ~madeShift                          % pHigh too small → enlarge
        pHigh = 2*pHigh;
    else                                       % pLow already shifts → shrink
        pLow  = 0.5*pLow;
    end
end

% ---------- bisection --------------------------------------------
for k = 1:maxIter
    pMid = 0.5*(pLow + pHigh);

    [~,xMid] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),pMid,pt0,pt1]);

    if isShifted(xMid(end),xmax(j))            % still tips → move left
        pHigh = pMid;
    else                                       % still recovers → move right
        pLow  = pMid;
    end

    if (pHigh - pLow) <= tol
        pCrit = pHigh;                         % reached desired accuracy
        break
    end
end

% ------ record result ---------------------------------------------
CriticalPulse(j,i) = pCrit;   % pre-allocate this matrix outside the loops
fprintf('Order %.2f   model %4d   pCrit = %.4f\n',Order(i),j,pCrit);

toc
    end
end
%%
%% ---------------------------------------------------------------
%%  SECOND PASS  –– verify / update column 2   (Order = 0.2)
%% ---------------------------------------------------------------
Tlong   = 150;           % longer integration time
pt0     = 10;            % pulse start
pt1     = 20;            % pulse end
tol     = 5e-4;          % accuracy for new search
maxIter = 30;

col      = 2;            % <<— column that corresponds to α = 0.2
alpha    = Order(col);   % (just for readability)
nModels  = length(PAR);

isShifted = @(xEnd,xUnst)  xEnd < xUnst - 1e-8;

% storage for updated pulses (start with first-pass values)
CriticalPulse_long      = CriticalPulse(:,col);
modelsNeedingAdjustment = false(nModels,1);

fprintf('\n===== second pass (α = %.2f,  T = %d) =====\n',alpha,Tlong);

for j = 1:nModels
    pTest = CriticalPulse_long(j);
    if isnan(pTest) || pTest==0
        modelsNeedingAdjustment(j) = true;
        continue
    end

    [~,x] = FDE_PI2_IM(alpha,F,JF,0,Tlong,xmin(j),h,[PAR(j,:) pTest pt0 pt1]);
    if ~isShifted(x(end),xmax(j))         % still recovers → needs update
        modelsNeedingAdjustment(j) = true;
    end
end

%% ---------- re-search only the troublesome models ----------------
idxBad = find(modelsNeedingAdjustment);

for jj = idxBad'
    % 1) start bracket around previous pCrit
    pLow  = max(CriticalPulse_long(jj)*0.5, 1e-4);
    pHigh = max(CriticalPulse_long(jj)*2  , 1e-3);

    % make sure upper end really shifts
    while true
        [~,xH] = FDE_PI2_IM(alpha,F,JF,0,Tlong,xmin(jj),h,[PAR(jj,:) pHigh pt0 pt1]);
        if isShifted(xH(end),xmax(jj)), break, end
        pHigh = 2*pHigh;                % escalate if still no shift
        if pHigh > 1e6                  % absolute ceiling – give up
            warning('Model %d never tipped up to p = %g',jj,pHigh);
            pHigh = NaN;  break
        end
    end
    if isnan(pHigh),  CriticalPulse_long(jj) = NaN;  continue, end

    % 2) bisection with longer horizon
    for k = 1:maxIter
        pMid      = 0.5*(pLow + pHigh);
        [~,xMid]  = FDE_PI2_IM(alpha,F,JF,0,Tlong,xmin(jj),h,[PAR(jj,:) pMid pt0 pt1]);
        if isShifted(xMid(end),xmax(jj))
            pHigh = pMid;
        else
            pLow  = pMid;
        end
        if (pHigh - pLow) <= tol, break, end
    end
    CriticalPulse_long(jj) = pHigh;
    fprintf('model %4d  updated pCrit = %.4f\n',jj,pHigh);
end

%% ---------------------------------------------------------------
%%  RESULTS
%% ---------------------------------------------------------------
save('CriticalPulse_long.mat','CriticalPulse_long','idxBad')

%% plot
dt=xmin-xmax;
EcoRL=[CriticalPulse(:,1), CriticalPulse_long]';

%colors
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];
figure

hold on

indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);

indP1 = repmat({'Memory increases \newline potential depth'},length(indP),1);
indP2 = repmat({'Memory decreases \newline potential depth'},length(indN),1);

% abs(Rdp(indP,2)-Rdp(indP,1))./Rdp(indP,1);
EcoRdp1=(EcoRL(2,indP)-EcoRL(1,indP))./(EcoRL(2,indP)+EcoRL(1,indP));
% scatter(abs(Rdp(indN,2)-Rdp(indN,1))./Rdp(indN,1),(EcoRL(2,indN)-EcoRL(1,indN))./EcoRL(1,indN),[],Cgr)
EcoRdp2=(EcoRL(2,indN)-EcoRL(1,indN))./(EcoRL(2,indN)+EcoRL(1,indN));

v=violinplot([EcoRdp1';EcoRdp2'],[indP1;indP2],...
'Width',0.5,'Bandwidth',0.02,'ViolinColor',Cgr,'ViolinAlpha',0.004,...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);

% v.ShowData=true; v.ViolinAlpha=0.1;
 h=v.ScatterPlot; h.MarkerFaceAlpha=1; h.MarkerFaceColor='flat';
 
 h.CData=repmat(h.CData,numel(h.YData),1);
 index=h.YData<3; h.CData(index,1)=Cpr(1); h.CData(index,2)=Cpr(2); h.CData(index,3)=Cpr(3);

ax = gca;
ax.YAxis.Scale ="log";

axis tight

ylabel('Relative effect of memory on ecological resilience')
% xlabel('Memory effects on potential depth')


%% Eco----Curv and pd
figure
p1=scatter(Flat(:,1),Rdp(:,1),60,(EcoRL(1,:)),'filled');

% ax = gca;
% ax.XDir = 'reverse';

xlabel('Flatness')
ylabel('Potential depth')

set(gca,'ColorScale','log') 
% ax.YAxis.Scale ="log";
% ax.XAxis.Scale ="log";
p1.MarkerEdgeColor=[.7,.7,.7];

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Engineering resilience';

set(gca,'CLim',[0.0164   2.7019]);

colormap(flipud(gray))