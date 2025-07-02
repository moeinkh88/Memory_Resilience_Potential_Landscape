% To provide resilience, first we need to find the minimum perturbation (resistance) to
% use less perturbation, ensuring the recovery happens
clear
clc

% Define the path to the data folder
dataPath = fullfile('..','..', 'data');

% Load all the .mat files from the data folder
load(fullfile(dataPath, 'PAR.mat'));
load(fullfile(dataPath, 'CurvR.mat'));
load(fullfile(dataPath, 'Rdp.mat'));
load(fullfile(dataPath, 'xmax.mat'));
load(fullfile(dataPath, 'xmin.mat'));
load(fullfile(dataPath, 'Flat.mat'));
load(fullfile(dataPath, 'EcoRL.mat'));
PAR(:,4)=[];

%% functions and conditions

F=@funPoly; % polynomial functions
JF=@JfunPoly; % Jacobian of functions

Order=1:-.2:.8; % order of derivatives
t0=0; % initial time

p=EcoRL(1,:)-0.0005; %  determine perturbations
dt=xmin-xmax; % distance between the bottom of valley (stable state) and hill top (unstable state)
%% evaluate resilience

for i=1:length(Order)
    for j=1:length(PAR)
h=.01; % step size
T=100; % final time (which could be updated if it is needed)
pt0=20;pt1=30; % pt0: time when  perturbation starts, pt1: perturbation ends
IndxP=pt1/h; % index when perturbation ends

Eps=dt(j)*1e-4; % scale of recovery threshold 
[t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(j),pt0,pt1]); % fist computing for the dynamics

while 1
Xmin=x(20/h-1); % find more accurate initial value close to stable state before perturbation
ind=find(abs(x(IndxP+1:end)-Xmin)<Eps); % find the index when the state is close to the vicinity of the original stable state based on Eps (recovery threshold)

if isempty(ind) % if it needs more time to recover, run it for longer time
    T=T+200;
    if i==1
    h=h/2;
    p(j)=p(j)/2;
    end
    [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(j),pt0,pt1]);
else
EngRL(i,j)=(2*(abs(Xmin-min(x)))/(abs(Xmin-min(x))+Eps)-1)/(t(ind(1))); % evaluating resilience
break
end
end 
disp(j);
    end
end

%%  RelEff=diff(EngRL)/sum(EngRL)
figure
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];

hold on

indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);

indP1 = repmat({'Memory increases \newline potential depth'},length(indP),1);
indP2 = repmat({'Memory decreases \newline potential depth'},length(indN),1);

% abs(Rdp(indP,2)-Rdp(indP,1))./Rdp(indP,1);
EngRdp1=diff(EngRL(:,indP))./sum(EngRL(:,indP));
% scatter(abs(Rdp(indN,2)-Rdp(indN,1))./Rdp(indN,1),(EcoRL(2,indN)-EcoRL(1,indN))./EcoRL(1,indN),[],Cgr)
EngRdp2=diff(EngRL(:,indN))./sum(EngRL(:,indN));

v=violinplot([EngRdp1';EngRdp2'],[indP1;indP2],...
'Width',0.5,'Bandwidth',0.02,'ViolinColor',Cgr,'ViolinAlpha',0.004,...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);

% v.ShowData=true; v.ViolinAlpha=0.1;
 h=v.ScatterPlot; h.MarkerFaceAlpha=1; h.MarkerFaceColor='flat';
 
 h.CData=repmat(h.CData,numel(h.YData),1);
 index=h.YData<3; h.CData(index,1)=Cpr(1); h.CData(index,2)=Cpr(2); h.CData(index,3)=Cpr(3);

ax = gca;
% ax.YAxis.Scale ="log";

ylabel('Relative impact of memory on engineering resilience')
% xlabel('Memory effects on potential depth')


%% Eng----Curv and pd
figure
p1=scatter(Flat(:,1),Rdp(:,1),60,(EngRL(1,:)),'filled');

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