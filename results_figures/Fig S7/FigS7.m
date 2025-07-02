clc
clear

% Define the path to the data folder
dataPath = fullfile('..', '..', 'data');

% Get a list of all .mat files in the data folder
matFiles = dir(fullfile(dataPath, '*.mat'));

% Loop through each .mat file and load it
for i = 1:length(matFiles)
    load(fullfile(dataPath, matFiles(i).name));
end

xt=0:.00001:5; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 


Nsample=1000;
for i=1:Nsample


%FLatness
% Flat(i,:)=ceil(max(max(abs(CurvR))))+CurvR(i,:);

[R,pv2]=corrcoef(Rdp(i,:),1-Order); % correlation between Depth and memory
corrRdp(i)=R(2);
[R,pv7]=corrcoef(Flat(i,:),1-Order); % correlation between Flatness and memory
corrCurvR(i)=R(2);

[RHO,pv9] = corr(Flat(i,:)',1-Order','Type','Spearman'); % correlation between Depth and memory (Spearman)
corrCurvRSpear(i)=RHO;


pv(i)=mean([pv2(2),pv7(2),pv9]); %average of p values

end

%colors
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];

indPc=find(corrRdp>0); % memory increases the depth
indNc=find(corrRdp<0); % memory decreases the depth


%% Fig S4a_b: Memory effects on potential depth in terms of potential depth and curvature and initial slope of potential


[CRcurv,iaRcurv] = unique(Flat(:,1));
figure

indPc=find(corrRdp(iaRcurv)>0);
indNc=find(corrRdp(iaRcurv)<0);

plot(CRcurv(indNc),corrRdp(iaRcurv(indNc)),'o','color',Cpr)
hold on
plot(CRcurv(indPc),corrRdp(iaRcurv(indPc)),'o','color',Cgr)
xlabel('Flatness')
ylabel('Correlation: memory and potential depth')


[CRdp,iaRdp] = unique(Rdp(:,1));
figure

indPc1=find(corrRdp(iaRdp)>0);
indNc1=find(corrRdp(iaRdp)<0);

semilogx(CRdp(indNc1),corrRdp(iaRdp(indNc1)),'o','color',Cpr)
hold on
semilogx(CRdp(indPc1),corrRdp(iaRdp(indPc1)),'o','color',Cgr)
xlabel('Potential depth')
ylabel('Correlation: memory and potential depth')

%%delta curv vs dp

figure

hold on

indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);
scatter(Rdp(indP,1),abs((Flat(indP,1)-Flat(indP,11))./(Flat(indP,1)+Flat(indP,11))),[],Cgr)
scatter(Rdp(indN,1),abs((Flat(indN,1)-Flat(indN,11))./(Flat(indN,1)+Flat(indN,11))),[],Cpr)

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Potential depth')
ylabel('Relative impact of memory on flatness')

%%delta curv vs curv

figure

hold on

scatter(Flat(indP,1),abs((Flat(indP,1)-Flat(indP,11))./(Flat(indP,1)+Flat(indP,11))),[],Cgr)
scatter(Flat(indN,1),abs((Flat(indN,1)-Flat(indN,11))./(Flat(indN,1)+Flat(indN,11))),[],Cpr)

% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Flatness')
ylabel('Relative impact of memory on flatness')

% set(gca,'FontSize',14)





%% plot: Xaxis determined by new metric combining both flatness and depth 

% Example weights
wf = 0.5; % weight for flatness
wd = 0.5; % weight for depth

% Normalize flatness and depth
normalized_flatness = (abs(CurvR(:,1)) - min(abs(CurvR(:,1)))) / (max(abs(CurvR(:,1))) - min(abs(CurvR(:,1))));
normalized_depth = (Rdp(:,1) - min(Rdp(:,1))) / (max(Rdp(:,1)) - min(Rdp(:,1)));

% Calculate the new metric
M = wf *  normalized_flatness + wd * normalized_depth;

%%imapct of memory on Flatness
figure

hold on
indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);

scatter(M(indP,1),abs((Flat(indP,1)-Flat(indP,11))./(Flat(indP,1)+Flat(indP,11))),[],Cgr)
scatter(M(indN,1),abs((Flat(indN,1)-Flat(indN,11))./(Flat(indN,1)+Flat(indN,11))),[],Cpr)

% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Basin sharpness index')
ylabel('Relative impact of memory on flatness')

%%imapct of memory on Depth

figure

hold on

scatter(M(indP,1),((Rdp(indP,11)-Rdp(indP,1))./(Rdp(indP,1)+Rdp(indP,11))),[],Cgr)
scatter(M(indN,1),((Rdp(indN,11)-Rdp(indN,1))./(Rdp(indN,1)+Rdp(indN,11))),[],Cpr)

% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Basin sharpness index')
ylabel('Relative impact of memory on drpth')
