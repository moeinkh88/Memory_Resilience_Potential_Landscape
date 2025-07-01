clear
clc

load('Rdp.mat')
load('CurvR.mat')
load('EcoRL.mat')
load('EngRL.mat')
load('Flat.mat')

%% Relative impact of memory on resilience VS Flatness

figure
p1=scatter(Flat(:,1),(diff(EngRL)./sum(EngRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
ax = gca;
xlabel('Flatness')
ylabel('Relative impact of memory on resilience')
p1.MarkerFaceColor = 'k';
p1.MarkerEdgeColor='w';

ax.YAxis.Scale ="log";

corEffMEngCurv=corrcoef(Flat(:,1),diff(EngRL)./sum(EngRL));
corEffMEngcurv=corEffMEngCurv(2);
display(corEffMEngcurv)

%% Relative impact of memory on resilience VS Depth
figure
p11=scatter(Rdp(:,1),(diff(EngRL)./sum(EngRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);


xlabel('Potential depth')
ylabel('Relative impact of memory on resilience')
p11.MarkerFaceColor = 'k';
p11.MarkerEdgeColor='w';

ax = gca;
ax.XAxis.Scale ="log";
ax.YAxis.Scale ="log";

corEffMEngDp=corrcoef(Rdp(:,1),diff(EngRL)./sum(EngRL));
corEffMEngdp=corEffMEngDp(2);
display(corEffMEngdp)

%% Relative impact of memory on resistance VS Depth

figure
p11=scatter(Rdp(:,1),(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);

ylabel('Relative effect of memory on resistance')
xlabel('Potential depth')
p11.MarkerFaceColor = 'k';
p11.MarkerEdgeColor='w';

ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";

colormap cool

[corEffMEcoDp,p1]=corrcoef(Rdp(:,1),(diff(EcoRL)./sum(EcoRL)));
corEffMEcodp=corEffMEcoDp(2);
display(corEffMEcodp)

%% Relative impact of memory on resistance VS Flatness

figure
p1=scatter(Flat(:,1),(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
ylabel('Relative effect of memory on EcoRL')
xlabel('Flatness')
set(gca,'yscale','log')
p1.MarkerFaceColor = 'k';
p1.MarkerEdgeColor='w';
% 

corEffMEcoCurv=corrcoef(Flat(:,1),(diff(EcoRL)./sum(EcoRL)));
corEffMEcocurv=corEffMEcoCurv(2);
display(corEffMEcocurv)
%% Relative impact of memory on resilience VS Basin sharpness index
% plot: Xaxis determined by new metric combining both flatness and depth 

% Example weights
wf = 0.5; % weight for flatness
wd = 0.5; % weight for depth

% Normalize flatness and depth
normalized_flatness = (abs(CurvR(:,1)) - min(abs(CurvR(:,1)))) / (max(abs(CurvR(:,1))) - min(abs(CurvR(:,1))));
normalized_depth = (Rdp(:,1) - min(Rdp(:,1))) / (max(Rdp(:,1)) - min(Rdp(:,1)));
% normalized_distance = (dt - min(dt)) / (max(dt) - min(dt));

% Calculate the new metric
% M = 1 - (wf *  normalized_flatness + wd * normalized_depth);
M = wf *  normalized_flatness + wd * normalized_depth;
% M = 2/4 *  normalized_flatness + 1/4 * normalized_depth +1/4 * normalized_distance;
% M = 1/2 *  normalized_flatness + 1/2 * sqrt( normalized_depth.^2 + normalized_distance.^2 );


figure
p11=scatter(M,(diff(EngRL)./sum(EngRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);


xlabel('Basin sharpness index')
ylabel('Relative impact of memory on resilience')
p11.MarkerFaceColor = 'k';
p11.MarkerEdgeColor='w';


%% Relative impact of memory on resistance VS Basin sharpness index
figure
p12=scatter(M,abs(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);


xlabel('Basin sharpness index')
ylabel('Relative impact of memory on resistance')
p12.MarkerFaceColor = 'k';
p12.MarkerEdgeColor='w';

ax = gca;
ax.YAxis.Scale ="log";
