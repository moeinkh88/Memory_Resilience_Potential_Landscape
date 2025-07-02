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

%% position of right and left depth and curvatures
%% Fig S18 (a) left and right scatter depth

%colors
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];
figure

hold on

p1=scatter(Ldp(:,1),Rdp(:,1),color='#0072BD');
p11=scatter(Ldp(:,end),Rdp(:,end),'+');

set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Memory=0','Memory=0.2')

p11.MarkerFaceColor = '#D95319';

set(p1,'LineWidth',1)
set(p11,'LineWidth',1)
ylabel('Potential depth of right valley')
xlabel('Potential depth of left valley')


%% Fig S18 (b) left and right curvature

figure

hold on
p1=scatter((CurvL(:,1)),(CurvR(:,1)),color='#0072BD');
p11=scatter((CurvL(:,end)),(CurvR(:,end)),'+');

set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Memory=0','Memory=0.2', location='southeast')

p11.MarkerFaceColor = '#D95319';

ylabel('Curvatures on the right valley')
xlabel('Curvatures on the left valley')

axis tight


%% Fig S18 (c) Parameters statistics

figure
% subplot(2,2,3)

g1 = repmat({'a3'},length(PAR(:,1)),1);
g2 = repmat({'a2'},length(PAR(:,2)),1);
g3 = repmat({'a1'},length(PAR(:,3)),1);

v=boxplot([PAR(:,1);PAR(:,2);PAR(:,3)],[g1;g2;g3],'Colors', ...
    [.5,.5,.5; .5,.5,.5;.5,.5,.5],'symbol','','orientation','horizontal');

h = findobj(gca,'Tag','Box');

for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.7);
end

ylabel('Parameters')
xlabel('Values')

%% Fig S18 (d) Critical points statistics (extreme points, stable and unstable points)
figure

g1 = repmat({'Peak'},length(xmax),1);
g2 = repmat({'Right valley'},length(xmin),1);

v=boxplot([xmax';xmin'],[g1;g2],'Colors', ...
    [.5,.5,.5; .5,.5,.5],'symbol','','orientation','horizontal');

 h = findobj(gca,'Tag','Box');

for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.7);
end

ylabel('Critical points')
xlabel('States')

