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

%% Fig a: dp vs curve

figure

hold on
t = tiledlayout(1,1);

ax1 = axes(t);
plot(ax1,Rdp(:,1),(CurvR(:,1)),'o',color='#0072BD');
ylabel('Potential curvature when memory=0')
xlabel('Potential depth when when memory=0')
set(gca,'xscale','log')
set(gca,'yscale','log')
ax1.XColor = '#0072BD';
ax1.YColor = '#0072BD';
% legend('Memory=0')


ax2 = axes(t);
plot(ax2,Rdp(:,11),(CurvR(:,11)),'+',color='#D95319');

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XColor = '#D95319';
ax2.YColor = '#D95319';
ax1.Box = 'off';
ax2.Box = 'off';
ax2.Color = 'none';


set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel('Potential curvature when memory=0.2')
xlabel('Potential depth when when memory=0.2')

axis tight

%% Fig b: box plot
figure
g3 = repmat({'Memory=0'},length(CurvR(:,1)),1);
g4 = repmat({'Memory=0.2'},length(CurvR(:,end)),1);

boxplot([(CurvR(:,1));(CurvR(:,end))],[g3;g4],'Colors',[0 0.4470 0.7410; 0.8500 0.3250 0.0980],'symbol','');

h = findobj(gca,'Tag','Box');
% Ch = findobj(gcf,'tag','Outliers');
% set(Ch,'Color',[1 1 1; 1 1 1]);

for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.7);
end

ax = gca;
ax.YAxis.Scale ="log";

ylabel('Potential curvature')

% plot(1:10);
% ha = annotation('arrow');  % store the arrow information in ha
% ha.Parent = gca;           % associate the arrow the the current axes
% ha.X = [2 1];          % the location in data units
% ha.Y = [10 -1];   
% 
% ha.LineWidth  = 3;          % make the arrow bolder for the picture
% ha.HeadWidth  = 30;
% ha.HeadLength = 30;
