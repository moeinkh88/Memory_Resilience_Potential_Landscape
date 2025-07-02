clear
clc


% Define the path to the data folder
dataPath = fullfile('..', '..', 'data');

% Get a list of all .mat files in the data folder
matFiles = dir(fullfile(dataPath, '*.mat'));

% Loop through each .mat file and load it
for i = 1:length(matFiles)
    load(fullfile(dataPath, matFiles(i).name));
end


%% Resilience heatmap VS depth and flatness scatter points
figure
p1=scatter(Flat(:,1),Rdp(:,1),60,(EngRL(1,:)),'filled');

xlabel('Flatness')
ylabel('Potential depth')

set(gca,'ColorScale','log') 
p1.MarkerEdgeColor=[.7,.7,.7];

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Resilience';

set(gca,'CLim',[0.0164   2.7019]);

colormap(flipud(gray))

%% Resistance heatmap VS depth and flatness scatter points
figure
p1=scatter(Flat(:,1),Rdp(:,1),60,(EcoRL(1,:)),'filled');
ax = gca;
xlabel('Flatness')
ylabel('Potential depth')
set(gca,'ColorScale','log') 
p1.MarkerEdgeColor=[.7,.7,.7];

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Resistance';

set(gca,'CLim',[0.0520 7.1130]);
colormap(flipud(gray))
