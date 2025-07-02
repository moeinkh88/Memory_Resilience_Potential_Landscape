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


%% Spearman Correlations

% Initialize correlation matrix (3 rows, 2 columns)
r = nan(3,2);

% Memory on Resistance
r(1,1) = corr(Flat(:,1), (diff(EcoRL)./sum(EcoRL))', 'Type', 'Spearman');    % Flatness
r(2,1) = corr(Rdp(:,1), (diff(EcoRL)./sum(EcoRL))', 'Type', 'Spearman');     % Depth

% Memory on Resilience
r(1,2) = corr(Flat(:,1), abs(diff(EngRL)./sum(EngRL))', 'Type', 'Spearman'); % Flatness
r(2,2) = corr(Rdp(:,1), abs(diff(EngRL)./sum(EngRL))', 'Type', 'Spearman');  % Depth

% Basin Index
norm_flat = (abs(CurvR(:,1)) - min(abs(CurvR(:,1)))) / (max(abs(CurvR(:,1))) - min(abs(CurvR(:,1))));
norm_depth = (Rdp(:,1) - min(Rdp(:,1))) / (max(Rdp(:,1)) - min(Rdp(:,1)));
M = 0.5 * norm_flat + 0.5 * norm_depth;

r(3,1) = corr(M, (diff(EcoRL)./sum(EcoRL))', 'Type', 'Spearman');    % Basin Index vs Resistance
r(3,2) = corr(M, abs(diff(EngRL)./sum(EngRL))', 'Type', 'Spearman'); % Basin Index vs Resilience


%% Plot as heatmap
figure;
imagesc(r);
axis equal tight;

% Define colormap (blue-white-red)
nColors = 256;
softness = 0.3;
blue = [linspace(softness,1,nColors/2)', linspace(softness,1,nColors/2)', ones(nColors/2,1)];
red = [ones(nColors/2,1), linspace(1,softness,nColors/2)', linspace(1,softness,nColors/2)'];
colormap([blue; red]);
colorbar;

% Labels
xticks(1:2);
yticks(1:3);
xticklabels({'Relative impact of \newline memory on resistance', 'Relative impact of \newline memory on resilience'});
yticklabels({'Flatness', 'Depth', 'Basin Index'});
xtickangle(0);  % Force horizontal tick labels
% Set font and appearance
set(gca, 'FontSize', 12, 'TickLength', [0 0]);

% Show correlation values
for i = 1:3
    for j = 1:2
        text(j, i, sprintf('%.2f', r(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', 'k', 'FontSize', 11);
    end
end

% Title
title('(b) Spearman correlation heatmap', 'FontSize', 13);
