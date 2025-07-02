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

%% plot: determining Basin sharpness index combining both flatness and depth 

dt=xmin-xmax; % distance from bottom of valleys to hill top thresholds

% Example weights
wf = 0.5; % weight for flatness
wd = 0.5; % weight for depth

% Normalize flatness and depth
normalized_flatness = (abs(CurvR(:,1)) - min(abs(CurvR(:,1)))) / (max(abs(CurvR(:,1))) - min(abs(CurvR(:,1))));
normalized_depth = (Rdp(:,1) - min(Rdp(:,1))) / (max(Rdp(:,1)) - min(Rdp(:,1)));
normalized_distance = (dt - min(dt)) / (max(dt) - min(dt));

% Calculate the new metric
M = wf *  normalized_flatness + wd * normalized_depth;


%% plot heatmap of values for first and last basin

[minM, idxMin] = min(M);
[maxM, idxMax] = max(M);
Im_M_RL = diff(EngRL) ./ sum(EngRL); % Relative impact of memory on resilience
Im_M_RS = diff(EcoRL) ./ sum(EcoRL); % Relative impact of memory on resistance
Im_M_D=((Rdp(:,11)-Rdp(:,1))./(Rdp(:,1)+Rdp(:,11))); % Relative impact of memory on depth
Im_M_F=abs((Flat(:,1)-Flat(:,11))./(Flat(:,1)+Flat(:,11))); % Relative impact of memory on flatness

Im_M_RL_AtMinM = Im_M_RL(idxMin); 
Im_M_RL_AtMaxM = Im_M_RL(idxMax);
Im_M_RS_AtMinM = Im_M_RS(idxMin);
Im_M_RS_AtMaxM = Im_M_RS(idxMax);
Im_M_D_AtMinM = Im_M_D(idxMin);
Im_M_D_AtMaxM = Im_M_D(idxMax);
Im_M_F_AtMinM = Im_M_F(idxMin);
Im_M_F_AtMaxM = Im_M_F(idxMax);


dataMatrix = [
    Im_M_F_AtMinM, Im_M_F_AtMaxM;
    Im_M_D_AtMinM, Im_M_D_AtMaxM;
    Im_M_RL_AtMinM, Im_M_RL_AtMaxM;
    Im_M_RS_AtMinM, Im_M_RS_AtMaxM
];

rowLabels = {'Flatness', 'Depth', 'Resilience', 'Resistance'};
colLabels = {'flat and shallow', 'steep and deep'};

figure;
h = heatmap(colLabels, rowLabels, dataMatrix);

% Optional: Additional customization
h.Title = 'Relative impact of memory on different metrics';
h.XLabel = 'Basin Type';
h.YLabel = 'Metric';
% h.ColorScaling = 'scaledcolumns'; % This option scales each column independently

% Define colors at specific points
colors = [0.5 0.5 0.5; % Gray (for -1)
          1 1 1;       % White (for 0)
          1 0.5 0];    % Orange (for 1)

% Create a custom colormap
customColormap = interp1([-1, 0, 1], colors, linspace(-1, 1, 256), 'linear');


h.ColorLimits = [-1, 1];  % Fix the color scale to range from -1 to 1
h.Colormap = customColormap;

% h.Colormap = jet;  % Using jet for better visual distinction, but you can choose any
% h.Colormap = parula;  % You can use other colormaps like jet, hot, cool, etc.

% Additional customization
% h.CellLabelColor = 'none';  % Optionally hide cell labels to avoid clutter
