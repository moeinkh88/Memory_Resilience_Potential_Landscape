%% Plotting: 4 panels (one per perturbation duration), 3 memory curves per panel

% plot comparison of model dynamics with different memory and perturbation
% duration
clc;
clear;

global r K A B

% coefficients
r = 0.8;
K = 3;
A = 0.2;
B = 0.6;

% Computing the dynamics
F = @funCP;
JF = @JfunCP;
h = 0.01;

al = [1 0.9 0.85]; % Memory parameter
t0 = 0;
T = 100;

xmin = [0, 1.9568];
xmax = 0.8432;

RX0 = xmin(2);

% Define multiple perturbation times
pulse_duration = [5, 7, 7.65, 7.7];

P = 0.3; % Perturbation strengths
pt0 = 10; % Starting time of the first perturbation
pt1 = pt0 + pulse_duration;


% Run simulations for each combination of P and al
for i = 1:length(pt1)
    for j = 1:length(al)
        % Run the simulation with the updated perturbation duration
        [Rt, Rx] = FDE_PI2_IM(al(j), F, JF, t0, T, RX0, h, [P, pt0, pt1(i)]);

        % Store results
        RX(i, j, :) = Rx;
    end
end

%% plotting
mem = 1 - al;
memLabels = arrayfun(@(m) sprintf('Memory = %.2f', m), mem, 'UniformOutput', false);

% Colors for the 3 memory levels (same order as al)
memColors = [0 0.4470 0.7410;      % memory 0.00
             0.8500 0.3250 0.0980; % memory 0.10
             0.5500 0.0000 0.0000];% memory 0.15

% Shade color for perturbation window
shadeColor = [0.9290 0.6940 0.1250];

% Global y limits for consistent scaling across panels
yMin = min(RX(:));
yMax = max(RX(:));
pad  = 0.06 * (yMax - yMin + eps);
yL   = [yMin - pad, yMax + pad];

figure('Color','w','Units','pixels','Position',[80 80 1200 750]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for i = 1:length(pt1)
    nexttile; hold on; box on;
    set(gca,'FontSize',11,'LineWidth',1);

    % Shade perturbation interval (do not show in legend)
    vb = [pt0    yL(1);
          pt1(i) yL(1);
          pt1(i) yL(2);
          pt0    yL(2)];
    patch('Faces',[1 2 3 4], 'Vertices',vb, ...
          'FaceColor', shadeColor, 'EdgeColor','none', ...
          'FaceAlpha',0.14, 'HandleVisibility','off');

    % Plot the 3 memory curves and store handles
    hLine = gobjects(length(al),1);
    for j = 1:length(al)
        ls = '-';
        if j == length(al)   % memory = 0.15 (last one)
            ls = '--';
        end

        hLine(j) = plot(Rt, squeeze(RX(i,j,:)), ...
                        'Color', memColors(j,:), ...
                        'LineStyle', ls, ...
                        'LineWidth', 1.8);
    end

    % Reference line (do not show in legend)
    yl = yline(xmax, ':', 'LineWidth', 1.9);
    yl.HandleVisibility = 'off';
    tx = text(pt0 - 5, xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');
    tx.HandleVisibility = 'off';

    title(sprintf('Perturbation duration = %.2f', pt1(i) - pt0));
    xlabel('Time');
    ylabel('State');
    ylim(yL);

    % Only one legend, and only for the lines (first panel only)
    if i == 1
        legend(hLine, memLabels, 'Location','best');
    end

    hold off;
end
