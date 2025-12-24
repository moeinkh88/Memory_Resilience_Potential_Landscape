%% Better visualization: 2x2 panels (one per P), overlay memory curves
% plot comparison of model dynamics with different memory and perturbation strength
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
% pulse_duration = 10;
% 
% P = [0.18 0.223 0.225 0.23]; % Perturbation strengths
pulse_duration = 5;

P = [0.31 0.38 0.42 0.43]; % Perturbation strengths
pt0(1) = 10; % Starting time of the first perturbation
pt1(1) = pt0(1) + pulse_duration;


% Run simulations for each combination of P and al
for i = 1:length(P)
    for j = 1:length(al)
        % Run the simulation with the updated perturbation strength
        [Rt, Rx] = FDE_PI2_IM(al(j), F, JF, t0, T, RX0, h, [P(i), pt0, pt1]);

        % Store results
        RX(i, j, :) = Rx;
    end
end

%% plotting

mem = 1 - al;  % memory amount (what you show in title)

% Force plotting order and styling by memory value
memTarget       = [0, 0.10, 0.15];
memLabelsTarget = arrayfun(@(m) sprintf('Memory = %.2f', m), memTarget, 'UniformOutput', false);

% Colors: mem 0 blue, mem 0.10 red, mem 0.15 dark red dashed
memColorMap = [
    0.0000 0.4470 0.7410   % blue
    1.0000 0.0000 0.0000   % red
    0.5500 0.0000 0.0000   % dark red
];
memStyleMap = {'-','-','--'};

% Global y limits for consistent scaling across panels
yMin = min(RX(:));
yMax = max(RX(:));
pad  = 0.06 * (yMax - yMin + eps);
yL   = [yMin - pad, yMax + pad];

figure('Color','w','Units','pixels','Position',[80 80 1400 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for i = 1:numel(P)
    nexttile; hold on; box on;
    set(gca,'FontSize',12,'LineWidth',1);

    xlim([t0 T]);
    ylim(yL);

    yl = ylim;
    patch([pt0 pt1 pt1 pt0], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.75 0.2], 'EdgeColor','none', 'FaceAlpha', 0.1 + 0.2*(i-1));

    % Plot curves in the forced order: 0, 0.10, 0.15
    h = gobjects(numel(memTarget),1);
    for k = 1:numel(memTarget)
        j = find(abs(mem - memTarget(k)) < 1e-9, 1, 'first');
        if isempty(j)
            error('Could not find memory value %.2f in mem = 1 - al.', memTarget(k));
        end
        h(k) = plot(Rt, squeeze(RX(i,j,:)), ...
            'LineWidth', 2.0, ...
            'Color', memColorMap(k,:), ...
            'LineStyle', memStyleMap{k});
    end

    yline(xmax, ':', 'LineWidth', 1.9);
    text(t0 + 0.02*(T-t0), xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');

    % grid on;
    title(sprintf('Perturbation = B + %.2f', P(i)), 'FontWeight','bold');
    xlabel('Time');
    ylabel('State');

    if i == 1
        legend(h, memLabelsTarget, 'Location','best', 'Box','off');
    end
end
