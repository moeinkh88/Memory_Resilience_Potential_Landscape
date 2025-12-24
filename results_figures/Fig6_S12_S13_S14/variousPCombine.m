%% Combined figure: 4 panel comparison + 1 sweep panel in the same figure
% Assumes funCP, JfunCP, and FDE_PI2_IM are on the MATLAB path.

clc; clear;
global r K A B

% coefficients
r = 0.8;
K = 3;
A = 0.2;
B = 0.6;

F  = @funCP;
JF = @JfunCP;

dt = 0.01;
t0 = 0;
T  = 70;

xmin = [0, 1.9568];
xmax = 0.8432;
RX0  = xmin(2);

%% Part 1: 2x2 panels, fixed P list, multiple memory levels
al1 = [1 0.9 0.85];          % memory parameter
mem1 = 1 - al1;              % memory amount

pulse_duration1 = 5;
P1 = [0.31 0.38 0.42 0.43];
pt0_1 = 10;
pt1_1 = pt0_1 + pulse_duration1;

RX4 = [];
Rt1 = [];

for i = 1:numel(P1)
    for j = 1:numel(al1)
        [Rt_tmp, Rx_tmp] = FDE_PI2_IM(al1(j), F, JF, t0, T, RX0, dt, [P1(i), pt0_1, pt1_1]);

        if isempty(Rt1), Rt1 = Rt_tmp; end
        RX4(i, j, :) = Rx_tmp;
    end
end

% Forced order and styling by memory value
memTarget = [0, 0.10, 0.15];
memLabelsTarget = arrayfun(@(m) sprintf('Memory = %.2f', m), memTarget, 'UniformOutput', false);

memColorMap = [
    0.0000 0.4470 0.7410   % blue
    1.0000 0.0000 0.0000   % red
    0.5500 0.0000 0.0000   % dark red
];
memStyleMap = {'-','-','--'};

%% Part 2: one panel sweep, fixed memory, many perturbations
al2  = 0.8;                  % memory parameter
mem2 = 1 - al2;              % should be 0.2

pulse_duration2 = 10;
P2 = 0.27:0.005:0.9;
pt0_2 = 10;
pt1_2 = pt0_2 + pulse_duration2;

RXsweep = [];
Rt2 = [];

for i = 1:numel(P2)
    [Rt_tmp, Rx_tmp] = FDE_PI2_IM(al2, F, JF, t0, T, RX0, dt, [P2(i), pt0_2, pt1_2]);

    if isempty(Rt2), Rt2 = Rt_tmp; end
    RXsweep(i, :) = Rx_tmp(:).';
end

% Detect first post pulse sign change toward decreasing state
tPoint = [];
xPoint = [];
pPoint = [];
dX = [];
dT = [];
speed = [];

for i = 1:numel(P2)
    idx = find(Rt2 > pt1_2, 1, 'first');
    if isempty(idx) || idx >= numel(Rt2), continue; end

    seg = RXsweep(i, idx:end);
    dSeg = diff(seg);

    if all(dSeg > 0)
        continue
    end

    zero_idx = find(dSeg < 0, 1, 'first');
    if isempty(zero_idx)
        continue
    end

    t_i = Rt2(idx + zero_idx);
    x_i = RXsweep(i, idx + zero_idx);

    dt_i = t_i - pt1_2;
    dx_i = abs(RXsweep(i, idx) - x_i);

    if dt_i <= 0
        continue
    end

    tPoint(end+1) = t_i;
    xPoint(end+1) = x_i;
    pPoint(end+1) = P2(i);
    dT(end+1) = dt_i;
    dX(end+1) = dx_i;
    speed(end+1) = dx_i / dt_i;
end

%% Global limits for consistent scaling across all tiles
allY = [RX4(:); RXsweep(:)];
yMin = min(allY);
yMax = max(allY);
pad  = 0.06 * (yMax - yMin + eps);
yL   = [yMin - pad, yMax + pad];

%% Plot: 3x2 tiled layout, bottom tile spans both columns
figure('Color','w','Units','pixels','Position',[60 60 1000 900]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

panelLetters = {'(a)','(b)','(c)','(d)','(e)'};

% First four tiles
for i = 1:numel(P1)
    ax = nexttile; hold(ax,'on'); box(ax,'on');
    set(ax,'FontSize',12,'LineWidth',1);

    xlim(ax,[t0 T]);
    ylim(ax,yL);

    yl = ylim(ax);
    patch(ax, [pt0_1 pt1_1 pt1_1 pt0_1], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.75 0.2], 'EdgeColor','none', 'FaceAlpha', 0.12);

    hLine = gobjects(numel(memTarget),1);
    for k = 1:numel(memTarget)
        j = find(abs(mem1 - memTarget(k)) < 1e-9, 1, 'first');
        if isempty(j)
            error('Could not find memory value %.2f in mem = 1 - al1.', memTarget(k));
        end

        hLine(k) = plot(ax, Rt1, squeeze(RX4(i,j,:)), ...
            'LineWidth', 2.0, 'Color', memColorMap(k,:), 'LineStyle', memStyleMap{k});
    end

    yline(ax, xmax, ':', 'LineWidth', 1.9);
    text(ax, t0 + 0.02*(T-t0), xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');

    title(ax, sprintf('Perturbation = B + %.2f', P1(i)), 'FontWeight','normal');
    xlabel(ax,'Time');
    ylabel(ax,'State');

    text(ax, -0.13, 1, panelLetters{i}, 'Units','normalized');

    if i == 1
        legend(ax, hLine, memLabelsTarget, 'Location','best', 'Box','off');
    end
end

% Bottom tile spanning both columns
ax5 = nexttile([1 2]); hold(ax5,'on'); box(ax5,'on');
set(ax5,'FontSize',11,'LineWidth',1);

xlim(ax5,[t0 T]);
ylim(ax5,yL);

yl = ylim(ax5);
patch(ax5, [pt0_2 pt1_2 pt1_2 pt0_2], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9290 0.6940 0.1250], 'EdgeColor','none', 'FaceAlpha', 0.14);

% Many trajectories, keep lines thin to reduce clutter
for i = 1:numel(P2)
    plot(ax5, Rt2, RXsweep(i,:), 'LineWidth', 0.6, 'Color', [0.85 0.3250 0.0980 0.35]);
end

% Overlay detected points
if ~isempty(tPoint)
    plot(ax5, tPoint, xPoint, 'k.', 'MarkerSize', 4.5);
end

yline(ax5, xmax, ':', 'LineWidth', 1.9);
text(ax5, t0 + 0.02*(T-t0), xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');

title(ax5, sprintf(['Sweep of perturbations, Memory = %.2f, Perturbations = B + %.2f:%.3f:%.2f\n',...
    'Black dots: first post pulse time where trajectory starts decreasing'], mem2, P2(1), P2(2)-P2(1), P2(end)), ...
    'FontWeight','normal');

xlabel(ax5,'Time');
ylabel(ax5,'State');

text(ax5, -.06, 1, panelLetters{5}, 'Units','normalized');

% Optional: overall title for the whole figure
% title(tl, 'Model dynamics: fixed P comparison (top) and perturbation sweep (bottom)', ...
    % 'FontSize', 14, 'FontWeight','bold');


%%
%% Plot: left 2 by 2 for (a)(b)(c)(d), right column for (e)
figure('Color','w','Units','pixels','Position',[60 60 1300 650]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

panelLetters = {'(a)','(b)','(c)','(d)','(e)'};

% Put the four small panels on the left side
tileSmall = [1 2 4 5];   % left two columns in a 2 by 3 layout

for i = 1:numel(P1)
    ax = nexttile(tl, tileSmall(i));
    hold(ax,'on'); box(ax,'on');
    set(ax,'FontSize',12,'LineWidth',1);

    xlim(ax,[t0 T]);
    ylim(ax,yL);

    yl = ylim(ax);
    patch(ax, [pt0_1 pt1_1 pt1_1 pt0_1], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.75 0.2], 'EdgeColor','none', 'FaceAlpha', 0.12);

    hLine = gobjects(numel(memTarget),1);
    for k = 1:numel(memTarget)
        j = find(abs(mem1 - memTarget(k)) < 1e-9, 1, 'first');
        if isempty(j)
            error('Could not find memory value %.2f in mem = 1 - al1.', memTarget(k));
        end

        hLine(k) = plot(ax, Rt1, squeeze(RX4(i,j,:)), ...
            'LineWidth', 2.0, 'Color', memColorMap(k,:), 'LineStyle', memStyleMap{k});
    end

    yline(ax, xmax, ':', 'LineWidth', 1.9);
    text(ax, t0 + 0.02*(T-t0), xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');

    title(ax, sprintf('Perturbation = B + %.2f', P1(i)), 'FontWeight','normal');
    xlabel(ax,'Time');
    ylabel(ax,'State');

    text(ax, -0.13, 1, panelLetters{i}, 'Units','normalized');

    if i == 1
        legend(ax, hLine, memLabelsTarget, 'Location','best', 'Box','off');
    end
end

% Put panel (e) on the right side spanning both rows
ax5 = nexttile(tl, 3, [2 1]);
hold(ax5,'on'); box(ax5,'on');
set(ax5,'FontSize',11,'LineWidth',1);

xlim(ax5,[t0 T]);
ylim(ax5,yL);

yl = ylim(ax5);
patch(ax5, [pt0_2 pt1_2 pt1_2 pt0_2], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9290 0.6940 0.1250], 'EdgeColor','none', 'FaceAlpha', 0.14);

for i = 1:numel(P2)
    plot(ax5, Rt2, RXsweep(i,:), 'LineWidth', 0.6, 'Color', [0.85 0.3250 0.0980 0.35]);
end

if ~isempty(tPoint)
    plot(ax5, tPoint, xPoint, 'k.', 'MarkerSize', 4.5);
end

yline(ax5, xmax, ':', 'LineWidth', 1.9);
text(ax5, t0 + 0.02*(T-t0), xmax, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');

title(ax5, sprintf('Sweep of perturbations, Memory = %.2f, \n Perturbations = B + %.2f:%.3f:%.2f', ...
    mem2, P2(1), P2(2)-P2(1), P2(end)), ...
    'FontWeight','normal');

xlabel(ax5,'Time');
ylabel(ax5,'State');

text(ax5, -.06, 1, panelLetters{5}, 'Units','normalized');
