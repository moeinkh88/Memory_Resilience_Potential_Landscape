%% Combined script: 4 duration panels plus 1 large period sweep panel in one figure
clc
clear

global r K A B

%% Parameters
r = 0.8;
K = 3;
A = 0.2;
B = 0.6;

F  = @funCP;
JF = @JfunCP;

h  = 0.01;
t0 = 0;
T  =80;

pt0 = 10;

% Equilibrium values you already use in your first script
xu  = 0.8432;
RX0 = 1.9568;

%% Part 1: 4 panels, fixed P, different durations, 3 memory curves
alList = [1 0.9 0.85];
mem    = minus(1, alList);
memLabels = arrayfun(@(m) sprintf('Memory = %.2f', m), mem, 'UniformOutput', false);

pulseDurList = [5 7 7.65 7.7];
pt1List = plus(pt0, pulseDurList);

P1 = 0.3;

% Run one sim to get time length, then preallocate
[RtA, RxA] = FDE_PI2_IM(alList(1), F, JF, t0, T, RX0, h, [P1, pt0, pt1List(1)]);
nT = numel(RtA);

RXdur = zeros(numel(pt1List), numel(alList), nT);
RXdur(1,1,:) = RxA;

for i = 1:numel(pt1List)
    for j = 1:numel(alList)
        if i == 1 && j == 1
            continue
        end
        [RtTmp, RxTmp] = FDE_PI2_IM(alList(j), F, JF, t0, T, RX0, h, [P1, pt0, pt1List(i)]);
        RXdur(i,j,:) = RxTmp;
    end
end

%% Part 2: 1 panel, fixed memory, fixed P, duration sweep
al2 = 0.8;     % memory = 0.2
P2  = 0.4;

pulseDurSweep = 5:0.1:10;
pt1Sweep = plus(pt0, pulseDurSweep);

[RtB, RxB] = FDE_PI2_IM(al2, F, JF, t0, T, RX0, h, [P2, pt0, pt1Sweep(1)]);
RXper = zeros(numel(pt1Sweep), numel(RxB));
RXper(1,:) = RxB;

for i = 2:numel(pt1Sweep)
    [RtTmp, RxTmp] = FDE_PI2_IM(al2, F, JF, t0, T, RX0, h, [P2, pt0, pt1Sweep(i)]);
    RXper(i,:) = RxTmp;
end

%% Shared plot styling
memColors = [0 0.4470 0.7410;
             0.8500 0.3250 0.0980;
             0.5500 0.0000 0.0000];

shadeColor = [0.9290 0.6940 0.1250];

yMinAll = min([RXdur(:); RXper(:)]);
yMaxAll = max([RXdur(:); RXper(:)]);
pad = 0.02 * (abs(minus(yMaxAll, yMinAll)) + eps);
yL  = [minus(yMinAll, pad)  plus(yMaxAll, pad)];

%% One figure with a 2 by 3 layout
% Left side: 2 by 2 duration panels
% Right side: one large panel spanning 2 rows
figure('Color','w','Units','pixels','Position',[80 80 1450 700]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax = gobjects(5,1);

% Tile mapping for the 4 small panels
tileIds = [1 2 4 5];
panelTags = {'(a)','(b)','(c)','(d)','(e)'};

for k = 1:4
    i = k;
    ax(k) = nexttile(tileIds(k));
    hold(ax(k),'on'); box(ax(k),'on');
    set(ax(k),'FontSize',11,'LineWidth',1);

    % Shade perturbation window
    vb = [pt0        yL(1);
          pt1List(i) yL(1);
          pt1List(i) yL(2);
          pt0        yL(2)];
    patch('Faces',[1 2 3 4], 'Vertices',vb, ...
          'FaceColor', shadeColor, 'EdgeColor','none', ...
          'FaceAlpha',0.14, 'HandleVisibility','off', 'Parent',ax(k));

    % Memory curves
    hLine = gobjects(numel(alList),1);
    for j = 1:numel(alList)
        hLine(j) = plot(ax(k), RtA, squeeze(RXdur(i,j,:)), ...
                        'Color', memColors(j,:), ...
                        'LineWidth', 1.8);
        if j == numel(alList)
            set(hLine(j),'LineStyle',':')
        end
    end

    yl = yline(ax(k), xu, ':', 'LineWidth', 1.9);
    yl.HandleVisibility = 'off';
    tx = text(ax(k), 5, xu, 'X_{u}', 'FontSize',12, 'VerticalAlignment','bottom');
    tx.HandleVisibility = 'off';

    title(ax(k), sprintf('Perturbation duration = %.2f', pulseDurList(i)), 'FontWeight','normal');
    xlabel(ax(k), 'Time');
    ylabel(ax(k), 'State');
    ylim(ax(k), yL);

    text(ax(k), -0.13, 1, panelTags{k}, 'Units','normalized', ...
         'FontSize',12, 'FontWeight','normal', 'VerticalAlignment','top');

    if k == 1
        legend(ax(k), hLine, memLabels, 'Location','best');
    end
end

%% Big right panel
ax(5) = nexttile(3,[2 1]);
hold(ax(5),'on'); box(ax(5),'on');
set(ax(5),'FontSize',10,'LineWidth',1);

epsSlope = 10^(minus(0,8));

for i = 1:numel(pt1Sweep)
    t = RtB;
    x = RXper(i,:);

    plot(ax(5), t, x, 'LineWidth', 0.8, 'Color', [0.8500 0.3250 0.0980]);

    idxEnd = find(t >= pt1Sweep(i), 1, 'first');
    if ~isempty(idxEnd)
        plot(ax(5), [pt1Sweep(i) pt1Sweep(i)], [yL(1) x(idxEnd)], ...
             'Color', [0.5 0.5 0.5], 'LineWidth', 0.1);
    end

    if isempty(idxEnd) || idxEnd >= numel(x) - 2
        continue
    end

    xSeg = x(idxEnd:end);
    xSeg = movmean(xSeg, 5);

    dx = diff(xSeg);

    up   = (dx > 0) & (abs(dx) > epsSlope);
    down = (dx < 0) & (abs(dx) > epsSlope);

    m = find(up(1:end-1) & down(2:end), 1, 'first');
    if ~isempty(m)
        dotIdx = idxEnd + m;
        plot(ax(5), t(dotIdx), x(dotIdx), 'k.', 'MarkerSize', 10);
    end
end

yline(ax(5), xu, ':', 'LineWidth', 1.3);
text(ax(5), 5, xu, 'X_{u}', 'FontSize', 12, 'VerticalAlignment','bottom');

title(ax(5), sprintf(['Duration sweep 5:0.1:10,  P = %.1f,  Memory = %.1f\n' ...
                      'Black dots show first peak after perturbation end'], ...
                      P2, 1 - al2), ...
     'FontWeight','normal');

xlabel(ax(5), 'Time');
ylabel(ax(5), 'State');
ylim(ax(5), yL);

text(ax(5), -0.13, 1, panelTags{5}, 'Units','normalized', ...
     'FontSize',12, 'FontWeight','normal', 'VerticalAlignment','top');
