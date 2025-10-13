clc; clear;

% -------- list the parameter sets you want to overlay ----------------
paramSets = [
    struct('r',80,'K',3.00000,'A',0.200000,'label','baseline'), ...
    struct('r',135.7950,'K',2.7799,'A',0.0112,'label','new set')
];

% -------- common c‑range ------------------------------------------------
c_min = 0;         % left edge for parameter B
c_max = 100;         % right edge (make sure it’s large enough)
dc    = 0.001;     % c step (resolution)

figure; hold on; box on;
col = lines(numel(paramSets));     % MATLAB’s default qualitative palette

lgdStr = {};      % legend labels

for k = 1:numel(paramSets)
    p = paramSets(k);

    % fetch branches for this parameter set
    [c1,Eq1,c2,Eq2,c3,Eq3] = bifurcation_branches( ...
                        p.r, p.K, p.A, c_min, c_max, dc);

    % ---- plot: solid = stable (choose your convention), dashed = unstable
    plot(c1,Eq1,'-','Color',col(k,:),'LineWidth',2);           % zero branch
    plot(c2,Eq2,'-','Color',col(k,:),'LineWidth',2);           % upper (stable)
    plot(c3,Eq3,'--','Color',col(k,:),'LineWidth',2);          % lower (unstable)

    % collect label for legend (only once per parameter set)
    lgdStr{end+1} = p.label;
end

xlabel('Parameter B');   ylabel('Equilibrium density');
set(gca,'FontSize',13);
legend(lgdStr,'Location','best');
title('Bifurcation diagram for different (r,K,A)');
hold off;


%% ------------------------------------------------------------------------
%  FUNCTION THAT RETURNS THE BRANCHES FOR ONE PARAMETER SET
%  ------------------------------------------------------------------------
function [c1,Eq1,c2,Eq2,c3,Eq3,c_int,eq_int] = ...
            bifurcation_branches(r,K,A,c_min,c_max,dc)
    % r,K,A : model parameters
    % c_min,c_max,dc : range and step for parameter c (B in your caption)
    %
    % Returns three branches (Eq1 0‑branch, Eq2 upper, Eq3 lower)
    % together with the transcritical / saddle‑node intersection point
    % (c_int,eq_int) so you still have that information if you need it.
    %
    syms c xx
    eqn = r*xx*(1 - xx/K) - c*xx/(xx + A);        % equilibrium equation
    sc  = solve(eqn,xx);                          % 0‑branch + two radicals
    
    sc2_func = matlabFunction(sc(2),'vars',c);
    sc3_func = matlabFunction(sc(3),'vars',c);

    % ---- find intersection of sc(2) and sc(3) ----
    c_rng     = c_min:dc:c_max;
    sc2_val   = sc2_func(c_rng);
    sc3_val   = sc3_func(c_rng);
    [~,idx]   = min(abs(sc2_val - sc3_val));
    c_int     = c_rng(idx);
    eq_int    = sc2_val(idx);

    % ---- find where upper branch first becomes positive -------------
    pos_idx   = find(sc2_val > 0,1,'first');
    c0        = c_rng(pos_idx);                   % your “c0”

    % ---- assemble the three branches --------------------------------
    c1 = c0:dc:1.5*c_int;     Eq1 = zeros(size(c1));      % 0‑branch
    c2 = c0:dc:c_int;         Eq2 = sc2_func(c2);         % upper
    c3 = c_min:dc:c_int;      Eq3 = sc3_func(c3);         % lower
end
