%% This code find the speed of unstable points shifting when there are 
% perturbations with the same strength but different periods
clc
clear

global r K A B

% coefficients
% r=80;
r=.8;
K=3;
A=.2;
% B=60;
B=.6;

%% evaluation of bifurcation

b = 0.4:0.01:2;

% Define the symbolic variables
syms c xx

% Define the equilibrium equation
eqn = (r * xx * (1 - xx / K) - c * xx / (xx + A));
sc = solve(eqn, xx); % sc will contain multiple solutions
disp(sc)

% Convert sc(2) and sc(3) to function handles
sc2_func = matlabFunction(sc(2));
sc3_func = matlabFunction(sc(3));

% Intersection func2 func3 
% Define the range of c values
c_range = 0:0.000001:5;

% Evaluate sc2_func and sc3_func over the range
sc2_values = sc2_func(c_range);
sc3_values = sc3_func(c_range);

% Compute the absolute difference between the two functions
difference_values = abs(sc2_values - sc3_values);

% Find the index where the difference is minimized (i.e., closest to intersection)
[min_difference, min_index] = min(difference_values);

% Get the value of c and the corresponding equilibrium density at the intersection
c_intersection = c_range(min_index);
eq_intersection = sc2_values(min_index); % or sc3_values(min_index), they are approximately equal


% Evaluate sc2_func over the range
sc2_values = sc2_func(c_range);

% Find the indices where sc2_func is positive
positive_indices = find(sc2_values > 0);

% Get the corresponding c0 values where sc2_func start being positive
c0= c_range(positive_indices(1));
%

% Define the ranges for c1, c2, c3
c1 = c0:0.001:2*c_intersection;
c2 = c0:0.001:c_intersection;
c3 = 0:0.001:c_intersection;

% Handle the solutions
% The first solution is constant, so handle it separately
Eq1 = zeros(size(c1)); % Since sc(1) is 0 for all values of c1

% Evaluate the equilibrium points by substituting the values of c
Eq2 = sc2_func(c2);
Eq3 = sc3_func(c3);

%% computing the dynamics


F = @funCP;
JF = @JfunCP;
h=0.01;

al = .8;
t0 = 0;
T = 100;

xmin = [0, sc3_func(B)];
xmax = sc2_func(B);

RX0 = xmin(2);

% Define multiple perturbation times
pulse_duration = 5:.1:10; 

P = 0.4; %perturbation
pt0 = 10; % Starting time of the first perturbation
pt1 = pt0 + pulse_duration;

% Run the simulation for each perturbation period
for i = 1:length(pt1)
    % Run the simulation with the updated perturbation strength
    [Rt, Rx] = FDE_PI2_IM(al, F, JF, t0, T, RX0, h, [P, pt0, pt1(i)]);
    RX(i, :) = Rx;
end


%% Plotting
colors = [0 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.4660 0.6740 0.1880;
          0.9290 0.6940 0.1250];

figure; hold on; box on;

% small threshold to ignore tiny numerical wiggles
epsSlope = 1e-8;

for i = 1:length(pt1)

    x = RX(i,:);
    t = Rt;

    % plot trajectory
    plot(t, x, 'Color', colors(2,:), 'LineWidth', 1.0);

    % vertical line at perturbation end
    idxEnd = find(t >= pt1(i), 1, 'first');
    if ~isempty(idxEnd)
        plot([pt1(i) pt1(i)], [0 x(idxEnd)], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1);
    end

    % find the dot: first peak after pt1(i)
    if isempty(idxEnd) || idxEnd >= numel(x)-2
        continue
    end

    % optional smoothing to avoid false sign flips
    xSeg = x(idxEnd:end);
    xSeg = movmean(xSeg, 5);

    dx = diff(xSeg);                 % slope sign (dt is constant in your sim)
    up   = dx >  epsSlope;
    down = dx < -epsSlope;

    % first index where it was going up then starts going down
    m = find(up(1:end-1) & down(2:end), 1, 'first');

    if ~isempty(m)
        dotIdx = idxEnd + m;         % local max is at the turning sample
        plot(t(dotIdx), x(dotIdx), 'k.', 'MarkerSize', 10);
    end
end

yline(0.8432, ':');
text(5, 0.8432, 'X_{u}', 'FontSize', 14);
title(['Dynamics for Different Perturbation Periods = 5:.1:10, Perturbation Strength=0.4, and Memory = 0.2 ' ...
    '\newline Black Dots= when states start converging to the zero stable state.']);

set(gca, 'FontSize', 10);
xlabel('Time');
ylabel('States');
