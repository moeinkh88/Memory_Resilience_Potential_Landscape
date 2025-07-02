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
h=0.01;
%%

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
c1 = c0:0.001:1.3*c_intersection;
c2 = c0:0.001:c_intersection;
c3 = 0:0.001:c_intersection;

% Handle the solutions
% The first solution is constant, so handle it separately
Eq1 = zeros(size(c1)); % Since sc(1) is 0 for all values of c1

% Evaluate the equilibrium points by substituting the values of c
Eq2 = sc2_func(c2);
Eq3 = sc3_func(c3);
%%


F = @funC_multiP;
JF = @JfunC_multiP;

al = 1:-0.2:0.8;
t0 = 0;
T = 100;

RX0 = 1.9568;

% Define multiple perturbation times
num_pulses = 1; % Number of pulses
pulse_interval = 50; % Time interval between pulses

P = 0.275; %perturbation

pt0 = zeros(2, num_pulses);
pt1 = zeros(2, num_pulses);

% c values
c_v =B;
c_v_2 =B+P;

% Evaluate sc2_func and sc3_func over the range
XU = sc2_func(c_v);
% XU2 = sc2_func(c_v_2);

for i = 1:length(al)
%     % Modify pulse duration for specific cases: Duration of each pulse (pt1 - pt0)
%     if i==1
        pulse_duration = 10; % Set pulse duration for al = 1
    % else 
    %     pulse_duration = 10; % Set pulse duration for al = 0.8
    % end

    % Update perturbation times
    pt0(i,1) = 10; % Starting time of the first perturbation
    pt1(i,1) = pt0(1) + pulse_duration;
    for j = 2:num_pulses
        pt0(i,j) = pt0(i,j-1) + pulse_interval;
        pt1(i,j) = pt0(i,j) + pulse_duration;
    end

    % Run the simulation with the updated perturbation times
    [Rt, Rx] = FDE_PI2_IM(al(i), F, JF, t0, T, RX0, h, [P, pt0(i,:), pt1(i,:), num_pulses]);

    Df = diff(Rx(:, end-1:end)');

    RX(i, :) = Rx;
end

%% Potential Landscape Calculation

h1 = h;
% Rate of change; right side of the main equation
% before perturb
dxR1 = diff(RX(1, :)') / h1;
J1 = 3:length(RX(1, :)) - 2;
dxR1(J1 - 1, 1) = 1 / (12 * h1) * (RX(1, J1 - 2)' - 8 * RX(1, J1 - 1)' + 8 * RX(1, J1 + 1)' - RX(1, J1 + 2)');

dxR2 = diff(RX(2, :)') / h1;
J2 = 3:length(RX(2, :)) - 2;
dxR2(J2 - 1, 1) = 1 / (12 * h1) * (RX(2, J2 - 2)' - 8 * RX(2, J2 - 1)' + 8 * RX(2, J2 + 1)' - RX(2, J2 + 2)');

% Sorting the states in the increasing order
DX1 = flip(dxR1);
Xall1 = flip(RX(1, 1:end-1)');

Q1 = -cumtrapz(Xall1(:, 1), DX1); % potential values

DX2 = flip(dxR2);
Xall2 = flip(RX(2, 1:end-1)');

Q2 = -cumtrapz(Xall2(:, 1), DX2); % potential values

Q22=Q2+(Q1(end)-Q2(end));
%% ploting

figure;

% Subplot 1: Bifurcation Plot
subplot(2, 2, 1);
p = plot(c1, Eq1, 'k', c2, Eq2, 'k--', c3, Eq3, 'k');
set(p, 'LineWidth', 3);
set(gca, 'FontSize', 14);
ylabel('Equilibrium density');
xlabel('Parameter B');
hold on;
xline(B, 'color', [0.5 0.5 0.5], 'LineWidth', 2);
xline(B + P, 'color', [0.5 0.5 0.5], 'LineWidth', 2);
title('(a) Bifurcation Plot');

% Define colors for dynamics plots
colors = [0 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.4660 0.6740 0.1880;
          0.9290 0.6940 0.1250];

% Subplot 2: Dynamics for Memory = 0
subplot(2, 2, 2);
% Highlight background as perturbation for Memory = 0
for j = 1:num_pulses
    vb = [pt0(1, j) 0; pt1(1, j) 0; pt1(1, j) max(max(RX(:))); pt0(1, j) max(max(RX(:)))];
    f = [1 2 3 4];
    patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);
    hold on;
end

p = plot(Rt, RX(1, :), 'color', colors(1, :));
hold on;
yline(XU, ':');
text(5, XU, 'X_{u}', 'FontSize', 14);
% yline(XU2, ':');
% text(5, XU2, 'X_{u2}', 'FontSize', 14);
title('(b) Memory = 0');
set(p, 'LineWidth', 3);
set(gca, 'FontSize', 14);
ylabel('States');

% Subplot 3: Dynamics for Memory = 0.2
subplot(2, 2, 3);
% Highlight background as perturbation for Memory = 0.2
for j = 1:num_pulses
    vb = [pt0(2, j) 0; pt1(2, j) 0; pt1(2, j) max(max(RX(:))); pt0(2, j) max(max(RX(:)))];
    f = [1 2 3 4];
    patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);
    hold on;
end

p = plot(Rt, RX(2, :), 'color', colors(2, :));
hold on;
yline(XU, ':');
text(5, XU, 'X_{u}', 'FontSize', 14);
% yline(XU2, ':');
% text(5, XU2, 'X_{u2}', 'FontSize', 14);
title('(c) Memory = 0.2');
set(p, 'LineWidth', 3);
set(gca, 'FontSize', 14);
xlabel('Time');
ylabel('States');

% Subplot 4: Potential Landscape
subplot(2, 2, 4);
p1 = plot(Xall1, Q1);
hold on;
p2 = plot(Xall2, Q22);
set(p1, 'LineWidth', 2);
set(p2, 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylabel('Potential energy');
xlabel('States');
xline(XU, '--k', 'LineWidth', 2);
% xline(XU2, '--k', 'LineWidth', 2);
title('(d) Potential Landscape');

%% Animation for panels (a), (b), and (c)
video = VideoWriter('dynamic_plots.avi'); % Create a video file named 'dynamic_plots.avi'
video.FrameRate = 10; % Set the frame rate (adjust as needed)
open(video); % Open the video for writing

% Create a figure with a fixed size
figure('Position', [200, 200, 1000, 800]); % Set position and size for consistency

Flip_Xall1=flip(Xall1); Flip_Q1=flip(Q1);
Flip_Xall2=flip(Xall2); Flip_Q2=flip(Q22);

step=50;
for i = 2:step:length(Rt)
    % Clear the figure to prevent overlapping
    clf;
    
    % Panel (a): Bifurcation Plot
    subplot(2, 2, 1);
    plot(c1, Eq1, 'k', c2, Eq2, 'k--', c3, Eq3, 'k', 'LineWidth', 2);
    set(gca, 'FontSize', 14);
    ylabel('Equilibrium density');
    xlabel('Parameter B');
    hold on;
    xline(B, 'color', [0.5 0.5 0.5], 'LineWidth', 2); 
    xline(B + P, 'color', [0.5 0.5 0.5], 'LineWidth', 2); 
    
   % Conditional plotting with an empty face and colored edge
    if Rt(i) >= pt0(1, 1) && Rt(i) <= pt1(1, 1)
        plot(B + P, Flip_Xall1(i), 'o', 'MarkerSize', 11, ...
             'MarkerEdgeColor', colors(1, :), 'MarkerFaceColor', colors(1, :), 'LineWidth', 2); % Blue outline if condition is met
    else    
        plot(B, Flip_Xall1(i), 'o', 'MarkerSize', 11, ...
             'MarkerEdgeColor', colors(1, :), 'MarkerFaceColor', colors(1, :), 'LineWidth', 2); % Blue outline, no fill
    end
    
    % Second conditional plotting with different coordinates
    if Rt(i) >= pt0(2, 1) && Rt(i) <= pt1(2, 1)
        plot(B + P, Flip_Xall2(i), 'o', 'MarkerSize', 10, ...
             'MarkerEdgeColor', colors(2, :), 'MarkerFaceColor', colors(2, :), 'LineWidth', 2); % Red outline if condition is met
    else
        plot(B, Flip_Xall2(i), 'o', 'MarkerSize', 10, ...
             'MarkerEdgeColor', colors(2, :), 'MarkerFaceColor', colors(2, :),'LineWidth', 2); % Red outline, no fill
    end
    
    title('(a) Bifurcation Plot');

    % Panel (b): Dynamics for Memory = 0
    subplot(2, 2, 2);
    for j = 1:num_pulses
    vb = [pt0(1, j) 0; pt1(1, j) 0; pt1(1, j) max(max(RX(:))); pt0(1, j) max(max(RX(:)))];
    f = [1 2 3 4];
    patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);
    hold on;
    end
    plot(Rt(1:i), RX(1, 1:i), 'color', colors(1, :), 'LineWidth', 2);
    set(gca, 'FontSize', 14);
    ylabel('States');
    title('(b) Memory = 0');
    hold on;
    yline(XU, ':');
    text(5, XU, 'X_{u}', 'FontSize', 14);
    % yline(XU2, ':');
    % text(5, XU2, 'X_{u2}', 'FontSize', 14);

    % Panel (c): Dynamics for Memory = 0.2
    subplot(2, 2, 3);
    for j = 1:num_pulses
    vb = [pt0(2, j) 0; pt1(2, j) 0; pt1(2, j) max(max(RX(:))); pt0(2, j) max(max(RX(:)))];
    f = [1 2 3 4];
    patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);
    hold on;
    end
    plot(Rt(1:i), RX(2, 1:i), 'color', colors(2, :), 'LineWidth', 2);
    set(gca, 'FontSize', 14);
    xlabel('Time');
    ylabel('States');
    title('(c) Memory = 0.2');
    hold on;
    yline(XU, ':');
    text(5, XU, 'X_{u}', 'FontSize', 14);
    % yline(XU2, ':');
    % text(5, XU2, 'X_{u2}', 'FontSize', 14);

    % Panel (d): Potential Energy vs States
    subplot(2, 2, 4);
    plot(Flip_Xall1(1:i), Flip_Q1(1:i), 'color', colors(1, :), 'LineWidth', 2);
    hold on;
    plot(Flip_Xall2(1:i), Flip_Q2(1:i), 'color', colors(2, :), 'LineWidth', 2);
    set(gca, 'FontSize', 14);
    ylabel('Potential energy');
    xlabel('States');
    xline(XU, '--k', 'LineWidth', 2);
    % xline(XU2, '--k', 'LineWidth', 2);
    title('(d) Potential Landscape');
    xlim([min(min([Xall1, Xall2])), max(max([Xall1, Xall2]))]); % Set fixed x-axis limits
    ylim([min(min(([Q1, Q22]))), max(max(([Q1, Q22])))]); % Set fixed y-axis limits


    % Add a ball (circle marker) and animate its movement along Flip_Xall1 vs Flip_Q1
    plot(Flip_Xall1(i), Flip_Q1(i)+.01, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(1, :), 'Color',colors(1, :)); % Initialize the ball
    plot(Flip_Xall2(i), Flip_Q2(i)+.01, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(2, :), 'Color', colors(2, :)); % Initialize the ball

    % Capture the frame for the video
    drawnow; % Make sure all graphics are updated before capturing
    frame = getframe(gcf);
    writeVideo(video, frame); % Write the current frame to the video
end

% Close the video writer object
close(video);
implay('dynamic_plots.avi');
