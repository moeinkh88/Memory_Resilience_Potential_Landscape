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

%% computing the dynamics


F = @funCP;
JF = @JfunCP;
h=0.01;

al = .8;
t0 = 0;
T = 100;

xmin = [0, 1.9568];
xmax = 0.8432;

RX0 = xmin(2);

% Define multiple perturbation times
pulse_duration = 10; 

P = 0.27:.005:0.9; %perturbation
pt0(1) = 10; % Starting time of the first perturbation
pt1(1) = pt0(1) + pulse_duration;

for i=1:length(P)
   
    % Run the simulation with the updated perturbation strength
    [Rt, Rx] = FDE_PI2_IM(al, F, JF, t0, T, RX0, h, [P(i), pt0, pt1]);
   
    RX(i, :) = Rx;
end

%% Calculate derivatives after perturbation time and filter results
Diff_T_XU = [];
Diff_XU = [];
Speed = [];
count = 1;

for i = 1:length(P)
    % Find the index where time is greater than pt1(i)
    idx = find(Rt > pt1, 1);
    
    % Calculate derivative after pt1(i)
    dRx = diff(RX(i, idx:end)) ;
    
    % Check if the derivative is always positive
    if all(dRx > 0)
        continue;
    end
    
    % Find the first time derivative is zero or close to zero
    % zero_idx = find(abs(dRx) < 1e-5, 1);
    zero_idx = find(dRx < 0, 1);
    if isempty(zero_idx)
        continue;
    end
    
    % Calculate T_XU and XU
    Diff_T_XU(count) = Rt(zero_idx);
    XU(count) = RX(i, idx + zero_idx);
    
    % Calculate Diff_XU
    Diff_XU(count) = abs(RX(i, idx) - XU(count));
    
    % Calculate speed
    Speed(count) = Diff_XU(count) / (Diff_T_XU(count));
    
    count = count + 1;
end

% Remove unused preallocate entries
Diff_T_XU = Diff_T_XU(1:count-1);
Diff_XU = Diff_XU(1:count-1);
Speed = Speed(1:count-1);


%% Plotting
% Define colors for dynamics plots
colors = [0 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.4660 0.6740 0.1880;
          0.9290 0.6940 0.1250];
figure;

% Highlight background as perturbation for Memory = 0.2
vb = [pt0 0; pt1 0; pt1 max(max(RX(:))); pt0 max(max(RX(:)))];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);
hold on;

I=length(P) - count + 1;

% Plot the dynamics for each perturbation value
for i = 1:length(P)
    pl = plot(Rt, RX(i, :), 'color', colors(2, :));
    % xline(pt1(i), 'color', [0.5 0.5 0.5]);
   if i>I
        plot(pt1 + Diff_T_XU(i-I), XU(i-I),'k.')
   end
end

% Add reference line and labels
yline(0.8432, ':');
text(5, 0.8432, 'X_{u}', 'FontSize', 14);
title(['Dynamics for Different Perturbation Strengths Memory = 0.2, ' ...
    'Perturbations = 0.27:0.005:0.9 \newline Black Dots= when states start converging to the zero stable state.']);
set(gca, 'FontSize', 10);
xlabel('Time');
ylabel('States');
