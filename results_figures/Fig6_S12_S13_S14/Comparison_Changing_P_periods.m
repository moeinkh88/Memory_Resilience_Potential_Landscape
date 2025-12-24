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
pulse_duration = [5, 7, 7.5, 7.9];

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

% plotting

% Prepare colors for plotting
colors = [0 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.4660 0.6740 0.1880;
          0.9290 0.6940 0.1250];

figure;

for i = 1:length(pt1)
    for j = 1:length(al)
        % Create subplot for each combination of P and al
        subplot(length(pt1), length(al), (i - 1) * length(al) + j);
        hold on;

        % Highlight background as perturbation for each Memory level
        vb = [pt0 0; pt1(i) 0; pt1(i) max(max(max((RX(:))))); pt0 max(max(max((RX(:)))))];
        f = [1 2 3 4];
        patch('Faces', f, 'Vertices', vb, 'FaceColor', colors(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.16);

        % Plot the dynamics
        if j==1
                    pl = plot(Rt, squeeze(RX(i, j, :)), 'color', colors(1, :),'LineWidth', 1.5);
        else
                    pl = plot(Rt, squeeze(RX(i, j, :)), 'color', colors(2, :),'LineWidth', 1.5);
        end

        % Add reference line and labels
        yline(xmax, ':');
        text(5, xmax, 'X_{u}', 'FontSize', 10);
        
        % Set titles and labels
        title(['Perturbation period = ' num2str(pt1(i)-pt0) ', Memory = ' num2str(1 - al(j))]);
        set(gca, 'FontSize', 8);
        xlabel('Time');
        ylabel('States');

        hold off;
    end
end

% Adjust the figure layout for better readability
% sgtitle('Dynamics for Different Perturbation Length and Memory Levels');

