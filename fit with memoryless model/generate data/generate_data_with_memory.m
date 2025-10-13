clc
clear 

global r K A B

% coefficients
r=.8;
K=6;
A=.2;
% B=1;

t0=0; % initial time
T=200; % final time
h=0.01; % computation step size

al=.8; % order of derivative (memory = 0 and 0.2)

F=@fun; % herbivory model equation function
JF=@Jfun; % Jacobian of function

%% bifurcation diagram; to see in which initial conditions run the solver
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
c1 = c0:0.001:1.5*c_intersection;
c2 = c0:0.001:c_intersection;
c3 = 0:0.001:c_intersection;

% Handle the solutions
% The first solution is constant, so handle it separately
Eq1 = zeros(size(c1)); % Since sc(1) is 0 for all values of c1

% Evaluate the equilibrium points by substituting the values of c
Eq2 = sc2_func(c2);
Eq3 = sc3_func(c3);

%let's plot the bifurcation diagram
figure
p=plot(c1,Eq1,'k',c2,Eq2,'k--',c3,Eq3,'k'); % Bifurcation diagram

set(p,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Equilibrium density')
xlabel('Paramter B')


%% Providing the dynamics with solution of the equation for different initial values

bb=.2:.1:1.2; % define some values within multistability region to generate data with initial conditions around the unstable points

X0_i=sc2_func(bb);


nb  = numel(bb);            % how many b–values you sweep
sol = struct( ...           % one “record” per parameter value
    'b',  cell(nb,1), ...   %   parameter B
    'IC', cell(nb,1), ...   %   1×2 vector of initial states
    't',  cell(nb,1), ...   %   time vector
    'x',  cell(nb,1));      %   Nt×2 matrix of trajectories

%
for i = 1:nb
    b      = bb(i);
    X01    = X0_i(i) - 1e-3;
    [t,x1] = FDE_PI2_IM(al,F,JF,t0,T,X01,h,b);

    X02    = X0_i(i) + 1e-3;
    [~,x2] = FDE_PI2_IM(al,F,JF,t0,T,X02,h,b);

    % store everything
    NtFull = numel(t);                  % full length (≈ T/h + 1)
    idx    = round(linspace(1, NtFull, T+1));   % T+1 evenly‑spaced indices
    
    sol(i).b  = b;
    sol(i).IC = [X01 X02];
    sol(i).t  = t(idx).';                      % column vector
    sol(i).x  =  [x1(idx)'  x2(idx)'];             % Nt × 2 matrix
end


%% plot the dynamics
figure;  hold on
for i = 1:nb
    plot(sol(i).t, sol(i).x(:,1), '--',  'DisplayName', ...
        sprintf('b = %.2f, IC = %.3f', sol(i).b, sol(i).IC(1)));
    plot(sol(i).t, sol(i).x(:,2), '-',   'DisplayName', ...
        sprintf('b = %.2f, IC = %.3f', sol(i).b, sol(i).IC(2)));
end
xlabel('Time');  ylabel('Trajectories');  legend show;  box on

%%
save('herbivore_dynamics_memory6.mat','sol')
%% 
tblAll = table();          % start empty

for i = 1:numel(sol)       % each parameter value (run)
    Nt = numel(sol(i).t);  % number of time steps
    
    % ---- first trajectory (perturbed –Δ) ----
    T1 = table( ...
        repmat(i,Nt,1), ...           % run index
        repmat(sol(i).b,Nt,1), ...    % parameter b
        repmat(1,Nt,1), ...           % ic_id = 1
        repmat(sol(i).IC(1),Nt,1), ...% initial condition value
        sol(i).t(:), ...              % time column
        sol(i).x(:,1));               % trajectory
    % ---- second trajectory (perturbed +Δ) ----
    T2 = table( ...
        repmat(i,Nt,1), ...
        repmat(sol(i).b,Nt,1), ...
        repmat(2,Nt,1), ...
        repmat(sol(i).IC(2),Nt,1), ...
        sol(i).t(:), ...
        sol(i).x(:,2));

    tblAll = [tblAll ; T1 ; T2];      %#ok<AGROW>
end

tblAll.Properties.VariableNames = ...
    {'run','b','ic_id','ic_val','t','x'};

% write to CSV
writetable(tblAll,'herbivore_dynamics_memory6.csv')

fprintf('CSV written: %s (%d rows)\n', ...
        'herbivore_dynamics_memory6.csv',height(tblAll));
