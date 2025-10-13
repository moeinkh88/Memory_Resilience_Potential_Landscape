clc
clear 


% coefficients
r=.8;
K=6;
A=.2; 
b=.4;

F=@fun; % herbivory model equation function
JF=@Jfun; % Jacobian of function

alpha=0.8;
t0=0; % initial time
T=200; % final time

h=0.01; % computation step size

xmin=[0,5.460549503]; % stable states of the model
xmax=.8432; % unstable state of the model

X0=xmin(2);

% perturbation
p1=1; % P is strength of perturbation (b+p)
pt0=10; % the time when perturbation starts
pt1=20; % the time when perturbation ends
p2=0.15; % P is strength of perturbation (b+p, second time)
pt2=90; % the time when perturbation starts (second time)
pt3=130; % the time when perturbation ends (second time)
%%
par=[p1,p2,pt0,pt1,pt2,pt3,r,K,A,b];
[t,x] = FDE_PI2_IM(alpha,F,JF,t0,T,X0,h,par);

%% plotting the dynamics
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980;
0.4660 0.6740 0.1880
0.9290 0.6940 0.1250];

figure

 % Highlight background as perturbation
    vb = [pt0 xmax; pt1 xmax; pt1 max(max(x(:,:))); pt0 max(max(x(:,:)))];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor',colors(4,:),'EdgeColor','non', 'FaceAlpha',.16);
%     set(gca,'YScale','log')
    hold on
 % Highlight background as perturbation
    vb = [pt2 xmax; pt3 xmax; pt3 max(max(x(:,:))); pt2 max(max(x(:,:)))];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor',colors(4,:),'EdgeColor','non', 'FaceAlpha',.16);


p=plot(t,x(1,:),'color',colors(1,:));
hold on

% set(p,'LineWidth',4)
% set(gca,'FontSize',14)



xlabel('Time')
ylabel('States')

% axis([0 T xmax max(max(x(:,:)))])


%%
function dx=fun(t,x,P)

p1=P(1);
p2=P(2);
pt0=P(3);pt1=P(4);
pt2=P(5);pt3=P(6);
r=P(7);
K=P(8);
A=P(9);
B=P(10);

%%%pertubation
if t>pt0 && t<pt1

    b=B+p1;

elseif t>pt2 && t<pt3
    b=B+p2;
else
    b=B;
end
%%%

dx= (r*x*(1-x/K)-b*x/(x+A));
end

function dx=Jfun(t,x,P)
p1=P(1);
p2=P(2);
pt0=P(3);pt1=P(4);
pt2=P(5);pt3=P(6);
r=P(7);
K=P(8);
A=P(9);
B=P(10);

%%%pertubation
if t>pt0 && t<pt1

    b=B+p1;

elseif t>pt2 && t<pt3
    b=B+p2;
else
    b=B;
end
%%%

dx= r*(1-x/K)-r/K*x-b/(x+A)+b*x/(x+A)^2;
end