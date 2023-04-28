clear
clc

load('PAR.mat')
PAR(:,4)=[];
load('Rdp.mat')
load('CurvR.mat')
% load('p.mat')
load('xmin.mat')
load('xmax.mat')

%% find Xmin and Xmax

% for i=1:length(PAR)
%     
% xt=0:.0001:10; % steps for polynomial function
% U=-polyint(PAR(i,:));
% U1=polyval(U,xt);
% TF=islocalmin(U1);
% TF1=islocalmax(U1);
% 
% xmin(i)=xt(TF); % critical points: bottom of valley
% xmax(i)=xt(TF1); % critical point: threshold 
% 
% end
%% EcoRL
load('EcoRL.mat')
% I1=find(EcoRL(1,:)==9.999);
% I2=find(EcoRL(1,:)==4.999);
% I3=find(EcoRL(1,:)==2.499);
% II=cat(I1,I2,I3);

p=EcoRL;
% p(1,II)=EcoRL(1,II)/2-1;

F=@funPoly;
JF=@JfunPoly;

Order=1:-.2:.8;
t0=0;
T=30;
h=0.01;

J=1;

for i=1:length(Order)
    for j=1:length(PAR)
tic
    p2=p(2,j); p1=.001; %initial guess
    pt0=10;pt1=20;
while 1

[t1,x1] = FDE_PI2_IM(Order(i),F,JF,t0,30,xmin(j),h,[PAR(j,:),p1,pt0,pt1]);
[t2,x2] = FDE_PI2_IM(Order(i),F,JF,t0,30,xmin(j),h,[PAR(j,:),p2,pt0,pt1]);

xmin(j)=x1(9.99/h);
if (x1(end)>x1(end-1) && x2(end)>x2(end-1)) || (x1(end)>xmax(j) && x2(end)>xmax(j)) || (abs(x1(end)-xmin(j))<1e-6 && abs(x2(end)-xmin(j))<1e-6)
    p(i,j)=p2;
    break
    elseif (x1(end)>xmax(j) && x2(end)<xmax(j)) || (x1(end)>x1(end-1) && x2(end)<x2(end-1)) ||  (abs(x1(end)-xmin(j))<1e-6 && x2(end)<1e-6)
    p1=(p2+p1)/2; p2=2*(p2+p1)/3; 
    if p1>p2
        p(i,j)=p1;
    break
    end
    elseif (x1(end)<x1(end-1) && x2(end)<x2(end-1)) || (x1(end)<xmax(j) && x2(end)<xmax(j)) ||  (x1(end)<1e-6 && x2(end)<1e-6)
%     else
    p1=p1/2; p2=p2/2;
%     if p2<.01 
% %         pt0=pt0/2;pt1=pt1/2;
%        I(J)=i; J=1+J; break
%     end
end

end
toc
    end
end
