clc
clear
global rr KK lambda aa Sigma dW W h noise
rr=1; 
KK=10;
lambda=2.75;
aa=1.6;
X01=3;
F1=@fun2;
T=1000;
t0=0;
h=0.01;
% W=OUP(T/h+1,1);
% dW=diff(W);
noise= 2*randn(1,T/h+1)+0;
al=1;
h=0.01;
    [t1,x1] = FDE_PI12_PC(al,F1,t0,T,X01,h);
    aa=1.5;
        [t08,x08] = FDE_PI12_PC(al,F1,t0,T,X01,h);

     
%%        
colors=[0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.6350 0.0780 0.1840];
figure

subplot(2,2,1)

p1=plot(x1,t1,'color',colors(1,:));
axis([0,7,0,T])
subplot(2,2,2)
p08=plot(x08,t08,'color',colors(2,:));
axis([0,7,0,T])


% subplot(2,4,4)
% % plot(W,t1)
% plot(noise,t1)

subplot(2,2,3)
edges = [0 0:0.01:6 7];
h1 = histogram(x1,edges);
h1.Normalization = 'countdensity';
subplot(2,2,4)
edges = [0 0:0.01:6 7];
h8 = histogram(x08,edges);
h8.Normalization = 'countdensity';

