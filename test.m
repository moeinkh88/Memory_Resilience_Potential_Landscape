xt=-2:.001:2;
for i=1:500
true=0;


while true==0
par=[-rand(1),randi([-1,1],1,3).*rand(1,3)]; % Randomly chosen coefficients of the polynomial     
U=-polyint(par); % potential of polynomial fun; -int(poly)
U1=polyval(U,xt);
TF=islocalmin(U1);
TF1=islocalmax(U1);
Depth=abs(U1(TF1)-U1(TF));

if length(Depth)==2  
    if Depth(1)>.05 && Depth(2)>.05
%     p1=plot(xt,polyval(U,xt));
    true=1;
    
    LRH(i,:)=U1(TF);
    end
end
end

end
%%
figure
subplot(2,2,[1,2])
g1 = repmat({'Left'},length(LRH(:,1)),1);
g2 = repmat({'right'},length(LRH(:,2)),1);
boxplot([LRH(:,1);LRH(:,2)],[g1;g2])

subplot(2,2,3)

histogram(LRH(:,1))

subplot(2,2,4)

histogram(LRH(:,2))