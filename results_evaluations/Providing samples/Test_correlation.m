xt=0:.00001:5; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

Il=1;Jl=1;Ir=1;Jr=1;

Nsample=1000;

for i=1:Nsample

[R,pv1]=corrcoef(Ldp(i,:),1-Order);
corrLdp(i)=R(2);
[R,pv2]=corrcoef(Rdp(i,:),1-Order);
corrRdp(i)=R(2);
[R,pv3]=corrcoef(LLp(i,:),1-Order);
corrLLp(i)=R(2);
[R,pv4]=corrcoef(RLp(i,:),1-Order);
corrRLp(i)=R(2);
[R,pv5]=corrcoef(Hpp(i,:),1-Order);
corrHpp(i)=R(2);
[R,pv6]=corrcoef(CurvL(i,:),1-Order);
corrCurvL(i)=R(2);
[R,pv7]=corrcoef(CurvR(i,:),1-Order);
corrCurvR(i)=R(2);


[RHO,pv8] = corr(CurvL(i,:)',1-Order','Type','Spearman');
corrCurvLSpear(i)=RHO;
[RHO,pv9] = corr(CurvR(i,:)',1-Order','Type','Spearman');
corrCurvRSpear(i)=RHO;

if corrRdp(i)>0
    pRelRdp(Ir)=abs(CurvR(i,1)-CurvR(i,end))/abs(CurvR(i,1));
    Ir=Ir+1;
else
    nRelRdp(Jr)=abs(CurvR(i,1)-CurvR(i,end))/abs(CurvR(i,1));
    Jr=Jr+1;
end

if corrLdp(i)>0
    pRelLdp(Il)=abs(CurvL(i,1)-CurvL(i,end))/abs(CurvL(i,1));
    Il=Il+1;
else
    nRelLdp(Jl)=abs(CurvL(i,1)-CurvL(i,end))/abs(CurvL(i,1));
    Jl=Jl+1;
end

pv(i)=mean([pv1(2),pv2(2),pv3(2),pv4(2),pv5(2),pv8,pv9]); %avarage of pvalues

end
