stock=table2array(ProblemSet1);
portfolio=[mean(stock(:,1:5),2),mean(stock(:,1:10),2),mean(stock(:,1:25),2),mean(stock(:,1:50),2)];
expc=mean(portfolio);
stdp=std(portfolio);
num=[5,10,25,50];
plot(num,stdp);

vartotal=var(portfolio);
varall=var(stock);
varind=[(1/25)*sum(varall(1:5)),(1/100)*sum(varall(1:10)),(1/(25*25))*sum(varall(1:25)),(1/(50*50))*sum(varall(1:50))];
stdcov=stdp-stdind;
varind./vartotal;
plot(num,stdind./stdp);

% tstat=expc./(stdp/180^0.5)
[h,p,ci,stats] = ttest(portfolio);

ttest(stock(:,1));
%random normal
x=randn(1,18000);
(max(x)-min(x))/(std(x));
skewness(x);
kurtosis(x)-3;
[(max(x)-min(x))/(std(x)),skewness(x),kurtosis(x)];

%stock
(max(stock(:,1))-min(stock(:,1)))/(stdall(1));
skewness(stock(:,1));
kurtosis(stock(:,1))-3;

%portfolio
(max(portfolio(:,4))-min(stock(:,4)))/(stdp(4));
skewness(portfolio(:,4));
kurtosis(portfolio(:,4))-3;

%index
index=table2array(ProblemSet2);
(max(index)-min(index))/(std(index));
skewness(index);
kurtosis(index)-3;

%Jarque–Bera test
%(180/6)*(skewness(x)^2+0.25*(kurtosis(x)-3)^2);
jbtest(stock(:,1));
jbtest(portfolio(:,4));
jbtest(index);


A=[];rs=[];

for i=1:10
    mdl=fitlm(index,stock(:,i))
    A=[A;(mdl.Coefficients.Estimate).']
    rs=[rs;mdl.Rsquared.Ordinary]
end

time =table2array(ProblemSet3);

volt1=[];volt2=[];
for i=1:168
    id=i+12
    volt1=[volt1;[std(t(1:id,1)),std(t(1:id,2))]]
    volt2=[volt2;[std(t(id-12:id,1)),std(t(id-12:id,2))]]
end

volt=[volt1,volt2];
tp=time(13:180);
tp=tp*100+1;
tdn=datenum(num2str(tp),'yyyymmdd');
plot(tdn,volt);datetick('x','yyyy-mm');
legend('TXN by method1','Market by method1','TXN by method2','Market by method2')

beta=[];
for i=1:168
    id=i+12
    [b1,bint1]=regress(t(1:id,1),t(1:id,2))
    [b2,bint2]=regress(t(id-11:id,1),t(id-11:id,2))
    beta=[beta;[b1,bint1,b2,bint2]]
end
plot(tdn,beta);datetick('x','yyyy-mm');
legend('beta method1','lowerbound method1','upperbound method1','beta method2','lowerbound method2','upperbound method2')