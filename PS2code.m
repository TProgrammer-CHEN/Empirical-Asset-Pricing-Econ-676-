clear;

raw = importdata('Problem_Set2.xls');
%% Means, Standard deviations and Covariance Matrix of Industries
industry_ret = raw.Industry_returns(:,2:end-1);
Means = mean((industry_ret));
Means = Means';
St_deviations = std((industry_ret));
Covar_matrix = cov((industry_ret));
%% Construction of Minimum Variance Portfolio
Ones = ones(10,1); 
weight_MVP = inv(Covar_matrix)*Ones/(Ones'*inv(Covar_matrix)*Ones);
Means_MVP = weight_MVP'*Means;
Stdev_MVP = sqrt(1/(Ones'*inv(Covar_matrix)*Ones));

%% Construction of Tangency Portfolio

rf = mean(raw.Industry_returns(:,end));
weight_TangP = inv(Covar_matrix)*(Means - Ones*rf)/(Ones'*inv(Covar_matrix)*(Means-Ones*rf));
Means_TangP = weight_TangP'*Means;
Stdev_TangP = sqrt(weight_TangP'*Covar_matrix*weight_TangP); 

%% Plot of Efficient Frontier with MVP and Tangent Portfolios

NumPorts = 10;
AssetBounds = [-0.7, 1];

LowerBound = AssetBounds(1);
UpperBound = AssetBounds(2);

p = Portfolio;
p = setAssetMoments(p, Means, Covar_matrix);
p = setDefaultConstraints(p);
p = setBounds(p, LowerBound, UpperBound);
plotFrontier(p, NumPorts);
hold on
labels = {'NoDur','Durbl','Manuf','Enrgy','HiTec','Telcm','Shops','Hlth','Utils','Other'};
plot(St_deviations, Means,'*r');
text(St_deviations, Means,labels,'VerticalAlignment','top','HorizontalAlignment','left')
hold on
sz = 100;
scatter(Stdev_MVP,Means_MVP,sz,'b','filled');
scatter(Stdev_TangP,Means_TangP,sz,'r','d','filled');
legend('Efficient Frontier','10 Inmdustries','Minimum Variance Portfolio','Tangency Portfolio', 'Location','northwest');
lgd.Fontsize = 16;
title('Efficient Frontier with MVP and Tangent Portfolios');
hold off
%% Question 1B: Increasing Mean Returns by 1 Standard error
Means2 = zeros(10,1);
for i=1:10
    Means2(i) = Means(i) + St_deviations(i)/sqrt(t); 
end

weight_TangP2 = inv(Covar_matrix)*(Means2 - Ones*rf)/(Ones'*inv(Covar_matrix)*(Means2-Ones*rf));
Means_TangP2 = weight_TangP2'*Means2;
Stdev_TangP2 = sqrt(weight_TangP2'*Covar_matrix*weight_TangP2); 

Means_MVP2 = weight_MVP'*Means2;
%% New Plot of Efficient Frontier with MVP and Tangent Portfolios

NumPorts = 10;
AssetBounds = [-0.7, 1];

LowerBound = AssetBounds(1);
UpperBound = AssetBounds(2);

p2 = Portfolio;
p2 = setAssetMoments(p2, Means2, Covar_matrix);
p2 = setDefaultConstraints(p2);
p2 = setBounds(p2, LowerBound, UpperBound);
plotFrontier(p2, NumPorts);
hold on
labels = {'NoDur','Durbl','Manuf','Enrgy','HiTec','Telcm','Shops','Hlth','Utils','Other'};
plot(St_deviations, Means2,'*r');
text(St_deviations, Means2,labels,'VerticalAlignment','top','HorizontalAlignment','left')
hold on
sz = 100;
scatter(Stdev_MVP,Means_MVP2,sz,'b','filled');
scatter(Stdev_TangP2,Means_TangP2,sz,'r','d','filled');
legend('Efficient Frontier','10 Industries','Minimum Variance Portfolio','Tangency Portfolio', 'Location','northwest');
lgd.Fontsize = 16;
title('Efficient Frontier with MVP and Tangent Portfolios (Mean Returns increased by 1 Standard Error)');
hold off

%% Question 1C: Change in Covariance Matrix changes MVP, Tangency portfolios and Efficient Frontier
Covar_matrix1C = diag(diag(Covar_matrix));

weight_MVP1C = inv(Covar_matrix1C)*Ones/(Ones'*inv(Covar_matrix1C)*Ones);
Means_MVP1C = weight_MVP1C'*Means;
Stdev_MVP1C = sqrt(1/(Ones'*inv(Covar_matrix1C)*Ones));


weight_TangP1C = inv(Covar_matrix1C)*(Means - Ones*rf)/(Ones'*inv(Covar_matrix1C)*(Means-Ones*rf));
Means_TangP1C = weight_TangP1C'*Means;
Stdev_TangP1C = sqrt(weight_TangP1C'*Covar_matrix1C*weight_TangP1C); 


%% Question 1C: Plot of Efficient Frontier with MVP and Tangent Portfolios
NumPorts = 10
AssetBounds = [0, 1];

LowerBound = AssetBounds(1);
UpperBound = AssetBounds(2);

p3 = Portfolio;
p3 = setAssetMoments(p3, Means, Covar_matrix1C);
p3 = setDefaultConstraints(p3);
p3 = setBounds(p3, LowerBound, UpperBound);
plotFrontier(p3, NumPorts);
hold on
labels = {'NoDur','Durbl','Manuf','Enrgy','HiTec','Telcm','Shops','Hlth','Utils','Other'};
plot(St_deviations, Means,'*r');
text(St_deviations, Means,labels,'VerticalAlignment','top','HorizontalAlignment','left')
hold on
sz = 100;
scatter(Stdev_MVP1C,Means_MVP1C,sz,'b','filled');
scatter(Stdev_TangP1C,Means_TangP1C,sz,'r','d','filled');
legend('Efficient Frontier','10 Industries','Minimum Variance Portfolio','Tangency Portfolio', 'Location','southwest');
lgd.Fontsize = 16;
title('Efficient Frontier with MVP and Tangent Portfolios (Diagonal Matrix of Variances as the Covariance Matrix)');
hold off

%% Identity Matrix as a Covariance Matrix
Covar_matrix1C_2 = diag(Ones);

weight_MVP1C_2 = inv(Covar_matrix1C_2)*Ones/(Ones'*inv(Covar_matrix1C_2)*Ones);
Means_MVP1C_2 = weight_MVP1C_2'*Means;
Stdev_MVP1C_2 = sqrt(1/(Ones'*inv(Covar_matrix1C_2)*Ones));


weight_TangP1C_2 = inv(Covar_matrix1C_2)*(Means - Ones*rf)/(Ones'*inv(Covar_matrix1C_2)*(Means-Ones*rf));
Means_TangP1C_2 = weight_TangP1C_2'*Means;
Stdev_TangP1C_2 = sqrt(weight_TangP1C_2'*Covar_matrix1C_2*weight_TangP1C_2); 

%% Question 1C part 2: Plot of Efficient Frontier with MVP and Tangent Portfolios

NumPorts = 10;
AssetBounds = [0, 1];

LowerBound = AssetBounds(1);
UpperBound = AssetBounds(2);

p4 = Portfolio;
p4 = setAssetMoments(p4, Means, Covar_matrix1C_2);
p4 = setDefaultConstraints(p4);
plotFrontier(p4, NumPorts);
hold on
labels = {'NoDur','Durbl','Manuf','Enrgy','HiTec','Telcm','Shops','Hlth','Utils','Other'};
plot(St_deviations, Means,'*r');
text(St_deviations, Means,labels,'VerticalAlignment','top','HorizontalAlignment','left')
hold on
sz = 100;
scatter(Stdev_MVP1C_2,Means_MVP1C_2,sz,'b','filled');
scatter(Stdev_TangP1C_2,Means_TangP1C_2,sz,'r','d','filled');
legend('Efficient Frontier','10 Industries','Minimum Variance Portfolio','Tangency Portfolio', 'Location','southwest');
lgd.Fontsize = 16;
title('Efficient Frontier with MVP and Tangent Portfolios (Identity Matrix as the Covariance Matrix)');
hold off

%%poblem d computation
rmvplist=[Means_MVP];smvplist=[Stdev_MVP];
rtglist=[Means_TangP];stglist=[Stdev_TangP];
t=size(industry_ret,1);
for i = 1:1000
    rets = mvnrnd(Means,Covar_matrix,t);
    rmean=mean(rets);rsd=std((rets));rcov = cov((rets));
    rmean=rmean';
    wmvp=inv(rcov)*Ones/(Ones'*inv(rcov)*Ones);
    rwmvp=wmvp'*Means;
    stdwmvp=sqrt(wmvp'*Covar_matrix*wmvp);
    rmvplist=[rmvplist;rwmvp];smvplist=[smvplist;stdwmvp];
    
    wtg=inv(rcov)*(rmean-Ones*rf)/(Ones'*inv(rcov)*(rmean-Ones*rf));
    rwtg=wtg'*Means;
    stdwtg=sqrt(wtg'*Covar_matrix*wtg);
    rtglist=[rtglist;rwtg];stglist=[stglist;stdwtg]; 
end
%% problem d graph plot 1
hold on
scatter(smvplist,rmvplist);
scatter(smvplist(1),rmvplist(1),'red','filled');
title('STD-MEAN Plot of MVP with simulation of Normal Distribution');
legend('Simulation','Real');
xlabel('Standard Deviation');ylabel('Expected Return');
hold off
%% problem d graph plot 2
hold on
scatter(stglist,rtglist);
scatter(stglist(1),rtglist(1),'red','filled');
title('STD-MEAN Plot of tangency portfolios with simulation of Normal Distribution');
legend('Simulation','Real');
xlabel('Standard Deviation');ylabel('Expected Return');
hold off

%% problem e computation
rmvplist=[Means_MVP];smvplist=[Stdev_MVP];
rtglist=[Means_TangP];stglist=[Stdev_TangP];
t=size(industry_ret,1);
for i = 1:1000
    rets =datasample(industry_ret,t);
    rmean=mean(rets);rsd=std((rets));rcov = cov((rets));
    rmean=rmean';
    wmvp=inv(rcov)*Ones/(Ones'*inv(rcov)*Ones);
    rwmvp=wmvp'*Means;
    stdwmvp=sqrt(wmvp'*Covar_matrix*wmvp);
    rmvplist=[rmvplist;rwmvp];smvplist=[smvplist;stdwmvp];
    
    wtg=inv(rcov)*(rmean-Ones*rf)/(Ones'*inv(rcov)*(rmean-Ones*rf));
    rwtg=wtg'*Means;
    stdwtg=sqrt(wtg'*Covar_matrix*wtg);
    rtglist=[rtglist;rwtg];stglist=[stglist;stdwtg]; 
end

%% problem e graph plot 1
hold on
scatter(smvplist,rmvplist);
scatter(smvplist(1),rmvplist(1),'red','filled');
title('STD-MEAN Plot of MVP with simulation of Empirical Distribution');
legend('Simulation','Real');
xlabel('Standard Deviation');ylabel('Expected Return');
hold off
%% problem e graph plot 2
hold on
scatter(stglist,rtglist);
scatter(stglist(1),rtglist(1),'red','filled');
title('STD-MEAN Plot of tangency portfolios with simulation of Empirical Distribution');
legend('Simulation','Real');
xlabel('Standard Deviation');ylabel('Expected Return');
hold off