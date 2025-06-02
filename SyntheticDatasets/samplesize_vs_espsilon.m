clear;
close all;

q = 0.05;
d = 30;
Nd = 1000000

eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*q) '_p_15.mat']);
ac_MoR15 = ac_MoR3;
eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*q) '_p_30.mat']);
ac_MoR30 = ac_MoR3;
eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*q) '_p_50.mat']);
ac_MoR50 = ac_MoR3;
eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*q) '_p_100.mat']);
ac_MoR100 = ac_MoR3;


eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*q) '_p_15.mat']);
ac_RoS15 = ac_RoS3;
eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*q) '_p_30.mat']);
ac_RoS30 = ac_RoS3;
eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*q) '_p_50.mat']);
ac_RoS50 = ac_RoS3;
eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*q) '_p_100.mat']);
ac_RoS100 = ac_RoS3;

% eval(['load simulation_results/mat_files/Pebound_q_' num2str(100*q) '_p_100.mat']);

X = [100 200 300 500 800 1000 2000 3000 5000 8000 10000 20000 30000 50000 80000 100000 200000 300000 500000];

for i = 1:20
    ac_t(:, i) = max(0, 1 - (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) - (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) - X*(1 - (d/(Nd-1)))^(Nd-1));
    ac_s(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) + X*(1 - (d/(Nd-1)))^(Nd-1));

    ac_s15(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) + X*(1 - (15/(Nd-1)))^(Nd-1));
    ac_s30(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) + X*(1 - (30/(Nd-1)))^(Nd-1));
    ac_s50(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) + X*(1 - (50/(Nd-1)))^(Nd-1));
    ac_s100(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q) + X*(1 - (100/(Nd-1)))^(Nd-1));
    
    
    d = 15;
    ac_u15(:, i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q*d/2) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q*d/2));
    d = 30;
    ac_u30(:, i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q*d/2) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q*d/2));
    d = 50;
    ac_u50(:, i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q*d/2) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q*d/2));
    d = 100;
    ac_u100(:, i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q*d/2) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q*d/2));
    
    beta = 1 + (i/100);
    ssize(:, i) = ones(1,length(X))*(log(2) + 0.5 * log(Nd)) / (q * (1 - (1/beta)*(log(beta) + 1)));
end



for i = 2:10
    yQuery = 0.05;
    [Y,ind]=unique(ac_MoR15(:,i), 'stable');
    X1 = X(ind);
    y1(i-1) = interp1(Y,X1,yQuery,'pchip');
    x1(i-1) = i/100;

    [Y,ind]=unique(ac_MoR30(:,i), 'stable');
    X1 = X(ind);
    y2(i-1) = interp1(Y,X1,yQuery,'pchip');
    x2(i-1) = i/100;

    [Y,ind]=unique(ac_MoR50(:,i), 'stable');
    X1 = X(ind);
    y3(i-1) = interp1(Y,X1,yQuery,'pchip');
    x3(i-1) = i/100;

    [Y,ind]=unique(ac_MoR100(:,i), 'stable');
    X1 = X(ind);
    y4(i-1) = interp1(Y,X1,yQuery,'pchip');
    x4(i-1) = i/100;

    [Y,ind]=unique(ac_s15(:,i), 'stable');
    X1 = X(ind);
    y7(i-1) = interp1(Y,X1,yQuery,'pchip');
    x7(i-1) = i/100;

    [Y,ind]=unique(ac_s30(:,i), 'stable');
    X1 = X(ind);
    y8(i-1) = interp1(Y,X1,yQuery,'pchip');
    x8(i-1) = i/100;

    [Y,ind]=unique(ac_s50(:,i), 'stable');
    X1 = X(ind);
    y9(i-1) = interp1(Y,X1,yQuery,'pchip');
    x9(i-1) = i/100;

    [Y,ind]=unique(ac_s100(:,i), 'stable');
    X1 = X(ind);
    y10(i-1) = interp1(Y,X1,yQuery,'pchip');
    x10(i-1) = i/100;
end


epsl = linspace(0.02,0.1,9)';
for i = 1:length(epsl)
    beta = 1 + epsl(i);
    x11(i) = epsl(i);
    y11(i) = (log(2) + 0.5 * log(Nd)) / (q * (1 - (1/beta)*(log(beta) + 1)));
end


figure;
semilogy(x1,y1), hold on;
semilogy(x2,y2), hold on;
semilogy(x3,y3), hold on;
semilogy(x4,y4), hold on;
% semilogy(x7,y7), hold on;
semilogy(x8,y8), hold on;
semilogy(x9,y9), hold on;
semilogy(x10,y10), hold on;
semilogy(x11,y11,'-.'), hold on;
xlabel('\epsilon')
ylabel('|S|')


data2save1(:,1) = linspace(0.02,0.1,9)';
data2save1(:,2) = y1;
data2save1(:,3) = y2;
data2save1(:,4) = y3;
data2save1(:,5) = y4;
data2save1(:,6) = y7;
data2save1(:,7) = y8;
data2save1(:,8) = y9;
data2save1(:,9) = y10;
data2save1(:,10) = y11;
eval(['save simulation_results/size_vs_epsilon_q_' num2str(q*100) '_P_' num2str(yQuery*100) '.dat data2save1 -ascii'])


for i = 2:10
    yQuery = 0.05;
    [Y,ind]=unique(ac_RoS15(:,i), 'stable');
    X1 = X(ind);
    y1(i-1) = interp1(Y,X1,yQuery,'pchip');
    x1(i-1) = i/100;

    [Y,ind]=unique(ac_RoS30(:,i), 'stable');
    X1 = X(ind);
    y2(i-1) = interp1(Y,X1,yQuery,'pchip');
    x2(i-1) = i/100;

    [Y,ind]=unique(ac_RoS50(:,i), 'stable');
    X1 = X(ind);
    y3(i-1) = interp1(Y,X1,yQuery,'pchip');
    x3(i-1) = i/100;

    [Y,ind]=unique(ac_RoS100(:,i), 'stable');
    X1 = X(ind);
    y4(i-1) = interp1(Y,X1,yQuery,'pchip');
    x4(i-1) = i/100;

    [Y,ind]=unique(ac_u15(:,i), 'stable');
    X1 = X(ind);
    y7(i-1) = interp1(Y,X1,yQuery,'pchip');
    x7(i-1) = i/100;

    [Y,ind]=unique(ac_u30(:,i), 'stable');
    X1 = X(ind);
    y8(i-1) = interp1(Y,X1,yQuery,'pchip');
    x8(i-1) = i/100;

    [Y,ind]=unique(ac_u50(:,i), 'stable');
    X1 = X(ind);
    y9(i-1) = interp1(Y,X1,yQuery,'pchip');
    x9(i-1) = i/100;

    [Y,ind]=unique(ac_u100(:,i), 'stable');
    X1 = X(ind);
    y10(i-1) = interp1(Y,X1,yQuery,'pchip');
    x10(i-1) = i/100;
end

figure;
semilogy(x1,y1), hold on;
semilogy(x2,y2), hold on;
semilogy(x3,y3), hold on;
semilogy(x4,y4), hold on;
semilogy(x7,y7,'--'), hold on;
semilogy(x8,y8,'--'), hold on;
semilogy(x9,y9,'--'), hold on;
semilogy(x10,y10,'--'), hold on;
semilogy(x11,y11,'-.'), hold on;
xlabel('\epsilon')
ylabel('|S|')


data2save1(:,1) = linspace(0.02,0.1,9)';
data2save1(:,2) = y1;
data2save1(:,3) = y2;
data2save1(:,4) = y3;
data2save1(:,5) = y4;
data2save1(:,6) = y7;
data2save1(:,7) = y8;
data2save1(:,8) = y9;
data2save1(:,9) = y10;
data2save1(:,10) = y11;
eval(['save simulation_results/size_vs_epsilon_q_' num2str(q*100) '_P_' num2str(yQuery*100) 'RoS.dat data2save1 -ascii'])