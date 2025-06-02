clear;
close all;

Nd = 1000000;
qs = [0.01 0.05 0.10 0.20];
X = [100 200 300 500 800 1000 2000 3000 5000 8000 10000 20000 30000 50000 80000 100000 200000 300000 500000];
eps_ind = 5;

for i = 1:4
    %     eval(['load simulation_results/mat_files/PeMoR1_q_' num2str(100*qs(i)) '.mat']);
    %     eval(['load simulation_results/mat_files/PeMoR2_q_' num2str(100*qs(i)) '.mat']);
    %     eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*qs(i)) '.mat']);
    %
    %     eval(['load simulation_results/mat_files/PeRoS1_q_' num2str(100*qs(i)) '.mat']);
    %     eval(['load simulation_results/mat_files/PeRoS2_q_' num2str(100*qs(i)) '.mat']);
    %     eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*qs(i)) '.mat']);
    %
    %     eval(['load simulation_results/mat_files/Pebound_q_' num2str(100*qs(i)) '.mat']);

    eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*qs(i)) '_p_15.mat']);
    ac_MoR15 = ac_MoR3;
    eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*qs(i)) '_p_30.mat']);
    ac_MoR30 = ac_MoR3;
    eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*qs(i)) '_p_50.mat']);
    ac_MoR50 = ac_MoR3;
    eval(['load simulation_results/mat_files/PeMoR3_q_' num2str(100*qs(i)) '_p_100.mat']);
    ac_MoR100 = ac_MoR3;

    eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*qs(i)) '_p_15.mat']);
    ac_RoS15 = ac_RoS3;
    eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*qs(i)) '_p_30.mat']);
    ac_RoS30 = ac_RoS3;
    eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*qs(i)) '_p_50.mat']);
    ac_RoS50 = ac_RoS3;
    eval(['load simulation_results/mat_files/PeRoS3_q_' num2str(100*qs(i)) '_p_100.mat']);
    ac_RoS100 = ac_RoS3;


    eval(['load simulation_results/mat_files/Pebound_q_' num2str(100*qs(i)) '_p_100.mat']);
    
    beta = 1 + eps_ind/100;
    d = 15;
    at_Mor15(:,i) = min(1, (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)) + X*(1 - (d/(Nd-1)))^(Nd-1));
    at_RoS15(:,i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)*d/2) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)*d/2));
    
    d = 30;
    at_Mor30(:,i) = min(1, (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)) + X*(1 - (d/(Nd-1)))^(Nd-1));
    at_RoS30(:,i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)*d/2) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)*d/2));

    d = 50;
    at_Mor50(:,i) = min(1, (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)) + X*(1 - (d/(Nd-1)))^(Nd-1));
    at_RoS50(:,i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)*d/2) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)*d/2));
    
    d = 100;
    at_Mor100(:,i) = min(1, (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)) + X*(1 - (d/(Nd-1)))^(Nd-1));
    at_RoS100(:,i) = min(1, (exp(1)/3.375).^(X*d/2) + (exp(beta - 1) / ((beta)^(beta))).^(X*qs(i)*d/2) + (exp((1/beta) - 1) / ((beta)^(-1/beta))).^(X*qs(i)*d/2));
 

    ssize(i) = (log(2) + 0.5*log(Nd)) / (qs(i) * (1 - (1/beta)*(log(beta) + 1)));

    MoR1(:,i) = ac_MoR15(:, eps_ind);
    MoR2(:,i) = ac_MoR30(:, eps_ind);
    MoR3(:,i) = ac_MoR50(:, eps_ind);
    MoR4(:,i) = ac_MoR100(:, eps_ind);

    RoS1(:,i) = ac_RoS15(:, eps_ind);
    RoS2(:,i) = ac_RoS30(:, eps_ind);
    RoS3(:,i) = ac_RoS50(:, eps_ind);
    RoS4(:,i) = ac_RoS100(:, eps_ind);


    s(:,i) = ac_s(:, eps_ind);
end

for i = 1:4
    yQuery = 0.05;
    [Y,ind]=unique(MoR1(:,i), 'stable');
    X1 = X(ind);
    y1(i) = interp1(Y,X1,yQuery,'pchip');
    x1(i) = qs(i);

    [Y,ind]=unique(MoR2(:,i), 'stable');
    X1 = X(ind);
    y2(i) = interp1(Y,X1,yQuery,'pchip');
    x2(i) = qs(i);

    [Y,ind]=unique(MoR3(:,i), 'stable');
    X1 = X(ind);
    y3(i) = interp1(Y,X1,yQuery,'pchip');
    x3(i) = qs(i);

    [Y,ind]=unique(MoR4(:,i), 'stable');
    X1 = X(ind);
    y4(i) = interp1(Y,X1,yQuery,'pchip');
    x4(i) = qs(i);



    [Y,ind]=unique(at_Mor15(:,i), 'stable');
    X1 = X(ind);
    y5(i) = interp1(Y,X1,yQuery,'pchip');
    x5(i) = qs(i);

    [Y,ind]=unique(at_Mor30(:,i), 'stable');
    X1 = X(ind);
    y6(i) = interp1(Y,X1,yQuery,'pchip');
    x6(i) = qs(i);

    [Y,ind]=unique(at_Mor50(:,i), 'stable');
    X1 = X(ind);
    y7(i) = interp1(Y,X1,yQuery,'pchip');
    x7(i) = qs(i);

    [Y,ind]=unique(at_Mor100(:,i), 'stable');
    X1 = X(ind);
    y8(i) = interp1(Y,X1,yQuery,'pchip');
    x8(i) = qs(i);


%     [Y,ind]=unique(s(:,i), 'stable');
%     X1 = X(ind);
%     y7(i) = interp1(Y,X1,yQuery,'pchip');
%     x7(i) = qs(i);
end

figure;
semilogy(x1,y1), hold on;
semilogy(x2,y2), hold on;
semilogy(x3,y3), hold on;
semilogy(x4,y4), hold on;
% semilogy(x5,y5), hold on;
semilogy(x6,y6), hold on;
semilogy(x7,y7), hold on;
semilogy(x8,y8), hold on;
semilogy(qs, ssize, '-.'), hold on;

data2save1(:,1) = [0.01 0.05 0.10 0.20]';
data2save1(:,2) = y1;
data2save1(:,3) = y2;
data2save1(:,4) = y3;
data2save1(:,5) = y4;
data2save1(:,6) = y5;
data2save1(:,7) = y6;
data2save1(:,8) = y7;
data2save1(:,9) = y8;
data2save1(:,10) = ssize';
eval(['save simulation_results/size_vs_rho_P_' num2str(yQuery*100) '.dat data2save1 -ascii'])



for i = 1:4
    yQuery = 0.05;
    [Y,ind]=unique(RoS1(:,i), 'stable');
    X1 = X(ind);
    y1(i) = interp1(Y,X1,yQuery,'pchip');
    x1(i) = qs(i);

    [Y,ind]=unique(RoS2(:,i), 'stable');
    X1 = X(ind);
    y2(i) = interp1(Y,X1,yQuery,'pchip');
    x2(i) = qs(i);

    [Y,ind]=unique(RoS3(:,i), 'stable');
    X1 = X(ind);
    y3(i) = interp1(Y,X1,yQuery,'pchip');
    x3(i) = qs(i);

    [Y,ind]=unique(RoS4(:,i), 'stable');
    X1 = X(ind);
    y4(i) = interp1(Y,X1,yQuery,'pchip');
    x4(i) = qs(i);



    [Y,ind]=unique(at_RoS15(:,i), 'stable');
    X1 = X(ind);
    y5(i) = interp1(Y,X1,yQuery,'pchip');
    x5(i) = qs(i);

    [Y,ind]=unique(at_RoS30(:,i), 'stable');
    X1 = X(ind);
    y6(i) = interp1(Y,X1,yQuery,'pchip');
    x6(i) = qs(i);

    [Y,ind]=unique(at_RoS50(:,i), 'stable');
    X1 = X(ind);
    y7(i) = interp1(Y,X1,yQuery,'pchip');
    x7(i) = qs(i);

    [Y,ind]=unique(at_RoS100(:,i), 'stable');
    X1 = X(ind);
    y8(i) = interp1(Y,X1,yQuery,'pchip');
    x8(i) = qs(i);


%     [Y,ind]=unique(s(:,i), 'stable');
%     X1 = X(ind);
%     y7(i) = interp1(Y,X1,yQuery,'pchip');
%     x7(i) = qs(i);
end

figure;
semilogy(x1,y1), hold on;
semilogy(x2,y2), hold on;
semilogy(x3,y3), hold on;
semilogy(x4,y4), hold on;
semilogy(x5,y5,'--'), hold on;
semilogy(x6,y6,'--'), hold on;
semilogy(x7,y7,'--'), hold on;
semilogy(x8,y8,'--'), hold on;
% semilogy(qs, ssize), hold on;

data2save1(:,1) = [0.01 0.05 0.10 0.20]';
data2save1(:,2) = y1;
data2save1(:,3) = y2;
data2save1(:,4) = y3;
data2save1(:,5) = y4;
data2save1(:,6) = y5;
data2save1(:,7) = y6;
data2save1(:,8) = y7;
data2save1(:,9) = y8;
data2save1(:,10) = ssize';
eval(['save simulation_results/size_vs_rho_P_' num2str(yQuery*100) 'RoS.dat data2save1 -ascii'])
