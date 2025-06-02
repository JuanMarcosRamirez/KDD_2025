clear;
close all;
trials = 20;
Nd = 250000; % no. of individuals the participants can cover

q = 0.20;
num_mat = 50;
num_hidpop = 10;
X = [100 200 300 500 800 1000 2000 3000 5000 8000 10000 20000 30000 50000 80000 100000 125000];

qRoS1 = zeros(length(X), num_mat*num_hidpop*trials);
qRoS2 = zeros(length(X), num_mat*num_hidpop*trials);
qRoS3 = zeros(length(X), num_mat*num_hidpop*trials);
qMoR1 = zeros(length(X), num_mat*num_hidpop*trials);
qMoR2 = zeros(length(X), num_mat*num_hidpop*trials);
qMoR3 = zeros(length(X), num_mat*num_hidpop*trials);

for ii = 1:num_mat
    ii
    eval(['load synthetic_data/q_' num2str(q*100) '/Gwtsr_v_' num2str(ii) '.mat'])
    eval(['load synthetic_data/q_' num2str(q*100) '/Gwosr_v_' num2str(ii) '.mat'])
    for jj = 1:num_hidpop
        eval(['load synthetic_data/q_' num2str(q*100) '/H1_v_' num2str(ii) '_t_' num2str(jj) '.mat'])
        eval(['load synthetic_data/q_' num2str(q*100) '/H2_v_' num2str(ii) '_t_' num2str(jj) '.mat'])
        eval(['load synthetic_data/q_' num2str(q*100) '/H3_v_' num2str(ii) '_t_' num2str(jj) '.mat'])
        H1 = full(H1);
        H2 = full(H2);
        H3 = full(H3);
        for n = 1:length(X)
            N = X(n);
            for i = 1:trials
                rp = randperm(Nd,N);
                G1s = Gwtsr(rp);
                G2s = Gwosr(rp);
                
%                 gamma = 1 + length(Gwosr)/sum(log(Gwosr/(min(Gwosr)-0.5)))

                H1s = H1(rp);
                H2s = H2(rp);
                H3s = H3(rp);

                qRoS1(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = sum(H1s(:)) / sum(G1s(:));
                qRoS2(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = sum(H2s(:)) / sum(G1s(:));
                qRoS3(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = sum(H3s(:)) / sum(G2s(:));

                qMoR1(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = nanmean(H1s./G1s);
                qMoR2(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = nanmean(H2s./G1s);
                qMoR3(n, num_hidpop*trials*(ii - 1) + trials*(jj - 1) + i) = nanmean(H3s./G2s);
            end
        end
    end
end

for n = 1:length(X)
    eM_RoS1(n,:) = max([q./full(qRoS1(n,:)); full(qRoS1(n,:))/q]);
    eM_RoS2(n,:) = max([q./full(qRoS2(n,:)); full(qRoS2(n,:))/q]);
    eM_RoS3(n,:) = max([q./full(qRoS3(n,:)); full(qRoS3(n,:))/q]);
    
    eM_MoR1(n,:) = max([q./full(qMoR1(n,:)); full(qMoR1(n,:))/q]);
    eM_MoR2(n,:) = max([q./full(qMoR2(n,:)); full(qMoR2(n,:))/q]);
    eM_MoR3(n,:) = max([q./full(qMoR3(n,:)); full(qMoR3(n,:))/q]);
end

errorM_RoS1 = mean(eM_RoS1,2);
errorM_RoS2 = mean(eM_RoS2,2);
errorM_RoS3 = mean(eM_RoS3,2);
errorM_MoR1 = mean(eM_MoR1,2);
errorM_MoR2 = mean(eM_MoR2,2);
errorM_MoR3 = mean(eM_MoR3,2);


indices=[1,4,6,9,11,14];
colors = {[0.3010 0.7450 0.9330] [0.4660 0.6740 0.1880] };
GroupedData = {eM_MoR3(indices,:)' eM_RoS3(indices,:)'};
figure;
hh1 = boxplot(GroupedData{1}, X(indices),'Colors',colors{1}, 'Widths',0.3, 'Positions',(1:length(X(indices)))-0.2);
hold on;
hh2 = boxplot(GroupedData{2}, 'Colors',colors{2}, 'Widths',0.3, 'Positions',(1:length(X(indices)))+0.2, 'Labels',X(indices), 'Symbol','+m');
set(hh1,'LineWidth',2, 'LineStyle','-', 'MarkerEdgeColor',"#7E2F8E")
set(hh2,'LineWidth',2, 'LineStyle','-', 'MarkerEdgeColor',"#7E2F8E")
% legend('MoR', 'RoS')
legend([hh1(5,1),hh2(5,1)], {'MoR','RoS'},'FontName','Times','FontSize',18)
set(gca,'LineWidth',2)
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times';
ax.YGrid = 'on'
pbaspect([1.25 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))
ylim([0.98,1.5])
ylabel(strcat("$", "\mathcal{E}_{\mathcal{M}}", "$"),'FontSize', 22, 'Interpreter', 'latex')
xlabel("Sample Size $|S|$",'FontSize', 22, 'Interpreter','latex')
title(strcat('(a) $\rho =~$', num2str(q)), 'FontSize', 22, 'Interpreter','latex')
ytickformat('%.2f')
print('-depsc', strcat('EpsFigs/SF_hist_rho_', num2str(100*q)))
exportgraphics(gcf,strcat('EpsFigs/SF_hist_rho_', num2str(100*q),'.png'),"Resolution",300)


data2save1(:,1) = X';
data2save1(:,2) = errorM_RoS1';
data2save1(:,3) = errorM_RoS2';
data2save1(:,4) = errorM_RoS3';
data2save1(:,5) = errorM_MoR1';
data2save1(:,6) = errorM_MoR2';
data2save1(:,7) = errorM_MoR3';
eval(['save simulation_results/error_trials_' num2str(trials*num_mat*num_hidpop) '_q_' num2str(100*q) '.dat data2save1 -ascii'])


figure;
subplot(131)
semilogx(X, errorM_MoR1)
hold on
semilogx(X, errorM_RoS1)
legend('MoR', 'RoS')
xlabel('Sample Size |S|')
ylabel('\epsilon')
title(['With replacement and with self-reporting. q = ' num2str(q)])

subplot(132)
semilogx(X, errorM_MoR2)
hold on
semilogx(X, errorM_RoS2)
legend('MoR', 'RoS')
xlabel('Sample Size |S|')
ylabel('\epsilon')
title(['Without replacement and with self-reporting. q = ' num2str(q)])

subplot(133)
semilogx(X, errorM_MoR3)
hold on
semilogx(X, errorM_RoS3)
legend('MoR', 'RoS')
xlabel('Sample Size |S|')
ylabel('\epsilon')
title(['Without replacement and without self-reporting. q = ' num2str(q)])


for n = 1:length(X)
    for i = 1:20
        ac_MoR1(n, i) = mean(eM_MoR1(n,:) > 1 + (i/100));
        ac_MoR2(n, i) = mean(eM_MoR2(n,:) > 1 + (i/100));
        ac_MoR3(n, i) = mean(eM_MoR3(n,:) > 1 + (i/100));
        
        ac_RoS1(n, i) = mean(eM_RoS1(n,:) > 1 + (i/100));
        ac_RoS2(n, i) = mean(eM_RoS2(n,:) > 1 + (i/100));
        ac_RoS3(n, i) = mean(eM_RoS3(n,:) > 1 + (i/100));
    end
end

gamma = 1.9375;
nu = 1 / sum((1:Nd).^(-gamma));
[gamma nu]

for i = 1:20
    ac_t(:, i) = max(0, 1 - (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) - (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q));
    ac_s(:, i) = min(1, (exp(i/100) / ((1+i/100)^(1+i/100))).^(X*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(X*q));
    gamma=2.2;
    delta = 0.5;
    mu = X * ((1-gamma)*(1-(Nd - 1)^(2-gamma))) / ((2-gamma)*(1-(Nd - 1)^(1-gamma)))
    ac_u(:, i) = min(1, (exp(-delta)/((1-delta)^(1-delta))).^(mu) + (exp(i/100) / ((1+i/100)^(1+i/100))).^(mu*q/2) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(mu*q/2));
    
%     for j = i:length(X)
%         R = X(j):2*X(j);
%         [X(j) min(R) max(R)]
%         for k = 1:length(R)
%             vcn(k) = nchoosek(R(k)-1, X(j)-1);
%         end
%         ac_u(j,i) = min(1, sum(((exp(i/100) / ((1+i/100)^(1+i/100))).^(R*q) ...
%             + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(R*q)) ...
%             .*(nu^(X(j))).*((R/X(j)).^(-gamma*X(j))).*vcn ...
%             ));
% %         ac_u(j,i) = min(1, sum(((exp(i/100) / ((1+i/100)^(1+i/100))).^(R*q) + (exp(-(i/100)/(1+i/100)) / ((1/(1+i/100))^(1/(1+i/100)))).^(R*q)).*(nu^(X(j))).*((R/X(j))^(-gamma * X(j))).*nchoosek(R-1,X(j)-1)));
%     end
end


figure;
subplot(131)
semilogx(X, ac_MoR1(:,5), 'LineWidth',2); hold on;
semilogx(X, ac_RoS1(:,5), 'LineWidth',2);
semilogx(X, ac_s(:,5), 'LineWidth',2);
xlabel('Sample Size')
ylabel('P(eM > 1 + \epsilon )')
legend('MoR', 'RoS', 'Bound')
title(['W-replace. W-self-rep. \epsilon = 0.05. q = ' num2str(q)])
subplot(132)
semilogx(X, ac_MoR2(:,5), 'LineWidth',2); hold on;
semilogx(X, ac_RoS2(:,5), 'LineWidth',2);
semilogx(X, ac_s(:,5), 'LineWidth',2);
xlabel('Sample Size')
ylabel('P(eM > 1 + \epsilon )')
legend('MoR', 'RoS', 'Bound')
title(['W/O-replace. W-self-rep.\epsilon = 0.05. q = ' num2str(q)])
subplot(133)
semilogx(X, ac_MoR3(:,5), 'LineWidth',2); hold on;
semilogx(X, ac_RoS3(:,5), 'LineWidth',2);
semilogx(X, ac_s(:,5), 'LineWidth',2);
xlabel('Sample Size')
ylabel('P(eM > 1 + \epsilon )')
legend('MoR', 'RoS', 'Bound')
title(['W/O-replace. W/O-self-rep.\epsilon = 0.05. q = ' num2str(q)])


data2save2(:,1) = X';
data2save2(:,2) = ac_MoR1(:,5);
data2save2(:,3) = ac_MoR2(:,5);
data2save2(:,4) = ac_MoR3(:,5);
data2save2(:,5) = ac_RoS1(:,5);
data2save2(:,6) = ac_RoS2(:,5);
data2save2(:,7) = ac_RoS3(:,5);
data2save2(:,8) = ac_s(:,5);
data2save2(:,9) = ac_u(:,5);
eval(['save simulation_results/Pe_trials_' num2str(trials*num_mat*num_hidpop) '_q_' num2str(100*q) 'eps_5.dat data2save2 -ascii'])

data2save2(:,1) = X';
data2save2(:,2) = ac_MoR1(:,10);
data2save2(:,3) = ac_MoR2(:,10);
data2save2(:,4) = ac_MoR3(:,10);
data2save2(:,5) = ac_RoS1(:,10);
data2save2(:,6) = ac_RoS2(:,10);
data2save2(:,7) = ac_RoS3(:,10);
data2save2(:,8) = ac_s(:,10);
data2save2(:,9) = ac_u(:,10);
eval(['save simulation_results/Pe_trials_' num2str(trials*num_mat*num_hidpop) '_q_' num2str(100*q) '_eps_10.dat data2save2 -ascii'])

eval(['save simulation_results/mat_files/sample_size_q_' num2str(100*q) '.mat X'])
eval(['save simulation_results/mat_files/PeMoR1_q_' num2str(100*q) '.mat ac_MoR1'])
eval(['save simulation_results/mat_files/PeMoR2_q_' num2str(100*q) '.mat ac_MoR2'])
eval(['save simulation_results/mat_files/PeMoR3_q_' num2str(100*q) '.mat ac_MoR3'])
eval(['save simulation_results/mat_files/PeRoS1_q_' num2str(100*q) '.mat ac_RoS1'])
eval(['save simulation_results/mat_files/PeRoS2_q_' num2str(100*q) '.mat ac_RoS2'])
eval(['save simulation_results/mat_files/PeRoS3_q_' num2str(100*q) '.mat ac_RoS3'])
eval(['save simulation_results/mat_files/Pebound_q_' num2str(100*q) '.mat ac_s'])