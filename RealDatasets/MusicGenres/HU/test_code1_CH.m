clear
filename = 'HU_genres.json';
str = fileread(filename);
data = jsondecode(str);
M = csvread("HU_edges.csv", 1);
N = numel(fieldnames(data));

load all_genres.mat
genre = all_genres{32}

Nv = round([0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.1 0.2 0.5] * N);

%Nv = round([0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.075 0.1] * N);
trials = 1000;

h = zeros(N,1);
for i = 1:N
%     eval(['h(' num2str(i) ') = sum(strcmp(data.x' num2str(i-1) ',"Electro"));']);
    eval(['B = data.x' num2str(i-1) ';']);
    h(i) = sum(strcmp(B, genre));
%     eval(['h(' num2str(i) ') = sum(strcmp(data.x' num2str(i-1) ',' all_genres{1} '));']);
end


real_rate = sum(h)/length(h)

epsilon = [0.05, 0.1, 0.2];
sample_size = 0.1 * (log(N)/real_rate) * (1 ./ (epsilon.^2));

% graph1 = graph(M(:,1)+1, M(:,2)+1);
% adjacency_mat = adjacency(graph1);
adjacency_mat1 = sparse(M(:,1)+1, M(:,2)+1, ones(length(M(:,1)),1), N, N);
adjacency_mat2 = sparse(M(:,2)+1, M(:,1)+1, ones(length(M(:,1)),1), N, N);

adjacency_mat = (adjacency_mat1)|(adjacency_mat2);

RoS = zeros(length(Nv), trials);
MoR = zeros(length(Nv), trials);

ERoS = zeros(length(Nv), trials);
EMoR = zeros(length(Nv), trials);

Ht = sum(adjacency_mat, 2);
Gt = sum(adjacency_mat(:,h==1), 2);

for j = 1:length(Nv)
    j
    for i = 1:trials
        indices = randperm(N, Nv(j));
        H = Ht(indices);
        G = Gt(indices);

        RoS(j, i) = sum(G) / sum(H);
        MoR(j, i) = nanmean(G./H);

        ERoS(j, i) = max([RoS(j, i)/real_rate, real_rate/RoS(j, i)]);
        EMoR(j, i) = max([MoR(j, i)/real_rate, real_rate/MoR(j, i)]);
        
    end
end

close all;
SEM_RoS = std(RoS, 0, 2)/ sqrt(trials);
ts_RoS = tinv([0.025  0.975], trials - 1);
ci_RoS = ts_RoS(1)*SEM_RoS;

SEM_MoR = std(MoR, 0, 2)/ sqrt(trials);
ts_MoR = tinv([0.025  0.975], trials - 1);
ci_MoR = ts_MoR(1)*SEM_MoR;


SEM_ERoS = std(ERoS, 0, 2)/ sqrt(trials);
ts_ERoS = tinv([0.025  0.975], trials - 1);
ci_ERoS = ts_ERoS(1)*SEM_ERoS;

SEM_EMoR = std(EMoR, 0, 2)/ sqrt(trials);
ts_EMoR = tinv([0.025  0.975], trials - 1);
ci_EMoR = ts_EMoR(1)*SEM_EMoR;



for eps = 1:20
    beta = 1 + eps/100;
    delta = 0.5;
    gamma=2.5;

    PeMoR_Bound(:,eps) = min(1, (exp(beta-1)/(beta^beta)).^(Nv*real_rate) ...
        + (exp((1/beta) - 1)/(beta^(-1/beta))).^(Nv*real_rate));

    mu = Nv * ((1-gamma)*(1-(N - 1)^(2-gamma))) / ((2-gamma)*(1-(N - 1)^(1-gamma)));
    PeRoS_Bound(:,eps) = min(1, (exp(-delta)/((1-delta)^(1-delta))).^(mu) ...
        + (exp(beta-1) / (beta^(beta))).^(mu*real_rate/2) ...
        + (exp((1/beta)-1) / ((beta)^(-1/beta))).^(mu*real_rate/2));

    PeRoS(:,eps) = mean(ERoS > 1 + eps/100, 2);
    PeMoR(:,eps) = mean(EMoR > 1 + eps/100, 2);
end

figure;
semilogx(Nv, PeMoR(:,10), 'LineWidth',2), hold on;
semilogx(Nv, PeRoS(:,10), 'LineWidth',2);
semilogx(Nv, PeMoR_Bound(:,10),'--','LineWidth',2);
semilogx(Nv, PeRoS_Bound(:,10),'--','LineWidth',2);
xlabel('Sample Size |S|')
ylabel('Perror > 1 + epsilon')
legend('MoR', 'RoS','Bound MoR', 'Bound RoS')


figure;
errorbar(Nv, mean(RoS,2), ci_RoS(:,1), '--o','LineWidth', 2);
hold on;
errorbar(Nv, mean(MoR,2), ci_MoR(:,1), 'k--o','LineWidth', 2);
hold on;
plot(Nv, ones(length(Nv))*real_rate, 'r--','LineWidth', 2);
ax = gca;
ax.FontSize = 14;
% ylim([0 0.1])
xlabel('Number of  of surveys','FontSize',14);
ylabel('Estimation','FontSize',14)
title(['|G|:' num2str(N) '. Trials: ' num2str(trials)],'FontSize',14)
legend('Ratio of sums', 'Mean of ratios', 'True value','FontSize',14)
grid on;
hold off;


figure;
errorbar(Nv, mean(ERoS,2), ci_ERoS(:,1), '--o','LineWidth', 2);
hold on;
errorbar(Nv, mean(EMoR,2), ci_EMoR(:,1), 'k--o','LineWidth', 2);
ax = gca;
ax.FontSize = 14;
% ylim([0 0.1])
xlabel('Number of  of surveys','FontSize',14);
ylabel('Estimation','FontSize',14)
title(['|G|:' num2str(N) '. Trials: ' num2str(trials)],'FontSize',14)
legend('Ratio of sums', 'Mean of Ratios', 'FontSize',14)
grid on;
hold off;