clear
close all;
filename = 'RO_genres.json';
str = fileread(filename);
data = jsondecode(str);
M = csvread("RO_edges.csv", 1);
N = numel(fieldnames(data));
load all_genres.mat

genre_ind = [80, 30, 40, 43, 44, 63, 32, 22, 62];
Nv = round([0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.1 0.2 0.5] * N);
trials = 10000;

for qq = 1:length(genre_ind)
    genre = all_genres{genre_ind(qq)};
    h = zeros(N,1);
    for i = 1:N
        eval(['B = data.x' num2str(i-1) ';']);
        h(i) = sum(strcmp(B, genre));
    end

    real_rate = sum(h)/length(h)
    rr_vec(qq) = real_rate;

    adjacency_mat = sparse(M(:,1)+1, M(:,2)+1, ones(length(M(:,1)),1), N, N);
    RoS = zeros(length(Nv), trials);
    MoR = zeros(length(Nv), trials);

    ERoS = zeros(length(Nv), trials);
    EMoR = zeros(length(Nv), trials);
    Rs = zeros(length(Nv), trials);

    Ht = sum(adjacency_mat, 2);
    Gt = sum(adjacency_mat(:,h==1), 2);

    gamma = 1 + N / sum(log((Ht+1)/2))

    for j = 1:length(Nv)
        [qq, j]
        for i = 1:trials
            indices = randperm(N, Nv(j));
            H = Ht(indices);
            G = Gt(indices);

            Rs(j, i) = sum(H);

            RoS(j, i) = sum(G) / sum(H);
            MoR(j, i) = nanmean(G./H);

            ERoS(j, i) = max([RoS(j, i)/real_rate, real_rate/RoS(j, i)]);
            EMoR(j, i) = max([MoR(j, i)/real_rate, real_rate/MoR(j, i)]);
        end
    end

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

    for j = 1:length(Nv)
        [ht, bins] = histcounts(Rs(j,:), linspace(min(Rs(j,:)), max(Rs(j,:)), max(Rs(j,:)) - min(Rs(j,:)) + 1));
        Prob{j} = ht/trials;
        Rv{j} = bins(end-1);
    end

    estimation_plot_data(:, 4*(qq-1) + 1) = mean(MoR, 2);
    estimation_plot_data(:, 4*(qq-1) + 2) = abs(ci_MoR);
    estimation_plot_data(:, 4*(qq-1) + 3) = mean(RoS, 2);
    estimation_plot_data(:, 4*(qq-1) + 4) = abs(ci_RoS);

    error_plot_data(:, 4*(qq-1) + 1) = mean(EMoR, 2);
    error_plot_data(:, 4*(qq-1) + 2) = abs(ci_EMoR);
    error_plot_data(:, 4*(qq-1) + 3) = mean(ERoS, 2);
    error_plot_data(:, 4*(qq-1) + 4) = abs(ci_ERoS);
    


    for eps = 1:20
        beta = 1 + eps/100;
        delta = 0.5;
%         gamma=2.5;

        PeMoR_Bound(:,eps) = min(1, (exp(beta-1)/(beta^beta)).^(Nv*real_rate) ...
            + (exp((1/beta) - 1)/(beta^(-1/beta))).^(Nv*real_rate));

        mu = Nv * ((1-gamma)*(1-(N - 1)^(2-gamma))) / ((2-gamma)*(1-(N - 1)^(1-gamma)));
        PeRoS_Bound2(:,eps) = min(1, (exp(-delta)/((1-delta)^(1-delta))).^(mu) ...
            + (exp(beta-1) / (beta^(beta))).^(mu*real_rate/2) ...
            + (exp((1/beta)-1) / ((beta)^(-1/beta))).^(mu*real_rate/2));

        for uu = 1:length(Nv)
            BRoS = sum( ( (exp(beta-1)/(beta^beta)).^(Rv{uu}*real_rate) + (exp((1/beta)-1)/(beta^(-1/beta))).^(Rv{uu}*real_rate) ).*Prob{uu} );
            PeRoS_Bound(uu, eps) = min(1, BRoS);
        end

        PeRoS(:,eps) = mean(ERoS > 1 + eps/100, 2);
        PeMoR(:,eps) = mean(EMoR > 1 + eps/100, 2);

        alpha = 0.50;
        SampSize(qq, eps) = (log(2) + alpha * log(N)) / (real_rate * (1 - (1/beta)*(log(beta) + 1)));
    end

    PMoR_plot_data(:, 2*(qq-1) + 1) = PeMoR(:,10);
    PMoR_plot_data(:, 2*(qq-1) + 2) = PeMoR_Bound(:,10);
    data2save1(:,1) = Nv';
    data2save1(:,2) = PeMoR(:,10);
    data2save1(:,3) = PeMoR_Bound(:,10);
    eval(['save Perror_plot/PeMoR_' num2str(qq) '_RO.dat data2save1 -ascii'])


    PRoS_plot_data(:, 3*(qq-1) + 1) = PeRoS(:,10);
    PRoS_plot_data(:, 3*(qq-1) + 2) = PeRoS_Bound(:,10);
    PRoS_plot_data(:, 3*(qq-1) + 3) = PeRoS_Bound2(:,10);
    data2save2(:,1) = Nv';
    data2save2(:,2) = PeRoS(:,10);
    data2save2(:,3) = PeRoS_Bound(:,10);
    data2save2(:,4) = PeRoS_Bound2(:,10);
    eval(['save Perror_plot/PeRoS_' num2str(qq) '_RO.dat data2save2 -ascii'])

end

figure;
for qq = 1:length(genre_ind)
    subplot(3,3,qq)
    errorbar(Nv, ...
        estimation_plot_data(:,4*(qq-1) + 1), ...
        estimation_plot_data(:,4*(qq-1) + 2), ...
        'LineWidth', 2),
    hold on;
    errorbar(Nv, ...
        estimation_plot_data(:,4*(qq-1) + 3), ...
        estimation_plot_data(:,4*(qq-1) + 4), ...
        'LineWidth', 2),
    plot(Nv, ...
        ones(length(Nv), 1)*rr_vec(qq), ...
        '--', ...
        'LineWidth', 2);
    ylabel('Estimation')
    xlabel('Samples Size |S|')
    title(['Genre: ' all_genres{genre_ind(qq)} '. Real rate = ' num2str(rr_vec(qq))])
    legend('MoR', 'RoS', 'Real Rate')
    axis('tight')
    set(gca, 'XScale', 'log')
end


figure;
for qq = 1:length(genre_ind)
    subplot(3,3,qq)
    errorbar(Nv, ...
        error_plot_data(:,4*(qq-1) + 1), ...
        error_plot_data(:,4*(qq-1) + 2), ...
        'LineWidth', 2),
    hold on;
    errorbar(Nv, ...
        error_plot_data(:,4*(qq-1) + 3), ...
        error_plot_data(:,4*(qq-1) + 4), ...
        'LineWidth', 2),
%     plot(Nv, ...
%         ones(length(Nv), 1)*rr_vec(qq), ...
%         '--', ...
%         'LineWidth', 2);
    ylabel('ErrorM')
    xlabel('Samples Size |S|')
    title(['Genre: ' all_genres{genre_ind(qq)} '. Real rate = ' num2str(rr_vec(qq))])
    legend('MoR', 'RoS')
    axis('tight')
    set(gca, 'XScale', 'log')
end



figure;
for qq = 1:length(genre_ind)
    subplot(3,3,qq)
    semilogx(Nv, ...
        PMoR_plot_data(:,2*(qq-1) + 1), ...
        'LineWidth', 2),
    hold on;
    semilogx(Nv, ...
        PMoR_plot_data(:,2*(qq-1) + 2), ...
        'LineWidth', 2),
%     plot(Nv, ...
%         ones(length(Nv), 1)*rr_vec(qq), ...
%         '--', ...
%         'LineWidth', 2);
    ylabel('PeMoR > 1 + epsilon')
    xlabel('Samples Size |S|')
    title(['Genre: ' all_genres{genre_ind(qq)} '. Real rate = ' num2str(rr_vec(qq))])
    legend('MoR', 'Bound MoR')
    axis('tight')
end

figure;
for qq = 1:length(genre_ind)
    subplot(3,3,qq)
    semilogx(Nv, ...
        PRoS_plot_data(:,3*(qq-1) + 1), ...
        'LineWidth', 2),
    hold on;
    semilogx(Nv, ...
        PRoS_plot_data(:,3*(qq-1) + 2), ...
        'LineWidth', 2),
    semilogx(Nv, ...
        PRoS_plot_data(:,3*(qq-1) + 3), ...
        'LineWidth', 2),
%     plot(Nv, ...
%         ones(length(Nv), 1)*rr_vec(qq), ...
%         '--', ...
%         'LineWidth', 2);
    ylabel('PeRoS > 1 + epsilon')
    xlabel('Samples Size |S|')
    title(['Genre: ' all_genres{genre_ind(qq)} '. Real rate = ' num2str(rr_vec(qq))])
    legend('RoS', 'Bound RoS 1', 'Bound RoS 2')
    axis('tight')
end


% figure;
% semilogx(Nv, PeMoR(:,10), 'LineWidth',2), hold on;
% semilogx(Nv, PeRoS(:,10), 'LineWidth',2);
% semilogx(Nv, PeMoR_Bound(:,10),'--','LineWidth',2);
% semilogx(Nv, PeRoS_Bound(:,10),'--','LineWidth',2);
% xlabel('Sample Size |S|')
% ylabel('Perror > 1 + epsilon')
% legend('MoR', 'RoS','Bound MoR', 'Bound RoS')
% 
% 
% figure;
% errorbar(Nv, mean(RoS,2), ci_RoS(:,1), '--o','LineWidth', 2);
% hold on;
% errorbar(Nv, mean(MoR,2), ci_MoR(:,1), 'k--o','LineWidth', 2);
% hold on;
% plot(Nv, ones(length(Nv))*real_rate, 'r--','LineWidth', 2);
% ax = gca;
% ax.FontSize = 14;
% % ylim([0 0.1])
% xlabel('Number of  of surveys','FontSize',14);
% ylabel('Estimation','FontSize',14)
% title(['|G|:' num2str(N) '. Trials: ' num2str(trials)],'FontSize',14)
% legend('Ratio of sums', 'Mean of ratios', 'True value','FontSize',14)
% grid on;
% hold off;
% 
% 
% figure;
% errorbar(Nv, mean(ERoS,2), ci_ERoS(:,1), '--o','LineWidth', 2);
% hold on;
% errorbar(Nv, mean(EMoR,2), ci_EMoR(:,1), 'k--o','LineWidth', 2);
% ax = gca;
% ax.FontSize = 14;
% % ylim([0 0.1])
% xlabel('Number of  of surveys','FontSize',14);
% ylabel('Estimation','FontSize',14)
% title(['|G|:' num2str(N) '. Trials: ' num2str(trials)],'FontSize',14)
% legend('Ratio of sums', 'Mean of Ratios', 'FontSize',14)
% grid on;
% hold off;