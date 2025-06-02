clear
close all;
filename = 'HR_genres.json';
str = fileread(filename);
data = jsondecode(str);
M = csvread("HR_edges.csv", 1);
N = numel(fieldnames(data));
load all_genres.mat

% genre_ind = [27, 80, 50, 66, 52, 71, 3, 68, 62];
genre_ind = [5, 50, 42, 44, 3, 68, 3, 68, 62];
Nv = round([0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.1 0.2 0.5] * N);
trials = 10000;

un_net = graph(M(:,1)+1, M(:,2)+1);
adj_mat= adjacency(un_net);
Ht     = sum(adj_mat, 2);

eps2t = 5;

for qq = 1:length(genre_ind)
    genre = all_genres{genre_ind(qq)};
    h = zeros(N,1);
    for i = 1:N
        eval(['B = data.x' num2str(i-1) ';']);
        h(i) = sum(strcmp(B, genre));
    end

    real_rate = sum(h)/length(h);
    rr_vec(qq) = real_rate;

    Gt = sum(adj_mat(:,h==1), 2);

    RoS = zeros(length(Nv), trials);
    MoR = zeros(length(Nv), trials);

    ERoS = zeros(length(Nv), trials);
    EMoR = zeros(length(Nv), trials);
    Rs = zeros(length(Nv), trials);

    for j = 1:length(Nv)
        [qq, j]
        for i = 1:trials
            indices = randperm(N, Nv(j));
            H = Ht(indices);
            G = Gt(indices);

            Rs(j, i) = sum(H);

            RoS(j, i) = sum(G) / sum(H);
            MoR(j, i) = mean(G./H);

            ERoS(j, i) = max([RoS(j, i)/real_rate, real_rate/RoS(j, i)]);
            EMoR(j, i) = max([MoR(j, i)/real_rate, real_rate/MoR(j, i)]);
        end
    end

    for j = 1:length(Nv)
        [ht, bins] = histcounts(Rs(j,:), linspace(min(Rs(j,:)), max(Rs(j,:)), max(Rs(j,:)) - min(Rs(j,:)) + 1));
        Prob{j} = ht/trials;
        Rv{j} = bins(1:end-1);

        [ht2, bins2] = histcounts(Rs(j,:), 100);
        Prob2{j} = ht2/trials;
        Rv2{j} = bins2(1:end-1);
    end
    
    for eps = 1:20
        beta = 1 + eps/100;

        % Real Probability
        PeRoS(:,eps) = mean(ERoS > 1 + eps/100, 2);
        PeMoR(:,eps) = mean(EMoR > 1 + eps/100, 2);

        PeMoR_Bound(:,eps) = min(1, (exp(beta-1)/(beta^beta)).^(Nv*real_rate) ...
            + (exp((1/beta) - 1)/(beta^(-1/beta))).^(Nv*real_rate));

        for uu = 1:length(Nv)
            BRoS = sum( ( (exp(beta-1)/(beta^beta)).^(Rv{uu}*real_rate) + (exp((1/beta)-1)/(beta^(-1/beta))).^(Rv{uu}*real_rate) ).*Prob{uu} );
            PeRoS_Bound(uu, eps) = min(1, BRoS);
        end

        alpha = 0.50;
        SampSize(qq, eps) = (log(2) + alpha * log(N)) / (real_rate * (1 - (1/beta)*(log(beta) + 1)));
    end

    PMoR_plot_data(:, qq) = PeMoR(:,eps2t);
    PRoS_plot_data(:, qq) = PeRoS(:,eps2t); 

    PRoS_Bound_plot_data(:, qq) = PeRoS_Bound(:,eps2t); 
    PMoR_Bound_plot_data(:, qq) = PeMoR_Bound(:,eps2t); 

    data2save1(:,1) = Nv';
    data2save1(:,2) = PeMoR(:,eps2t);
    data2save1(:,3) = PeMoR_Bound(:,eps2t);
    eval(['save Perror/PeMoR_' num2str(qq) '_eps_' num2str(eps2t) '_HR.dat data2save1 -ascii'])

    data2save2(:,1) = Nv';
    data2save2(:,2) = PeRoS(:,10);
    data2save2(:,3) = PeRoS_Bound(:,10);
    eval(['save Perror/PeRoS_' num2str(qq) '_eps_' num2str(eps2t) '_HR.dat data2save2 -ascii'])

end

figure;
for j = 1:length(Nv)
    subplot(2,5,j)
    bar(Rv2{j}, Prob2{j})
    ylabel('Probability', 'FontSize',16,'Interpreter','latex')
    xlabel('$R_s$', 'FontSize',16,'Interpreter','latex')
    title(['$|S| = $' num2str(Nv(j))], 'FontSize',16,'Interpreter','latex')
end
f = gcf ;
f.WindowState = 'maximized';
print(gcf, '-dpng', 'HR_hist_Rs.png','-r600');
exportgraphics(gcf,"myplot.png")

figure;
for j = 1 : length(genre_ind)
    subplot(3,3,j);
    semilogx(Nv, PMoR_plot_data(:,j),'LineWidth',2), hold on;
    semilogx(Nv, PMoR_Bound_plot_data(:,j),'LineWidth',2);
    plot([SampSize(j, eps2t) SampSize(j, eps2t)], [0 1], '--', 'LineWidth',2), hold off;
    xlim([min(Nv), max(Nv)])
    xlabel('Sample Size |S|')
    ylabel('PeMoR > 1 + 0.1');
    legend('Simulation', 'Theoretical')
    title(['Genre: ' all_genres{genre_ind(j)} '. Real rate = ' num2str(rr_vec(j))])
end

figure;
for j = 1 : length(genre_ind)
    subplot(3,3,j);
    semilogx(Nv, PRoS_plot_data(:,j),'LineWidth',2), hold on;
    semilogx(Nv, PRoS_Bound_plot_data(:,j),'LineWidth',2);
    plot([SampSize(j, eps2t) SampSize(j, eps2t)], [0 1], '--', 'LineWidth',2), hold off;
    xlim([min(Nv), max(Nv)])
    xlabel('Sample Size |S|')
    ylabel('PeMoR > 1 + 0.1');
    legend('Simulation', 'Theoretical')
    title(['Genre: ' all_genres{genre_ind(j)} '. Real rate = ' num2str(rr_vec(j))])
end