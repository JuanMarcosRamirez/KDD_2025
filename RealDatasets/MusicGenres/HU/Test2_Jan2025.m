clear
close all;
filename = 'HU_genres.json';
str = fileread(filename);
data = jsondecode(str);
M = csvread("HU_edges.csv", 1);
N = numel(fieldnames(data));
load all_genres.mat

% genre_ind = [50, 42, 44, 42, 44, 68];
Nv = round([0.001 0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.1] * N);
genre_ind = 74;
trials = 10000;

un_net = graph(M(:,1)+1, M(:,2)+1);
adj_mat= adjacency(un_net);
Ht     = sum(adj_mat, 2);
genre = all_genres{genre_ind};

h = zeros(N,1);
for i = 1:N
    eval(['B = data.x' num2str(i-1) ';']);
    h(i) = sum(strcmp(B, genre));
end
real_rate = sum(h)/length(h);
Gt = sum(adj_mat(:,h==1), 2);

t = -0.20:0.005:0.20;

EX1 = zeros(length(t), length(Nv));
EX2 = zeros(length(t), length(Nv));
for k = 1:length(Nv)
    Z = zeros(Nv(k), trials);
    for i = 1:trials
        indices = randperm(N, Nv(k));
        H = Ht(indices);
        G = Gt(indices);
        Z(:,i)= G./H;
    end
    
    Z1 = sum(Z);
    Z2 = Z1;
    for i = 1:length(t)
        disp(['Size: ', num2str(Nv(k)),'. t = ', num2str(t(i)), '.'])
        X1 = exp(t(i) * Z2);
        EX1(i,k) = mean(X1);
        Z_hat = binornd(Nv(k),real_rate,trials,1);
        X2 = exp(t(i) * Z_hat);
        EX2(i,k) = mean(X2);
    end
end

alpha = 1.0;
% contourf(log10(Nv),t,log(alpha*EX2-EX1))
contourf(Nv,t,log(alpha*EX2-EX1))
set(gca,'xscale','log')
xlabel('|S|')
ylabel('t')
colorbar
title(['MoR ($\rho = $' num2str(100*real_rate,4) '\%): $\log(E[e^{t\hat{Z}}] - E[e^{tZ}])$'],'Interpreter','latex')

% semilogy(t,EX1), hold on;
% semilogy(t,EX2);
% xlabel('t')
% legend('From Real Network','Binomial Distribution')
