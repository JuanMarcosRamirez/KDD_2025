clear
close all;
filename = 'HR_genres.json';
str = fileread(filename);
data = jsondecode(str);
M = csvread("HR_edges.csv", 1);
N = numel(fieldnames(data));
load all_genres.mat

genre_ind = [27, 80, 50, 66, 52, 71, 3, 68, 62];
% Nv = round([0.002 0.0035 0.005 0.01 0.02 0.035 0.05 0.1 0.2 0.5] * N);
% Nv = round(0.02 * N);
Nv = 1000;
genre_ind = 50;
trials = 100000;

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

Z = zeros(Nv, trials);
% prod_set = zeros(2^Nv-1,trials);
for i = 1:trials
    indices = randperm(N, Nv);
    H = Ht(indices);
    G = Gt(indices);
    Z(:,i)= G./H;

    % vec = Z(:,i)';
    % products = []; 
    % n = length(vec);
    % for k = 1:n
    %     combinations_k = nchoosek(vec, k);
    %     products_k = prod(combinations_k, 2); % Row-wise product
    %     products = [products; products_k];
    % end
    % prod_set(:,i) = products; 
end

% EP = mean(prod_set, 2);
% 
evec = mean(Z, 2);
% PE = [];
% n = length(evec);
% for k = 1:n
%     combinations_k = nchoosek(evec, k);
%     products_k = prod(combinations_k, 2); % Row-wise product
%     PE = [PE; products_k];
% end
% mean(EP < PE)

% % evec = mean(Z,2);
% 
% evec = zeros(Nv,1);
% num_bins = 100; % Number of bins
% for i = 1:Nv
%     [counts, edges] = histcounts(Z(i,:), num_bins); % Counts and bin edges
%     bin_centers = edges(1:end-1) + diff(edges)/2;
%     probabilities = counts / sum(counts);
%     evec(i) = sum(bin_centers .* probabilities);
% end
% 
% PE = [];
% n = length(evec);
% for k = 1:n
%     combinations_k = nchoosek(evec, k);
%     products_k = prod(combinations_k, 2); % Row-wise product
%     PE = [PE; products_k];
% end
% 
% EP = zeros(2^Nv-1,1);
% for i = 1:2^Nv-1
%     [counts, edges] = histcounts(prod_set(i,:), num_bins); % Counts and bin edges
%     bin_centers = edges(1:end-1) + diff(edges)/2;
%     probabilities = counts / sum(counts);
%     EP(i) = sum(bin_centers .* probabilities);
% end
% % EP = mean(prod_set,2);
% mean(EP < PE)