function [infected_nodes,recovered_nodes,vis] = SIR_net(graph,beta,gamma,infected_nodes,recovered_nodes,min_inf)

% inf_matrix = zeros(length(infected_nodes),nsteps);
% rec_matrix = zeros(length(infected_nodes),nsteps);
% vis_matrix = zeros(length(infected_nodes),nsteps);
%for iteration = 1:nsteps
prm = 0;
while prm < min_inf 
    [infected_nodes,recovered_nodes,vis] = SIR_step(graph,beta,gamma,infected_nodes,recovered_nodes);
%     inf_matrix(:,iteration) = infected_nodes;
%     rec_matrix(:,iteration) = recovered_nodes;
%     vis_matrix(:,iteration) = vis;
    prm = sum(infected_nodes)
end

end
    