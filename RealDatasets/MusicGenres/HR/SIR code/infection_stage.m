function h = infection_stage(graph,beta,infected_nodes,recovered_nodes)
% graph, the underlying network
% beta, infection probability
% infected nodes, binary vector, the ones represent infected individuals

new_infected_nodes = infected_nodes;
for i = 1:length(infected_nodes)
    if infected_nodes(i) == 1
%         nei = graph.neighbors(infected_nodes(i));
        nei = graph.neighbors(i);
        for j =  1:length(nei)
            if unifrnd(0,1) < beta && new_infected_nodes(nei(j))==0 && recovered_nodes(nei(j))==0
                new_infected_nodes(nei(j)) = 1;
            end
        end
    end
end

h = new_infected_nodes;
end