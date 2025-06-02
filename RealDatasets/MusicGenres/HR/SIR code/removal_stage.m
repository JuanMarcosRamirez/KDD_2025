function [inf,rec] = removal_stage(gamma,infected_nodes,recovered_nodes)
% gamma: removal rate
% infected_nodes: binary vector, the ones represent the infected nodes
% recovered_nodes: binary vector, the ones represent the recovered nodes

new_infected_nodes = infected_nodes;

for i=1:length(infected_nodes)
    if infected_nodes(i)==1
        if unifrnd(0,1) < gamma
            recovered_nodes(i) = 1;
            new_infected_nodes(i)=0;
        end
    end
end
inf = new_infected_nodes;
rec = recovered_nodes;
end

            
