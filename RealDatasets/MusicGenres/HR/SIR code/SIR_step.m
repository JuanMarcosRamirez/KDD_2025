function [inf,rec,vis] = SIR_step(graph,beta,gamma,infected_nodes,recovered_nodes)

inf = infection_stage(graph,beta,infected_nodes,recovered_nodes);

[inf,rec] = removal_stage(gamma,inf,recovered_nodes);

vis = infection_visibility(graph,inf);
end
