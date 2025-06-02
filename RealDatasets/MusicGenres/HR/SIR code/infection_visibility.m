function vis = infection_visibility(graph,infected_nodes)

vis = zeros(length(infected_nodes),1);
for i=1:length(infected_nodes)
    nei = graph.neighbors(i);
    for j=1:length(nei)
        if infected_nodes(nei(j))==1
            vis(i)=vis(i)+1;
        end
    end
end
end
