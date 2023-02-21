%this function will store the intersection loops after a surface intersection
%to be used with boolean of stl files

%will be changed to only need to be supplied an edge list, instead of
%taking both surface structures and combining them together, will minimize
%tolerance headaches

function [edge_loops,edge_loops_directed,pinch_ids]  = form_intersection_loops_boolean(all_edges)

%Input
% all_edges : an edge list, each edge is a row
%Output
% edge_loops : cell list of connected edge loops that are unique
% edge_loops_directed : cell list groups of connected edges maintaning direction of the original edge
%
%combine the edges together 
%all_edges = unique([new_edges_tmpA; new_edges_tmpB],'rows','stable');
all_edges = unique(all_edges,'rows','stable');
%remove any edges with the same start and end id
all_edges(all_edges(:,1) == all_edges(:,2) ,:)=[];
%test all_edges
%all_edges = [1 3; 3 2; 1 4; 4 2; 1 5; 5 2; 1 6; 6 2]; %for soft loops (works)
%all_edges = [1 2; 2 3; 3 4; 4 5; 5 1; 1 9; 9 8; 8 7; 7 6; 6 1]; %for when 2 or more closed loops touch each other (works)
%all_edges = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 5; 5 9; 9 3; 3 10; 10 1]; %when there are 3 loops touching each other (closed loops)
%all_edges = [1 2; 2 3; 3 4; 4 5; 5 6; 6 4; 4 1; 7 8; 8 12; 12 9; 9 10; 10 11; 11 12; 12 13; 13 7]; %for 2 sets of closed loops toching each other
%check if any node is shared more than twice
[a,b]=hist(all_edges(:),unique(all_edges));
pinch_ids = b(a>2);
free_ends = b(a<2);
if numel(free_ends) > 1
    fprintf('Free Edges detected!\n')
%     hold on
%     scatter3(intSurface_tmp.vertices(free_ends(:),1),intSurface_tmp.vertices(free_ends(:),2),intSurface_tmp.vertices(free_ends(:),3),100,'filled')
%     scatter3(surfaceB.vertices(pinch_ids(:),1),surfaceB.vertices(pinch_ids(:),2),surfaceB.vertices(pinch_ids(:),3),100,'filled')
% 
%     stop
    %need to fill in the gaps before proceeding
end
%scatter3(nodes_master_new_FV_global(pinch_ids(:),1),nodes_master_new_FV_global(pinch_ids(:),2),nodes_master_new_FV_global(pinch_ids(:),3),'filled')
%all_edges(any(all_edges==389,2),:)

%create the list of closed loops for each of the surfaces
%connect head to tail of all edges
all_edges_tmp = all_edges;
edge_loops = cell(1);
edge_loops_directed = cell(1);
edge_node_loop_ids = unique(all_edges_tmp);
start_id = edge_node_loop_ids(1);
current_id = start_id;
loop_num = 1;
edge_collector = [];
edge_collector_directed = [];
soft_loop_on = 0;
% edge_collector_re_directed = all_edges_tmp*0;
% direct_count = 0;
%assuming no open loops exist
%counter = 0;
while size(all_edges_tmp,1)>=1
    %direct_count = direct_count + 1;
    %still need to evaulate if there is a pinch point!!!
    %fprintf('%f\n',size(all_edges_tmp,1));
    %get the edge that contains the starting node
    edge_bool_w1 = ismember(all_edges_tmp, current_id);
        
    %select the first edge (if more than 2 exist, then the node is a pinch point)
    selected_edges_tmp = all_edges_tmp(any(edge_bool_w1,2),:);
    if (ismember(current_id, pinch_ids) || ismember(current_id, free_ends)) && ~isempty(edge_collector)

        %need to create mulitple open loops, then combine them together at the end to form closed loops
        
        %reverse the loop and continue searching until you hit another pinch-point
        
        %continue until you reach another critical marker
        
        %turn on the pinch switch if off
        if soft_loop_on == 0
            soft_loop_on = 1;
            %update the starting id
            start_id=current_id;
            %invert the global collector and switch directions
            edge_collector = flipud(edge_collector);
            edge_collector_directed = rot90(edge_collector_directed,2);
            current_id = edge_collector_directed(end,2);
            %current_edge = edge_collector_directed(end,2);
            %current_edge = edge_collector(end,:);
            next_id=current_id;
            
        else
           %it is assumed that the other pinch point has been selected
           %this algortihm will not work if there are multiple pinch points
           %along a path
           soft_loop_on = 0;
           next_id = current_id;
           start_id = current_id; %end the loop (will be an open loop)
        end    

    else
        %create the edges to be appended to the edge collectors
        current_edge = selected_edges_tmp(1,:);
        %select the next edge node id
        next_id = current_edge(~ismember(current_edge,current_id));
        modified_edge = [current_id next_id];
        
        %remove the selected edge from the database
        row_query = find(any(edge_bool_w1,2));
        all_edges_tmp(row_query(1),:) = [];
        
        %store the edge into the collector
        edge_collector = [edge_collector; current_edge ]; %#ok<AGROW>
        edge_collector_directed = [edge_collector_directed; modified_edge]; %#ok<AGROW>
        
        %evaluate the numberd of edges eliminated from the selection versus
        %the number of edges supplied to the collectors, if there is a
        %diffirence, then stop the program
        
        
    end

    %asses if the next id is the same as the starting id
    if next_id == start_id %|| ismember(next_id,pinch_ids)
        %the loop is closed
        
        %add the collector the cell
        edge_loops{loop_num} = edge_collector;
        edge_loops_directed{loop_num} = edge_collector_directed;
        
        %increment the loop number
        loop_num = loop_num + 1;

        %select a new starting id
        if isempty(all_edges_tmp)
            %if all line segments or edges have been tested, then exit the
            %while loop
            break
        else
            edge_node_loop_ids = unique(all_edges_tmp);
            start_id = edge_node_loop_ids(1);
            current_id = start_id;
            %empty the collector(s)
            edge_collector = []; 
            edge_collector_directed = [];
            soft_loop_on = 0;
        end
        
    elseif isempty(all_edges_tmp)
        
        %the process ended with an open loop?
        
        edge_loops{loop_num} = edge_collector;
        edge_loops_directed{loop_num} = edge_collector_directed;
        
        break
    else
        %continue searching 
        current_id = next_id;
    end
    
    
end


end
