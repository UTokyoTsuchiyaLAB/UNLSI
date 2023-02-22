%this is a sub function for the boolean operations

%this organizes all the sub-surfaces together given an edge loop list and a
%surface with the edge list imprinted onto it
%https://arxiv.org/pdf/1308.4434.pdf

function [newSF,num_loops_matrix,orientation_loops_matrix] = create_sub_surfaces_boolean(FV, edge_loops,visual)

        %input, set of closed loops and a triangular surface S
        %output, a set of sub-surfaces
        tria_surf = FV.faces;
        %divide up the surfaces into sub surfaces
        
        %need to stop growing or selecting triangles when it reaches other edge loops
        
        %create triangle to edge connectiivty list
        all_edges_list = unique([tria_surf(:,1:2); tria_surf(:,2:3); tria_surf(:,[3 1])],'rows');
        [~,ids_A1] = ismember(tria_surf(:,1:2),all_edges_list,'rows');
        [~,ids_A2] = ismember(tria_surf(:,2:3),all_edges_list,'rows');
        [~,ids_A3] = ismember(tria_surf(:,[3 1]),all_edges_list,'rows');
        tria_edge_connectivity_list = [ids_A1 ids_A2 ids_A3];
        
        %create another connectivity list that does not take into account the edge directions
        all_edges_list_unq = unique(sort([tria_surf(:,1:2); tria_surf(:,2:3); tria_surf(:,[3 1])],2),'rows');
        [~,ids_A1_unq] = ismember(sort(tria_surf(:,1:2),2),all_edges_list_unq,'rows');
        [~,ids_A2_unq] = ismember(sort(tria_surf(:,2:3),2),all_edges_list_unq,'rows');
        [~,ids_A3_unq] = ismember(sort(tria_surf(:,[3 1]),2),all_edges_list_unq,'rows');
        tria_edge_connectivity_list_unq = [ids_A1_unq ids_A2_unq ids_A3_unq];
        
        trias_FV_tmp = tria_surf;
        newSF =cell(length(edge_loops)+1,1);
        loop_selection = 1:length(edge_loops);
        %create a matrix to distinguish if the contour was flipped or not
        flipped = 0;
        for hj = 1:length(edge_loops)
            
            selected_loop = edge_loops{hj};
            %while loop if no trias selected but trias still remain,
            %reverse the loop direction
            max_iter = 0; %to avoid infinite while loops
            while max_iter<3
                
                max_iter = max_iter +1;
                
                %create a list of triangles to omit from the selection as they
                %are part of other loops
                
                edge_omits = edge_loops(loop_selection(~ismember(loop_selection,hj)));
                edge_omits = unique(sort(cell2mat(edge_omits(:)),2),'rows'); %flip the direction of the edge? 
               
                %how to ensure that the loop is counter clockwise though? make sure that if a triangle has a another edge loop edge,
                
                %that no other tirangle can share those edges afterwards... (do not need directions in this case...
                
                %store edges that can only be shared once by any triangle...
                if ~isempty(edge_omits)
                    [~,omit_loop_edge_ids] = ismember(sort(edge_omits,2), sort(all_edges_list_unq,2),'rows');
                end
                
                %create loop connectivity list
                [~,selected_loop_edge_ids] = ismember(selected_loop, all_edges_list,'rows');
                [~,selected_loop_edge_ids_unq] = ismember(sort(selected_loop,2), all_edges_list_unq,'rows'); %there are zeros..?
    
                %grab triangles with edge directions matching the loop direction
                
                ta_bool_tot = ismember(tria_edge_connectivity_list,selected_loop_edge_ids); %needs to be updated
                grown_trias = trias_FV_tmp(any(ta_bool_tot,2),:); %tria to nodes
                num_trias = size(grown_trias,1);
    
                %get non-unique edges of the grown trias (direction does not matter anymore)
                grown_trias_edges = tria_edge_connectivity_list_unq(any(ta_bool_tot,2),:); %tria to edges
                
                %get the directed trias to edge list as well
                grown_trias_edges_directed =  tria_edge_connectivity_list(any(ta_bool_tot,2),:); %tria to directed edges
    
                %select the next triangle edges (a list)
                next_edges_to_find_unq = unique(grown_trias_edges(~ismember(grown_trias_edges,selected_loop_edge_ids_unq)));
                
                %query out edges that are attached to other contour loops
                if ~isempty(edge_omits)
                    edges_to_other_loops_bool = ismember(next_edges_to_find_unq, omit_loop_edge_ids);
                    next_edges_to_find_unq(any(edges_to_other_loops_bool,2),:) = [];
                end
                
                %select triangles which share edges (not with the closed loop)
                next_bool = ismember(tria_edge_connectivity_list_unq,next_edges_to_find_unq);
                
                %must ensure that all the edges are selected before passing to the next
                %step, assuming that the correct orientations are set before
                
                
                if any(any(next_bool))
                    break;
                elseif ~any(any(next_bool)) && ~isempty(trias_FV_tmp)
                    selected_loop = fliplr(edge_loops{hj});
                    flipped = 1;
                else
                    break;
                end
            %if no triangles are found for a given loop, but triangles exist, then this while loop will go on for ever...
            end %end while loop
            diff_trias = 1;
            while diff_trias ~= 0 %while the number of triangles is increasing

                %find the new triangles (with nodal ids)
                triangles_added = trias_FV_tmp(any(next_bool,2),:);
                triangles_added = triangles_added(~ismember(sort(triangles_added,2),sort(grown_trias,2),'rows'),:); %remove previously selected triangles (nodal ids)
                %get the unique edge ids of the new triangles
                edges_ids_bulk = tria_edge_connectivity_list_unq(any(next_bool,2),:);
                edges_ids_bulk = edges_ids_bulk(~ismember(sort(edges_ids_bulk,2),sort(grown_trias_edges,2),'rows'),:); %remove previously selected triangles (edges ids)
                %get the directed edge ids of the new triangles
                directed_ids_bulk = tria_edge_connectivity_list(any(next_bool,2),:);
                directed_ids_bulk = directed_ids_bulk(~ismember(sort(directed_ids_bulk,2),sort(grown_trias_edges_directed,2),'rows'),:); %remove previously selected triangles (edges ids)
                
                %add the new triangles to the surface
                grown_trias = unique([grown_trias; triangles_added],'rows','stable'); %new trias list (triangles in terms of nodal ids)
                grown_trias_edges = unique(sort([grown_trias_edges; edges_ids_bulk],2),'rows','stable'); %new tris to edge id list
                
                %evaluate the size
                diff_trias = num_trias - size(grown_trias,1);
                num_trias = size(grown_trias,1);

                %remove the triangles from the database (tris with node ids)
                rem_tri_bool1 = ismember(sort(trias_FV_tmp,2),sort(grown_trias,2),'rows');
                trias_FV_tmp(rem_tri_bool1,:) = [];
                
                %update the edge connectivity list, remove unneded trias
                rem_tri_bool2 = ismember(sort(tria_edge_connectivity_list_unq,2),sort(grown_trias_edges,2),'rows');
                tria_edge_connectivity_list_unq(rem_tri_bool2,:) = [];
                
                %add triangles with directed edges into another collector
                grown_trias_edges_directed = [grown_trias_edges_directed; directed_ids_bulk];
                
                %remove triangles from the database (tris with edge ids {with directions})
                tria_edge_connectivity_list(rem_tri_bool1,:) =[];
                
                %select the next criteria for finding more triangles
                %get the edge list of the added triangles
                
                current_edge_ids = unique(grown_trias_edges);
                
                %get the edges that are refrenced once (they are the free edges)
                [~, ~, edgeNos] = unique(grown_trias_edges);
                hm = accumarray(edgeNos,1);
                edge_list_mat = [current_edge_ids hm];
                edge_search = edge_list_mat(hm==1,1);
                
                %query out the edges that include the closed loop contour
                next_edges_to_find_unq = edge_search(~ismember(edge_search, selected_loop_edge_ids_unq));
                
                %query out edges that are attached to other contour loops
                if ~isempty(edge_omits)
                    edges_to_other_loops_bool = ismember(next_edges_to_find_unq, omit_loop_edge_ids);
                    next_edges_to_find_unq(any(edges_to_other_loops_bool,2),:) = [];
                end
                %create the boolean to select the next triangles
                next_bool = ismember(tria_edge_connectivity_list_unq,next_edges_to_find_unq);          
             
                if visual == 1 && any(any(next_bool))
                    %debug_plots
                    figure(9)
                    clf;
                    view(3);
                    grid on
                    hold on
                    %grown_trias = trias_FV_tmpA(ta_bool_tot,:);
                    %TRA = triangulation(grown_trias ,nodes_master_new_FV_global);    
                    %trisurf(TRA)
                    FV_test.faces = grown_trias;
                    FV_test.vertices = FV.vertices;
                    %patch(FV_test,'FaceColor','b','FaceAlpha',0.25,'Displayname','Design Space')
                    TGH = triangulation(FV_test.faces,FV_test.vertices);
                    trisurf(TGH,'edgecolor','none')
                    %trisurf(TGH,'edgecolor','black')
                    axis equal
                    %plot the sharp features aswell
                    F = featureEdges(TGH,pi/10);
                    %hold off
                    if 1==2
                        %plot the normal vectors
                        F2 = faceNormal(TGH);
                        centroids = (TGH.Points(TGH.ConnectivityList(:,1),:)+TGH.Points(TGH.ConnectivityList(:,2),:)+TGH.Points(TGH.ConnectivityList(:,3),:))/3;
                        quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'Color','b'); %check the face normals, must be consistent

                    end
                end
                if ~any(any(next_bool))
                    %no more trias :(
                    break
                end
                
            end
            %add the triangles list into a cell structure
            newSF{hj}.trias = grown_trias;
            newSF{hj}.edges = grown_trias_edges;
            newSF{hj}.directed_edges = grown_trias_edges_directed;
            newSF{hj}.orientation = flipped; %the opposing contours must be opposite (reverse the loop orientation)
        end


%add the unsearched traingles into the structure at the end
newSF{end}.trias = trias_FV_tmp;
newSF{end}.edges = tria_edge_connectivity_list_unq;
newSF{end}.directed_edges = tria_edge_connectivity_list;
newSF{end}.orientation = 1;

%determine the number of loops each surface shares and form a decision matrix
num_loops_matrix = zeros(length(newSF),length(edge_loops));
orientation_loops_matrix = zeros(length(newSF),length(edge_loops));
for uuj = 1:length(newSF)
    
    current_surface_edges = newSF{uuj}.edges;
    current_directed_edges = newSF{uuj}.directed_edges;
    %check which edges are present in the surface
    for uit = 1:length(edge_loops)
         %regular connectivity matrix creation
         selected_loop = edge_loops{uit};
         %evaluate edges that are not part of the global edge list in the
         %mesh and omit them (this is a patch due to an error in previous sections, not sure how to fix it)
         [~,selected_loop_edge_ids_unq] = ismember(sort(selected_loop,2), all_edges_list_unq,'rows');
         all_boolean = ismember(selected_loop_edge_ids_unq(selected_loop_edge_ids_unq~=0),current_surface_edges); %must update this line when debugging for finding the root cause of the problem
         num_loops_matrix(uuj,uit) = (all(all_boolean));
         %contour orientation along the surface connectivity matrix
         [~,selected_loop_edge_ids] = ismember(selected_loop, all_edges_list,'rows');
         contained_boolean = ismember(selected_loop_edge_ids,current_directed_edges);

         if all(all(contained_boolean))==1 && num_loops_matrix(uuj,uit)==1
             %clockwise (or couter-clockwise??)
             orientation_loops_matrix(uuj,uit) = 1;
         elseif num_loops_matrix(uuj,uit)==1
             %reversed
             orientation_loops_matrix(uuj,uit) = -1;
         end
         
    end
    
    
end

%remove empty feilds
idx_rm = zeros(length(newSF),1);
for w = 1:length(newSF)

    empt_bool = isempty(newSF{w}.trias);
    if empt_bool
        idx_rm(w) = 1;
    end
    
end

%save the fields correctly
newSF = newSF(find(~idx_rm));
num_loops_matrix = num_loops_matrix(find(~idx_rm),:);
orientation_loops_matrix = orientation_loops_matrix(find(~idx_rm),:);

end