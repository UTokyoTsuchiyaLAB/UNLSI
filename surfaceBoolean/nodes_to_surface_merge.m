%this function will combine the surface nodes from the truncation algotihm
%with the surface(iso) mesh of the design surface
function [trias,nodes_master_new_tmp,edges_tmp,nodes_at_edges] = nodes_to_surface_merge(nodes_master_new,surface_nodes,FV , edges,tolerance)
    %it must be remeshed in hypermesh before being sent to the fem many bad
    %trias will be created
        
    FV_tmp = FV;
    
    TRf = triangulation(FV_tmp.faces,FV_tmp.vertices);    
    normals_all = faceNormal(TRf);

    %set up required information
    num_trias = size(FV_tmp.faces,1);
    trias_cell = cell(num_trias,1);
    %update the global nodes list
    %tolerance = 10^-8;
    %surface nodes mmust already be a part of the master list
    nodes_master_new_tmp = [nan(3); nodes_master_new; FV_tmp.vertices]; %to not re-order, or will break other functions
    
    %renumber the nodes in the iso surface
    iso_nodes = FV_tmp.vertices(FV_tmp.faces(:),:);
    [~,iso_nodes2] = ismembertol(iso_nodes,nodes_master_new_tmp,tolerance,'Byrows',true);
    FV_tmp.faces(:) = iso_nodes2;
    %FV_tmp.faces = FV_tmp.faces + size(nodes_master_new,1)+3
   
    %evaluate if the edges are supplied
    constraint_simple = [1 2; 2 3; 3 1]; %simple constraint
    if isempty(edges)
        edges_on = 0;
        edges_tmp = [];
    else
        edges_tmp = edges;
        edges_on = 1;

        %renumber the edge ids
        [~,edges_idnew] = ismembertol(nodes_master_new(edges_tmp(:),:),nodes_master_new_tmp,tolerance,'Byrows',true);
        edges_tmp(:) = edges_idnew; 
                
    end

    %get the surfacenode ids
    %surface_nodes = uniquetol(surface_nodes,tolerance,'Byrows',true);
    [~,surface_nodes_ids] = ismembertol(surface_nodes,nodes_master_new_tmp,tolerance,'Byrows',true);
    num_surf_nodes = size(surface_nodes,1);
    
    %deal with duplicate surface nodes
    surf_nodes_unq_ids = unique(surface_nodes_ids);
    if numel(surf_nodes_unq_ids)~=numel(unique(surface_nodes_ids))
        %need to remove some nodes so that each surface node has a unique
        %id identifier
        stop %error
    end
    
    %initilize empty cells if additional nodes are required (when edges overlap)
    new_nodes_cell = cell(num_trias,1);
    new_edges_cell = cell(num_trias,1);
    new_tris_cell = cell(num_trias,1);
    nodal_ids_at_edges = cell(num_trias,1);
    %evaluate the largest area among the triangles
    %[A]= areaIsosurface(FV.faces,FV.vertices);
    %max_area = max(A);
    
    %calculate the normals before entering the loop
    %normals_all = Tnorm(nodes_master_new_tmp,FV_tmp.faces);
    
    %debug_nodes_cell = cell(num_trias,1);
    
    tic
    
    
    parfor i = 1:num_trias
    
        %evaluate if there are any of the surface nodes inside it
        triangle = FV_tmp.faces(i,:);
        tri_nodes = nodes_master_new_tmp(triangle(:),:);
        normal = normals_all(i,:);
        centroid_t = mean(tri_nodes);
        
        %directions from centroid to all other surface nodes
        directions_to_surf_nodes = centroid_t - surface_nodes;
        %normalize
        directions_to_surf_nodes =  directions_to_surf_nodes./sqrt(sum(directions_to_surf_nodes.^2,2));
        %dot product
        ortho = dot(repmat(normal,num_surf_nodes,1),directions_to_surf_nodes,2);
        %boolean_orthogonal = ortho > -0.001 & ortho <0.001;
        boolean_orthogonal = ismembertol(ortho,0,tolerance);
        %selected_nodes = surface_nodes(boolean_orthogonal,:);
        selected_nodes_ids = surface_nodes_ids(boolean_orthogonal,:);
        
        %%% Causes some edges to no longer be a valid selection in proceeding boolean operations...?
        %filter out nodal ids that are selected twice
        selected_nodes_ids = unique(selected_nodes_ids,'stable');
        %get node ids that are the same as the current triangle nodal ids
        is_in_tri_bool = ismember(selected_nodes_ids,triangle);
        
        
        %store the nodes that are causing problemsinto a seperate structure for debugging purposes
        %if ~isempty(selected_nodes_ids)
        %    debug_nodes_cell{i} = nodes_master_new_tmp(selected_nodes_ids(is_in_tri_bool),:);
        %end
        
        %remove surface nodes which are equal to nodes in the triangle
        selected_nodes_ids(is_in_tri_bool) = [];
        %reselect the nodes again
        selected_nodes = nodes_master_new_tmp(selected_nodes_ids,:);
        
        %%%
        %%scatter3(selected_nodes(:,1),selected_nodes(:,2),selected_nodes(:,3))
        if any(any(boolean_orthogonal))

            M = rotate_vect_by_another(normal,[0 0 1]); 
            %ROTATE THE NODES%
            %the selected surface nodes
            tem_surf = selected_nodes - centroid_t;
            tem_surf = tem_surf*M; %test diffirent orientations
            tem_surf = tem_surf + centroid_t;
            %the FV nodes
            tri_nodestm = tri_nodes - centroid_t;
            tri_nodestm = tri_nodestm*M;
            tri_nodestm = tri_nodestm + centroid_t;
            
            if edges_on
                
                %turn on edge constraints (assuming that the edges lie
                %inside the triangles with intersections at the triangle
                %edges
                
                %evaluate if any of the nodes are edge nodes (nodes on the same plane as the triangle)
                edge_boolean1 = ismember(edges_tmp, selected_nodes_ids);
                bad_selection = 0; %trias with these edges, will be omitted in the final surface
                if (any(any(edge_boolean1)))
                    %an edge boolean is within the plane
                    %find out if any of the edge nodes are touching the
                    %boundary of the triangle.
                    
                    %if there are, then collect them and link them to each
                    %other edge by distance from a specified node
                    
                    nodeA = tri_nodestm([1 2 3],1:2);
                    nodeB = tri_nodestm([2 3 1],1:2);
                    
                    %original distances for edges of triangle
                    edge_distance = sqrt(sum((nodeA-nodeB).^2,2));
                    
                    %x_distances
                    dxA = nodeA(:,1) - tem_surf(:,1)';
                    dxB = nodeB(:,1) - tem_surf(:,1)';
                    %y_distance
                    dyA = nodeA(:,2) - tem_surf(:,2)';
                    dyB = nodeB(:,2) - tem_surf(:,2)';
                    
                    %compile total distances among all nodes
                    dist_node_A = sqrt(dxA.^2 + dyA.^2);
                    dist_node_B = sqrt(dxB.^2 + dyB.^2);
                    
                    %evaluate the boolean
                    dist_bool1 = ismembertol(dist_node_A + dist_node_B,edge_distance,tolerance);
                    
                    if any(any(dist_bool1))
                        %create the constraints
                        %3X N where N is the number of points
                        
                        %distances from node A
                        ratios = dist_node_A./edge_distance;
                        constraint = [];
                        edge_constraint = cell(3,1);
                        for iq = 1:3
                            
                            %the original edge constraint
                            edge_original = constraint_simple(iq,:); %node A is the first
                            
                            %organize based on distance
                            dist_curr = ratios(iq,dist_bool1(iq,:));
                            surf_ids_curr = selected_nodes_ids(dist_bool1(iq,:),:);
                            
                            %create the bad selection edges
                            if ~isempty(surf_ids_curr)
                                bad_selection = 1;
                            end
                            
                            %combine all nodes by ids together then sort
                            org_mat = uniquetol(sortrows([dist_curr' surf_ids_curr; 0 edge_original(1); 1 edge_original(2)],1),tolerance,'Byrows',true);
                            num_o = size(org_mat,1);
                            cont_tmp = zeros(num_o-1,2);
                            for hj = 1:(num_o-1)
                               sel_rows = hj:(hj+1);
                               cont_tmp(hj,:) = org_mat(sel_rows,2)'; 
                            end
                            constraint = [constraint; cont_tmp];
                            edge_constraint{iq} = cont_tmp;
                        end
                    
                       %debug plot for nodes and creating constraints
%                       figure(3)
%                       clf
%                       hold on
%                       scatter(tem_surf(:,1),tem_surf(:,2))
%                       text(tem_surf(:,1),tem_surf(:,2),cellstr(num2str(selected_nodes_ids(:))))
%                       axis equal
%                       scatter(tri_nodestm(:,1),tri_nodestm(:,2))
%                       text(tri_nodestm(:,1),tri_nodestm(:,2),{'1','2','3'})
                    
                        %generate a polyshape from the contour, need to
                        %link up heads to tails of the constraint array
                        %this code is not needed but is calculated
                        %optionally (for polyshape)
                        
                        constraint_tmp = constraint;
                        id_start = constraint_tmp(1);
                        contour_constraint = constraint;
                        for ert = 1:size(constraint,1)
                            dfb = ismember(constraint_tmp,id_start);
                            if any(any(dfb))
                                tip_bool = find(any(dfb,2));
                                sel_c = constraint_tmp(tip_bool(1),:);
                                contour_constraint(ert,:) = sel_c;
                                id_start = sel_c(~ismember(sel_c,id_start));
                                constraint_tmp(tip_bool(1),:) = [];
                            else
                                break
                            end
                        end
                        %create the polyshape
                        cont_path = unique(contour_constraint,'stable');
                         nodal_ids_at_edges{i} = cont_path(:); %for closing gaps in the mesh after remeshing (if errors occur)                        
                        %evaluate the rest of the selected nodes and use
                        %the ones which are completely inside the triangles
                        %only!!
                        tem_bool = zeros(size(tem_surf,1),1);
                        P1 = tri_nodestm(1,1:2); P2 = tri_nodestm(2,1:2); P3 = tri_nodestm(3,1:2);
                        P12 = P1-P2; P23 = P2-P3; P31 = P3-P1;
                        for kl = 1:size(tem_surf(:,1:2),1)
                               Pn=tem_surf(kl,1:2);
                               Pn_id = selected_nodes_ids(kl,:);
                               if ismember(Pn_id,constraint)
                                   t = 1;
                               else
                                   t = sign(det([P31;P23]))*sign(det([P3-Pn;P23])) >= 0 & sign(det([P12;P31]))*sign(det([P1-Pn;P31])) >= 0 & sign(det([P23;P12]))*sign(det([P2-Pn;P12])) >= 0 ;
                               end
                               tem_bool(kl) = t;     
                        end
                        tem_bool = logical(tem_bool);
                        %remove bad interior edges
                        interior_edges = edges_tmp(all(edge_boolean1,2),:);
                        if any(~tem_bool)
                            %remove contraints.. with the selected node
                            rem_node_id = selected_nodes_ids(~tem_bool,:);
                            rem_bool = ismember(interior_edges,rem_node_id);
                            interior_edges(any(rem_bool,2),:)=[];
                            %remove bad nodes outside of the traingle (but
                            %not those on the edges)
                            tem_surf(~tem_bool,:)=[];
                            selected_nodes_ids(~tem_bool)=[];
                        end
                        
                        %interpolate the ids so that the highest id is
                        %equal to the number of surface nodes found + 3
                        
                        constraint = [constraint; interior_edges];
                        %interp
                        interp_surf_nodes = [1:3 selected_nodes_ids']; %some inccorrect ids are being selected..?
                        lib = [1:3 4:(length(selected_nodes_ids)+3)];
                        constraint = unique(sort(interp1(interp_surf_nodes,lib, constraint),2),'rows');
                        %update edge constraint cells
                        for kkpp = 1:3
                            edge_constraint{kkpp} = interp1(interp_surf_nodes,lib,unique(edge_constraint{kkpp}));
                        end
                        cont_path = interp1(interp_surf_nodes,lib,cont_path);
                        P = [tri_nodestm(:,1:2); tem_surf(:,1:2)];
                       
                        
                        %pgon = polyshape(P(interp1(interp_surf_nodes,lib, cont_path),:));
                                                
                        %debugging plots
                        %figure(2)
                        %clf;
                        %hold on
                        %axis equal
                        %g2 = LATTICEfcn.plot_many_line_segments([P zeros(size(P,1),1)], constraint ) ;
                        %set(g2, 'Color','blue','LineWidth',0.8,'Displayname','Lattice Net')
                        %scatter(P(:,1),P(:,2))
                        
                    else
                        constraint = constraint_simple;
                        %remove nodes that are not in the triangle (same as before)
                        
                        %Step 1 Bonding Box
                        x_side_bool = tem_surf(:,1) < min(tri_nodestm(:,1)) | tem_surf(:,1) > max(tri_nodestm(:,1));
                        y_side_bool = tem_surf(:,2) < min(tri_nodestm(:,2)) | tem_surf(:,2) > max(tri_nodestm(:,2));
                        tem_surf(x_side_bool | y_side_bool,:)=[];
                        if ~isempty(tem_surf)
                            %Step 2 Robust calculation of in-out of triangle
                            tem_bool = zeros(size(tem_surf,1),1);
                            P1 = tri_nodestm(1,1:2); P2 = tri_nodestm(2,1:2); P3 = tri_nodestm(3,1:2);
                            P12 = P1-P2; P23 = P2-P3; P31 = P3-P1;
                            for kl = 1:size(tem_surf(:,1:2),1)
                                   Pn=tem_surf(kl,1:2);
                                   Pn_id = selected_nodes_ids(kl,:);
                                   if ismember(Pn_id,constraint)
                                       t = 1;
                                   else
                                       t = sign(det([P31;P23]))*sign(det([P3-Pn;P23])) >= 0 & sign(det([P12;P31]))*sign(det([P1-Pn;P31])) >= 0 & sign(det([P23;P12]))*sign(det([P2-Pn;P12])) >= 0 ;
                                   end
                                   tem_bool(kl) = t;     
                            end
                            tem_bool = logical(tem_bool);
                            tem_surf(~tem_bool,:)=[];
                        end
                        %create the points list and the lib array
                        P = [tri_nodestm(:,1:2); tem_surf(:,1:2)]; 
                        lib = 1:(3+size(selected_nodes,1));
                        
                    end
                else
                    constraint = constraint_simple;
                    P = [tri_nodestm(:,1:2); tem_surf(:,1:2)]; 
                    lib = 1:(3+size(selected_nodes,1));
                end
                
            else
                constraint = constraint_simple;
                P = [tri_nodestm(:,1:2); tem_surf(:,1:2)];
                lib = 1:(3+size(selected_nodes,1));
            end
            
            %coincident nodes error occur with edge constraints turn on...
            
            DTP = delaunayTriangulation(P, constraint); %round to 6 to avoid slim triangles
                        
            if edges_on   
                %so that all the triangle normals face outwards, flip them
                new_trias = fliplr(DTP.ConnectivityList);
                %check for bad trias with "bad selection"
                if bad_selection %perform triangle validity checks

                    %need to also check if a triangle has all 3 nodes along
                    % an edge (Matlab Glitch)
                    glitch_bool_1 = ismember(new_trias,unique(edge_constraint{1}));
                    glitch_bool_2 = ismember(new_trias,unique(edge_constraint{2}));
                    glitch_bool_3 = ismember(new_trias,unique(edge_constraint{3}));
                    new_trias(all(glitch_bool_1,2) | all(glitch_bool_2,2) | all(glitch_bool_3,2),:)=[];
                    
                end
            else
                %select the interior triangles only
                IO = isInterior(DTP);            %needed?
                new_trias = fliplr(DTP(IO,:));
            end
            
            %remap the nodal ids then store the triangles in the cell structure
            new_triangles = interp1(lib,[triangle(:); selected_nodes_ids],new_trias);
            
            %evaluate in case of very large areas (improper id allocation)
            %temporary patch because I am unable to determine what the
            %error is
              
            %[A]= areaIsosurface(new_triangles,nodes_master_new_tmp);
             
%              if any(A < tolerance )  
%                  stop
%                 figure(1)
%                 clf
%                 hold on
%                 triplot(DTP)
%                 centroid_trias = DTP.Points(DTP.ConnectivityList(:,1),:) + DTP.Points(DTP.ConnectivityList(:,2),:) + DTP.Points(DTP.ConnectivityList(:,3),:);
%                 scatter(centroid_trias(:,1)/3,centroid_trias(:,2)/3)
%                 text(P(:,1),P(:,2),cellstr(num2str((1:size(P,1))')))
%              end
            
            if any(any(isnan(new_triangles))) || any(any(new_triangles==0))
                %add new point and edges to the databases?
                %find the new nodes
                nan_bool = ismembertol(DTP.Points,P,tolerance,'Byrows',true);
                new_nodes = DTP.Points(~nan_bool,:);
                %update the nodes with a height + rotation
                avg_height = mean(tem_surf(:,3));
                node_list = [new_nodes repmat(avg_height,size(new_nodes,1),1)];
                %rotate nodes
                M2 = rotate_vect_by_another([0 0 1],normal); 
                node_list2  = node_list  - centroid_t;
                node_list2  = node_list2*M2;
                node_list2  = node_list2  + centroid_t;
                new_nodes_cell{i} = node_list2;
                
                %find the new constraints
                edge_bool_2 = ismember(sort(DTP.Constraints,2),sort(constraint,2),'rows');
                new_edges = DTP.Constraints(~edge_bool_2,:);
                new_edges_interpd = interp1(lib,[triangle(:); selected_nodes_ids],new_edges);
                new_edges_cell{i} = [new_edges new_edges_interpd];
                
                %remove the triangles with nans
                ttr = new_trias(any(isnan(new_triangles),2),:);
                new_tris_cell{i} = [new_triangles(any(isnan(new_triangles),2),:) ttr];
                new_triangles(any(isnan(new_triangles),2),:) = [];
            end
            
            trias_cell{i} = new_triangles;
        
        else
            trias_cell{i} = triangle;
        end
                
    end
    
    clear DTP
    clear IO
    
    %check the normal orientations of the newly created triangles
    %flip normals that are incorrect
    nodes_master_new_tmp2 = nodes_master_new_tmp;
    nodes_master_new_tmp2(1:3,:) = ones(3);
    %repeat this step until no more triangles were flipped?
    parfor r = 1:length(trias_cell)
        
        old_normal = normals_all(r,:);
        current_tris_an = trias_cell{r};
        TRf_tmp = triangulation(current_tris_an,nodes_master_new_tmp2);    
        new_normals = faceNormal(TRf_tmp);
        
        is_ok = ismembertol(round(new_normals,3),round(old_normal,3),tolerance,'Byrows',true );

        if any(~is_ok)
           current_tris_an(~is_ok,:) = fliplr(current_tris_an(~is_ok,:));
           trias_cell{r} = current_tris_an; 
        end
        
    end
    
    %create new triangles
    trias = cell2mat(trias_cell);
    
    %query out unused cells
    empties = find(cellfun(@isempty,new_nodes_cell));
    new_nodes_cell(empties)=[];
    new_edges_cell(empties)=[];
    new_tris_cell(empties)=[];
    
    num_new_trias_append = size(cell2mat(new_tris_cell),1);
    trias_append = zeros(num_new_trias_append,3);
    edges_append = zeros(size(cell2mat(new_edges_cell),1),2);
    
    if ~isempty(new_nodes_cell)
        %append new triangles with new nodes
        nodal_push = size(nodes_master_new_tmp,1);
        bump = 0;
        bump2 = 0;
        for uy = 1:size(new_nodes_cell,1)
   
            %nodes_tmp_new = new_nodes_cell{uy};
            edges_tmp_new = new_edges_cell{uy};
            trias_tmp_new = new_tris_cell{uy};
            
            %create interpolation library
            tmp_mat_edges_new = edges_tmp_new(:,3:4);
            nan_pos = isnan(tmp_mat_edges_new);
            tmp_mat_edges = edges_tmp_new(:,1:2);
            old_ids = tmp_mat_edges(nan_pos);
            
            trias_old = trias_tmp_new(:,4:6);
            trias_new_curr = trias_tmp_new(:,1:3); 
            %nan_pos_tris = isnan(trias_new_curr);
            
            %get number of unique ids
            total_old_ids = unique(old_ids);
            
            %for each new id, update the traingles and the new edges
            for qw = 1:length(total_old_ids)
                
                nodal_push = nodal_push +1; %id for the next node
                updated_id = nodal_push;
                
                %update the edges
                current_id = total_old_ids(qw);
                sell_bool_id_curr = ismember(tmp_mat_edges, current_id);
                tmp_mat_edges_new(sell_bool_id_curr) = repmat(updated_id,sum(sum(sell_bool_id_curr)),1);
                
                %update the triangles
                sell_bool_id_curr_tris = ismember(trias_old, current_id);
                trias_new_curr(sell_bool_id_curr_tris) = repmat(updated_id,sum(sum(sell_bool_id_curr_tris)),1);

            end
            %append the new triangles to the new_trias matrix
            %update trias_append
            
            rows_append_sel = (1:size(trias_new_curr,1)) + bump;
            bump=bump+length(rows_append_sel);
            trias_append(rows_append_sel,:) = trias_new_curr;
            
            %append the new edges
            rows_append_sel2 = (1:size(tmp_mat_edges_new,1)) + bump2;
            bump2=bump2+length(rows_append_sel2);
            edges_append(rows_append_sel2,:) = tmp_mat_edges_new;
            
        end
        %nodes to be appended at end
        nodes_master_new_tmp = [nodes_master_new_tmp; cell2mat(new_nodes_cell) ];
        
        if edges_on
            %edges to be appedned at the end
            edges_tmp = [edges_tmp; edges_append];
        end
        %new_triangles appended
        trias = [trias; trias_append];
        
    end

    %nodes_master_new_tmp2 = uniquetol(nodes_master_new_tmp(4:end,:),'Byrows',true);
    nodes_master_new_tmp2 = nodes_master_new_tmp(4:end,:);
    
    %renumber the trias
    [~,trias_ids_2] = ismembertol(nodes_master_new_tmp(trias(:),:),nodes_master_new_tmp2,tolerance,'Byrows',true);
    trias(:) = trias_ids_2;
    trias = fliplr(trias);
    
    %renumber the edges and edge list nodes
    if edges_on
        [~,edges_ids_2] = ismembertol(nodes_master_new_tmp(edges_tmp(:),:),nodes_master_new_tmp2,tolerance,'Byrows',true);
        edges_tmp(:) = edges_ids_2;      
        
        nodes_at_edges = unique(cell2mat(nodal_ids_at_edges));
        [~,nodes_at_edges] = ismembertol(nodes_master_new_tmp(nodes_at_edges(:),:),nodes_master_new_tmp2,tolerance,'Byrows',true);
        nodes_at_edges(nodes_at_edges==0)=[]; %remove the ids that could not be found?
    else
        edges_tmp=[];
        nodes_at_edges = [];
    end
        
    %update the master nodes list
    %nodes_master_new_tmp(1:3,:) = [];
    nodes_master_new_tmp = nodes_master_new_tmp2;
    
    
%     figure('Name','Surface Merge Results')
%     clf;
%     view(3)
%     hold on
%     FV_new.faces=trias;
%     FV_new.vertices = nodes_master_new_tmp;
%     fr=triangulation(FV_new.faces,FV_new.vertices);
%     %patch(FV_new,'FaceColor','b','FaceAlpha',0.25)
%     trisurf(fr)
%     g3 = LATTICEfcn.plot_many_line_segments(nodes_master_new_tmp, edges_tmp ) ;
%     set(g3, 'Color','red','LineWidth',3,'Displayname','Lattice Net')
%     axis equal

%debug surface plot for assesing the normal directions of the new mesh
%     figure(2)
% clf;
% view(3)
% hold on
% axis equal
% tmp_tris = cell2mat(trias_cell);
% nodes_master_new_tmp(1:3,1:3) = eye(3);
% TR = triangulation(tmp_tris,nodes_master_new_tmp);    
% trisurf(TR,'edgecolor','black')
% F = faceNormal(TR);
% centroids = (TR.Points(TR.ConnectivityList(:,1),:)+TR.Points(TR.ConnectivityList(:,2),:)+TR.Points(TR.ConnectivityList(:,3),:))/3;
% quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F(:,1),F(:,2),F(:,3),0.5,'Color','b'); %check the face normals, must be consistent
% g6 = LATTICEfcn.plot_many_line_segments(nodes_master_new_tmp, edges_tmp ) ; %cell2mat(edge_loops') %selected_loop 
% set(g6, 'Color','r','LineWidth',3)

    fprintf('Solid Skin Meshing Time : %.3f min\n',toc/60)
    %drip
    %drop
end


