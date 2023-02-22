
function [surfaceA, surfaceB] = pre_process_mesh(intSurface_tmp,FV_A, FV_A_tmp ,FV_B, FV_B_tmp,sur_A_int,sur_B_int, edge_list_init,edges_co_planar_tmp,tolerance)

%this function will perform post processing on the intersected surface to
%ensure a water tight mesh

%Input
% intSurface_tmp : temp intersection structure from SurfaceIntersection()
% FV_A_tmp : temp surface A
% FV_B_tmp : temp surface B
% edge_list_init : edge list for intersections between surfaces
% edges_co_planar_tmp : coplanar edges
% tolerance : point coordinate tolerance (floating point correction)

%Output
% surfaceA : cleaned surface FV_A_tmp
% surfaceB : cleaned surface FV_B_tmp

fprintf('Collecting Information on Intersections\n');
%to remesh both boundaries

edge_nodes_ids_tmp = unique(edge_list_init);
surface_nodes = intSurface_tmp.vertices(edge_nodes_ids_tmp ,:);

fprintf('Remeshing\n');

%this function will create an erronous edge list when soft contours are supplied
[trias_FV_tmpA,nodes_master_new_FV_tmpA,new_edges_tmpA,nodes_at_edges_A] = nodes_to_surface_merge(intSurface_tmp.vertices,surface_nodes, FV_A_tmp , edge_list_init,tolerance); %check for errors here
[trias_FV_tmpB,nodes_master_new_FV_tmpB,new_edges_tmpB,nodes_at_edges_B] = nodes_to_surface_merge(intSurface_tmp.vertices,surface_nodes, FV_B_tmp , edge_list_init,tolerance);

%create new structures for each of the surfaces
surfaceA.faces = trias_FV_tmpA;
surfaceA.vertices = nodes_master_new_FV_tmpA;
surfaceA.edges = new_edges_tmpA;

surfaceB.faces = trias_FV_tmpB;
surfaceB.vertices = nodes_master_new_FV_tmpB;
surfaceB.edges = new_edges_tmpB;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combine the old triangles with the new ones%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Combining and Unifying Meshes\n');

%The first surface
FV_A_tmp = FV_A;
FV_A_tmp.faces(find(sur_A_int),:) = [];
surfaceA.faces = fliplr(surfaceA.faces);%flip faces to return to normal orientations
surfaceA = combine_surfaces(surfaceA,FV_A_tmp); 
%the second surface
FV_B_tmp = FV_B;
FV_B_tmp.faces(find(sur_B_int),:) = [];
surfaceB.faces = fliplr(surfaceB.faces);
surfaceB = combine_surfaces(surfaceB,FV_B_tmp);

%unify the surface meshes together now
all_vertices = uniquetol([surfaceA.vertices; surfaceB.vertices],tolerance,'Byrows',true);

%renumber the faces and edges
[~,loc_face_A] = ismembertol(surfaceA.vertices(surfaceA.faces(:),:),all_vertices,tolerance,'Byrows',true);
surfaceA.faces(:) = loc_face_A;
[~,loc_face_B] = ismembertol(surfaceB.vertices(surfaceB.faces(:),:),all_vertices,tolerance,'Byrows',true);
surfaceB.faces(:) = loc_face_B;
[~,loc_edge_A] = ismembertol(surfaceA.vertices(surfaceA.edges(:),:),all_vertices,tolerance,'Byrows',true);
surfaceA.edges(:) = loc_edge_A;
[~,loc_edge_B] = ismembertol(surfaceB.vertices(surfaceB.edges(:),:),all_vertices,tolerance,'Byrows',true);
surfaceB.edges(:) = loc_edge_B;

%update the co-planar edge vertices ids
[~,loc_edge_co_planar] = ismembertol(intSurface_tmp.vertices(edges_co_planar_tmp(:),:),all_vertices,tolerance,'Byrows',true);
edges_co_planar_tmp(:) = loc_edge_co_planar; 

%scatter3(all_vertices(edges_co_planar_tmp(:),1),all_vertices(edges_co_planar_tmp(:),2),all_vertices(edges_co_planar_tmp(:),3))

%update the ids of the nodes located at the original triangle edges
[~,nodes_at_edges_A_updated] = ismembertol(surfaceA.vertices(nodes_at_edges_A,:),all_vertices,tolerance,'Byrows',true);
[~,nodes_at_edges_B_updated] = ismembertol(surfaceB.vertices(nodes_at_edges_B,:),all_vertices,tolerance,'Byrows',true);

nodes_ids_at_edges_total = unique([surfaceA.edges(:); surfaceB.edges(:)]);

%scatter3(all_vertices(nodes_ids_at_edges_total(:),1),all_vertices(nodes_ids_at_edges_total(:),2),all_vertices(nodes_ids_at_edges_total(:),3))

%allocate the new nodes to the surfaces
surfaceA.vertices = all_vertices;
surfaceB.vertices = all_vertices;


%checks
%test if the new triangles are watertight (no holes)
%all edges must be shared at least twice
%no infinitesimally small area triangles
%[AA]= areaIsosurface(surfaceA.faces,surfaceA.vertices);
%[AB]= areaIsosurface(surfaceB.faces,surfaceB.vertices);

if  test_mesh_watertight(surfaceA,1) && test_mesh_watertight(surfaceB,1) 
    fprintf('Meshes OK!\n');
else
    fprintf('Error in Meshes!\n');

    if ~test_mesh_watertight(surfaceA,1) || ~test_mesh_watertight(surfaceB,1)
        fprintf('Surfaces are not Water Tight\n');
        %stop
        %need to only consider nodes that are on the original edges of the
        %original mesh
        
        %stop
        
        %close element gaps if they exist
        [~,list_edges_A] =  test_mesh_watertight(surfaceA,1);
        old_error_node_ids_A = unique(list_edges_A(list_edges_A(:,end)~=2,1:2));
        error_nodes_A = surfaceA.vertices(old_error_node_ids_A,:);
        [~,list_edges_B] =  test_mesh_watertight(surfaceB,1);
        old_error_node_ids_B = unique(list_edges_B(list_edges_B(:,end)~=2,1:2));
        error_nodes_B = surfaceB.vertices(old_error_node_ids_B,:);
        error_nodes_total = [error_nodes_A; error_nodes_B];
        error_nodes_ids_total = [old_error_node_ids_A; old_error_node_ids_B];
        %get the edge list of a that contain error nodes
        surf_a_bad_edges_bool = ismember(surfaceA.edges,error_nodes_ids_total);
        surf_a_bad_edges_list = surfaceA.edges(all(surf_a_bad_edges_bool,2),:);
        %filter out uneeded nodes
        filter_bool = ~ismember(error_nodes_ids_total,nodes_ids_at_edges_total); %nodes_ids_at_edges_total
        error_nodes_ids_total(filter_bool,:)=[];
        error_nodes_total(filter_bool,:) = [];
        %get a list of all the free edges and group them toegether as loops
        if ~isempty(old_error_node_ids_A)
            loops_of_A = form_intersection_loops_boolean(list_edges_A(list_edges_A(:,end)<2,1:2));
        else
           loops_of_A = []; 
        end
        if ~isempty(old_error_node_ids_B)
            loops_of_B = form_intersection_loops_boolean(list_edges_B(list_edges_B(:,end)<2,1:2));
        else
            loops_of_B = [];
        end
        total_loops = [loops_of_A loops_of_B];
        %add triangles or combine???
        library = zeros(length(total_loops),2);
        for jjuy = 1:length(total_loops)
            
            current_edge_group = total_loops{jjuy};            
            %equivalence nodal ids for the edge nodes then update the surface structure
            part_g = ismember(current_edge_group,error_nodes_ids_total);
            
            id_collect = unique(current_edge_group(part_g));
            min_id = min(id_collect);
            
            other_id = id_collect(id_collect~=min_id);
            if numel(other_id) > 1
                fprintf('Warning Errors in Mesh Detected!\n')
            end
            if ~isempty(other_id)
                surfaceA.faces(ismember(surfaceA.faces,other_id)) = min_id;
                surfaceB.faces(ismember(surfaceB.faces,other_id)) = min_id;
                surfaceA.edges(ismember(surfaceA.edges,other_id)) = min_id;
                surfaceB.edges(ismember(surfaceB.edges,other_id)) = min_id;
                library(jjuy,:) = [min_id other_id(1)];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %remove bad triangles.....
        
        %remove triangles and edges that share the same node twice (needed?)
        bad_tri_bool_A = surfaceA.faces(:,1)==surfaceA.faces(:,2) | surfaceA.faces(:,2)==surfaceA.faces(:,3) | surfaceA.faces(:,3)==surfaceA.faces(:,1);
        surfaceA.faces(bad_tri_bool_A,:) = [];
        
        bad_tri_bool_B = surfaceB.faces(:,1)==surfaceB.faces(:,2) | surfaceB.faces(:,2)==surfaceB.faces(:,3) | surfaceB.faces(:,3)==surfaceB.faces(:,1);
        surfaceB.faces(bad_tri_bool_B,:) = [];
        
        %surf A edges
        surfaceA.edges(surfaceA.edges(:,1)==surfaceA.edges(:,2),:)=[];
        %surf B edges
        surfaceB.edges(surfaceB.edges(:,1)==surfaceB.edges(:,2),:)=[];
        
        %evaluate if the surfaces are closed or not
        is_A_ok = test_mesh_watertight(surfaceA,1);
        is_B_ok = test_mesh_watertight(surfaceB,1);
        
        if ~(is_A_ok || is_B_ok)
            error('Error in Mesh stiching')
            %error in the mesh
        end
    end
end


end