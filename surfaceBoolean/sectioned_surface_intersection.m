%this function will split up the input surfaces to save time during the
%surface intersection and speed up the calcualtions

function [intSurface_tmp,FV_A_tmp, FV_B_tmp,sur_A_int,sur_B_int,edges_co_planar_tmp ,co_planar_switch] = sectioned_surface_intersection(FV_A,FV_B,tolerance)

%this function will calculate the intersection of the intersecting surfaces
%it uses SurfaceIntersection() with additional post-processsing to ensure a
%watertight mesh

%Input
% FV_A : Surface A
% FV_B : Surface B
% tolerance: numeric tolerance for combining nodes/point coordiantes

%Output
% intSurface_tmp : intersected surface structure from SurfaceIntersection()
% FV_A_tmp : temp structure for surface A
% FV_B_tmp : temp structure for surface B
% sur_A_int : intersection boolean for surface A
% sur_B_int : intersection boolean for suface B
% edges_co_planar_tmp : co-planar triangles
% co_planar_switch : if co-planar triangles are found
%%
%Collect intersecting triangles by bounding box to save memory in the
%intersection calculations step
fprintf('Locating Intersection Triangles\n');
[sur_A_int_A , sur_B_int_A] = find_intersecting_triangles_by_sphere(FV_A, FV_B);
[sur_B_int_B , sur_A_int_B] = find_intersecting_triangles_by_sphere(FV_B, FV_A);

%combine the booleans (not the greatest method to do this)
sur_A_int = sur_A_int_A | sur_A_int_B;
sur_B_int = sur_B_int_A | sur_B_int_B;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the surface intersections%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Performing Surface intersection calculations\n');
%reverse = 0;

FV_A_tmp = FV_A;
FV_A_tmp.faces(~sur_A_int,:) = []; 
FV_B_tmp = FV_B;
FV_B_tmp.faces(~sur_B_int,:) = [];
 
%need to split up the triangles into sections so that the amount of memory used is minimized
%(can be parfor'd fto speed up the process)

%might be faster to just loop thorugh each triangle indivudually...(parfor)
try

    %old method, using the entire surface for the surface to surface intersection
    [~, intSurface_tmp] = SurfaceIntersection(FV_A_tmp,FV_B_tmp,'pointroundingtol',1/10e-6);  

    %clean up the intersection once more to avoid duplicate nodes
    nodes_global = uniquetol(intSurface_tmp.vertices,tolerance,'Byrows',true);
    
    %clean up edges
    [~,new_ids_edges] = ismembertol(intSurface_tmp.vertices(intSurface_tmp.edges(:),:),nodes_global,tolerance,'Byrows',true);
    intSurface_tmp.edges(:) = new_ids_edges;
    %remove edges wich have the same node attached twice
    intSurface_tmp.edges(intSurface_tmp.edges(:,1)==intSurface_tmp.edges(:,2),:)=[];
    %clean up faces
    [~,new_ids_faces] = ismembertol(intSurface_tmp.vertices(intSurface_tmp.faces(:),:),nodes_global,tolerance,'Byrows',true);
    intSurface_tmp.faces(:) = new_ids_faces;
    face_removal_bool = intSurface_tmp.faces(:,1)==intSurface_tmp.faces(:,2) | intSurface_tmp.faces(:,2)==intSurface_tmp.faces(:,3) | intSurface_tmp.faces(:,3)==intSurface_tmp.faces(:,1);
    intSurface_tmp.faces(face_removal_bool,:) = [];
    %append the new nodes
    intSurface_tmp.vertices = nodes_global;
        
catch
    %memory limit
    %get limits (x,y,z)
    fprintf('Memory Limit Reached...Splitting Mesh\n');
    vert_tot = [FV_A_tmp.vertices; FV_B_tmp.vertices];
    lims_max_tot = max(vert_tot);
    lims_min_tot = min(vert_tot);
    
    struct_list_init = {FV_A_tmp,FV_B_tmp};
    
    %split along the longest dimension
    limits_tot = (lims_max_tot-lims_min_tot);
    %get position of largest gap
    max_tot_pos = find(limits_tot==max(limits_tot));
    pos_dim = max_tot_pos(1);
    divisions_mesh = 12;
    
    vert_select_array = lims_min_tot(pos_dim):((lims_max_tot(pos_dim)-lims_min_tot(pos_dim))/divisions_mesh):lims_max_tot(pos_dim);
    over_lap_percent = 0.1;
    over_lap = ((lims_max_tot(pos_dim)-lims_min_tot(pos_dim))/divisions_mesh)*(1+over_lap_percent);
    mesh_secs = cell(divisions_mesh,1);
    mesh_nodes = cell(divisions_mesh,1);
    parfor ut = 1:divisions_mesh
        
        vert_select_min = vert_select_array(ut)-over_lap;
        vert_select_max = vert_select_array(ut+1)+over_lap;
        %select the subset of triangles between both meshes
        tmp_struct = struct_list_init;
        for kl = 1:length(struct_list_init)
            current_struct = struct_list_init{kl}; 
            vert_bool = current_struct.vertices(:,pos_dim) <= vert_select_max & current_struct.vertices(:,pos_dim) >= vert_select_min;
            %selct all triangles which share the chosen nodes
            face_tmp = current_struct.faces;
            face_tmp(:) = vert_bool(face_tmp(:));
            tri_bool = any(face_tmp,2);
            tmp_struct{kl}.faces = current_struct.faces(tri_bool,:);
        end
        %perform the intersection calculations if triangles are selected
        if ~any(cell2mat(cellfun(@(d) isempty(d.faces) ,tmp_struct,'UniformOutput',false)))
            
            [~, intSurface_tmp_tmp] = SurfaceIntersection(tmp_struct{1},tmp_struct{2},'pointroundingtol',1/tolerance);  
            
            %save the output
            mesh_secs{ut} = intSurface_tmp_tmp;
            mesh_nodes{ut} = intSurface_tmp_tmp.vertices;
        end
        
    end
    
    %create the global nodes list
    nodes_global = uniquetol([vert_tot; cell2mat(mesh_nodes)],tolerance,'Byrows',true);
    %filter out empty cells
    filt_bool = cell2mat(cellfun(@isempty,mesh_secs,'UniformOutput',false));
    mesh_secs = mesh_secs(find(~filt_bool));
    filt_bool = cell2mat(cellfun(@(d) isempty(d.faces) ,mesh_secs,'UniformOutput',false));
    mesh_secs = mesh_secs(find(~filt_bool));
    
    %re-assemble the intersection surfaces
    edge_list_tot = cell(length(mesh_secs),1);
    trian_list_tot = cell(length(mesh_secs),1);
    parfor h = 1:length(mesh_secs)
    
        mesh_curr = mesh_secs{h};
        
        %update the structure so that all the edges have the same vertices
        %numbering as before..
        %for the faces
        [~,new_ids_faces] = ismembertol(mesh_curr.vertices(mesh_curr.faces(:),:),nodes_global,tolerance,'Byrows',true);
        mesh_curr.faces(:) = new_ids_faces;
        trian_list_tot{h} = reshape(mesh_curr.faces,numel(mesh_curr.faces)/3,3);
        %for the edges
        [~,new_ids_edges] = ismembertol(mesh_curr.vertices(mesh_curr.edges(:),:),nodes_global,tolerance,'Byrows',true);
        mesh_curr.edges(:) = new_ids_edges;
        %make sure that the number of columns are consisttent
        edge_list_tot{h} =  reshape(mesh_curr.edges,numel(mesh_curr.edges)/2,2);
        
    end
    
    %finalize calculations
    edge_list_tot = unique(cell2mat(edge_list_tot),'rows');
    trian_list_tot = unique(cell2mat(trian_list_tot),'rows');
    intSurface_tmp.vertices = nodes_global;
    intSurface_tmp.edges = edge_list_tot;
    intSurface_tmp.faces = trian_list_tot;
end



%% determine if co-planarity is found between both objects
int_faces = intSurface_tmp.faces;
co_p_faces = ~any (int_faces(:,1) == int_faces(:,2) | int_faces(:,2) == int_faces(:,3) | int_faces(:,3) == int_faces(:,1) ,2);
co_planar_switch = 0;
if any(co_p_faces)

    %there is co-planerity
    fprintf('Co-Panar Surfaces Detected! Results may be erronous.\n');

    co_planar_switch = 1;
    
    %extract the edges from the faces and use those as a closed loop
    critical_faces = int_faces(co_p_faces,:);
    
    %get the free edges in the critical_faces matrix
    crit_test.faces = critical_faces;
    [~,edge_list_crit_test] = test_mesh_watertight(crit_test,1);
    edges_co_planar_tmp = edge_list_crit_test(edge_list_crit_test(:,end)==1,1:2);
    
    [~,loc_of_true_edge] = ismember(sort(edges_co_planar_tmp,2),sort(intSurface_tmp.edges,2),'rows');
    edges_co_planar_tmp = unique(intSurface_tmp.edges(loc_of_true_edge,:),'rows','stable');
    %get the edge loops for the co_planaer areas to determine the correct iso-surfaces
    [~,edge_loops_co_planar]  = form_intersection_loops_boolean(edges_co_planar_tmp);
    
    %remove edge loops which have less than 3 connected edges
    loop_sizes = cellfun(@(dd) size(dd,1),edge_loops_co_planar,'UniformOutput',false); 
    cell_bools_1 = cellfun(@(dc) dc<3,loop_sizes,'UniformOutput',false);
    edge_loops_co_planar(find(cell2mat(cell_bools_1'))) = [];  
    edges_co_planar = cell2mat(edge_loops_co_planar');
        
    %remove trias that are not needed, as well as any of there edges
    non_crit_faces = intSurface_tmp.faces(~co_p_faces,:);
    %extract the new edges from the faces
    edge_kb1 = non_crit_faces(:,1)==non_crit_faces(:,2);
    new_edges_part1 = non_crit_faces(edge_kb1,[1 2]);
    edge_kb2 = non_crit_faces(:,2)==non_crit_faces(:,3);
    new_edges_part2 = non_crit_faces(edge_kb2,[2 3]);
    edge_kb3 = non_crit_faces(:,1)==non_crit_faces(:,3);
    new_edges_part3 = non_crit_faces(edge_kb3,[3 1]);
    %group the triangle edges together
    new_tri_edge_from_faces = [new_edges_part1; new_edges_part2; new_edges_part3];
    %remove bad triangles
    intSurface_tmp.faces(~co_p_faces,:) = [];
    %update the edges field in the structure
    old_edges = intSurface_tmp.edges;
    %remove edges associated with the co-planaer triangles
    co_p_edges = unique([critical_faces(:,[1 2]); critical_faces(:,[2 3]); critical_faces(:,[3 1])],'rows','stable');
    rem_bool_ed1 = ismember(sort(old_edges,2),sort(co_p_edges,2),'rows');
    old_edges(rem_bool_ed1 ,:)=[];
    %remove any old edges that form bad loops (same as with the free-edge co-planar loops)
    try
        old_edges_loops = form_intersection_loops_boolean(old_edges);
        old_loop_sizes = cellfun(@(dd) size(dd,1),old_edges_loops,'UniformOutput',false); 
        cell_bools_2 = cellfun(@(dc) dc<3,old_loop_sizes,'UniformOutput',false);
        old_edges_loops(find(cell2mat(cell_bools_2'))) = [];  
        old_edges = cell2mat(old_edges_loops');
    catch
       %old edges are not update
       fprintf('Additional Errors Detected\n');
    end
    
    new_edge_list = unique([new_tri_edge_from_faces; old_edges; edges_co_planar_tmp],'rows');
    intSurface_tmp.edges = new_edge_list;
else

    edges_co_planar_tmp = [];
    
end

end