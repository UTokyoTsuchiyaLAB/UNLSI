%this function is for the bounding box on triagnles between 2 surfaces
%needs to be updated...to be spheres?

function [sur_A_int , sur_B_int] = find_intersecting_triangles_by_sphere(FV_A, FV_B)

%surface B intersecting surface A
tolerance = 10^-9;

surfaceA_tmp = FV_A;
surfaceB_tmp = FV_B;
triA_tmp = surfaceA_tmp.faces;
triB_tmp = surfaceB_tmp.faces;
nodes_surfA = surfaceA_tmp.vertices;
nodes_surfB = surfaceB_tmp.vertices;
sur_A_int = repmat({0},size(triA_tmp,1),1);
sur_B_int = repmat({0},size(triB_tmp,1),1);

centroids_A = (nodes_surfA(triA_tmp(:,1),:) + nodes_surfA(triA_tmp(:,2),:) +nodes_surfA(triA_tmp(:,3),:)  )/3;
centroids_B = (nodes_surfB(triB_tmp(:,1),:) + nodes_surfB(triB_tmp(:,2),:) +nodes_surfB(triB_tmp(:,3),:)  )/3;

%get the sides of all the triangles in A
AA = sqrt(sum((nodes_surfA(triA_tmp(:,1),:) - nodes_surfA(triA_tmp(:,2),:)).^2,2));
BB = sqrt(sum((nodes_surfA(triA_tmp(:,2),:) - nodes_surfA(triA_tmp(:,3),:)).^2,2));
CC = sqrt(sum((nodes_surfA(triA_tmp(:,3),:) - nodes_surfA(triA_tmp(:,1),:)).^2,2));
circumradius_A = (AA.*BB.*CC)./sqrt((AA+BB+CC).*(BB+CC-AA).*(AA+BB-CC));
sur_b_tri_ids = repmat({[]},size(triB_tmp,1),1);

parfor ux = 1:size(triA_tmp,1)

    %test if there are any nodes from surf B inside the triangles bounding
    %box
    %cur_tri = triA_tmp(ux,:);
    %nodes_curr = nodes_surfA(cur_tri,:);
    center_a = centroids_A(ux,:);
    cur_rad = circumradius_A(ux,:)*3; %radius is is the only vairable that determines the robustness...
    max_vals = center_a+cur_rad;
    min_vals = center_a-cur_rad;
    
    nb1 = nodes_surfB(:,1) <= max_vals(1);
    nb2 = nodes_surfB(:,2) <= max_vals(2);
    nb3 = nodes_surfB(:,3) <= max_vals(3);
    
    nb4 = nodes_surfB(:,1) >= min_vals(1);
    nb5 = nodes_surfB(:,2) >= min_vals(2);
    nb6 = nodes_surfB(:,3) >= min_vals(3);
    
    nb7 = centroids_B(:,1) <= max_vals(1);
    nb8 = centroids_B(:,2) <= max_vals(2);
    nb9 = centroids_B(:,3) <= max_vals(3);
    
    nb10 = centroids_B(:,1) >= min_vals(1);
    nb11 = centroids_B(:,2) >= min_vals(2);
    nb12 = centroids_B(:,3) >= min_vals(3);
    
    %get nodal indeces using a logical boolean
    glo_b = nb1 & nb2 & nb3 & nb4 & nb5 & nb6;
    centers_b = nb7 & nb8 & nb9 & nb10 & nb11 & nb12;
        
    if any(any(glo_b) | any(centers_b))
        %link up the nodal ids with the triangles
        nodal_ids_tmp = find(glo_b);
        tri_sel_bool = any(ismember(triB_tmp,nodal_ids_tmp),2) | centers_b;
        sur_A_int{ux} = 1;
        nums_tmp = find(tri_sel_bool);
        sur_b_tri_ids{ux,1} = nums_tmp(:);
        %(find(tri_sel_bool))  = repmat({1},sum(tri_sel_bool),1);
    end
    
end

triBs_ids = unique(cell2mat(sur_b_tri_ids));

sur_A_int = cell2mat(sur_A_int);
sur_B_int(triBs_ids) = repmat({1},length(triBs_ids),1);
sur_B_int = cell2mat(sur_B_int);

end