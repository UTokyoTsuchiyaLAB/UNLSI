%this will perform a boolean on cube and sphere
close all
clear all
clear path
clc
%turn off warnings
warning('off','all')
%tolerance for nodes
tolerance = 10^-8; %floating point tolerance
visual=1; %plot the floodfill algorithm? (your choice 1=yes,0=no)

%Surface boolean between volumes A and B

%create a cube
side_length = 10;
side_array = -(side_length/2):(side_length/30):(side_length/2);
[xc,yc,zc] = meshgrid(side_array,side_array,side_array);
xc = xc(:);
yc = yc(:);
zc = zc(:);

[k1,av1] = convhull(xc,yc,zc,'Simplify',false);

%create a sphere
[XS,YS,ZS] = sphere(21);
R = 6; %radius of the sphere
XS = XS(:)*R;
YS = YS(:)*R;
ZS = ZS(:)*R;
%scatter3(x,y,z)
[k2,av2] = convhull(XS,YS,ZS,'Simplify',false);


%Plot the cube
figure(1)
clf
hold on
surf_c = trisurf(k1,xc,yc,zc,'FaceColor','cyan');
axis equal
view(3)
hold off

%plot the sphere
figure(2)
clf
hold on
surf_s = trisurf(k2,XS,YS,ZS,'FaceColor','cyan');
axis equal
view(3)
hold off


%contruct the structure to be used for surface booleans
FV_A.faces = surf_c.Faces;
FV_A.vertices = surf_c.Vertices;

FV_B.faces = surf_s.Faces;
FV_B.vertices = surf_s.Vertices;

%plot the shapes together, ensure normals are pointing in the same direction
figure(3)
clf
hold on
view(3)
axis equal
Triang_A = triangulation(FV_A.faces,FV_A.vertices);
trisurf(Triang_A,'edgecolor','k','Facecolor','y')
F = featureEdges(Triang_A,pi/10);
F2 = faceNormal(Triang_A);
centroids = (Triang_A.Points(Triang_A.ConnectivityList(:,1),:)+Triang_A.Points(Triang_A.ConnectivityList(:,2),:)+Triang_A.Points(Triang_A.ConnectivityList(:,3),:))/3;
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'Color','b'); %check the face normals, must be consistent

% axis equal
Triang_B = triangulation(FV_B.faces,FV_B.vertices);
trisurf(Triang_B,'edgecolor','k','Facecolor','g')
F = featureEdges(Triang_B,pi/10);
F2 = faceNormal(Triang_B);
centroids = (Triang_B.Points(Triang_B.ConnectivityList(:,1),:)+Triang_B.Points(Triang_B.ConnectivityList(:,2),:)+Triang_B.Points(Triang_B.ConnectivityList(:,3),:))/3;
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'Color','b'); %check the face normals, must be consistent

xlabel('x')
ylabel('y')
zlabel('z')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Mesh Quality %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%test if the input meshes are closed or not

logC = test_mesh_watertight(FV_B,1);
logS = test_mesh_watertight(FV_A,1);

%must loop for each free volume in "Sphere"
if ~logC || ~logS
    fprintf('Warning! Initial Meshes are not Watertight!\n');
end

[intSurface_tmp,FV_A_tmp, FV_B_tmp,sur_A_int,sur_B_int,edges_co_planar_tmp ,co_planar_switch] = sectioned_surface_intersection(FV_A,FV_B,tolerance);

%cleaned intersection line segments
edge_list_init = intSurface_tmp.edges; %edges seem to be directed properly if all the triangles have the same normals

%% plot the contours

hold on; 
view(3)
%plot the directions of the edges
start_points = intSurface_tmp.vertices(edge_list_init(:,1),:);
directions = intSurface_tmp.vertices(edge_list_init(:,2),:) - start_points;
directions = directions./sqrt(sum(directions.^2,2));
quiver3(start_points(:,1),start_points(:,2),start_points(:,3),directions(:,1),directions(:,2),directions(:,3),1,'Color','red'); 
axis equal

%%
%%%%%%%%%%%%%%%%%%%
%%%Clean Meshes %%%
%%%%%%%%%%%%%%%%%%%

%clean meshes before connecting intersection edges together
[surfaceA, surfaceB] = pre_process_mesh(intSurface_tmp,FV_A, FV_A_tmp ,FV_B, FV_B_tmp,sur_A_int,sur_B_int, edge_list_init,edges_co_planar_tmp,tolerance);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FROMING INTERSECTION LOOPS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Forming Intersection Loops\n');
%for soft contour loops (intersecting contours)
%the loops must be split into two parts
%project lines to have and 'x' shape
all_edges = unique([surfaceA.edges; surfaceB.edges],'rows','stable');
[edge_loops,edge_loops_directed]  = form_intersection_loops_boolean(all_edges);


%%
%%%%%%%%%%%%%%%%%%%%%%%
%CREATING SUB SURFACES%
%%%%%%%%%%%%%%%%%%%%%%%
[m_aBlocks, orientation_matrix_tot ] = create_sub_blocks(surfaceA,surfaceB,edge_loops,co_planar_switch,visual);

%%
%%%%%%%
%Update the boolean when co-planar is turned on
%%%%%%%
if co_planar_switch
     
     m_aBlocks{1}.iso_bool(sum(orientation_matrix_tot,2)==0,:) = 0; 
     m_aBlocks{2}.iso_bool(sum(orientation_matrix_tot,2)==0,:) = 0;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine if the subblocks are union, subtraction or intersection%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%union = A+B-intersection
%subtraction = union - 1st or 2nd Shape
%intersection = seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Solutions and Ensure Proper Orientations and Connectivity%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_titles(1) = "A-B"; 
fig_titles(2) = "B-A"; 
fig_titles(3) = "Intersection"; 
fig_titles(4) = "Union"; 

new_nodes = uniquetol(surfaceA.vertices,tolerance,'Byrows',true);

%plot the new Shapes
    figure(4)
    clf;
    t = tiledlayout(2,2);
    
for yyu = 1:length(m_aBlocks)
    tris_block = m_aBlocks{yyu}.trias;
    nodes_block = surfaceA.vertices(tris_block(:),:);
    [~,new_ids] = ismembertol(nodes_block,new_nodes,tolerance,'Byrows',true );
    tris_block(:) = new_ids;
    m_aBlocks{yyu}.trias = tris_block;
    
    nexttile
    
    title(fig_titles(yyu))
    axis equal
    hold on
    upper = triangulation(m_aBlocks{yyu}.trias,surfaceA.vertices);
    trisurf(upper,'edgecolor','k')
    F = featureEdges(upper,pi/10);
    view(3)
    axis equal
end

