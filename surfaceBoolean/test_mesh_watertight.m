%this function will test if a mesh is water_tight

function [answer, edge_matrix, answer2, edge_matrix_directed] = test_mesh_watertight(FV,option)
%Input
%FV: a structure with fields faces and vertices
%option
%if 1, then analyze the triangles only
%if 2, the analyze the edge list only
%else then both edges and triangles(faces) will be analyzed

%extract the faces from the triangles if requested

%Output
% answer : boolean indicator if the mesh is watertight (no free edges)
% edge_matrix : edge connectivity table and number of occurances of that edge in the mesh
% answer2 : if all edges are shared only twice
% edge_matrix_directed :same as edge_matrix but for dge connectivity table that does not sort the edges (direction maintained)

%determine if a triangulation class is used or not
type = class(FV);
switch type
    case 'triangulation'
        FV2.faces = FV.ConnectivityList;
        FV2.vertices = FV.Points;
        FV=FV2;
        clear FV2
    case 'struct'
        %proceed as normal
    otherwise
        error('Must input a Structure with fields faces and vertices or an object of class triangulation.')
end

if option~=2
    trias_FV_tmpA = FV.faces;
    %tria_edges = sort([trias_FV_tmpA(:,1:2); trias_FV_tmpA(:,2:3); trias_FV_tmpA(:,[3 1])],2);
    tria_edges_directed = [trias_FV_tmpA(:,1:2); trias_FV_tmpA(:,2:3); trias_FV_tmpA(:,[3 1])];
    tria_edges_directed_unique = unique(tria_edges_directed,'rows','stable');
    tria_edges_sort = sort(tria_edges_directed,2);
    tria_edges_sort_unique = unique(tria_edges_sort,'rows'); 
else
    tria_edges_sort = [];
    tria_edges_sort_unique = [];
    tria_edges_directed = [];
    tria_edges_directed_unique = [];
end
        
if option~=1
    edge_list_directed  = FV.edges;
    edge_list_directed_unique = unique(edge_list_directed,'rows');
    edge_list_sort = sort(edge_list_directed,2);
    edge_list_sort_unique = unique(edge_list_sort,'rows'); 
else
    edge_list_sort =[];
    edge_list_sort_unique = [];
    edge_list_directed = [];
    edge_list_directed_unique = [];
end

%add up the edge lists for the unique edges (non-directed) (is a soft
%constraint on manifold meshes)
total_edge_list_sort_unique = [tria_edges_sort_unique; edge_list_sort_unique];
total_edge_list_sort = [tria_edges_sort; edge_list_sort];
[~, ~, pNos] = unique(total_edge_list_sort,'rows');
hm = accumarray(pNos,1);
answer = min(hm) ==2;
edge_matrix = [total_edge_list_sort_unique hm];

%for the directed case
total_directed_edge_list = [tria_edges_directed; edge_list_directed];
total_directed_edge_list_unique = [tria_edges_directed_unique; edge_list_directed_unique];
[~, ~, pNos2] = unique(total_directed_edge_list,'rows');
hm2 = accumarray(pNos2,1);
answer2 = min(hm2) ==2;
edge_matrix_directed = [total_directed_edge_list_unique hm2];

end