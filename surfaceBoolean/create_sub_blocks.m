
function [m_aBlocks, orientation_matrix_tot ] = create_sub_blocks(surfaceA,surfaceB,edge_loops,co_planar_switch,visual)


%this surface will organize the faces into volumes
%intersection information is used to determine a seed to begin a floodfill
%algorithm (visual=1 to see it)
%an orientation matrix is created to determine which type of boolean
%(intersetion, subtraction, etc) is asscociated with each new volume.
%all volume are placed in m_aBlocks in a particular order (A-B, B-A, Intersection, Union)

%Input
% surfaceA : surface A
% surfaceB : surface B
% edge_loops : edge loops for intersections
% co_planar_switch : coplanar triangles are found
% visual : show floodfill algortihm

%Output
% m_aBlocks: structure for each boolean operation
% orientation_matrix_tot : boolean matrix for the intersections

%%%%%%%%%%%%%%%%%%%%%%%
%CREATING SUB SURFACES%
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Creating Sub Surfaces\n');
%FV = surfaceB;
%this code will not pick up interior holes? - > will be lumped into the final structure...

%this is a flood fill algortimh, it is very sensitive to erros in the mesh
[newSF_B, decision_matrix_B,orientation_array_B] = create_sub_surfaces_boolean(surfaceB, edge_loops,visual); 
[newSF_A, decision_matrix_A,orientation_array_A] = create_sub_surfaces_boolean(surfaceA, edge_loops,visual);

%not picking up soft loops...
if sum(sum(abs(orientation_array_A)))==0
    error('Triagulated Mesh Boolean Failed\n');   
end
%begin assembling the sub blocks into union, subtration and intersection
%matrix is isosurface vs. contour present inside it

%[decision_matrix_A orientation_array_A]
%[decision_matrix_B orientation_array_B]

%%

fprintf('Distinguishing Iso-Surfaces\n');

%there can be multiple public sub-surfaces (3) which have more than one
%close loop...
%need to code soft-boundaries (shapes with the same thicknesses (e.g; 2 cylinders))

%a private loop can only be attached to a public of another loop (for hollow cases)
%how to deal with surfaces that are not affected by the intersection?

%combine newSF_A with newSF_B befre this step?
newSF_tot = cell(length(newSF_A)+length(newSF_B),1);
newSF_tot(1:length(newSF_A)) = newSF_A;
for ws = 1:length(newSF_A) 
    newSF_tot{ws}.surf_origin = 0;
end
countb=0;
for qqw = (length(newSF_A)+1):length(newSF_A)+length(newSF_B)
    countb=countb+1;
    newSF_tot{qqw}.trias = newSF_B{countb}.trias;
    newSF_tot{qqw}.edges = newSF_B{countb}.edges;
    newSF_tot{qqw}.orientation = newSF_B{countb}.orientation;
    newSF_tot{qqw}.surf_origin = 1;
end

%combine the dcision matrices
decision_matrix_tot = [decision_matrix_A; decision_matrix_B];

%combine the orienation matrices
orientation_matrix_tot = [orientation_array_A; orientation_array_B];

%surface_origin_matrix
surface_origin_matrix = [ones(size(orientation_array_A,1),1);zeros(size(orientation_array_B,1),1)];

%for this case, a public can only be selected once?
if sum(sum(orientation_matrix_tot))~=0
    error('Boolean Failed, Iso Surfaces can not be distinguished');
    stop %error?
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%SUBTRACTIONS FOR A AND B%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%combine surfaces together to obtain both of the possible subtractions
%lump together all surfaces with equal orientations and equal number of shared loops
%fill up gaps with private surfaces during construction (when there are
%single open loops during matching surfaces with equal number of closed loops)

%A-B = all isosurfaces with at least a negative orientation closed loop
%B-A = all isosurface with at least a positive orientaiton closed loop
%both subtractions do not share any iso-surfaces together

%The first isosurface created is A - B (flip the normals of B) (B=0)
%The second Isosurface created is B - A (Flip the normals of A) (A=1)

block_store = cell(2,1);
prev_s_sum = [];
pos_rows = any(orientation_matrix_tot==1,2);  
neg_rows = any(orientation_matrix_tot==-1,2);
common_orients = [pos_rows neg_rows];
s_row = sum(abs(orientation_matrix_tot),2);
num_loop_groups = unique(s_row);
num_isos = size(orientation_matrix_tot,1);
for fk = 1:length(block_store)
    m_aBlocks{fk}.trias = [];
    %look in the orientation matrices

    if fk==1
        current_row = neg_rows;
    else
        current_row = pos_rows;
    end
    
    select_r_num = find(current_row);
    
    %evaluate if the matrix has all loops used and two surface sharing
    %each loop
    eval_matrix = decision_matrix_tot(select_r_num,:);
    
    all_loops_closed = sum(eval_matrix); %each columns must be equal to 2
    
    %find faces to flip
    flip_case = surface_origin_matrix(select_r_num);
    
    %assemble the sub block
    for huk = 1:length(select_r_num)
        %check if the selected iso_surface is a co-planar surface, if it
        %is, omit it, test if the suraface shares all edges with it
        %edges_co_planar_tmp
        pass=0;
        if co_planar_switch
            %current_tiangles_to_add = newSF_tot{select_r_num(huk)}.trias;
            %current_tiangles_to_add_edges = unique(sort([current_tiangles_to_add(:,[1 2]); current_tiangles_to_add(:,[2 3]); current_tiangles_to_add(:,[3 1])],2),'rows','stable');
            %co_planar_bool = ismember(current_tiangles_to_add_edges, sort(edges_co_planar_tmp,2),'rows');
            %test if all orientaions sum up to zero for that loop, if so,
            %then that loop needs to be omitted from the subtraction sub-blocks
            
            if sum(orientation_matrix_tot(select_r_num(huk),:)) ~=0 %~all(co_planar_bool)
                pass=1;
            end
        else
           pass = 1; 
        end
        if pass
            if (fk==1 && flip_case(huk)==0) || (fk==2 && flip_case(huk)==1)
                %flip all faces associated with B (B=0) or
                trias_add = fliplr(newSF_tot{select_r_num(huk)}.trias);
            else
                %do not flip the faces
                trias_add = newSF_tot{select_r_num(huk)}.trias;
            end
            
            m_aBlocks{fk}.trias =  [m_aBlocks{fk}.trias; trias_add];    
        end
    end
    %create and append the iso_bool array
    m_aBlocks{fk}.iso_bool = ismember(1:num_isos,select_r_num)';
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get The Union from Boolean Logic%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%using the subtraction
%The union can be obtained from either:
%Case 1: Union = A - subtraction iso_boolean_A_minus_B
%Case 2: Union = B - subtraction iso_boolean_B_minus_A
%then select the surfaces wich are turned off

%for case 1
fk=3;
m_aBlocks{fk}.trias = [];

%do surface A minus the subtraction
sub_bool_1 = m_aBlocks{1}.iso_bool;

%ISO_BOOL Suraface A
iso_bool_surf_A = [ones(size(orientation_array_A,1),1); zeros(size(orientation_array_B,1),1)];
%Iso_Bool Surface B
iso_bool_surf_B = [zeros(size(orientation_array_A,1),1); ones(size(orientation_array_B,1),1)];
%assemble the sub block
int_bool = iso_bool_surf_A - sub_bool_1;
select_r_num_int1 = find(int_bool~=-0);

for huk = 1:length(select_r_num_int1)
   trias_add = newSF_tot{select_r_num_int1(huk)}.trias;
   m_aBlocks{fk}.trias =  [m_aBlocks{fk}.trias; trias_add]; 
end
m_aBlocks{fk}.iso_bool = int_bool~=-0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the intersection as the inverse of the union...%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fk=4;
m_aBlocks{fk}.trias = [];
m_aBlocks{fk}.iso_bool = ~m_aBlocks{3}.iso_bool;
select_r_num_int2 = find(m_aBlocks{fk}.iso_bool);

for huk = 1:length(select_r_num_int2)
   trias_add = newSF_tot{select_r_num_int2(huk)}.trias;
   m_aBlocks{fk}.trias =  [m_aBlocks{fk}.trias; trias_add]; 
end



end