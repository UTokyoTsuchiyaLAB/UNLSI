%combine surfaces together
%used in the boolean code for triangulated surfaces

function FV_out = combine_surfaces(FV_A,FV_B)
%Input
% FV_A : surface structure with feilds vertices and faces to be combined with FV_B
% FV_B : surface structure with feilds vertices and faces to be combined with FV_A
%Output
% FV_out : new surface with unified nodes

    tolerance = 10^-10;
    %find out if the edge field exists
    fnA = fieldnames(FV_A);
    fnB = fieldnames(FV_B);
    
    trias_FV_tmpA = FV_A.faces;
    nodes_master_new_FV_tmpA = FV_A.vertices;
    if any(contains(upper(fnA),'EDGES'))
        new_edges__tmpA = FV_A.edges;
    end
    
    trias_FV_tmpB = FV_B.faces;
    nodes_master_new_FV_tmpB = FV_B.vertices;
    
    if any(contains(upper(fnB),'EDGES'))
        new_edges__tmpB = FV_B.edges;
    end
    
    %combine both meshes together so that they share the same nodal list
    nodes_master_new_FV_global = uniquetol([nodes_master_new_FV_tmpA; nodes_master_new_FV_tmpB],tolerance,'Byrows',true);
    
    FV_out.vertices = nodes_master_new_FV_global;
    
    %update the triangle connectivity tables for surface A
    if ~isempty(trias_FV_tmpA)
        [~,locins1] = ismembertol(nodes_master_new_FV_tmpA(trias_FV_tmpA(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
        trias_FV_tmpA(:) = locins1;
    end
    %for surface B
    if ~isempty(trias_FV_tmpB)
        [~,locins2] = ismembertol(nodes_master_new_FV_tmpB(trias_FV_tmpB(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
        trias_FV_tmpB(:) = locins2;
    end
    
    FV_out.faces = unique([trias_FV_tmpA; trias_FV_tmpB],'rows');
        
    if any(contains(upper(fnA),'EDGES')) || any(contains(upper(fnB),'EDGES'))
        if any(contains(upper(fnA),'EDGES')) && any(contains(upper(fnB),'EDGES'))
            %update the new edges for surface A
            [~,locins3] =ismembertol(nodes_master_new_FV_tmpA(new_edges__tmpA(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
            new_edges__tmpA(:) = locins3;
            [~,locins4] =ismembertol(nodes_master_new_FV_tmpB(new_edges__tmpB(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
            new_edges__tmpB(:) = locins4;
        elseif any(contains(upper(fnA),'EDGES')) && ~any(contains(upper(fnB),'EDGES'))
            [~,locins3] =ismembertol(nodes_master_new_FV_tmpA(new_edges__tmpA(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
            new_edges__tmpA(:) = locins3;
            new_edges__tmpB = [];
        elseif ~any(contains(upper(fnA),'EDGES')) && any(contains(upper(fnB),'EDGES'))
            new_edges__tmpA =[];
            [~,locins4] =ismembertol(nodes_master_new_FV_tmpB(new_edges__tmpB(:),:),nodes_master_new_FV_global,tolerance,'Byrows',true);
            new_edges__tmpB(:) = locins4;
        end
    
        FV_out.edges = unique([new_edges__tmpA; new_edges__tmpB],'rows');
    
    end

end