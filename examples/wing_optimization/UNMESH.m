classdef UNMESH

    properties
        orgMesh
        orgSurf
        designVariables
    end

    methods(Access = public)
        function obj = UNMESH(orgVerts,orgCon,orgSurfVerts,orgSurfCon)

        end

        function [modVerts,Con] = meshDeformation(obj,modSurfVerts,modSurfCon)

        end
        
        function [modVerts,Con] = makeMeshGradient
        end

    end

end