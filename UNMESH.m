classdef UNMESH

    properties
        orgMesh
        orgSurf
        designVariables
        lb
        ub
        designScale
        gradSurf
    end

    methods(Access = public)
        function obj = UNMESH(orgVerts,orgCon,orgSurfVerts,orgSurfCon,designVariables,lb,ub)
            obj.orgMesh = triangulation(orgCon,orgVerts);
            obj.orgSurf =  triangulation(orgSurfCon,orgSurfVerts);
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            if nargin == 4
                designVariables = [];
            end
            obj.designVariables = (designVariables(:)'-obj.lb)./obj.designScale;
            
        end

        function [modVerts,con] = meshDeformation(obj,modSurfVerts,modSurfCon)
            if nargin == 3
                %connectivityが一致しているか確認
                if any(obj.orgSurf.ConnectivityList(:) == modSurfCon(:))
                    error("Surf connectivity is not match")
                end
            end

            md.x = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,1)-obj.orgSurf.Points(:,1),'linear','linear');
            md.y = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,2)-obj.orgSurf.Points(:,2),'linear','linear');
            md.z = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,3)-obj.orgSurf.Points(:,3),'linear','linear');
            dVerts(:,1) = md.x(obj.orgMesh.Points);
            dVerts(:,2) = md.y(obj.orgMesh.Points);
            dVerts(:,3) = md.z(obj.orgMesh.Points);
            modVerts = obj.orgMesh.Points+dVerts;
            con = obj.orgSurf.ConnectivityList;
        end
        
        function obj = makeMeshGradient(obj,surfGenFun)
            %設計変数勾配による表面近似
            pert = 0.001./obj.designScale;
            ndim = numel(obj.designVariables);
            for i = 1:ndim
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) + pert(i)).*obj.designScale(i)+obj.lb(i);
                modSurf = surfGenFun(sampleDes);
                dmodSurf = modSurf-obj.orgSurf.Points;
                sampleSurff = dmodSurf(:);
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                modSurf = surfGenFun(sampleDes);
                dmodSurf = modSurf-obj.orgSurf.Points;
                sampleSurfr = dmodSurf(:);
                obj.gradSurf(:,i) = (sampleSurff-sampleSurfr)./2./pert(i);
            end
        end

        function modSurf = makeSurffromVariables(obj,x)
            %スケーリングされてない
            scaledVar = (x(:)'- obj.lb) ./ obj.designScale; 
            modSurf = obj.orgSurf.Points + reshape(obj.gradSurf*(scaledVar(:)-obj.designVariables(:)),size(obj.orgSurf.Points));
        end

        function [obj0, dobj_dx, con0, dcons_dx] = calcObjandConsGradients(obj,objandConsFun)
            
            [obj0,con0] = objandConsFun(obj.designVariables.*obj.designScale+obj.lb);
            dobj_dx = zeros(1,numel(obj.designVariables));
            if not(isempty(con0))
                dcons_dx = zeros(numel(cons),numel(obj.designVariables));
            else
                dcons_dx = [];
            end
            pert = sqrt(eps);
            for i = 1:numel(obj.designVariables)
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) + pert).*obj.designScale(i)+obj.lb(i);
                [objf,conf] = objandConsFun(sampleDes);
                dobj_dx(i) = (objf-obj0)/pert;
                if not(isempty(con0))
                    dcons_dx(:,i) = (conf-con0)/pert;
                end
            end
        end
    end

end