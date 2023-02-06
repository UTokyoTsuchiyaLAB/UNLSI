classdef UNLSI

    properties
        tri
        surfID
        wakeline
        paneltype
        flow %各flowconditionが格納されたセル配列
        cluster
    end

    methods(Access = public)
        function obj = UNLSI(verts,connectivity,surfID,wakeline)
            %Constructor
            obj.tri = triangulation(connectivity,verts);
            obj.surfID = surfID;
            obj.wakeline = wakeline;
            obj.paneltype = ones(size(connectivity,1),1);
            obj.flow = {};
            obj.cluster = cell([1,numel(obj.paneltype)]);
        end

        function obj = setVerts(obj,verts)
            obj.tri.Points = verts;
        end

        function obj = setPanelType(obj,ID,typename)
            %paneltypeの指定関数
            %panelname : body 1 機体表面パネル
            %panelname : base 2 ベース面パネル
            %panelname : structure 3 構造パネル 
            if strcmpi(typename,"body")
                obj.paneltype(obj.surfID==ID,1) = 1;
            elseif strcmpi(typename,"base")
                obj.paneltype(obj.surfID==ID,1) = 2;
            elseif strcmpi(typename,"structure")
                obj.paneltype(obj.surfID==ID,1) = 3;    
            end
        end

        function obj = makeCluster(obj,nCluster)
           for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1 %表面パネルであれば
                    obj.cluster{i} = i;
                    index = 1;
                    while(1)
                        neighbor = obj.tri.neighbors(obj.cluster{i}(index));
                        diffcluster = setdiff(neighbor,obj.cluster{i});
                        if isempty(diffcluster) || numel(obj.cluster{i})>nCluster
                            break;
                        end
                        obj.cluster{i} = [obj.cluster{i},diffcluster];
                        index = index + 1;
                    end
                end
           end
        end
        
        function obj = addFlowCondition(obj,Mach,newtoniantype)
            if nargin < 3
                newtoniantype = "OldTangentCone";
            end
            nFlow = numel(obj.flow);
            obj.flow{nFlow+1}.Mach = Mach;
            obj.flow{nFlow+1}.newtoniantype = newtoniantype;
            kappa = 1.4;

            if obj.flow{nFlow+1}.Mach > 1
                delta = linspace(-pi,pi,500);
                Cp = zeros(size(delta));
                for j =1:size(delta,2)
                    if delta(j) >= 0
                        if strcmpi(newtoniantype,'TangentConeEdwards')
                            Mas = (0.87*Mach-0.554)*sin(delta(j))+0.53;
                            Cp(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(newtoniantype,'OldTangentCone')
                            Mas = 1.090909*Mach*sin(delta(j))+exp(-1.090909*Mach*sin(delta(j)));
                            Cp(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(newtoniantype,'TangentWedge')
                            if delta(j)>45.585*pi/180
                                Cp(j)=((1.2*Mach*sin(delta(j))+exp(-0.6*Mach*sin(delta(j))))^2-1.0)/(0.6*Mach^2);
                                %R=1/Mach^2+(kappa+1)/2*delta(j)/sqrt(Mach^2-1);
                            elseif delta(j)<0.035
                                Cp(j)=kappa*Mach^2*delta(j)/sqrt(Mach^4-1.0);
                            else
                                b = -((Mach^2+2)/Mach^2)-kappa*sin(delta(j))^2;
                                c = (2*Mach^2+1)/Mach^4+((kappa+1)^2/4+(kappa-1)/Mach^2)*sin(delta(j))^2;
                                d = -cos(delta(j))^2/Mach^4;
                                q = (b^2-3*c)/9;
                                rr = (b*(2*b^2-9*c)+27*d)/54;
                                disc = q^3-rr^2;
                                if disc<0
                                    Cp(j)=((1.2*Mach*sin(delta(j))+exp(-0.6*Mach*sin(delta(j))))^2-1.0)/(0.6*Mach^2);
                                else
                                    r = roots([1,b,c,d]);
                                    %ts = asin(sqrt(r(2)));
                                    R=r(2);
                                    Cp(j) = 4*(Mach^2*R-1)/((kappa+1)*Mach^2);
                                end
                            end
                        end
                        Pe_Pinf= Cp(j)*(kappa*Mach^2)/2+1;
                        Te_Tinf = Pe_Pinf^(1/(kappa/(kappa-1)));
                        Me = sqrt((2+(kappa+1)*Mach^2)/((kappa+1)*Te_Tinf)-2/(kappa+1));
                    else
                        Me = UNLSI.bisection(@(M2)UNLSI.pmsolve(Mach,M2,-delta(j),kappa),Mach,300);
                        Te_Tinf = (2+(kappa+1)*Mach^2)/(2+(kappa+1)*Me^2);
                        Pe_Pinf=(Te_Tinf)^(kappa/(kappa-1));
                        Cp(j) = 2/(kappa*Mach^2)*(Pe_Pinf-1);
                    end
                end
                obj.flow{nFlow+1}.pp = griddedInterpolant(delta,Cp,'spline');
            end
        end
    end

    methods(Static)
        function res = pmsolve(M1,M2,nu,kappa)
            M1nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M1^2-1)))-atan(sqrt(M1^2-1));
            M2nu = M1nu+nu;
            res = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M2^2-1)))-atan(sqrt(M2^2-1))-M2nu;
        end

        function [p info] = bisection(fun,a,b)
            % provide the equation you want to solve with R.H.S = 0 form. 
            % Write the L.H.S by using inline function
            % Give initial guesses.
            % Solves it by method of bisection.
            % A very simple code. But may come handy
            i = 1;
            if fun(a)*fun(b)>0 
                p = (a + b)/2;
                info = -1;
            else
                p = (a + b)/2;
                err = abs(fun(p));
                while(err > sqrt(eps)) 
                   if fun(a)*fun(p)<0 
                       b = p;
                   else
                       a = p;          
                   end
                    p = (a + b)/2; 
                   err = abs(fun(p));
                   i = i+1;
                end
                info = 1;
            end
        end
    end

end