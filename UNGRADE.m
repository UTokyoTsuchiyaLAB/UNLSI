classdef UNGRADE < UNLSI

    properties
        orgMesh
        orgSurf
        designVariables
        surfGenFun
        meshGenFun
        lb
        ub
        unlsiParam
        designScale
        gradSurf
        approxMat %ソルバーの近似行列
        hessianUpdate
    end

    methods(Access = public)

        function obj = UNGRADE(surfGenFun,designVariables,lb,ub,orgVerts,orgCon,surfID,wakelineID,halfmesh)
            %%%%%%%%%%%%メッシュの登録と基準サーフェス生成関数の登録%%%%%%%%%%%%%%%
            %orgVerts,orgCon　実際に解析を行うメッシュ（openVSPのCFDツール等で生成した（基本的に）オーバーラップのない非構造メッシュ）
            %meshGenFun　設計変数（スケールされていない）から基準サーフェス（解析メッシュと同じ表面を持つ、トリムされていないメッシュ）
            %designVariables : 設計変数
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj = obj@UNLSI(orgVerts,orgCon,surfID,wakelineID,halfmesh)
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            obj.orgMesh = triangulation(orgCon,orgVerts);
            obj.surfGenFun = surfGenFun;
            [orgSurfVerts,orgSurfCon,desOrg] = obj.surfGenFun(designVariables(:)');
            obj.designVariables = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgSurf =  triangulation(orgSurfCon,orgSurfVerts);
        end

        function obj = setMeshGenFun(obj,meshGenFun)
            obj.meshGenFun = meshGenFun;
        end

        function obj = modifyMesh(obj,designVariables,orgVerts,orgCon,surfID,wakelineID)
            obj = obj.setMesh(orgVerts,orgCon,surfID,wakelineID);
            obj.orgMesh = triangulation(orgCon,orgVerts);
            [orgSurfVerts,orgSurfCon,desOrg] = obj.surfGenFun(designVariables(:)');
            obj.designVariables = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgSurf =  triangulation(orgSurfCon,orgSurfVerts);
        end

        function obj = modifyMeshfromVariables(obj,designVariables)
            [orgVerts,orgCon,surfID,wakelineID,designVariables] = obj.meshGenFun(designVariables);
            obj = obj.modifyMesh(designVariables,orgVerts,orgCon,surfID,wakelineID);
        end

        function [modVerts,con] = meshDeformation(obj,modSurfVerts,modSurfCon)
            %%%%%%%%%%%%非構造メッシュのメッシュ変形%%%%%%%
            %補間ベースのメッシュ変形
            %modSurfVers : 変形後の基準サーフェスの接点情報
            %modSurfCon : 入力すると一応conが一致しているかチェックする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 3
                %connectivityが一致しているか確認
                if any(obj.orgSurf.ConnectivityList(:) ~= modSurfCon(:))
                    error("Surf connectivity is not match")
                end
            end
            md.x = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,1)-obj.orgSurf.Points(:,1),'natural','linear');
            md.y = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,2)-obj.orgSurf.Points(:,2),'natural','linear');
            md.z = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,3)-obj.orgSurf.Points(:,3),'natural','linear');
            dVerts(:,1) = md.x(obj.orgMesh.Points);
            dVerts(:,2) = md.y(obj.orgMesh.Points);
            dVerts(:,3) = md.z(obj.orgMesh.Points);
            modVerts = obj.orgMesh.Points+dVerts;
            con = obj.orgSurf.ConnectivityList;
        end
        
        function obj = makeMeshGradient(obj)
            %%%%%%%%%%%設計変数勾配による標本表面の近似関数の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pert = 0.01./obj.designScale;
            ndim = numel(obj.designVariables);
            desOrg = obj.designVariables.*obj.designScale+obj.lb;
            [surforg,~,desOrg] = obj.surfGenFun(desOrg);
            for i = 1:ndim
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) + pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,desBuff] = obj.surfGenFun(sampleDes);
                pertf = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurff = dmodSurf(:);
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,desBuff] = obj.surfGenFun(sampleDes);
                pertr = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurfr = dmodSurf(:);
                obj.gradSurf(:,i) = (sampleSurff-sampleSurfr)./(pertf-pertr);
            end
        end

        function [modSurf,modMesh] = variables2Mesh(obj,x,method)
           %%%%%%%%%%%%設計変数を変更した場合の標本近似表面の節点値の計算


           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if nargin == 2
               method = "raw";
           end
           if strcmpi(method,"raw")
               modSurf = obj.surfGenFun(x(:)');
           elseif strcmpi(method,"linear")
               scaledVar = (x(:)'- obj.lb) ./ obj.designScale; 
               modSurf = obj.orgSurf.Points + reshape(obj.gradSurf*(scaledVar(:)-obj.designVariables(:)),size(obj.orgSurf.Points));
           end
           if nargout == 2
               modMesh = obj.meshDeformation(modSurf);
           end
        end
        
        function obj = calcApproximatedEquation(obj)
            if obj.approximated == 1
                error("This instance is approximated. Please execute obj.makeEquation()");
            end
            nbPanel = sum(obj.paneltype == 1);
            %接点に繋がるIDを決定
            vertAttach = obj.tri.vertexAttachments();
            pert = sqrt(eps);
            for i = 1:numel(vertAttach)
                if mod(i,floor(numel(vertAttach)/10))==0 || i == 1
                    fprintf("%d/%d ",i,numel(vertAttach));
                end
                %paneltype == 1以外を削除する
                vertAttach{i}(obj.paneltype(vertAttach{i}) ~=1) = [];
                obj.approxMat.calcIndex{i} = sort(obj.IndexPanel2Solver(vertAttach{i}));
                %このvertsがwakeに含まれているか
                for j = 1:3
                    newVerts = obj.tri.Points;
                    newVerts(i,j) = obj.tri.Points(i,j)+pert;
                    obj2 = obj.setVerts(newVerts);
                    [VortexAr,VortexBr,VortexAc,VortexBc] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});
                    %
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakePanelLength,obj.nWake);
                                VortexAr(:,interpID(1)) = VortexAr(:,interpID(1)) - influence;
                                VortexAr(:,interpID(2)) = VortexAr(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakePanelLength,obj.nWake);
                                VortexAr(:,interpID(1)) = VortexAr(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexAr(:,interpID(2)) = VortexAr(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAc(:,b) = VortexAc(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAc(:,b) = VortexAc(:,b) + influence;
                                end
                            end
                        end
                    end
                    obj.approxMat.dVAr{i,j} = (VortexAr-obj.LHS(obj.approxMat.calcIndex{i},:))./pert;
                    obj.approxMat.dVBr{i,j} = (VortexBr-obj.RHS(obj.approxMat.calcIndex{i},:))./pert;
                    obj.approxMat.dVAc{i,j} = (VortexAc-obj.LHS(:,obj.approxMat.calcIndex{i}))./pert;
                    obj.approxMat.dVBc{i,j} = (VortexBc-obj.RHS(:,obj.approxMat.calcIndex{i}))./pert;
                end
            end
            fprintf("\n");
        end

        function approxmatedObj = makeAproximatedInstance(obj,modifiedVerts)
                approxmatedObj = obj.setVerts(modifiedVerts);
                approxmatedObj.approxMat = [];
                approxmatedObj.approximated = 1;
                for i = 1:size(modifiedVerts,1)
                    for j = 1:3
                        dv = (approxmatedObj.tri.Points(i,j)-obj.tri.Points(i,j));
                        obj.approxMat.dVAc{i,j}(obj.approxMat.calcIndex{i},:) = 0; 
                        obj.approxMat.dVBc{i,j}(obj.approxMat.calcIndex{i},:) = 0; 
                        approxmatedObj.LHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.LHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVAr{i,j}.*dv;
                        approxmatedObj.RHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.RHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVBr{i,j}.*dv;
                        approxmatedObj.LHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.LHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAc{i,j}.*dv;
                        approxmatedObj.RHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.RHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVBc{i,j}.*dv;
                    end
                end
        end

        function viewMesh(obj,modMesh,fig)
            %%%%%%%%%%%設計変数勾配による標本表面の近似関数の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            viewtri = triangulation(obj.orgMesh.ConnectivityList,modMesh);
            figure(fig);
            trisurf(viewtri, 'FaceAlpha', 0, 'EdgeColor', 'black');
            axis equal;drawnow();
        end

        function viewSurf(obj,modSurf,fig)
            %%%%%%%%%%%設計変数勾配による標本表面の近似関数の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            viewtri = triangulation(obj.orgSurf.ConnectivityList,modSurf);
            figure(fig);
            trisurf(viewtri, 'FaceAlpha', 0, 'EdgeColor', 'black');
            axis equal;drawnow();
        end

        function checkSurfGenWork(obj,pert,fig)
            ndim = numel(obj.designVariables);
            randparam = rand(1).*ones(1,ndim);
            randDes = randparam.*obj.designScale+obj.lb;
            modSurf = obj.surfGenFun(randDes);
            viewtri = triangulation(obj.orgSurf.ConnectivityList,modSurf);
            figure(fig);clf
            trisurf(viewtri);
            axis equal;drawnow();
            fprintf("randamized parameter : ");
            disp(randDes);
            pause(1);
            for i = 1:ndim
                sampleDes = randparam.*obj.designScale+obj.lb;
                sampleDes(i) = (randparam(i) + pert).*obj.designScale(i)+obj.lb(i);
                modSurf = obj.surfGenFun(sampleDes);
                viewtri = triangulation(obj.orgSurf.ConnectivityList,modSurf);
                figure(fig);clf
                trisurf(viewtri);
                axis equal;drawnow();
                pause(1)
            end
        end

        function obj = setOptCondition(obj,Mach,alpha,beta,wakeLength,n_wake,n_divide,nCluster,edgeAngleThreshold,Re,Lch,k,LTratio,CfeCoefficient)
            obj = obj.flowCondition(1,Mach);
            obj.unlsiParam.alpha = alpha;
            obj.unlsiParam.beta = beta;        
            obj.unlsiParam.n_wake = n_wake;
            obj.unlsiParam.wakeLength = wakeLength;
            obj.unlsiParam.n_divide = n_divide;
            obj.unlsiParam.nCluster = nCluster;
            obj.unlsiParam.edgeAngleThreshold = edgeAngleThreshold;
            obj.unlsiParam.Re = Re;
            obj.unlsiParam.Lch = Lch;
            obj.unlsiParam.k = k;
            obj.unlsiParam.LTratio = LTratio;
            obj.unlsiParam.coefficient = CfeCoefficient;

        end

        function obj = setHessianUpdate(obj,H0,TR,method,nMemory)
            obj.hessianUpdate.H0 = H0;
            obj.hessianUpdate.H = H0;
            obj.hessianUpdate.nMemory = nMemory;
            obj.hessianUpdate.TR = TR;
            obj.hessianUpdate.xScaled = [];
            obj.hessianUpdate.dL_dx = [];
            if strcmpi(method,"SR1")
                obj.hessianUpdate.updateFunction = @obj.SR1;
            elseif strcmpi(method,"SSR1")
                obj.hessianUpdate.updateFunction = @obj.SSR1;
            elseif strcmpi(method,"BFGS")
                obj.hessianUpdate.updateFunction = @obj.BFGS;
            elseif strcmpi(method,"MBFGS")
                obj.hessianUpdate.updateFunction = @obj.MBFGS;
            elseif strcmpi(method,"DBFGS")
                obj.hessianUpdate.updateFunction = @obj.DBFGS;
            elseif strcmpi(method,"SR1_BFGS")
                obj.hessianUpdate.updateFunction = @(s,y,H)obj.SR1_BFGS(obj,s,y,H);
            elseif strcmpi(method,"SSR1_MBFGS")
                obj.hessianUpdate.updateFunction = @(s,y,H)obj.SSR1_MBFGS(obj,s,y,H);
            end
        end

        function [dx,obj] = finddx(obj,objandConsFun,cmin,cmax)
            %%%%%%%%%%%%指定した評価関数と制約条件における設計変数勾配の計算%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %初期点の解析
            obj = obj.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
            desOrg = obj.designVariables.*obj.designScale+obj.lb;
            [u0,~] = obj.solvePertPotential(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
            [AERODATA0,Cp0,Cfe0,R0,obj] = obj.solveFlowForAdjoint(u0,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
            [I0,con0] = objandConsFun(desOrg,AERODATA0,Cp0,Cfe0);
            disp([I0,con0(:)']);
            %u微分の計算
            pert = sqrt(eps);
            for i = 1:numel(u0)
                u = u0;
                u(i) = u(i)+pert;
                [AERODATA,Cp,Cfe] = obj.solveFlowForAdjoint(u,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                [I,con] = objandConsFun(desOrg,AERODATA,Cp,Cfe);
                dI_du(i) = (I-I0)/pert;
                if not(isempty(con))
                    dcon_du(:,i) = (con-con0)/pert;
                end
            end
            dR_du = obj.LHS;
            %x微分の計算
            %メッシュの節点勾配を作成
            obj = obj.makeMeshGradient();

            for i= 1:numel(obj.designVariables)
                x = obj.designVariables;
                x(i) = obj.designVariables(i)+pert;
                des = x.*obj.designScale+obj.lb;
                [~,modMesh] = obj.variables2Mesh(des,'linear');
                %変数が少ないときは直接作成
                obj2 = obj.setVerts(modMesh);
                obj2 = obj2.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
                [AERODATA,Cp,Cfe,R] = obj2.solveFlowForAdjoint(u0,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                [I,con] = objandConsFun(des,AERODATA,Cp,Cfe);
                dI_dx(i) = (I-I0)/pert;
                if not(isempty(con))
                    dcon_dx(:,i) = (con-con0)/pert;
                end
                dR_dx(:,i) = (R-R0)./pert;
            end
            %全微分でdxからduを求める行列
            duMatVec = -(dR_du)\([R0,dR_dx]);
            duMat = duMatVec(:,2:end);
            duVec = duMatVec(:,1);
            
            objTotalGrad = dI_dx+dI_du*duMat;
            if not(isempty(con))
                conTotalGrad = dcon_dx+dcon_du*duMat;
            end
            %線形計画問題に変換
            if not(isempty(con))
                alin = [-conTotalGrad;conTotalGrad];
                blin = [con0(:)-cmin(:);cmax(:)-con0(:)];
            else
                alin = [];
                blin = [];
            end
            lbf = -obj.designVariables;
            ubf = 1-obj.designVariables;
            %[dxscaled,fval,exitflag,output,lambda] = quadprog(obj.hessianUpdate.H,objTotalGrad,alin,blin,[],[],lbf,ubf);
            options = optimoptions(@fmincon,'Algorithm','sqp');
            [dxscaled,fval,exitflag,output,lambda] = fmincon(@(dx)obj.fminconObj(dx,obj.hessianUpdate.H,objTotalGrad),zeros(numel(obj.designVariables),1),alin,blin,[],[],lbf,ubf,@(dx)obj.fminconNlc(dx,obj.hessianUpdate.TR),options);
            if not(isempty(con))
                lambdaR = -(dR_du)\(dI_du+lambda.ineqlin'*[-dcon_du;dcon_du])';
                dL_dx = dI_dx+lambdaR'*dR_dx+lambda.ineqlin'*[-dcon_dx;dcon_dx];
            else
                lambdaR = -(dR_du)\(dI_du)';

                dL_dx = dI_dx+lambdaR'*dR_dx;
            end
            dx = dxscaled(:)'.*obj.designScale;
            
            %HessianUpdate
            obj.hessianUpdate.xScaled = [obj.hessianUpdate.xScaled;obj.designVariables];
            obj.hessianUpdate.dL_dx = [obj.hessianUpdate.dL_dx;dL_dx];
            if size(obj.hessianUpdate.xScaled,1) > 1
                n_iter = size(obj.hessianUpdate.xScaled,1)-1;
                if n_iter > obj.hessianUpdate.nMemory
                    n_iter = obj.hessianUpdate.nMemory;
                end
                obj.hessianUpdate.H = obj.hessianUpdate.H0;
                for i = 1:n_iter
                    s = obj.hessianUpdate.xScaled(end-(i-1),:)-obj.hessianUpdate.xScaled(end-i,:);
                    y = obj.hessianUpdate.dL_dx(end-(i-1),:)-obj.hessianUpdate.dL_dx(end-i,:);
                    obj.hessianUpdate.H = obj.hessianUpdate.updateFunction(s,y,obj.hessianUpdate.H);
                end
            end
        end

        function obj = updateVariables(obj,FcnObjandCon,cmin,cmax)
            [dx,obj] = obj.finddx(FcnObjandCon,cmin,cmax);
            obj.plotGeometry(1,obj.Cp,[-2,1]);
            newDes = obj.designVariables.*obj.designScale+obj.lb+dx;
            disp(newDes);
            obj = obj.modifyMeshfromVariables(newDes);
            obj = obj.setCf(1,obj.unlsiParam.Re,obj.unlsiParam.Lch,obj.unlsiParam.k,obj.unlsiParam.LTratio,obj.unlsiParam.coefficient);
            obj = obj.makeCluster(obj.unlsiParam.nCluster,obj.unlsiParam.edgeAngleThreshold);
        end
    end

    methods(Static)
        function res = fminconObj(dx,H,g)
            res = 0.5 * dx(:)'*H*dx(:) + g*dx(:);
        end

        function [c,ceq] = fminconNlc(dx,TR)
            ceq = [];
            c = sum(dx.^2)-(TR)^2;
        end

        function Bkp1 = SR1_BFGS(obj,s,y,Bk)
            b = (s(:)'*Bk*s(:))/(s(:)'*y(:));
            eta2 = 10^-7;
            if b < 1-eta2 && not(isinf(b))
                fprintf('Hessian update mode is SR1/BFGS. This time is SR1\n')
                Bkp1 = obj.SR1(s,y,Bk);
            else
                fprintf('Hessian update mode is SR1/BFGS. This time is BFGS\n')
                Bkp1 = obj.BFGS(s,y,Bk);
            end
        end
        function Bkp1 = SSR1_MBFGS(obj,s,y,Bk)
            b = (s(:)'*Bk*s(:))/(s(:)'*y(:));
            eta2 = 10^-7;
            if b < 1-eta2 && not(isinf(b))
                fprintf('Hessian update mode is SSR1/MBFGS. This time is SSR1\n')
                Bkp1 = obj.SSR1(s,y,Bk);
            else
                fprintf('Hessian update mode is SSR1/MBFGS. This time is MBFGS\n')
                Bkp1 = obj.MBFGS(s,y,Bk);
            end
        end

        function Bkp1 = BFGS(s,y,Bk)
            %%%%%%%%%%%%%%%%%BFGS%%%%
            %準ニュートン法のBFGSアップデートの手実装
            %s 設計変数の変化量
            %y ヤコビアン
            %%%%%%%%%%%%%%%%%%%%%%%%
            s = s(:);
            y=y(:);
            if s'*y<=0
                Bkp1 = Bk;
            else
                Bkp1 = Bk-(Bk*s*(Bk*s)')/(s.'*Bk*s)+(y*y.')/(s.'*y);
                coeB=(s'*y)/(s'*Bkp1*s);
                if coeB>0 && coeB<1
                    Bkp1=coeB.*Bkp1;
                end
            end
        end
        function Bkp1 = MBFGS(s,y,Bk)
            s = s(:);
            y=y(:);
            if s'*(Bk*s-y)==0
                Bkp1 = Bk;
            else
                if s'*y>=0.2*(s'*Bk*s)
                    phi = 1;
                else
                    phi = 0.8*(s'*Bk*s)/(s'*(Bk*s-y));
                end
                yhat = phi.*y+(1-phi).*(Bk*s);
                Bkp1 = Bk-(Bk*s*(Bk*s)')/(s.'*Bk*s)+(yhat*yhat.')/(s.'*yhat);
                coeB=(s'*y)/(s'*Bkp1*s);
                if coeB>0 && coeB<1
                    Bkp1=coeB.*Bkp1;
                end
            end
        end

        function Bkp1 = DBFGS(s,y,Bk)
            s = s(:);
            y=y(:);
            if y'*s > 0
                pk = y'*s/(s'*Bk*s);
                bk = 1./pk;
                hk = y'*s/(y'*(Bk\y));
                ak = bk*hk-1;
                if hk<1
                    theta = 1/(1-bk);
                else
                    theta = 0;
                end
                sigma2 = max(1-1/ak,0.5);
                sigma3 = Inf;
                if pk<1-sigma2
                    phi = sigma2/(1-pk);
                elseif pk>1+sigma3
                    phi = sigma3/(pk-1);
                else
                    phi = 1;
                end
                yhat = phi.*y+(1-phi).*(Bk*s);
                w = sqrt(s'*Bk*s).*(yhat./(y'*s)-Bk*s/(s'*Bk*s));
                Bkp1 = Bk-(Bk*s*(Bk*s)')/(s.'*Bk*s)+(yhat*yhat.')/(s.'*yhat)+theta*(w*w');
            else
                Bkp1 = Bk;
            end
    
        end

        function Bkp1 = SR1(s,y,Bk)
            s = s(:);
            y=y(:);
            if abs(s'*(y-Bk*s)) < 10^-8*norm(s)*norm(y-Bk*s) || ((y-Bk*s)'*s)==0
                Bkp1 = Bk;
            else
                Bkp1 = Bk+(y-Bk*s)*(y-Bk*s)'/((y-Bk*s)'*s);
                coeB =(s'*y)/(s'*Bkp1*s);
                if coeB>0 && coeB<1
                    Bkp1=coeB.*Bkp1;
                end
            end
        end

        function Bkp1 = SSR1(s,y,Bk)
            s = s(:);
            y=y(:);
            if abs(s'*(y-Bk*s)) < 10^-8*norm(s)*norm(y-Bk*s) || ((y-Bk*s)'*s)==0 || y'*s-s'*Bk*s<=0
                lambda_k = y'*y/(y'*s)-sqrt((y'*y)^2/(y'*s)^2-(y'*y)/(s'*s));
                if isnan(lambda_k)||not(isreal(lambda_k))
                    Bkp1 = Bk;
                else
                    Bk = 1./lambda_k.*eye(size(Bk,1));
                    Bkp1 = Bk+(y-Bk*s)*(y-Bk*s)'/((y-Bk*s)'*s);
                    coeB =(s'*y)/(s'*Bkp1*s);
                    if coeB>0 && coeB<1
                        Bkp1=coeB.*Bkp1;
                    end
                end
            else
                Bkp1 = Bk+(y-Bk*s)*(y-Bk*s)'/((y-Bk*s)'*s);
                coeB =(s'*y)/(s'*Bkp1*s);
                if coeB>0 && coeB<1
                    Bkp1=coeB.*Bkp1;
                end
            end
        end

    end

end