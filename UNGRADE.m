 classdef UNGRADE < UNLSI

    properties
        orgMesh
        orgGeom
        unscaledVar
        scaledVar
        geomGenFun
        meshGenFun
        lb
        ub
        setting
        designScale
        gradSurf
        approxMat
        approximated        
        gradSREF
        gradBREF
        gradCREF
        gradXYZREF
        gradArginx
        argin_x
        iteration
        Hessian
        history
        flowNoList
        stabbbScaling
    end

    methods(Access = public)

        function obj = UNGRADE(meshGenFun,geomGenFun,unscaledVariables,lb,ub,halfmesh,Mach,alpha,beta)
            %%%%%%%%%%%%メッシュの登録と基準サーフェス生成関数の登録%%%%%%%%%%%%%%%
            %orgVerts,orgCon　実際に解析を行うメッシュ（openVSPのCFDツール等で生成した（基本的に）オーバーラップのない非構造メッシュ）
            %meshGenFun　設計変数（スケールされていない）から基準サーフェス（解析メッシュと同じ表面を持つ、トリムされていないメッシュ）
            %designVariables : 設計変数
            %lb：設計変数の下限
            %ub：設計変数の上限
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [orgMeshVerts, orgMeshCon,surfID,wakeLineID, desOrg] = meshGenFun(unscaledVariables(:)');
            obj = obj@UNLSI(orgMeshVerts,orgMeshCon,surfID,wakeLineID,halfmesh)
            
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            if any(obj.designScale<0)
                error("check the ub and lb setting");
            end
            obj.geomGenFun = geomGenFun;
            obj.meshGenFun = meshGenFun;
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(desOrg);
            obj.orgMesh = triangulation(orgMeshCon,orgMeshVerts);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
            obj.unscaledVar = desOrg(:)';
            obj.scaledVar = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);

            obj.setting.Mach = Mach;
            obj.setting.alpha = alpha;
            obj.setting.beta = beta;
            if numel(obj.setting.Mach) ~= numel(obj.setting.alpha) || numel(obj.setting.Mach) ~= numel(obj.setting.beta) || numel(obj.setting.beta) ~= numel(obj.setting.alpha)
                error("No. of case is not match");
            end
            obj.iteration = 0;
            
            %%%%%%%%%%%オプションのデフォルト設定
            obj.setting.nMemory = 20;
            obj.setting.H0 = eye(numel(obj.lb));
            obj.setting.n_wake = 5;
            obj.setting.wakeLength = 20;
            obj.setting.n_divide = 10;
            obj.setting.nCluster = 50;
            obj.setting.edgeAngleThreshold = 50;
            obj.setting.Re = 500000;
            obj.setting.Lch = 1;
            obj.setting.k = 0.052*(10^-5);
            obj.setting.LTratio = 0;
            obj.setting.coefficient = 1;
            obj.setting.gradientCalcMethod = "direct";
            obj.setting.HessianUpdateMethod = "SR1";
            obj.setting.betaLM = 0.5;
            obj.setting.TrustRegion = 0.1;
            %%%%%%%%%%%%%

            obj.Hessian = obj.setting.H0;

            uniMach = unique(Mach);
            for i = 1:numel(uniMach)
                obj = obj.flowCondition(i,uniMach(i));
                obj = obj.setCf(i,obj.setting.Re,obj.setting.Lch,obj.setting.k,obj.setting.LTratio,obj.setting.coefficient);
            end
            for i = 1:numel(Mach)
                obj.flowNoList(i,1) = find(uniMach == Mach(i));
                obj.flowNoList(i,2) = uniMach(obj.flowNoList(i,1));
                obj.flowNoList(i,3) = obj.flowNoList(i,2)<1;
            end
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeCluster(obj.setting.nCluster,obj.setting.edgeAngleThreshold);
            end

        end
        
        function obj = updateMeshGeomfromVariables(obj,unscaledVariables)
            %%%%%%%設計変数からメッシュや解析条件を更新する%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [orgMeshVerts, orgMeshCon,surfID,wakeLineID, unscaledVariables] = obj.meshGenFun(unscaledVariables(:)');
            obj = obj.setMesh(orgMeshVerts,orgMeshCon,surfID,wakeLineID);
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
            obj.orgMesh = triangulation(orgMeshCon,orgMeshVerts);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
            obj.unscaledVar = desOrg(:)';
            obj.scaledVar = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeCluster(obj.setting.nCluster,obj.setting.edgeAngleThreshold);
            end
            for i = 1:numel(obj.flow)
                obj = obj.setCf(i,obj.setting.Re,obj.setting.Lch,obj.setting.k,obj.setting.LTratio,obj.setting.coefficient);
            end
        end

        function [modGeomVerts,modSREF,modBREF,modCREF,modXYZREF,modargin_x,unscaledVariables] = calcGeomfromVariables(obj,unscaledVariables)
            [modGeomVerts,~,modSREF,modBREF,modCREF,modXYZREF,modargin_x,unscaledVariables] = obj.geomGenFun(unscaledVariables);
        end
        
   
        function [modGeomVerts,modMeshVerts] = calcApproximatedMeshGeom(obj,unscaledVariables)
            %%%%%%%%%%%設計変数勾配によるMeshとGeomの一次近似の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           scaledVariables = (unscaledVariables(:)'- obj.lb) ./ obj.designScale; 
           modGeomVerts = obj.orgGeom.Points + reshape(obj.gradSurf*(scaledVariables(:)-obj.scaledVar(:)),size(obj.orgGeom.Points));
           if nargout == 2
               modMeshVerts = obj.meshDeformation(obj,modGeomVerts);
           end
        end
        
    
        function viewMesh(obj,modMesh,fig,~)
            %%%%%%%%%%%Meshの表示%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            viewtri = triangulation(obj.orgMesh.ConnectivityList,modMesh);
            figure(fig);
            if nargin < 4
                trisurf(viewtri, 'FaceAlpha', 0, 'EdgeColor', 'black');
            else
                deform = vecnorm(modMesh-obj.orgMesh.Points,2,2);
                trisurf(viewtri, deform ,'FaceAlpha', 0.8,'EdgeAlpha',0.15);
            end
            axis equal;drawnow();
        end

        function viewGeom(obj,modGeom,fig,~)
            %%%%%%%%%%%Geomの表示%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            viewtri = triangulation(obj.orgGeom.ConnectivityList,modGeom);
            figure(fig);
            if nargin < 4
                trisurf(viewtri, 'FaceAlpha', 0, 'EdgeColor', 'black');
            else
                deform = vecnorm(modGeom-obj.orgGeom.Points,2,2);
                trisurf(viewtri, deform ,'FaceAlpha', 0.8,'EdgeAlpha',0.15);
            end
            axis equal;drawnow();
        end

        function checkGeomGenWork(obj,pert,fig)
            %%%%%%%%%%%設計変数が機能しているかチェックする%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = numel(obj.scaledVar);
            randparam = rand(1).*ones(1,ndim);
            randDes = randparam.*obj.designScale+obj.lb;
            modSurforg = obj.geomGenFun(randDes);
            if nargin == 3
                viewtri = triangulation(obj.orgGeom.ConnectivityList,modSurforg);
                figure(fig);clf
                trisurf(viewtri);
                axis equal;drawnow();
                fprintf("randamized parameter : ");
                disp(randDes);
                pause(1);
            end
            for i = 1:ndim
                sampleDes = randparam.*obj.designScale+obj.lb;
                sampleDes(i) = (randparam(i) + pert).*obj.designScale(i)+obj.lb(i);
                modSurf = obj.geomGenFun(sampleDes);
                if all(abs(modSurforg(:)-modSurf(:)) <sqrt(eps))
                    error("Variables No. %d : NOT MOVED",i);
                end
                if nargin == 3
                    viewtri = triangulation(obj.orgGeom.ConnectivityList,modSurf);
                    figure(fig);clf
                    trisurf(viewtri);
                    axis equal;drawnow();
                    pause(1);
                end
            end
        end

        function obj = setCfParameter(obj,Re,Lch,k,LTratio,CfeCoefficient)       
            obj.setting.Re = Re;
            obj.setting.Lch = Lch;
            obj.setting.k = k;
            obj.setting.LTratio = LTratio;
            obj.setting.coefficient = CfeCoefficient;
            for i = 1:numel(obj.flow)
                obj = obj.setCf(i,obj.setting.Re,obj.setting.Lch,obj.setting.k,obj.setting.LTratio,obj.setting.coefficient);
            end
        end
        
        function obj = setOptions(obj,varargin)
            if mod(numel(varargin),1) == 1
                error("input style not match");
            end
            for iter = 1:numel(varargin)/2
                if isfield(obj.setting,varargin{2*iter-1})
                    obj.setting = setfield(obj.setting,varargin{2*iter-1},varargin{2*iter});
                    if strcmp(varargin{2*iter-1},"nCluster")
                        if any(obj.flowNoList(:,3) == 1)
                            obj = obj.makeCluster(obj.setting.nCluster,obj.setting.edgeAngleThreshold);
                        end
                    end
                    if strcmp(varargin{2*iter-1},"Mach")
                        if numel(obj.setting.Mach) ~= numel(obj.setting.alpha) || numel(obj.setting.Mach) ~= numel(obj.setting.beta) || numel(obj.setting.beta) ~= numel(obj.setting.alpha)
                            error("No. of case is not match");
                        end
                        uniMach = unique(obj.setting.Mach);
                        for i = 1:numel(uniMach)
                            obj = obj.flowCondition(i,uniMach(i));
                            obj = obj.setCf(i,obj.setting.Re,obj.setting.Lch,obj.setting.k,obj.setting.LTratio,obj.setting.coefficient);
                        end
                        for i = 1:numel(obj.setting.Mach)
                            obj.flowNoList(i,1) = find(uniMach == obj.setting.Mach(i));
                            obj.flowNoList(i,2) = uniMach(obj.flowNoList(i,1));
                            obj.flowNoList(i,3) = obj.flowNoList(i,2)<1;
                        end
                    end
                else
                    error("Field name is not match");
                end
            end


        end

        function [nextUnscaledVar,obj] = calcNextVariables(obj,objandConsFun,cmin,cmax,varargin)
            if nargin>4
                obj = obj.setOptions(varargin);
            end
            obj.iteration = obj.iteration + 1;
            fprintf("Itaration No. %d : Started\n -- GradientCalcMethod : %s\n -- HessianUpdateMethod : %s\n",obj.iteration,obj.setting.gradientCalcMethod,obj.setting.HessianUpdateMethod);
            %%%%%%%%%%%%指定した評価関数と制約条件における設計変数勾配の計算%%%%%%%%
            %method：設計変数に関する偏微分の計算方法
            % 'direct'-->メッシュを再生成，'chain'-->基準メッシュの変位を補間,
            % 'nonlin'-->approxmatによる近似パネル法行列を使った直接最適化
            % 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %初期点の解析
            nbPanel = sum(obj.paneltype == 1);
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeEquation(obj.setting.wakeLength,obj.setting.n_wake,obj.setting.n_divide);
            end
            obj.approximated = 0;
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            u0 = zeros(nbPanel*size(obj.flowNoList,1),1);
            R0 = zeros(nbPanel*size(obj.flowNoList,1),1);
            for iter = 1:numel(obj.flow)
                alphabuff = obj.setting.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.setting.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                [u0solve,~] = obj.solvePertPotential(iter,alphabuff,betabuff);%ポテンシャルを求める
                [AERODATA0,Cp0,Cfe0,R0solve,obj] = obj.solveFlowForAdjoint(u0solve,iter,alphabuff,betabuff);%ポテンシャルから空力係数を計算
                %結果をマッピング
                lktable = find(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                for i = 1:numel(lktable)
                    u0((lktable(i)-1)*nbPanel+1:lktable(i)*nbPanel,1) = u0solve(nbPanel*(i-1)+1:nbPanel*i,1);
                    R0((lktable(i)-1)*nbPanel+1:lktable(i)*nbPanel,1) = R0solve(nbPanel*(i-1)+1:nbPanel*i,1);
                end

            end
            obj.plotGeometry(1,obj.Cp{1}(:,1),[-2,1]);
            [I0,con0] = objandConsFun(desOrg,AERODATA0,Cp0,Cfe0,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x);
            fprintf("Orginal Objective Value and Constraints:\n")
            disp([I0,con0(:)']);
            obj.history.objVal(obj.iteration) = I0;
            obj.history.conVal(:,obj.iteration) = con0(:);
            fprintf("AERODATA of Itaration No.%d ->\n",obj.iteration);
            disp(vertcat(AERODATA0{:}));
            if not(strcmpi(obj.setting.gradientCalcMethod,"nonlin"))
                %%%
                %明示・非明示随伴方程式法の実装
                %u微分の計算
                pert = sqrt(eps);
                for i = 1:numel(u0)
                    u = u0;
                    u(i) = u(i)+pert;
                    for iter = 1:numel(obj.flow)
                        lktable = find(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        usolve = zeros(nbPanel*numel(lktable),1);
                        for k = 1:numel(lktable)
                            usolve(nbPanel*(k-1)+1:nbPanel*k,1) = u((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1);
                        end
                        alphabuff = obj.setting.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        betabuff = obj.setting.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        [AERODATA,Cp,Cfe,~,obj] = obj.solveFlowForAdjoint(usolve,iter,alphabuff,betabuff);
                    end
                    [I,con] = objandConsFun(desOrg,AERODATA,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x);
                    dI_du(i) = (I-I0)/pert;%評価関数のポテンシャルに関する偏微分
                    if not(isempty(con))
                        dcon_du(:,i) = (con-con0)/pert;
                    end
                end
                for i = 1:size(obj.flowNoList,1)
                    if i == 1
                        if obj.flowNoList(i,3) == 1
                            dR_du = obj.LHS; %亜音速
                        else
                            dR_du = eye(nbPanel);
                        end
                    else
                        if obj.flowNoList(i,3) == 1
                            dR_du = blkdiag(dR_du,obj.LHS);
                        else
                            dR_du = blkdiag(dR_du,eye(nbPanel));
                        end
                    end
                end
                %x微分の計算
                %メッシュの節点勾配を作成
                obj = obj.makeMeshGradient(obj);
                if strcmpi(obj.setting.gradientCalcMethod,'chain')
                    if any(obj.flowNoList(:,3) == 1)
                        obj = obj.calcApproximatedEquation(obj);
                    end
                end
                for i= 1:numel(obj.scaledVar)
                    x = obj.scaledVar;
                    x(i) = obj.scaledVar(i)+pert;
                    des = x.*obj.designScale+obj.lb;
                    [~,modMesh] = obj.calcApproximatedMeshGeom(des);
                    if strcmpi(obj.setting.gradientCalcMethod,'direct')
                        %変数が少ないときは直接作成
                        obj2 = obj.setVerts(modMesh);
                        if any(obj.flowNoList(:,3) == 1)
                            obj2 = obj2.makeEquation(obj.setting.wakeLength,obj.setting.n_wake,obj.setting.n_divide);
                        end
                        obj2.approximated = 0;
                    elseif strcmpi(obj.setting.gradientCalcMethod,'chain')
                        if any(obj.flowNoList(:,3) == 1)
                            obj2 = obj.makeApproximatedInstance(modMesh);
                        else
                            obj2 = obj.setVerts(modMesh);
                        end
                    else
                        error("Supported method is 'direct' or 'chain'.");
                    end
                        
                    %基準面積等の設計変数変化
                    SREF2 = obj.SREF+obj.gradSREF*(x(:)-obj.scaledVar(:));
                    BREF2 = obj.BREF+obj.gradBREF*(x(:)-obj.scaledVar(:));
                    CREF2 = obj.CREF+obj.gradCREF*(x(:)-obj.scaledVar(:));
                    XYZREF2 = obj.XYZREF+(obj.gradXYZREF*(x(:)-obj.scaledVar(:)))';
                    argin_x2 = obj.argin_x+obj.gradArginx*(x(:)-obj.scaledVar(:));
                    obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                    obj2 = obj2.setRotationCenter(XYZREF2);
                    R = zeros(nbPanel*size(obj.flowNoList,1),1);
                    for iter = 1:numel(obj.flow)
                        lktable = find(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        u0solve = zeros(nbPanel*numel(lktable),1);
                        for k = 1:numel(lktable)
                            u0solve(nbPanel*(k-1)+1:nbPanel*k,1) = u0((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1);
                        end
                        alphabuff = obj.setting.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        betabuff = obj.setting.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                        [AERODATA,Cp,Cfe,Rsolve,obj2] = obj2.solveFlowForAdjoint(u0solve,iter,alphabuff,betabuff);
                        for k = 1:numel(lktable)
                            R((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1) = Rsolve(nbPanel*(k-1)+1:nbPanel*k,1);
                        end
                    end
                    
                    [I,con] = objandConsFun(des,AERODATA,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                    dI_dx(i) = (I-I0)/pert;%評価関数の設計変数に関する偏微分
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
                lbf = -obj.scaledVar;
                ubf = 1-obj.scaledVar;
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                [~,~,~,~,lambda] = linprog(objTotalGrad,alin,blin,[],[],lbf,ubf);
                
                rho = 1000;
                if not(isempty(con0))
                    lambdaR = -(dR_du)\(dI_du+lambda.ineqlin'*[-dcon_du;dcon_du])';
                    Lorg = I0 + lambda.ineqlin'*[-con0;con0];
                    penaltyorg = rho*sum(max(max(0,cmin-con0),max(con0-cmax)));
                    dL_dx = dI_dx+lambdaR'*dR_dx+lambda.ineqlin'*[-dcon_dx;dcon_dx];
                else
                    lambdaR = -(dR_du)\(dI_du)';
                    Lorg = I0 ;
                    penaltyorg = 0;
                    dL_dx = dI_dx+lambdaR'*dR_dx;
                end
                obj.history.LagrangianVal(obj.iteration) = Lorg;
                obj.history.penaltyVal(obj.iteration) = penaltyorg;
                [dxscaled] = fmincon(@(x)obj.fminconObj(x,(obj.setting.betaLM)*obj.Hessian+(1-obj.setting.betaLM)*diag(diag(obj.Hessian)),dL_dx),zeros(numel(obj.scaledVar),1),alin,blin,[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.setting.TrustRegion),options);
                dx = dxscaled(:)'.*obj.designScale;
                
                %精度評価
                desdx = obj.scaledVar.*obj.designScale+obj.lb+dx(:)';
                [modSurfdx,~,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx,desdx] = obj.geomGenFun(desdx);
                modMeshdx = obj.meshDeformation(obj,modSurfdx);
                objdx = obj.setVerts(modMeshdx);
                if any(obj.flowNoList(:,3) == 1)
                    objdx = objdx.makeEquation(obj.setting.wakeLength,obj.setting.n_wake,obj.setting.n_divide);
                end
                objdx = objdx.setREFS(SREFdx,BREFdx,CREFdx);
                objdx = objdx.setRotationCenter(XYZREFdx);
                for iter = 1:numel(obj.flow)
                    alphabuff = obj.setting.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.setting.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    [udx,~] = objdx.solvePertPotential(iter,alphabuff,betabuff);%ポテンシャルを求める
                    [AERODATA,Cp,Cfe,~,objdx] = objdx.solveFlowForAdjoint(udx,iter,alphabuff,betabuff);%ポテンシャルから空力係数を計算
                end
               
                [Idx,condx] = objandConsFun(desdx,AERODATA,Cp,Cfe,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx);
                if not(isempty(con0))
                    Ldx = Idx + lambda.ineqlin'*[-condx;condx];
                    if penaltyorg < sqrt(eps)
                        penaltydx = 0;
                    else
                       penaltydx = rho*sum(max(max(0,cmin-condx),max(condx-cmax))); 
                    end
                else
                    Ldx = Idx ;
                    penaltydx = 0;
                end
                acc = (Ldx-Lorg)/(0.5*(dxscaled(:)'*obj.Hessian*dxscaled(:))+dL_dx*dxscaled(:));
               
                obj.history.LpredictVal(obj.iteration) = Ldx;
                obj.history.trAccuracyVal(obj.iteration) = acc;
                obj.history.dxNormVal(obj.iteration) = norm(dxscaled);
                fprintf("Variables:\n")
                disp(desOrg);
                fprintf("------>\n")
                disp(desdx);
                fprintf("Objective and Constraints:\n")
                disp([I0,con0(:)']);
                fprintf("------>\n");
                disp([Idx,condx(:)']);
                fprintf("dx norm :%f\nLagrangian Value : %f -> %f\nPenalty Value : %f -> %f\nHessian Approximation Accuracy:%f\n",norm(dxscaled),Lorg,Ldx,penaltyorg,penaltydx,acc);
                dx = desdx-desOrg;


                if strcmpi(obj.setting.HessianUpdateMethod,"SR1")
                   updateFunction = @obj.SR1;
                elseif strcmpi(obj.setting.HessianUpdateMethod,"BFGS")
                    updateFunction = @obj.BFGS;
                elseif strcmpi(obj.setting.HessianUpdateMethod,"DFP")
                    updateFunction = @obj.DFP;
                elseif strcmpi(obj.setting.HessianUpdateMethod,"Broyden")
                    updateFunction = @obj.Broyden;
                end
                %Hessianの更新
                obj.history.scaledVar(obj.iteration,:) = obj.scaledVar;
                obj.history.dL_dx(obj.iteration,:) = dL_dx;
                if size(obj.history.scaledVar,1) > 1
                    n_iter = size(obj.history.scaledVar,1)-1;
                    if n_iter > obj.setting.nMemory
                        n_iter = obj.setting.nMemory;
                    end
                    obj.Hessian = obj.setting.H0;
                    for i = n_iter:-1:1
                        s = obj.history.scaledVar(end-(i-1),:)-obj.history.scaledVar(end-i,:);
                        y = obj.history.dL_dx(end-(i-1),:)-obj.history.dL_dx(end-i,:);
                        if norm(s)>sqrt(eps) && norm(y)>sqrt(eps)
                            obj.Hessian = updateFunction(s,y,obj.Hessian);
                        end
                    end
                    s = obj.history.scaledVar(end,:)-obj.history.scaledVar(end-1,:);
                    y = obj.history.dL_dx(end,:)-obj.history.dL_dx(end-1,:);
                end
            else
                %近似行列による直接最適化
                obj = obj.makeMeshGradient();
                %%%%TODO : 超音速のみの時の実装を考える
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.calcApproximatedEquation();
                end
                lbf = -obj.scaledVar;
                ubf = 1-obj.scaledVar;
                
                if isempty(con0)
                    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');
                    [dxscaled,fval,exitflag,output,lambda] = fmincon(@(x)obj.nonlnObj(x,obj,objandConsFun,cmin,cmax),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.nonlnNormCon(x,obj.optimization.TRmax),options);
                else
                    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');
                    [dxscaled,fval,exitflag,output,lambda] = fmincon(@(x)obj.nonlnObj(x,obj,objandConsFun,cmin,cmax),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.nonlnCon(x,obj,objandConsFun,cmin,cmax,obj.optimization.TRmax),options);                
                end
                dx = dxscaled(:)'.*obj.designScale;
                %精度評価
                desdx = desOrg+dx(:)';
                [modSurfdx,~,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx,desdx] = obj.geomGenFun(desdx);
                modMeshdx = obj.meshDeformation(obj,modSurfdx);
                objdx = obj.setVerts(modMeshdx);
                if any(obj.flowNoList(:,3) == 1)
                    objdx = objdx.makeEquation(obj.setting.wakeLength,obj.setting.n_wake,obj.setting.n_divide);
                end
                objdx = objdx.setREFS(SREFdx,BREFdx,CREFdx);
                objdx = objdx.setRotationCenter(XYZREFdx);
                objdx = objdx.solveFlow(1,obj.setting.alpha,obj.setting.beta);
                %[udx,~] = obj.solvePertPotential(1,obj.setting.alpha,obj.setting.beta);%ポテンシャルを求める
                %[AERODATA,Cp,Cfe,Rdx,~] = objdx.solveFlowForAdjoint(udx,1,obj.setting.alpha,obj.setting.beta);%ポテンシャルから空力係数を計算
                [Idx,condx] = objandConsFun(desdx,objdx.AERODATA,objdx.Cp,objdx.Cfe,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx);
                fprintf("Variables:\n")
                disp(desOrg);
                fprintf("------>\n")
                disp(desdx);
                fprintf("Objective and Constraints:\n")
                disp([I0,con0(:)']);
                fprintf("------>\n");
                disp([Idx,condx(:)']);
                fprintf("dx norm :%f\n",norm((desdx-desOrg)./obj.designScale));
            end
            nextUnscaledVar = obj.scaledVar.*obj.designScale+obj.lb+dx;
            fprintf("Itaration No. %d : Completed\n",obj.iteration);
        end


        function plotOptimizationState(obj,fig)
            figure(fig);clf;
            plotiter = 1:obj.iteration;
            subplot(5,1,1);grid on;
            plot(plotiter,obj.history.objVal,"-o");
            subplot(5,1,2);grid on;
            plot(plotiter,obj.history.conVal,"-o");
            subplot(5,1,3);grid on;
            plot(plotiter,obj.history.LagrangianVal,"-o");
            subplot(5,1,4);grid on;
            plot(plotiter,obj.history.trAccuracyVal,"-o");
            subplot(5,1,5);grid on;
            plot(plotiter,obj.history.dxNormVal,"-o");
            drawnow();
        end
    end

    methods(Static)
        function res = nonlnObj(dx,obj,objandConsFun,cmin,cmax)
                x = dx(:)'+obj.scaledVar;
                des = x.*obj.designScale+obj.lb;
                [~,modMesh] = obj.variables2Mesh(des,'linear');
                if any(obj.flowNoList(:,3) == 1)
                    obj2 = obj.makeApproximatedInstance(modMesh);
                else
                    obj2 = obj.setVerts(modMesh);
                end
                SREF2 = obj.optSREF+obj.gradSREF*dx(:);
                BREF2 = obj.optBREF+obj.gradBREF*dx(:);
                CREF2 = obj.optCREF+obj.gradCREF*dx(:);
                XYZREF2 = obj.optXYZREF+(obj.gradXYZREF*dx(:))';
                argin_x2 = obj.argin_x+obj.gradArginx*dx(:);
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                obj2 = obj2.setRotationCenter(XYZREF2);
                obj2 = obj2.solveFlow(1,obj.setting.alpha,obj.setting.beta);
                [I,con] = objandConsFun(des,obj2.AERODATA,obj2.Cp,obj2.Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                res = I;
        end

        function [c,ceq] = nonlnCon(dx,obj,objandConsFun,cmin,cmax,TR)
                x = dx(:)'+obj.scaledVar;
                des = x.*obj.designScale+obj.lb;
                [~,modMesh] = obj.variables2Mesh(des,'linear');
                if any(obj.flowNoList(:,3) == 1)
                    obj2 = obj.makeApproximatedInstance(modMesh);
                else
                    obj2 = obj.setVerts(modMesh);
                end
                SREF2 = obj.optSREF+obj.gradSREF*dx(:);
                BREF2 = obj.optBREF+obj.gradBREF*dx(:);
                CREF2 = obj.optCREF+obj.gradCREF*dx(:);
                XYZREF2 = obj.optXYZREF+(obj.gradXYZREF*dx(:))';
                argin_x2 = obj.argin_x+obj.gradArginx*dx(:);
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                obj2 = obj2.setRotationCenter(XYZREF2);
                obj2 = obj2.solveFlow(1,obj.setting.alpha,obj.setting.beta);
                [I,con] = objandConsFun(des,obj2.AERODATA,obj2.Cp,obj2.Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                ceq = [];
                c = [-con+cmin;con-cmax];
                c(end+1) = sum(dx.^2)-(TR)^2;   
        end
        
        function [c,ceq] = nonlnNormCon(dx,TR)
            c = sum(dx.^2)-(TR)^2;   
            ceq = [];
        end

        function res = fminconObj(dx,H,g)
            %%%SQPの更新ベクトルの目的関数%%%%%
            %dx：更新ベクトル
            %H：ヘシアン
            %g：勾配
            %%%%%%%%%%%%%%%%%%%%%%
            res = 0.5 * dx(:)'*H*dx(:) + g*dx(:);
        end

        function [c,ceq] = fminconNlc(dx,TR)
            ceq = [];
            c = sum(dx.^2)-(TR)^2;
        end

        function [c,ceq] = fminconNlc2(dx,TR,H,dx1scaled)
            c = sum(dx.^2)-(TR)^2;
            ceq = dx1scaled'*H*dx;
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
            md.x = scatteredInterpolant(obj.orgGeom.Points,modSurfVerts(:,1)-obj.orgGeom.Points(:,1),'linear','linear');
            md.y = scatteredInterpolant(obj.orgGeom.Points,modSurfVerts(:,2)-obj.orgGeom.Points(:,2),'linear','linear');
            md.z = scatteredInterpolant(obj.orgGeom.Points,modSurfVerts(:,3)-obj.orgGeom.Points(:,3),'linear','linear');
            dVerts(:,1) = md.x(obj.orgMesh.Points);
            dVerts(:,2) = md.y(obj.orgMesh.Points);
            dVerts(:,3) = md.z(obj.orgMesh.Points);
            modVerts = obj.orgMesh.Points+dVerts;
            con = obj.orgGeom.ConnectivityList;
        end
        
        function obj = makeMeshGradient(obj)
            %%%%%%%%%%%設計変数勾配による標本表面の近似関数の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pert = 0.01./obj.designScale;
            ndim = numel(obj.scaledVar);
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            [surforg,~,~,~,~,~,~,desOrg] = obj.geomGenFun(desOrg);
            for i = 1:ndim
                sampleDes = obj.scaledVar.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.scaledVar(i) + pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFf,BREFf,CREFf,XYZREFf,argin_xf,desBuff] = obj.geomGenFun(sampleDes);
                pertf = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurff = dmodSurf(:);
                sampleDes = obj.scaledVar.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.scaledVar(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFr,BREFr,CREFr,XYZREFr,argin_xr,desBuff] = obj.geomGenFun(sampleDes);
                pertr = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurfr = dmodSurf(:);
                obj.gradSurf(:,i) = (sampleSurff-sampleSurfr)./(pertf-pertr);
                obj.gradSREF(i) = (SREFf-SREFr)./(pertf-pertr);
                obj.gradBREF(i) = (BREFf-BREFr)./(pertf-pertr);
                obj.gradCREF(i) = (CREFf-CREFr)./(pertf-pertr);
                obj.gradXYZREF(:,i) = (XYZREFf(:)-XYZREFr(:))./(pertf-pertr);
                obj.gradArginx(:,i) = (argin_xf(:)-argin_xr(:))./(pertf-pertr);

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

        function approxmatedObj = makeApproximatedInstance(obj,modifiedVerts)
                nPanel = numel(obj.paneltype);
                nbPanel = sum(obj.paneltype == 1);
                approxmatedObj = obj.setVerts(modifiedVerts);
                approxmatedObj.mu2v{1} = sparse(nPanel,nbPanel);
                approxmatedObj.mu2v{2} = sparse(nPanel,nbPanel);
                approxmatedObj.mu2v{3} = sparse(nPanel,nbPanel);
            
                for i = 1:nPanel
                    if approxmatedObj.paneltype(i) == 1
                        CPmat =approxmatedObj.center(approxmatedObj.cluster{i},1:3);
                        pnt = approxmatedObj.center(i,:);
                        m = approxmatedObj.tri.Points(approxmatedObj.tri.ConnectivityList(i,1),:)'-pnt(:);
                        m = m./norm(m);
                        l = cross(m,approxmatedObj.normal(i,:)');
                        Minv = [l,m,approxmatedObj.normal(i,:)'];
                        lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                        bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                        Bmat=pinv(bb,sqrt(eps));
                        Vnmat = Minv(:,[1,2])*[1,0,0,0;0,1,0,0]*Bmat;
                        for iter = 1:3
                            approxmatedObj.mu2v{iter}(i,approxmatedObj.IndexPanel2Solver(approxmatedObj.cluster{i})) = Vnmat(iter,:);
                        end
                    end
                end
            
                
                for wakeNo = 1:numel(obj.wakeline)
                    theta = linspace(pi,0,numel(obj.wakeline{wakeNo}.edge)*obj.LLT.n_interp+1);
                    iter = 1;
                    obj.LLT.sp{wakeNo} = [];
                    obj.LLT.calcMu{wakeNo} = zeros(1,nbPanel);
                    s = zeros(1,numel(obj.wakeline{wakeNo}.edge));
                    for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                        s(edgeNo+1) = s(edgeNo) + norm([obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),2:3)-obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),2:3)]);
                    end
                    for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                        obj.LLT.sp{wakeNo} =  [obj.LLT.sp{wakeNo},(s(edgeNo)+s(edgeNo+1))./2];
                    end
                    sd = (s(end)-s(1))*(cos(theta)./2+0.5)+s(1);
                    obj.LLT.sinterp{wakeNo} = (sd(2:end)+sd(1:end-1))./2;
                    obj.LLT.yinterp{wakeNo} = interp1(s,obj.tri.Points(obj.wakeline{wakeNo}.edge(:),2),obj.LLT.sinterp{wakeNo},'linear','extrap');
                    obj.LLT.zinterp{wakeNo} = interp1(s,obj.tri.Points(obj.wakeline{wakeNo}.edge(:),3),obj.LLT.sinterp{wakeNo},'linear','extrap');
                    yd = interp1(s,obj.tri.Points(obj.wakeline{wakeNo}.edge(:),2),sd,'linear','extrap');
                    zd = interp1(s,obj.tri.Points(obj.wakeline{wakeNo}.edge(:),3),sd,'linear','extrap');
                    obj.LLT.phiinterp{wakeNo} = atan((zd(2:end)-zd(1:end-1))./(yd(2:end)-yd(1:end-1)));
                    obj.LLT.spanel{wakeNo} = (sd(2:end)-sd(1:end-1))./2;
                    
                end
                obj.LLT.Qij = obj.Calc_Q(horzcat(obj.LLT.yinterp{:}),horzcat(obj.LLT.zinterp{:}),horzcat(obj.LLT.phiinterp{:}),horzcat(obj.LLT.spanel{:}),obj.halfmesh);


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

        %%%%%%%%%%%%%%%%%%%%%%%%
        %準ニュートン法のアップデートの手実装
        %s 設計変数の変化量
        %y ヤコビアンの変化
        %%%%%%%%%%%%%%%%%%%%%%%%
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
        function Bkp1 = DFP(s,y,Bk)
            s = s(:);
            y=y(:);
            ndim = numel(s);
            if s'*y<=0
                Bkp1 = Bk;
            else
                Bkp1 = (eye(ndim)-(y*s')/(y'*s))*Bk*(eye(ndim)-(s*y')/(y'*s))+(y*y')/(y'*s);
                coeB=(s'*y)/(s'*Bkp1*s);
                if coeB>0 && coeB<1
                    Bkp1=coeB.*Bkp1;
                end
            end
        end
        function Bkp1 = BFGS(s,y,Bk)
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
        function Bkp1 = Broyden(s,y,Bk)
            s = s(:);
            y=y(:);
            if s'*y<=0
                Bkp1 = Bk;
            else
                Bkp1 = Bk+((y-Bk*s)/(s'*s))*s';
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
                end
            else
                Bkp1 = Bk+(y-Bk*s)*(y-Bk*s)'/((y-Bk*s)'*s);
            end
        end
        
        function Bkp1 = Barzilai_Borwein(s,y,Bk)
            s = s(:);
            y=y(:);
            delta = 100.0;
            alpha = (s'*y)/(s'*s);
            if alpha <= 0
                alpha = norm(s)/norm(y);
            end
            alphak = min(alpha, delta / norm(y));
            Bk = eye(numel(s));
            Bkp1 = alphak * Bk;
        end

    end

end