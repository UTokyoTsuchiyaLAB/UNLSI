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
        approximated
        optimization
        optSREF
        optBREF
        optCREF
        optXYZREF
        gradSREF
        gradBREF
        gradCREF
        gradXYZREF
        gradArginx
        argin_x
    end

    methods(Access = public)

        function obj = UNGRADE(surfGenFun,designVariables,lb,ub,orgVerts,orgCon,surfID,wakelineID,halfmesh)
            %%%%%%%%%%%%メッシュの登録と基準サーフェス生成関数の登録%%%%%%%%%%%%%%%
            %orgVerts,orgCon　実際に解析を行うメッシュ（openVSPのCFDツール等で生成した（基本的に）オーバーラップのない非構造メッシュ）
            %meshGenFun　設計変数（スケールされていない）から基準サーフェス（解析メッシュと同じ表面を持つ、トリムされていないメッシュ）
            %designVariables : 設計変数
            %lb：設計変数の下限
            %ub：設計変数の上限
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj = obj@UNLSI(orgVerts,orgCon,surfID,wakelineID,halfmesh)
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            obj.orgMesh = triangulation(orgCon,orgVerts);
            obj.surfGenFun = surfGenFun;
            [orgSurfVerts,orgSurfCon,obj.optSREF,obj.optBREF,obj.optCREF,obj.optXYZREF,obj.argin_x,desOrg] = obj.surfGenFun(designVariables(:)');
            obj = obj.setREFS(obj.optSREF,obj.optBREF,obj.optCREF);
            obj = obj.setRotationCenter(obj.optXYZREF);
            obj.designVariables = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgSurf =  triangulation(orgSurfCon,orgSurfVerts);
        end

        function obj = setMeshGenFun(obj,meshGenFun)
            obj.meshGenFun = meshGenFun;
        end

        function obj = modifyMesh(obj,designVariables,orgVerts,orgCon,surfID,wakelineID)
            obj = obj.setMesh(orgVerts,orgCon,surfID,wakelineID);
            obj.orgMesh = triangulation(orgCon,orgVerts);
            [orgSurfVerts,orgSurfCon,obj.optSREF,obj.optBREF,obj.optCREF,obj.optXYZREF,obj.argin_x,desOrg] = obj.surfGenFun(designVariables(:)');
            obj = obj.setREFS(obj.optSREF,obj.optBREF,obj.optCREF);
            obj = obj.setRotationCenter(obj.optXYZREF);
            obj.designVariables = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgSurf =  triangulation(orgSurfCon,orgSurfVerts);
            obj = obj.makeCluster(obj.unlsiParam.nCluster,obj.unlsiParam.edgeAngleThreshold);
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
            md.x = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,1)-obj.orgSurf.Points(:,1),'linear','linear');
            md.y = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,2)-obj.orgSurf.Points(:,2),'linear','linear');
            md.z = scatteredInterpolant(obj.orgSurf.Points,modSurfVerts(:,3)-obj.orgSurf.Points(:,3),'linear','linear');
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
            [surforg,~,obj.optSREF,obj.optBREF,obj.optCREF,obj.optXYZREF,obj.argin_x,desOrg] = obj.surfGenFun(desOrg);
            for i = 1:ndim
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) + pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFf,BREFf,CREFf,XYZREFf,argin_xf,desBuff] = obj.surfGenFun(sampleDes);
                pertf = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurff = dmodSurf(:);
                sampleDes = obj.designVariables.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.designVariables(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFr,BREFr,CREFr,XYZREFr,argin_xr,desBuff] = obj.surfGenFun(sampleDes);
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
            modSurforg = obj.surfGenFun(randDes);
            viewtri = triangulation(obj.orgSurf.ConnectivityList,modSurforg);
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
                if all(abs(modSurforg(:)-modSurf(:)) <sqrt(eps))
                    warning("Variables No. %d : NOT MOVED",i);
                end
                figure(fig);clf
                trisurf(viewtri);
                axis equal;drawnow();
                pause(1)
            end
        end

        function obj = setOptFlowCondition(obj,Mach,alpha,beta,wakeLength,n_wake,n_divide,nCluster,edgeAngleThreshold,Re,Lch,k,LTratio,CfeCoefficient)
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
            obj = obj.setCf(1,obj.unlsiParam.Re,obj.unlsiParam.Lch,obj.unlsiParam.k,obj.unlsiParam.LTratio,obj.unlsiParam.coefficient);
            obj = obj.makeCluster(obj.unlsiParam.nCluster,obj.unlsiParam.edgeAngleThreshold);
        end

        function obj = setOptimization(obj,H0,TR,TRmin,TRmax,method,nMemory)
            obj.optimization.H0 = H0;
            obj.optimization.H = H0;
            obj.optimization.nMemory = nMemory;
            obj.optimization.TR = TR;
            obj.optimization.TRmin = TRmin;
            obj.optimization.TRmax = TRmax;
            obj.optimization.xScaled = [];
            obj.optimization.dL_dx = [];
            if strcmpi(method,"SR1")
                obj.optimization.updateFunction = @obj.SR1;
            elseif strcmpi(method,"SSR1")
                obj.optimization.updateFunction = @obj.SSR1;
            elseif strcmpi(method,"BFGS")
                obj.optimization.updateFunction = @obj.BFGS;
            elseif strcmpi(method,"MBFGS")
                obj.optimization.updateFunction = @obj.MBFGS;
            elseif strcmpi(method,"DBFGS")
                obj.optimization.updateFunction = @obj.DBFGS;
            elseif strcmpi(method,"SR1_BFGS")
                obj.optimization.updateFunction = @(s,y,H)obj.SR1_BFGS(obj,s,y,H);
            elseif strcmpi(method,"SSR1_MBFGS")
                obj.optimization.updateFunction = @(s,y,H)obj.SSR1_MBFGS(obj,s,y,H);
            end
        end

        function [dx,obj] = finddx(obj,objandConsFun,method,cmin,cmax)
            %%%%%%%%%%%%指定した評価関数と制約条件における設計変数勾配の計算%%%%%%%%
            %method：設計変数に関する偏微分の計算方法
            % 'direct'-->メッシュを再生成，'chain'-->基準メッシュの変位を補間,
            % 'nonlin'-->approxmatによる近似パネル法行列を使った直接最適化
            % 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %初期点の解析
            obj = obj.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
            obj.approximated = 0;
            desOrg = obj.designVariables.*obj.designScale+obj.lb;
            obj = obj.setREFS(obj.optSREF,obj.optBREF,obj.optCREF);
            obj = obj.setRotationCenter(obj.optXYZREF);
            [u0,~] = obj.solvePertPotential(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルを求める
            [AERODATA0,Cp0,Cfe0,R0,obj] = obj.solveFlowForAdjoint(u0,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルから空力係数を計算
            obj.plotGeometry(1,obj.Cp,[-2,1]);
            [I0,con0] = objandConsFun(desOrg,AERODATA0,Cp0,Cfe0,obj.optSREF,obj.optBREF,obj.optCREF,obj.optXYZREF,obj.argin_x);

            if not(strcmpi(method,"nonlin"))
                %%%
                %明示・非明示随伴方程式法の実装
                %u微分の計算
                pert = sqrt(eps);
                for i = 1:numel(u0)
                    u = u0;
                    u(i) = u(i)+pert;
                    [AERODATA,Cp,Cfe] = obj.solveFlowForAdjoint(u,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                    [I,con] = objandConsFun(desOrg,AERODATA,Cp,Cfe,obj.optSREF,obj.optBREF,obj.optCREF,obj.optXYZREF,obj.argin_x);
                    dI_du(i) = (I-I0)/pert;%評価関数のポテンシャルに関する偏微分
                    if not(isempty(con))
                        dcon_du(:,i) = (con-con0)/pert;
                    end
                end
                dR_du = obj.LHS;
                %x微分の計算
                %メッシュの節点勾配を作成
                obj = obj.makeMeshGradient();
                if strcmpi(method,'chain')
                    obj = obj.calcApproximatedEquation();
                end
                for i= 1:numel(obj.designVariables)
                    x = obj.designVariables;
                    x(i) = obj.designVariables(i)+pert;
                    des = x.*obj.designScale+obj.lb;
                    [~,modMesh] = obj.variables2Mesh(des,'linear');
                    if strcmpi(method,'direct')
                        %変数が少ないときは直接作成
                        obj2 = obj.setVerts(modMesh);
                        obj2 = obj2.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
                        obj2.approximated = 0;
                    elseif strcmpi(method,'chain')
                        obj2 = obj.makeApproximatedInstance(modMesh);
                    else
                        error("Supported method is 'direct' or 'chain'.");
                    end
                        
                    %基準面積等の設計変数変化
                    SREF2 = obj.optSREF+obj.gradSREF*(x(:)-obj.designVariables(:));
                    BREF2 = obj.optBREF+obj.gradBREF*(x(:)-obj.designVariables(:));
                    CREF2 = obj.optCREF+obj.gradCREF*(x(:)-obj.designVariables(:));
                    XYZREF2 = obj.optXYZREF+(obj.gradXYZREF*(x(:)-obj.designVariables(:)))';
                    argin_x2 = obj.argin_x+obj.gradArginx*(x(:)-obj.designVariables(:));
                    obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                    obj2 = obj2.setRotationCenter(XYZREF2);
                    [AERODATA,Cp,Cfe,R] = obj2.solveFlowForAdjoint(u0,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
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
                lbf = -obj.designVariables;
                ubf = 1-obj.designVariables;
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed');
                [dxscaled,fval,exitflag,output,lambda] = fmincon(@(dx)obj.fminconObj(dx,obj.optimization.H,objTotalGrad),zeros(numel(obj.designVariables),1),alin,blin,[],[],lbf,ubf,@(dx)obj.fminconNlc(dx,obj.optimization.TR),options);
                if not(isempty(con0))
                    lambdaR = -(dR_du)\(dI_du+lambda.ineqlin'*[-dcon_du;dcon_du])';
                    Lorg = I0 + lambda.ineqlin'*[-con0;con0];
                    dL_dx = dI_dx+lambdaR'*dR_dx+lambda.ineqlin'*[-dcon_dx;dcon_dx];
                else
                    lambdaR = -(dR_du)\(dI_du)';
                    Lorg = I0 ;
                    dL_dx = dI_dx+lambdaR'*dR_dx;
                end
                dx = dxscaled(:)'.*obj.designScale;
                
                %精度評価
                desdx = x.*obj.designScale+obj.lb+dx(:)';
                [modSurfdx,~,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx,desBuff] = obj.surfGenFun(desdx);
                modMeshdx = obj.meshDeformation(modSurfdx);
                objdx = obj.setVerts(modMeshdx);
                objdx = objdx.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
                objdx = objdx.setREFS(SREFdx,BREFdx,CREFdx);
                objdx = objdx.setRotationCenter(XYZREFdx);
                %objdx = objdx.solveFlow(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                [udx,~] = obj.solvePertPotential(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルを求める
                [AERODATA,Cp,Cfe,Rdx,~] = objdx.solveFlowForAdjoint(udx,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルから空力係数を計算
                [Idx,condx] = objandConsFun(des,AERODATA,Cp,Cfe,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx);
                if not(isempty(con0))
                    Ldx = Idx + lambda.ineqlin'*[-condx;condx];
                else
                    Ldx = Idx ;
                end
                if not(isempty(con0))
                    acc = (Ldx-Lorg)/(fval+(lambda.ineqlin'*[-dcon_dx;dcon_dx])*dxscaled(:));
                else
                    acc = (Ldx-Lorg)/(fval);
                end
                fprintf("Variables:\n")
                disp(desOrg);
                fprintf("------>\n")
                disp(desdx);
                fprintf("Objective and Constraints:\n")
                disp([I0,con0(:)']);
                fprintf("------>\n");
                disp([Idx,condx(:)']);
                fprintf("dx norm :%f\nLagrangian Value : %f -> %f\nHessian Approximation Accuracy:%f\n",norm(dxscaled),Lorg,Ldx,acc);
                if acc < 0.15
                    obj.optimization.TR = obj.optimization.TR * 0.9;
                elseif acc > 0.5
                    obj.optimization.TR = obj.optimization.TR / 0.9;
                end
                if obj.optimization.TR > obj.optimization.TRmax
                    obj.optimization.TR = obj.optimization.TRmax;
                elseif obj.optimization.TR < obj.optimization.TRmin
                    obj.optimization.TR = obj.optimization.TRmin;
                end

                %Hessianの更新
                obj.optimization.xScaled = [obj.optimization.xScaled;obj.designVariables];
                obj.optimization.dL_dx = [obj.optimization.dL_dx;dL_dx];
                if size(obj.optimization.xScaled,1) > 1
                    n_iter = size(obj.optimization.xScaled,1)-1;
                    if n_iter > obj.optimization.nMemory
                        n_iter = obj.optimization.nMemory;
                    end
                    obj.optimization.H = obj.optimization.H0;
                    for i = 1:n_iter
                        s = obj.optimization.xScaled(end-(i-1),:)-obj.optimization.xScaled(end-i,:);
                        y = obj.optimization.dL_dx(end-(i-1),:)-obj.optimization.dL_dx(end-i,:);
                        obj.optimization.H = obj.optimization.updateFunction(s,y,obj.optimization.H);
                    end
                end
            else
                %近似行列による直接最適化
                obj = obj.makeMeshGradient();
                obj = obj.calcApproximatedEquation();
                lbf = -obj.designVariables;
                ubf = 1-obj.designVariables;
                
                if isempty(con0)
                    %options = optimoptions(@patternsearch,'MaxIterations',100,'Display','iter');
                    %dxscaled = patternsearch(@(dx)obj.nonlnObj(dx,obj,objandConsFun,cmin,cmax),zeros(numel(obj.designVariables),1),[],[],[],[],lbf,ubf,[],options);
                    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');
                    [dxscaled,fval,exitflag,output,lambda] = fmincon(@(dx)obj.nonlnObj(dx,obj,objandConsFun,cmin,cmax),zeros(numel(obj.designVariables),1),[],[],[],[],lbf,ubf,@(dx)obj.nonlnNormCon(dx,obj.optimization.TRmax),options);
                else
                    %options = optimoptions(@patternsearch,'MaxIterations',100,'Display','iter');
                    %dxscaled = patternsearch(@(dx)obj.nonlnObj(dx,obj,objandConsFun,cmin,cmax),zeros(numel(obj.designVariables),1),[],[],[],[],lbf,ubf,@(dx)obj.nonlnCon(dx,obj,objandConsFun,cmin,cmax),options);
                    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');
                    [dxscaled,fval,exitflag,output,lambda] = fmincon(@(dx)obj.nonlnObj(dx,obj,objandConsFun,cmin,cmax),zeros(numel(obj.designVariables),1),[],[],[],[],lbf,ubf,@(dx)obj.nonlnCon(dx,obj,objandConsFun,cmin,cmax,obj.optimization.TRmax),options);                
                end
                %dxscaled = ga(@(dx)obj.nonlnObj(dx,obj,objandConsFun,cmin,cmax),numel(obj.designVariables),[],[],[],[],lbf,ubf,@(dx)obj.nonlnCon(dx,obj,objandConsFun,cmin,cmax,1),options);
                dx = dxscaled(:)'.*obj.designScale;
                %精度評価
                desdx = desOrg+dx(:)';
                [modSurfdx,~,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx,desBuff] = obj.surfGenFun(desdx);
                modMeshdx = obj.meshDeformation(modSurfdx);
                objdx = obj.setVerts(modMeshdx);
                objdx = objdx.makeEquation(obj.unlsiParam.wakeLength,obj.unlsiParam.n_wake,obj.unlsiParam.n_divide);
                objdx = objdx.setREFS(SREFdx,BREFdx,CREFdx);
                objdx = objdx.setRotationCenter(XYZREFdx);
                %objdx = objdx.solveFlow(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                [udx,~] = obj.solvePertPotential(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルを求める
                [AERODATA,Cp,Cfe,Rdx,~] = objdx.solveFlowForAdjoint(udx,1,obj.unlsiParam.alpha,obj.unlsiParam.beta);%ポテンシャルから空力係数を計算
                [Idx,condx] = objandConsFun(desdx,AERODATA,Cp,Cfe,SREFdx,BREFdx,CREFdx,XYZREFdx,argin_xdx);
                fprintf("Variables:\n")
                disp(desOrg);
                fprintf("------>\n")
                disp(desdx);
                fprintf("Objective and Constraints:\n")
                disp([I0,con0(:)']);
                fprintf("------>\n");
                disp([Idx,condx(:)']);
                fprintf("dx norm :%f\n",norm(dxscaled));
            end

           
        end

        function [dx,obj] = updateVariables(obj,FcnObjandCon,method,cmin,cmax)
            [dx,obj] = obj.finddx(FcnObjandCon,method,cmin,cmax);
            newDes = obj.designVariables.*obj.designScale+obj.lb+dx;
            obj = obj.modifyMeshfromVariables(newDes);
            obj = obj.setCf(1,obj.unlsiParam.Re,obj.unlsiParam.Lch,obj.unlsiParam.k,obj.unlsiParam.LTratio,obj.unlsiParam.coefficient);
            obj = obj.makeCluster(obj.unlsiParam.nCluster,obj.unlsiParam.edgeAngleThreshold);
        end
    end

    methods(Static)
        function res = nonlnObj(dx,obj,objandConsFun,cmin,cmax)
                x = dx(:)'+obj.designVariables;
                des = x.*obj.designScale+obj.lb;
                [~,modMesh] = obj.variables2Mesh(des,'linear');
                obj2 = obj.makeApproximatedInstance(modMesh);
                SREF2 = obj.optSREF+obj.gradSREF*dx(:);
                BREF2 = obj.optBREF+obj.gradBREF*dx(:);
                CREF2 = obj.optCREF+obj.gradCREF*dx(:);
                XYZREF2 = obj.optXYZREF+(obj.gradXYZREF*dx(:))';
                argin_x2 = obj.argin_x+obj.gradArginx*dx(:);
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                obj2 = obj2.setRotationCenter(XYZREF2);
                obj2 = obj2.solveFlow(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
                [I,con] = objandConsFun(des,obj2.AERODATA,obj2.Cp,obj2.Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                res = I;
        end

        function [c,ceq] = nonlnCon(dx,obj,objandConsFun,cmin,cmax,TR)
                x = dx(:)'+obj.designVariables;
                des = x.*obj.designScale+obj.lb;
                [~,modMesh] = obj.variables2Mesh(des,'linear');
                obj2 = obj.makeApproximatedInstance(modMesh);
                SREF2 = obj.optSREF+obj.gradSREF*dx(:);
                BREF2 = obj.optBREF+obj.gradBREF*dx(:);
                CREF2 = obj.optCREF+obj.gradCREF*dx(:);
                XYZREF2 = obj.optXYZREF+(obj.gradXYZREF*dx(:))';
                argin_x2 = obj.argin_x+obj.gradArginx*dx(:);
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2);
                obj2 = obj2.setRotationCenter(XYZREF2);
                obj2 = obj2.solveFlow(1,obj.unlsiParam.alpha,obj.unlsiParam.beta);
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
            res = 0.5 * dx(:)'*H*dx(:) + g*dx(:);
        end

        function [c,ceq] = fminconNlc(dx,TR)
            ceq = [];
            c = sum(dx.^2)-(TR)^2;
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

    end

end