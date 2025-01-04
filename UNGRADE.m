classdef UNGRADE < UNLSI

    %%%%%%%%%%%%UNGRDE%%%%%%%%%%%%%%%%%%%
    %The MIT License Copyright (c) 2023 Naoto Morita
    %以下に定める条件に従い、本ソフトウェアおよび関連文書のファイル（以下「ソフトウェア」）の複製を取得するすべての人に対し、ソフトウェアを無制限に扱うことを無償で許可します。
    %これには、ソフトウェアの複製を使用、複写、変更、結合、掲載、頒布、サブライセンス、および/または販売する権利、およびソフトウェアを提供する相手に同じことを許可する権利も無制限に含まれます。
    %上記の著作権表示および本許諾表示を、ソフトウェアのすべての複製または重要な部分に記載するものとします。
    %ソフトウェアは「現状のまま」で、明示であるか暗黙であるかを問わず、何らの保証もなく提供されます。
    %ここでいう保証とは、商品性、特定の目的への適合性、および権利非侵害についての保証も含みますが、それに限定されるものではありません。
    %作者または著作権者は、契約行為、不法行為、またはそれ以外であろうと、ソフトウェアに起因または関連し、あるいはソフトウェアの使用またはその他の扱いによって生じる一切の請求、損害、その他の義務について何らの責任も負わないものとします。

    %また、本ソフトウェアに付帯するソフトウェアとしてOpenVSPおよびdistanceVertex2Mesh.mを用いています。
    %上記ソフトウェアのライセンス表記を以下に示します。
    %openVSP
    %Copyright (c) 2012 United States Government as represented by the Administrator for The National Aeronautics and Space Administration. All Rights Reserved.
    %
    %DISTANCEVERTEX2MESH - calculate the distance between vertices and a mesh
    %Author: Christopher Haccius
    %Telecommunications Lab, Saarland University, Germany
    %email: haccius@nt.uni-saarland.de
    %March 2015; Last revision: 26-March-2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        orgMesh %解析を行うメッシュ
        orgGeom %メッシュ変形に用いる基準形状
        unscaledVar %スケーリングされていない設計変数
        scaledVar %ub-lbでスケーリングした設計変数
        geomGenFun %基準形状を作成する関数
        meshGenFun %解析メッシュを作成する関数
        lb %設計変数の下限
        ub %設計変数の上限
        settingUNGRADE %最適化・解析の設定値
        geom2MeshMat %メッシュ変形行列
        phiRBF%メッシュ行列作成に使うRBF関数
        designScale %ub-lbのスケーリングパラメータ
        gradMesh %メッシュ節点の変化を設計変数で一次近似したもの
        gradSREF %基準面積の設計変数に対する一次近似
        gradBREF %横基準長の設計変数に対する一次近似
        gradCREF %縦基準長の設計変数に対する一次近似
        gradXYZREF %回転中心の設計変数に対する一次近似
        gradArginx %任意変数（重量など）の設計変数に対する一次近似
        argin_x %任意変数（重量やバッテリー搭載量などのユーザーで自由に設定する変数）
        iteration %設計更新を行った階数
        Hessian %準ニュートン法で近似されたヘッシアン
        BBMatrix
        history %設計更新の履歴
        flowNoList %マッハ数とflowNoの関連付け行列
        LagrangianInfo
        eqnSol %構造解析の解を一時的に保存
    end

    methods(Access = public)

        function obj = UNGRADE(meshGenFun,geomGenFun,unscaledVariables,lb,ub,halfmesh,vertsMergeTol)
            %%%%%%%%%%%%メッシュの登録と基準サーフェス生成関数の登録%%%%%%%%%%%%%%%
            %orgVerts,orgCon　実際に解析を行うメッシュ（openVSPのCFDツール等で生成した（基本的に）オーバーラップのない非構造メッシュ）
            %meshGenFun　設計変数（スケールされていない）から基準サーフェス（解析メッシュと同じ表面を持つ、トリムされていないメッシュ）
            %designVariables : 設計変数
            %lb：設計変数の下限
            %ub：設計変数の上限
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            [orgMeshVerts, orgMeshCon,surfID,wakeLineID, desOrg] = meshGenFun(unscaledVariables(:)');
            if ~exist('vertsMergeTol', 'var')
                vertsTol = 0.001;
            else
                vertsTol = vertsMergeTol;
            end
            obj = obj@UNLSI(orgMeshVerts,orgMeshCon,surfID,wakeLineID,halfmesh,vertsTol);
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            if any(obj.designScale<0)
                error("check the ub and lb setting");
            end

            obj.geomGenFun = geomGenFun;
            obj.meshGenFun = meshGenFun;


            %%%%%%%%%%%オプションのデフォルト設定
            % gradientCalcMethodex
            % 'direct'-->直接パネル法行列を計算し差分を作成
            % 'chain'-->各節点位置の変化に対する近似パネル法行列と、設計変数変化に対する節点位置の変化を用いてチェインルールで勾配計算する。
            % 'nonlin'-->approxmatによる近似パネル法行列を使った直接最適化
            %
            obj.settingUNGRADE.nMemory = 20; %記憶制限準ニュートン法のさかのぼり回数（厳密な記憶制限法は実装していない）
            obj.settingUNGRADE.H0 = eye(numel(obj.lb)); %初期ヘッシアン
            obj.settingUNGRADE.meshGradientPerturbation = 0.001;%メッシュ勾配をとるときのスケールされてない設計変数の摂動量
            obj.settingUNGRADE.gradientCalcMethod = "chain"; %設計変数偏微分の取得方法："direct", "chain", "nonlin"
            obj.settingUNGRADE.HessianUpdateMethod = "BFGS"; %ヘッシアン更新は "BFGS","DFP","Broyden","SR1"から選択
            obj.settingUNGRADE.updateMethod = "Quasi-Newton";%解を更新する方法 Quasi-Newton, Steepest-Descent,Levenberg–Marquardt,
            obj.settingUNGRADE.betaLM = 0.1; %レーベンバーグマッカートの重み係数
            obj.settingUNGRADE.alphaSD = 0.5; %最急降下法の重み
            obj.settingUNGRADE.TrustRegion = 0.1; %設計更新を行う際のスケーリングされた設計変数更新量の最大値
            obj.settingUNGRADE.dynCoefFlag = 0;
            obj.settingUNGRADE.r0RBF = 1;
            obj.settingUNGRADE.femCouplingFlag = 2; %0:カップリングしない 1:弱錬成（変形後の空力を考慮しない） 2:強連成（変形後の空力を考慮する）
            obj.settingUNGRADE.strongCouplingTol = 1e-10;
            obj.settingUNGRADE.selfLoadFlag = 1;
            %%%%%%%%%%%%%
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol,"delete");
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(desOrg);
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF,optXYZREF);
            fprintf("Reference Value\nSREF : %f BREF : %f CREF : %f\nXREF : %f YREF : %f ZREF : %f\n",optSREF,optBREF,optCREF,optXYZREF(1),optXYZREF(2),optXYZREF(3));
            disp("Argin Value");
            disp(obj.argin_x);
            obj.unscaledVar = desOrg(:)';
            obj.scaledVar = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);

            %Geom to Meshのメッシュ変形行列を作成
            %{
            Rgg = obj.calcRMat(obj.orgGeom.Points,obj.orgGeom.Points);
            obj.phiRBF = @(r)obj.phi3(r,obj.settingUNLSI.r0RBF);
            Mgg = obj.phiRBF(Rgg);
            Pgg = [ones(size(obj.orgGeom.Points,1),1),obj.orgGeom.Points]';
            Rmg = obj.calcRMat(obj.orgMesh.Points,obj.orgGeom.Points);
            Amg = [ones(size(obj.orgMesh.Points,1),1),obj.orgMesh.Points,obj.phiRBF(Rmg)];
            invM = pinv(Mgg);
            Mp = pinv(Pgg*invM*Pgg');
            obj.geom2MeshMat = Amg * [Mp*Pgg*invM;invM-invM*Pgg'*Mp*Pgg*invM];
            %}
            obj.phiRBF = @(r)obj.phi3(r,obj.settingUNGRADE.r0RBF);
            meshScaling = abs(max(obj.orgGeom.Points,[],1)-min(obj.orgGeom.Points,[],1));
            meshScaling(meshScaling==0) = 1;
            Agg = obj.phiRBF(obj.calcRMat(obj.orgGeom.Points./meshScaling,obj.orgGeom.Points./meshScaling));
            Amg = obj.phiRBF(obj.calcRMat(obj.orgMesh.Points./meshScaling,obj.orgGeom.Points./meshScaling));
            obj.geom2MeshMat = Amg*pinv(Agg);
            %}

            obj.iteration = 0;
            fprintf("iteration No. %d ---> Variables:\n",obj.iteration)
            disp(desOrg);
            obj.LagrangianInfo.Lorg = [];
            obj.LagrangianInfo.dL_dx = [];
            obj.LagrangianInfo.alin = [];
            obj.LagrangianInfo.blin = [];

            obj.Hessian = obj.settingUNGRADE.H0;
            obj.BBMatrix = eye(numel(obj.lb));

            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        end

        function obj = setFlowCondition(obj,alpha,beta,Mach,Re,dynPress)
            if nargin<6
                obj.settingUNGRADE.Mach = Mach;
                obj.settingUNGRADE.alpha = alpha;
                obj.settingUNGRADE.beta = beta;
                obj.settingUNGRADE.Re = Re;
                if any([numel(obj.settingUNGRADE.Mach),numel(obj.settingUNGRADE.alpha),numel(obj.settingUNGRADE.beta)]~=numel(obj.settingUNGRADE.Mach))
                    error("No. of case is not match");
                end
            else
                obj.settingUNGRADE.Mach = Mach;
                obj.settingUNGRADE.alpha = alpha;
                obj.settingUNGRADE.beta = beta;
                obj.settingUNGRADE.Re = Re;
                obj.settingUNGRADE.dynPress = dynPress;
                if any([numel(obj.settingUNGRADE.Mach),numel(obj.settingUNGRADE.alpha),numel(obj.settingUNGRADE.beta)]~=numel(obj.settingUNGRADE.dynPress))
                    error("No. of case is not match");
                end
            end

            uniMach = unique(Mach);
            for i = 1:numel(uniMach)
                obj = obj.flowCondition(i,uniMach(i));
                obj = obj.setCf(i,obj.settingUNGRADE.Re);
            end
            for i = 1:numel(Mach)
                obj.flowNoList(i,1) = find(uniMach == Mach(i));
                obj.flowNoList(i,2) = uniMach(obj.flowNoList(i,1));
                obj.flowNoList(i,3) = obj.flowNoList(i,2)<1;
                obj.flowNoList(i,4) =  obj.settingUNGRADE.alpha(i);
                obj.flowNoList(i,5) =  obj.settingUNGRADE.beta(i);
                obj.flowNoList(i,6) =  obj.settingUNGRADE.Re;
                if nargin>5
                    obj.flowNoList(i,7) = obj.settingUNGRADE.dynPress(i);
                end
            end
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeCluster();
            end
        end

        function [con,verts,femID] = calcFemMeshfromVariables(obj,unscaledVariables)
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
            [orgMeshVerts,orgMeshCon] = obj.meshDeformation(obj,orgGeomVerts,orgGeomCon);
            verts = obj.femMeshDeformation(obj,orgMeshVerts,orgMeshCon);
            con = obj.femtri.ConnectivityList;
            femID = obj.femID;
        end

        function obj = updateMeshGeomfromVariables(obj,unscaledVariables,meshDeformFlag)
            %%%%%%%設計変数からメッシュや解析条件を更新する%%%%%%%%%%%%%%
            %unscaledVariables : スケーリングされていない設計変数
            %meshDeformFlag : 1:新しくメッシュを作らずにmeshdeformationによって変形したものを新しいメッシュとする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            %if nargin<3
            %    meshDeformFlag = 0;
            %end
            if meshDeformFlag == 1
                [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
                [orgMeshVerts,orgMeshCon] = obj.meshDeformation(obj,orgGeomVerts,orgGeomCon);
                unscaledVariables = desOrg(:)';
                surfID = obj.surfID;
                wakelineID = obj.wakelineID;
                fprintf("iteration No. %d set by Mesh Deformation ---> Variables:\n",obj.iteration)
                disp(desOrg);
                obj = obj.setMesh(orgMeshVerts,orgMeshCon,surfID,wakelineID);
                obj.history.setMeshbyDeform(obj.iteration) = 1;
            else
                [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
                [orgMeshVerts, orgMeshCon,surfID,wakelineID, unscaledVariables] = obj.meshGenFun(desOrg(:)');
                fprintf("iteration No. %d set by Mesh Generation ---> Variables:\n",obj.iteration)
                obj = obj.setMesh(orgMeshVerts,orgMeshCon,surfID,wakelineID);
                disp(unscaledVariables);
                obj.history.setMeshbyDeform(obj.iteration) = 0;
            end
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol,"delete");
            fprintf("Reference Value\nSREF : %f BREF : %f CREF : %f\nXREF : %f YREF : %f ZREF : %f\n",optSREF,optBREF,optCREF,optXYZREF(1),optXYZREF(2),optXYZREF(3));
            disp("Argin Value");
            disp(obj.argin_x);
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF,optXYZREF);
            obj.unscaledVar = unscaledVariables(:)';
            obj.scaledVar = (unscaledVariables(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);

            %Geom to Meshのメッシュ変形行列を作成
            %{
            Rgg = obj.calcRMat(obj.orgGeom.Points,obj.orgGeom.Points);
            obj.phiRBF = @(r)obj.phi3(r,obj.settingUNLSI.r0RBF);
            Mgg = obj.phiRBF(Rgg);
            Pgg = [ones(size(obj.orgGeom.Points,1),1),obj.orgGeom.Points]';
            Rmg = obj.calcRMat(obj.orgMesh.Points,obj.orgGeom.Points);
            Amg = [ones(size(obj.orgMesh.Points,1),1),obj.orgMesh.Points,obj.phiRBF(Rmg)];
            invM = pinv(Mgg);
            Mp = pinv(Pgg*invM*Pgg');
            obj.geom2MeshMat = Amg * [Mp*Pgg*invM;invM-invM*Pgg'*Mp*Pgg*invM];
            %}
            meshScaling = max(abs(obj.orgGeom.Points),[],1)-min(abs(obj.orgGeom.Points),[],1);
            meshScaling(meshScaling==0) = 1;
            Agg = obj.phiRBF(obj.calcRMat(obj.orgGeom.Points./meshScaling,obj.orgGeom.Points./meshScaling));
            Amg = obj.phiRBF(obj.calcRMat(obj.orgMesh.Points./meshScaling,obj.orgGeom.Points./meshScaling));
            obj.geom2MeshMat = Amg*pinv(Agg);
            %}

            obj.LagrangianInfo.Lorg = [];
            obj.LagrangianInfo.dL_dx = [];
            obj.LagrangianInfo.alin = [];
            obj.LagrangianInfo.blin = [];
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeCluster();
            end
            for i = 1:numel(obj.flow)
                obj = obj.setCf(i,obj.settingUNGRADE.Re);
            end
            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        end

        function obj = solveAnalysis(obj,flowNo,alpha,beta)
            %%%%%%%%%%%現在の機体形状に対する解析の実行%%%%%
            %flowNo : UNLSIのflowNo
            %alpha : 解析する迎角
            %beta : 解析する横滑り角
            %omega : 回転角速度（任意）
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(obj.LHS)
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.makeEquation();
                end
            end
            if obj.settingUNGRADE.femCouplingFlag > 9  && isempty(obj.femLHS)
                obj = obj.makeFemEquation();
            end
            u0 = obj.solvePertPotential(flowNo,alpha,beta);
            [~,~,~,~,obj] = obj.solveFlowForAdjoint(u0,flowNo,alpha,beta);
            if obj.settingUNGRADE.dynCoefFlag == 1
                [udwf,udwr] = obj.calcDynCoefdu(flowNo,alpha,beta);
                [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,flowNo,alpha,beta);
            end
            if obj.settingUNGRADE.femCouplingFlag > 0
                %構造カップリングの求解
                for i = 1:numel(obj.Cp)
                    [deltabuff,deltadot,~] = obj.solveFem(obj.Cp{i}.*obj.settingUNGRADE.dynPress(1),obj.settingUNGRADE.selfLoadFlag);
                    for j = 1:numel(deltabuff)
                        delta{i,j} = deltabuff{j};
                    end
                end
                obj.eqnSol.delta = delta;
                obj.eqnSol.deltadot = deltadot;
            end
        end

        function [I,con,obj] = evaluateObjFun(obj,objandConsFun)
            %%%%%%%%%%%現在の機体形状に対する評価関数と制約条件の評価%%%%%
            %objandConsFun : 評価関数と制約条件を出す関数。calcNextVarとインプットは同じ

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(obj.LHS)
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.makeEquation();
                end
            end
            if obj.settingUNGRADE.femCouplingFlag > 0 && isempty(obj.femLHS)
                obj = obj.makeFemEquation();
            end
            obj.approximated = 0;
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            delta0 = {};
            u0 = [];
            for iter = 1:numel(obj.flow)
                alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                u0buff = obj.solvePertPotential(iter,alphabuff,betabuff);
                u0 = [u0;u0buff];
                [~,~,~,~,obj] = obj.solveFlowForAdjoint(u0buff,iter,alphabuff,betabuff);
                if obj.settingUNGRADE.dynCoefFlag == 1
                    [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff,betabuff);
                    [~,obj] =  obj.calcDynCoefforAdjoint(u0buff,udwf,udwr,iter,alphabuff,betabuff);
                end
                if obj.settingUNGRADE.femCouplingFlag > 0
                    [deltabuff,deltadot,~] = obj.solveFem(obj.Cp{iter}.*obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach),obj.settingUNGRADE.selfLoadFlag);
                    delta0 = [delta0,deltabuff];
                    obj.eqnSol.delta = delta0;
                    obj.eqnSol.deltadot = deltadot;
                end
            end
            if obj.settingUNGRADE.femCouplingFlag == 0
                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,obj.AERODATA,obj.DYNCOEF,obj.Cp,obj.Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,[]);
            else
                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,obj.AERODATA,obj.DYNCOEF,obj.Cp,obj.Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,delta0);
            end
            [I,con] = objandConsFun(geom,aero,fem);

            %             obj = obj.calcApproximatedEquation(1);
            %             nfem = size(obj.femLHS,1);
            %             for iter = 1:10
            %                 [jac,f] = obj.calcStrongCouplingJaccobian(u0,delta0);
            %                 disp([norm(f),max(abs(f))]);
            %                 duddelta = -((jac'*jac+0.1.*eye(size(jac)))\jac'*f);
            %                 u0 = u0 + full(duddelta(1:size(u0,1),1));
            %                 deltap0 = duddelta(size(u0,1)+1:end,1);
            %                 for i = 1:size(obj.flowNoList,1)
            %                     disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
            %                     disp_buff(obj.femutils.InvMatIndex,1) = deltap0(nfem*(i-1)+1:nfem*i);
            %                     ddelta{i} = zeros(size(obj.femtri.Points,1),6);
            %                     ddelta{i}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
            %                     ddelta{i}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
            %                     ddelta{i}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
            %                     ddelta{i}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
            %                     ddelta{i}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
            %                     ddelta{i}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
            %                     delta0{i} = delta0{i} +  ddelta{i};
            %                 end
            %             end

        end

        function [f,jac] = calcStrongCouplingJaccobian(obj,u,delta)
            %deltaから機体メッシュを再生成
            nflowVar = size(obj.LHS,1);
            dS_du = [];
            dR_ddelta = [];
            calcCount = 1;
            R0all = [];
            S0all = [];
            for iter = 1:numel(obj.flow)
                alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                for jter = 1:numel(alphabuff)
                    modVerts = obj.calcModifiedVerts(delta{calcCount});
                    obj2 = obj.makeApproximatedInstance(modVerts);
                    u0 = u(nflowVar*(calcCount-1)+1:nflowVar*calcCount,1);
                    [AERODATA,Cp,Cfe,R0,obj2] = obj2.solveFlowForAdjoint(u0,iter,alphabuff(jter),betabuff(jter));
                    R0all = [R0all;R0];
                    obj2.plotGeometry(3,Cp{iter},[-2,1]);
                    %[deltaOut,~,S0] = obj.solveFem(Cp{i}(:,j).*obj.settingUNGRADE.dynPress(1),0);
                    %delta{i,j} = deltaOut{1};
                    deltaIn{1} = delta{calcCount};
                    [~,~,S0] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                    S0all = [S0all;S0];
                    if nargout>1
                        %u微分
                        tic;
                        pert = eps^(1/3);
                        for k = 1:numel(u0)
                            uf = u0;
                            uf(k) = uf(k)+pert;
                            [AERODATA,Cp,Cfe,Rf,~] = obj2.solveFlowForAdjoint(uf,iter,alphabuff(jter),betabuff(jter));
                            deltaIn{1} = delta{calcCount};
                            [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                            dS_du_buff(:,k) = (Sf-S0)/pert;
                        end
                        if calcCount == 1
                            dS_du = dS_du_buff;
                        else
                            dS_du = blkdiag(dS_du,dS_du_buff);
                        end
                        toc;
                        tic;
                        %delta微分の計算
                        %dR_ddelta
                        deltap0_buff = delta{calcCount}(:);
                        deltap0 = deltap0_buff(obj.femutils.MatIndex==1,1);
                        nfem = size(obj.femLHS,1);
                        for k = 1:size(deltap0,1)
                            if k < sum(obj.femutils.MatIndex(1:sub2ind(size(delta{calcCount}),1,4)))
                                deltapf = deltap0;
                                deltapf(k,1) = deltap0(k,1)+pert;
                                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                                disp_buff(obj.femutils.InvMatIndex,1)=deltapf;
                                deltaf = delta;
                                deltaf{calcCount} = zeros(size(obj.femtri.Points,1),6);
                                deltaf{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                                modVerts = obj.calcModifiedVerts(deltaf{calcCount});
                                obj3 = obj2.makeApproximatedInstance(modVerts,1);
                                [AERODATA,Cp,Cfe,Rf,~] = obj3.solveFlowForAdjoint(u0,iter,alphabuff(jter),betabuff(jter));
                                deltaIn{1} = deltaf{calcCount};
                                [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                                dR_ddelta_buff(:,k) = (Rf-R0)./pert;
                                dS_ddelta_buff(:,k) = (Sf-S0)./pert;
                            else
                                dR_ddelta_buff(:,k) = 0;
                                dS_ddelta_buff(:,k) = obj.femLHS(:,k);
                            end
                        end
                        toc;
                        if calcCount == 1
                            dR_ddelta = dR_ddelta_buff;
                            dS_ddelta = dS_ddelta_buff;
                        else
                            dR_ddelta = blkdiag(dR_ddelta,dR_ddelta_buff);
                            dS_ddelta = blkdiag(dS_ddelta,dS_ddelta_buff);
                        end
                    end
                    calcCount = calcCount + 1;
                end
            end

            for i = 1:size(obj.flowNoList,1)
                if i == 1
                    if obj.flowNoList(i,3) == 1
                        dR_du = (obj2.LHS+obj2.wakeLHS); %亜音速
                    else
                        dR_du = eye(nbPanel);
                    end
                else
                    if obj.flowNoList(i,3) == 1
                        dR_du = blkdiag(dR_du,(obj2.LHS+obj2.wakeLHS));
                    else
                        nbPanel = sum(obj.paneltype == 1);
                        dR_du = blkdiag(dR_du,eye(nbPanel));
                    end
                end
            end
            if nargout>1
                jac= [dR_du,dR_ddelta;dS_du./10000,dS_ddelta./10000];
            end
            f = [R0all;S0all./10000];

        end

        function viewMesh(obj,fig,modMesh,~)
            %%%%%%%%%%%Meshの表示%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %modMesh:表示するメッシュ
            %fig : 描画するfigure
            %~ : 適当な引数を入れるとメッシュ変化量を表示する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                modMesh = obj.orgMesh.Points;
            end
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

        function viewGeom(obj,fig,modGeom,~)
            %%%%%%%%%%%Geomの表示%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %modMesh:表示するメッシュ
            %fig : 描画するfigure
            %~ : 適当な引数を入れるとメッシュ変化量を表示する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                modGeom = obj.orgGeom.Points;
            end
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
        %{
        function [errMax,errMin] = calcErrMeshGeom(obj)
            [errMax,errMin] = obj.dist2geom(obj.orgGeom,obj.orgMesh);
        end
        %}
        function checkGeomGenWork(obj,pert,checkVar,fig)
            %%%%%%%%%%%設計変数が機能しているかチェックする%%%%%%%%%%%
            %最後の引数によって
            %メッシュの変化量を色に表示するか動かを設定できる
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = numel(obj.scaledVar);
            if nargin == 2
                checkVar = 1:ndim;
            end
            randparam = rand(1).*ones(1,ndim);
            randDes = randparam.*obj.designScale+obj.lb;
            modSurforg = obj.geomGenFun(randDes);
            if nargin == 4
                viewtri = triangulation(obj.orgGeom.ConnectivityList,modSurforg);
                figure(fig);clf
                trisurf(viewtri);
                axis equal;drawnow();
                fprintf("randamized parameter : ");
                disp(randDes);
                pause(1);
            end
            for i = checkVar
                sampleDes = randparam.*obj.designScale+obj.lb;
                sampleDes(i) = (randparam(i) + pert).*obj.designScale(i)+obj.lb(i);
                modSurf = obj.geomGenFun(sampleDes);
                if all(abs(modSurforg(:)-modSurf(:)) <sqrt(eps))
                    error("Variables No. %d : NOT MOVED",i);
                end
                if nargin == 4
                    viewtri = triangulation(obj.orgGeom.ConnectivityList,modSurf);
                    figure(fig);clf
                    trisurf(viewtri);
                    axis equal;drawnow();
                    pause(1);
                end
            end
        end


        function obj = setUNGRADESettings(obj,varargin)
            %%%%%%%%%設定を変更する%%%%%%%%%%%%
            %varargin : 変数名と値のペアで入力
            %obj.settingの構造体の内部を変更する
            %構造体の内部変数の名前と一致していないとエラーが出るようになっている
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mod(numel(varargin),1) == 1
                error("input style not match");
            end
            for iter = 1:numel(varargin)/2
                if isfield(obj.settingUNGRADE,varargin{2*iter-1})
                    obj.settingUNGRADE = setfield(obj.settingUNGRADE,varargin{2*iter-1},varargin{2*iter});
                    if strcmpi(varargin{2*iter-1},"Mach")
                        if numel(obj.settingUNGRADE.Mach) ~= numel(obj.settingUNGRADE.alpha) || numel(obj.settingUNGRADE.Mach) ~= numel(obj.settingUNGRADE.beta) || numel(obj.settingUNGRADE.beta) ~= numel(obj.settingUNGRADE.alpha)
                            error("No. of case is not match");
                        end
                        uniMach = unique(obj.settingUNGRADE.Mach);
                        for i = 1:numel(uniMach)
                            obj = obj.flowCondition(i,uniMach(i));
                            obj = obj.setCf(i,obj.settingUNGRADE.Re,obj.settingUNGRADE.Lch,obj.settingUNGRADE.k,obj.settingUNGRADE.LTratio,obj.settingUNGRADE.coefficient);
                        end
                        for i = 1:numel(obj.settingUNGRADE.Mach)
                            obj.flowNoList(i,1) = find(uniMach == obj.settingUNGRADE.Mach(i));
                            obj.flowNoList(i,2) = uniMach(obj.flowNoList(i,1));
                            obj.flowNoList(i,3) = obj.flowNoList(i,2)<1;
                        end
                    end
                else
                    error("Field name is not match");
                end
            end
            if obj.iteration == 0
                obj.Hessian = obj.settingUNGRADE.H0;
            end
        end

        function [obj] = calcLagrangianGradient(obj,objandConsFun,cmin,cmax,varargin)
            %%%%%%%%%%%%指定した評価関数と制約条件における設計変数勾配の計算%%%%%%%%
            %objandConFun : [res con] 評価関数値と制約条件値を計算する関数
            %cmin,cmax : 制約の上下限値
            %varargin : 追加するとsetOptionsをその変数で実行してくれる
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(cmax-cmin<0)
                error("check the cmin and cmax setting");
            end
            cScale = cmax-cmin;        
            if nargin>4
                obj = obj.setUNGRADESettings(varargin{:});
            end
            obj.iteration = obj.iteration + 1;
            fprintf("iteration No. %d : Gradient Caluculation Started\n -- GradientCalcMethod : %s\n",obj.iteration,obj.settingUNGRADE.gradientCalcMethod);
            %初期点の解析
            nbPanel = sum(obj.paneltype == 1);
            if isempty(obj.LHS)
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.makeEquation();
                end
            end
            if obj.settingUNGRADE.femCouplingFlag > 0
                if isempty(obj.femLHS)
                    obj = obj.makeFemEquation();
                end
            end
            if obj.settingUNGRADE.femCouplingFlag == 2 || strcmpi(obj.settingUNGRADE.gradientCalcMethod,"chain")
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.calcApproximatedEquation(1);
                end      
                if obj.settingUNGRADE.femCouplingFlag > 0
                    obj = obj.calcApproximatedFemEquation();
                end
            end
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            u0 = zeros(nbPanel*size(obj.flowNoList,1),1);
            R0 = zeros(nbPanel*size(obj.flowNoList,1),1);
            delta0 = {};
            S0 = [];
            for iter = 1:numel(obj.flow)
                alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                [u0solve,~] = obj.solvePertPotential(iter,alphabuff,betabuff);%ポテンシャルを求める
                if obj.settingUNGRADE.dynCoefFlag == 1
                    [udwf{iter},udwr{iter}] = obj.calcDynCoefdu(iter,alphabuff,betabuff);
                end
                [AERODATA0,Cp0,Cfe0,R0solve,obj] = obj.solveFlowForAdjoint(u0solve,iter,alphabuff,betabuff);%ポテンシャルから空力係数を計算
                if obj.settingUNGRADE.dynCoefFlag == 1
                    [DYNCOEF0,obj] =  obj.calcDynCoefforAdjoint(u0solve,udwf{iter},udwr{iter},iter,alphabuff,betabuff);
                else
                    DYNCOEF0 = [];
                end
                %結果をマッピング
                lktable = find(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                for i = 1:numel(lktable)
                    u0((lktable(i)-1)*nbPanel+1:lktable(i)*nbPanel,1) = u0solve(nbPanel*(i-1)+1:nbPanel*i,1);
                    R0((lktable(i)-1)*nbPanel+1:lktable(i)*nbPanel,1) = R0solve(nbPanel*(i-1)+1:nbPanel*i,1);
                end
                if obj.settingUNGRADE.femCouplingFlag > 0 %構造カップリングの求解
                    [deltabuff,deltadot,S0solve] = obj.solveFem(Cp0{iter}.*obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach),obj.settingUNGRADE.selfLoadFlag);
                    S0 = [S0;S0solve];
                    delta0 = [delta0,deltabuff];
                end
            end
            if obj.settingUNGRADE.femCouplingFlag > 1
                %連成解を求める
                nfem = size(obj.femLHS,1);
                for iter = 1:3
                    [f,jac] = obj.calcStrongCouplingJaccobian(u0,delta0);
                    jacupdateflag = 1;
                    for j = 1:20
                        duddelta = -(lsqminnorm(jac'*jac,jac',obj.settingUNLSI.lsqminnormTol)*f);
                        u1 = u0 + full(duddelta(1:size(u0,1),1));
                        deltap1 = duddelta(size(u0,1)+1:end,1);
                        for i = 1:size(obj.flowNoList,1)
                            disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                            disp_buff(obj.femutils.InvMatIndex,1) = deltap1(nfem*(i-1)+1:nfem*i);
                            ddelta{i} = zeros(size(obj.femtri.Points,1),6);
                            ddelta{i}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                            ddelta{i}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                            ddelta{i}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                            ddelta{i}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                            ddelta{i}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                            ddelta{i}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                            delta1{i} = delta0{i} +  ddelta{i};
                        end
                        fnext = obj.calcStrongCouplingJaccobian(u1,delta1);
                        disp(full([norm(f),max(abs(f))]))
                        if norm(fnext)<norm(f) && max(abs(fnext))<max(abs(f))
                            jacupdateflag = 0;
                            f = fnext;
                            u0 = u1;
                            delta0 = delta1;
                            if norm(f)<obj.settingUNGRADE.strongCouplingTol ||  max(abs(f))<obj.settingUNGRADE.strongCouplingTol
                                break;
                            end
                        else
                            break;
                        end
                    end
                    if norm(f)<obj.settingUNGRADE.strongCouplingTol ||  max(abs(f))<obj.settingUNGRADE.strongCouplingTol || jacupdateflag == 1
                        break;
                    end     
                end
                disp("Strong Coupling Residial : (norm) (max)")
                disp(full([norm(f),max(abs(f))]))
                calcCount = 1;
                S0 = [];
                R0 = [];
                for iter = 1:numel(obj.flow)
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    for jter = 1:numel(alphabuff)
                        usolve = u0((calcCount-1)*nbPanel+1:calcCount*nbPanel,1);
                        modVerts = obj.calcModifiedVerts(delta0{calcCount});
                        %aeroObj{calcCount} = obj.makeApproximatedInstance(modVerts);
                        aeroObj{calcCount} = obj.setVerts(modVerts);
                        aeroObj{calcCount} = aeroObj{calcCount}.makeEquation();
                        aeroObj{calcCount} = aeroObj{calcCount}.calcApproximatedEquation(1);

                        if obj.settingUNGRADE.dynCoefFlag == 1
                            [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff(jter),betabuff(jter));
                        end
                        [AERODATAbuff,Cpbuff,Cfebuff,Rsolve,aeroObj{calcCount}] = aeroObj{calcCount}.solveFlowForAdjoint(usolve,iter,alphabuff(jter),betabuff(jter));
                        AERODATA0{iter}(jter,:) = AERODATAbuff{iter};
                        Cp0{iter}(:,jter) = Cpbuff{iter};
                        Cfe0{iter} = Cfebuff{iter};
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                            deltaIn{1} = delta0{calcCount};
                            [~,~,Ssolve] = obj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                        end
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            DYNCOEFbuff =  aeroObj{calcCount}.calcDynCoefforAdjoint(usolve,udwf,udwr,iter,alphabuff(jter),betabuff(jter));
                            DYNCOEF0{iter,jter} = DYNCOEFbuff{iter,1};
                        else
                            DYNCOEF = [];
                        end
                        S0 = [S0;Ssolve];
                        R0 = [R0;Rsolve];
                        calcCount = calcCount + 1;
                    end
                end
            else
                for iter = 1:numel(obj.flow)
                    calcCount = 1;
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    for jter = 1:numel(alphabuff)
                        aeroObj{calcCount} = obj;
                        calcCount = calcCount + 1;
                    end
                end
            end

            if obj.settingUNGRADE.femCouplingFlag > 0
                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA0,DYNCOEF0,Cp0,Cfe0,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,delta0);
            else
                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA0,DYNCOEF0,Cp0,Cfe0,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,[]);
            end

            [I0,con0] = objandConsFun(geom,aero,fem);
            con0 = (con0-cmin)./cScale;
            fprintf("Orginal Objective Value and Constraints:\n")
            disp([I0,((con0.*cScale)+cmin)']);
            obj.history.objVal(obj.iteration) = I0;
            obj.history.conVal(:,obj.iteration) = con0(:);
            fprintf("AERODATA of iteration No.%d ->\n",obj.iteration);
            disp(vertcat(AERODATA0{:}));
            if obj.settingUNGRADE.dynCoefFlag == 1
                fprintf("DYNCOEF of iteration No.%d ->\n",obj.iteration);
                disp(vertcat(DYNCOEF0{:}));
            end
            %%%
            %明示・非明示随伴方程式法の実装
            %u微分の計算
            invR.center = obj.RbfinvRMake(obj,obj.center,1,0.01);
            invR.verts = obj.RbfinvRMake(obj,obj.tri.Points,1,0.01);
            pert = eps^(1/3);
            calcCount = 1;
            dI_du = [];
            dcon_du = [];
            dS_du = [];
            for iter = 1:numel(obj.flow)
                alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                for jter = 1:numel(alphabuff)
                    usolve = u0((calcCount-1)*nbPanel+1:calcCount*nbPanel,1);
                    if obj.settingUNGRADE.dynCoefFlag == 1
                        [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff(jter),betabuff(jter));
                    end
                    for i = 1:numel(usolve)
                        uf = usolve;
                        ur = usolve;
                        uf(i) = uf(i)+pert;
                        ur(i) = ur(i)-pert;
                        [AERODATAbuff,Cpbuff,Cfebuff,~,aeroObj{calcCount}] = aeroObj{calcCount}.solveFlowForAdjoint(uf,iter,alphabuff(jter),betabuff(jter));
                        AERODATA = AERODATA0;
                        AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                        Cp = Cp0;
                        Cp{iter}(:,jter) = Cpbuff{iter};
                        Cfe = Cfe0;
                        Cfe{iter} = Cfebuff{iter};
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                            deltaIn{1} = delta0{calcCount};
                            [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                        end
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            DYNCOEFbuff =  aeroObj{calcCount}.calcDynCoefforAdjoint(uf,udwf,udwr,iter,alphabuff(jter),betabuff(jter));
                            DYNCOEF = DYNCOEF0;
                            DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                        else
                            DYNCOEF = [];
                        end
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,delta0,invR);
                        else
                            [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,[],invR);
                        end
                        [If,conf] = objandConsFun(geom,aero,fem);
                        conf = (conf-cmin)./cScale;

                        [AERODATAbuff,Cpbuff,Cfebuff,~,aeroObj{calcCount}] = aeroObj{calcCount}.solveFlowForAdjoint(ur,iter,alphabuff(jter),betabuff(jter));
                        AERODATA = AERODATA0;
                        AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                        Cp = Cp0;
                        Cp{iter}(:,jter) = Cpbuff{iter};
                        Cfe = Cfe0;
                        Cfe{iter} = Cfebuff{iter};
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                            deltaIn{1} = delta0{calcCount};
                            [~,~,Sr] = obj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                        end
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            DYNCOEFbuff =  aeroObj{calcCount}.calcDynCoefforAdjoint(ur,udwf,udwr,iter,alphabuff(jter),betabuff(jter));
                            DYNCOEF = DYNCOEF0;
                            DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                        else
                            DYNCOEF = [];
                        end
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,delta0,invR);
                        else
                            [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,[],invR);
                        end
                        [Ir,conr] = objandConsFun(geom,aero,fem);
                        conr = (conr-cmin)./cScale;

                        dI_du_buff(i) = (If-Ir)/pert/2;%評価関数のポテンシャルに関する偏微分
                        if not(isempty(con0))
                            dcon_du_buff(:,i) = (conf-conr)/pert/2;
                        end
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            dS_du_buff(:,i) = (Sf-Sr)/pert/2;
                        end
                    end
                    dI_du = [dI_du,dI_du_buff];
                    dcon_du = [dcon_du,dcon_du_buff];
                    if obj.settingUNGRADE.femCouplingFlag > 0
                        dS_du = blkdiag(dS_du,dS_du_buff);
                    end
                    if calcCount == 1
                        if obj.flowNoList(calcCount,3) == 1
                            dR_du = (aeroObj{calcCount}.LHS+aeroObj{calcCount}.wakeLHS); %亜音速
                        else
                            dR_du = eye(nbPanel);
                        end
                    else
                        if obj.flowNoList(calcCount,3) == 1
                            dR_du = blkdiag(dR_du,(aeroObj{calcCount}.LHS+aeroObj{calcCount}.wakeLHS));
                        else
                            dR_du = blkdiag(dR_du,eye(nbPanel));
                        end
                    end
                    calcCount = calcCount + 1;
                end
            end


            %メッシュの節点勾配を作成
            obj = obj.makeMeshGradient(obj,obj.settingUNGRADE.meshGradientPerturbation./obj.designScale);
            if obj.settingUNGRADE.femCouplingFlag > 0
                %delta微分の計算
                %dR_ddelta
                nfem = size(obj.femLHS,1);
                invR.center = obj.RbfinvRMake(obj,obj.center,1,0.01);
                invR.verts = obj.RbfinvRMake(obj,obj.tri.Points,1,0.01);
                calcCount = 1;
                dR_ddelta = [];
                dS_ddelta = [];
                dI_ddelta = [];
                dcon_ddelta = [];
                for iter = 1:numel(obj.flow)
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    for jter = 1:numel(alphabuff)
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff(jter),betabuff(jter));
                        end
                        deltap0_buff = delta0{calcCount}(:);
                        deltap0 = deltap0_buff(obj.femutils.MatIndex==1,1);
                        usolve = u0((calcCount-1)*nbPanel+1:calcCount*nbPanel,1);
                        for i = 1:size(deltap0,1)
                            if i < sum(obj.femutils.MatIndex(1:sub2ind(size(delta0{calcCount}),1,4)))
                                deltapf = deltap0;
                                deltapf(i,1) = deltap0(i,1)+pert;
                                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                                disp_buff(obj.femutils.InvMatIndex,1)=deltapf;
                                deltaf = delta0;
                                deltaf{calcCount} = zeros(size(obj.femtri.Points,1),6);
                                deltaf{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                                deltaf{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                                if obj.settingUNGRADE.femCouplingFlag > 1 %強連成
                                    modVerts = obj.calcModifiedVerts(deltaf{calcCount});
                                    aeroObj2 = aeroObj{calcCount}.makeApproximatedInstance(modVerts);
                                    deltaIn{1} = deltaf{calcCount};
                                    [AERODATAbuff,Cpbuff,Cfebuff,Rf,aeroObj2] = aeroObj2.solveFlowForAdjoint(usolve,iter,alphabuff(jter),betabuff(jter));
                                    %%%%%%%%%%%%%DYNCOEF
                                    if obj.settingUNGRADE.dynCoefFlag == 1
                                        DYNCOEFbuff =  aeroObj2.calcDynCoefforAdjoint(usolve,udwf,udwr,iter,alphabuff(jter),betabuff(jter));
                                        DYNCOEF = DYNCOEF0;
                                        DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                                    else
                                        DYNCOEF = [];
                                    end
                                    [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                                    AERODATA = AERODATA0;
                                    AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                                    Cp = Cp0;
                                    Cp{iter}(:,jter) = Cpbuff{iter};
                                    Cfe = Cfe0;
                                    Cfe{iter} = Cfebuff{iter};
                                else
                                    Cp = Cp0;
                                    Cfe = Cfe0;
                                    AERODATA = AERODATA0;
                                    DYNCOEF = DYNCOEF0;
                                end
                                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,deltaf,invR);
                                [If,conf] = objandConsFun(geom,aero,fem);
                                conf = (conf-cmin)./cScale;
                                clear aeroObj2

                                deltapr = deltap0;
                                deltapr(i,1) = deltap0(i,1)-pert;
                                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                                disp_buff(obj.femutils.InvMatIndex,1)=deltapr;
                                deltar = delta0;
                                deltar{calcCount} = zeros(size(obj.femtri.Points,1),6);
                                deltar{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                                deltar{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                                deltar{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                                deltar{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                                deltar{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                                deltar{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                                if obj.settingUNGRADE.femCouplingFlag > 1 %強連成
                                    modVerts = obj.calcModifiedVerts(deltar{calcCount});
                                    aeroObj2 = aeroObj{calcCount}.makeApproximatedInstance(modVerts);
                                    deltaIn{1} = deltar{calcCount};
                                    [AERODATAbuff,Cpbuff,Cfebuff,Rr,aeroObj2] = aeroObj2.solveFlowForAdjoint(usolve,iter,alphabuff(jter),betabuff(jter));
                                    [~,~,Sr] = obj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                                    if obj.settingUNGRADE.dynCoefFlag == 1
                                        DYNCOEFbuff =  aeroObj2.calcDynCoefforAdjoint(usolve,udwf,udwr,iter,alphabuff(jter),betabuff(jter));
                                        DYNCOEF = DYNCOEF0;
                                        DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                                    else
                                        DYNCOEF = [];
                                    end
                                    AERODATA = AERODATA0;
                                    AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                                    Cp = Cp0;
                                    Cp{iter}(:,jter) = Cpbuff{iter};
                                    Cfe = Cfe0;
                                    Cfe{iter} = Cfebuff{iter};
                                    dS_ddelta_buff(:,i) = (Sf-Sr)/pert/2;
                                    dR_ddelta_buff(:,i) = (Rf-Rr)/pert/2;
                                else
                                    Cp = Cp0;
                                    Cfe = Cfe0;
                                    AERODATA = AERODATA0;
                                    DYNCOEF = DYNCOEF0;
                                end
                                [geom,aero,fem] = obj.makeObjArguments(obj,desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x,deltar,invR);
                                [Ir,conr] = objandConsFun(geom,aero,fem);
                                conr = (conr-cmin)./cScale;
                                clear aeroObj2

                                dI_ddelta_buff(i) = (If-Ir)/pert/2;%評価関数のポテンシャルに関する偏微分
                                if not(isempty(con0))
                                    dcon_ddelta_buff(:,i) = (conf-conr)/pert/2;
                                end
                            else
                                dS_ddelta_buff(:,i) = obj.femLHS(:,i);
                                dR_ddelta_buff(:,i) = 0;
                                dI_ddelta_buff(i) = 0;
                                if not(isempty(con0))
                                    dcon_ddelta_buff(:,i) = 0;
                                end
                            end
                        end
                        if obj.settingUNGRADE.femCouplingFlag < 2 %強連成
                            dS_ddelta_buff = obj.femLHS;
                            dR_ddelta_buff = zeros(size(aeroObj{calcCount}.LHS,1),size(deltap0,1));
                        end
                        dR_ddelta = blkdiag(dR_ddelta,dR_ddelta_buff);
                        dS_ddelta = blkdiag(dS_ddelta,dS_ddelta_buff);
                        dI_ddelta = [dI_ddelta,dI_ddelta_buff];
                        if not(isempty(con0))
                            dcon_ddelta = [dcon_ddelta,dcon_ddelta_buff];
                        end
                        calcCount = calcCount + 1;
                    end

                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %x微分の計算
            pert = (eps)^(1/3);
            for i= 1:numel(obj.scaledVar)
                x = obj.scaledVar;
                x(i) = obj.scaledVar(i)+pert;
                des = x.*obj.designScale+obj.lb;
                modMesh = obj.orgMesh.Points + reshape(obj.gradMesh*(x(:)-obj.scaledVar(:)),size(obj.orgMesh.Points));
                if obj.settingUNGRADE.femCouplingFlag == 2 || strcmpi(obj.settingUNGRADE.gradientCalcMethod,"chain")
                    if any(obj.flowNoList(:,3) == 1)
                        obj2 = obj.makeApproximatedInstance(modMesh);
                    else
                        obj2 = obj.setVerts(modMesh);
                    end
                else
                    if any(obj.flowNoList(:,3) == 1)
                        obj2 = obj.setVerts(modMesh);
                        obj2 = obj2.makeEquation();
                    else
                        obj2 = obj.setVerts(modMesh);
                    end
                end


                %基準面積等の設計変数変化
                SREF2 = obj.SREF+obj.gradSREF*(x(:)-obj.scaledVar(:));
                BREF2 = obj.BREF+obj.gradBREF*(x(:)-obj.scaledVar(:));
                CREF2 = obj.CREF+obj.gradCREF*(x(:)-obj.scaledVar(:));
                XYZREF2 = obj.XYZREF+(obj.gradXYZREF*(x(:)-obj.scaledVar(:)))';
                argin_x2 = obj.argin_x+obj.gradArginx*(x(:)-obj.scaledVar(:));
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2,XYZREF2);

                if obj.settingUNGRADE.femCouplingFlag > 0
                    %基準メッシュから構造メッシュへの投影
                    modFemMesh = obj.femMeshDeformation(obj,modMesh,obj.orgMesh.ConnectivityList);
                    if obj.settingUNGRADE.femCouplingFlag == 2 || strcmpi(obj.settingUNGRADE.gradientCalcMethod,"chain")
                        femObj = obj.makeApproximatedFemInstance(modFemMesh);
                    else
                        femobj = obj.setFemVerts(modFemMesh,0);
                        femObj = femobj.makeFemEquation();
                    end
                end

                calcCount = 1;
                Rf = [];
                Sf = [];
                AERODATA = AERODATA0;
                Cp = Cp0;
                Cfe = Cfe0;
                DYNCOEF = DYNCOEF0;
                for iter = 1:numel(obj.flow)
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    if obj.settingUNGRADE.femCouplingFlag > 0
                        dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    end
                    for jter = 1:numel(alphabuff)
                        u0solve = u0((calcCount-1)*nbPanel+1:calcCount*nbPanel,1);
                        if obj.settingUNGRADE.femCouplingFlag > 1
                            modVerts = obj2.calcModifiedVerts(delta0{calcCount});
                            aeroDesObj = aeroObj{calcCount}.makeApproximatedInstance(modVerts);
                        else
                            aeroDesObj = obj2;
                        end
                        [AERODATAbuff,Cpbuff,Cfebuff,Rbuff,aeroDesObj] = aeroDesObj.solveFlowForAdjoint(u0solve,iter,alphabuff(jter),betabuff(jter));
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            [udwfdx,udwrdx] = aeroDesObj.calcDynCoefdu(iter,alphabuff(jter),betabuff(jter));
                            DYNCOEFbuff =  aeroDesObj.calcDynCoefforAdjoint(u0solve,udwfdx,udwrdx,iter,alphabuff(jter),betabuff(jter));
                            DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                        else
                            DYNCOEF = [];
                        end
                        AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                        Cp{iter}(:,jter) = Cpbuff{iter};
                        Cfe{iter} = Cfebuff{iter};
                        Rf = [Rf;Rbuff];
                        %%%%%%%%%%%%%%%
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            deltaIn{1} = delta0{calcCount};
                            [~,~,Sbuff] = femObj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                            Sf = [Sf;Sbuff];
                        end
                        calcCount = calcCount + 1;
                    end
                end
                if obj.settingUNGRADE.femCouplingFlag > 0
                    [geom,aero,fem] = obj2.makeObjArguments(obj2,des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2,delta0);
                else
                    [geom,aero,fem] = obj2.makeObjArguments(obj2,des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2,[]);
                end
                [If,conf] = objandConsFun(geom,aero,fem);
                conf = (conf-cmin)./cScale;
                clear obj2 aeroDesObj femObj


                x = obj.scaledVar;
                x(i) = obj.scaledVar(i)-pert;
                des = x.*obj.designScale+obj.lb;
                %[~,modMesh] = obj.calcApproximatedMeshGeom(des);
                modMesh = obj.orgMesh.Points + reshape(obj.gradMesh*(x(:)-obj.scaledVar(:)),size(obj.orgMesh.Points));
                if obj.settingUNGRADE.femCouplingFlag == 2 || strcmpi(obj.settingUNGRADE.gradientCalcMethod,"chain")
                    if any(obj.flowNoList(:,3) == 1)
                        obj2 = obj.makeApproximatedInstance(modMesh);
                    else
                        obj2 = obj.setVerts(modMesh);
                    end
                else
                    if any(obj.flowNoList(:,3) == 1)
                        obj2 = obj.setVerts(modMesh);
                        obj2 = obj2.makeEquation();
                    else
                        obj2 = obj.setVerts(modMesh);
                    end
                end

                %基準面積等の設計変数変化
                SREF2 = obj.SREF+obj.gradSREF*(x(:)-obj.scaledVar(:));
                BREF2 = obj.BREF+obj.gradBREF*(x(:)-obj.scaledVar(:));
                CREF2 = obj.CREF+obj.gradCREF*(x(:)-obj.scaledVar(:));
                XYZREF2 = obj.XYZREF+(obj.gradXYZREF*(x(:)-obj.scaledVar(:)))';
                argin_x2 = obj.argin_x+obj.gradArginx*(x(:)-obj.scaledVar(:));
                obj2 = obj2.setREFS(SREF2,BREF2,CREF2,XYZREF2);

                if obj.settingUNGRADE.femCouplingFlag > 0
                    %基準メッシュから構造メッシュへの投影
                    modFemMesh = obj.femMeshDeformation(obj,modMesh,obj.orgMesh.ConnectivityList);
                    femObj = obj.makeApproximatedFemInstance(modFemMesh);
                end

                calcCount = 1;
                Rr = [];
                Sr = [];
                AERODATA = AERODATA0;
                Cp = Cp0;
                Cfe = Cfe0;
                DYNCOEF = DYNCOEF0;
                for iter = 1:numel(obj.flow)
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    if obj.settingUNGRADE.femCouplingFlag > 0
                        dynbuff = obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    end

                    for jter = 1:numel(alphabuff)
                        u0solve = u0((calcCount-1)*nbPanel+1:calcCount*nbPanel,1);
                        if obj.settingUNGRADE.femCouplingFlag > 1
                            modVerts = obj2.calcModifiedVerts(delta0{calcCount});
                            aeroDesObj = aeroObj{calcCount}.makeApproximatedInstance(modVerts);
                        else
                            aeroDesObj = obj2;
                        end
                        [AERODATAbuff,Cpbuff,Cfebuff,Rbuff,aeroDesObj] = aeroDesObj.solveFlowForAdjoint(u0solve,iter,alphabuff(jter),betabuff(jter));
                        if obj.settingUNGRADE.dynCoefFlag == 1
                            [udwfdx,udwrdx] = aeroDesObj.calcDynCoefdu(iter,alphabuff,betabuff);
                            DYNCOEFbuff =  aeroDesObj.calcDynCoefforAdjoint(u0solve,udwfdx,udwrdx,iter,alphabuff(jter),betabuff(jter));
                            DYNCOEF{iter,jter} = DYNCOEFbuff{iter,1};
                        else
                            DYNCOEF = [];
                        end
                        AERODATA{iter}(jter,:) = AERODATAbuff{iter};
                        Cp{iter}(:,jter) = Cpbuff{iter};
                        Cfe{iter} = Cfebuff{iter};
                        Rr = [Rr;Rbuff];
                        %%%%%%%%%%%%%%%
                        if obj.settingUNGRADE.femCouplingFlag > 0
                            deltaIn{1} = delta0{calcCount};
                            [~,~,Sbuff] = femObj.solveFemForAdjoint(deltaIn,Cpbuff{iter}.*dynbuff(jter),obj.settingUNGRADE.selfLoadFlag);
                            Sr = [Sr;Sbuff];
                        end
                        calcCount = calcCount + 1;
                    end
                end
                if obj.settingUNGRADE.femCouplingFlag > 0
                    [geom,aero,fem] = obj2.makeObjArguments(obj2,des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2,delta0);
                else
                    [geom,aero,fem] = obj2.makeObjArguments(obj2,des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2,[]);
                end
                [Ir,conr] = objandConsFun(geom,aero,fem);
                conr = (conr-cmin)./cScale;
                clear obj2 aeroDesObj femObj

                dI_dx(i) = (If-Ir)./pert/2;%評価関数の設計変数に関する偏微分
                if not(isempty(conr))
                    dcon_dx(:,i) = (conf-conr)./pert/2;
                end
                dR_dx(:,i) = (Rf-Rr)./pert/2;
                if obj.settingUNGRADE.femCouplingFlag > 0
                    dS_dx(:,i) = (Sf-Sr)./pert/2;
                end
            end

            %全微分でdxからduを求める行列
            if obj.settingUNGRADE.femCouplingFlag > 0
                duddeltaMatVec = lsqminnorm(-([dR_du,dR_ddelta;dS_du./10000,dS_ddelta./10000]),([R0,dR_dx;S0./10000,dS_dx./10000]),obj.settingUNLSI.lsqminnormTol);
                duddeltaVec = duddeltaMatVec(:,1);
                duddeltaMat = duddeltaMatVec(:,2:end);
                objTotalGrad = dI_dx+[dI_du,dI_ddelta]*duddeltaMat;
                obj.LagrangianInfo.objTotalGrad = objTotalGrad;
                obj.LagrangianInfo.I0 = I0;
                if not(isempty(conf))
                    conTotalGrad = dcon_dx+[dcon_du,dcon_ddelta]*duddeltaMat;
                    errorMod = [dcon_du,dcon_ddelta]*duddeltaVec;
                end
            else
                duMatVec = -(dR_du)\([R0,dR_dx]);
                duMat = duMatVec(:,2:end);
                objTotalGrad = dI_dx+dI_du*duMat;
                obj.LagrangianInfo.objTotalGrad = objTotalGrad;
                obj.LagrangianInfo.I0 = I0;
                if not(isempty(conf))
                    conTotalGrad = dcon_dx+dcon_du*duMat;
                    errorMod = 0;
                end
            end
            %線形計画問題に変換
            lbf = -obj.scaledVar;
            ubf = 1-obj.scaledVar;
            obj.LagrangianInfo.lactiveBound = obj.scaledVar'<=0;
            obj.LagrangianInfo.uactiveBound = obj.scaledVar'>=1;
            lbflin = lbf;
            lbflin(obj.LagrangianInfo.lactiveBound==0) = -Inf;
            ubflin = ubf;
            ubflin(obj.LagrangianInfo.uactiveBound==0) =  Inf;
            if not(isempty(conf))
                obj.LagrangianInfo.clactiveBound = con0 <= 0-errorMod;
                obj.LagrangianInfo.cuactiveBound = con0 >= 1-errorMod;
                bothInacitve = and(obj.LagrangianInfo.clactiveBound==0,obj.LagrangianInfo.cuactiveBound==0);
                obj.LagrangianInfo.clactiveBound(bothInacitve) = 1;
                obj.LagrangianInfo.cuactiveBound(bothInacitve) = 1;

            else
                obj.LagrangianInfo.clactiveBound = [];
                obj.LagrangianInfo.cuactiveBound = [];
            end
            if not(isempty(conf))
                obj.LagrangianInfo.alin = [-conTotalGrad;conTotalGrad];
                obj.LagrangianInfo.blin = [con0(:)-0+errorMod;1-con0(:)-errorMod];
                alinlin = obj.LagrangianInfo.alin;
                blinlin = obj.LagrangianInfo.blin;
                alinlin([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]==0,:) = [];
                blinlin([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]==0,:) = [];
            else
                obj.LagrangianInfo.alin = [];
                obj.LagrangianInfo.blin = [];
            end

            options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','MaxFunctionEvaluations',10000);
            [~,~,exitFlag,~,obj.LagrangianInfo.lambda] = fmincon(@(x)obj.fminconObj(x,eye(numel(obj.scaledVar)),objTotalGrad),zeros(numel(obj.scaledVar),1),alinlin,blinlin,[],[],lbflin,ubflin,[],options);
            if exitFlag<1
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','MaxFunctionEvaluations',10000);
                [~,~,~,~,obj.LagrangianInfo.lambda] = fmincon(@(x)obj.fminconObj(x,eye(numel(obj.scaledVar)),objTotalGrad),zeros(numel(obj.scaledVar),1),alinlin,blinlin,[],[],lbflin,ubflin,[],options);
            end
            %lambdaRS = (sparse([designUpdate.dR_du,zeros(size(designUpdate.dR_du,1),size(designUpdate.dS_ddelta,2));designUpdate.dS_du,designUpdate.dS_ddelta]'))\(sparse(-[designUpdate.dI_du';designUpdate.dI_ddelta']-(lambdaC'*([designUpdate.dC_du,designUpdate.dC_ddelta]))'));
            if obj.settingUNGRADE.femCouplingFlag > 0
                if not(isempty(conf))
                    lambdaC = zeros(size(obj.LagrangianInfo.blin,1),1);
                    if not(isempty(obj.LagrangianInfo.lambda.ineqlin))
                        lambdaC([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]==1,1) = obj.LagrangianInfo.lambda.ineqlin;
                    end
                    lambdaRS = lsqminnorm(-([dR_du,dR_ddelta;dS_du./10000,dS_ddelta./10000]'),([dI_du,dI_ddelta]+(lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-[dcon_du,dcon_ddelta];[dcon_du,dcon_ddelta]])',obj.settingUNLSI.lsqminnormTol) ;
                    lambdaR = lambdaRS(1:size(dR_du,1),1);
                    lambdaS = lambdaRS(size(dR_du,1)+1:end,1);
                    obj.LagrangianInfo.Lorg = I0 + (lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[0-con0;con0-1;]+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                    obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+lambdaS'*dS_dx./10000+(lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-dcon_dx;dcon_dx]+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
                else
                    lambdaRS = lsqminnorm(-([dR_du,dR_ddelta;dS_du./10000,dS_ddelta./10000]'),[dI_du,dI_ddelta]',obj.settingUNLSI.lsqminnormTol);
                    lambdaR = lambdaRS(1:size(dR_du,1),1);
                    lambdaS = lambdaRS(size(dR_du,1)+1:end,1);
                    obj.LagrangianInfo.Lorg = I0+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                    obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+lambdaS'*dS_dx./10000+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
                end
            else
                if not(isempty(conf))
                    lambdaC = zeros(size(obj.LagrangianInfo.blin,1),1);
                    if not(isempty(obj.LagrangianInfo.lambda.ineqlin))
                        lambdaC([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]==1,1) = obj.LagrangianInfo.lambda.ineqlin;
                    end
                    lambdaR = lsqminnorm(-(dR_du'),(dI_du+(lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-dcon_du;dcon_du])',obj.settingUNLSI.lsqminnormTol);
                    obj.LagrangianInfo.Lorg = I0 + (lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[0-con0;con0-1;]+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                    obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+(lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-dcon_dx;dcon_dx]+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
                else
                    lambdaR = lsqminnrom(-(dR_du'),dI_du',obj.settingUNLSI.lsqminnormTol);
                    obj.LagrangianInfo.Lorg = I0+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                    obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
                end
            end
            %obj.LagrangianInfo.dLorg = obj.LagrangianInfo.Lorg;
            %obj.LagrangianInfo.dL_dx = obj.LagrangianInfo.dL_dx;

            obj.history.LagrangianVal(obj.iteration) = obj.LagrangianInfo.Lorg;
            obj.history.scaledVar(obj.iteration,:) = obj.scaledVar;
            obj.history.dL_dx(obj.iteration,:) = obj.LagrangianInfo.dL_dx;
            fprintf("Lagrangian Value : %f \nGradient of Lagrangian : \n",obj.LagrangianInfo.Lorg);
            disp(obj.LagrangianInfo.dL_dx);
            fprintf("iteration No. %d : Gradient Caluculation Completed\n",obj.iteration);
        end

        function [nextUnscaledVar,scaleddxNorm] = descentLagrangian(obj,varargin)
            %%%%%%%%%%%%%Lagrangianの降下方向に向かう設計変数を計算する%%%%%%%%%
            %
            %dog-legとsteepest-descentはアルゴリズム対応が途中
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin>1
                obj = obj.setUNGRADESettings(varargin{:});
            end
            fprintf("--updateMethod : %s\n",obj.settingUNGRADE.updateMethod)

            if strcmpi(obj.settingUNGRADE.updateMethod,"Quasi-Newton")
                fprintf("--TrustRegion : %f\n",obj.settingUNGRADE.TrustRegion)
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"Levenberg–Marquardt")
                fprintf("--TrustRegion : %f\n",obj.settingUNGRADE.TrustRegion)
                fprintf("--Beta(Levenberg-Marquardt) : %f\n",obj.settingUNGRADE.betaLM)
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"dogleg")
                fprintf("--TrustRegion : %f\n",obj.settingUNGRADE.TrustRegion)
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"Steepest-Descent")
                fprintf("--Alpha(Steepest-Descent) : %f\n",obj.settingUNGRADE.alphaSD)
                fprintf("--TrustRegion : %f\n",obj.settingUNGRADE.TrustRegion)
            end
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            lbf = -obj.scaledVar;
            ubf = 1-obj.scaledVar;
            %lbf(lbf<-obj.settingUNGRADE.dxMax) = -obj.settingUNGRADE.dxMax;
            %ubf(ubf>obj.settingUNGRADE.dxMax) = obj.settingUNGRADE.dxMax;
            ndim = numel(lbf);
            ncon = numel(obj.LagrangianInfo.blin);
            if strcmpi(obj.settingUNGRADE.updateMethod,"Quasi-Newton")
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,obj.Hessian,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion,ndim),options);
                dxscaled = dxscaled(1:ndim,1);
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"Levenberg–Marquardt")
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,obj.Hessian+obj.settingUNGRADE.betaLM*diag(diag(obj.Hessian)),obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion,ndim),options);
                dxscaled = dxscaled(1:ndim,1);
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"dogleg")
                %%%%%%%%%%%%%%%%%%%%%%%%%%%dog-legとsteepest-descentはアルゴリズム対応が途中
                %sd方向を求める行列
                Aeq = [eye(ndim),obj.LagrangianInfo.dL_dx(:)];
                beq = zeros(ndim,1);
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                xsd_t = fmincon(@(x)obj.fminconObj(x,obj.Hessian,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar)+1,1),[],[],Aeq,beq,[lbf,-100],[ubf,100],[],options);
                xsd = xsd_t(1:ndim);
                xqn = fmincon(@(x)obj.fminconObj(x,obj.Hessian,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,[],options);
                if norm(xqn) <= obj.settingUNGRADE.TrustRegion
                    dxscaled = xqn;
                elseif norm(xqn) > obj.settingUNGRADE.TrustRegion && norm(xsd) > obj.settingUNGRADE.TrustRegion
                    dxscaled = (obj.settingUNGRADE.TrustRegion / norm(xsd) ).* xsd;
                else
                    sdl = fsolve(@(s)obj.doglegsearch(s,xsd,xqn,obj.settingUNGRADE.TrustRegion),1);
                    dxscaled = xsd + sdl .* (xqn-xsd);
                end
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"Steepest-Descent")
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,2.*eye(ndim),2.*obj.settingUNGRADE.alphaSD.*obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion,ndim),options);
                %dxscaled = -obj.settingUNGRADE.alphaSD.*obj.LagrangianInfo.dL_dx';
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"bb")
                options = optimoptions(@fmincon,'Algorithm','sqp','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,obj.BBMatrix,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion,ndim),options);
            end

            %線形化した厳密な評価関数を用いて直線探索をする。
            %             rho = 1000;
            %             dx0ls = dxscaled;
            %             fvalOrg = rho*sum(max([-obj.LagrangianInfo.blin,zeros(size(obj.LagrangianInfo.blin))],[],2));
            %             for i = 1:300
            %                 conPred = obj.LagrangianInfo.alin * dxscaled(:)-obj.LagrangianInfo.blin;
            %                 fvaldx = obj.LagrangianInfo.dL_dx*dxscaled(:) + rho*sum(max([conPred,zeros(size(obj.LagrangianInfo.blin))],[],2));
            %                 if fvaldx>fvalOrg
            %                     dxscaled = 0.99.*dxscaled;
            %                     if all(obj.settingUNGRADE.meshGradientPerturbation > abs(dxscaled(:)).*obj.designScale(:))
            %                         dxscaled = dxscaled./0.99;
            %                         break;
            %                     end
            %                 else
            %                     break;
            %                 end
            %             end
            %             fprintf("\nSimplified Line Search : alpha %f \n",norm(dxscaled)/norm(dx0ls));
            dx = dxscaled(:)'.*obj.designScale;
            fprintf("dx Norm : %f  unscaled : %f \n",norm(dxscaled),norm(dx));
            nextUnscaledVar = desOrg + dx(:)';
            scaleddxNorm = norm(dxscaled);
        end

        function [Lval,I,con,obj] = calcLagrangian(obj,objandConsFun,cmin,cmax)
            [I,con,obj] = evaluateObjFun(obj,objandConsFun);
            cScale = cmax-cmin;
            con = (con-cmin)./cScale;
            if not(isempty(con))
                lambdaC = zeros(numel(con)*2,1);
                if not(isempty(obj.LagrangianInfo.lambda.ineqlin))
                    lambdaC([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]==1,1) = obj.LagrangianInfo.lambda.ineqlin;
                end
                Lval = I + (lambdaC'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[0-con;con-1;]+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
            else
                Lval = I+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
            end
        end

        function obj = updateHessian(obj,varargin)
            %%%%%%%%%%%%%%%%ヘッシアン更新%%%%%%%%%%%
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin>1
                obj = obj.setUNGRADESettings(varargin{:});
            end
            fprintf("-- HessianUpdateMethod : %s\n",obj.settingUNGRADE.HessianUpdateMethod);
            if strcmpi(obj.settingUNGRADE.HessianUpdateMethod,"SR1")
                updateFunction = @obj.SR1;
            elseif strcmpi(obj.settingUNGRADE.HessianUpdateMethod,"BFGS")
                updateFunction = @obj.BFGS;
            elseif strcmpi(obj.settingUNGRADE.HessianUpdateMethod,"DFP")
                updateFunction = @obj.DFP;
            elseif strcmpi(obj.settingUNGRADE.HessianUpdateMethod,"Broyden")
                updateFunction = @obj.Broyden;
            end
            %Hessianの更新
            if size(obj.history.scaledVar,1) > 1
                n_iter = size(obj.history.scaledVar,1)-1;
                %                 if n_iter > obj.settingUNGRADE.nMemory
                %                     n_iter = obj.settingUNGRADE.nMemory;
                %                 end
                for i = 1:n_iter
                    s = obj.history.scaledVar(end,:)-obj.history.scaledVar(end-i,:);
                    y = obj.history.dL_dx(end,:)-obj.history.dL_dx(end-i,:);
                    if norm(s)>sqrt(eps) && norm(y)>sqrt(eps)
                        break;
                    end
                end
                obj.Hessian = updateFunction(s,y,obj.Hessian);
                obj.BBMatrix = obj.Barzilai_Borwein(s,y);
            end
            if issymmetric(obj.Hessian)
                fprintf("Hessian is symmetric ");
                d = eig(obj.Hessian);
                if all(d > 0)
                    fprintf("and positive-define\n");
                elseif all(d >= 0)
                    fprintf("and semi-positive-define\n");
                else
                    fprintf("and not positive-define\n");
                end
            else
                fprintf("Hessian is Unsymmetric\n");
            end
            fprintf("Hessian diag value -->\n");
            disp(diag(obj.Hessian)');
        end

        function [nextUnscaledVar,obj] = calcNextVariables(obj,objandConsFun,cmin,cmax,varargin)
            %%%%%%%%%%%%後方互換性用、とりあえず次の変数を計算する%%%%%%%%
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin>4
                obj = obj.setUNGRADESettings(varargin{:});
            end
            obj = obj.calcLagrangianGradient(objandConsFun,cmin,cmax);%次の設計変数を計算する
            obj = obj.updateHessian();%次の設計変数を計算する
            nextUnscaledVar = obj.descentLagrangian();%次の設計変数を計算する
        end

        function plotOptimizationState(obj,fig,filename)
            %%%%%%%%%%%%最適化の履歴をプロットする%%%%%%%%
            %fig : 描画するFigure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(obj.history.conVal)
                wNum = 4;
            else
                wNum = 5;
            end
            figure(fig);clf;
            plotiter = 1:obj.iteration;

            subplot(wNum,1,1)
            plot(plotiter,obj.history.scaledVar,"-o","LineWidth",1);
            ylabel("Var(scaled)");
            set(gca,"FontSize",14,'xticklabel',[]);grid on;

            subplot(wNum,1,2);grid on;
            plot(plotiter,obj.history.dL_dx,"-o","LineWidth",1);
            %plot(plotiter,obj.history.dL_dxLP,"--o","LineWidth",1);hold off;
            ylabel("dL/dx");
            set(gca,"FontSize",14,'xticklabel',[]);grid on;
            subplot(wNum,1,3);grid on;
            plot(plotiter,obj.history.LagrangianVal,"-o","LineWidth",1);
            %plot(plotiter,obj.history.LagrangianValLP,"--o","LineWidth",1);hold off;
            ylabel("Lag val");
            set(gca,"FontSize",14,'xticklabel',[]);grid on;
            subplot(wNum,1,4);grid on;
            plot(plotiter,obj.history.objVal,"-o","LineWidth",1);
            ylabel("Obj val");

            if not(isempty(obj.history.conVal))
                set(gca,"FontSize",14,'xticklabel',[]);grid on;
                subplot(wNum,1,5);grid on;hold on;
                for i = 1:size(obj.history.conVal,1)
                    plot(plotiter,obj.history.conVal(i,:),"-o","LineWidth",1);
                end
                hold off;
                ylabel("Cons val");
                set(gca,"FontSize",14);grid on;
            else
                set(gca,"FontSize",14);grid on;
            end
            xlabel("Iteration")
            set(gca,"FontSize",14);

            drawnow();
            if nargin > 2
                print(filename,'-dpng');
            end
        end
    end

    methods(Static)

        function [geom,aero,fem] = makeObjArguments(obj,x,AERODATA,DYNCOEF,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x,delta,invR)
            %geom:機体構造関連
            %   geom.x:設計変数
            %   geom.aeroMesh:空力解析の解析メッシュ
            %   geom.surfID:ID
            %   geom.SREF,BREF,CREF,XYZREF,argin_x
            %aero:空力解析関連
            %   aero.floadata:関数の流れの状態
            %   aero.DATA;空力解析結果
            %       1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt
            %   aero.DYNCOEF:動安定微係数（フラグを1にした場合のみ）
            %   aero.Cp([x,y,z]):CpをRBFN補間する関数
            %   aero.Cfe([x,y,z]):CfeをRBFN補間する関数
            %fem:構造解析関連
            %   fem.mesh:構造解析メッシュ
            %   fem.delta([x,y,z]):変位をRBFN補間する関数。空力メッシュに投影したものなので注意
            geom.x = x;
            geom.aeroMesh = obj.tri;
            geom.surfID = obj.surfID;
            geom.SREF = SREF;geom.CREF = CREF;geom.BREF = BREF;geom.XYZREF = XYZREF;geom.argin_x = argin_x;
            aero.DATA = AERODATA;
            aero.DYNCOEF = DYNCOEF;

            iter = 1;
            for i = 1:numel(obj.flow)
                aero.flowdata{i} = AERODATA{i}(:,1:4);
                if nargin < 13
                    ppCp = obj.RbfppMake(obj,obj.center,Cp{i},1,0.01);
                    ppCfe = obj.RbfppMake(obj,obj.center,Cfe{i},1,0.01);
                else
                    ppCp = obj.invRppMake(invR.center,Cp{i});
                    ppCfe = obj.invRppMake(invR.center,Cfe{i});
                end
                aero.Cp{i} = @(val_interp)obj.execRbfInterp(obj,ppCp,val_interp);
                aero.Cfe{i} = @(val_interp)obj.execRbfInterp(obj,ppCfe,val_interp);

                if not(isempty(delta))
                    for jter = 1:numel(obj.settingUNGRADE.dynPress(obj.flowNoList(:,2)==obj.flow{i}.Mach))
                        aerodelta = zeros(size(obj.tri.Points));
                        aerodelta(obj.femutils.usedAeroVerts,1) = obj.fem2aeroMat*delta{iter}(obj.femutils.usedVerts,1);
                        aerodelta(obj.femutils.usedAeroVerts,2) = obj.fem2aeroMat*delta{iter}(obj.femutils.usedVerts,2);
                        aerodelta(obj.femutils.usedAeroVerts,3) = obj.fem2aeroMat*delta{iter}(obj.femutils.usedVerts,3);
                        if nargin < 13
                            ppDelta = obj.RbfppMake(obj,obj.tri.Points,aerodelta,1,0.01);
                        else
                            ppDelta = obj.invRppMake(invR.verts,aerodelta);
                        end
                        fem.delta{iter} = @(val_interp)obj.execRbfInterp(obj,ppDelta,val_interp);
                        iter = iter+1;
                    end
                else
                    iter = iter+1;
                end
            end
            if not(isempty(delta))
                fem.femmesh = obj.femtri;
            else
                fem = [];
            end

            %fem.delta = @(val_interp)obj.execRbfInterp(obj,pp,val_interp);
        end

        function res = doglegsearch(s,xsd,xqn,TR)
            vec = xsd+s*(xqn-xsd);
            res = sqrt(sum(vec.^2))-TR;
        end

        function res = fminconObj(dx,H,g)
            %%%SQPの更新ベクトルの目的関数%%%%%
            %dx：更新ベクトル
            %H：ヘシアン
            %g：勾配
            %%%%%%%%%%%%%%%%%%%%%%
            ndim = numel(g);
            res = 0.5 * dx(1:ndim)'*H*dx(1:ndim) + g*dx(1:ndim);
        end

        function [c,ceq] = fminconNlc(dx,TR,nvar)
            ceq = [];
            c = sum(dx(1:nvar).^2)/TR^2-1;
        end

        function [modVerts,con] = meshDeformation(obj,modGeomVerts,modGeomCon)
            %%%%%%%%%%%%非構造メッシュのメッシュ変形%%%%%%%
            %補間ベースのメッシュ変形
            %modSurfVers : 変形後の基準サーフェスの接点情報
            %modSurfCon : 入力すると一応conが一致しているかチェックする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 3
                %connectivityが一致しているか確認
                if any(obj.orgGeom.ConnectivityList(:) ~= modGeomCon(:))
                    error("Surf connectivity is not match")
                end
            end

            modGeomVerts = modGeomVerts-obj.orgGeom.Points;
            dVerts(:,1) = obj.geom2MeshMat*modGeomVerts(:,1);
            dVerts(:,2) = obj.geom2MeshMat*modGeomVerts(:,2);
            dVerts(:,3) = obj.geom2MeshMat*modGeomVerts(:,3);
            modVerts = obj.tri.Points+dVerts;
            con = obj.tri.ConnectivityList;
        end
        function [modVerts,con] = femMeshDeformation(obj,modMeshVerts,modMeshCon)
            %%%%%%%%%%%%非構造メッシュのメッシュ変形%%%%%%%
            %補間ベースのメッシュ変形
            %modSurfVers : 変形後の基準サーフェスの接点情報
            %modSurfCon : 入力すると一応conが一致しているかチェックする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 3
                %connectivityが一致しているか確認
                if any(obj.orgMesh.ConnectivityList(:) ~= modMeshCon(:))
                    error("Surf connectivity is not match")
                end
            end

            modMeshVerts = modMeshVerts-obj.orgMesh.Points;
            dVerts(:,1) = obj.aero2femMat*modMeshVerts(obj.femutils.usedAeroVerts,1);
            dVerts(:,2) = obj.aero2femMat*modMeshVerts(obj.femutils.usedAeroVerts,2);
            dVerts(:,3) = obj.aero2femMat*modMeshVerts(obj.femutils.usedAeroVerts,3);
            modVerts = obj.femtri.Points;
            modVerts(obj.femutils.usedVerts,:) = obj.femtri.Points(obj.femutils.usedVerts,:)+dVerts;
            con = obj.femtri.ConnectivityList;
        end
        function obj = makeMeshGradient(obj,pert)
            %%%%%%%%%%%設計変数勾配による標本表面の近似関数の作成%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = numel(obj.scaledVar);
            obj.gradMesh = zeros(size(obj.orgMesh.Points(:),1),ndim);
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            [surforg,~,~,~,~,~,~,desOrg] = obj.geomGenFun(desOrg);
            for i = 1:ndim
                sampleDes = obj.scaledVar.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.scaledVar(i) + pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFf,BREFf,CREFf,XYZREFf,argin_xf,desBuff] = obj.geomGenFun(sampleDes);
                pertf = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                modMeshf(:,1) = obj.geom2MeshMat*dmodSurf(:,1);
                modMeshf(:,2) = obj.geom2MeshMat*dmodSurf(:,2);
                modMeshf(:,3) = obj.geom2MeshMat*dmodSurf(:,3);
                sampleDes = obj.scaledVar.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.scaledVar(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFr,BREFr,CREFr,XYZREFr,argin_xr,desBuff] = obj.geomGenFun(sampleDes);
                pertr = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                modMeshr(:,1) = obj.geom2MeshMat*dmodSurf(:,1);
                modMeshr(:,2) = obj.geom2MeshMat*dmodSurf(:,2);
                modMeshr(:,3) = obj.geom2MeshMat*dmodSurf(:,3);

                obj.gradSREF(i) = (SREFf-SREFr)./(pertf-pertr);
                obj.gradBREF(i) = (BREFf-BREFr)./(pertf-pertr);
                obj.gradCREF(i) = (CREFf-CREFr)./(pertf-pertr);
                obj.gradXYZREF(:,i) = (XYZREFf(:)-XYZREFr(:))./(pertf-pertr);
                obj.gradArginx(:,i) = (argin_xf(:)-argin_xr(:))./(pertf-pertr);
                obj.gradMesh(:,i) = (modMeshf(:)-modMeshr(:))./(pertf-pertr);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        %準ニュートン法のアップデートの手実装
        %s 設計変数の変化量
        %y ヤコビアンの変化
        %%%%%%%%%%%%%%%%%%%%%%%%
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
            Bkp1 = Bk+((y-Bk*s)/(s'*s))*s';
            %Bkp1 = (Bkp1+Bkp1')./2;
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
        function Bkp1 = Barzilai_Borwein(s,y)
            s = s(:);
            y=y(:);
            %delta = 100.0;
            alphak = (s'*y)/(s'*s);
            %alpha = (s'*y)/(y'*y);
            %if alpha <= 0
            %    alpha = norm(s)/norm(y);
            %end
            %alphak = min(alpha, delta / norm(y));
            Bk = eye(numel(s));
            Bkp1 = alphak * Bk;
        end

    end

end