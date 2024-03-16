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
        meshScaling%メッシュ行列のスケーリング変数
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
            obj.settingUNGRADE.gradientCalcMethod = "direct"; %設計変数偏微分の取得方法："direct", "chain", "nonlin"
            obj.settingUNGRADE.HessianUpdateMethod = "BFGS"; %ヘッシアン更新は "BFGS","DFP","Broyden","SR1"から選択
            obj.settingUNGRADE.updateMethod = "Levenberg–Marquardt";%解を更新する方法 Steepest-Descent,Levenberg–Marquardt,dogleg
            obj.settingUNGRADE.betaLM = 0.5; %レーベンバーグマッカートの重み係数
            obj.settingUNGRADE.alphaSD = 0.5; %最急降下法の重み 
            obj.settingUNGRADE.TrustRegion = 0.1; %設計更新を行う際のスケーリングされた設計変数更新量の最大値
            obj.settingUNGRADE.dynCoefFlag = 0;
            obj.settingUNGRADE.r0RBF = 1;
            %%%%%%%%%%%%%
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol,"delete");
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(desOrg);
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
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

        function obj = setFlowCondition(obj,alpha,beta,Mach,Re)
            obj.settingUNGRADE.Mach = Mach;
            obj.settingUNGRADE.alpha = alpha;
            obj.settingUNGRADE.beta = beta;
            obj.settingUNGRADE.Re = Re;
            if numel(obj.settingUNGRADE.Mach) ~= numel(obj.settingUNGRADE.alpha) || numel(obj.settingUNGRADE.Mach) ~= numel(obj.settingUNGRADE.beta) || numel(obj.settingUNGRADE.beta) ~= numel(obj.settingUNGRADE.alpha)
                error("No. of case is not match");
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
            end
            if any(obj.flowNoList(:,3) == 1)
                obj = obj.makeCluster();
            end
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
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
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
            u0 = obj.solvePertPotential(flowNo,alpha,beta);
            [~,~,~,~,obj] = obj.solveFlowForAdjoint(u0,flowNo,alpha,beta);
            if obj.settingUNGRADE.dynCoefFlag == 1
                [udwf,udwr] = obj.calcDynCoefdu(flowNo,alpha,beta);
                [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,flowNo,alpha,beta);
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
            obj.approximated = 0;
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            for iter = 1:numel(obj.flow)
                alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                u0 = obj.solvePertPotential(iter,alphabuff,betabuff);
                [~,~,~,~,obj] = obj.solveFlowForAdjoint(u0,iter,alphabuff,betabuff);
                if obj.settingUNGRADE.dynCoefFlag == 1
                    [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff,betabuff);
                    [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,iter,alphabuff,betabuff);
                end
            end
            [I,con] = objandConsFun(desOrg,obj.AERODATA,obj.DYNCOEF,obj.Cp,obj.Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x);
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
            obj.approximated = 0;
            desOrg = obj.scaledVar.*obj.designScale+obj.lb;
            u0 = zeros(nbPanel*size(obj.flowNoList,1),1);
            R0 = zeros(nbPanel*size(obj.flowNoList,1),1);
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
            end
            [I0,con0] = objandConsFun(desOrg,AERODATA0,DYNCOEF0,Cp0,Cfe0,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x);
            fprintf("Orginal Objective Value and Constraints:\n")
            disp([I0,con0(:)']);
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
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    [AERODATA,Cp,Cfe,~,obj] = obj.solveFlowForAdjoint(usolve,iter,alphabuff,betabuff);
                    if obj.settingUNGRADE.dynCoefFlag == 1
                        DYNCOEF =  obj.calcDynCoefforAdjoint(usolve,udwf{iter},udwr{iter},iter,alphabuff,betabuff);
                    else
                        DYNCOEF = [];
                    end
                end
                [I,con] = objandConsFun(desOrg,AERODATA,DYNCOEF,Cp,Cfe,obj.SREF,obj.BREF,obj.CREF,obj.XYZREF,obj.argin_x);
                dI_du(i) = (I-I0)/pert;%評価関数のポテンシャルに関する偏微分
                if not(isempty(con))
                    dcon_du(:,i) = (con-con0)/pert;
                end
            end
            for i = 1:size(obj.flowNoList,1)
                if i == 1
                    if obj.flowNoList(i,3) == 1
                        dR_du = -obj.LHS; %亜音速
                    else
                        dR_du = eye(nbPanel);
                    end
                else
                    if obj.flowNoList(i,3) == 1
                        dR_du = blkdiag(dR_du,-obj.LHS);
                    else
                        dR_du = blkdiag(dR_du,eye(nbPanel));
                    end
                end
            end
            %x微分の計算
            %メッシュの節点勾配を作成
            obj = obj.makeMeshGradient(obj,obj.settingUNGRADE.meshGradientPerturbation./obj.designScale);
            if strcmpi(obj.settingUNGRADE.gradientCalcMethod,'chain')
                if any(obj.flowNoList(:,3) == 1)
                    obj = obj.calcApproximatedEquation(obj);
                end
            end

            pert = sqrt(eps);
            for i= 1:numel(obj.scaledVar)
                x = obj.scaledVar;
                x(i) = obj.scaledVar(i)+pert;
                des = x.*obj.designScale+obj.lb;
                %[~,modMesh] = obj.calcApproximatedMeshGeom(des);
                modMesh = obj.orgMesh.Points + reshape(obj.gradMesh*(x(:)-obj.scaledVar(:)),size(obj.orgMesh.Points));
                if strcmpi(obj.settingUNGRADE.gradientCalcMethod,'direct')
                    %変数が少ないときは直接作成
                    obj2 = obj.setVerts(modMesh);
                elseif strcmpi(obj.settingUNGRADE.gradientCalcMethod,'chain')
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
                if strcmpi(obj.settingUNGRADE.gradientCalcMethod,'direct')
                    %変数が少ないときは直接作成
                    if any(obj.flowNoList(:,3) == 1)
                        obj2 = obj2.makeEquation();
                    end
                    obj2.approximated = 0;
                end
                R = zeros(nbPanel*size(obj.flowNoList,1),1);
                for iter = 1:numel(obj.flow)
                    lktable = find(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    u0solve = zeros(nbPanel*numel(lktable),1);
                    for k = 1:numel(lktable)
                        u0solve(nbPanel*(k-1)+1:nbPanel*k,1) = u0((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1);
                    end
                    alphabuff = obj.settingUNGRADE.alpha(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    betabuff = obj.settingUNGRADE.beta(obj.flowNoList(:,2)==obj.flow{iter}.Mach);
                    [AERODATA,Cp,Cfe,Rsolve,obj2] = obj2.solveFlowForAdjoint(u0solve,iter,alphabuff,betabuff);
                    if obj.settingUNGRADE.dynCoefFlag == 1
                        [udwfdx,udwrdx] = obj2.calcDynCoefdu(iter,alphabuff,betabuff);
                        DYNCOEF =  obj2.calcDynCoefforAdjoint(u0solve,udwfdx,udwrdx,iter,alphabuff,betabuff);
                    else
                        DYNCOEF = [];
                    end
                    for k = 1:numel(lktable)
                        R((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1) = Rsolve(nbPanel*(k-1)+1:nbPanel*k,1);
                    end
                end
                
                [I,con] = objandConsFun(des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                dI_dx(i) = (I-I0)/pert;%評価関数の設計変数に関する偏微分
                if not(isempty(con))
                    dcon_dx(:,i) = (con-con0)/pert;
                end
                dR_dx(:,i) = -(R-R0)./pert;
                clear obj2;
            end
            %全微分でdxからduを求める行列
            %invdR_du = pinv(dR_du);
            duMatVec = -(dR_du)\([R0,dR_dx]);
            %duMatVec = -(dR_du)\([R0,dR_dx]);
            duMat = duMatVec(:,2:end);
            %dxMat = -dR_dx'*dR_du;
            %duVec = duMatVec(:,1);
            objTotalGrad = dI_dx+dI_du*duMat;
            obj.LagrangianInfo.objTotalGrad = objTotalGrad;
            obj.LagrangianInfo.I0 = I0;
            if not(isempty(con))
                conTotalGrad = dcon_dx+dcon_du*duMat;
            end
            %線形計画問題に変換
            lbf = -obj.scaledVar;
            ubf = 1-obj.scaledVar;
            obj.LagrangianInfo.lactiveBound = obj.scaledVar'<=0;
            obj.LagrangianInfo.uactiveBound = obj.scaledVar'>=1;
            if not(isempty(con))
                obj.LagrangianInfo.clactiveBound = con0 <= cmin;
                obj.LagrangianInfo.cuactiveBound = con0 >= cmax;
            else
                obj.LagrangianInfo.clactiveBound = [];
                obj.LagrangianInfo.cuactiveBound = [];
            end
            if not(isempty(con))
                obj.LagrangianInfo.alin = [-conTotalGrad;conTotalGrad];
                obj.LagrangianInfo.blin = [con0(:)-cmin(:);cmax(:)-con0(:)];
            else
                obj.LagrangianInfo.alin = [];
                obj.LagrangianInfo.blin = [];
            end

            options = optimoptions(@linprog,'Algorithm','interior-point','Display','final','MaxIter',10000);
            [~,~,exitFlag,~,obj.LagrangianInfo.lambda] = linprog(objTotalGrad,obj.LagrangianInfo.alin,obj.LagrangianInfo.blin,[],[],lbf,ubf,options);
            if exitFlag<1
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed',"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                [~,~,~,~,obj.LagrangianInfo.lambda] = fmincon(@(x)obj.fminconObj(x,obj.Hessian,objTotalGrad),zeros(numel(obj.scaledVar),1),obj.LagrangianInfo.alin,obj.LagrangianInfo.blin,[],[],lbf,ubf,[],options);
            end
      
            if not(isempty(con))
                lambdaR = -(dR_du')\(dI_du+(obj.LagrangianInfo.lambda.ineqlin'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-dcon_du;dcon_du])';
                obj.LagrangianInfo.Lorg = I0 + (obj.LagrangianInfo.lambda.ineqlin'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[cmin-con0;con0-cmax;]+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+(obj.LagrangianInfo.lambda.ineqlin'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[-dcon_dx;dcon_dx]+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
            else
                lambdaR = -(dR_du')\dI_du';
                obj.LagrangianInfo.Lorg = I0+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
                obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+(obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'-(obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)';
            end
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

            if strcmpi(obj.settingUNGRADE.updateMethod,"Levenberg–Marquardt")
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
            lbf(lbf<-obj.settingUNGRADE.TrustRegion) = -obj.settingUNGRADE.TrustRegion;
            ubf(ubf>obj.settingUNGRADE.TrustRegion) = obj.settingUNGRADE.TrustRegion;
            ndim = numel(lbf);
            ncon = numel(obj.LagrangianInfo.blin);
            if strcmpi(obj.settingUNGRADE.updateMethod,"Levenberg–Marquardt")
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);     
                if any([obj.LagrangianInfo.cuactiveBound;obj.LagrangianInfo.clactiveBound]==1)
                    H = blkdiag(obj.Hessian+obj.settingUNGRADE.betaLM*diag(diag(obj.Hessian)),eye(ncon).*10);
                    dxscaled = fmincon(@(x)obj.fminconObj(x,H,[obj.LagrangianInfo.dL_dx,ones(1,ncon).*10]),zeros(numel(obj.scaledVar)+numel(obj.LagrangianInfo.blin),1),[obj.LagrangianInfo.alin,-eye(ncon)],obj.LagrangianInfo.blin,[],[],[lbf,ones(1,ncon).*0],[ubf,ones(1,ncon).*Inf],[],options);
                else
                    dxscaled = fmincon(@(x)obj.fminconObj(x,obj.Hessian+obj.settingUNGRADE.betaLM*diag(diag(obj.Hessian)),obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],lbf,ubf,[],options);
                end
                dxscaled = dxscaled(1:ndim,1);
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"dogleg")

                %%%%%%%%%%%%%%%%%%%%%%%%%%%dog-legとsteepest-descentはアルゴリズム対応が途中

                %sd方向を求める行列
                Aeq = [eye(ndim),obj.LagrangianInfo.dL_dx(:)];
                beq = zeros(ndim,1);
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
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
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,2.*eye(ndim),2.*obj.settingUNGRADE.alphaSD.*obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),[],[],[],[],[],[],@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion),options);
                %dxscaled = -obj.settingUNGRADE.alphaSD.*obj.LagrangianInfo.dL_dx';
            end
            dx = dxscaled(:)'.*obj.designScale;
            nextUnscaledVar = desOrg + dx(:)';
            scaleddxNorm = norm(dxscaled);
        end

        function [Lval,I,con,obj] = calcLagrangian(obj,objandConsFun,cmin,cmax)
            [I,con,obj] = evaluateObjFun(obj,objandConsFun);
            if not(isempty(con))
                Lval = I + (obj.LagrangianInfo.lambda.ineqlin'*diag([obj.LagrangianInfo.clactiveBound;obj.LagrangianInfo.cuactiveBound]))*[cmin-con;con-cmax;]+((obj.LagrangianInfo.lambda.upper.*obj.LagrangianInfo.uactiveBound)'*(obj.scaledVar'-1))'+((obj.LagrangianInfo.lambda.lower.*obj.LagrangianInfo.lactiveBound)'*(-obj.scaledVar'))';
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
                if n_iter > obj.settingUNGRADE.nMemory
                    n_iter = obj.settingUNGRADE.nMemory;
                end
                obj.Hessian = obj.settingUNGRADE.H0;
                for i = n_iter:-1:1
                    s = obj.history.scaledVar(end-(i-1),:)-obj.history.scaledVar(end-i,:);
                    y = obj.history.dL_dx(end-(i-1),:)-obj.history.dL_dx(end-i,:);
                    if norm(s)>sqrt(eps) && norm(y)>sqrt(eps)
                        obj.Hessian = updateFunction(s,y,obj.Hessian);
                    end
                end
            end
            if size(obj.history.scaledVar,1) > 1
                s = obj.history.scaledVar(end,:)-obj.history.scaledVar(end-1,:);
                y = obj.history.dL_dx(end,:)-obj.history.dL_dx(end-1,:);
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

        function plotOptimizationState(obj,fig)
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
            ylabel("dL/dx");
            set(gca,"FontSize",14,'xticklabel',[]);grid on;
            subplot(wNum,1,3);grid on;
            plot(plotiter,obj.history.LagrangianVal,"-o","LineWidth",1);
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
        end
    end

    methods(Static)
        
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