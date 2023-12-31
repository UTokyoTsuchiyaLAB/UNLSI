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
        designScale %ub-lbのスケーリングパラメータ
        gradMesh %メッシュ節点の変化を設計変数で一次近似したもの
        approxMat %パネル法行列の近似行列を作成するためのセル
        approximated %パネル法行列が近似されたものかどうか    
        gradSREF %基準面積の設計変数に対する一次近似
        gradBREF %横基準長の設計変数に対する一次近似
        gradCREF %縦基準長の設計変数に対する一次近似
        gradXYZREF %回転中心の設計変数に対する一次近似
        gradArginx %任意変数（重量など）の設計変数に対する一次近似
        argin_x %任意変数（重量やバッテリー搭載量などのユーザーで自由に設定する変数）
        iteration %設計更新を行った階数
        Hessian %準ニュートン法で近似されたヘッシアン
        history %設計更新の履歴
        flowNoList %マッハ数とflowNoの関連付け行列
        LagrangianInfo
    end

    methods(Access = public)

        function obj = UNGRADE(meshGenFun,geomGenFun,unscaledVariables,lb,ub,halfmesh)
            %%%%%%%%%%%%メッシュの登録と基準サーフェス生成関数の登録%%%%%%%%%%%%%%%
            %orgVerts,orgCon　実際に解析を行うメッシュ（openVSPのCFDツール等で生成した（基本的に）オーバーラップのない非構造メッシュ）
            %meshGenFun　設計変数（スケールされていない）から基準サーフェス（解析メッシュと同じ表面を持つ、トリムされていないメッシュ）
            %designVariables : 設計変数
            %lb：設計変数の下限
            %ub：設計変数の上限
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            [orgMeshVerts, orgMeshCon,surfID,wakeLineID, desOrg] = meshGenFun(unscaledVariables(:)');
            obj = obj@UNLSI(orgMeshVerts,orgMeshCon,surfID,wakeLineID,halfmesh);
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.designScale = (ub-lb);
            if any(obj.designScale<0)
                error("check the ub and lb setting");
            end

            obj.geomGenFun = geomGenFun;
            obj.meshGenFun = meshGenFun;


            %%%%%%%%%%%オプションのデフォルト設定
            % gradientCalcMethod
            % 'direct'-->直接パネル法行列を計算し差分を作成
            % 'chain'-->各節点位置の変化に対する近似パネル法行列と、設計変数変化に対する節点位置の変化を用いてチェインルールで勾配計算する。
            % 'nonlin'-->approxmatによる近似パネル法行列を使った直接最適化
            % 
            obj.settingUNGRADE.nMemory = 6; %記憶制限準ニュートン法のさかのぼり回数（厳密な記憶制限法は実装していない）
            obj.settingUNGRADE.H0 = eye(numel(obj.lb)); %初期ヘッシアン
            obj.settingUNGRADE.wakeLength = 20;
            obj.settingUNGRADE.meshGradientPerturbation = 0.01;%メッシュ勾配をとるときのスケールされてない設計変数の摂動量
            obj.settingUNGRADE.gradientCalcMethod = "direct"; %設計変数偏微分の取得方法："direct", "chain", "nonlin"
            obj.settingUNGRADE.HessianUpdateMethod = "BFGS"; %ヘッシアン更新は "BFGS","DFP","Broyden","SR1"から選択
            obj.settingUNGRADE.updateMethod = "Levenberg–Marquardt";%解を更新する方法 Steepest-Descent,Levenberg–Marquardt,dogleg
            obj.settingUNGRADE.betaLM = 0.5; %レーベンバーグマッカートの重み係数
            obj.settingUNGRADE.alphaSD = 0.5; %最急降下法の重み 
            obj.settingUNGRADE.TrustRegion = 0.1; %設計更新を行う際のスケーリングされた設計変数更新量の最大値
            %%%%%%%%%%%%%
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol,obj.settingUNLSI.checkMeshMethod);
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(desOrg);
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
            obj.unscaledVar = desOrg(:)';
            obj.scaledVar = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);
            
            obj.iteration = 0;
            fprintf("iteration No. %d ---> Variables:\n",obj.iteration)
            disp(desOrg);
            obj.LagrangianInfo.Lorg = [];
            obj.LagrangianInfo.dL_dx = [];
            obj.LagrangianInfo.alin = [];
            obj.LagrangianInfo.blin = [];  

            obj.Hessian = obj.settingUNGRADE.H0;

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
            if nargin<3
                meshDeformFlag = 0;
            end
            if meshDeformFlag == 1
                [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
                [orgMeshVerts,orgMeshCon] = obj.meshDeformation(obj,orgGeomVerts,orgGeomCon);
                surfID = obj.surfID;
                wakelineID = obj.wakelineID;
                fprintf("iteration No. %d set by Mesh Deformation ---> Variables:\n",obj.iteration)
                disp(unscaledVariables);
                obj = obj.setMesh(orgMeshVerts,orgMeshCon,surfID,wakelineID);
                obj.history.setMeshbyDeform(obj.iteration) = 1;
            else
                [orgMeshVerts, orgMeshCon,surfID,wakelineID, unscaledVariables] = obj.meshGenFun(unscaledVariables(:)');
                fprintf("iteration No. %d set by Mesh Generation ---> Variables:\n",obj.iteration)
                disp(unscaledVariables);
                obj = obj.setMesh(orgMeshVerts,orgMeshCon,surfID,wakelineID);
                [orgGeomVerts,orgGeomCon,optSREF,optBREF,optCREF,optXYZREF,obj.argin_x,desOrg] = obj.geomGenFun(unscaledVariables);
                obj.history.setMeshbyDeform(obj.iteration) = 0;
            end
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            obj.orgMesh = triangulation(obj.tri.ConnectivityList,obj.tri.Points);
            obj = obj.setREFS(optSREF,optBREF,optCREF);
            obj = obj.setRotationCenter(optXYZREF);
            obj.unscaledVar = desOrg(:)';
            obj.scaledVar = (desOrg(:)'-obj.lb)./obj.designScale;
            obj.orgGeom =  triangulation(orgGeomCon,orgGeomVerts);
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
            [udwf,udwr] = obj.calcDynCoefdu(flowNo,alpha,beta);
            [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,flowNo,alpha,beta);
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
                [udwf,udwr] = obj.calcDynCoefdu(iter,alphabuff,betabuff);
                [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,iter,alphabuff,betabuff);
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
                [udwf{iter},udwr{iter}] = obj.calcDynCoefdu(iter,alphabuff,betabuff);
                [AERODATA0,Cp0,Cfe0,R0solve,obj] = obj.solveFlowForAdjoint(u0solve,iter,alphabuff,betabuff);%ポテンシャルから空力係数を計算
                [DYNCOEF0,obj] =  obj.calcDynCoefforAdjoint(u0solve,udwf{iter},udwr{iter},iter,alphabuff,betabuff);
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
            fprintf("DYNCOEF of iteration No.%d ->\n",obj.iteration);
            disp(vertcat(DYNCOEF0{:}));
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
                    DYNCOEF =  obj.calcDynCoefforAdjoint(usolve,udwf{iter},udwr{iter},iter,alphabuff,betabuff);
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
                        dR_du = obj.LHS'; %亜音速
                    else
                        dR_du = eye(nbPanel);
                    end
                else
                    if obj.flowNoList(i,3) == 1
                        dR_du = blkdiag(dR_du,obj.LHS');
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
                    [udwfdx,udwrdx] = obj2.calcDynCoefdu(iter,alphabuff,betabuff);
                    DYNCOEF =  obj2.calcDynCoefforAdjoint(u0solve,udwfdx,udwrdx,iter,alphabuff,betabuff);
                    for k = 1:numel(lktable)
                        R((lktable(k)-1)*nbPanel+1:lktable(k)*nbPanel,1) = Rsolve(nbPanel*(k-1)+1:nbPanel*k,1);
                    end
                end
                
                [I,con] = objandConsFun(des,AERODATA,DYNCOEF,Cp,Cfe,SREF2,BREF2,CREF2,XYZREF2,argin_x2);
                dI_dx(i) = (I-I0)/pert;%評価関数の設計変数に関する偏微分
                if not(isempty(con))
                    dcon_dx(:,i) = (con-con0)/pert;
                end
                dR_dx(:,i) = (R-R0)./pert;
                clear obj2;
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
            gradientCalcdx = 1e-5;
            lbf = -obj.scaledVar;
            ubf = 1-obj.scaledVar;
            lactiveBound = obj.scaledVar'-gradientCalcdx<0;
            uactiveBound = obj.scaledVar'+gradientCalcdx>1;
            lbfGrad = -ones(size(lbf)).*gradientCalcdx;
            ubfGrad =  ones(size(lbf)).*gradientCalcdx;
            lbfGrad(lactiveBound) = lbf(lactiveBound);
            ubfGrad(uactiveBound) = ubf(uactiveBound);
            if not(isempty(con))
                obj.LagrangianInfo.alin = [-conTotalGrad;conTotalGrad;eye(numel(ubf));-eye(numel(lbf))];
                obj.LagrangianInfo.blin = [con0(:)+dcon_du*duVec-cmin(:);cmax(:)-con0(:)-dcon_du*duVec;ubf';-lbf'];
                blinGrad = [con0(:)+dcon_du*duVec-cmin(:);cmax(:)-con0(:)-dcon_du*duVec;ubfGrad';-lbfGrad'];
            else
                obj.LagrangianInfo.alin = [eye(numel(ubf));-eye(numel(lbf))];
                obj.LagrangianInfo.blin = [ubf';-lbf'];
                blinGrad = [ubfGrad';-lbfGrad'];
            end

            options = optimoptions(@linprog,'Algorithm','interior-point','Display','final','MaxIter',10000);
            [~,~,exitFlag,~,lambda] = linprog(objTotalGrad,obj.LagrangianInfo.alin,blinGrad,[],[],[],[],options);
            if exitFlag < 1
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed',"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                [~,~,~,~,lambda] = fmincon(@(x)obj.fminconObj(x,zeros(numel(ubf),numel(ubf)),objTotalGrad),zeros(numel(obj.scaledVar),1),obj.LagrangianInfo.alin,blinGrad,[],[],[],[],[],options);
            end
            rho = 1000;
            if not(isempty(con0))
                lambdaR = -(dR_du)\(dI_du+lambda.ineqlin(1:size(con,1)*2,1)'*[-dcon_du;dcon_du])';
                obj.LagrangianInfo.Lorg = I0 + lambda.ineqlin'*[-con0;con0;obj.scaledVar'.*uactiveBound;-obj.scaledVar'.*lactiveBound];
                penaltyorg = rho*sum(max(max(0,cmin-con0),max(con0-cmax)));
                obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+lambda.ineqlin'*[-dcon_dx;dcon_dx;diag(uactiveBound);-diag(lactiveBound)];
            else
                lambdaR = -(dR_du)\(dI_du)';
                obj.LagrangianInfo.Lorg = I0 + lambda.ineqlin'*[obj.scaledVar'.*uactiveBound;-obj.scaledVar'.*lactiveBound];
                penaltyorg = 0;
                obj.LagrangianInfo.dL_dx = dI_dx+lambdaR'*dR_dx+lambda.ineqlin'*[diag(uactiveBound);-diag(lactiveBound)];
            end
            obj.history.LagrangianVal(obj.iteration) = obj.LagrangianInfo.Lorg;
            obj.history.penaltyVal(obj.iteration) = penaltyorg;
            obj.history.scaledVar(obj.iteration,:) = obj.scaledVar;
            obj.history.dL_dx(obj.iteration,:) = obj.LagrangianInfo.dL_dx;
            fprintf("Lagrangian Value : %f \nPenalty Value : %f\nGradient of Lagrangian : \n",obj.LagrangianInfo.Lorg,penaltyorg);
            disp(obj.LagrangianInfo.dL_dx);
            fprintf("iteration No. %d : Gradient Caluculation Completed\n",obj.iteration);
        end

        function [nextUnscaledVar,scaleddxNorm] = descentLagrangian(obj,varargin)
            %%%%%%%%%%%%%Lagrangianの降下方向に向かう設計変数を計算する%%%%%%%%%
            %
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
            ndim = numel(lbf);
            if strcmpi(obj.settingUNGRADE.updateMethod,"Levenberg–Marquardt")
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                dxscaled = fmincon(@(x)obj.fminconObj(x,obj.Hessian+obj.settingUNGRADE.betaLM*diag(diag(obj.Hessian)),obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),obj.LagrangianInfo.alin,obj.LagrangianInfo.blin,[],[],lbf,ubf,@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion),options);
            elseif strcmpi(obj.settingUNGRADE.updateMethod,"dogleg")
                %sd方向を求める行列
                Aeq = [eye(ndim),obj.LagrangianInfo.dL_dx(:)];
                beq = zeros(ndim,1);
                options = optimoptions(@fmincon,'Algorithm','interior-point','Display','final-detailed','EnableFeasibilityMode',true,"SubproblemAlgorithm","cg",'MaxFunctionEvaluations',10000);
                xsd_t = fmincon(@(x)obj.fminconObj(x,obj.Hessian,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar)+1,1),[obj.LagrangianInfo.alin,zeros(size(obj.LagrangianInfo.alin,1),1)],obj.LagrangianInfo.blin,Aeq,beq,[lbf,-100],[ubf,100],[],options);
                xsd = xsd_t(1:ndim);
                xqn = fmincon(@(x)obj.fminconObj(x,obj.Hessian,obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),obj.LagrangianInfo.alin,obj.LagrangianInfo.blin,[],[],lbf,ubf,[],options);
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
                dxscaled = fmincon(@(x)obj.fminconObj(x,2.*eye(ndim),2.*obj.settingUNGRADE.alphaSD.*obj.LagrangianInfo.dL_dx),zeros(numel(obj.scaledVar),1),obj.LagrangianInfo.alin,obj.LagrangianInfo.blin,[],[],[],[],@(x)obj.fminconNlc(x,obj.settingUNGRADE.TrustRegion),options);
            end
            dx = dxscaled(:)'.*obj.designScale;
            nextUnscaledVar = desOrg + dx(:)';
            scaleddxNorm = norm(dxscaled);
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

        function [c,ceq] = fminconNlc(dx,TR)
            ceq = [];
            c = sum(dx.^2)-(TR)^2;
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
            md.x = scatteredInterpolant(obj.orgGeom.Points,modGeomVerts(:,1)-obj.orgGeom.Points(:,1),'linear','linear');
            md.y = scatteredInterpolant(obj.orgGeom.Points,modGeomVerts(:,2)-obj.orgGeom.Points(:,2),'linear','linear');
            md.z = scatteredInterpolant(obj.orgGeom.Points,modGeomVerts(:,3)-obj.orgGeom.Points(:,3),'linear','linear');
            dVerts(:,1) = md.x(obj.tri.Points);
            dVerts(:,2) = md.y(obj.tri.Points);
            dVerts(:,3) = md.z(obj.tri.Points);
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
                sampleSurff = dmodSurf(:);
                sampleDes = obj.scaledVar.*obj.designScale+obj.lb;
                sampleDes(i) = (obj.scaledVar(i) - pert(i)).*obj.designScale(i)+obj.lb(i);
                [modSurf,~,SREFr,BREFr,CREFr,XYZREFr,argin_xr,desBuff] = obj.geomGenFun(sampleDes);
                pertr = (desBuff(i)-desOrg(i))/obj.designScale(i);
                dmodSurf = modSurf-surforg;
                sampleSurfr = dmodSurf(:);
                gradSurf = (sampleSurff-sampleSurfr)./(pertf-pertr);
                obj.gradSREF(i) = (SREFf-SREFr)./(pertf-pertr);
                obj.gradBREF(i) = (BREFf-BREFr)./(pertf-pertr);
                obj.gradCREF(i) = (CREFf-CREFr)./(pertf-pertr);
                obj.gradXYZREF(:,i) = (XYZREFf(:)-XYZREFr(:))./(pertf-pertr);
                obj.gradArginx(:,i) = (argin_xf(:)-argin_xr(:))./(pertf-pertr);
                reshapeGrad = reshape(gradSurf,size(modSurf));
                md.x = scatteredInterpolant(obj.orgGeom.Points,reshapeGrad(:,1),'linear','linear');
                md.y = scatteredInterpolant(obj.orgGeom.Points,reshapeGrad(:,2),'linear','linear');
                md.z = scatteredInterpolant(obj.orgGeom.Points,reshapeGrad(:,3),'linear','linear');
                dVerts(:,1) = md.x(obj.orgMesh.Points);
                dVerts(:,2) = md.y(obj.orgMesh.Points);
                dVerts(:,3) = md.z(obj.orgMesh.Points);
                obj.gradMesh(:,i) = dVerts(:);
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
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
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
                if obj.paneltype(i) == 1
                    CPmat =obj.center(obj.cluster{i},1:3);
                    pnt = obj.center(i,:);
                    m = obj.tri.Points(obj.tri.ConnectivityList(i,1),:)'-pnt(:);
                    m = m./norm(m);
                    l = cross(m,obj.normal(i,:)');
                    Minv = [l,m,obj.normal(i,:)'];
                    lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                    bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                    Bmat=pinv(bb,sqrt(eps));
                    %Vnmat = Minv(:,[1,2])*[1,0,0,0;0,1,0,0]*Bmat;
                    Vnmat = Minv*[1,0,0,0;0,1,0,0;0,0,1,0;]*Bmat;
                    for iter = 1:3
                        approxmatedObj.mu2v{iter}(i,obj.IndexPanel2Solver(obj.cluster{i})) = Vnmat(iter,:);
                    end
                end
            end
            
                
            for wakeNo = 1:numel(obj.wakeline)
                theta = linspace(pi,0,size(obj.wakeline{wakeNo}.validedge,2)*obj.LLT.n_interp+1);
                obj.LLT.sp{wakeNo} = [];
                obj.LLT.calcMu{wakeNo} = zeros(1,nbPanel);
                s = zeros(1,numel(obj.wakeline{wakeNo}.edge));
                jter = 1;
                for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                    s(edgeNo+1) = s(edgeNo) + norm(obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),2:3)-obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),2:3));
                    if obj.wakeline{wakeNo}.valid(edgeNo) == 1
                        obj.LLT.sp{wakeNo} =  [obj.LLT.sp{wakeNo},(s(edgeNo)+s(edgeNo+1))./2];
                        jter = jter+1;
                    end
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
        %{

        function [errMax,errMean] = dist2geom(geom,mesh)
            numV = size(geom.Points,1);
            numF = size(geom.ConnectivityList,1);
            tV = size(mesh.Points,1);
            dMin = Inf(tV,1);
            for c1 = 1:tV 
                nn = geom.nearestNeighbor(mesh.Points(c1,:));
                dMin(c1) = norm(geom.Points(nn,:)-mesh.Points(c1,:));

                ne = geom.vertexAttachments(nn);
                for c2 = ne{1}
                    for c3 = 1:2 
                        for c4 = c3+1:3
                            v1 = mesh.Points(c1,:);
                            v2 = geom.Points(geom.ConnectivityList(c2,c3),:);
                            v3 = geom.Points(geom.ConnectivityList(c2,c4),:);
                            if ( all(min(v2,v3) < (v1 + dMin(c1))) && all(max(v2,v3) > (v1 - dMin(c1))) )
                                d = norm(cross((v3-v2),(v2-v1)))/norm(v3 - v2);
                                s = - (v2-v1)*(v3-v2)' / (norm(v3-v2))^2;
                                if (s>=0 && s<=1)
                                    if d < dMin(c1) 
                                        dMin(c1) = d;
                                    end
                                end
                            end
                        end
                    end
                end

                for c2 = ne{1}
                    v1 = mesh.Points(c1,:);
                    v2 = geom.Points(geom.ConnectivityList(c2,1),:);
                    v3 = geom.Points(geom.ConnectivityList(c2,2),:);
                    v4 = geom.Points(geom.ConnectivityList(c2,3),:);
                    if ( all(min([v2;v3;v4]) < (v1 + dMin(c1))) && all(max([v2;v3;v4]) > (v1 - dMin(c1))) )
                        n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
                        d = abs(n * (v1 - v2)');
                        n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
                        f1 = v1 + d * n;
                        f2 = v1 - d * n;
                        m = [v3-v2;v4-v2]';
                        if(abs(f1-v2)>0)
                            r1 = m\(f1-v2)';
                        else
                            r1 = [inf;inf];
                        end
                        if(abs(f2-v2)>0)
                            r2 = m\(f2-v2)';
                        else
                            r2 = [inf;inf];
                        end
                        if ((sum(r1)<=1 && sum(r1)>=0 && all(r1 >=0) && all(r1 <=1)) || (sum(r2)<=1 && sum(r2)>=0 && all(r2 >=0) && all(r2 <=1)))
                            if d < dMin(c1)
                                dMin(c1) = d;
                            end
                        end
                    end
                end
            end
            
            pVError = dMin; 
            errMean = mean(pVError);
            errMax = max(pVError);
        end

        function dist = edgeDist(v1,v2,v3)
            d = norm(cross((v3-v2),(v2-v1)))/norm(v3 - v2);
            s = - (v2-v1)*(v3-v2)' / (norm(v3-v2))^2;
            if (s>=0 && s<=1)
                dist = d;
            else
                dist = inf;
            end
        end

        function dist = faceDist(v1,v2,v3,v4)
            n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
            d = abs(n * (v1 - v2)');
            n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
            f1 = v1 + d * n;
            f2 = v1 - d * n;
            m = [v3-v2;v4-v2]';
            if(abs(f1-v2)>0)
                r1 = m\(f1-v2)';
            else
                r1 = [inf;inf];
            end
            if(abs(f2-v2)>0)
                r2 = m\(f2-v2)';
            else
                r2 = [inf;inf];
            end
            if ((sum(r1)<=1 && sum(r1)>=0 && all(r1 >=0) && all(r1 <=1)) || (sum(r2)<=1 && sum(r2)>=0 && all(r2 >=0) && all(r2 <=1)))
                dist = d;
            else
                dist = inf;
            end
        end
        %}
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

    end

end