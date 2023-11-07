classdef UNLSI

    properties
        tri %MATLAB triangulationクラス
        surfID %各パネルのタグ番号
        wakeline %wakeパネルの設定
        wakelineID %生データ
        halfmesh %半裁かどうか(半裁:1)
        SREF %基準面積
        CREF %縦方向基準長
        BREF %横方向基準長
        XYZREF %回転中心
        paneltype %各パネルのタイプ 1:body 2:base 3:structure
        cpcalctype %Cpを計算する方法の選択 %1:1-V^2 2: -2u 3: -2u-v^2-w^2
        IndexPanel2Solver %パネルのインデックス⇒ソルバー上でのインデックス
        VindWakes
        settingUNLSI
        orgNormal %各パネルの法線ベクトル
        modNormal %舵角等によって変更した法線ベクトル
        center %各パネルの中心
        area %各パネルの面積
        flowNoTable %flowNoが整理されたMatrix
        flow %各flowconditionが格納されたセル配列
        prop %各プロペラ方程式情報が格納されたセル配列
        cluster %パネルクラスターの情報
        LHS %パネル法連立方程式の左辺行列
        RHS %パネル法連立方程式の右辺行列
        mu2v %ポテンシャル⇒機体表面速度への変換行列
        Cp %圧力係数
        Cfe %表面摩擦係数
        deflAngle
        deflGroup
        AERODATA %空力解析結果の格納
        DYNCOEF %同安定微係数結果の格納
        LLT
        ppCoef
        ppDyn
    end

    methods(Access = public)
        function obj = UNLSI(verts,connectivity,surfID,wakelineID,halfmesh)
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            %%%%%%%%%%%Constructor%%%%%%%%%%%%%
            %verts:頂点座標
            %conectivity:各パネルの頂点ID
            %surfID:パネルのID
            %wakelineID:wakeをつけるエッジ上の頂点ID配列を要素にもつセル配列
            %halfmesh:半裁メッシュの場合に1を指定 引数を省略すれば自動設定
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %データの格納、パネル法線、面積の計算

            %halfmeshの判定
            if all(verts(:,2) >=0) ||  all(verts(:,2) <=0)
                if nargin > 4 
                    if halfmesh == 0
                        warning("Signatures of verts' y cordinates are all same, but the half mesh settingUNLSI is 0 (not halfmesh)");
                    end
                else
                    halfmesh = 1;
                end
            else
               if nargin > 4 
                    if halfmesh == 1
                        warning("Signatures of verts' y cordinates are NOT all same, but the half mesh settingUNLSI is 1 (halfmesh)");
                    end
                else
                    halfmesh = 0;
               end
            end
                     
            %%%%%%%%%%%settingUNLSIの初期設定%%%%%%%%%%%%%
            obj.settingUNLSI.checkMeshMethod = "delete";
            obj.settingUNLSI.checkMeshTol = sqrt(eps);
            obj.settingUNLSI.LLTnInterp = 10;
            obj.settingUNLSI.nCluster = 50;%クラスターの目標数（達成できなければそこで打ち切り）
            obj.settingUNLSI.edgeAngleThreshold = 50;%この角度以下であれば近隣パネルとして登録する（角度差が大きすぎるモノをはじくため）
            obj.settingUNLSI.wingWakeLength = 100; %各wakeパネルの長さ(機体基準長基準）
            obj.settingUNLSI.nCalcDivide = 5;%makeEquationは各パネル数×各パネル数サイズの行列を扱うため、莫大なメモリが必要となる。一度に計算する列をobj.settingUNLSI.nCalcDivide分割してメモリの消費量を抑えている。
            obj.settingUNLSI.angularVelocity = [];
            obj.settingUNLSI.propCalcFlag = 0;
            obj.settingUNLSI.deflDerivFlag = 1;
            obj.settingUNLSI.propWakeLength = 3;
            obj.settingUNLSI.nPropWake = 101; %propWake円周上の点の数
            obj.settingUNLSI.kCf = 1.015*(10^-5);
            obj.settingUNLSI.laminarRatio = 0;
            obj.settingUNLSI.coefCf = 1;
            obj.settingUNLSI.newtoniantype = "OldTangentCone";
            obj.settingUNLSI.resultSearchMethod = "and";
            obj.settingUNLSI.Vinf = 15;
            obj.settingUNLSI.rho = 1.225;
            obj.settingUNLSI.kappa = 1.4;
            obj.settingUNLSI.nGriddedInterp = 90;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.halfmesh = halfmesh;
            obj.flowNoTable = [];
            obj.flow = {};
            obj.deflAngle = [];
            obj.SREF = 1;
            obj.CREF = 1;
            obj.BREF = 1;
            obj.XYZREF = [0,0,0];

            obj = obj.setMesh(verts,connectivity,surfID,wakelineID);
        end

        function obj = setOptions(obj,varargin)
           %%%%%%%%%設定を変更する%%%%%%%%%%%%
           %varargin : 変数名と値のペアで入力
           %obj.settingの構造体の内部を変更する
           %構造体の内部変数の名前と一致していないとエラーが出るようになっている
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mod(numel(varargin),1) == 1
                error("input style not match");
            end
            for iter = 1:numel(varargin)/2
                if isfield(obj.settingUNLSI,varargin{2*iter-1})
                    obj.settingUNLSI = setfield(obj.settingUNLSI,varargin{2*iter-1},varargin{2*iter});
                else
                    error("Field name is not match");
                end
            end
        end

        function obj = setREFS(obj,SREF,BREF,CREF)
            %%%%%%%%%%%Referentials settingUNLSI%%%%%%%%%%%%%
            %SREF:基準面積
            %BREF:横方向基準長(eg.翼幅)
            %CREF:縦方向基準長(eg.MAC,機体全長
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.SREF = SREF;
            obj.CREF = CREF;
            obj.BREF = BREF;
        end

        function obj = setRotationCenter(obj,XYZREF)
            %%%%%%%%%%%Rotation Center settingUNLSI%%%%%%%%%%%%%
            %回転中心を設定する。
            %回転中心：モーメント計算の基準位置,主流角速度の回転中心
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.XYZREF = XYZREF(:)';
        end

        function obj = setDeflGroup(obj,groupNo,groupName,groupID,deflGain)
            obj.deflGroup{groupNo}.name = groupName;
            obj.deflGroup{groupNo}.ID = groupID;
            obj.deflGroup{groupNo}.gain = deflGain;
        end

        function obj = setDeflAngle(obj,ID,rotAxis,angle)
            %%%%%%%%%%%Rotation Center settingUNLSI%%%%%%%%%%%%%
            %疑似的な舵角を設定する。
            %そのID上のパネルの法線ベクトルをロドリゲスの回転ベクトルによって回転する
            % 誤差については要検討。
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for iter = 1:numel(ID)
                i = find(obj.deflAngle(:,1)==ID(iter));
                if isempty(i)
                    error("Invalid ID");
                end
                dcm = obj.rod2dcm(rotAxis(iter,:),angle(iter));
                obj.deflAngle(i,2:end) = [rotAxis(iter,:),angle(iter),dcm(:)'];
            end
            for i = 1:numel(obj.paneltype)
                %[obj.area(i,1),~ , obj.orgNormal(i,:)] = obj.vertex(obj.tri.Points(obj.tri.ConnectivityList(i,1),:),obj.tri.Points(obj.tri.ConnectivityList(i,2),:),obj.tri.Points(obj.tri.ConnectivityList(i,3),:));
                %obj.center(i,:) = [mean(obj.tri.Points(obj.tri.ConnectivityList(i,:),1)),mean(obj.tri.Points(obj.tri.ConnectivityList(i,:),2)),mean(obj.tri.Points(obj.tri.ConnectivityList(i,:),3))];
                obj.modNormal(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])*obj.orgNormal(i,:)')';
            end
        end

        function plotGeometry(obj,figureNo,triColor,colorlim,method,extrapmethod)
            %%%%%%%%%%%plotting%%%%%%%%%%%%%
            %機体メッシュもしくは状態量をプロットする。
            %カラー値は各パネルごとの値（例えばCp）を入れると、sccatteredInterpを用いて接点の値に補間しプロットする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(figureNo);clf;
            if nargin<3
                trisurf(obj.tri,"FaceAlpha",0);
                if obj.halfmesh == 1
                    hold on
                    trisurf(obj.tri.ConnectivityList,obj.tri.Points(:,1),-obj.tri.Points(:,2),obj.tri.Points(:,3),"FaceAlpha",0);
                    hold off;
                end
            else
                if nargin == 4
                    method = "linear";
                    extrapmethod = "linear";
                end
                if strcmpi(method,"exact")
                    %指定がexactの場合はパネル一枚一枚描画する
                    hold on;
                    for i = 1:numel(obj.surfID)
                        if obj.paneltype(i) == 1
                            trisurf([1,2,3],obj.tri.Points(obj.tri(i,:),1),obj.tri.Points(obj.tri(i,:),2),obj.tri.Points(obj.tri(i,:),3),triColor(i),'EdgeAlpha',0.15);
                            if obj.halfmesh == 1
                                trisurf([1,2,3],obj.tri.Points(obj.tri(i,:),1),-obj.tri.Points(obj.tri(i,:),2),obj.tri.Points(obj.tri(i,:),3),triColor(i),'EdgeAlpha',0.15);
                            end
                        else
                            trisurf([1,2,3],obj.tri.Points(obj.tri(i,:),1),obj.tri.Points(obj.tri(i,:),2),obj.tri.Points(obj.tri(i,:),3),0,'EdgeAlpha',0.15);
                            if obj.halfmesh == 1
                                trisurf([1,2,3],obj.tri.Points(obj.tri(i,:),1),-obj.tri.Points(obj.tri(i,:),2),obj.tri.Points(obj.tri(i,:),3),0,'EdgeAlpha',0.15);
                            end
                        end
                    end
                    colormap jet;
                    hold off;
                    caxis(colorlim);
                else
                    F = scatteredInterpolant(obj.center,triColor,method,extrapmethod);
                    c = F(obj.tri.Points);
                    trisurf(obj.tri,c,'FaceColor','interp','EdgeAlpha',0.15);
                    colormap jet;
                    if obj.halfmesh == 1
                        hold on;
                        trisurf(obj.tri.ConnectivityList,obj.tri.Points(:,1),-obj.tri.Points(:,2),obj.tri.Points(:,3),c,'FaceColor','interp','EdgeAlpha',0.15);
                        colormap jet;
                        hold off
                    end
                    caxis(colorlim);
                end
            end
            axis equal;xlabel("x");ylabel("y");zlabel("z");
            drawnow();pause(0.1);

        end

        function obj = setMesh(obj,verts,connectivity,surfID,wakelineID)
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            obj.tri = triangulation(connectivity,verts);
            obj.surfID = surfID;
            obj.wakelineID = wakelineID;
            obj.wakeline = [];
            for i = 1:numel(wakelineID)
                obj.wakeline{i}.edge = wakelineID{i};
            end
            if isempty(obj.deflAngle)
                eyebuff = eye(3);
                obj.deflAngle = [unique(obj.surfID),zeros(numel(unique(obj.surfID)),4),repmat(eyebuff(:)',[numel(unique(obj.surfID)),1])];
            end
            obj.paneltype = ones(size(connectivity,1),1);
            obj.cpcalctype = ones(size(connectivity,1),1);
            obj.IndexPanel2Solver = 1:numel(obj.paneltype);
            obj.cluster = cell([1,numel(obj.paneltype)]);
            obj.area = zeros(numel(obj.paneltype),1);
            obj.orgNormal = zeros(numel(obj.paneltype),3);
            obj.modNormal = zeros(numel(obj.paneltype),3);
            obj.center = zeros(numel(obj.paneltype),3);
            obj.Cp = {};
            obj.Cfe = {};
            obj.LHS = [];
            obj.RHS = [];
            lltninterp = obj.settingUNLSI.LLTnInterp;
            obj.LLT = [];
            obj.settingUNLSI.LLTnInterp = lltninterp;
            for i = 1:numel(obj.paneltype)
                [obj.area(i,1),~ , obj.orgNormal(i,:)] = obj.vertex(verts(connectivity(i,1),:),verts(connectivity(i,2),:),verts(connectivity(i,3),:));
                obj.center(i,:) = [mean(verts(obj.tri.ConnectivityList(i,:),1)),mean(verts(obj.tri.ConnectivityList(i,:),2)),mean(verts(obj.tri.ConnectivityList(i,:),3))];
                obj.modNormal(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])*obj.orgNormal(i,:)')';
            end

            
            %半裁メッシュの境界表面上のwakeを削除
            deleteiter = [];
            for iter = 1:numel(obj.wakeline)
                deletelist = [];
                for j = 1:numel(obj.wakeline{iter}.edge)-1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{iter}.edge(j),obj.wakeline{iter}.edge(j+1));
                    if numel(attachpanel{1}) < 2
                        deletelist = [deletelist,j];
                    end
                end
                obj.wakeline{iter}.edge(deletelist) = [];
                if numel(obj.wakeline{iter}.edge) < 3
                    deleteiter = [deleteiter,iter];
                end
            end
            obj.wakeline(deleteiter) = [];
            %wakeのつくパネルIDを特定
            wakesurfID = [];
            for wakeNo = 1:numel(obj.wakeline)
                panelNo = 1;
                obj.wakeline{wakeNo}.valid = zeros(1,numel(obj.wakeline{wakeNo}.edge)-1);
                for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{wakeNo}.edge(edgeNo),obj.wakeline{wakeNo}.edge(edgeNo+1));
                    wakesurfID = [wakesurfID,setdiff(obj.surfID(attachpanel{1}),wakesurfID)];
                    if numel(attachpanel{1})==2
                        obj.wakeline{wakeNo}.validedge(1,panelNo) = obj.wakeline{wakeNo}.edge(edgeNo);
                        obj.wakeline{wakeNo}.validedge(2,panelNo) = obj.wakeline{wakeNo}.edge(edgeNo+1);
                        obj.wakeline{wakeNo}.valid(edgeNo) = 1;
                        phivert = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),2:3)-obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),2:3);
                        phiwake = atan2d(phivert(2),phivert(1));
                        if abs(phiwake)<90
                            obj.wakeline{wakeNo}.upperID(panelNo) = attachpanel{1}(1);
                            obj.wakeline{wakeNo}.lowerID(panelNo) = attachpanel{1}(2);
                        else
                            obj.wakeline{wakeNo}.upperID(panelNo) = attachpanel{1}(2);
                            obj.wakeline{wakeNo}.lowerID(panelNo) = attachpanel{1}(1);
                        end
                        panelNo = panelNo + 1;
                    else
                        warning("invalid wake edge detected");
                    end
                end
            end
            if numel(obj.prop)>0
                for propNo = 1:numel(obj.prop)
                    obj = obj.setProp(propNo,obj.prop{propNo}.ID,obj.prop{propNo}.diameter,obj.prop{propNo}.XZsliced);
                end
            end
            %wakeSurfIDにあるもの以外はcpcalctypeをlinearに
            for i = setdiff(unique(obj.surfID)',wakesurfID)
                obj = obj.setCpCalcType(i,"linear");
            end
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol,obj.settingUNLSI.checkMeshMethod);
            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        end

        function obj = setVerts(obj,verts)
            %%%%%%%%%%%%%%%%%接点座標の変更%%%%%%%%%%%%%%%
            %三角形の接続関係を変えずに、接点の座標のみ変更する
            %変更時したときに代わる値は全てここで計算しなおす。
            %最適化を行うときに、設計変数の勾配を計算するとき等に使用する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %triangulationはreadonlyなので、新しく作り直す
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            con = obj.tri.ConnectivityList;
            obj.tri = triangulation(con,verts);
            
            for i = 1:numel(obj.paneltype)
                [obj.area(i,1),~ , obj.orgNormal(i,:)] = obj.vertex(verts(obj.tri(i,1),:),verts(obj.tri(i,2),:),verts(obj.tri(i,3),:));
                obj.center(i,:) = [mean(verts(obj.tri.ConnectivityList(i,:),1)),mean(verts(obj.tri.ConnectivityList(i,:),2)),mean(verts(obj.tri.ConnectivityList(i,:),3))];
                obj.modNormal(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])*obj.orgNormal(i,:)')';
            end
            obj.checkMesh(obj.settingUNLSI.checkMeshTol,obj.settingUNLSI.checkMeshMethod);
            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        end

        function obj = checkMesh(obj,tol,outType)
            %%%%%%%%%%%%%%%%パネル面積の確認%%%%%%%%%%%%%%%
            %tol:最小許容面積
            %outType:"warning"を指定するとエラーの代わりに警告を表示する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                outType = "error";
            end
            if any(obj.area<tol)
                if strcmpi(outType,'warning')
                    warning("some panel areas are too small");
                elseif strcmpi(outType,'delete')
                    
                    deleteIndex = obj.area<tol;
                    if any(deleteIndex)
                        warning("some panel areas are too small and deleted");
                    end
                    newCon = obj.tri.ConnectivityList;
                    newCon(deleteIndex,:) = [];
                    obj.surfID(deleteIndex,:) = [];
                    obj = obj.setMesh(obj.tri.Points,newCon,obj.surfID,obj.wakelineID,"warning",0);
                else
                    error("some panel areas are too small");
                end
            end
        end

        function obj = mergeVerts(obj,tol)
            %%%%%%%%%%%%%%%%近接しているvertsをマージする%%%%%%%%%%%%%%%
            %%%%%%%%%%%%CAUTISON BUGGY%%%%%%%%%%%%%%%%
            %まだきちんとできていない
            %tol:マージする最小距離
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nVerts = size(obj.tri.Points,1);
            mergePair = [];
            newCon = obj.tri.ConnectivityList;
            newVerts = obj.tri.Points;
            for i = 1:nVerts-1
                sindex = i+1:nVerts;
                nDist = vecnorm(obj.tri.Points(sindex,:)-repmat(obj.tri.Points(i,:),[numel(sindex),1]),2,2);
                if any(nDist<tol)
                    mergePoint = find(nDist<0.01)+i;
                    mergePair = [mergePair;repmat(i,numel(mergePoint),1),mergePoint];
                end
            end
            for i = 1:size(mergePair,1)
                %ID1 = obj.tri.vertexAttachments(mergePair(i,1));
                ID2 = obj.tri.vertexAttachments(mergePair(i,2));
                for j = 1:numel(ID2{1})
                    if numel(setdiff(obj.tri(ID2{1}(j),:),mergePair(i,:))) > 1
                        newCon(ID2{1}(j),obj.tri(ID2{1}(j),:)==mergePair(i,2)) = mergePair(i,1);
                    end
                end
                newVerts(mergePair(i,1),:) = (obj.tri.Points(mergePair(i,1),:)+obj.tri.Points(mergePair(i,2),:))./2;
            end
            for i = 1:numel(obj.wakeline)
                wakelineID{i} = obj.wakeline{i}.edge;
            end
            obj = obj.setMesh(newVerts,newCon,obj.surfID,wakelineID);
        end
        
        function obj = setPanelType(obj,ID,typename)
            %%%%%%%%%%%%%%%%paneltypeの指定関数%%%%%%%%%%%%%%%%
            %panelname : body 1 機体表面パネル
            %panelname : base 2 ベース面パネル
            %panelname : structure 3 構造パネル 
            %panelname : prop 4 アクチュエータディスク（プロペラ）パネル 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmpi(typename,"body")
                obj.paneltype(obj.surfID==ID,1) = 1;
            elseif strcmpi(typename,"base")
                obj.paneltype(obj.surfID==ID,1) = 2;
            elseif strcmpi(typename,"structure")
                obj.paneltype(obj.surfID==ID,1) = 3;
            elseif strcmpi(typename,"prop")
                obj.paneltype(obj.surfID==ID,1) = 4;  
            end
            index = 1;
            %パネルタイプによってはパネル法連立方程式に含まれないので、連立方程式上でのインデックスと全体のインデックスを関連付ける
            obj.IndexPanel2Solver = zeros(numel(obj.paneltype),1);
            for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1
                    obj.IndexPanel2Solver(i) = index;
                    index = index+1;
                end
            end
        end

        function obj = setBasewithAngle(obj,alpha,beta,threshold,ID)
            %%%%%%%%%%%%%%%%角度によるベース面指定%%%%%%%%%%%%%%%%
            %alphabetaで指定した角度に対してthreshold(deg)角度以下のパネルをベース面として設定する
            %alpha,beta : AoA Sideslip(deg)
            %threshold : deg 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flowVec(1) = cosd(alpha)*cosd(beta);
            flowVec(2) = -sind(beta);
            flowVec(3) = sind(alpha)*cosd(beta);
            for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1 && obj.surfID(i) == ID
                    if(abs(acosd(dot(flowVec,obj.orgNormal(i,:))))<threshold)
                        obj.paneltype(i,1) = 2;
                    end
                end
            end
            index = 1;
            %パネルタイプによってはパネル法連立方程式に含まれないので、連立方程式上でのインデックスと全体のインデックスを関連付ける
            obj.IndexPanel2Solver = zeros(numel(obj.paneltype),1);
            for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1
                    obj.IndexPanel2Solver(i) = index;
                    index = index+1;
                end
            end
        end

        function obj = setCpCalcType(obj,ID,typename)
            %%%%%%%%%%%%%%%%Cpの計算方法の指定関数%%%%%%%%%%%%%%%%
            %panelname : incompressible 1 一般的な圧力推算 1-V^2
            %panelname : linear         2 線形化 -2u
            %panelname : slender body   3 -2u-w^2-v^2 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmpi(typename,"incompressible")
                obj.cpcalctype(obj.surfID==ID,1) = 1;
            elseif strcmpi(typename,"linear")
                obj.cpcalctype(obj.surfID==ID,1) = 2;
            elseif strcmpi(typename,"slender-body")
                obj.cpcalctype(obj.surfID==ID,1) = 3;
            else
                error("invalid calc type name")
            end
        end

        function obj = makeCluster(obj)
            %%%%%%%%%%%%%パネルクラスターの作成%%%%%%%%%%%%%%%%%%
            %各パネルの近隣パネルを指定の数集める
            %クラスター内のIDは統一（TODO：統一しないオプションをつけるか検討）
            %パネル法連立方程式を解いて得られるポテンシャルから機体表面に沿った微分を計算するための準備
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1 %表面パネルであれば
                    obj.cluster{i} = i;
                    index = 1;
                    while(1)
                        neighbor = obj.tri.neighbors(obj.cluster{i}(index));
                        neighbor(isnan(neighbor)) = [];
                        %角度が厳しいところは削除
                        deletelist = [];
                        for j = 1:numel(neighbor)
                            edgeAngle = acos(dot(obj.orgNormal(obj.cluster{i}(index),:),obj.orgNormal(neighbor(j),:)));
                            if abs(edgeAngle)> obj.settingUNLSI.edgeAngleThreshold*pi/180 || obj.paneltype(neighbor(j)) ~= 1
                                deletelist = [deletelist,j];
                            end
                        end
                        neighbor(deletelist) = [];
                        diffcluster = setdiff(neighbor,obj.cluster{i});

                        if numel(obj.cluster{i})> obj.settingUNLSI.nCluster
                            break;
                        end
                        obj.cluster{i} = [obj.cluster{i},diffcluster];
                        index = index + 1;
                        if numel(obj.cluster{i})<index
                            break;
                        end
                    end
                end
           end
        end

        function obj = makeEquation(obj)
            %%%%%%%%%%%%%パネル法連立方程式行列の作成%%%%%%%%%%%%
            %パネル法の根幹
            %表面微分行列も併せて作成している
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            %パネル微分行列の作成
            obj.mu2v{1} = sparse(nPanel,nbPanel);
            obj.mu2v{2} = sparse(nPanel,nbPanel);
            obj.mu2v{3} = sparse(nPanel,nbPanel);
            
            for i = 1:nPanel
                if obj.paneltype(i) == 1
                    CPmat =obj.center(obj.cluster{i},1:3);
                    pnt = obj.center(i,:);
                    m = obj.tri.Points(obj.tri.ConnectivityList(i,1),:)'-pnt(:);
                    m = m./norm(m);
                    l = cross(m,obj.orgNormal(i,:)');
                    Minv = [l,m,obj.orgNormal(i,:)'];
                    lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                    bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                    Bmat=pinv(bb,sqrt(eps));
                    %Vnmat = Minv(:,[1,2])*[1,0,0,0;0,1,0,0]*Bmat;
                    Vnmat = Minv*[1,0,0,0;0,1,0,0;0,0,1,0;]*Bmat;
                    for iter = 1:3
                        obj.mu2v{iter}(i,obj.IndexPanel2Solver(obj.cluster{i})) = Vnmat(iter,:);
                    end
                end
            end

            %誘導抗力計算用の行列の作成
            
            %パネル方連立方程式行列の作成
            %機体パネル⇒機体パネルへの影響
            
            si = floor(nbPanel/obj.settingUNLSI.nCalcDivide).*(0:obj.settingUNLSI.nCalcDivide-1)+1;
            ei = [floor(nbPanel/obj.settingUNLSI.nCalcDivide).*(1:obj.settingUNLSI.nCalcDivide-1),nbPanel];
            obj.LHS = zeros(nbPanel);
            obj.RHS = zeros(nbPanel);

            for i= 1:obj.settingUNLSI.nCalcDivide
                [~,~,VortexAc,VortexBc] = obj.influenceMatrix(obj,[],si(i):ei(i));
                obj.LHS(:,si(i):ei(i)) = VortexAc; 
                obj.RHS(:,si(i):ei(i)) = VortexBc;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %wakeパネル⇒機体パネルへの影響
            %
            
            for wakeNo = 1:numel(obj.wakeline)
                theta = linspace(pi,0,size(obj.wakeline{wakeNo}.validedge,2)*obj.settingUNLSI.LLTnInterp+1);
                iter = 1;
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
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                    interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                    [influence] = obj.wakeInfluenceMatrix(obj,wakeNo,edgeNo,1:nbPanel);
                    obj.LHS(:,interpID(1)) = obj.LHS(:,interpID(1)) - influence;
                    obj.LHS(:,interpID(2)) = obj.LHS(:,interpID(2)) + influence;
                    obj.LLT.calcMu{wakeNo}(iter,interpID(1)) = -1;
                    obj.LLT.calcMu{wakeNo}(iter,interpID(2)) = 1;
                    iter = iter+1;
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
            if numel(obj.wakeline)>0
                obj.LLT.Qij = obj.Calc_Q(horzcat(obj.LLT.yinterp{:}),horzcat(obj.LLT.zinterp{:}),horzcat(obj.LLT.phiinterp{:}),horzcat(obj.LLT.spanel{:}),obj.halfmesh);
            end
        end

        function obj = setProp(obj,propNo,ID,diameter,XZsliced)
            %IDのパネルタイプをプロペラ==4に変更
            obj = obj.setPanelType(ID,'prop');
            obj.prop{propNo}.ID = ID;
            nPanel = numel(obj.paneltype);
            %プロペラディスク⇒bodyパネルへの影響係数
            obj.prop{propNo}.diameter = diameter;
            obj.prop{propNo}.area = pi*(obj.prop{propNo}.diameter/2)^2;
            obj.prop{propNo}.XZsliced  = XZsliced;
            %ペラパネルの重心位置を計算
            mom = [0,0,0];
            propArea = 0;
            for i = 1:nPanel
                if obj.surfID(i) == ID
                    darea = obj.vertex(obj.tri.Points(obj.tri(i,1),:),obj.tri.Points(obj.tri(i,2),:),obj.tri.Points(obj.tri(i,3),:));
                    propArea = propArea + darea;
                    mom = mom + darea.*obj.center(i,:);
                    if XZsliced == 1
                        mom = mom + darea.*[obj.center(i,1),-obj.center(i,2),obj.center(i,3)];
                    end
                end
            end
            if XZsliced == 1
                propArea = 2*propArea;
            end
            obj.prop{propNo}.normal = mean(obj.orgNormal(obj.surfID == ID,:),1);
            obj.prop{propNo}.center = mom./propArea;

            obj = obj.setPropState(propNo,0,0,0);

            obj = makePropEquation(obj,propNo);


        end

        function obj = makePropEquation(obj,propNo)
            %とりあえずpaneltypeで計算しているが、IDで管理したい
            %muとsigma両方
            [VortexAi,VortexBi,VortexAo,VortexBo,propWakeA] = obj.propellerInfluenceMatrix(obj,propNo,obj.settingUNLSI.propWakeLength);
            obj.prop{propNo}.LHSi = VortexAi;
            obj.prop{propNo}.RHSi = VortexBi;
            obj.prop{propNo}.LHSo = VortexAo;
            obj.prop{propNo}.RHSo = VortexBo;
            obj.prop{propNo}.wakeLHS = propWakeA;
        end

        function obj = flowCondition(obj,flowNo,Mach,Re)
            %%%%%%%%%%%%%%%%主流の設定%%%%%%%%%%%%%%%%%%%%%%
            %UNLSIでは流れの状態の変化によって不連続的な変化を起こさないよう（設計変数の微分が連続になるよう）考慮されている。
            %特に超音速流解析では主流の状態によって各パネルの膨張と圧縮が切り替わるため、特別な考慮が必要となる。
            %UNLSIでは各マッハ数に応じてパネル角度-180deg~180degにて事前に圧力係数を計算し、実際の計算ではgriddedInterpolantを用いて値を求めている。
            %したがって、主流状態と補間関数を事前につくる必要がある。
            
            %flowNo:作成する流れのID
            %Mach:マッハ数
            %超音速解析(修正ニュートン流理論)の手法:"OldTangentCone"(デフォルト),"TangentConeEdwards","TangentWedge"
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.flow{flowNo}.Mach = Mach;
            if Mach < 1
                obj.flow{flowNo}.flowtype = "Souce-Doublet 1st-order Panel";
            else
                obj.flow{flowNo}.flowtype = strcat(obj.settingUNLSI.newtoniantype,"+Prandtl-Meyer");
            end

            if obj.flow{flowNo}.Mach > 1
                delta = linspace(-pi,pi,500);
                Cpfc = zeros(size(delta));
                for j =1:size(delta,2)
                    if delta(j) >= 0
                        if strcmpi(obj.settingUNLSI.newtoniantype,'TangentConeEdwards')
                            Mas = (0.87*Mach-0.554)*sin(delta(j))+0.53;
                            Cpfc(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(obj.settingUNLSI.newtoniantype,'OldTangentCone')
                            Mas = 1.090909*Mach*sin(delta(j))+exp(-1.090909*Mach*sin(delta(j)));
                            Cpfc(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(obj.settingUNLSI.newtoniantype,'TangentWedge')
                            if delta(j)>45.585*pi/180
                                Cpfc(j)=((1.2*Mach*sin(delta(j))+exp(-0.6*Mach*sin(delta(j))))^2-1.0)/(0.6*Mach^2);
                                %R=1/Mach^2+(obj.settingUNLSI.kappa+1)/2*delta(j)/sqrt(Mach^2-1);
                            elseif delta(j)<0.035
                                Cpfc(j)=obj.settingUNLSI.kappa*Mach^2*delta(j)/sqrt(Mach^4-1.0);
                            else
                                b = -((Mach^2+2)/Mach^2)-obj.settingUNLSI.kappa*sin(delta(j))^2;
                                c = (2*Mach^2+1)/Mach^4+((obj.settingUNLSI.kappa+1)^2/4+(obj.settingUNLSI.kappa-1)/Mach^2)*sin(delta(j))^2;
                                d = -cos(delta(j))^2/Mach^4;
                                q = (b^2-3*c)/9;
                                rr = (b*(2*b^2-9*c)+27*d)/54;
                                disc = q^3-rr^2;
                                if disc<0
                                    Cpfc(j)=((1.2*Mach*sin(delta(j))+exp(-0.6*Mach*sin(delta(j))))^2-1.0)/(0.6*Mach^2);
                                else
                                    r = roots([1,b,c,d]);
                                    %ts = asin(sqrt(r(2)));
                                    R=r(2);
                                    Cpfc(j) = 4*(Mach^2*R-1)/((obj.settingUNLSI.kappa+1)*Mach^2);
                                end
                            end
                        end
                        Pe_Pinf= Cpfc(j)*(obj.settingUNLSI.kappa*Mach^2)/2+1;
                        Te_Tinf = Pe_Pinf^(1/(obj.settingUNLSI.kappa/(obj.settingUNLSI.kappa-1)));
                        Me = sqrt((2+(obj.settingUNLSI.kappa+1)*Mach^2)/((obj.settingUNLSI.kappa+1)*Te_Tinf)-2/(obj.settingUNLSI.kappa+1));
                    else
                        Me = UNLSI.bisection(@(M2)UNLSI.pmsolve(Mach,M2,-delta(j),obj.settingUNLSI.kappa),Mach,300);
                        Te_Tinf = (2+(obj.settingUNLSI.kappa+1)*Mach^2)/(2+(obj.settingUNLSI.kappa+1)*Me^2);
                        Pe_Pinf=(Te_Tinf)^(obj.settingUNLSI.kappa/(obj.settingUNLSI.kappa-1));
                        Cpfc(j) = 2/(obj.settingUNLSI.kappa*Mach^2)*(Pe_Pinf-1);
                    end
                end
                obj.flow{flowNo}.pp = griddedInterpolant(delta,Cpfc,'spline');
            end
            if nargin == 3
                Re = 0;
            else
                obj = obj.setCf(flowNo,Re);
            end
            obj.flowNoTable(flowNo,:) = [Mach,Re]; 
        end

        function obj = setCf(obj,flowNo,Re)
            %%%%%%%%%%%%%%%%%%%%%Cf計算における設定%%%%%%%%%%%%%%%%%%%
            %Cf計算の詳細はhttps://openvsp.org/wiki/doku.php?id=parasitedragが詳しい
            %flowNo:流れのID
            %Re:レイノルズ数
            %Lch:代表長さ
            %k:skin roughness value
                %Camouflage paint on aluminium
                %k = 1.015*(10^-5); 
                %Smooth paint
                %k = 0.634*(10^-5);
                %Produciton sheet metal
                %k = 0.405*(10^-5);
                %Polished sheet metal
                %k = 0.152*(10^-5);
                %Smooth molded composite
                %k = 0.052*(10^-5);
            %obj.settingUNLSI.laminarRatio:層流の割合
            %obj.settingUNLSI.coefCf:調整用の係数(デフォルト:1)

            %TODO:パネルIDによってCfの値を切り替えられるようにする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{

            %}
            obj.flowNoTable(flowNo,:) = [obj.flow{flowNo}.Mach,Re]; 
            Lch = obj.CREF;
            obj.flow{flowNo}.Re = Re;
            if(obj.flow{flowNo}.Mach<0.9)
                Re_Cut = 38.21*((Lch/obj.settingUNLSI.kCf)^1.053);
            else
                Re_Cut = 44.62*((Lch/obj.settingUNLSI.kCf)^1.053).*((obj.flow{flowNo}.Mach)^1.16);
            end
            if Re>Re_Cut
                Re = Re_Cut;
            end
            Cf_L = 1.328./sqrt(Re);
            Cf_T = 0.455/((log10(Re).^(2.58))*((1+0.144*obj.flow{flowNo}.Mach*obj.flow{flowNo}.Mach)^0.65));
            obj.Cfe{flowNo} = zeros(numel(obj.paneltype),1);
            obj.Cfe{flowNo}(obj.paneltype==1,1) =(obj.settingUNLSI.laminarRatio*Cf_L + (1-obj.settingUNLSI.laminarRatio)*Cf_T)*obj.settingUNLSI.coefCf;
        end

        function [obj,thrust,power,Jref] = setPropState(obj,propNo,CT,CP,rpm)

            obj.prop{propNo}.CT = CT;
            obj.prop{propNo}.CP = CP;
            obj.prop{propNo}.rpm = rpm;
            thrust(propNo) = obj.prop{propNo}.CT * obj.settingUNLSI.rho * (obj.prop{propNo}.rpm/60)^2 * obj.prop{propNo}.diameter^4;
            power(propNo) =  obj.prop{propNo}.CP * obj.settingUNLSI.rho * (obj.prop{propNo}.rpm/60)^3 * obj.prop{propNo}.diameter^5;
            obj.prop{propNo}.thrust = thrust(propNo);
            obj.prop{propNo}.power = power(propNo);
            Jref(propNo) =  obj.settingUNLSI.Vinf / ( 2*obj.prop{propNo}.rpm*obj.prop{propNo}.diameter/2/60);
            if nargout>1
                fprintf("Prop No %d (ID;%d)\n",propNo,obj.prop{propNo}.ID);
                fprintf("Prop Thrust\n");
                disp(thrust);
                fprintf("Prop Power\n");
                disp(power);
                fprintf("Prop Advance Ratio\n");
                disp(Jref);
                fprintf("Prop Ref Efficiency\n");
                disp(Jref.*obj.prop{propNo}.CT./obj.prop{propNo}.CP);
            end
        end

        function obj = solveFlow(obj,alpha,beta,Mach,Re)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %結果はobj.AERODATAに格納される。
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt 7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14:E(翼効率) 15:CFx 16:CFy 17:CFz 18:CMx 19:CMy 20:CMz 21:CMl 22:CMm 23:CMn 24:FOpt 
            %上記で求めていないものは0が代入される
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %omega:主流の回転角速度(deg/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.flowNoTable)
                if nargin == 5
                    obj = obj.flowCondition(1,Mach,Re);
                else
                    obj = obj.flowCondition(1,Mach);
                end
                flowNo = 1;
            else
                flowNo = 0;
                for i = 1:size(obj.flowNoTable,1)
                    if obj.flowNoTable(i,1) == Mach && obj.flowNoTable(i,2) == Re
                        flowNo = i;
                    end
                end
                if flowNo == 0
                    if nargin == 5
                        obj = obj.flowCondition(size(obj.flowNoTable,1)+1,Mach,Re);
                    else
                        obj = obj.flowCondition(size(obj.flowNoTable,1)+1,Mach);
                    end
                    flowNo = size(obj.flowNoTable,1);
                end
            end
            if any(size(alpha) ~= size(beta))
                error("analysis points are not match");
            end
            obj.AERODATA{flowNo} = [];
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            obj.Cp{flowNo} = zeros(nPanel,numel(alpha));
            for iterflow = 1:numel(alpha)
                T(1,1) = cosd(alpha(iterflow))*cosd(beta(iterflow));
                T(1,2) = cosd(alpha(iterflow))*sind(beta(iterflow));
                T(1,3) = -sind(alpha(iterflow));
                T(2,1) = -sind(beta(iterflow));
                T(2,2) = cosd(beta(iterflow));
                T(2,3) = 0;
                T(3,1) = sind(alpha(iterflow))*cosd(beta(iterflow));
                T(3,2) = sind(alpha(iterflow))*sind(beta(iterflow));
                T(3,3) = cosd(alpha(iterflow));

                if isempty(obj.settingUNLSI.angularVelocity)
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0]))';
                    end
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.settingUNLSI.angularVelocity(iterflow,:)./180.*pi,rvec(:)')'))';
                    end
                end
                
                Tvec(:,1) = obj.modNormal(:,2).* Vinf(:,3)-obj.modNormal(:,3).* Vinf(:,2);
                Tvec(:,2) = obj.modNormal(:,3).* Vinf(:,1)-obj.modNormal(:,1).* Vinf(:,3);
                Tvec(:,3) = obj.modNormal(:,1).* Vinf(:,2)-obj.modNormal(:,2).* Vinf(:,1);
                s(:,1) = Tvec(:,2).*obj.modNormal(:,3)-Tvec(:,3).*obj.modNormal(:,2);
                s(:,2) = Tvec(:,3).*obj.modNormal(:,1)-Tvec(:,1).*obj.modNormal(:,3);
                s(:,3) = Tvec(:,1).*obj.modNormal(:,2)-Tvec(:,2).*obj.modNormal(:,1);
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                                                    
                    end
                    
                    if obj.settingUNLSI.propCalcFlag == 1
                        RHV = obj.RHS*sigmas;
                        RHVprop = obj.calcPropRHV(obj,T);
                        RHV = RHV + RHVprop;
                    else
                        RHV = obj.RHS*sigmas;
                    end
                    u =  -obj.LHS\RHV;
                    %figure(2);clf;
                    %plot(u);
                    dv = zeros(nPanel,3);
                    for i = 1:3
                        dv(:,i) = obj.mu2v{i}*u;
                    end
                    dv = Vinf + dv;
                    dv = dv - obj.orgNormal.*(dot(obj.orgNormal,dv,2));
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 1,iterflow) = (1-(sqrt(sum(dv(obj.paneltype==1 & obj.cpcalctype == 1,:).^2,2))).^2)./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 3,iterflow) = (-2*(dv(obj.paneltype==1 & obj.cpcalctype == 3,1))-dv(obj.paneltype==1 & obj.cpcalctype == 3,2).^2-dv(obj.paneltype==1 & obj.cpcalctype == 3,3).^2)./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 2,iterflow) = (-2*(dv(obj.paneltype==1 & obj.cpcalctype == 2,1)-Vinf(obj.paneltype==1 & obj.cpcalctype == 2,1)))./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-0.139-0.419.*(obj.flow{flowNo}.Mach-0.161).^2);
                    uinterp = [];
                    for i = 1:numel(obj.wakeline)
                        uinterp = [uinterp,interp1(obj.LLT.sp{i},obj.LLT.calcMu{i}*u,obj.LLT.sinterp{i},'linear','extrap')];
                    end
                    if numel(obj.wakeline)>0
                        Vind = obj.LLT.Qij*uinterp';
                        CLt = (2.*horzcat(obj.LLT.spanel{:}).*cos(horzcat(obj.LLT.phiinterp{:}))*uinterp')/(0.5*obj.SREF)/(1-obj.flow{flowNo}.Mach^2);
                        CDt = ((uinterp.*horzcat(obj.LLT.spanel{:}))*Vind)/(0.5*obj.SREF)/norm(1-obj.flow{flowNo}.Mach^2)^3;
                    else
                        CLt = 0;
                        CDt = 0;
                    end
                    else
                    %超音速
                    delta = zeros(nbPanel,1);
                    iter = 1;
                    %各パネルが主流となす角度を求める
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            delta(iter,1) = acos(dot(obj.orgNormal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
                            iter = iter+1;
                        end
                    end
                    %用意された応答曲面をもちいてパネルの角度からCpを求める
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = obj.flow{flowNo}.pp(delta);%Cp
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-obj.flow{flowNo}.Mach.^(-2)+0.57.*obj.flow{flowNo}.Mach.^(-4));
                end
                %Cp⇒力への変換
                dCA_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,3)).*obj.area./obj.SREF;
                dCA_f = (+obj.Cfe{flowNo}.*s(:,1)).*obj.area./obj.SREF;
                dCY_f = (+obj.Cfe{flowNo}.*s(:,2)).*obj.area./obj.SREF;
                dCN_f = (+obj.Cfe{flowNo}.*s(:,3)).*obj.area./obj.SREF;
                dCM = cross(obj.center-repmat(obj.XYZREF,[size(obj.center,1),1]),[dCA_p+dCA_f,dCY_p+dCY_f,dCN_p+dCN_f]);
                dCMX = dCM(:,1)./obj.BREF;
                dCMY = dCM(:,2)./obj.CREF;
                dCMZ = dCM(:,3)./obj.BREF;

                if obj.halfmesh == 1
                    %半裁
                    CNp = sum(dCN_p)*2;
                    CAp = sum(dCA_p)*2;
                    CYp = 0;
                    CNf = sum(dCN_f)*2;
                    CAf = sum(dCA_f)*2;
                    CYf = 0;
                    CMX = 0;
                    CMY = sum(dCMY)*2;
                    CMZ = 0;
                else
                    CNp = sum(dCN_p);
                    CAp = sum(dCA_p);
                    CYp = sum(dCY_p);
                    CNf = sum(dCN_f);
                    CAf = sum(dCA_f);
                    CYf = sum(dCY_f);
                    CMX = sum(dCMX);
                    CMY = sum(dCMY);
                    CMZ = sum(dCMZ);
                end
                CL = T(:,3)'*[CAp+CAf;CYp+CYf;CNp+CNf];
                CDi = T(:,1)'*[CAp;CYp;CNp];
                if obj.flow{flowNo}.Mach < 1
                    if obj.halfmesh == 1
                        CLt = 2*CLt;
                        CDt = 2*CDt;
                    end
                else
                    CLt = CL;
                    CDt = CDi;
                end
                CDo = T(:,1)'*[CAf;CYf;CNf];
                CDtot = CDi+CDo;
                CDtott = CDt+CDo;
                if obj.halfmesh == 1
                    CY = 0;
                else
                    CY = T(:,2)'*[CAp+CAf;CYp+CYf;CNp+CNf];
                end
                AR = obj.BREF^2/obj.SREF;
                if isfield(obj.flow{flowNo},"Re")
                    ReOut = obj.flow{flowNo}.Re/1000000;
                else
                    ReOut = 0;
                end
                CMl = -CMX;
                CMm = CMY;
                CMn = -CMZ;
                obj.AERODATA{flowNo}(iterflow,:) = [beta(iterflow),obj.flow{flowNo}.Mach,alpha(iterflow),ReOut,CL,CLt,CDo,CDi,CDtot,CDt,CDtott,CY,CLt/CDtott,CLt^2/pi/AR/CDt,CAp+CAf,CYp+CYf,CNp+CNf,CMX,CMY,CMZ,CMl,CMm,CMn,0];
                
            end
        end

        function [u,R] = solvePertPotential(obj,flowNo,alpha,beta)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %Adjoint法実装のためポテンシャルの変動値のみ求める。

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(size(alpha) ~= size(beta))
                error("analysis points are not match");
            end
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            R = [];
            u = [];
            for iterflow = 1:numel(alpha)
                T(1,1) = cosd(alpha(iterflow))*cosd(beta(iterflow));
                T(1,2) = cosd(alpha(iterflow))*sind(beta(iterflow));
                T(1,3) = -sind(alpha(iterflow));
                T(2,1) = -sind(beta(iterflow));
                T(2,2) = cosd(beta(iterflow));
                T(2,3) = 0;
                T(3,1) = sind(alpha(iterflow))*cosd(beta(iterflow));
                T(3,2) = sind(alpha(iterflow))*sind(beta(iterflow));
                T(3,3) = cosd(alpha(iterflow));
                if isempty(obj.settingUNLSI.angularVelocity)
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0]))';
                    end
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.settingUNLSI.angularVelocity(iterflow,:)./180.*pi,rvec(:)')'))';
                    end
                end
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                    end
                    if obj.settingUNLSI.propCalcFlag == 1
                        RHV = obj.RHS*sigmas;
                        RHVprop = obj.calcPropRHV(obj,T);
                        RHV = RHV + RHVprop;
                    else
                       RHV = obj.RHS*sigmas;
                    end
                    usolve =  -obj.LHS\RHV;
                    Rsolve = obj.LHS*usolve+RHV;
                    u = [u;usolve];
                    R = [R;Rsolve];
                else
                    %超音速
                    delta = zeros(nbPanel,1);
                    iter = 1;
                    %各パネルが主流となす角度を求める
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            delta(iter,1) = acos(dot(obj.orgNormal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
                            iter = iter+1;
                        end
                    end
                    %用意された応答曲面をもちいてパネルの角度からCpを求める
                    usolve = obj.flow{flowNo}.pp(delta);%Cp
                    Rsolve = usolve-obj.flow{flowNo}.pp(delta);%Cp
                    u = [u;usolve];
                    R = [R;Rsolve];
                end
            end
        end
        
        function [AERODATA,Cp,Cfe,R,obj] = solveFlowForAdjoint(obj,u,flowNo,alpha,beta)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %ポテンシャルから力を求める
            %結果は配列に出力される
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CDo 7:CDi 8:CDtot 9:CDt 10:CDtot_t 11:CS 12:L/D E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
            %上記で求めていないものは0が代入される
            %u: doubletの強さ
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %omega:主流の回転角速度(deg/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(size(alpha) ~= size(beta))
                error("analysis points are not match");
            end
            obj.AERODATA{flowNo} = [];
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            obj.Cp{flowNo} = zeros(nPanel,numel(alpha));
            for iterflow = 1:numel(alpha)
                T(1,1) = cosd(alpha(iterflow))*cosd(beta(iterflow));
                T(1,2) = cosd(alpha(iterflow))*sind(beta(iterflow));
                T(1,3) = -sind(alpha(iterflow));
                T(2,1) = -sind(beta(iterflow));
                T(2,2) = cosd(beta(iterflow));
                T(2,3) = 0;
                T(3,1) = sind(alpha(iterflow))*cosd(beta(iterflow));
                T(3,2) = sind(alpha(iterflow))*sind(beta(iterflow));
                T(3,3) = cosd(alpha(iterflow));
                if isempty(obj.settingUNLSI.angularVelocity)
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0]))';
                    end
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.settingUNLSI.angularVelocity(iterflow,:)./180.*pi,rvec(:)')'))';
                    end
                end
                Tvec(:,1) = obj.modNormal(:,2).* Vinf(:,3)-obj.modNormal(:,3).* Vinf(:,2);
                Tvec(:,2) = obj.modNormal(:,3).* Vinf(:,1)-obj.modNormal(:,1).* Vinf(:,3);
                Tvec(:,3) = obj.modNormal(:,1).* Vinf(:,2)-obj.modNormal(:,2).* Vinf(:,1);
                s(:,1) = Tvec(:,2).*obj.modNormal(:,3)-Tvec(:,3).*obj.modNormal(:,2);
                s(:,2) = Tvec(:,3).*obj.modNormal(:,1)-Tvec(:,1).*obj.modNormal(:,3);
                s(:,3) = Tvec(:,1).*obj.modNormal(:,2)-Tvec(:,2).*obj.modNormal(:,1);
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                    end
                    if obj.settingUNLSI.propCalcFlag == 1
                       RHVprop = obj.calcPropRHV(obj,T);
                    end
                    usolve = u(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1);
                    dv = zeros(nPanel,3);
                    for i = 1:3
                        dv(:,i) = obj.mu2v{i}*usolve;
                    end
                    dv = Vinf + dv;
                    dv = dv - obj.orgNormal.*(dot(obj.orgNormal,dv,2));
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 1,iterflow) = (1-(sqrt(sum(dv(obj.paneltype==1 & obj.cpcalctype == 1,:).^2,2))).^2)./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 3,iterflow) = (-2*(dv(obj.paneltype==1 & obj.cpcalctype == 3,1)+-1)-dv(obj.paneltype==1 & obj.cpcalctype == 3,3).^2)./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==1 & obj.cpcalctype == 2,iterflow) = (-2*(dv(obj.paneltype==1 & obj.cpcalctype == 2,1)-Vinf(obj.paneltype==1 & obj.cpcalctype == 2,1)))./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-0.139-0.419.*(obj.flow{flowNo}.Mach-0.161).^2);
                    uinterp = [];
                    for i = 1:numel(obj.wakeline)
                        uinterp = [uinterp,interp1(obj.LLT.sp{i},obj.LLT.calcMu{i}*(usolve),obj.LLT.sinterp{i},'linear','extrap')];
                    end
                    Vind = obj.LLT.Qij*uinterp';

                    CLt = (2.*horzcat(obj.LLT.spanel{:}).*cos(horzcat(obj.LLT.phiinterp{:}))*uinterp')/(0.5*obj.SREF)/(1-obj.flow{flowNo}.Mach^2);
                    CDt = ((uinterp.*horzcat(obj.LLT.spanel{:}))*Vind)/(0.5*obj.SREF)/norm(1-obj.flow{flowNo}.Mach^2)^3;
                else
                    %超音速
                    delta = zeros(nbPanel,1);
                    iter = 1;
                    %各パネルが主流となす角度を求める
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            delta(iter,1) = acos(dot(obj.orgNormal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
                            iter = iter+1;
                        end
                    end
                    usolve = u(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1);
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = usolve;%Cp
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-obj.flow{flowNo}.Mach.^(-2)+0.57.*obj.flow{flowNo}.Mach.^(-4));
                end
                %Cp⇒力への変換
                dCA_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-obj.Cp{flowNo}(:,iterflow).*obj.modNormal(:,3)).*obj.area./obj.SREF;
                dCA_f = (+obj.Cfe{flowNo}.*s(:,1)).*obj.area./obj.SREF;
                dCY_f = (+obj.Cfe{flowNo}.*s(:,2)).*obj.area./obj.SREF;
                dCN_f = (+obj.Cfe{flowNo}.*s(:,3)).*obj.area./obj.SREF;
                dCM = cross(obj.center-repmat(obj.XYZREF,[size(obj.center,1),1]),[dCA_p+dCA_f,dCY_p+dCY_f,dCN_p+dCN_f]);
                dCMX = dCM(:,1)./obj.BREF;
                dCMY = dCM(:,2)./obj.CREF;
                dCMZ = dCM(:,3)./obj.BREF;
                if obj.halfmesh == 1
                    %半裁
                    CNp = sum(dCN_p)*2;
                    CAp = sum(dCA_p)*2;
                    CYp = 0;
                    CNf = sum(dCN_f)*2;
                    CAf = sum(dCA_f)*2;
                    CYf = 0;
                    CMX = 0;
                    CMY = sum(dCMY)*2;
                    CMZ = 0;
                else
                    CNp = sum(dCN_p);
                    CAp = sum(dCA_p);
                    CYp = sum(dCY_p);
                    CNf = sum(dCN_f);
                    CAf = sum(dCA_f);
                    CYf = sum(dCY_f);
                    CMX = sum(dCMX);
                    CMY = sum(dCMY);
                    CMZ = sum(dCMZ);
                end
    
                CL = T(:,3)'*[CAp+CAf;CYp+CYf;CNp+CNf];
                CDi = T(:,1)'*[CAp;CYp;CNp];
                if obj.flow{flowNo}.Mach < 1
                    if obj.halfmesh == 1
                        CLt = 2*CLt;
                        CDt = 2*CDt;
                    end
                else
                    CLt = CL;
                    CDt = CDi;
                end
                CDo = T(:,1)'*[CAf;CYf;CNf];
                CDtot = CDi+CDo;
                CDtott = CDt+CDo;
                if obj.halfmesh == 1
                    CY = 0;
                else
                    CY = T(:,2)'*[CAp+CAf;CYp+CYf;CNp+CNf];
                end
                AR = obj.BREF^2/obj.SREF;                
                if isfield(obj.flow{flowNo},"Re")
                    ReOut = obj.flow{flowNo}.Re/1000000;
                else
                    ReOut = 0;
                end
                CMl = -CMX;
                CMm = CMY;
                CMn = -CMZ;
                obj.AERODATA{flowNo}(iterflow,:) = [beta(iterflow),obj.flow{flowNo}.Mach,alpha(iterflow),ReOut,CL,CLt,CDo,CDi,CDtot,CDt,CDtott,CY,CLt/CDtott,CLt^2/pi/AR/CDt,CAp+CAf,CYp+CYf,CNp+CNf,CMX,CMY,CMZ,CMl,CMm,CMn,0];
                AERODATA = obj.AERODATA;
                Cp = obj.Cp;
                Cfe = obj.Cfe;
                if obj.flow{flowNo}.Mach < 1
                    RHV = obj.RHS*sigmas;
                    if obj.settingUNLSI.propCalcFlag == 1
                        RHV = RHV+RHVprop;
                    end
                    R(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1) = obj.LHS*usolve+RHV;
                else
                    R(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1) = usolve-obj.flow{flowNo}.pp(delta);
                end
            end
            %disp([CL,CDo,CDi,CDtot,CMY]);
        end

        function AERODATA = getAERODATA(obj,alpha,beta,Mach,Re)
            if nargin == 1
                AERODATA = obj.AERODATA;
            else
                %結果を検索して引っかかったrowのものを出力
                alldata = vertcat(obj.AERODATA{:});
                if obj.settingUNLSI.resultSearchMethod == "or"
                    if nargin == 2
                        outindex = [find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 3
                        outindex = [find(any(alldata(:,1) == beta,2));find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 4
                        outindex =[find(any(alldata(:,1) == beta,2));find(any(alldata(:,2) == Mach,2));find(any(alldata(:,3) == alpha,2))];
                    else
                       outindex =[find(any(alldata(:,1) == beta,2));find(any(alldata(:,2) == Mach,2));find(any(alldata(:,3) == alpha,2));find(any(alldata(:,4) == Re/1000000,2))];
                    end
                    outindex = unique(outindex);
                    AERODATA = alldata(outindex,:);
                elseif obj.settingUNLSI.resultSearchMethod == "and"
                    if nargin == 2
                        outindex = [find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 3
                        outindex = intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,3) == alpha,2)));
                    elseif nargin == 4
                        outindex =intersect(intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,2) == Mach,2))),find(any(alldata(:,3) == alpha,2)));
                    else
                       outindex =intersect(intersect(intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,2) == Mach,2))),find(any(alldata(:,3) == alpha,2))),find(any(alldata(:,4) == Re/1000000,2)));
                    end
                    AERODATA = alldata(outindex,:);
                else
                    AERODATA = [];
                end
            end
        end

        function Cp = getCp(obj,alpha,beta,Mach,Re)
            if nargin == 1
                Cp = obj.Cp;
            else
                %結果を検索して引っかかったrowのものを出力
                alldata = vertcat(obj.AERODATA{:});
                allCp = horzcat(obj.Cp{:});
                if obj.settingUNLSI.resultSearchMethod == "or"
                    if nargin == 2
                        outindex = [find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 3
                        outindex = [find(any(alldata(:,1) == beta,2));find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 4
                        outindex =[find(any(alldata(:,1) == beta,2));find(any(alldata(:,2) == Mach,2));find(any(alldata(:,3) == alpha,2))];
                    else
                       outindex =[find(any(alldata(:,1) == beta,2));find(any(alldata(:,2) == Mach,2));find(any(alldata(:,3) == alpha,2));find(any(alldata(:,4) == Re/1000000,2))];
                    end
                    outindex = unique(outindex);
                    Cp = allCp(:,outindex);s
                elseif obj.settingUNLSI.resultSearchMethod == "and"
                    if nargin == 2
                        outindex = [find(any(alldata(:,3) == alpha,2))];
                    elseif nargin == 3
                        outindex = intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,3) == alpha,2)));
                    elseif nargin == 4
                        outindex =intersect(intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,2) == Mach,2))),find(any(alldata(:,3) == alpha,2)));
                    else
                       outindex =intersect(intersect(intersect(find(any(alldata(:,1) == beta,2)),find(any(alldata(:,2) == Mach,2))),find(any(alldata(:,3) == alpha,2))),find(any(alldata(:,4) == Re/1000000,2)));
                    end
                    Cp = allCp(:,outindex);
                else
                    Cp = [];
                end
            end
        end

        function [DYNCOEF,dynCoefStruct] = getDYNCOEF(obj,alpha,beta,Mach,Re)
            if obj.settingUNLSI.resultSearchMethod == "and"
                outindex = [];
                for i = 1:size(obj.DYNCOEF,1)
                    for j = 1:size(obj.DYNCOEF,2)
                        if not(isempty(obj.DYNCOEF{i,j}))
                            if any(obj.DYNCOEF{i,j}(3,1) == alpha) && any(obj.DYNCOEF{i,j}(1,1) == beta) && any(obj.DYNCOEF{i,j}(2,1) == Mach) && any(obj.DYNCOEF{i,j}(4,1) == Re/1000000)
                                outindex = [outindex;[i,j]];
                            end
                        end
                    end
                end
            elseif obj.settingUNLSI.resultSearchMethod == "or"
                outindex = [];
                for i = 1:size(obj.DYNCOEF,1)
                    for j = 1:size(obj.DYNCOEF,2)
                        if not(isempty(obj.DYNCOEF{i,j}))
                            if any(obj.DYNCOEF{i,j}(3,1) == alpha) || any(obj.DYNCOEF{i,j}(1,1) == beta) || any(obj.DYNCOEF{i,j}(2,1) == Mach) || any(obj.DYNCOEF{i,j}(4,1) == Re/1000000)
                                outindex = [outindex;[i,j]];
                            end
                        end
                    end
                end
            end
            iter = 1;
            for i = 1:size(outindex,1)
               DYNCOEF{iter} = obj.DYNCOEF{outindex(iter,1),outindex(iter,2)};
               if(nargout>1)
                    dynCoefStruct{iter}.beta = DYNCOEF{iter}(1,1);
                    dynCoefStruct{iter}.Mach = DYNCOEF{iter}(2,1);
                    dynCoefStruct{iter}.alpha = DYNCOEF{iter}(3,1);
                    dynCoefStruct{iter}.Re = DYNCOEF{iter}(4,1)*1000000;
                    dynCoefStruct{iter}.Cyb = DYNCOEF{iter}(2,2);
                    dynCoefStruct{iter}.Clb = DYNCOEF{iter}(4,2);
                    dynCoefStruct{iter}.Cnb = DYNCOEF{iter}(6,2);
                    dynCoefStruct{iter}.Cxa = DYNCOEF{iter}(1,3);
                    dynCoefStruct{iter}.Cza = DYNCOEF{iter}(3,3);
                    dynCoefStruct{iter}.Cma = DYNCOEF{iter}(5,3);
                    dynCoefStruct{iter}.Cyp = DYNCOEF{iter}(2,4);
                    dynCoefStruct{iter}.Clp = DYNCOEF{iter}(4,4);
                    dynCoefStruct{iter}.Cnp = DYNCOEF{iter}(6,4);
                    dynCoefStruct{iter}.Cxq = DYNCOEF{iter}(1,5);
                    dynCoefStruct{iter}.Czq = DYNCOEF{iter}(3,5);
                    dynCoefStruct{iter}.Cmq = DYNCOEF{iter}(5,5);
                    dynCoefStruct{iter}.Cyr = DYNCOEF{iter}(2,6);
                    dynCoefStruct{iter}.Clr = DYNCOEF{iter}(4,6);
                    dynCoefStruct{iter}.Cnr = DYNCOEF{iter}(6,6);
                end
               iter = iter+1;
            end 
        end

        function  [modal,DYNCOEF,dynCoefStruct] = getModal(obj,alpha,beta,Mach,Re,UREF,mass,Inatia)
            [DYNCOEF,dynCoefStruct] = obj.getDYNCOEF(alpha,beta,Mach,Re);
            for i = 1:numel(DYNCOEF)
	            %Mq = obj.settingUNLSI.rho*UREF*SREF*CREF*CREF/(4*Inatia(2,2))*dynCoefStruct{i}.Cmq;
	            Yb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF/(2*mass)*dynCoefStruct{i}.Cyb;
	            Lb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF*obj.BREF/(2*Inatia(1,1))*dynCoefStruct{i}.Clb;
	            Nb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF*obj.BREF/(2*Inatia(3,3))*dynCoefStruct{i}.Cnb;
	            %Yp = obj.settingUNLSI.rho*UREF*SREF*BREF/(4*mass)*dynCoefStruct{i}.Cyp;
	            Lp = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(1,1))*dynCoefStruct{i}.Clp;
	            Np = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(3,3))*dynCoefStruct{i}.Cnp;
	            %Yr = obj.settingUNLSI.rho*UREF*SREF*BREF/(4*mass)*dynCoefStruct{i}.Cyr;
	            Lr = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(1,1))*dynCoefStruct{i}.Clr;
	            Nr = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(3,3))*dynCoefStruct{i}.Cnr;
                %Xu = obj.settingUNLSI.rho*UREF*SREF/(2*mass)*(dynCoefStruct{i}.Cxu);
                %Zu = obj.settingUNLSI.rho*UREF*SREF/(2*mass)*(dynCoefStruct{i}.CZu-2*dynCoefStruct{i}.CL);
                %Xa = obj.settingUNLSI.rho*UREF^2*SREF/(2*mass)*dynCoefStruct{i}.Cxa;
                Za = obj.settingUNLSI.rho*UREF^2*obj.SREF/(2*mass)*dynCoefStruct{i}.Cza;
                %Zq = obj.settingUNLSI.rho*UREF*SREF*CREF/(4*mass)*(dynCoefStruct{i}.Czq);
                %Mu = obj.settingUNLSI.rho*UREF*SREF*CREF/(2*Inatia(2,2))*dynCoefStruct{i}.Cmu;
                Ma = obj.settingUNLSI.rho*UREF^2*obj.SREF*obj.CREF/(2*Inatia(2,2))*dynCoefStruct{i}.Cma;
                Mq = obj.settingUNLSI.rho*UREF*obj.SREF*obj.CREF^2/(4*Inatia(2,2))*dynCoefStruct{i}.Cmq;
        
                
                modal.shortPriod.omegan = sqrt(-Ma+Za/UREF*Mq);
                modal.shortPriod.eta =(-Za/UREF-Mq)/2/modal.shortPriod.omegan;
                
                if not(Yb == 0 && Lp == 0)
                    modal.roll.T = -1/(Lp);
                    lambdas_u = (Lb*Nr-Nb*Lr);
                    lambdas_l = (Lp*Nr-Np*Lr)*Yb+(Lb*Nr-Nb*Lr)*UREF-Lb*9.8;
                    modal.spiral.T = lambdas_l/lambdas_u;
            
                    modal.datchRoll.omegan = sqrt(Nb-(Np/Lp)*Lb);
                    modal.datchRoll.eta = -(Nr-(Np/Lp)*Lr+(Np/Lp^2)*Lb)/2/modal.datchRoll.omegan;
                else
                    modal.roll = [];
                    modal.spiral = [];
                    modal.datchRoll = [];
                end
            end
        end

        function obj = calcDynCoef(obj,alpha,beta,Mach,Re)
                %%%動微係数の計算%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%TO DO alpha とbetaも
                %結果はDyncoefに格納される
                % 横軸    beta, alpha, p, q, r
                %        --------------------
                %    Cx  |
                % 縦 Cz  |
                % 軸 Cmx |
                %    Cmy |
                %    Cmz |
                %dynCoefStructは構造体で出力する。
                %flowNo:解きたい流れのID
                %alpha:迎角[deg]
                %beta:横滑り角[deg]
                %difference:有限差分の方法 "forward"(デフォルト)-前進差分 "central"-中心差分
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(obj.flowNoTable)
                    if nargin == 5
                        obj = obj.flowCondition(1,Mach,Re);
                    else
                        obj = obj.flowCondition(1,Mach);
                    end
                    flowNo = 1;
                else
                    flowNo = 0;
                    for i = 1:size(obj.flowNoTable,1)
                        if obj.flowNoTable(i,1) == Mach && obj.flowNoTable(i,2) == Re
                            flowNo = i;
                        end
                    end
                    if flowNo == 0
                        if nargin == 5
                            obj = obj.flowCondition(size(obj.flowNoTable,1)+1,Mach,Re);
                        else
                            obj = obj.flowCondition(size(obj.flowNoTable,1)+1,Mach);
                        end
                        flowNo = size(obj.flowNoTable,1);
                    end
                end
                u0 = obj.solvePertPotential(flowNo,alpha,beta);
                [udwf,udwr] = obj.calcDynCoefdu(flowNo,alpha,beta);
                [~,obj] =  obj.calcDynCoefforAdjoint(u0,udwf,udwr,flowNo,alpha,beta);
        end

        function [udwf,udwr] = calcDynCoefdu(obj,flowNo,alpha,beta)
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            dw = eps^(1/3);%rad
            if obj.flow{flowNo}.Mach < 1
                udw = zeros(nbPanel*numel(alpha),5);
                for iterflow = 1:numel(alpha)
                    %betaについて
                    Vinfdw = repmat([-cosd(alpha(iterflow))*sind(beta(iterflow))*dw,-cosd(beta(iterflow))*dw,-sind(alpha(iterflow))*sind(beta(iterflow))*dw],[nPanel,1]);
                    vdw = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            vdw(iter,1) = dot(Vinfdw(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                    end
                    rhsdw = obj.RHS*vdw;
                    udw(1+nbPanel*(iterflow-1):nbPanel*iterflow,1) = -obj.LHS\rhsdw;
        
                    %alphaについて
                    Vinfdw = repmat([-sind(alpha(iterflow))*cosd(beta(iterflow))*dw,0,cosd(alpha(iterflow))*cosd(beta(iterflow))*dw],[nPanel,1]);
                    vdw = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            vdw(iter,1) = dot(Vinfdw(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                    end
                    rhsdw = obj.RHS*vdw;
                    udw(1+nbPanel*(iterflow-1):nbPanel*iterflow,2) = -obj.LHS\rhsdw;
                end

                for jter = 1:3%p,q,rについて差分をとる
                    if isempty(obj.settingUNLSI.angularVelocity)
                        omegadw = zeros(1,3);
                    else
                        omegadw = obj.settingUNLSI.angularVelocity(iterflow,:);
                    end
                    omegadw(jter) = omegadw(jter)+dw;
                    Vinfdw = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinfdw(i,:) = -(cross(omegadw,rvec(:)))';
                    end
                    vdw = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            vdw(iter,1) = dot(Vinfdw(i,:)',obj.orgNormal(i,:)');
                            iter = iter+1;
                        end
                    end
                    rhsdw = obj.RHS*vdw;
                    buff = -obj.LHS\rhsdw;
                    for iterflow = 1:numel(alpha)
                        udw(1+nbPanel*(iterflow-1):nbPanel*iterflow,2+jter) = buff;
                    end
                end

                udwf = udw;
                udwr = udw;
                if obj.settingUNLSI.deflDerivFlag == 1
                    orgDeflAngle = obj.deflAngle;
                    for iterflow = 1:numel(alpha)
                        u0 = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                        for i = 1:numel(obj.deflGroup)
                            [deflID,~] = find(obj.deflGroup{i}.ID(:)'==orgDeflAngle(:,1));
                            if any(abs(vecnorm(orgDeflAngle(deflID,2:4),2,2)-1)>sqrt(eps))
                                error("Defl Axis Not Set");
                            end
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)+dw*180/pi*obj.deflGroup{i}.gain(:));
                            uf = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                            udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,5+i) = uf-u0;
                            [deflID,~] = find(obj.deflGroup{i}.ID==orgDeflAngle(:,1));
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)-dw*180/pi*obj.deflGroup{i}.gain(:));
                            ur = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                            udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,5+i) = u0-ur;
                            [deflID,~] = find(obj.deflGroup{i}.ID==orgDeflAngle(:,1));
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5));
                        end
                    end
                end
            else
                for iterflow = 1:numel(alpha)
                    u0 = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                    %betaについて
                    uf = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow)+dw/pi*180);
                    ur = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow)-dw/pi*180);
                    udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,1) = uf-u0;
                    udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,1) = u0-ur;
                    %alphaについて
                    uf = obj.solvePertPotential(flowNo,alpha(iterflow)+dw/pi*180,beta(iterflow));
                    ur = obj.solvePertPotential(flowNo,alpha(iterflow)-dw/pi*180,beta(iterflow));
                    udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,2) = uf-u0;
                    udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,2) = u0-ur;
                    avorg = obj.settingUNLSI.angularVelocity;
                    for jter = 1:3%p,q,rについて差分をとる
                        if isempty(obj.settingUNLSI.angularVelocity)
                            omegadwf = zeros(1,3);
                            omegadwr = zeros(1,3);
                        else
                            omegadwf = obj.settingUNLSI.angularVelocity(iterflow,:);
                            omegadwr = obj.settingUNLSI.angularVelocity(iterflow,:);
                        end
                        omegadwf(jter) = omegadwf(jter)+dw/pi*180;
                        omegadwr(jter) = omegadwr(jter)-dw/pi*180;
                        obj.settingUNLSI.angularVelocity = omegadwf;
                        uf = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow),omegadwf);
                        obj.settingUNLSI.angularVelocity = omegadwr;
                        ur = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow),omegadwr);
                        obj.settingUNLSI.angularVelocity = avorg;
                        udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,2+jter) = uf-u0;
                        udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,2+jter) = u0-ur;
                    end
                end
                if obj.settingUNLSI.deflDerivFlag == 1
                    orgDeflAngle = obj.deflAngle;
                    for iterflow = 1:numel(alpha)
                        u0 = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow),0);
                        for i = 1:numel(obj.deflGroup)
                            [deflID,~] = find(obj.deflGroup{i}.ID(:)'==orgDeflAngle(:,1));
                            if any(abs(vecnorm(orgDeflAngle(deflID,2:4),2,2)-1)>sqrt(eps))
                                error("Defl Axis Not Set");
                            end
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)+dw*180/pi*obj.deflGroup{i}.gain(:));
                            uf = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow),0);
                            udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,5+i) = uf-u0;
                            [deflID,~] = find(obj.deflGroup{i}.ID==orgDeflAngle(:,1));
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)-dw*180/pi*obj.deflGroup{i}.gain(:));
                            ur = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow),0);
                            udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,5+i) = u0-ur;
                            [deflID,~] = find(obj.deflGroup{i}.ID==orgDeflAngle(:,1));
                            obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5));
                        end
                    end
                end
            end
        end

        function [DYNCOEF,obj] = calcDynCoefforAdjoint(obj,u0,udwf,udwr,flowNo,alpha,beta)
                %%%動微係数の計算%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%TO DO alpha とbetaも
                %結果はDyncoefに格納される
                % 横軸     flowdata  beta, alpha, p, q, r
                %        --------------------
                %    Cx  |beta
                % 縦 Cy  |Mach
                %    Cz  |alpha
                % 軸 Cmx |Re/e6
                %    Cmy | 0
                %    Cmz | 0
                %dynCoefStructは構造体で出力する。
                %flowNo:解きたい流れのID
                %alpha:迎角[deg]
                %beta:横滑り角[deg]
                %difference:有限差分の方法 "forward"(デフォルト)-前進差分 "central"-中心差分
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %ノミナルについて流れ変数を取得する。
                dw = eps^(1/3);%rad
                %有限差分による計算
                ind = 15:20;%AERODATA = [18:CMX, 19:CMY, 20: CMZ]
                nondim = [obj.BREF/(2*1),obj.CREF/(2*1),obj.BREF/(2*1)];
                for i = 1:numel(alpha)
                    obj.DYNCOEF{flowNo,i} = zeros(6,6);
                end
                %betaについて
                AERODATAf = obj.solveFlowForAdjoint(u0+udwf(:,1),flowNo,alpha,beta+dw*180/pi);
                AERODATAr = obj.solveFlowForAdjoint(u0-udwr(:,1),flowNo,alpha,beta-dw*180/pi);
                tmp = (AERODATAf{flowNo} - AERODATAr{flowNo})./(dw*2);
                for i = 1:size(AERODATAf{flowNo},1)
                    obj.DYNCOEF{flowNo,i}(:,2) = tmp(i,ind)'; % beta
                end

                %alphaについて
                AERODATAf = obj.solveFlowForAdjoint(u0+udwf(:,2),flowNo,alpha+dw*180/pi,beta);
                AERODATAr = obj.solveFlowForAdjoint(u0-udwr(:,2),flowNo,alpha-dw*180/pi,beta);
                tmp = (AERODATAf{flowNo} - AERODATAr{flowNo})./(dw*2);
                for i = 1:size(AERODATAf{flowNo},1)
                    obj.DYNCOEF{flowNo,i}(:,3) = tmp(i,ind)'; % alpha
                end
                avorg = obj.settingUNLSI.angularVelocity;
                for jter = 1:3%p,q,rについて差分をとる
                    if isempty(obj.settingUNLSI.angularVelocity)
                        omegadwf = zeros(1,3);
                        omegadwr = zeros(1,3);
                    else
                        omegadwf = obj.settingUNLSI.angularVelocity(iterflow,:);
                        omegadwr = obj.settingUNLSI.angularVelocity(iterflow,:);
                    end
                    omegadwf(jter) = omegadwf(jter)+dw.*180./pi;
                    omegadwr(jter) = omegadwr(jter)-dw.*180./pi;
                    obj.settingUNLSI.angularVelocity = repmat(omegadwf,[numel(alpha),1]);
                    AERODATAf = obj.solveFlowForAdjoint(u0+udwf(:,2+jter),flowNo,alpha,beta);
                    obj.settingUNLSI.angularVelocity = repmat(omegadwr,[numel(alpha),1]);
                    AERODATAr = obj.solveFlowForAdjoint(u0-udwr(:,2+jter),flowNo,alpha,beta);
                    obj.settingUNLSI.angularVelocity = avorg;
                    tmp = (AERODATAf{flowNo} - AERODATAr{flowNo})./(2*dw*nondim(jter));
                    for i = 1:size(AERODATAf{flowNo},1)
                        obj.DYNCOEF{flowNo,i}(:,3+jter) = tmp(i,ind)'; % pqr
                    end
                end
                if size(udwf,2)>5
                    orgDeflAngle = obj.deflAngle;
                    for jter = 1:numel(obj.deflGroup)
                        [deflID,~] = find(obj.deflGroup{jter}.ID(:)'==orgDeflAngle(:,1));
                        if any(abs(vecnorm(orgDeflAngle(deflID,2:4),2,2)-1)>sqrt(eps))
                            error("Defl Axis Not Set");
                        end
                        obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)+dw*180/pi*obj.deflGroup{jter}.gain(:));
                        AERODATAf = obj.solveFlowForAdjoint(u0+udwf(:,5+jter),flowNo,alpha,beta);
                        obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5)-dw*180/pi*obj.deflGroup{jter}.gain(:));
                        AERODATAr = obj.solveFlowForAdjoint(u0-udwr(:,5+jter),flowNo,alpha,beta);
                        obj = obj.setDeflAngle(orgDeflAngle(deflID,1),orgDeflAngle(deflID,2:4),orgDeflAngle(deflID,5));
                        tmp = (AERODATAf{flowNo} - AERODATAr{flowNo})./(2*dw);
                        for i = 1:size(AERODATAf{flowNo},1)
                            obj.DYNCOEF{flowNo,i}(:,6+jter) = tmp(i,ind)'; % defl
                        end
                    end
                end

                axisRot = [ 0, 1, 0,-1, 0,-1;
                           -1, 0,-1, 0, 1, 0;
                            0,-1, 0, 1, 0, 1;
                           -1, 0,-1, 0, 1, 0;
                            0,-1, 0, 1, 0, 1]';

                for i = 1:numel(obj.deflGroup)
                    axisRot(:,5+i) = [-1, 1, -1,-1, 1,-1]';
                end
                for iterflow = 1:numel(alpha)
                    obj.DYNCOEF{flowNo,iterflow}(:,2:end) = (axisRot.*obj.DYNCOEF{flowNo,iterflow}(:,2:end));
                    %安定軸に変換
                    T(1,1) = cosd(alpha(iterflow))*cosd(beta(iterflow));
                    T(1,2) = cosd(alpha(iterflow))*sind(beta(iterflow));
                    T(1,3) = -sind(alpha(iterflow));
                    T(2,1) = -sind(beta(iterflow));
                    T(2,2) = cosd(beta(iterflow));
                    T(2,3) = 0;
                    T(3,1) = sind(alpha(iterflow))*cosd(beta(iterflow));
                    T(3,2) = sind(alpha(iterflow))*sind(beta(iterflow));
                    T(3,3) = cosd(alpha(iterflow));
                    for i = 1:size(obj.DYNCOEF{flowNo,iterflow},2)-1
                        obj.DYNCOEF{flowNo,iterflow}(1:3,i+1) = T'*obj.DYNCOEF{flowNo,iterflow}(1:3,i+1);
                    end
                    %解析状況を付加
                    obj.DYNCOEF{flowNo,iterflow}(1,1) = beta(iterflow);
                    obj.DYNCOEF{flowNo,iterflow}(2,1) = obj.flow{flowNo}.Mach;
                    obj.DYNCOEF{flowNo,iterflow}(3,1) = alpha(iterflow);
                    if isfield(obj.flow{flowNo},"Re")
                        ReOut = obj.flow{flowNo}.Re/1000000;
                    else
                        ReOut = 0;
                    end
                    obj.DYNCOEF{flowNo,iterflow}(4,1) = ReOut;
                end
                

                DYNCOEF = obj.DYNCOEF;
            end
        
        function obj = makeSurrogateModel(obj,alphaRange,betaRange,MachRange,Re)
            index = [3,1,11,12,6,18:20];
            sign = [1,1,1,1,1,-1,1,-1];
            [x1,x2] = ndgrid(alphaRange,betaRange);
            for iter = 1:numel(MachRange)
                obj = obj.solveFlow(x1(:),x2(:),MachRange(iter),Re);
                aerodata = obj.getAERODATA(alphaRange,betaRange,MachRange(iter),Re);
                coefdata = aerodata(:,index).*repmat(sign,[size(aerodata,1),1]);
                for i = 1:6
                    coefrbf{iter,i} = obj.RbfppMake(obj,coefdata(:,1:2),coefdata(:,i+2),1,0.01);
                end
                obj = obj.calcDynCoef(x1(:),x2(:),MachRange(iter),Re);
                outDyn = obj.getDYNCOEF(x1(:),x2(:),MachRange(iter),Re);
                dyndata = [];
                for a = 1:numel(outDyn)
                    dyndata = [dyndata;outDyn{a}(:)'];
                end
                for i = 1:size(dyndata(:,7:end),2)
                    dynrbf{iter,i} = obj.RbfppMake(obj,[dyndata(:,3),dyndata(:,1)],dyndata(:,6+i),1,0.01);
                end
            end
            alpharange2 = linspace(min(alphaRange),max(alphaRange),obj.settingUNLSI.nGriddedInterp);
            betarange2 = linspace(min(betaRange),max(betaRange),obj.settingUNLSI.nGriddedInterp);
            if numel(MachRange)>1
                machrange2 = linspace(min(MachRange),max(MachRange),obj.settingUNLSI.nGriddedInterp);
                [X1,X2,X3] = ndgrid(alpharange2,betarange2,machrange2);
                for i = 1:6
                    for a = 1:numel(alpharange2)
                        for b = 1:numel(betarange2)
                            for m = 1:numel(MachRange)
                                coefInt(m) = obj.execRbfInterp(obj,coefrbf{m,i},[X1(a,b,1),X2(a,b,1)]);
                            end
                            for c = 1:numel(machrange2)
                                data(a,b,c) = interp1(MachRange,coefInt,X3(a,b,c),'linear','extrap');
                            end
                        end
                    end
                    obj.ppCoef{i} = griddedInterpolant(X1,X2,X3,data,"linear","linear");
                end
                for a = 1:numel(alpharange2)
                    for b = 1:numel(betarange2)
                        for i = 1:size(dyndata(:,7:end),2)
                            for m = 1:numel(MachRange)
                                dyncoefInt(m,i) = obj.execRbfInterp(obj,dynrbf{m,i},[X1(a,b,1),X2(a,b,1)]);
                            end
                            for c = 1:numel(machrange2)
                                dyndata2{i}(a,b,c) = interp1(MachRange,dyncoefInt(:,i),X3(a,b,c),'linear','extrap');
                            end
                        end
                    end
                end
                for i = 1:size(dyncoefInt,2)
                    obj.ppDyn{i} = griddedInterpolant(X1,X2,X3,dyndata2{i},"linear","linear");
                end
            else
                [X1,X2] = ndgrid(alpharange2,betarange2);
                for i = 1:6
                    for a = 1:numel(alpharange2)
                        for b = 1:numel(betarange2)
                            data(a,b) = obj.execRbfInterp(obj,coefrbf{1,i},[X1(a,b),X2(a,b)]);
                        end
                    end
                    obj.ppCoef{i} = griddedInterpolant(X1,X2,data,"linear","linear");
                end
                for a = 1:numel(alpharange2)
                    for b = 1:numel(betarange2)
                        for i = 1:size(dyndata(:,7:end),2)
                            dyndata2{i}(a,b) = obj.execRbfInterp(obj,dynrbf{1,i},[X1(a,b),X2(a,b)]);
                        end
                    end
                end
                for i = 1:size(dyndata(:,7:end),2)
                    obj.ppDyn{i} = griddedInterpolant(X1,X2,dyndata2{i},"linear","linear");
                end
            end
        end
        
        function [ppCoef,ppDyn,testData] = getSurrogateModel(obj,alphaRange,betaRange,Mach,dispCoeforDyn,dispIndex,figureNo)
            ppCoef = obj.ppCoef;
            ppDyn = obj.ppDyn;
            testData = [];
            if nargin>1
                if numel(obj.ppCoef{1}.GridVectors) == 2
                    [x1,x2] = ndgrid(alphaRange,betaRange);
                    for a = 1:numel(alphaRange)
                        for b = 1:numel(betaRange)
                            if strcmpi(dispCoeforDyn,"coef")
                                testData(a,b) = ppCoef{dispIndex}(x1(a,b),x2(a,b));
                            else
                                testData(a,b) = ppDyn{dispIndex}(x1(a,b),x2(a,b));
                            end
                        end
                    end
                else
                    [x1,x2] = ndgrid(alphaRange,betaRange);
                    for a = 1:numel(alphaRange)
                        for b = 1:numel(betaRange)
                            if strcmpi(dispCoeforDyn,"coef")
                                testData(a,b) = ppCoef{dispIndex}(x1(a,b),x2(a,b),Mach);
                            else
                                testData(a,b) = ppDyn{dispIndex}(x1(a,b),x2(a,b),Mach);
                            end
                        end
                    end
                end
                figure(figureNo);clf;
                mesh(x1,x2,testData);
                xlabel("alpha(deg)");
                ylabel("beta(deg)");
                zlabel("dispCoeforDyn");
            end
        end
    end         
    methods(Static)

        function res = softplus(x,beta)
            res = 1./beta * log(1+exp(beta.*x));
        end

        function res = sigmoid(x,c,a)
            res = 1./(1 + exp(-a.*(x-c)));
        end
        
        function dcm = rod2dcm(ra,angle)
            %angleはdegreeなので注意
            dcm = zeros(3,3);
            cd  = cosd(angle);
            sd = sind(angle);
            dcm(1,1) = ra(1)^2*(1-cd)+cd;
            dcm(1,2) = ra(1)*ra(2)*(1-cd)-ra(3)*sd;
            dcm(1,3) = ra(1)*ra(3)*(1-cd)+ra(2)*sd;
            dcm(2,1) = ra(1)*ra(2)*(1-cd)+ra(3)*sd;
            dcm(2,2) = ra(2)^2*(1-cd)+cd;
            dcm(2,3) = ra(2)*ra(3)*(1-cd)-ra(1)*sd;
            dcm(3,1) = ra(1)*ra(3)*(1-cd)-ra(2)*sd;
            dcm(3,2) = ra(2)*ra(3)*(1-cd)+ra(1)*sd;
            dcm(3,3) = ra(3)^2*(1-cd)+cd;
        end

        function RHVprop = calcPropRHV(obj,T)
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            muprop =  zeros(nPanel,1);
            RHVprop = zeros(nbPanel,1);
            for propNo = 1:numel(obj.prop)
                VinfMag = dot(-obj.prop{propNo}.normal,obj.settingUNLSI.Vinf*(T*[1;0;0])');%角速度は後で
                Jratio = VinfMag / ( 2*obj.prop{propNo}.rpm*(obj.prop{propNo}.diameter/2)/60);
                Vh = -0.5*VinfMag + sqrt((0.5*VinfMag)^2 + obj.prop{propNo}.thrust/(2*obj.settingUNLSI.rho*obj.prop{propNo}.area));
                omegaProp = obj.prop{propNo}.rpm*2*pi/60;
                CT_h = obj.prop{propNo}.thrust/ ( obj.settingUNLSI.rho * obj.prop{propNo}.area * (omegaProp*obj.prop{propNo}.diameter/2)^2 );
                Vo = Vh/sqrt(1+ CT_h*log(CT_h/2)+CT_h/2); 
                eta_mom = 2/(1+sqrt(1+obj.prop{propNo}.CT));
                eta_prop = Jratio*obj.prop{propNo}.CT/obj.prop{propNo}.CP;
                dCpwakeval = (2*obj.settingUNLSI.rho*Vh^2*(omegaProp*(obj.prop{propNo}.diameter/2))^4/((omegaProp*(obj.prop{propNo}.diameter/2))^2+Vh^2)^2)/(0.5*obj.settingUNLSI.rho*VinfMag.^2)*eta_prop / eta_mom;
                Vnprop = Vo/obj.settingUNLSI.Vinf;
                for i = 1:nPanel
                    if obj.surfID(i) == obj.prop{propNo}.ID 
                        %プロペラ軸との最短距離を求める
                        r = obj.center(i,:)-obj.prop{propNo}.center;
                        dotCoef = dot(r,-obj.prop{propNo}.normal);
                        shortestPoint = obj.prop{propNo}.center - dotCoef.*obj.prop{propNo}.normal;
                        rs = norm(obj.center(i,:)-shortestPoint);
                        muprop(i,1)= (2*obj.settingUNLSI.rho*Vh^2*(omegaProp*rs)^4/((omegaProp*rs)^2+Vh^2)^2)/(0.5*obj.settingUNLSI.rho*VinfMag.^2)*eta_prop / eta_mom;
                    end
                end
                RHVprop = RHVprop + (obj.prop{propNo}.LHSi-obj.prop{propNo}.LHSo)*muprop(obj.surfID == obj.prop{propNo}.ID,1) - 2*obj.prop{propNo}.wakeLHS*dCpwakeval + (obj.prop{propNo}.RHSi)*(Vnprop*ones(size(obj.prop{propNo}.RHSi,2),1));
            end
        end

        function [VortexAr,VortexBr,VortexAc,VortexBc] = influenceMatrix(obj,rowIndex,colIndex)
            %%%%%%%%%%%%%%%%%%%影響係数の計算　パネル⇒パネル%%%%%%%%%%%%%%%
            %rowIndex : 計算する行
            %colIndex : 計算する列
            %Vortex~r : numel(rowIndex)×nbPanelの影響係数
            %Vortex~c : nbPanel×numel(rowIndex)の影響係数
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            verts = obj.tri.Points;
            con  = obj.tri.ConnectivityList(obj.paneltype==1,:);
            center = obj.center(obj.paneltype==1,:);
            normal = obj.orgNormal(obj.paneltype==1,:);

            nbPanel = size(con,1);
            nrow = numel(rowIndex);
            ncol = numel(colIndex);
            VortexAr = [];
            VortexBr = [];
            VortexAc = [];
            VortexBc = [];
            if not(isempty(rowIndex))
    
                POI.X(:,1) = center(rowIndex,1);
                POI.Y(:,1) = center(rowIndex,2);
                POI.Z(:,1) = center(rowIndex,3);
                
                c.X(1,:) = center(:,1)';
                c.Y(1,:) = center(:,2)';
                c.Z(1,:) = center(:,3)';
                n.X(1,:) = normal(:,1)';
                n.Y(1,:) = normal(:,2)';
                n.Z(1,:) = normal(:,3)';
                N1.X(1,:) = verts(con(:,1),1)';
                N1.Y(1,:) = verts(con(:,1),2)';
                N1.Z(1,:) = verts(con(:,1),3)';
                N2.X(1,:) = verts(con(:,2),1)';
                N2.Y(1,:) = verts(con(:,2),2)';
                N2.Z(1,:) = verts(con(:,2),3)';
                N3.X(1,:) = verts(con(:,3),1)';
                N3.Y(1,:) = verts(con(:,3),2)';
                N3.Z(1,:) = verts(con(:,3),3)';
                POI.X = repmat(POI.X,[1,nbPanel]);
                POI.Y = repmat(POI.Y,[1,nbPanel]);
                POI.Z = repmat(POI.Z,[1,nbPanel]);
    
                c.X = repmat(c.X,[nrow,1]);
                c.Y = repmat(c.Y,[nrow,1]);
                c.Z = repmat(c.Z,[nrow,1]);
                n.X = repmat(n.X,[nrow,1]);
                n.Y = repmat(n.Y,[nrow,1]);
                n.Z = repmat(n.Z,[nrow,1]);
                N1.X = repmat(N1.X,[nrow,1]);
                N1.Y = repmat(N1.Y,[nrow,1]);
                N1.Z = repmat(N1.Z,[nrow,1]);
                N2.X = repmat(N2.X,[nrow,1]);
                N2.Y = repmat(N2.Y,[nrow,1]);
                N2.Z = repmat(N2.Z,[nrow,1]);
                N3.X = repmat(N3.X,[nrow,1]);
                N3.Y = repmat(N3.Y,[nrow,1]);
                N3.Z = repmat(N3.Z,[nrow,1]);
    
                n12.X = (N1.X+N2.X)./2;
                n12.Y = (N1.Y+N2.Y)./2;
                n12.Z = (N1.Z+N2.Z)./2;
    
                pjk.X = POI.X-c.X;
                pjk.Y = POI.Y-c.Y;
                pjk.Z = POI.Z-c.Z;
                PN = obj.matrix_dot(pjk,n);
    
                %1回目
                a.X = POI.X-N1.X;
                a.Y = POI.Y-N1.Y;
                a.Z = POI.Z-N1.Z;
                b.X = POI.X-N2.X;
                b.Y = POI.Y-N2.Y;
                b.Z = POI.Z-N2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N2.X-N1.X;
                s.Y = N2.Y-N1.Y;
                s.Z = N2.Z-N1.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAr = phiV;
                VortexBr = srcV;
                %2回目
                a.X = POI.X-N2.X;
                a.Y = POI.Y-N2.Y;
                a.Z = POI.Z-N2.Z;
                b.X = POI.X-N3.X;
                b.Y = POI.Y-N3.Y;
                b.Z = POI.Z-N3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N3.X-N2.X;
                s.Y = N3.Y-N2.Y;
                s.Z = N3.Z-N2.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAr = VortexAr+phiV;
                VortexBr = VortexBr+srcV;
                %3回目
                a.X = POI.X-N3.X;
                a.Y = POI.Y-N3.Y;
                a.Z = POI.Z-N3.Z;
                b.X = POI.X-N1.X;
                b.Y = POI.Y-N1.Y;
                b.Z = POI.Z-N1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N1.X-N3.X;
                s.Y = N1.Y-N3.Y;
                s.Z = N1.Z-N3.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAr = VortexAr+phiV;
                VortexBr = VortexBr+srcV;
    
                %半裁の考慮
                if obj.halfmesh == 1
                    POI.Y = -POI.Y;
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    PN = obj.matrix_dot(pjk,n);
                    %1回目
                    a.X = POI.X-N1.X;
                    a.Y = POI.Y-N1.Y;
                    a.Z = POI.Z-N1.Z;
                    b.X = POI.X-N2.X;
                    b.Y = POI.Y-N2.Y;
                    b.Z = POI.Z-N2.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N2.X-N1.X;
                    s.Y = N2.Y-N1.Y;
                    s.Z = N2.Z-N1.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAr = VortexAr+phiV;
                    VortexBr = VortexBr+srcV;
                    %2回目
                    a.X = POI.X-N2.X;
                    a.Y = POI.Y-N2.Y;
                    a.Z = POI.Z-N2.Z;
                    b.X = POI.X-N3.X;
                    b.Y = POI.Y-N3.Y;
                    b.Z = POI.Z-N3.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N3.X-N2.X;
                    s.Y = N3.Y-N2.Y;
                    s.Z = N3.Z-N2.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAr = VortexAr+phiV;
                    VortexBr = VortexBr+srcV;
                    %3回目
                    a.X = POI.X-N3.X;
                    a.Y = POI.Y-N3.Y;
                    a.Z = POI.Z-N3.Z;
                    b.X = POI.X-N1.X;
                    b.Y = POI.Y-N1.Y;
                    b.Z = POI.Z-N1.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N1.X-N3.X;
                    s.Y = N1.Y-N3.Y;
                    s.Z = N1.Z-N3.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAr = VortexAr+phiV;
                    VortexBr = VortexBr+srcV;

                end
                clear POI c n N1 N2 N3
            end
            
            newCol = setdiff(1:nbPanel,rowIndex);
            %上で計算したところはコピーする
            
            if not(isempty(colIndex))
                if not(isempty(newCol))
                    POI.X(:,1) = center(newCol,1);
                    POI.Y(:,1) = center(newCol,2);
                    POI.Z(:,1) = center(newCol,3);
                    
                    c.X(1,:) = center(colIndex,1)';
                    c.Y(1,:) = center(colIndex,2)';
                    c.Z(1,:) = center(colIndex,3)';
                    n.X(1,:) = normal(colIndex,1)';
                    n.Y(1,:) = normal(colIndex,2)';
                    n.Z(1,:) = normal(colIndex,3)';
                    N1.X(1,:) = verts(con(colIndex,1),1)';
                    N1.Y(1,:) = verts(con(colIndex,1),2)';
                    N1.Z(1,:) = verts(con(colIndex,1),3)';
                    N2.X(1,:) = verts(con(colIndex,2),1)';
                    N2.Y(1,:) = verts(con(colIndex,2),2)';
                    N2.Z(1,:) = verts(con(colIndex,2),3)';
                    N3.X(1,:) = verts(con(colIndex,3),1)';
                    N3.Y(1,:) = verts(con(colIndex,3),2)';
                    N3.Z(1,:) = verts(con(colIndex,3),3)';
                    POI.X = repmat(POI.X,[1,ncol]);
                    POI.Y = repmat(POI.Y,[1,ncol]);
                    POI.Z = repmat(POI.Z,[1,ncol]);
        
                    c.X = repmat(c.X,[nbPanel-nrow,1]);
                    c.Y = repmat(c.Y,[nbPanel-nrow,1]);
                    c.Z = repmat(c.Z,[nbPanel-nrow,1]);
                    n.X = repmat(n.X,[nbPanel-nrow,1]);
                    n.Y = repmat(n.Y,[nbPanel-nrow,1]);
                    n.Z = repmat(n.Z,[nbPanel-nrow,1]);
                    N1.X = repmat(N1.X,[nbPanel-nrow,1]);
                    N1.Y = repmat(N1.Y,[nbPanel-nrow,1]);
                    N1.Z = repmat(N1.Z,[nbPanel-nrow,1]);
                    N2.X = repmat(N2.X,[nbPanel-nrow,1]);
                    N2.Y = repmat(N2.Y,[nbPanel-nrow,1]);
                    N2.Z = repmat(N2.Z,[nbPanel-nrow,1]);
                    N3.X = repmat(N3.X,[nbPanel-nrow,1]);
                    N3.Y = repmat(N3.Y,[nbPanel-nrow,1]);
                    N3.Z = repmat(N3.Z,[nbPanel-nrow,1]);
        
                    n12.X = (N1.X+N2.X)./2;
                    n12.Y = (N1.Y+N2.Y)./2;
                    n12.Z = (N1.Z+N2.Z)./2;
        
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    PN = obj.matrix_dot(pjk,n);
        
                    %1回目
                    a.X = POI.X-N1.X;
                    a.Y = POI.Y-N1.Y;
                    a.Z = POI.Z-N1.Z;
                    b.X = POI.X-N2.X;
                    b.Y = POI.Y-N2.Y;
                    b.Z = POI.Z-N2.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N2.X-N1.X;
                    s.Y = N2.Y-N1.Y;
                    s.Z = N2.Z-N1.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAcb = phiV;
                    VortexBcb = srcV;
                    %2回目
                    a.X = POI.X-N2.X;
                    a.Y = POI.Y-N2.Y;
                    a.Z = POI.Z-N2.Z;
                    b.X = POI.X-N3.X;
                    b.Y = POI.Y-N3.Y;
                    b.Z = POI.Z-N3.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N3.X-N2.X;
                    s.Y = N3.Y-N2.Y;
                    s.Z = N3.Z-N2.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAcb = VortexAcb+phiV;
                    VortexBcb = VortexBcb+srcV;
                    %3回目
                    a.X = POI.X-N3.X;
                    a.Y = POI.Y-N3.Y;
                    a.Z = POI.Z-N3.Z;
                    b.X = POI.X-N1.X;
                    b.Y = POI.Y-N1.Y;
                    b.Z = POI.Z-N1.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    s.X = N1.X-N3.X;
                    s.Y = N1.Y-N3.Y;
                    s.Z = N1.Z-N3.Z;
                    smdot = obj.matrix_dot(s,m);
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                    VortexAcb = VortexAcb+phiV;
                    VortexBcb = VortexBcb+srcV;
        
                    %半裁の考慮
                    if obj.halfmesh == 1
                        POI.Y = -POI.Y;
                        pjk.X = POI.X-c.X;
                        pjk.Y = POI.Y-c.Y;
                        pjk.Z = POI.Z-c.Z;
                        PN = obj.matrix_dot(pjk,n);
                        %1回目
                        a.X = POI.X-N1.X;
                        a.Y = POI.Y-N1.Y;
                        a.Z = POI.Z-N1.Z;
                        b.X = POI.X-N2.X;
                        b.Y = POI.Y-N2.Y;
                        b.Z = POI.Z-N2.Z;
                        anorm = obj.matrix_norm(a);
                        bnorm = obj.matrix_norm(b);
                        m = obj.getUnitVector(c,n12);
                        l = obj.matrix_cross(m,n);
                        s.X = N2.X-N1.X;
                        s.Y = N2.Y-N1.Y;
                        s.Z = N2.Z-N1.Z;
                        smdot = obj.matrix_dot(s,m);
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                        srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                        VortexAcb = VortexAcb+phiV;
                        VortexBcb = VortexBcb+srcV;
                        %2回目
                        a.X = POI.X-N2.X;
                        a.Y = POI.Y-N2.Y;
                        a.Z = POI.Z-N2.Z;
                        b.X = POI.X-N3.X;
                        b.Y = POI.Y-N3.Y;
                        b.Z = POI.Z-N3.Z;
                        anorm = obj.matrix_norm(a);
                        bnorm = obj.matrix_norm(b);
                        m = obj.getUnitVector(c,n12);
                        l = obj.matrix_cross(m,n);
                        s.X = N3.X-N2.X;
                        s.Y = N3.Y-N2.Y;
                        s.Z = N3.Z-N2.Z;
                        smdot = obj.matrix_dot(s,m);
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                        srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                        VortexAcb = VortexAcb+phiV;
                        VortexBcb = VortexBcb+srcV;
                        %3回目
                        a.X = POI.X-N3.X;
                        a.Y = POI.Y-N3.Y;
                        a.Z = POI.Z-N3.Z;
                        b.X = POI.X-N1.X;
                        b.Y = POI.Y-N1.Y;
                        b.Z = POI.Z-N1.Z;
                        anorm = obj.matrix_norm(a);
                        bnorm = obj.matrix_norm(b);
                        m = obj.getUnitVector(c,n12);
                        l = obj.matrix_cross(m,n);
                        s.X = N1.X-N3.X;
                        s.Y = N1.Y-N3.Y;
                        s.Z = N1.Z-N3.Z;
                        smdot = obj.matrix_dot(s,m);
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                        srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                        VortexAcb = VortexAcb+phiV;
                        VortexBcb = VortexBcb+srcV;
                    end
                    VortexAc(newCol,:) = VortexAcb;
                    VortexBc(newCol,:) = VortexBcb;
                end
                if not(isempty(rowIndex))
                    VortexAc(rowIndex,:) = VortexAr(:,colIndex);
                    VortexBc(rowIndex,:) = VortexBr(:,colIndex);
                end
            end
            if not(isempty(rowIndex))
                for i = 1:numel(rowIndex)
                    VortexAr(i,rowIndex(i))=-2*pi;
                end
            end
            if not(isempty(colIndex))
                for i = 1:numel(colIndex)
                    VortexAc(colIndex(i),i)=-2*pi;
                end
            end
        end

        function [VortexA] = wakeInfluenceMatrix(obj,wakeNo,edgeNo,rowIndex)
            %%%%%%%%%%%%wakeからパネルへの影響関数%%%%%%%%%%%%%%%%%
            %wakeNoとedgeNoを指定して全パネルへの影響を計算
            %
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nbPanel = size(obj.tri.ConnectivityList(obj.paneltype==1,:),1);
            center = obj.center(obj.paneltype==1,:);
            POI.X = center(rowIndex,1);
            POI.Y = center(rowIndex,2);
            POI.Z = center(rowIndex,3);
            nrow = numel(rowIndex);
            VortexA = zeros(nrow,1);
            for i = 1:1
                wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+[obj.settingUNLSI.wingWakeLength*obj.CREF*(i-1),0,0];
                wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+[obj.settingUNLSI.wingWakeLength*obj.CREF*(i-1),0,0];
                wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+[obj.settingUNLSI.wingWakeLength*obj.CREF*i,0,0];
                wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+[obj.settingUNLSI.wingWakeLength*obj.CREF*i,0,0];
                [~, ~, nbuff] = obj.vertex(wakepos(1,:),wakepos(2,:),wakepos(3,:));
                
                ulvec = obj.center(obj.wakeline{wakeNo}.lowerID(edgeNo),:)-obj.center(obj.wakeline{wakeNo}.upperID(edgeNo),:);
                uldist = dot(ulvec,nbuff);
                if uldist >= 0
                    Nw1.X = repmat(wakepos(1,1),[nrow,1]);
                    Nw1.Y = repmat(wakepos(1,2),[nrow,1]);
                    Nw1.Z = repmat(wakepos(1,3),[nrow,1]);
                    Nw2.X = repmat(wakepos(2,1),[nrow,1]);
                    Nw2.Y = repmat(wakepos(2,2),[nrow,1]);
                    Nw2.Z = repmat(wakepos(2,3),[nrow,1]);
                    Nw3.X = repmat(wakepos(3,1),[nrow,1]);
                    Nw3.Y = repmat(wakepos(3,2),[nrow,1]);
                    Nw3.Z = repmat(wakepos(3,3),[nrow,1]);
                    Nw4.X = repmat(wakepos(4,1),[nrow,1]);
                    Nw4.Y = repmat(wakepos(4,2),[nrow,1]);
                    Nw4.Z = repmat(wakepos(4,3),[nrow,1]);
                else
                    Nw1.X = repmat(wakepos(2,1),[nrow,1]);
                    Nw1.Y = repmat(wakepos(2,2),[nrow,1]);
                    Nw1.Z = repmat(wakepos(2,3),[nrow,1]);
                    Nw2.X = repmat(wakepos(1,1),[nrow,1]);
                    Nw2.Y = repmat(wakepos(1,2),[nrow,1]);
                    Nw2.Z = repmat(wakepos(1,3),[nrow,1]);
                    Nw3.X = repmat(wakepos(4,1),[nrow,1]);
                    Nw3.Y = repmat(wakepos(4,2),[nrow,1]);
                    Nw3.Z = repmat(wakepos(4,3),[nrow,1]);
                    Nw4.X = repmat(wakepos(3,1),[nrow,1]);
                    Nw4.Y = repmat(wakepos(3,2),[nrow,1]);
                    Nw4.Z = repmat(wakepos(3,3),[nrow,1]);
                end
                [b1, b2, nw] = obj.vertex([Nw1.X(1),Nw1.Y(1),Nw1.Z(1)],[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)],[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)]);
                cw = mean([[Nw1.X(1),Nw1.Y(1),Nw1.Z(1)];[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)];[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)];[Nw4.X(1),Nw4.Y(1),Nw4.Z(1)]],1);
                n.X = repmat(nw(1),[nrow,1]);
                n.Y = repmat(nw(2),[nrow,1]);
                n.Z = repmat(nw(3),[nrow,1]);
                %Amat = repmat(b1*2,[nbPanel,1]);
                c.X = repmat(cw(1),[nrow,1]);
                c.Y = repmat(cw(2),[nrow,1]);
                c.Z = repmat(cw(3),[nrow,1]);
                n12.X = (Nw3.X+Nw4.X)./2;
                n12.Y = (Nw3.Y+Nw4.Y)./2;
                n12.Z = (Nw3.Z+Nw4.Z)./2;

                pjk.X = POI.X-c.X;
                pjk.Y = POI.Y-c.Y;
                pjk.Z = POI.Z-c.Z;
                PN = obj.matrix_dot(pjk,n);
                
                %1回目
                a.X = POI.X-Nw1.X;
                a.Y = POI.Y-Nw1.Y;
                a.Z = POI.Z-Nw1.Z;
                b.X = POI.X-Nw2.X;
                b.Y = POI.Y-Nw2.Y;
                b.Z = POI.Z-Nw2.Z;
                s.X = Nw2.X-Nw1.X;
                s.Y = Nw2.Y-Nw1.Y;
                s.Z = Nw2.Z-Nw1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                VortexA = VortexA+phiV;
                %2回目
                a.X = POI.X-Nw2.X;
                a.Y = POI.Y-Nw2.Y;
                a.Z = POI.Z-Nw2.Z;
                b.X = POI.X-Nw3.X;
                b.Y = POI.Y-Nw3.Y;
                b.Z = POI.Z-Nw3.Z;
                s.X = Nw3.X-Nw2.X;
                s.Y = Nw3.Y-Nw2.Y;
                s.Z = Nw3.Z-Nw2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                VortexA = VortexA+phiV;
                %3回目
                a.X = POI.X-Nw3.X;
                a.Y = POI.Y-Nw3.Y;
                a.Z = POI.Z-Nw3.Z;
                b.X = POI.X-Nw4.X;
                b.Y = POI.Y-Nw4.Y;
                b.Z = POI.Z-Nw4.Z;
                s.X = Nw4.X-Nw3.X;
                s.Y = Nw4.Y-Nw3.Y;
                s.Z = Nw4.Z-Nw3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                VortexA = VortexA+phiV;
                %4回目
                a.X = POI.X-Nw4.X;
                a.Y = POI.Y-Nw4.Y;
                a.Z = POI.Z-Nw4.Z;
                b.X = POI.X-Nw1.X;
                b.Y = POI.Y-Nw1.Y;
                b.Z = POI.Z-Nw1.Z;
                s.X = Nw1.X-Nw4.X;
                s.Y = Nw1.Y-Nw4.Y;
                s.Z = Nw1.Z-Nw4.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                VortexA = VortexA+phiV;

                %半裁
                if obj.halfmesh == 1
                    POI.Y = -POI.Y;
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    PN = obj.matrix_dot(pjk,n);
                    
                    %1回目
                    a.X = POI.X-Nw1.X;
                    a.Y = POI.Y-Nw1.Y;
                    a.Z = POI.Z-Nw1.Z;
                    b.X = POI.X-Nw2.X;
                    b.Y = POI.Y-Nw2.Y;
                    b.Z = POI.Z-Nw2.Z;
                    s.X = Nw2.X-Nw1.X;
                    s.Y = Nw2.Y-Nw1.Y;
                    s.Z = Nw2.Z-Nw1.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    
                    VortexA = VortexA+phiV;
                    %2回目
                    a.X = POI.X-Nw2.X;
                    a.Y = POI.Y-Nw2.Y;
                    a.Z = POI.Z-Nw2.Z;
                    b.X = POI.X-Nw3.X;
                    b.Y = POI.Y-Nw3.Y;
                    b.Z = POI.Z-Nw3.Z;
                    s.X = Nw3.X-Nw2.X;
                    s.Y = Nw3.Y-Nw2.Y;
                    s.Z = Nw3.Z-Nw2.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    VortexA = VortexA+phiV;
                    %3回目
                    a.X = POI.X-Nw3.X;
                    a.Y = POI.Y-Nw3.Y;
                    a.Z = POI.Z-Nw3.Z;
                    b.X = POI.X-Nw4.X;
                    b.Y = POI.Y-Nw4.Y;
                    b.Z = POI.Z-Nw4.Z;
                    s.X = Nw4.X-Nw3.X;
                    s.Y = Nw4.Y-Nw3.Y;
                    s.Z = Nw4.Z-Nw3.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    VortexA = VortexA+phiV;
                    %4回目
                    a.X = POI.X-Nw4.X;
                    a.Y = POI.Y-Nw4.Y;
                    a.Z = POI.Z-Nw4.Z;
                    b.X = POI.X-Nw1.X;
                    b.Y = POI.Y-Nw1.Y;
                    b.Z = POI.Z-Nw1.Z;
                    s.X = Nw1.X-Nw4.X;
                    s.Y = Nw1.Y-Nw4.Y;
                    s.Z = Nw1.Z-Nw4.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    VortexA = VortexA+phiV;
                end
            end

        end
        
        function [VortexAi,VortexBi,VortexAo,VortexBo,wakeVA] = propellerInfluenceMatrix(obj,propNo,propWakeLength)
            ID = obj.prop{propNo}.ID;
            verts = obj.tri.Points;
            con  = obj.tri.ConnectivityList(obj.paneltype==1,:);
            %normal = obj.orgNormal(obj.paneltype==1,:);

            nbPanel = size(con,1);
            ncol = sum(obj.surfID == ID);
    
            POI.X(:,1) = obj.center(obj.paneltype==1,1);
            POI.Y(:,1) = obj.center(obj.paneltype==1,2);
            POI.Z(:,1) = obj.center(obj.paneltype==1,3);
            
            c.X(1,:) = obj.center(obj.surfID == ID,1)';
            c.Y(1,:) = obj.center(obj.surfID == ID,2)';
            c.Z(1,:) = obj.center(obj.surfID == ID,3)';
            n.X(1,:) = obj.orgNormal(obj.surfID == ID,1)';
            n.Y(1,:) = obj.orgNormal(obj.surfID == ID,2)';
            n.Z(1,:) = obj.orgNormal(obj.surfID == ID,3)';
            N1.X(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,1),1)';
            N1.Y(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,1),2)';
            N1.Z(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,1),3)';
            N2.X(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,2),1)';
            N2.Y(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,2),2)';
            N2.Z(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,2),3)';
            N3.X(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,3),1)';
            N3.Y(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,3),2)';
            N3.Z(1,:) = verts(obj.tri.ConnectivityList(obj.surfID == ID,3),3)';
            POI.X = repmat(POI.X,[1,ncol]);
            POI.Y = repmat(POI.Y,[1,ncol]);
            POI.Z = repmat(POI.Z,[1,ncol]);
            
            %流入側
            c.X = repmat(c.X,[nbPanel,1]);
            c.Y = repmat(c.Y,[nbPanel,1]);
            c.Z = repmat(c.Z,[nbPanel,1]);
            n.X = repmat(n.X,[nbPanel,1]);
            n.Y = repmat(n.Y,[nbPanel,1]);
            n.Z = repmat(n.Z,[nbPanel,1]);
            N1.X = repmat(N1.X,[nbPanel,1]);
            N1.Y = repmat(N1.Y,[nbPanel,1]);
            N1.Z = repmat(N1.Z,[nbPanel,1]);
            N2.X = repmat(N2.X,[nbPanel,1]);
            N2.Y = repmat(N2.Y,[nbPanel,1]);
            N2.Z = repmat(N2.Z,[nbPanel,1]);
            N3.X = repmat(N3.X,[nbPanel,1]);
            N3.Y = repmat(N3.Y,[nbPanel,1]);
            N3.Z = repmat(N3.Z,[nbPanel,1]);

            n12.X = (N1.X+N2.X)./2;
            n12.Y = (N1.Y+N2.Y)./2;
            n12.Z = (N1.Z+N2.Z)./2;

            pjk.X = POI.X-c.X;
            pjk.Y = POI.Y-c.Y;
            pjk.Z = POI.Z-c.Z;
            PN = obj.matrix_dot(pjk,n);

            %1回目
            a.X = POI.X-N1.X;
            a.Y = POI.Y-N1.Y;
            a.Z = POI.Z-N1.Z;
            b.X = POI.X-N2.X;
            b.Y = POI.Y-N2.Y;
            b.Z = POI.Z-N2.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N2.X-N1.X;
            s.Y = N2.Y-N1.Y;
            s.Z = N2.Z-N1.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAi = phiV;
            VortexBi = srcV;
            %2回目
            a.X = POI.X-N2.X;
            a.Y = POI.Y-N2.Y;
            a.Z = POI.Z-N2.Z;
            b.X = POI.X-N3.X;
            b.Y = POI.Y-N3.Y;
            b.Z = POI.Z-N3.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N3.X-N2.X;
            s.Y = N3.Y-N2.Y;
            s.Z = N3.Z-N2.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAi = VortexAi+phiV;
            VortexBi = VortexBi+srcV;
            %3回目
            a.X = POI.X-N3.X;
            a.Y = POI.Y-N3.Y;
            a.Z = POI.Z-N3.Z;
            b.X = POI.X-N1.X;
            b.Y = POI.Y-N1.Y;
            b.Z = POI.Z-N1.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N1.X-N3.X;
            s.Y = N1.Y-N3.Y;
            s.Z = N1.Z-N3.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAi = VortexAi+phiV;
            VortexBi = VortexBi+srcV;

            %流出側
            clear n;
            n.X(1,:) = -obj.orgNormal(obj.surfID == ID,1)';
            n.Y(1,:) = -obj.orgNormal(obj.surfID == ID,2)';
            n.Z(1,:) = -obj.orgNormal(obj.surfID == ID,3)';
            n.X = repmat(n.X,[nbPanel,1]);
            n.Y = repmat(n.Y,[nbPanel,1]);
            n.Z = repmat(n.Z,[nbPanel,1]);
            PN = obj.matrix_dot(pjk,n);

            %1回目
            a.X = POI.X-N2.X;
            a.Y = POI.Y-N2.Y;
            a.Z = POI.Z-N2.Z;
            b.X = POI.X-N1.X;
            b.Y = POI.Y-N1.Y;
            b.Z = POI.Z-N1.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N1.X-N2.X;
            s.Y = N1.Y-N2.Y;
            s.Z = N1.Z-N2.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAo = phiV;
            VortexBo = srcV;
            %2回目
            a.X = POI.X-N3.X;
            a.Y = POI.Y-N3.Y;
            a.Z = POI.Z-N3.Z;
            b.X = POI.X-N2.X;
            b.Y = POI.Y-N2.Y;
            b.Z = POI.Z-N2.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N2.X-N3.X;
            s.Y = N2.Y-N3.Y;
            s.Z = N2.Z-N3.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAo = VortexAo+phiV;
            VortexBo = VortexBo+srcV;
            %3回目
            a.X = POI.X-N1.X;
            a.Y = POI.Y-N1.Y;
            a.Z = POI.Z-N1.Z;
            b.X = POI.X-N3.X;
            b.Y = POI.Y-N3.Y;
            b.Z = POI.Z-N3.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            m = obj.getUnitVector(c,n12);
            l = obj.matrix_cross(m,n);
            s.X = N3.X-N1.X;
            s.Y = N3.Y-N1.Y;
            s.Z = N3.Z-N1.Z;
            smdot = obj.matrix_dot(s,m);
            snorm = obj.matrix_norm(s);
            Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            PB = PA-Al.*smdot;
            phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
            VortexAo = VortexAo+phiV;
            VortexBo = VortexBo+srcV;

            %半裁の考慮
            if obj.halfmesh == 1
                clear n;
                POI.Y = -POI.Y;
                pjk.X = POI.X-c.X;
                pjk.Y = POI.Y-c.Y;
                pjk.Z = POI.Z-c.Z;
                %流入
                n.X(1,:) = obj.orgNormal(obj.surfID == ID,1)';
                n.Y(1,:) = obj.orgNormal(obj.surfID == ID,2)';
                n.Z(1,:) = obj.orgNormal(obj.surfID == ID,3)';
                n.X = repmat(n.X,[nbPanel,1]);
                n.Y = repmat(n.Y,[nbPanel,1]);
                n.Z = repmat(n.Z,[nbPanel,1]);
                PN = obj.matrix_dot(pjk,n);
                %1回目
                a.X = POI.X-N1.X;
                a.Y = POI.Y-N1.Y;
                a.Z = POI.Z-N1.Z;
                b.X = POI.X-N2.X;
                b.Y = POI.Y-N2.Y;
                b.Z = POI.Z-N2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N2.X-N1.X;
                s.Y = N2.Y-N1.Y;
                s.Z = N2.Z-N1.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAi = VortexAi+phiV;
                VortexBi = VortexBi+srcV;
                %2回目
                a.X = POI.X-N2.X;
                a.Y = POI.Y-N2.Y;
                a.Z = POI.Z-N2.Z;
                b.X = POI.X-N3.X;
                b.Y = POI.Y-N3.Y;
                b.Z = POI.Z-N3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N3.X-N2.X;
                s.Y = N3.Y-N2.Y;
                s.Z = N3.Z-N2.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAi = VortexAi+phiV;
                VortexBi = VortexBi+srcV;
                %3回目
                a.X = POI.X-N3.X;
                a.Y = POI.Y-N3.Y;
                a.Z = POI.Z-N3.Z;
                b.X = POI.X-N1.X;
                b.Y = POI.Y-N1.Y;
                b.Z = POI.Z-N1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N1.X-N3.X;
                s.Y = N1.Y-N3.Y;
                s.Z = N1.Z-N3.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAi = VortexAi+phiV;
                VortexBi = VortexBi+srcV;

                %流出
                clear n
                n.X(1,:) = -obj.orgNormal(obj.surfID == ID,1)';
                n.Y(1,:) = -obj.orgNormal(obj.surfID == ID,2)';
                n.Z(1,:) = -obj.orgNormal(obj.surfID == ID,3)';
                n.X = repmat(n.X,[nbPanel,1]);
                n.Y = repmat(n.Y,[nbPanel,1]);
                n.Z = repmat(n.Z,[nbPanel,1]);
                PN = obj.matrix_dot(pjk,n);
                %1回目
                a.X = POI.X-N2.X;
                a.Y = POI.Y-N2.Y;
                a.Z = POI.Z-N2.Z;
                b.X = POI.X-N1.X;
                b.Y = POI.Y-N1.Y;
                b.Z = POI.Z-N1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N1.X-N2.X;
                s.Y = N1.Y-N2.Y;
                s.Z = N1.Z-N2.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAo = VortexAo+phiV;
                VortexBo = VortexBo+srcV;
                %2回目
                a.X = POI.X-N3.X;
                a.Y = POI.Y-N3.Y;
                a.Z = POI.Z-N3.Z;
                b.X = POI.X-N2.X;
                b.Y = POI.Y-N2.Y;
                b.Z = POI.Z-N2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N2.X-N3.X;
                s.Y = N2.Y-N3.Y;
                s.Z = N2.Z-N3.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAo = VortexAo+phiV;
                VortexBo = VortexBo+srcV;
                %3回目
                a.X = POI.X-N1.X;
                a.Y = POI.Y-N1.Y;
                a.Z = POI.Z-N1.Z;
                b.X = POI.X-N3.X;
                b.Y = POI.Y-N3.Y;
                b.Z = POI.Z-N3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                s.X = N3.X-N1.X;
                s.Y = N3.Y-N1.Y;
                s.Z = N3.Z-N1.Z;
                smdot = obj.matrix_dot(s,m);
                snorm = obj.matrix_norm(s);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                srcV = Al.*(log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm)-PN.*phiV;
                VortexAo = VortexAo+phiV;
                VortexBo = VortexBo+srcV;
            end

            %プロペラwakeの計算
            % 円周上の点の座標を計算する
            if obj.prop{propNo}.XZsliced == 1
                theta = linspace(pi, 0, obj.settingUNLSI.nPropWake);
            else
                theta = linspace(2*pi, 0, obj.settingUNLSI.nPropWake);
            end
            cirp = obj.prop{propNo}.diameter/2 * [cos(theta); sin(theta); zeros(1, obj.settingUNLSI.nPropWake)];
            
            % 円周上の点を中心座標と法線ベクトルを使用して回転・平行移動する
            rotation_matrix = vrrotvec2mat(vrrotvec([0, 0, 1], obj.prop{propNo}.normal)); % 回転行列を計算
            pellerEdge = (rotation_matrix * cirp + obj.prop{propNo}.center')';
                        
            POI.X = obj.center(obj.paneltype==1,1);
            POI.Y = obj.center(obj.paneltype==1,2);
            POI.Z = obj.center(obj.paneltype==1,3);
            pellerNormal = mean(obj.orgNormal(obj.surfID == ID,:),1);
            wakeVA = zeros(nbPanel,1);
            for wakeNo = 1:size(pellerEdge,1)-1
                wakepos(1,:) = pellerEdge(wakeNo,:)+[0,0,0];
                wakepos(2,:) = pellerEdge(wakeNo+1,:)+[0,0,0];
                wakepos(3,:) = pellerEdge(wakeNo+1,:)-propWakeLength*obj.prop{propNo}.diameter.*pellerNormal;
                wakepos(4,:) = pellerEdge(wakeNo,:)-propWakeLength*obj.prop{propNo}.diameter.*pellerNormal;
                [~, ~, nbuff] = obj.vertex(wakepos(1,:),wakepos(2,:),wakepos(3,:));
                Nw1.X = repmat(wakepos(2,1),[nbPanel,1]);
                Nw1.Y = repmat(wakepos(2,2),[nbPanel,1]);
                Nw1.Z = repmat(wakepos(2,3),[nbPanel,1]);
                Nw2.X = repmat(wakepos(1,1),[nbPanel,1]);
                Nw2.Y = repmat(wakepos(1,2),[nbPanel,1]);
                Nw2.Z = repmat(wakepos(1,3),[nbPanel,1]);
                Nw3.X = repmat(wakepos(4,1),[nbPanel,1]);
                Nw3.Y = repmat(wakepos(4,2),[nbPanel,1]);
                Nw3.Z = repmat(wakepos(4,3),[nbPanel,1]);
                Nw4.X = repmat(wakepos(3,1),[nbPanel,1]);
                Nw4.Y = repmat(wakepos(3,2),[nbPanel,1]);
                Nw4.Z = repmat(wakepos(3,3),[nbPanel,1]);
                [b1, b2, nw] = obj.vertex([Nw1.X(1),Nw1.Y(1),Nw1.Z(1)],[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)],[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)]);
                cw = mean([[Nw1.X(1),Nw1.Y(1),Nw1.Z(1)];[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)];[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)];[Nw4.X(1),Nw4.Y(1),Nw4.Z(1)]],1);
                n.X = repmat(nw(1),[nbPanel,1]);
                n.Y = repmat(nw(2),[nbPanel,1]);
                n.Z = repmat(nw(3),[nbPanel,1]);
                %Amat = repmat(b1*2,[nbPanel,1]);
                c.X = repmat(cw(1),[nbPanel,1]);
                c.Y = repmat(cw(2),[nbPanel,1]);
                c.Z = repmat(cw(3),[nbPanel,1]);
                n12.X = (Nw3.X+Nw4.X)./2;
                n12.Y = (Nw3.Y+Nw4.Y)./2;
                n12.Z = (Nw3.Z+Nw4.Z)./2;

                pjk.X = POI.X-c.X;
                pjk.Y = POI.Y-c.Y;
                pjk.Z = POI.Z-c.Z;
                PN = obj.matrix_dot(pjk,n);
                
                %1回目
                a.X = POI.X-Nw1.X;
                a.Y = POI.Y-Nw1.Y;
                a.Z = POI.Z-Nw1.Z;
                b.X = POI.X-Nw2.X;
                b.Y = POI.Y-Nw2.Y;
                b.Z = POI.Z-Nw2.Z;
                s.X = Nw2.X-Nw1.X;
                s.Y = Nw2.Y-Nw1.Y;
                s.Z = Nw2.Z-Nw1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                wakeVA = wakeVA+phiV;
                %2回目
                a.X = POI.X-Nw2.X;
                a.Y = POI.Y-Nw2.Y;
                a.Z = POI.Z-Nw2.Z;
                b.X = POI.X-Nw3.X;
                b.Y = POI.Y-Nw3.Y;
                b.Z = POI.Z-Nw3.Z;
                s.X = Nw3.X-Nw2.X;
                s.Y = Nw3.Y-Nw2.Y;
                s.Z = Nw3.Z-Nw2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                wakeVA = wakeVA+phiV;
                %3回目
                a.X = POI.X-Nw3.X;
                a.Y = POI.Y-Nw3.Y;
                a.Z = POI.Z-Nw3.Z;
                b.X = POI.X-Nw4.X;
                b.Y = POI.Y-Nw4.Y;
                b.Z = POI.Z-Nw4.Z;
                s.X = Nw4.X-Nw3.X;
                s.Y = Nw4.Y-Nw3.Y;
                s.Z = Nw4.Z-Nw3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                wakeVA = wakeVA+phiV;
                %4回目
                a.X = POI.X-Nw4.X;
                a.Y = POI.Y-Nw4.Y;
                a.Z = POI.Z-Nw4.Z;
                b.X = POI.X-Nw1.X;
                b.Y = POI.Y-Nw1.Y;
                b.Z = POI.Z-Nw1.Z;
                s.X = Nw1.X-Nw4.X;
                s.Y = Nw1.Y-Nw4.Y;
                s.Z = Nw1.Z-Nw4.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
                m = obj.getUnitVector(c,n12);
                l = obj.matrix_cross(m,n);
                smdot = obj.matrix_dot(s,m);
                Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                PB = PA-Al.*smdot;
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                wakeVA = wakeVA+phiV;

                %半裁
                if obj.halfmesh == 1
                    POI.Y = -POI.Y;
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    PN = obj.matrix_dot(pjk,n);
                    
                    %1回目
                    a.X = POI.X-Nw1.X;
                    a.Y = POI.Y-Nw1.Y;
                    a.Z = POI.Z-Nw1.Z;
                    b.X = POI.X-Nw2.X;
                    b.Y = POI.Y-Nw2.Y;
                    b.Z = POI.Z-Nw2.Z;
                    s.X = Nw2.X-Nw1.X;
                    s.Y = Nw2.Y-Nw1.Y;
                    s.Z = Nw2.Z-Nw1.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    
                    wakeVA = wakeVA+phiV;
                    %2回目
                    a.X = POI.X-Nw2.X;
                    a.Y = POI.Y-Nw2.Y;
                    a.Z = POI.Z-Nw2.Z;
                    b.X = POI.X-Nw3.X;
                    b.Y = POI.Y-Nw3.Y;
                    b.Z = POI.Z-Nw3.Z;
                    s.X = Nw3.X-Nw2.X;
                    s.Y = Nw3.Y-Nw2.Y;
                    s.Z = Nw3.Z-Nw2.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    wakeVA = wakeVA+phiV;
                    %3回目
                    a.X = POI.X-Nw3.X;
                    a.Y = POI.Y-Nw3.Y;
                    a.Z = POI.Z-Nw3.Z;
                    b.X = POI.X-Nw4.X;
                    b.Y = POI.Y-Nw4.Y;
                    b.Z = POI.Z-Nw4.Z;
                    s.X = Nw4.X-Nw3.X;
                    s.Y = Nw4.Y-Nw3.Y;
                    s.Z = Nw4.Z-Nw3.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    wakeVA = wakeVA+phiV;
                    %4回目
                    a.X = POI.X-Nw4.X;
                    a.Y = POI.Y-Nw4.Y;
                    a.Z = POI.Z-Nw4.Z;
                    b.X = POI.X-Nw1.X;
                    b.Y = POI.Y-Nw1.Y;
                    b.Z = POI.Z-Nw1.Z;
                    s.X = Nw1.X-Nw4.X;
                    s.Y = Nw1.Y-Nw4.Y;
                    s.Z = Nw1.Z-Nw4.Z;
                    anorm = obj.matrix_norm(a);
                    bnorm = obj.matrix_norm(b);
                    m = obj.getUnitVector(c,n12);
                    l = obj.matrix_cross(m,n);
                    smdot = obj.matrix_dot(s,m);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                    wakeVA = wakeVA+phiV;
                end
            end
        end
        
        function Q_ij=Calc_Q(y,z,phi,dS,halfmesh)
            yd_ij = zeros(numel(y),numel(y));
            zd_ij = zeros(numel(y),numel(y));
            ydd_ij = zeros(numel(y),numel(y));
            zdd_ij = zeros(numel(y),numel(y));
            R_2_Pij = zeros(numel(y),numel(y));
            R_2_Mij = zeros(numel(y),numel(y));
            Rd_2_Pij = zeros(numel(y),numel(y));
            Rd_2_Mij = zeros(numel(y),numel(y));
            Q_ij_1 = zeros(numel(y),numel(y));
            Q_ij_2 = zeros(numel(y),numel(y));
            Q_ij_3 = zeros(numel(y),numel(y));
            Q_ij_4 = zeros(numel(y),numel(y));
            Q_ij = zeros(numel(y),numel(y));
	        for i=1:numel(y)
		        for j=1:numel(y)
			        yd_ij(i,j)=(y(i)-y(j))*cos(phi(j))+(z(i)-z(j))*sin(phi(j));
			        zd_ij(i,j)=-(y(i)-y(j))*sin(phi(j))+(z(i)-z(j))*cos(phi(j));
			        R_2_Pij(i,j)=(yd_ij(i,j)-dS(j))^2+zd_ij(i,j)^2;
			        R_2_Mij(i,j)=(yd_ij(i,j)+dS(j))^2+zd_ij(i,j)^2;
			        Q_ij_1(i,j)=((yd_ij(i,j)-dS(j))/R_2_Pij(i,j)-(yd_ij(i,j)+dS(j))/R_2_Mij(i,j))*cos(phi(i)-phi(j));
			        Q_ij_2(i,j)=((zd_ij(i,j))/R_2_Pij(i,j)-(zd_ij(i,j))/R_2_Mij(i,j))*sin(phi(i)-phi(j));
			        
                    if halfmesh == 1
                        ydd_ij(i,j)=(y(i)+y(j))*cos(phi(j))-(z(i)-z(j))*sin(phi(j));
			            zdd_ij(i,j)=(y(i)+y(j))*sin(phi(j))+(z(i)-z(j))*cos(phi(j));
			            Rd_2_Pij(i,j)=(ydd_ij(i,j)+dS(j))^2+zdd_ij(i,j)^2;
			            Rd_2_Mij(i,j)=(ydd_ij(i,j)-dS(j))^2+zdd_ij(i,j)^2;
                        Q_ij_3(i,j)=((ydd_ij(i,j)-dS(j))/Rd_2_Mij(i,j)-(ydd_ij(i,j)+dS(j))/Rd_2_Pij(i,j))*cos(phi(i)+phi(j));
			            Q_ij_4(i,j)=((zdd_ij(i,j))/Rd_2_Mij(i,j)-(zdd_ij(i,j))/Rd_2_Pij(i,j))*sin(phi(i)+phi(j));
                    end
                    if halfmesh == 1
			            Q_ij(i,j)=-1/2/pi*(Q_ij_1(i,j)+Q_ij_2(i,j)+Q_ij_3(i,j)+Q_ij_4(i,j));
                    else
                        Q_ij(i,j)=-1/2/pi*(Q_ij_1(i,j)+Q_ij_2(i,j));
                    end
		        end
	        end
        end

        function res = pmsolve(M1,M2,nu,kappa)
            M1nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M1^2-1)))-atan(sqrt(M1^2-1));
            M2nu = M1nu+nu;
            res = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M2^2-1)))-atan(sqrt(M2^2-1))-M2nu;
        end

        function [p, info] = bisection(fun,a,b)
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

        function [darea, dvolume, n] = vertex(P0,P1,P2)
            %座標を与えてSTLの繰り返し部分を追加書き込みする
            x10 = P1(1)-P0(1);y10 = P1(2)-P0(2);z10 = P1(3)-P0(3);
            x20 = P2(1)-P0(1);y20 = P2(2)-P0(2);z20 = P2(3)-P0(3);
            nx=y10*z20-z10*y20;
            ny=z10*x20-x10*z20;
            nz=x10*y20-y10*x20;
            N = norm([nx,ny,nz]);
            darea = sqrt(nx*nx+ny*ny+nz*nz)/2;
            dvolume = (nx*P0(1)+ny*P0(2)+nz*P0(3))/6;
            n = [nx,ny,nz]./N;
        end

        function m = getUnitVector(p1,p2)
            m.X = p2.X-p1.X;
            m.Y = p2.Y-p1.Y;
            m.Z = p2.Z-p1.Z;
            no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
            index = find(no);
            m.X(index) = m.X(index) ./ no(index);
            m.Y(index) = m.Y(index) ./ no(index);
            m.Z(index) = m.Z(index) ./ no(index);
        end
        
        function no = matrix_norm(m)
            no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
        end
        
        function md = matrix_dot(A,B)
            md = A.X.*B.X+A.Y.*B.Y+A.Z.*B.Z;
        end
        
        function mc = matrix_cross(A,B)
            mc.X = A.Y.*B.Z-A.Z.*B.Y;
            mc.Y = A.Z.*B.X-A.X.*B.Z;
            mc.Z = A.X.*B.Y-A.Y.*B.X;
        end
        
        function [pp] = RbfppMake(obj,xd,fd,rbfMode,r0,colmax)
            if rbfMode ==4
                phi = @(r,r0)obj.phi4(r,r0);
            elseif rbfMode == 1
                phi = @(r,r0)obj.phi1(r,r0);
            elseif rbfMode == 2
                phi = @(r,r0)obj.phi2(r,r0);
            else
                error('未実装')
            end
        
            if nargin > 5
                colmax = colmax(:);
                scaleShift = zeros(size(colmax,1),1);
                scaleWeight = colmax-scaleShift;
            else
                for i = 1:size(xd,2)
                    scaleShift(i,1) = min(xd(:,i));
                    colmax(i,1) = max(xd(:,i));
                    scaleWeight(i,1) = colmax(i,1)-scaleShift(i,1);
                    if scaleWeight(i,1) == 0
                        scaleWeight(i,1) = 1;
                    end 
                end
            end
            %xdのスケールの変更
            for i = 1:size(xd,2)
                xd(:,i) = (xd(:,i)-scaleShift(i,1))./scaleWeight(i,1);
            end
            
            X = xd';
            H = sum(xd.^2,2);
            H = repmat(H,[1,size(X,2)]);
            r = sqrt(H'-2.*X'*X+H);
            a = phi(r,r0);
            invR = pinv(a);
            w = invR*fd;
            pp.w = w;
            pp.rbfMode = rbfMode;
            pp.nSample = size(xd,1);
            pp.nDesign = size(xd,2);
            pp.val_samp = xd;
            pp.res_samp = fd;
            pp.R0 = r0;
            pp.scaleShift = scaleShift;
            pp.scaleWeight = scaleWeight;
        end
    
        function  phi = phi1(r,r0)
            phi = sqrt(r.*r+r0.*r0);
        end
        
        function  phi = phi2(r,r0)
            phi = 1./sqrt(r.*r+r0.*r0);
        end
        
        function phi = phi4(r,r0)
            phi = exp(-0.5.*r.^2/r0.^2);
        end

       function fi = execRbfInterp(obj,pp,val_interp)
            nSamp = size(pp.val_samp,1);
            nInterp = size(val_interp,1);
            val_interp= (val_interp-repmat(pp.scaleShift(:)',[nInterp,1]))./repmat(pp.scaleWeight(:)',[nInterp,1]);
            Xi = sum(val_interp.^2,2);
            H1 = repmat(Xi,[1,nSamp]);
            Xs = sum(pp.val_samp.^2,2);
            H2 = repmat(Xs',[nInterp,1]);
            M = val_interp*pp.val_samp';
            r = sqrt(H1-2.*M+H2);
            fi = obj.phi1(r,pp.R0)*pp.w(:);
        end

    end

end