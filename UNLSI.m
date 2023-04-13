classdef UNLSI

    %%%%%残タスク
    %solveFlowForAdjointの実装
    %calcDyncoefの実装
    %%%%%

    properties
        tri %MATLAB triangulationクラス
        surfID %各パネルのタグ番号
        wakeline %wakeパネルの設定
        halfmesh %半裁かどうか(半裁:1)
        SREF %基準面積
        CREF %縦方向基準長
        BREF %横方向基準長
        XYZREF %回転中心
        paneltype %各パネルのタイプ 1:body 2:base 3:structure
        IndexPanel2Solver %パネルのインデックス⇒ソルバー上でのインデックス
        VindWake
        wakePanelLength
        nWake
        normal %各パネルの法線ベクトル
        center %各パネルの中心
        area %各パネルの面積
        flow %各flowconditionが格納されたセル配列
        cluster %パネルクラスターの情報
        LHS %パネル法連立方程式の左辺行列
        RHS %パネル法連立方程式の右辺行列
        mu2v %ポテンシャル⇒機体表面速度への変換行列
        Cp %圧力係数
        Cfe %表面摩擦係数
        AERODATA %結果の格納
        LLT
    end

    methods(Access = public)
        function obj = UNLSI(verts,connectivity,surfID,wakelineID,halfmesh)
            %%%%%%%%%%%Constructor%%%%%%%%%%%%%
            %verts:頂点座標
            %conectivity:各パネルの頂点ID
            %surfID:パネルのID
            %wakelineID:wakeをつけるエッジ上の頂点ID配列を要素にもつセル配列
            %halfmesh:半裁メッシュの場合に1を指定
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %データの格納、パネル法線、面積の計算
            obj.tri = triangulation(connectivity,verts);
            obj.surfID = surfID;
            for i = 1:numel(wakelineID)
                obj.wakeline{i}.edge = wakelineID{i};
            end
            obj.paneltype = ones(size(connectivity,1),1);
            obj.IndexPanel2Solver = 1:numel(obj.paneltype);
            obj.halfmesh = halfmesh;
            obj.LLT.n_interp = 10;
            obj.flow = {};
            obj.SREF = 1;
            obj.CREF = 1;
            obj.BREF = 1;
            obj.XYZREF = [0,0,0];
            obj.cluster = cell([1,numel(obj.paneltype)]);
            obj.area = zeros(numel(obj.paneltype),1);
            obj.normal = zeros(numel(obj.paneltype),3);
            obj.center = zeros(numel(obj.paneltype),3);
            obj.Cp = {};
            obj.Cfe = {};
            for i = 1:numel(obj.paneltype)
                [obj.area(i,1),~ , obj.normal(i,:)] = obj.vertex(verts(connectivity(i,1),:),verts(connectivity(i,2),:),verts(connectivity(i,3),:));
                obj.center(i,:) = [mean(verts(obj.tri.ConnectivityList(i,:),1)),mean(verts(obj.tri.ConnectivityList(i,:),2)),mean(verts(obj.tri.ConnectivityList(i,:),3))];
            end

            
            %半裁メッシュの境界表面上のwakeを削除
            if halfmesh == 1
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
                    if numel(obj.wakeline{iter}.edge) == 1
                        deleteiter = [deleteiter,iter];
                    end
                end
                obj.wakeline(deleteiter) = [];
            end
            %wakeのつくパネルIDを特定
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{wakeNo}.edge(edgeNo),obj.wakeline{wakeNo}.edge(edgeNo+1));
                    phivert = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),2:3)-obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),2:3);
                    phiwake = atan2d(phivert(2),phivert(1));
                    if abs(phiwake)<90
                        obj.wakeline{wakeNo}.upperID(edgeNo) = attachpanel{1}(1);
                        obj.wakeline{wakeNo}.lowerID(edgeNo) = attachpanel{1}(2);
                    else
                        obj.wakeline{wakeNo}.upperID(edgeNo) = attachpanel{1}(2);
                        obj.wakeline{wakeNo}.lowerID(edgeNo) = attachpanel{1}(1);
                    end
                end
            end
            obj.checkMesh(sqrt(eps),"warning");
        end
       

        function obj = setREFS(obj,SREF,BREF,CREF)
            %%%%%%%%%%%Referentials Setting%%%%%%%%%%%%%
            %SREF:基準面積
            %BREF:横方向基準長(eg.翼幅)
            %CREF:縦方向基準長(eg.MAC,機体全長
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.SREF = SREF;
            obj.CREF = CREF;
            obj.BREF = BREF;
        end

        function obj = setRotationCenter(obj,XYZREF)
            %%%%%%%%%%%Rotation Center Setting%%%%%%%%%%%%%
            %回転中心を設定する。
            %回転中心：モーメント計算の基準位置,主流角速度の回転中心
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.XYZREF = XYZREF(:)';
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
            axis equal;xlabel("x");ylabel("y");zlabel("z");
            drawnow();pause(0.1);

        end

        function obj = setMesh(obj,verts,connectivity,surfID,wakelineID)
            obj.tri = triangulation(connectivity,verts);
            obj.surfID = surfID;
            for i = 1:numel(wakelineID)
                obj.wakeline{i}.edge = wakelineID{i};
            end
            obj.paneltype = ones(size(connectivity,1),1);
            obj.IndexPanel2Solver = 1:numel(obj.paneltype);
            obj.cluster = cell([1,numel(obj.paneltype)]);
            obj.area = zeros(numel(obj.paneltype),1);
            obj.normal = zeros(numel(obj.paneltype),3);
            obj.center = zeros(numel(obj.paneltype),3);
            obj.Cp = {};
            obj.Cfe = {};
            obj.LHS = [];
            obj.RHS = [];
            for i = 1:numel(obj.paneltype)
                [obj.area(i,1),~ , obj.normal(i,:)] = obj.vertex(verts(connectivity(i,1),:),verts(connectivity(i,2),:),verts(connectivity(i,3),:));
                obj.center(i,:) = [mean(verts(obj.tri.ConnectivityList(i,:),1)),mean(verts(obj.tri.ConnectivityList(i,:),2)),mean(verts(obj.tri.ConnectivityList(i,:),3))];
            end

            
            %半裁メッシュの境界表面上のwakeを削除
            if obj.halfmesh == 1
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
                    if numel(obj.wakeline{iter}.edge) == 1
                        deleteiter = [deleteiter,iter];
                    end
                end
                obj.wakeline(deleteiter) = [];
            end
            %wakeのつくパネルIDを特定
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge)-1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{wakeNo}.edge(edgeNo),obj.wakeline{wakeNo}.edge(edgeNo+1));
                    phivert = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),2:3)-obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),2:3);
                    phiwake = atan2d(phivert(2),phivert(1));
                    if abs(phiwake)<90
                        obj.wakeline{wakeNo}.upperID(edgeNo) = attachpanel{1}(1);
                        obj.wakeline{wakeNo}.lowerID(edgeNo) = attachpanel{1}(2);
                    else
                        obj.wakeline{wakeNo}.upperID(edgeNo) = attachpanel{1}(2);
                        obj.wakeline{wakeNo}.lowerID(edgeNo) = attachpanel{1}(1);
                    end
                end
            end
            obj.checkMesh(sqrt(eps),"warning");
        end

        function obj = setVerts(obj,verts)
            %%%%%%%%%%%%%%%%%接点座標の変更%%%%%%%%%%%%%%%
            %三角形の接続関係を変えずに、接点の座標のみ変更する
            %変更時したときに代わる値は全てここで計算しなおす。
            %最適化を行うときに、設計変数の勾配を計算するとき等に使用する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %triangulationはreadonlyなので、新しく作り直す
            con = obj.tri.ConnectivityList;
            obj.tri = triangulation(con,verts);
            
            for i = 1:numel(obj.paneltype)
                [obj.area(i,1),~ , obj.normal(i,:)] = obj.vertex(verts(obj.tri(i,1),:),verts(obj.tri(i,2),:),verts(obj.tri(i,3),:));
                obj.center(i,:) = [mean(verts(obj.tri.ConnectivityList(i,:),1)),mean(verts(obj.tri.ConnectivityList(i,:),2)),mean(verts(obj.tri.ConnectivityList(i,:),3))];
            end
            obj.checkMesh(sqrt(eps),"warning");
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
                else
                    error("some panel areas are too small");
                end
            end
        end

        function obj = setPanelType(obj,ID,typename)
            %%%%%%%%%%%%%%%%paneltypeの指定関数%%%%%%%%%%%%%%%%
            %panelname : body 1 機体表面パネル
            %panelname : base 2 ベース面パネル
            %panelname : structure 3 構造パネル 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmpi(typename,"body")
                obj.paneltype(obj.surfID==ID,1) = 1;
            elseif strcmpi(typename,"base")
                obj.paneltype(obj.surfID==ID,1) = 2;
            elseif strcmpi(typename,"structure")
                obj.paneltype(obj.surfID==ID,1) = 3;    
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

        function obj = makeCluster(obj,nCluster,edgeAngleThreshold)
            %%%%%%%%%%%%%パネルクラスターの作成%%%%%%%%%%%%%%%%%%
            %各パネルの近隣パネルを指定の数集める
            %クラスター内のIDは統一（TODO：統一しないオプションをつけるか検討）
            %パネル法連立方程式を解いて得られるポテンシャルから機体表面に沿った微分を計算するための準備
            %nCluster:クラスターの目標数（達成できなければそこで打ち切り）
            %edgeAngleThreshold(deg):この角度以下であれば近隣パネルとして登録する（角度差が大きすぎるモノをはじくため）
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
                            edgeAngle = acos(dot(obj.normal(obj.cluster{i}(index),:),obj.normal(neighbor(j),:)));
                            if abs(edgeAngle)>edgeAngleThreshold*pi/180 || obj.surfID(obj.cluster{i}(index)) ~= obj.surfID(neighbor(j)) || obj.paneltype(neighbor(j)) ~= 1
                                deletelist = [deletelist,j];
                            end
                        end
                        neighbor(deletelist) = [];
                        diffcluster = setdiff(neighbor,obj.cluster{i});

                        if numel(obj.cluster{i})>nCluster
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

        function obj = makeEquation(obj,wakepanellength,nwake,n_divide)
            %%%%%%%%%%%%%パネル法連立方程式行列の作成%%%%%%%%%%%%
            %パネル法の根幹
            %表面微分行列も併せて作成している
            %wakepanellength:各wakeパネルの長さ
            %nwake:wakeパネルの数　つまりwakepanellength*nwakeがwakeの長さ、機体長の100倍以上あれば十分
            %n_divide:この関数は各パネル数×各パネル数サイズの行列を扱うため、莫大なメモリが必要となる。一度に計算する列をn_divide分割してメモリの消費量を抑えている。
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 3
                n_divide = 1;
            end
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            obj.wakePanelLength = wakepanellength;
            obj.nWake = nwake;
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
                    l = cross(m,obj.normal(i,:)');
                    Minv = [l,m,obj.normal(i,:)'];
                    lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                    bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                    Bmat=pinv(bb,sqrt(eps));
                    Vnmat = Minv(:,[1,2])*[1,0,0,0;0,1,0,0]*Bmat;
                    for iter = 1:3
                        obj.mu2v{iter}(i,obj.IndexPanel2Solver(obj.cluster{i})) = Vnmat(iter,:);
                    end
                end
            end

            %誘導抗力計算用の行列の作成
            
            %パネル方連立方程式行列の作成
            %機体パネル⇒機体パネルへの影響
            
            si = floor(nbPanel/n_divide).*(0:n_divide-1)+1;
            ei = [floor(nbPanel/n_divide).*(1:n_divide-1),nbPanel];
            obj.LHS = zeros(nbPanel);
            obj.RHS = zeros(nbPanel);

            for i= 1:n_divide
                [~,~,VortexAc,VortexBc] = obj.influenceMatrix(obj,[],si(i):ei(i));
                obj.LHS(:,si(i):ei(i)) = VortexAc; 
                obj.RHS(:,si(i):ei(i)) = VortexBc;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %wakeパネル⇒機体パネルへの影響
            %
            
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
                    interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                    interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                    [influence] = obj.wakeInfluenceMatrix(obj,wakeNo,edgeNo,1:nbPanel,wakepanellength,nwake);
                    obj.LHS(:,interpID(1)) = obj.LHS(:,interpID(1)) - influence;
                    obj.LHS(:,interpID(2)) = obj.LHS(:,interpID(2)) + influence;
                    obj.LLT.sp{wakeNo} =  [obj.LLT.sp{wakeNo},(s(edgeNo)+s(edgeNo+1))./2];
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
            obj.LLT.Qij = obj.Calc_Q(horzcat(obj.LLT.yinterp{:}),horzcat(obj.LLT.zinterp{:}),horzcat(obj.LLT.phiinterp{:}),horzcat(obj.LLT.spanel{:}),obj.halfmesh);
        end



        function obj = flowCondition(obj,flowNo,Mach,newtoniantype)
            %%%%%%%%%%%%%%%%主流の設定%%%%%%%%%%%%%%%%%%%%%%
            %UNLSIでは流れの状態の変化によって不連続的な変化を起こさないよう（設計変数の微分が連続になるよう）考慮されている。
            %特に超音速流解析では主流の状態によって各パネルの膨張と圧縮が切り替わるため、特別な考慮が必要となる。
            %UNLSIでは各マッハ数に応じてパネル角度-180deg~180degにて事前に圧力係数を計算し、実際の計算ではgriddedInterpolantを用いて値を求めている。
            %したがって、主流状態と補間関数を事前につくる必要がある。
            
            %flowNo:作成する流れのID
            %Mach:マッハ数
            %超音速解析(修正ニュートン流理論)の手法:"OldTangentCone"(デフォルト),"TangentConeEdwards","TangentWedge"
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 3
                newtoniantype = "OldTangentCone";
            end
            obj.flow{flowNo}.Mach = Mach;
            if Mach < 1
                obj.flow{flowNo}.flowtype = "Souce-Doublet 1st-order Panel";
            else
                obj.flow{flowNo}.flowtype = strcat(newtoniantype,"+Prandtl-Meyer");
            end
            kappa = 1.4;

            if obj.flow{flowNo}.Mach > 1
                delta = linspace(-pi,pi,500);
                Cpfc = zeros(size(delta));
                for j =1:size(delta,2)
                    if delta(j) >= 0
                        if strcmpi(newtoniantype,'TangentConeEdwards')
                            Mas = (0.87*Mach-0.554)*sin(delta(j))+0.53;
                            Cpfc(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(newtoniantype,'OldTangentCone')
                            Mas = 1.090909*Mach*sin(delta(j))+exp(-1.090909*Mach*sin(delta(j)));
                            Cpfc(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
                        elseif strcmpi(newtoniantype,'TangentWedge')
                            if delta(j)>45.585*pi/180
                                Cpfc(j)=((1.2*Mach*sin(delta(j))+exp(-0.6*Mach*sin(delta(j))))^2-1.0)/(0.6*Mach^2);
                                %R=1/Mach^2+(kappa+1)/2*delta(j)/sqrt(Mach^2-1);
                            elseif delta(j)<0.035
                                Cpfc(j)=kappa*Mach^2*delta(j)/sqrt(Mach^4-1.0);
                            else
                                b = -((Mach^2+2)/Mach^2)-kappa*sin(delta(j))^2;
                                c = (2*Mach^2+1)/Mach^4+((kappa+1)^2/4+(kappa-1)/Mach^2)*sin(delta(j))^2;
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
                                    Cpfc(j) = 4*(Mach^2*R-1)/((kappa+1)*Mach^2);
                                end
                            end
                        end
                        Pe_Pinf= Cpfc(j)*(kappa*Mach^2)/2+1;
                        Te_Tinf = Pe_Pinf^(1/(kappa/(kappa-1)));
                        Me = sqrt((2+(kappa+1)*Mach^2)/((kappa+1)*Te_Tinf)-2/(kappa+1));
                    else
                        Me = UNLSI.bisection(@(M2)UNLSI.pmsolve(Mach,M2,-delta(j),kappa),Mach,300);
                        Te_Tinf = (2+(kappa+1)*Mach^2)/(2+(kappa+1)*Me^2);
                        Pe_Pinf=(Te_Tinf)^(kappa/(kappa-1));
                        Cpfc(j) = 2/(kappa*Mach^2)*(Pe_Pinf-1);
                    end
                end
                obj.flow{flowNo}.pp = griddedInterpolant(delta,Cpfc,'spline');
            end
        end

        function obj = setCf(obj,flowNo,Re,Lch,k,LTratio,coefficient)
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
            %LTratio:層流の割合
            %coefficient:調整用の係数(デフォルト:1)

            %TODO:パネルIDによってCfの値を切り替えられるようにする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 6
                coefficient = 1;
            end
            %{

            %}
            if(obj.flow{flowNo}.Mach<0.9)
                Re_Cut = 38.21*((Lch/k)^1.053);
            else
                Re_Cut = 44.62*((Lch/k)^1.053).*((obj.flow{flowNo}.Mach)^1.16);
            end
            if Re>Re_Cut
                Re = Re_Cut;
            end
            Cf_L = 1.328./sqrt(Re);
            Cf_T = 0.455/((log10(Re).^(2.58))*((1+0.144*obj.flow{flowNo}.Mach*obj.flow{flowNo}.Mach)^0.65));
            obj.Cfe{flowNo}(obj.paneltype==1,1) =(LTratio*Cf_L + (1-LTratio)*Cf_T)*coefficient;
        end

        function obj = solveFlow(obj,flowNo,alpha,beta,omega)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %結果はobj.AERODATAに格納される。
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt 7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14:E(翼効率) 15:CFx 16:CFy 17:CFz 18:CMx 19:CMy 20:CMz 21:CMl 22:CMm 23:CMn 24:FOpt 
            %上記で求めていないものは0が代入される
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %omega:主流の回転角速度(deg/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(size(alpha) ~= size(beta))
                error("analysis points are not match");
            end
            obj.AERODATA{flowNo} = [];
            obj.Cp{flowNo} = [];
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
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

                if nargin < 5
                    Vinf = repmat((T*[1;0;0])',[nPanel,1]);
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (T*[1;0;0])'-(cross(omega(iterflow,:)./180.*pi,rvec(:)'));
                    end
                end
                Tvec(:,1) = obj.normal(:,2).* Vinf(:,3)-obj.normal(:,3).* Vinf(:,2);
                Tvec(:,2) = obj.normal(:,3).* Vinf(:,1)-obj.normal(:,1).* Vinf(:,3);
                Tvec(:,3) = obj.normal(:,1).* Vinf(:,2)-obj.normal(:,2).* Vinf(:,1);
                s(:,1) = Tvec(:,2).*obj.normal(:,3)-Tvec(:,3).*obj.normal(:,2);
                s(:,2) = Tvec(:,3).*obj.normal(:,1)-Tvec(:,1).*obj.normal(:,3);
                s(:,3) = Tvec(:,1).*obj.normal(:,2)-Tvec(:,2).*obj.normal(:,1);
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.normal(i,:)');
                            iter = iter+1;
                        end
                    end
                    RHV = obj.RHS*sigmas;
                    u =  -obj.LHS\RHV;
                    potential = u + sum(Vinf(obj.paneltype == 1,:).*obj.center(obj.paneltype == 1,:),2);
                    dv = zeros(nPanel,3);
                    for i = 1:3
                        dv(:,i) = obj.mu2v{i}*(potential);
                    end
                
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = (1-sum(dv(obj.paneltype==1,:).^2,2))./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-0.139-0.419.*(obj.flow{flowNo}.Mach-0.161).^2);
                    uinterp = [];
                    for i = 1:numel(obj.wakeline)
                        uinterp = [uinterp,interp1(obj.LLT.sp{i},obj.LLT.calcMu{i}*u,obj.LLT.sinterp{i},'linear','extrap')];
                    end
                    Vind = obj.LLT.Qij*uinterp';
                    %{
                    figure(3);clf;hold on;
                    plot(horzcat(obj.LLT.yinterp{:}),Vind');
                    plot(horzcat(obj.LLT.yinterp{:}),uinterp);
                    ylim([-0.5,3]);
                    %}
                    CLt = (2.*horzcat(obj.LLT.spanel{:}).*cos(horzcat(obj.LLT.phiinterp{:}))*uinterp')/(0.5*obj.SREF)/(1-obj.flow{flowNo}.Mach^2);
                    CDt = ((uinterp.*horzcat(obj.LLT.spanel{:}))*Vind)/(0.5*obj.SREF)/norm(1-obj.flow{flowNo}.Mach^2)^3;
                else
                    %超音速
                    delta = zeros(nbPanel,1);
                    iter = 1;
                    %各パネルが主流となす角度を求める
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            delta(iter,1) = acos(dot(obj.normal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
                            iter = iter+1;
                        end
                    end
                    %用意された応答曲面をもちいてパネルの角度からCpを求める
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = obj.flow{flowNo}.pp(delta);%Cp
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-obj.flow{flowNo}.Mach.^(-2)+0.57.*obj.flow{flowNo}.Mach.^(-4));
                end
                %Cp⇒力への変換
                dCA_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,3)).*obj.area./obj.SREF;
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
                
                obj.AERODATA{flowNo}(iterflow,:) = [beta(iterflow),obj.flow{flowNo}.Mach,alpha(iterflow),0,CL,CLt,CDo,CDi,CDtot,CDt,CDtott,CY,CLt/CDtott,CLt^2/pi/AR/CDt,CAp+CAf,CYp+CYf,CNp+CNf,CMX,CMY,CMZ,0,0,0,0];
                
            end
            %disp([CL,CDo,CDi,CDtot,CMY]);
        end

        function [u,R] = solvePertPotential(obj,flowNo,alpha,beta,omega)
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
                if nargin < 5
                    Vinf = repmat((T*[1;0;0])',[nPanel,1]);
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (T*[1;0;0])'-(cross(omega(iterflow,:)./180.*pi,rvec(:)))';
                    end
                end
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.normal(i,:)');
                            iter = iter+1;
                        end
                    end
                    RHV = obj.RHS*sigmas;
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
                            delta(iter,1) = acos(dot(obj.normal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
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
        
        function [AERODATA,Cp,Cfe,R,obj] = solveFlowForAdjoint(obj,u,flowNo,alpha,beta,omega)
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
            obj.Cp{flowNo} = [];
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
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
                if nargin < 6
                    Vinf = repmat((T*[1;0;0])',[nPanel,1]);
                else
                    Vinf = zeros(nPanel,3);
                    for i = 1:nPanel
                       rvec = obj.center(i,:)'-obj.XYZREF(:);
                       Vinf(i,:) = (T*[1;0;0])'-(cross(omega(iterflow,:)./180.*pi,rvec(:)))';
                    end
                end
                Tvec(:,1) = obj.normal(:,2).* Vinf(:,3)-obj.normal(:,3).* Vinf(:,2);
                Tvec(:,2) = obj.normal(:,3).* Vinf(:,1)-obj.normal(:,1).* Vinf(:,3);
                Tvec(:,3) = obj.normal(:,1).* Vinf(:,2)-obj.normal(:,2).* Vinf(:,1);
                s(:,1) = Tvec(:,2).*obj.normal(:,3)-Tvec(:,3).*obj.normal(:,2);
                s(:,2) = Tvec(:,3).*obj.normal(:,1)-Tvec(:,1).*obj.normal(:,3);
                s(:,3) = Tvec(:,1).*obj.normal(:,2)-Tvec(:,2).*obj.normal(:,1);
                if obj.flow{flowNo}.Mach < 1
                    %亜音速
                    sigmas = zeros(nbPanel,1);
                    iter = 1;
                    for i = 1:nPanel
                        if obj.paneltype(i) == 1
                            sigmas(iter,1) = dot(Vinf(i,:)',obj.normal(i,:)');
                            iter = iter+1;
                        end
                    end
                    usolve = u(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1);
                    potential = usolve + sum(Vinf(obj.paneltype == 1,:).*obj.center(obj.paneltype == 1,:),2);
                    dv = zeros(nPanel,3);
                    for i = 1:3
                        dv(:,i) = obj.mu2v{i}*(potential);
                    end
                
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = (1-sum(dv(obj.paneltype==1,:).^2,2))./(1-obj.flow{flowNo}.Mach^2);
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-0.139-0.419.*(obj.flow{flowNo}.Mach-0.161).^2);
                    uinterp = [];
                    for i = 1:numel(obj.wakeline)
                        uinterp = [uinterp,interp1(obj.LLT.sp{i},obj.LLT.calcMu{i}*usolve,obj.LLT.sinterp{i},'linear','extrap')];
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
                            delta(iter,1) = acos(dot(obj.normal(i,:)',Vinf(i,:)')/norm(Vinf(i,:)))-pi/2;%パネル角度
                            iter = iter+1;
                        end
                    end
                    usolve = u(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1);
                    obj.Cp{flowNo}(obj.paneltype==1,iterflow) = usolve;%Cp
                    obj.Cp{flowNo}(obj.paneltype==2,iterflow) = (-obj.flow{flowNo}.Mach.^(-2)+0.57.*obj.flow{flowNo}.Mach.^(-4));
                end
                %Cp⇒力への変換
                dCA_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-obj.Cp{flowNo}(:,iterflow).*obj.normal(:,3)).*obj.area./obj.SREF;
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
                obj.AERODATA{flowNo}(iterflow,:) = [beta(iterflow),obj.flow{flowNo}.Mach,alpha(iterflow),0,CL,CLt,CDo,CDi,CDtot,CDt,CDtott,CY,CLt/CDtott,CLt^2/pi/AR/CDt,CAp+CAf,CYp+CYf,CNp+CNf,CMX,CMY,CMZ,0,0,0,0];
                AERODATA = obj.AERODATA;
                Cp = obj.Cp;
                Cfe = obj.Cfe;
                if obj.flow{flowNo}.Mach < 1
                    RHV = obj.RHS*sigmas;
                    R(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1) = obj.LHS*usolve+RHV;
                else
                    R(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1) = usolve-obj.flow{flowNo}.pp(delta);
                end
            end
            %disp([CL,CDo,CDi,CDtot,CMY]);
        end
        
        function dynCoef = calcDynCoef(obj,flowNo,alpha,beta,difference)
            %%%動微係数の計算%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%TO DO alpha とbetaも
            %結果はDyncoefに格納される
            %[1:Clp, 2:Cmp, 3:Cnp;
            % 1:Clq, 2:Cmq, 3:Cnq;
            % 1:Clr, 2:Cmr, 3,Cnr]
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %difference:有限差分の方法 "forward"(デフォルト)-前進差分 "central"-中心差分
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist("difference","var"); difference = "forward"; end
            
            obj = obj.solveFlow(flowNo,alpha,beta);
            %有限差分による計算
            ind = [15:20];%AERODATA = [18:CMX, 19:CMY, 20: CMZ]
            nondim = [obj.BREF/(2*1),obj.CREF/(2*1),obj.BREF/(2*1)];
            dynCoef = zeros(5,6);
            if strcmp(difference,"forward")
                dw = sqrt(eps);
                delta = obj.solveFlow(flowNo,alpha,beta+dw);
                tmp = (delta.AERODATA{flowNo} - obj.AERODATA{flowNo})./(dw/180*pi);
                dynCoef(1,:) = tmp(ind); %beta
                delta = obj.solveFlow(flowNo,alpha+dw,beta);
                tmp = (delta.AERODATA{flowNo} - obj.AERODATA{flowNo})./(dw/180*pi);
                dynCoef(2,:) = tmp(ind); % alpha
                for i = 1:3%p,q,rについて差分をとる
                    omega = zeros(1,3);
                    omega(i) = dw;
                    delta = obj.solveFlow(flowNo,alpha,beta,omega);
                    tmp = (delta.AERODATA{flowNo} - obj.AERODATA{flowNo})./((dw/180*pi)*nondim(i));
                    dynCoef(2+i,:) = tmp(ind);
                end
            elseif strcmp(difference,"central")
                dw = eps^(1/3);
                delta1 = obj.solveFlow(flowNo,alpha,beta+dw);
                delta2 = obj.solveFlow(flowNo,alpha,beta-dw);
                tmp = (delta1.AERODATA{flowNo} - delta2.AERODATA{flowNo})./(2*(dw/180*pi));
                dynCoef(1,:) = tmp(ind);
                delta1 = obj.solveFlow(flowNo,alpha+dw,beta);
                delta2 = obj.solveFlow(flowNo,alpha-dw,beta);
                tmp = (delta1.AERODATA{flowNo} - delta2.AERODATA{flowNo})./(2*(dw/180*pi));
                dynCoef(2,:) = tmp(ind);
                for i = 1:3
                    omega = zeros(1,3);
                    omega(i) = dw;
                    delta1 = obj.solveFlow(flowNo,alpha,beta,omega);
                    delta2 = obj.solveFlow(flowNo,alpha,beta,-omega);
                    tmp = (delta1.AERODATA{flowNo} - delta2.AERODATA{flowNo})./(2*(dw/180*pi)*nondim(i));
                    dynCoef(2+i,:) = tmp(ind);
                end
            else 
                error("Supported difference methods are ""foraward"" and ""central"".")
            end
            axisRot = [ 0, 1, 0,-1, 0,-1;
                       -1, 0,-1, 0, 1, 0;
                        0,-1, 0, 1, 0, 1;
                        0, 0,-1, 0, 1, 0;
                        0,-1, 0, 1, 0, 1];
            dynCoef = (axisRot.*dynCoef)';

        end
    end

    methods(Static)
        
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
            normal = obj.normal(obj.paneltype==1,:);

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
    
                %1回目
                a.X = POI.X-N1.X;
                a.Y = POI.Y-N1.Y;
                a.Z = POI.Z-N1.Z;
                b.X = POI.X-N2.X;
                b.Y = POI.Y-N2.Y;
                b.Z = POI.Z-N2.Z;
                s.X = N2.X-N1.X;
                s.Y = N2.Y-N1.Y;
                s.Z = N2.Z-N1.Z;
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                VortexAr = phiV;
                VortexBr = srcV;
                %2回目
                a.X = POI.X-N2.X;
                a.Y = POI.Y-N2.Y;
                a.Z = POI.Z-N2.Z;
                b.X = POI.X-N3.X;
                b.Y = POI.Y-N3.Y;
                b.Z = POI.Z-N3.Z;
                s.X = N3.X-N2.X;
                s.Y = N3.Y-N2.Y;
                s.Z = N3.Z-N2.Z;
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                VortexAr = VortexAr+phiV;
                VortexBr = VortexBr+srcV;
                %3回目
                a.X = POI.X-N3.X;
                a.Y = POI.Y-N3.Y;
                a.Z = POI.Z-N3.Z;
                b.X = POI.X-N1.X;
                b.Y = POI.Y-N1.Y;
                b.Z = POI.Z-N1.Z;
                s.X = N1.X-N3.X;
                s.Y = N1.Y-N3.Y;
                s.Z = N1.Z-N3.Z; 
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                VortexAr = VortexAr+phiV;
                VortexBr = VortexBr+srcV;
    
                %半裁の考慮
                if obj.halfmesh == 1
                    POI.Y = -POI.Y;
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    %1回目
                    a.X = POI.X-N1.X;
                    a.Y = POI.Y-N1.Y;
                    a.Z = POI.Z-N1.Z;
                    b.X = POI.X-N2.X;
                    b.Y = POI.Y-N2.Y;
                    b.Z = POI.Z-N2.Z;
                    s.X = N2.X-N1.X;
                    s.Y = N2.Y-N1.Y;
                    s.Z = N2.Z-N1.Z;
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                    VortexAr = VortexAr+phiV;
                    VortexBr = VortexBr+srcV;
                    %2回目
                    a.X = POI.X-N2.X;
                    a.Y = POI.Y-N2.Y;
                    a.Z = POI.Z-N2.Z;
                    b.X = POI.X-N3.X;
                    b.Y = POI.Y-N3.Y;
                    b.Z = POI.Z-N3.Z;
                    s.X = N3.X-N2.X;
                    s.Y = N3.Y-N2.Y;
                    s.Z = N3.Z-N2.Z;
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                    VortexAr = VortexAr+phiV;
                    VortexBr = VortexBr+srcV;
                    %3回目
                    a.X = POI.X-N3.X;
                    a.Y = POI.Y-N3.Y;
                    a.Z = POI.Z-N3.Z;
                    b.X = POI.X-N1.X;
                    b.Y = POI.Y-N1.Y;
                    b.Z = POI.Z-N1.Z;
                    s.X = N1.X-N3.X;
                    s.Y = N1.Y-N3.Y;
                    s.Z = N1.Z-N3.Z; 
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
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
        
                    %1回目
                    a.X = POI.X-N1.X;
                    a.Y = POI.Y-N1.Y;
                    a.Z = POI.Z-N1.Z;
                    b.X = POI.X-N2.X;
                    b.Y = POI.Y-N2.Y;
                    b.Z = POI.Z-N2.Z;
                    s.X = N2.X-N1.X;
                    s.Y = N2.Y-N1.Y;
                    s.Z = N2.Z-N1.Z;
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                    VortexAcb = phiV;
                    VortexBcb = srcV;
                    %2回目
                    a.X = POI.X-N2.X;
                    a.Y = POI.Y-N2.Y;
                    a.Z = POI.Z-N2.Z;
                    b.X = POI.X-N3.X;
                    b.Y = POI.Y-N3.Y;
                    b.Z = POI.Z-N3.Z;
                    s.X = N3.X-N2.X;
                    s.Y = N3.Y-N2.Y;
                    s.Z = N3.Z-N2.Z;
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                    VortexAcb = VortexAcb+phiV;
                    VortexBcb = VortexBcb+srcV;
                    %3回目
                    a.X = POI.X-N3.X;
                    a.Y = POI.Y-N3.Y;
                    a.Z = POI.Z-N3.Z;
                    b.X = POI.X-N1.X;
                    b.Y = POI.Y-N1.Y;
                    b.Z = POI.Z-N1.Z;
                    s.X = N1.X-N3.X;
                    s.Y = N1.Y-N3.Y;
                    s.Z = N1.Z-N3.Z; 
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                    VortexAcb = VortexAcb+phiV;
                    VortexBcb = VortexBcb+srcV;
        
                    %半裁の考慮
                    if obj.halfmesh == 1
                        POI.Y = -POI.Y;
                        pjk.X = POI.X-c.X;
                        pjk.Y = POI.Y-c.Y;
                        pjk.Z = POI.Z-c.Z;
                        %1回目
                        a.X = POI.X-N1.X;
                        a.Y = POI.Y-N1.Y;
                        a.Z = POI.Z-N1.Z;
                        b.X = POI.X-N2.X;
                        b.Y = POI.Y-N2.Y;
                        b.Z = POI.Z-N2.Z;
                        s.X = N2.X-N1.X;
                        s.Y = N2.Y-N1.Y;
                        s.Z = N2.Z-N1.Z;
                        phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                        srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                        VortexAcb = VortexAcb+phiV;
                        VortexBcb = VortexBcb+srcV;
                        %2回目
                        a.X = POI.X-N2.X;
                        a.Y = POI.Y-N2.Y;
                        a.Z = POI.Z-N2.Z;
                        b.X = POI.X-N3.X;
                        b.Y = POI.Y-N3.Y;
                        b.Z = POI.Z-N3.Z;
                        s.X = N3.X-N2.X;
                        s.Y = N3.Y-N2.Y;
                        s.Z = N3.Z-N2.Z;
                        phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                        srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
                        VortexAcb = VortexAcb+phiV;
                        VortexBcb = VortexBcb+srcV;
                        %3回目
                        a.X = POI.X-N3.X;
                        a.Y = POI.Y-N3.Y;
                        a.Z = POI.Z-N3.Z;
                        b.X = POI.X-N1.X;
                        b.Y = POI.Y-N1.Y;
                        b.Z = POI.Z-N1.Z;
                        s.X = N1.X-N3.X;
                        s.Y = N1.Y-N3.Y;
                        s.Z = N1.Z-N3.Z; 
                        phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                        srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
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

        function [VortexA] = wakeInfluenceMatrix(obj,wakeNo,edgeNo,rowIndex,xwake,nwake)
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
            for i = 1:nwake
                wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),:)+[xwake*(i-1),0,0];
                wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),:)+[xwake*(i-1),0,0];
                wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo+1),:)+[xwake*i,0,0];
                wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo),:)+[xwake*i,0,0];
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
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
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
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
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
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
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
                phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                VortexA = VortexA+phiV;

                %半裁
                if obj.halfmesh == 1
                    POI.Y = -POI.Y;
                    pjk.X = POI.X-c.X;
                    pjk.Y = POI.Y-c.Y;
                    pjk.Z = POI.Z-c.Z;
                    
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
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    %srcV = (obj.matrix_dot(n,obj.matrix_cross(s,a))).*(1./obj.matrix_norm(s).*log(((obj.matrix_norm(a)+obj.matrix_norm(b)+obj.matrix_norm(s))./(obj.matrix_norm(a)+obj.matrix_norm(b)-obj.matrix_norm(s)))))-(obj.matrix_dot(pjk,n)).*phiV;
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
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
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
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
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
                    phiV = atan2(obj.matrix_dot(s,obj.getUnitVector(c,n12)).*obj.matrix_dot(pjk,n).*(obj.matrix_norm(b).*(obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-obj.matrix_norm(a).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))),((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s)))).*((obj.matrix_dot(a,obj.matrix_cross(obj.matrix_cross(obj.getUnitVector(c,n12),n),obj.matrix_cross(a,s))))-(obj.matrix_dot(n,obj.matrix_cross(s,a))).*obj.matrix_dot(s,obj.getUnitVector(c,n12)))+obj.matrix_dot(pjk,n).^2.*obj.matrix_norm(a).*obj.matrix_norm(b).*(obj.matrix_dot(s,obj.getUnitVector(c,n12))).^2));
                    VortexA = VortexA+phiV;
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
    end

end