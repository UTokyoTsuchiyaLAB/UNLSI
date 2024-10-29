classdef UNLSI

    %%%%%%%%%%%%UNLSI%%%%%%%%%%%%%%%%%%%
    %The MIT License Copyright (c) 2023 Naoto Morita
    %以下に定める条件に従い、本ソフトウェアおよび関連文書のファイル（以下「ソフトウェア」）の複製を取得するすべての人に対し、ソフトウェアを無制限に扱うことを無償で許可します。
    %これには、ソフトウェアの複製を使用、複写、変更、結合、掲載、頒布、サブライセンス、および/または販売する権利、およびソフトウェアを提供する相手に同じことを許可する権利も無制限に含まれます。
    %上記の著作権表示および本許諾表示を、ソフトウェアのすべての複製または重要な部分に記載するものとします。
    %ソフトウェアは「現状のまま」で、明示であるか暗黙であるかを問わず、何らの保証もなく提供されます。
    %ここでいう保証とは、商品性、特定の目的への適合性、および権利非侵害についての保証も含みますが、それに限定されるものではありません。 
    %作者または著作権者は、契約行為、不法行為、またはそれ以外であろうと、ソフトウェアに起因または関連し、あるいはソフトウェアの使用またはその他の扱いによって生じる一切の請求、損害、その他の義務について何らの責任も負わないものとします。
    
    %また、本ソフトウェアに付帯するソフトウェアとしてOpenVSPおよびdistanceVertex2Mesh.m, SurfaceIntersection.mを用いています。
    %上記ソフトウェアのライセンス表記を以下に示します。
    %openVSP
    %Copyright (c) 2012 United States Government as represented by the Administrator for The National Aeronautics and Space Administration. All Rights Reserved.
    %
    %DISTANCEVERTEX2MESH - calculate the distance between vertices and a mesh
    %Author: Christopher Haccius
    %Telecommunications Lab, Saarland University, Germany
    %email: haccius@nt.uni-saarland.de
    %March 2015; Last revision: 26-March-2015
    %SurfaceIntersection
    %Copyright (c) 2014, Jaroslaw Tuszynski
    %All rights reserved.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        settingUNLSI %パラメータセッティング
        tri %MATLAB triangulationクラス
        id %各パネルのID
        surfID %各パネルのタグ番号
        wakeline %wakeパネルの設定
        wakelineID %生データ
        halfmesh %半裁かどうか(半裁:1)
        SREF %基準面積
        CREF %縦方向基準長
        BREF %横方向基準長
        XYZREF %回転中心
        rotCenter %各パネルのRotationの中心
        paneltype %各パネルのタイプ 1:body 2:base 3:structure
        cpcalctype %Cpを計算する方法の選択 %1:1-V^2 2: -2u 3: -2u-v^2-w^2
        IndexPanel2Solver %パネルのインデックス⇒ソルバー上でのインデックス
        approxMat %パネル法行列の近似行列を作成するためのセル
        approxFemMat %パネル法行列の近似行列を作成するためのセル
        approximated %パネル法行列が近似されたものかどうか
        femapproximated %パネル法行列が近似されたものかどうか
        orgNormal %各パネルの法線ベクトル
        modNormal %舵角等によって変更した法線ベクトル
        center %各パネルの中心
        area %各パネルの面積
        angularVelocity %各パネルの回転角速度
        flowNoTable %flowNoが整理されたMatrix
        flow %各flowconditionが格納されたセル配列
        prop %各プロペラ方程式情報が格納されたセル配列
        cluster %パネルクラスターの情報
        LHS %パネル法連立方程式の左辺行列
        wakeLHS %パネル法のwake影響行列
        RHS %パネル法連立方程式の右辺行列
        mu2v %ポテンシャル⇒機体表面速度への変換行列
        verts2centerMat %節点での量⇛パネルセンターでの量への変換行列
        Cp %圧力係数
        Cfe %表面摩擦係数
        CpLimit%圧力係数の最大最小値
        deflAngle %舵角
        deflGroup %舵面グループ
        AERODATA %空力解析結果の格納
        DYNCOEF %同安定微係数結果の格納
        LLT %トレフツ面解析用の構造体
        ppCoef %空力補完局面(griddedInterp)
        ppDyn %動安定微係数補完局面(griddedInterp)
        femutils
        femLHS
        femLHS_L
        femMass
        femMass_L
        femDamp
        femtri
        femID
        femarea
        femNormal
        femcenter
        femThn
        femE
        femRho
        femVisc
        fem2aeroMat
        aero2femMat
        femverts2centerMat %femメッシュでの節点量⇛パネルセンター量への変換行列
        df_ddelta %FEMにおける力の変位に関する微分
        femEigenVec %Femの結果を固有値解析した固有ベクトル
        femEigenVal %固有値
        selfLoadFlag
    end

    methods(Access = public)
        % UNLSIクラスのコンストラクタ
        % verts: 頂点座標の行列
        % connectivity: 頂点の接続情報の行列
        % surfID: サーフェスIDのベクトル
        % wakelineID: ウェイクラインIDのベクトル
        % halfmesh: 半裁フラグ (0: 半裁でない, 1: 半裁)
        % vertsMergeTol: 頂点のマージ許容誤差 (省略可)
        function obj = UNLSI(verts,connectivity,surfID,wakelineID,halfmesh,vertsMergeTol)
            warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
            if ~exist('vertsMergeTol', 'var')
                obj.settingUNLSI.vertsTol = 0.001;
            else
                obj.settingUNLSI.vertsTol = vertsMergeTol;
            end

            % halfmeshの判定
            if all(verts(:,2) >=0) ||  all(verts(:,2) <=0)
                if exist('halfmesh', 'var') 
                    if halfmesh == 0
                        warning("Signatures of verts' y cordinates are all same, but the half mesh settingUNLSI is 0 (not halfmesh)");
                    end
                else
                    halfmesh = 1;
                end
            else
               if exist('halfmesh', 'var')
                    if halfmesh == 1
                        warning("Signatures of verts' y cordinates are NOT all same, but the half mesh settingUNLSI is 1 (halfmesh)");
                    end
                else
                    halfmesh = 0;
               end
            end

            if isempty(halfmesh)
                if all(verts(:,2) >=0) ||  all(verts(:,2) <=0)
                    halfmesh = 1;
                else
                    halfmesh = 0;
                end
            end

            %%%%%%%%%%%settingUNLSIの初期設定%%%%%%%%%%%%%
            obj.settingUNLSI.checkMeshTol = sqrt(eps); % メッシュの精度をチェックするための許容誤差
            obj.settingUNLSI.LLTnInterp = 20; % LLT法の補間点の数
            obj.settingUNLSI.nCluster = 50; % パネルクラスターの目標数
            obj.settingUNLSI.edgeAngleThreshold = 50; % 近隣パネルとして登録するための角度の閾値
            obj.settingUNLSI.nCalcDivide = 5; % makeEquationの計算を分割するための分割数
            %obj.settingUNLSI.angularVelocity = []; % 正規化された主流の角速度
            obj.settingUNLSI.propCalcFlag = 1; % プロペラの計算フラグ
            obj.settingUNLSI.deflDerivFlag = 1; % 舵角の微分フラグ
            obj.settingUNLSI.propWakeLength = 3; % プロペラwakeの長さ
            obj.settingUNLSI.nPropWake = 101; % propWake円周上の点の数
            obj.settingUNLSI.kCf = 1.015*(10^-5); % 摩擦係数に関するk
            obj.settingUNLSI.laminarRatio = 0; % 層流の割合
            obj.settingUNLSI.coefCf = 1; % Cfの係数(単純な倍数)
            obj.settingUNLSI.newtoniantype = "OldTangentCone"; % 圧縮におけるニュートン流のタイプ
            obj.settingUNLSI.resultSearchMethod = "and"; % 結果の検索方法
            obj.settingUNLSI.Vinf = 15; % 流速
            obj.settingUNLSI.rho = 1.225; % 空気密度
            obj.settingUNLSI.kappa = 1.4; % 比熱比
            obj.settingUNLSI.g0 = 9.8; % 重力加速度
            obj.settingUNLSI.nGriddedInterp = 90; % griddedInterpの補間点の数
            obj.settingUNLSI.ode23AbsTol = 1e-6; % ode23の絶対許容誤差
            obj.settingUNLSI.ode23RelTol = 1e-3; % ode23の相対許容誤差
            obj.settingUNLSI.lsqminnormTol = 1e-12;% 最小二乗法における正規方程式の解を求める際の許容誤差を設定します。
            obj.settingUNLSI.xWakeAttach = 0.1;
            obj.settingUNLSI.defaultWakeLength = 100; % デフォルトの各wakeパネルの長さ（機体基準長基準）
            obj.settingUNLSI.nWakeMax = 20;%marchwakeにおけるwakeの最大個数
            obj.settingUNLSI.deltaVortexCore = 0.1;
            obj.settingUNLSI.pnThreshold = sqrt(eps);
            obj.settingUNLSI.defaultCpLimit = [-1000,1000];
            obj.settingUNLSI.CpSlope = 5;
            obj.settingUNLSI.selfLoadFlag = 0;
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


  
        function obj = setUNLSISettings(obj,varargin)
            % setUNLSISettings(obj,varargin)
            %   設定を変更する関数です。
            %
            % 入力:
            %   - obj: UNLSI オブジェクト
            %   - varargin: 変数名と値のペアで入力
            %
            % 出力:
            %   - obj: 更新された UNLSI オブジェクト
            %
            % 説明:
            %   この関数は、UNLSI オブジェクトの設定を変更します。
            %   varargin には、変数名と値のペアを指定します。
            %   obj.settingUNLSI の構造体の内部を変更します。
            %   構造体の内部変数の名前と一致していない場合はエラーが発生します。
            %
            % 例:
            %   obj = setUNLSISettings(obj, 'kCf', 0.5, 'laminarRatio', 0.2);
            %   上記の例では、obj の kCf と laminarRatio の値が変更されます。
            %   変更後の obj が返されます。
            %
            % 注意:
            %   'kCf'、'laminarRatio'、'coefCf'、'newtoniantype'、'kappa' のいずれかの変数が変更された場合、
            %   obj の flowNoTable の各行に対して flowCondition を呼び出します。
            %
            %   例:
            %   obj = setUNLSISettings(obj, 'kCf', 0.5);
            %   上記の例では、obj の kCf の値が変更されます。
            %   obj の flowNoTable の各行に対して flowCondition を呼び出します。
            %
            %   もし、変数名が一致しない場合はエラーが発生します。
            %
            %   例:
            %   obj = setUNLSISettings(obj, 'invalidField', 0.5);
            %   上記の例では、変数名が一致しないためエラーが発生します。
            %
            %   注意事項を守り、正しい変数名と値のペアを指定してください。
            %
            if mod(numel(varargin),1) == 1
                error("input style not match");
            end
            for iter = 1:numel(varargin)/2
                if isfield(obj.settingUNLSI,varargin{2*iter-1})
                    obj.settingUNLSI = setfield(obj.settingUNLSI,varargin{2*iter-1},varargin{2*iter});
                    if strcmp(varargin{2*iter-1}, "kCf") || strcmp(varargin{2*iter-1},"laminarRatio") || strcmp(varargin{2*iter-1},"coefCf") || strcmp(varargin{2*iter-1}, "newtoniantype") || strcmp(varargin{2*iter-1}, "kappa")
                        for i = 1:size(obj.flowNoTable,1)
                            obj = obj.flowCondition(i,obj.flowNoTable(i,1),obj.flowNoTable(i,2));
                        end
                    end
                else
                    error("Field name is not match");
                end
            end
        end

        function obj = setREFS(obj,SREF,BREF,CREF,XYZREF)
            % setREFS - 参照面積と基準長を設定する関数
            %
            %   obj = setREFS(obj,SREF,BREF,CREF) は、参照面積(SREF)、横方向基準長(BREF)、
            %   縦方向基準長(CREF)を設定するメソッドです。
            %
            %   入力:
            %   obj - UNLSIオブジェクト
            %   SREF - 基準面積
            %   BREF - 横方向基準長(例: 翼幅)
            %   CREF - 縦方向基準長(例: MAC, 機体全長)
            %
            %   出力:
            %   obj - 更新されたUNLSIオブジェクト
            %
            %   例:
            %   obj = setREFS(obj, 10, 5, 20)
            %
            %   参考文献:
            %   [1] UNLSI Documentation, 2021
            %
            %   この関数は、UNLSIクラスのプロパティであるSREF、BREF、CREFを設定します。
            %   これらの値は、後続の計算や解析に使用されます。
            %
            %   この関数は、UNLSIオブジェクトを返します。更新されたオブジェクトは、
            %   後続の操作で使用することができます。
            %
            %   この関数は、UNLSIクラスのメソッドです。UNLSIクラスのインスタンスを作成し、
            %   そのインスタンスに対してこのメソッドを呼び出すことができます。
            %
            %   この関数は、UNLSIパッケージの一部です。


            obj.SREF = SREF;
            obj.CREF = CREF;
            obj.BREF = BREF;
            obj.XYZREF = XYZREF(:)';
        end

        function obj = setCpLimit(obj,limitVal,ID)
            if nargin == 2
                obj.CpLimit = repmat(limitVal,[size(obj.paneltype,1),1]);
            else
                if size(limitVal,1)==numel(ID)
                    for i = 1:numel(ID)
                        obj.CpLimit(obj.surfID == ID(i),:) = repmat(limitVal,[sum(obj.surfID == ID(i)),1]);
                    end
                elseif size(limitVal,1) == 1
                    obj.CpLimit(any(obj.surfID == ID,2),:) = repmat(limitVal,[sum(any(obj.surfID == ID,2)),1]);
                else
                    error("Invalid Input");
                end
            end
        end


        function obj = setRotation(obj,ID,rotCenter,angularVel)
            %%%%%%%%%%%Rotation Center settingUNLSI%%%%%%%%%%%%%
            %各パネルの回転を設定する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if numel(ID) == size(rotCenter,1) && numel(ID) == size(angularVel,1)
                for iter = 1:numel(ID)
                    obj.rotCenter(obj.surfID == ID(iter),1:3) = rotCenter(iter,1:3);
                    obj.angularVelocity(obj.surfID == ID(iter),1:3) = angularVel(iter,1:3);
                end
            elseif numel(ID) == size(rotCenter,1) 
                for iter = 1:numel(ID)
                    obj.rotCenter(obj.surfID == ID(iter),1:3) = rotCenter(iter,1:3);
                    obj.angularVelocity(obj.surfID == ID(iter),1:3) = repmat(angularVel(:)',[sum(obj.surfID == ID(iter)),1]);
                end
            else
                for iter = 1:numel(ID)
                    obj.rotCenter(obj.surfID == ID(iter),1:3) = repmat(rotCenter(:)',[sum(obj.surfID == ID(iter)),1]);
                    obj.angularVelocity(obj.surfID == ID(iter),1:3) = repmat(angularVel(:)',[sum(obj.surfID == ID(iter)),1]);
                end
            end
        end

        function obj = setDeflGroup(obj,groupNo,groupName,groupID,deflGain)
            % setDeflGroup - デフレクトグループの設定を行う関数
            %   obj = setDeflGroup(obj, groupNo, groupName, groupID, deflGain)は、指定されたデフレクトグループの設定を行います。
            %   obj: UNLSIオブジェクト
            %   groupNo: デフレクトグループの番号
            %   groupName: デフレクトグループの名前
            %   groupID: デフレクトグループのID
            %   deflGain: デフレクトゲイン
            %   返り値:
            %   obj: 更新されたUNLSIオブジェクト
            %
            %   例:
            %   obj = setDeflGroup(obj, 1, 'Group1', 123, 0.5)
            %
            %   注意:
            %   - groupNoは1から始まる番号である必要があります。
            %   - groupNameはデフレクトグループの名前です。
            %   - groupIDはデフレクトグループのIDです。
            %   - deflGainはデフレクトゲインです。
            %   - obj.deflGroup{groupNo}にデフレクトグループの情報が格納されます。
            %   - obj.deflGroup{groupNo}.nameにデフレクトグループの名前が格納されます。
            %   - obj.deflGroup{groupNo}.IDにデフレクトグループのIDが格納されます。
            %   - obj.deflGroup{groupNo}.gainにデフレクトゲインが格納されます。
            obj.deflGroup{groupNo}.name = groupName;
            obj.deflGroup{groupNo}.ID = groupID(:)';
            obj.deflGroup{groupNo}.gain = deflGain(:)';
        end

        function obj = setDeflAngle(obj,ID,rotAxis,angle)
            %%%%%%%%%%%Rotation Center settingUNLSI%%%%%%%%%%%%%
            %疑似的な舵角を設定する。
            %そのID上のパネルの法線ベクトルをロドリゲスの回転ベクトルによって回転する
            % 誤差については要検討。
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rotAxis = rotAxis./vecnorm(rotAxis,2,2);
            for iter = 1:numel(ID)
                i = find(obj.deflAngle(:,1)==ID(iter));
                if isempty(i)
                    error("Invalid ID");
                end
                dcm = obj.rod2dcm(rotAxis(iter,:),angle(iter));
                obj.deflAngle(i,2:end) = [rotAxis(iter,:),angle(iter),dcm(:)'];
            end
            indices = arrayfun(@(x) find(obj.deflAngle(:, 1) == x), obj.surfID(:, 1), 'UniformOutput', false);      
            deflMatrices = cell2mat(cellfun(@(idx) reshape(obj.deflAngle(idx, 6:14), [3, 3]), indices, 'UniformOutput', false));
            obj.modNormal = reshape(cell2mat(arrayfun(@(i) (deflMatrices((i-1)*3+1:i*3, :) * obj.orgNormal(i, :)')', 1:numel(obj.paneltype), 'UniformOutput', false)),3,[])';
        end


        function plotGeometry(obj,figureNo,triColor,colorlim,method,xyz,euler)
            % plotGeometry(obj,figureNo,triColor,colorlim,method,extrapmethod)
            %   機体メッシュもしくは状態量をプロットする。
            %   カラー値は各パネルごとの値（例えばCp）を入れると、接点の値に補間しプロットする
            %
            % Inputs:
            %   - obj: オブジェクト自体
            %   - figureNo: プロットするフィギュアの番号
            %   - triColor: 各パネルのカラー値（例えばCp）
            %   - colorlim: カラーマップの範囲
            %   - method: カラーマップ補間の方法（"exact"または"linear"）
            %   - extrapmethod: カラーマップ補間の外挿方法（"exact"または"linear"）
            %
            % Usage:
            %   obj.plotGeometry(figureNo,triColor,colorlim,method,extrapmethod)
            %
            % Description:
            %   この関数は、機体メッシュもしくは状態量をプロットするために使用されます。
            %   カラー値は各パネルごとの値（例えばCp）を入れると、接点の値に補間しプロットします。
            %   プロットするフィギュアの番号、カラーマップの範囲、カラーマップ補間の方法、カラーマップ補間の外挿方法を指定することができます。
            %   カラーマップ補間の方法は、"exact"または"linear"を選択できます。
            %   カラーマップ補間の外挿方法は、"exact"または"linear"を選択できます。
            %
            % Example:
            %   obj.plotGeometry(1, Cp, [-1, 1], "exact", "linear")
            %
            % See also: trisurf, colormap
            if not(isa(figureNo, 'matlab.graphics.axis.Axes'))
                figure(figureNo);clf;
                ax = axes(gcf);
            else
                ax = figureNo;
            end
            if nargin<3
                trisurf(obj.tri,"FaceAlpha",0 ,'Parent',ax);
                if obj.halfmesh == 1
                    hold(ax,"on");
                    trisurf(obj.tri.ConnectivityList,obj.tri.Points(:,1),-obj.tri.Points(:,2),obj.tri.Points(:,3),"FaceAlpha",0,'Parent',ax);
                    hold(ax,"off");
                end
            else
                if nargin == 4
                    method = "linear";
                    xyz = [0,0,0];
                    euler = [0,0,0];
                elseif nargin == 5
                    xyz = [0,0,0];
                    euler = [0,0,0];
                end
                
                roll = euler(1);
                pitch = euler(2);
                yaw = euler(3);
                R = [cosd(yaw)*cosd(pitch), cosd(yaw)*sind(pitch)*sind(roll) - sind(yaw)*cosd(roll), cosd(yaw)*sind(pitch)*cosd(roll) + sind(yaw)*sind(roll);
                    sind(yaw)*cosd(pitch), sind(yaw)*sind(pitch)*sind(roll) + cosd(yaw)*cosd(roll), sind(yaw)*sind(pitch)*cosd(roll) - cosd(yaw)*sind(roll);
                    -sind(pitch), cosd(pitch)*sind(roll), cosd(pitch)*cosd(roll)];
                    
                affineverts = obj.tri.Points*R'+repmat(xyz,size(obj.tri.Points,1),1);

                if strcmpi(method,"exact")
                    %指定がexactの場合はパネル一枚一枚描画する
                    hold(ax,"on");
                    for i = 1:numel(obj.surfID)
                        if obj.paneltype(i) == 1
                            trisurf([1,2,3],affineverts(obj.tri(i,:),1),affineverts(obj.tri(i,:),2),affineverts(obj.tri(i,:),3),triColor(i),'EdgeAlpha',0.15,'Parent',ax);
                            if obj.halfmesh == 1
                                trisurf([1,2,3],affineverts(obj.tri(i,:),1),-affineverts(obj.tri(i,:),2),affineverts(obj.tri(i,:),3),triColor(i),'EdgeAlpha',0.15,'Parent',ax);
                            end
                        else
                            trisurf([1,2,3],affineverts(obj.tri(i,:),1),affineverts(obj.tri(i,:),2),affineverts(obj.tri(i,:),3),0,'EdgeAlpha',0.15,'Parent',ax);
                            if obj.halfmesh == 1
                                trisurf([1,2,3],affineverts(obj.tri(i,:),1),-affineverts(obj.tri(i,:),2),affineverts(obj.tri(i,:),3),0,'EdgeAlpha',0.15,'Parent',ax);
                            end
                        end
                    end
                    colormap(ax,"jet");
                    hold(ax,"off");
                    caxis(ax,colorlim);
                else
                    c = obj.verts2centerMat'*triColor;
                    trisurf(obj.tri.ConnectivityList,affineverts(:,1),affineverts(:,2),affineverts(:,3),c,'FaceColor','interp','EdgeAlpha',0.15,'Parent',ax);
                    colormap(ax,"jet");
                    if obj.halfmesh == 1
                        hold(ax,"on");
                        trisurf(obj.tri.ConnectivityList,affineverts(:,1),-affineverts(:,2),affineverts(:,3),c,'FaceColor','interp','EdgeAlpha',0.15,'Parent',ax);
                        colormap(ax,"jet");
                        hold(ax,"off");
                    end
                    caxis(ax,colorlim);
                end
            end
            axis(ax,"equal");xlabel(ax,"x");ylabel(ax,"y");zlabel(ax,"z");grid(ax,"on")
            drawnow();pause(0.1);
        end

 

        function obj = setMesh(obj, verts, connectivity, surfID, wakelineID)
            % setMesh - メッシュを設定する関数
            %   obj = setMesh(obj, verts, connectivity, surfID, wakelineID) は、UNLSIオブジェクトのメッシュを設定するための関数です。
            %   vertsは頂点座標の行列、connectivityは頂点の接続情報の行列、surfIDはパネルの表面IDのベクトル、wakelineIDはwakeのエッジIDのセル配列です。
            %   関数は、メッシュの各パネルに関連する情報を設定し、必要に応じてメッシュの修正を行います。
            %   返り値は、メッシュが設定されたUNLSIオブジェクトです。
            % ワーニングメッセージをオフにする
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
            % vertsのチェックとマージ
            [verts, connectivity, wakelineID] = obj.mergeVerts(obj.settingUNLSI.vertsTol, verts, connectivity, wakelineID);
            
            % 三角形分割オブジェクトを作成
            obj.tri = triangulation(connectivity, verts);
            
            % surfIDとwakelineIDを設定
            obj.surfID = surfID;
            obj.CpLimit = repmat(obj.settingUNLSI.defaultCpLimit,[numel(obj.surfID),1]);
            obj.wakelineID = wakelineID;
            
            % wakelineを初期化
            obj.wakeline = [];
            
            % wakelineのエッジ情報を設定
            for i = 1:numel(wakelineID)
                obj.wakeline{i}.edge = wakelineID{i};
            end
            
            % deflAngleが空の場合、初期化
            if isempty(obj.deflAngle)
                eyebuff = obj.rod2dcm([0,1,0],0);
                obj.deflAngle = [unique(obj.surfID), zeros(numel(unique(obj.surfID)), 1),ones(numel(unique(obj.surfID)), 1),zeros(numel(unique(obj.surfID)), 2), repmat(eyebuff(:)', [numel(unique(obj.surfID)), 1])];
            end
            
            % paneltypeとcpcalctypeを初期化
            obj.paneltype = ones(size(connectivity, 1), 1);
            obj.cpcalctype = ones(size(connectivity, 1), 1);
            
            % IndexPanel2Solverを初期化
            obj.IndexPanel2Solver = 1:numel(obj.paneltype);
            
            % clusterを初期化
            obj.cluster = cell([1, numel(obj.paneltype)]);
            
            % area, orgNormal, modNormal, centerを初期化
            obj.area = zeros(numel(obj.paneltype), 1);
            obj.orgNormal = zeros(numel(obj.paneltype), 3);
            obj.modNormal = zeros(numel(obj.paneltype), 3);
            obj.center = zeros(numel(obj.paneltype), 3);
            obj.rotCenter = zeros(numel(obj.paneltype), 3);
            obj.angularVelocity = zeros(numel(obj.paneltype), 3);
            
            % Cp, Cfe, LHS, RHSを初期化
            obj.Cp = {};
            obj.Cfe = {};
            obj.LHS = [];
            obj.RHS = [];
            obj.femutils = [];
            obj.femLHS = [];
            obj.femMass = [];
            obj.femDamp = [];
            obj.femtri = [];
            obj.femID = [];
            obj.femarea = [];
            obj.femNormal = [];
            obj.femcenter = [];
            obj.femThn = [];
            obj.femE = [];
            obj.femRho = [];
            obj.femVisc = [];
            obj.fem2aeroMat = [];
            obj.aero2femMat = [];

            
            % LLTとLLTnInterpを初期化
            lltninterp = obj.settingUNLSI.LLTnInterp;
            obj.LLT = [];
            obj.settingUNLSI.LLTnInterp = lltninterp;
            obj.modNormal = zeros(numel(obj.paneltype),3);
            obj.center = (verts(obj.tri.ConnectivityList(:,1),:)+verts(obj.tri.ConnectivityList(:,2),:)+verts(obj.tri.ConnectivityList(:,3),:))./3;
            obj.area = 0.5 * sqrt(sum(cross(obj.tri.Points(obj.tri.ConnectivityList(:, 2), :) - obj.tri.Points(obj.tri.ConnectivityList(:, 1), :), obj.tri.Points(obj.tri.ConnectivityList(:, 3), :) - obj.tri.Points(obj.tri.ConnectivityList(:, 1), :), 2).^2, 2));
            obj.orgNormal = faceNormal(obj.tri);
            for i = 1:numel(obj.paneltype)
                obj.modNormal(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])*obj.orgNormal(i,:)')';
            end
%             indices = arrayfun(@(x) find(obj.deflAngle(:, 1) == x), obj.surfID(:, 1), 'UniformOutput', false);      
%             deflMatrices = cell2mat(cellfun(@(idx) reshape(obj.deflAngle(idx, 6:14), [3, 3]), indices, 'UniformOutput', false));
%             obj.modNormal = cell2mat(arrayfun(@(i) (deflMatrices((i-1)*3+1:i*3, :) * obj.orgNormal(i, :)')', 1:numel(obj.paneltype), 'UniformOutput', false));
            
            % 半裁メッシュの境界表面上のwakeを削除
            deleteiter = [];
            for iter = 1:numel(obj.wakeline)
                deletelist = [];
                for j = 1:numel(obj.wakeline{iter}.edge) - 1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{iter}.edge(j), obj.wakeline{iter}.edge(j + 1));
                    if numel(attachpanel{1}) < 2
                        deletelist = [deletelist, j];
                    end
                end
                obj.wakeline{iter}.edge(deletelist) = [];
                if numel(obj.wakeline{iter}.edge) < 3
                    deleteiter = [deleteiter, iter];
                end
            end
            obj.wakeline(deleteiter) = [];
            
            % wakeのつくパネルIDを特定
            wakesurfID = [];
            for wakeNo = 1:numel(obj.wakeline)
                panelNo = 1;
                obj.wakeline{wakeNo}.valid = zeros(1, numel(obj.wakeline{wakeNo}.edge) - 1);
                for edgeNo = 1:numel(obj.wakeline{wakeNo}.edge) - 1
                    attachpanel = obj.tri.edgeAttachments(obj.wakeline{wakeNo}.edge(edgeNo), obj.wakeline{wakeNo}.edge(edgeNo + 1));
                    wakesurfID = [wakesurfID, setdiff(obj.surfID(attachpanel{1}), wakesurfID)];
                    if numel(attachpanel{1}) == 2
                        obj.wakeline{wakeNo}.validedge(1, panelNo) = obj.wakeline{wakeNo}.edge(edgeNo);
                        obj.wakeline{wakeNo}.validedge(2, panelNo) = obj.wakeline{wakeNo}.edge(edgeNo + 1);
                        obj.wakeline{wakeNo}.valid(edgeNo) = 1;
                        phivert = obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo + 1), 2:3) - obj.tri.Points(obj.wakeline{wakeNo}.edge(edgeNo), 2:3);
                        phiwake = atan2d(phivert(2), phivert(1));
                        if abs(phiwake) < 90
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
            
            %デフォルトのwake形状を定義する
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    xwakeshape = linspace(0,obj.settingUNLSI.defaultWakeLength*obj.CREF,6);
                    obj.wakeline{wakeNo}.wakeShape{edgeNo} =  [xwakeshape(2:end)',zeros(5,2)];
                end
                obj.wakeline{wakeNo}.wakeOrg = obj.getWakePosition(wakeNo);
            end
 

            % propが存在する場合、設定を反映
            if numel(obj.prop) > 0
                for propNo = 1:numel(obj.prop)
                    propstate = [obj.prop{propNo}.CT, obj.prop{propNo}.CP, obj.prop{propNo}.rpm];
                    obj = obj.setProp(propNo, obj.prop{propNo}.ID, obj.prop{propNo}.diameter, obj.prop{propNo}.XZsliced);
                    obj = obj.setPropState(propNo, propstate(1), propstate(2), propstate(3));
                end
            end
            
            % wakeSurfIDにないパネルのcpcalctypeをlinearに設定
            for i = setdiff(unique(obj.surfID)', wakesurfID)
                obj = obj.setCpCalcType(i, "linear");
            end
            
            % メッシュのチェックと修正
            fprintf("\n Panel Check\n");
            fprintf("Minimum area = %d\n",min(obj.area));
            fprintf('Maximum area = %d\n',max(obj.area));
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol, "delete");
            % obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol, "warning");
            obj = obj.checkWakeIntersect();
            % ワーニングメッセージをオンにする
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');
            
            % 空力のパネル中心座標からメッシュノードへの投影変換行列を作成
            phi = @(r) obj.phi3(r, 1);
            Avv = phi(obj.calcRMat(obj.tri.Points, obj.tri.Points));
            Acv = phi(obj.calcRMat(obj.center, obj.tri.Points));
            obj.verts2centerMat = Acv * pinv(Avv);
            
            % メッシュのチェックと修正
            obj = obj.checkMesh(obj.settingUNLSI.checkMeshTol, "delete");
            obj = obj.checkWakeIntersect();
            % ワーニングメッセージをオンにする
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');
            
            % 空力のパネル中心座標からメッシュノードへの投影変換行列を作成
            phi = @(r) obj.phi3(r, 1);
            Avv = phi(obj.calcRMat(obj.tri.Points, obj.tri.Points));
            Acv = phi(obj.calcRMat(obj.center, obj.tri.Points));
            obj.verts2centerMat = Acv * pinv(Avv);
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
            modVerts = any((verts-obj.tri.Points)~=0,2);
            modTri = vertexAttachments(obj.tri,find(modVerts));
            obj.tri = triangulation(con,verts);
            index = unique(horzcat(modTri{:}));
            obj.center(index,:) = (verts(obj.tri.ConnectivityList(index,1),:)+verts(obj.tri.ConnectivityList(index,2),:)+verts(obj.tri.ConnectivityList(index,3),:))./3;
            obj.area(index,1) = 0.5 * sqrt(sum(cross(obj.tri.Points(obj.tri.ConnectivityList(index, 2), :) - obj.tri.Points(obj.tri.ConnectivityList(index, 1), :), obj.tri.Points(obj.tri.ConnectivityList(index, 3), :) - obj.tri.Points(obj.tri.ConnectivityList(index, 1), :), 2).^2, 2));
            obj.orgNormal(index,:) = faceNormal(obj.tri,index(:));
            for i = 1:numel(index)
                obj.modNormal(index(i),:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(index(i),1)),6:14),[3,3])*obj.orgNormal(index(i),:)')';
            end
            %indices = arrayfun(@(x) find(obj.deflAngle(:, 1) == x), obj.surfID(:, 1), 'UniformOutput', false);      
            %deflMatrices = cell2mat(cellfun(@(idx) reshape(obj.deflAngle(idx, 6:14), [3, 3]), indices, 'UniformOutput', false));
            %obj.modNormal = cell2mat(arrayfun(@(i) (deflMatrices((i-1)*3+1:i*3, :) * obj.orgNormal(i, :)')', 1:numel(obj.paneltype), 'UniformOutput', false));
        end

        function obj = exportMesh(obj,filetype,filename,uv1,uv2,uv3)
            %%%%%%%%%%%%%%%%%メッシュのエクスポート%%%%%%%%%%%%%%%
            %メッシュをtriフォーマットかvspgeomフォーマットでエクスポートする
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmpi(filetype,"tri")
                filename = strcat(filename,".tri");

                verts = obj.tri.Points;
                connectivity = obj.tri.ConnectivityList;
                surfID = obj.surfID;
                [prows, pcols] = size(verts);
                [crows, ccols] = size(connectivity);
                [srows, scols] = size(surfID);
                fid = fopen(filename, 'w');
                if fid == -1
                    error('ファイルを開けませんでした: %s', filename);
                end
                fprintf(fid, '%d %d \n', prows,crows);%heder
                for i = 1:prows
                    for j = 1:pcols
                        fprintf(fid, '%.6f ', verts(i, j)); % データをスペース区切りで書き込む
                    end
                    fprintf(fid, '\n'); % 行の終わりに改行を追加
                end
                
                
                for i = 1:crows
                    for j = 1:ccols
                        fprintf(fid, '%d ', connectivity(i, j)); % データをスペース区切りで書き込む
                    end
                    fprintf(fid, '\n'); % 行の終わりに改行を追加
                end

                for i = 1:srows
                    for j = 1:scols
                        fprintf(fid, '%d ', surfID(i, j)); % データをスペース区切りで書き込む
                    end
                    fprintf(fid, '\n'); % 行の終わりに改行を追加
                end
                fclose(fid);
            elseif strcmpi(filetype,"vspgeom")
                fid = fopen(filename, 'w');
                if fid == -1
                    error('ファイルを開けませんでした: %s', filename);
                end
                fprintf(fid,'%d\n',size(obj.tri.Points,1));%節点数       
                for i = 1:size(obj.tri.Points,1)
                    fprintf(fid,'%f %f %f\n',obj.tri.Points(i,1),obj.tri.Points(i,2),obj.tri.Points(i,3));
                end         
                fprintf(fid,'%d\n',size(obj.tri.ConnectivityList,1));%三角形要素数
                for i = 1:size(obj.tri.ConnectivityList,1)
                    fprintf(fid,'3 %d %d %d\n',obj.tri.ConnectivityList(i,1),obj.tri.ConnectivityList(i,2),obj.tri.ConnectivityList(i,3));
                end

                %surf ID
                for i = 1:size(obj.tri.ConnectivityList,1)
                    fprintf(fid,'%d %f %f %f %f %f %f\n',obj.surfID(i),uv1(i,1),uv1(i,2),uv2(i,1),uv2(i,2),uv3(i,1),uv3(i,2));
                end

                %wake
                fprintf(fid,'%d\n',size((obj.wakelineID),1));%number of wakes
                   %以下copilot 要確認 
                for i = 1:size(obj.wakelineID,1)
                    fprintf(fid,'%d\n',size(obj.wakelineID{i},2));%number of edges
                    for j = 1:size(obj.wakelineID{i},2)
                        fprintf(fid,'%d\n',obj.wakelineID{i}(j));
                    end
                end
                error("未作成");
            else
                error("Invalid filetype");
            end
        end
        

        
        function obj = checkWakeIntersect(obj)
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    for i = 1:size(obj.wakeline{wakeNo}.wakeShape{edgeNo},1)
                        if i == 1
                            surface1.vertices(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:);
                            surface1.vertices(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:);
                            surface1.vertices(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i,:);
                            surface1.vertices(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i,:);
                        else
                            surface1.vertices(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i-1,:);
                            surface1.vertices(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i-1,:);
                            surface1.vertices(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i,:);
                            surface1.vertices(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+obj.wakeline{wakeNo}.wakeShape{edgeNo}(i,:);
                        end
                        surface1.faces = [1,2,3;3,4,1];
                        surface2.vertices = obj.tri.Points;
                        surface2.faces = obj.tri.ConnectivityList(obj.paneltype == 1,:);
                        [~, Surf12] = obj.SurfaceIntersection(obj,surface1, surface2);
                        if i == 1
                            if size(Surf12.vertices,1) <= 2
                                obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) = 1;
                            else
                                obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) = 0;
                            end
                        else
                            if size(Surf12.vertices,1) < 2 
                                obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) = 1;
                            else
                                obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) = 0;
                            end
                        end
                    end
                    
                end
                
            end
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
                    obj = obj.setMesh(obj.tri.Points,newCon,obj.surfID,obj.wakelineID);
                else
                    error("some panel areas are too small");
                end
            end
        end
  

        function [verts, cons, wedata] = mergeVerts(obj, tol, verts, cons, wedata)
    
         % 
            % mergeVerts関数は、指定された許容誤差(tol)を使用して、頂点(verts)と接続(cons)をマージします。
            % マージされた結果の頂点と接続は、vertsとconsの出力引数として返されます。
            % 
            % 入力:
            %   - tol: マージするための許容誤差
            %   - verts: 頂点の座標行列
            %   - cons: 頂点の接続行列
            %   - wedata (オプション): 頂点の追加データ行列
            %
            % 出力:
            %   - verts: マージ後の頂点の座標行列
            %   - cons: マージ後の頂点の接続行列
            %   - wedata: マージ後の頂点の追加データ行列 (オプション)
            %
            % 例:
            %   tol = 0.1;
            %   verts = [0 0; 1 1; 2 2; 3 3];
            %   cons = [1 2; 2 3; 3 4; 4 1];
            %   [mergedVerts, mergedCons] = mergeVerts(tol, verts, cons);
            %
            %   出力:
            %   mergedVerts =
            %       0 0
            %       1 1
            %       3 3
            %   
            %   mergedCons =
            %       1 2
            %       2 3
            %       3 1
            %
            % この関数は、指定された許容誤差内で頂点をマージし、重複する頂点を削除します。
            % マージされた頂点に対応する接続も更新されます。
            % オプションのwedata引数が指定された場合、マージされた頂点に対応する追加データも更新されます。
            % マージが行われた場合、警告メッセージが表示されます。
            % マージが行われなかった場合、ループを終了します。
            % オプションのwedata引数が指定されていない場合、wedataは空の行列として返されます。
            % 頂点と接続から三角形を作成
            mergetri = triangulation(cons, verts);
            
            % 頂点のマージを繰り返す
            for iter = 1:size(verts, 1)
                modflag = 0;
                for i = 1:size(verts, 1)
                    % 許容誤差内で重複する頂点を検索
                    index = find(vecnorm(verts(i, :) - verts(i+1:end, :), 2, 2) < tol);
                    if ~isempty(index)
                        j = index(1) + i;
                        Vbuff = vertexAttachments(mergetri, j);
                        V = Vbuff{1};
                        
                        % 接続を更新
                        for k = 1:numel(V)
                            cons(V(k), cons(V(k), :) == j) = i;
                        end
                        
                        % 頂点を削除
                        verts(j, :) = [];
                        
                        % 頂点のインデックスを更新
                        cons(cons(:) > j) = cons(cons(:) > j) - 1;
                        
                        % オプションのwedata引数が指定された場合、追加データを更新
                        if nargin > 4
                            for k = 1:numel(wedata)
                                wedata{k}(wedata{k} == j) = i;
                                wedata{k}(wedata{k} > j) = wedata{k}(wedata{k} > j) - 1;
                                wedata{k} = unique(wedata{k}, 'stable');
                            end
                        end
                        
                        modflag = 1;
                        warning("points(%d and %d) are merged", i, j);
                        break;
                    else
                        modflag = 0;
                    end
                end
                
                % 三角形を更新
                mergetri = triangulation(cons, verts);
                
                % マージが行われなかった場合、ループを終了
                if modflag == 0
                    break;
                end
            end
            
            % オプションのwedata引数が指定されていない場合、wedataは空の行列とする
            if nargin < 4
                wedata = [];
            end
        end
   
        function obj = checkFemMesh(obj,tol,outType)
            %%%%%%%%%%%%%%%%パネル面積の確認%%%%%%%%%%%%%%%
            %tol:最小許容面積
            %outType:"warning"を指定するとエラーの代わりに警告を表示する
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                outType = "error";
            end
            if any(obj.femarea<tol)
                if strcmpi(outType,'warning')
                    warning("some panel areas are too small");
                elseif strcmpi(outType,'delete')
                    
                    deleteIndex = obj.femarea<tol;
                    if any(deleteIndex)
                        warning("some panel areas are too small and deleted");
                    end
                    newCon = obj.femtri.ConnectivityList;
                    newCon(deleteIndex,:) = [];
                    obj.femID(deleteIndex,:) = [];
                    obj = obj.setFemMesh(obj.femtri.Points,newCon,obj.femID,obj.femutils.aeroIDLink);
                else
                    error("some panel areas are too small");
                end
            end
            fprintf("\n Panel Check\n");
            fprintf("Minimum FEM area=%d\n",min(obj.femarea));
            fprintf("Maximum FEM area=%d\n",max(obj.femarea));
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
            if nargin == 4
                ID = unique(obj.surfID);
            end
            flowVec(1) = cosd(alpha)*cosd(beta);
            flowVec(2) = -sind(beta);
            flowVec(3) = sind(alpha)*cosd(beta);
            for i = 1:numel(obj.paneltype)
                if obj.paneltype(i) == 1 && any(obj.surfID(i) == ID(:)')
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
            % makeEquation関数はパネル法連立方程式行列を作成する関数です。
            % パネル法の根幹であり、表面微分行列も併せて作成します。
            % 
            % 入力:
            %   - obj: UNLSIオブジェクト
            %
            % 出力:
            %   - obj: パネル法連立方程式行列が作成されたUNLSIオブジェクト
            %
            % 使用例:
            %   obj = UNLSI();
            %   obj = obj.makeEquation();
            %
            
            % パネル数とバウンダリーパネル数を取得
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            
            % パネル微分行列の初期化
            obj.mu2v{1} = sparse(nPanel,nbPanel);
            obj.mu2v{2} = sparse(nPanel,nbPanel);
            obj.mu2v{3} = sparse(nPanel,nbPanel);
            
            % パネル法連立方程式行列の作成
            for i = 1:nPanel
                if obj.paneltype(i) == 1
                    CPmat = obj.center(obj.cluster{i},1:3);
                    pnt = obj.center(i,:);
                    m = obj.tri.Points(obj.tri.ConnectivityList(i,1),:)'-pnt(:);
                    m = m./norm(m);
                    l = cross(m,obj.orgNormal(i,:)');
                    Minv = [l,m,obj.orgNormal(i,:)'];
                    % Minvinv_temp = pinv(Minv);
                    % B = CPmat - repmat(pnt,[size(CPmat,1),1]);
                    % lmnMat = lsqminnorm(Minv',B')';
                    lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                    % lmnMat = (Minvinv_temp*(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                    bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                    Bmat=pinv(bb,sqrt(eps));
                    Vnmat = Minv*[1,0,0,0;0,1,0,0;0,0,1,0;]*Bmat;
                    for iter = 1:3
                        obj.mu2v{iter}(i,obj.IndexPanel2Solver(obj.cluster{i})) = Vnmat(iter,:);
                    end
                end
            end

            % 誘導抗力計算用の行列の作成
            
            % パネル方連立方程式行列の作成
            % 機体パネル⇒機体パネルへの影響
            si = floor(nbPanel/obj.settingUNLSI.nCalcDivide).*(0:obj.settingUNLSI.nCalcDivide-1)+1;
            ei = [floor(nbPanel/obj.settingUNLSI.nCalcDivide).*(1:obj.settingUNLSI.nCalcDivide-1),nbPanel];
            obj.LHS = zeros(nbPanel);
            obj.wakeLHS = zeros(nbPanel);
            obj.RHS = zeros(nbPanel);

            for i= 1:obj.settingUNLSI.nCalcDivide
                [~,~,VortexAc,VortexBc] = obj.influenceMatrix(obj,[],si(i):ei(i));
                obj.LHS(:,si(i):ei(i)) = VortexAc; 
                obj.RHS(:,si(i):ei(i)) = VortexBc;
            end
            obj.LHS(~isfinite(obj.LHS)) = 0;
            obj.RHS(~isfinite(obj.RHS)) = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [U, S, V] = svds(obj.LHS,size(obj.LHS,1));
            obj.LHS = U*S*V';
            %wakeパネル⇒機体パネルへの影響
            %
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                    interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                    [influence] = obj.wakeInfluenceMatrix(obj,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                    obj.wakeLHS(:,interpID(1)) = obj.wakeLHS(:,interpID(1)) - influence;
                    obj.wakeLHS(:,interpID(2)) = obj.wakeLHS(:,interpID(2)) + influence;
                end
            end

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

        function obj = setWakeShape(obj,wakeShape,ID)
            %%%%%%%%%%%%%%%%wake形状の指定関数%%%%%%%%%%%%%%%%
            %wakeNo : wake番号
            %edgeNo : edge番号
            %wakeShape : [x,y,z]の形状
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                for wakeNo = 1:numel(obj.wakeline)
                    for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                        obj.wakeline{wakeNo}.wakeShape{edgeNo} = wakeShape;
                    end
                end
            else
                for wakeNo = 1:numel(obj.wakeline)
                    for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                        if any(obj.surfID(obj.wakeline{wakeNo}.upperID(edgeNo))==ID) && any(obj.surfID(obj.wakeline{wakeNo}.lowerID(edgeNo))==ID) 
                            obj.wakeline{wakeNo}.wakeShape{edgeNo} = wakeShape;
                        end
                    end
                end
            end
            obj = obj.checkWakeIntersect();
        end


        function obj = setHelixWake(obj,ID,dt,rpm,rotAxis,rotOrigin)
            %%%%%%%%%%%%%%%%wake形状の指定関数%%%%%%%%%%%%%%%%
            %wakeShape : [x,y,z]の形状
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K = [0, -rotAxis(3), rotAxis(2);
                 rotAxis(3), 0, -rotAxis(1);
                -rotAxis(2), rotAxis(1), 0];
            rotAngle = -rpm*pi/30 * dt;
            dcm = eye(3) + sin(rotAngle) * K + (1 - cos(rotAngle)) * K^2;
            for wakeIter = 1:numel(obj.wakeline)
                wakeOrg = obj.getWakePosition(wakeIter);
                for edgeNo = 1:size(obj.wakeline{wakeIter}.validedge,2)
                    if any(obj.surfID(obj.wakeline{wakeIter}.upperID(edgeNo))==ID) && any(obj.surfID(obj.wakeline{wakeIter}.lowerID(edgeNo))==ID) 
                        for i = 1:obj.settingUNLSI.nWakeMax
                            if i == 1
                                nextR = (wakeOrg(edgeNo,:)-rotOrigin)*dcm-(wakeOrg(edgeNo,:)-rotOrigin);
                                obj.wakeline{wakeIter}.wakeShape{edgeNo}(i,1:3) = nextR + dt * obj.settingUNLSI.Vinf.*rotAxis(:)';
                            else
                                nextR = (obj.wakeline{wakeIter}.wakeShape{edgeNo}(i-1,1:3)+wakeOrg(edgeNo,:)-rotOrigin)*dcm-wakeOrg(edgeNo,:)+rotOrigin;
                                obj.wakeline{wakeIter}.wakeShape{edgeNo}(i,1:3) = nextR + dt * obj.settingUNLSI.Vinf.*rotAxis(:)';
                            end
                        end
                    end
                end
            end
            obj = obj.checkWakeIntersect();
            obj = obj.makeWakeEquation();
        end

        function obj = marchWake(obj,dt,alpha,beta,Mach,Re)
            %%%%%%%%%%%%%%%%wake形状の更新関数%%%%%%%%%%%%%%%%
            %flowVec : 流速ベクトル
            %dt : 時間ステップ
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(obj.flowNoTable)
                if nargin == 6
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
                    if nargin == 6
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
            
            wakePos = {};
            wakeOrgBuff = {};
            snorm = [];
            wakeNormal = [];
            for wakeNo = 1:numel(obj.wakeline)
                wakeOrg{wakeNo} = obj.getWakePosition(wakeNo);
                for i = 1:size(wakeOrg{wakeNo},1)
                    wakePos = [wakePos,obj.wakeline{wakeNo}.wakeOrg(i,:)+obj.wakeline{wakeNo}.wakeShape{i}];
                    wakeOrgBuff = [wakeOrgBuff,wakeOrg{wakeNo}(i,:)];
                    snorm = [snorm;norm(obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,i),:)-obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,i),:))];
                    normBuff = (obj.modNormal(obj.wakeline{wakeNo}.upperID(i),:)+obj.modNormal(obj.wakeline{wakeNo}.lowerID(i),:))./2;
                    wakeNormal = [wakeNormal;normBuff/norm(normBuff)];
                end
            end
            u0 = obj.solvePertPotential(flowNo,alpha,beta);
            initPos = [];
            for iter = 1:numel(wakePos)
                if size(wakePos{1},1) == obj.settingUNLSI.nWakeMax
                    nWake(iter) = size(wakePos{iter},1)-1;
                else
                    nWake(iter) = size(wakePos{iter},1);
                end
                for k = nWake(iter):-1:1
                    initPos= [initPos;wakePos{iter}(k,:)];
                end
            end
            flowVec = zeros(size(initPos,1),3);
            for i = 1:size(initPos,1)
               %rvec = obj.center(i,:)'-obj.rotCenter(i,:)';
               flowVec(i,:) = obj.settingUNLSI.Vinf.*[cosd(alpha)*cosd(beta),-sind(beta),sind(alpha)*cosd(beta)];%-(cross(obj.angularVelocity(i,:)./180.*pi,rvec(:)')')';
            end
           [VmuWake,VmuBody]  = obj.makeVelocityInfluence(initPos);
            vel = obj.calcVelocity(u0,VmuWake,VmuBody);
            vel = obj.settingUNLSI.Vinf.*vel+flowVec;
            k1 = dt .* vel(:);
            [VmuWake,VmuBody] = obj.makeVelocityInfluence(reshape(initPos(:)+0.5.*k1,[],3));
            vel = obj.calcVelocity(u0,VmuWake,VmuBody);
            vel = obj.settingUNLSI.Vinf.*vel+flowVec;
            k2 = dt .* vel(:);
            [VmuWake,VmuBody] = obj.makeVelocityInfluence(reshape(initPos(:)+0.5.*k2,[],3));
            vel = obj.calcVelocity(u0,VmuWake,VmuBody);
            vel = obj.settingUNLSI.Vinf.*vel+flowVec;
            k3 = dt .* vel(:);
            [VmuWake,VmuBody] = obj.makeVelocityInfluence(reshape(initPos(:)+k3,[],3));
            vel = obj.calcVelocity(u0,VmuWake,VmuBody);
            vel = obj.settingUNLSI.Vinf.*vel+flowVec;
            k4 = dt .* vel(:);
            newPos = reshape(initPos(:)+(k1+2.*k2+2.*k3+k4)./6,[],3);
%                 options = odeset("RelTol",1e-2,"AbsTol",1e-2);
%                 [~,yspan] = ode113(@(t,y)obj.marchWakeShapeODE(t,y,u0,flowVec,obj.settingUNLSI.Vinf),[0,dt],initPos(:),options);
%                 newPos = reshape(yspan(end,:)',[],3);
            
            for iter = 1:numel(wakePos)
                for k = nWake(iter):-1:1
                    wakePos{iter}(k+1,:) = newPos(1,:);
                    newPos(1,:) = [];
                end
            end
            for i = 1:numel(wakePos)
                wakePos{i}(1,:) = wakeOrgBuff{i}(1,:) + snorm(i,1)*obj.settingUNLSI.xWakeAttach.*wakeNormal(i,:);
            end

            id = 1;
            for wakeNo = 1:numel(obj.wakeline)
                for iter = 1:size(wakeOrg{wakeNo},1)
                    obj.wakeline{wakeNo}.wakeShape{iter}=wakePos{id}-wakeOrgBuff{id}(1,:);
                    id = id+1;
                end
                obj.wakeline{wakeNo}.wakeOrg = wakeOrg{wakeNo};
            end
            obj = obj.checkWakeIntersect();       
            obj = obj.makeWakeEquation();
        end

        function dy = marchWakeShapeODE(obj,t,y,u0,flowVec,Vnorm)
            %%%%%%%%%%%%%%%%wake形状の更新関数%%%%%%%%%%%%%%%%
            %flowVec : 流速ベクトル
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [VmuWake,VmuBody] = obj.makeVelocityInfluence(reshape(y,[],3));
            vel = obj.calcVelocity(u0,VmuWake,VmuBody);
            vel = Vnorm.*vel+flowVec;
            dy = vel(:);
        end

        function controlPoint = getWakePosition(obj,wakeNo)
            %%%%%%%%%%%%%%%%wake位置の取得関数%%%%%%%%%%%%%%%%
            %wakeNo : wake番号
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            controlPoint = zeros(size(obj.wakeline{wakeNo}.validedge,2),3);
            for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                controlPoint(edgeNo,:) = (obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:))./2;
            end
        end

        function plotWakeShape(obj,figNo,xyz,euler)
            %%%%%%%%%%%%%%%%wake形状のプロット関数%%%%%%%%%%%%%%%%
            %figNo : figure番号
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                xyz = [0,0,0];
                euler = [0,0,0];
            elseif nargin == 3
                euler = [0,0,0];
            end
            
            roll = euler(1);
            pitch = euler(2);
            yaw = euler(3);
            R = [cosd(yaw)*cosd(pitch), cosd(yaw)*sind(pitch)*sind(roll) - sind(yaw)*cosd(roll), cosd(yaw)*sind(pitch)*cosd(roll) + sind(yaw)*sind(roll);
                sind(yaw)*cosd(pitch), sind(yaw)*sind(pitch)*sind(roll) + cosd(yaw)*cosd(roll), sind(yaw)*sind(pitch)*cosd(roll) - cosd(yaw)*sind(roll);
                -sind(pitch), cosd(pitch)*sind(roll), cosd(pitch)*cosd(roll)];
                
                
            figure(figNo);
            hold on;
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    wakeLine1 = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+[0,0,0;obj.wakeline{wakeNo}.wakeShape{edgeNo}];
                    wakeLine1 = wakeLine1 * R' + xyz(:)';
                    wakeLine2 = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+[0,0,0;obj.wakeline{wakeNo}.wakeShape{edgeNo}];
                    wakeLine2 = wakeLine2 * R' + xyz(:)';
                    plot3(wakeLine1(:,1),wakeLine1(:,2),wakeLine1(:,3),'r');
                    plot3(wakeLine2(:,1),wakeLine2(:,2),wakeLine2(:,3),'r');
                    for i = 1:size(wakeLine1,1)
                        plot3([wakeLine1(i,1);wakeLine2(i,1)],[wakeLine1(i,2);wakeLine2(i,2)],[wakeLine1(i,3);wakeLine2(i,3)],'r');
                    end
                    if obj.halfmesh == 1
                        plot3(wakeLine1(:,1),-wakeLine1(:,2),wakeLine1(:,3),'r');
                        plot3(wakeLine2(:,1),-wakeLine2(:,2),wakeLine2(:,3),'r');
                        for i = 1:size(wakeLine1,1)
                            plot3([wakeLine1(i,1);wakeLine2(i,1)],-[wakeLine1(i,2);wakeLine2(i,2)],[wakeLine1(i,3);wakeLine2(i,3)],'r');
                        end
                    end
                end
            end
            hold off;drawnow();
        end


        function obj = makeWakeEquation(obj)
            nbPanel = sum(obj.paneltype == 1);
            obj.wakeLHS = zeros(nbPanel);
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                    interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                    [influence] = obj.wakeInfluenceMatrix(obj,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                    obj.wakeLHS(:,interpID(1)) = obj.wakeLHS(:,interpID(1)) - influence;
                    obj.wakeLHS(:,interpID(2)) = obj.wakeLHS(:,interpID(2)) + influence;
                end
            end
        end

        function [VmuWake,VmuBody] = makeVelocityInfluence(obj,controlPoint)
        
            bPanel = sum(obj.paneltype == 1);
            VmuWake.X = zeros(size(controlPoint,1),bPanel);
            VmuWake.Y = zeros(size(controlPoint,1),bPanel);
            VmuWake.Z = zeros(size(controlPoint,1),bPanel);

            %wakeパネル⇒機体パネルへの影響
            for wakeNo = 1:numel(obj.wakeline)
                for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                    interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                    interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                    [influence] = obj.wakeVelocityInfluence(obj,wakeNo,edgeNo,controlPoint,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                    
                    VmuWake.X(:,interpID(1)) = VmuWake.X(:,interpID(1)) - influence.X;
                    VmuWake.X(:,interpID(2)) = VmuWake.X(:,interpID(2)) + influence.X;
                    VmuWake.Y(:,interpID(1)) = VmuWake.Y(:,interpID(1)) - influence.Y;
                    VmuWake.Y(:,interpID(2)) = VmuWake.Y(:,interpID(2)) + influence.Y;
                    VmuWake.Z(:,interpID(1)) = VmuWake.Z(:,interpID(1)) - influence.Z;
                    VmuWake.Z(:,interpID(2)) = VmuWake.Z(:,interpID(2)) + influence.Z;
                end
            end

            if nargout > 1
                VmuBody = obj.velocityInfluence(obj,controlPoint);
            end
        end

        function obj = calcApproximatedFemEquation(obj)
            % calcApproximatedEquationメソッドは、近似方程式を計算するための関数です。
            % 近似フラグが立っている場合は、エラーが発生します。
            % nbPanelは、paneltypeが1であるパネルの数を表します。
            % 接点に繋がるIDを決定し、近似行列の計算インデックスを設定します。
            % また、各頂点がwakeに含まれているかどうかを判定し、影響行列を計算します。
            % 最後に、微小変位に対する微分係数を計算します。  

            if obj.femapproximated == 1
                error("This instance is approximated. Please execute obj.makeEquation()");
            end

            % 接点に繋がるIDを決定
            vertAttach = obj.femtri.vertexAttachments();
            pert = eps^(1/3);
            for i = 1:numel(vertAttach)
                
                if mod(i,floor(numel(vertAttach)/10))==0 || i == 1
                    fprintf("%d/%d \n",i,numel(vertAttach));
                end
                % paneltype == 1以外を削除する
                obj.approxFemMat.calcIndex{i} = vertAttach{i};
                matsize = size(obj.femLHS);
                obj.approxFemMat.dLHSX{i} = sparse(matsize(1),matsize(2));
                obj.approxFemMat.dLHSY{i} = sparse(matsize(1),matsize(2));
                obj.approxFemMat.dLHSZ{i} = sparse(matsize(1),matsize(2));
                j=1;%X
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)+pert;
                obj2f = obj.setFemVerts(newVerts,1);
                
                obj2f= obj2f.makeFemEquation(obj.approxFemMat.calcIndex{i});
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)-pert;
                obj2r = obj.setFemVerts(newVerts,1);
                obj2r.femRho = obj.femRho;obj2r.femE = obj.femE;obj2r.femThn = obj.femThn;obj2r.femVisc = obj.femVisc;
                obj2r= obj2r.makeFemEquation(obj.approxFemMat.calcIndex{i});
                obj.approxFemMat.dLHSX{i} = (obj2f.femLHS-obj2r.femLHS)./pert/2;
                j=2;%Y
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)+pert;
                obj2f = obj.setFemVerts(newVerts,1);
                
                obj2f= obj2f.makeFemEquation(obj.approxFemMat.calcIndex{i});
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)-pert;
                obj2r = obj.setFemVerts(newVerts,1);
                obj2r.femRho = obj.femRho;obj2r.femE = obj.femE;obj2r.femThn = obj.femThn;obj2r.femVisc = obj.femVisc;
                obj2r= obj2r.makeFemEquation(obj.approxFemMat.calcIndex{i});
                obj.approxFemMat.dLHSY{i} = (obj2f.femLHS-obj2r.femLHS)./pert/2;
                j=3;%Z
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)+pert;
                obj2f = obj.setFemVerts(newVerts,1);
                
                obj2f= obj2f.makeFemEquation(obj.approxFemMat.calcIndex{i});
                newVerts = obj.femtri.Points;
                newVerts(i,j) = obj.femtri.Points(i,j)-pert;
                obj2r = obj.setFemVerts(newVerts,1);
                obj2r.femRho = obj.femRho;obj2r.femE = obj.femE;obj2r.femThn = obj.femThn;obj2r.femVisc = obj.femVisc;
                obj2r= obj2r.makeFemEquation(obj.approxFemMat.calcIndex{i});
                obj.approxFemMat.dLHSZ{i} = (obj2f.femLHS-obj2r.femLHS)./pert/2;
            end
        end


        function obj = calcApproximatedEquation(obj,wakeGradFlag)
            % calcApproximatedEquationメソッドは、近似方程式を計算するための関数です。
            % 近似フラグが立っている場合は、エラーが発生します。
            % nbPanelは、paneltypeが1であるパネルの数を表します。
            % 接点に繋がるIDを決定し、近似行列の計算インデックスを設定します。
            % また、各頂点がwakeに含まれているかどうかを判定し、影響行列を計算します。
            % 最後に、微小変位に対する微分係数を計算します。
            if nargin == 1
                wakeGradFlag = 0;
            end
            obj.approxMat.wakeGradFlag = wakeGradFlag;    

            if obj.approximated == 1
                error("This instance is approximated. Please execute obj.makeEquation()");
            end
            nbPanel = sum(obj.paneltype == 1);
            % 接点に繋がるIDを決定
            vertAttach = obj.tri.vertexAttachments();
            pert = eps^(1/3);
            for i = 1:numel(vertAttach)
                if mod(i,floor(numel(vertAttach)/10))==0 || i == 1
                    fprintf("%d/%d \n",i,numel(vertAttach));
                end
                % paneltype == 1以外を削除する
                vertAttach{i}(obj.paneltype(vertAttach{i}) ~=1) = [];
                obj.approxMat.calcIndex{i} = sort(obj.IndexPanel2Solver(vertAttach{i}));
                j=1;%X
                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)+pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArf,VortexBrf,VortexAcf,VortexBcf] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});
                
                if obj.approxMat.wakeGradFlag == 1
                    VortexArwf = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwf = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence;
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) + influence;
                                end
                            end
                        end
                    end
                end

                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)-pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArr,VortexBrr,VortexAcr,VortexBcr] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});

                if obj.approxMat.wakeGradFlag == 1
                    VortexArwr = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwr = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence;
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) + influence;
                                end
                            end
                        end
                    end
                end

                obj.approxMat.dVArX{i} = (VortexArf-VortexArr)./pert/2;
                obj.approxMat.dVBrX{i} = (VortexBrf-VortexBrr)./pert/2;
                obj.approxMat.dVAcX{i} = (VortexAcf-VortexAcr)./pert/2;
                obj.approxMat.dVBcX{i} = (VortexBcf-VortexBcr)./pert/2;
                if obj.approxMat.wakeGradFlag == 1
                    obj.approxMat.dVArwX{i} = (VortexArwf-VortexArwr)./pert/2;
                    obj.approxMat.dVAcwX{i} = (VortexAcwf-VortexAcwr)./pert/2;
                end
                
                j=2;%Y
                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)+pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArf,VortexBrf,VortexAcf,VortexBcf] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});
                
                if obj.approxMat.wakeGradFlag == 1
                    VortexArwf = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwf = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence;
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) + influence;
                                end
                            end
                        end
                    end
                end

                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)-pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArr,VortexBrr,VortexAcr,VortexBcr] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});
                
                if obj.approxMat.wakeGradFlag == 1
                    VortexArwr = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwr = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence;
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) + influence;
                                end
                            end
                        end
                    end
                end

                obj.approxMat.dVArY{i} = (VortexArf-VortexArr)./pert/2;
                obj.approxMat.dVBrY{i} = (VortexBrf-VortexBrr)./pert/2;
                obj.approxMat.dVAcY{i} = (VortexAcf-VortexAcr)./pert/2;
                obj.approxMat.dVBcY{i} = (VortexBcf-VortexBcr)./pert/2;
                if obj.approxMat.wakeGradFlag == 1
                    obj.approxMat.dVArwY{i} = (VortexArwf-VortexArwr)./pert/2;
                    obj.approxMat.dVAcwY{i} = (VortexAcwf-VortexAcwr)./pert/2;
                end

                j=3;%Z
                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)+pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArf,VortexBrf,VortexAcf,VortexBcf] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});
                
                if obj.approxMat.wakeGradFlag == 1
                    VortexArwf = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwf = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence;
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwf(:,interpID(1)) = VortexArwf(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwf(:,interpID(2)) = VortexArwf(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwf(:,b) = VortexAcwf(:,b) + influence;
                                end
                            end
                        end
                    end
                end
                
                newVerts = obj.tri.Points;
                newVerts(i,j) = obj.tri.Points(i,j)-pert;
                obj2 = obj.setVerts(newVerts);
                [VortexArr,VortexBrr,VortexAcr,VortexBcr] = obj2.influenceMatrix(obj2,obj.approxMat.calcIndex{i},obj.approxMat.calcIndex{i});

                if obj.approxMat.wakeGradFlag == 1
                    VortexArwr = zeros(numel(obj.approxMat.calcIndex{i}),nbPanel);
                    VortexAcwr = zeros(nbPanel,numel(obj.approxMat.calcIndex{i}));
                    for wakeNo = 1:numel(obj.wakeline)
                        for edgeNo = 1:size(obj.wakeline{wakeNo}.validedge,2)
                            interpID(1) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.upperID(edgeNo));
                            interpID(2) = obj.IndexPanel2Solver(obj.wakeline{wakeNo}.lowerID(edgeNo));
                            if isempty(intersect(interpID,obj.approxMat.calcIndex{i}))
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,obj.approxMat.calcIndex{i},obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence;
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence;
                            else
                                influence = obj2.wakeInfluenceMatrix(obj2,wakeNo,edgeNo,1:nbPanel,obj.wakeline{wakeNo}.wakeShape{edgeNo});
                                VortexArwr(:,interpID(1)) = VortexArwr(:,interpID(1)) - influence(obj.approxMat.calcIndex{i},:);
                                VortexArwr(:,interpID(2)) = VortexArwr(:,interpID(2)) + influence(obj.approxMat.calcIndex{i},:);
                                if not(isempty(intersect(interpID(1),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(1));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) - influence;
                                end
                                if not(isempty(intersect(interpID(2),obj.approxMat.calcIndex{i})))
                                    [~,b] = find(obj.approxMat.calcIndex{i}==interpID(2));
                                    VortexAcwr(:,b) = VortexAcwr(:,b) + influence;
                                end
                            end
                        end
                    end
                end

                obj.approxMat.dVArZ{i} = (VortexArf-VortexArr)./pert/2;
                obj.approxMat.dVBrZ{i} = (VortexBrf-VortexBrr)./pert/2;
                obj.approxMat.dVAcZ{i} = (VortexAcf-VortexAcr)./pert/2;
                obj.approxMat.dVBcZ{i} = (VortexBcf-VortexBcr)./pert/2;

               if obj.approxMat.wakeGradFlag == 1
                    obj.approxMat.dVArwZ{i} = (VortexArwf-VortexArwr)./pert/2;
                    obj.approxMat.dVAcwZ{i} = (VortexAcwf-VortexAcwr)./pert/2;
                end

                
                cutoffVal = 0.0;
                obj.approxMat.dVArX{i}(abs(obj.approxMat.dVArX{i})<cutoffVal) = 0;
                obj.approxMat.dVBrX{i}(abs(obj.approxMat.dVBrX{i})<cutoffVal) = 0;
                obj.approxMat.dVAcX{i}(abs(obj.approxMat.dVAcX{i})<cutoffVal) = 0;
                obj.approxMat.dVBcX{i}(abs(obj.approxMat.dVBcX{i})<cutoffVal) = 0;
                obj.approxMat.dVArY{i}(abs(obj.approxMat.dVArY{i})<cutoffVal) = 0;
                obj.approxMat.dVBrY{i}(abs(obj.approxMat.dVBrY{i})<cutoffVal) = 0;
                obj.approxMat.dVAcY{i}(abs(obj.approxMat.dVAcY{i})<cutoffVal) = 0;
                obj.approxMat.dVBcY{i}(abs(obj.approxMat.dVBcY{i})<cutoffVal) = 0;
                obj.approxMat.dVArZ{i}(abs(obj.approxMat.dVArZ{i})<cutoffVal) = 0;
                obj.approxMat.dVBrZ{i}(abs(obj.approxMat.dVBrZ{i})<cutoffVal) = 0;
                obj.approxMat.dVAcZ{i}(abs(obj.approxMat.dVAcZ{i})<cutoffVal) = 0;
                obj.approxMat.dVBcZ{i}(abs(obj.approxMat.dVBcZ{i})<cutoffVal) = 0;


                if obj.approxMat.wakeGradFlag == 1
                    obj.approxMat.dVArX{i}(abs(obj.approxMat.dVArX{i})<cutoffVal) = 0;
                    obj.approxMat.dVArY{i}(abs(obj.approxMat.dVArY{i})<cutoffVal) = 0;
                    obj.approxMat.dVArZ{i}(abs(obj.approxMat.dVArZ{i})<cutoffVal) = 0;
                    obj.approxMat.dVAcX{i}(abs(obj.approxMat.dVAcX{i})<cutoffVal) = 0;
                    obj.approxMat.dVAcY{i}(abs(obj.approxMat.dVAcY{i})<cutoffVal) = 0;
                    obj.approxMat.dVAcZ{i}(abs(obj.approxMat.dVAcZ{i})<cutoffVal) = 0;

                end
            end
            fprintf("\n");
        end
        

        function approxmatedObj = makeApproximatedInstance(obj,modifiedVerts,skipCpMatandLLT)
            if nargin < 3
                skipCpMatandLLT = 0;
            end
            nPanel = numel(obj.paneltype);
            nbPanel = sum(obj.paneltype == 1);
            approxmatedObj = obj.setVerts(modifiedVerts);
            
            % パネル法連立方程式行列の作成
            if skipCpMatandLLT == 0
                for i = 1:nPanel
                    if approxmatedObj.paneltype(i) == 1
                        CPmat = approxmatedObj.center(approxmatedObj.cluster{i},1:3);
                        pnt = approxmatedObj.center(i,:);
                        m = approxmatedObj.tri.Points(approxmatedObj.tri.ConnectivityList(i,1),:)'-pnt(:);
                        m = m./norm(m);
                        l = cross(m,approxmatedObj.orgNormal(i,:)');
                        Minv = [l,m,approxmatedObj.orgNormal(i,:)'];
                        lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
                        bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
                        Bmat=pinv(bb,sqrt(eps));
                        Vnmat = Minv*[1,0,0,0;0,1,0,0;0,0,1,0;]*Bmat;
                        for iter = 1:3
                            approxmatedObj.mu2v{iter}(i,approxmatedObj.IndexPanel2Solver(obj.cluster{i})) = Vnmat(iter,:);
                        end
                    end
                end
            end



            %approxmatedObj.approxMat = [];
            approxmatedObj.approximated = 1;
            dv = approxmatedObj.tri.Points-obj.tri.Points;
            calcIndex = find(abs(dv(:,1))>sqrt(eps))';
            for i = calcIndex
                if ~isempty(obj.approxMat.calcIndex{i})
                    obj.approxMat.dVAcX{i}(obj.approxMat.calcIndex{i},:) = 0; 
                    obj.approxMat.dVBcX{i}(obj.approxMat.calcIndex{i},:) = 0; 
                    approxmatedObj.LHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.LHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArX{i}.*dv(i,1);
                    approxmatedObj.RHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.RHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVBrX{i}.*dv(i,1);
                    approxmatedObj.LHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.LHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcX{i}.*dv(i,1);
                    approxmatedObj.RHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.RHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVBcX{i}.*dv(i,1);
                end
            end
            calcIndex = find(abs(dv(:,2))>sqrt(eps))';
            for i = calcIndex
                if ~isempty(obj.approxMat.calcIndex{i})
                    obj.approxMat.dVAcY{i}(obj.approxMat.calcIndex{i},:) = 0; 
                    obj.approxMat.dVBcY{i}(obj.approxMat.calcIndex{i},:) = 0;
                    approxmatedObj.LHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.LHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArY{i}.*dv(i,2);
                    approxmatedObj.RHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.RHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVBrY{i}.*dv(i,2);
                    approxmatedObj.LHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.LHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcY{i}.*dv(i,2);
                    approxmatedObj.RHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.RHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVBcY{i}.*dv(i,2);
                end
            end
            calcIndex = find(abs(dv(:,3))>sqrt(eps))';
            for i = calcIndex
                if ~isempty(obj.approxMat.calcIndex{i})
                    obj.approxMat.dVAcZ{i}(obj.approxMat.calcIndex{i},:) = 0; 
                    obj.approxMat.dVBcZ{i}(obj.approxMat.calcIndex{i},:) = 0;
                    approxmatedObj.LHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.LHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArZ{i}.*dv(i,3);
                    approxmatedObj.RHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.RHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVBrZ{i}.*dv(i,3);
                    approxmatedObj.LHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.LHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcZ{i}.*dv(i,3);
                    approxmatedObj.RHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.RHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVBcZ{i}.*dv(i,3);
                end
            end
            if obj.approxMat.wakeGradFlag == 1
                calcIndex = find(abs(dv(:,1))>sqrt(eps))';
                for i = calcIndex
                    obj.approxMat.dVAcwX{i}(obj.approxMat.calcIndex{i},:) = 0;  
                    approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArwX{i}.*dv(i,1);
                    approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcwX{i}.*dv(i,1);
                end
                calcIndex = find(abs(dv(:,2))>sqrt(eps))';
                for i = calcIndex  
                    obj.approxMat.dVAcwY{i}(obj.approxMat.calcIndex{i},:) = 0; 
                    approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArwY{i}.*dv(i,2);
                    approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcwY{i}.*dv(i,2);
                end
                calcIndex = find(abs(dv(:,3))>sqrt(eps))';
                for i = calcIndex
                    obj.approxMat.dVAcwZ{i}(obj.approxMat.calcIndex{i},:) = 0;
                    approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:) = approxmatedObj.wakeLHS(obj.approxMat.calcIndex{i},:)+obj.approxMat.dVArwZ{i}.*dv(i,3);
                    approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i}) = approxmatedObj.wakeLHS(:,obj.approxMat.calcIndex{i})+obj.approxMat.dVAcwZ{i}.*dv(i,3);
                end
            else
                approxmatedObj = approxmatedObj.makeWakeEquation();
            end
            if skipCpMatandLLT == 0
                for wakeNo = 1:numel(approxmatedObj.wakeline)
                    theta = linspace(pi,0,size(approxmatedObj.wakeline{wakeNo}.validedge,2)*approxmatedObj.settingUNLSI.LLTnInterp+1);
                    iter = 1;
                    approxmatedObj.LLT.sp{wakeNo} = [];
                    approxmatedObj.LLT.calcMu{wakeNo} = zeros(1,nbPanel);
                    s = zeros(1,numel(approxmatedObj.wakeline{wakeNo}.edge));
                    jter = 1;
                    for edgeNo = 1:numel(approxmatedObj.wakeline{wakeNo}.edge)-1
                        s(edgeNo+1) = s(edgeNo) + norm(approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(edgeNo),2:3)-approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(edgeNo+1),2:3));
                        if approxmatedObj.wakeline{wakeNo}.valid(edgeNo) == 1
                            approxmatedObj.LLT.sp{wakeNo} =  [approxmatedObj.LLT.sp{wakeNo},(s(edgeNo)+s(edgeNo+1))./2];
                            jter = jter+1;
                        end
                    end
                    for edgeNo = 1:size(approxmatedObj.wakeline{wakeNo}.validedge,2)
                        interpID(1) = approxmatedObj.IndexPanel2Solver(approxmatedObj.wakeline{wakeNo}.upperID(edgeNo));
                        interpID(2) = approxmatedObj.IndexPanel2Solver(approxmatedObj.wakeline{wakeNo}.lowerID(edgeNo));
                        approxmatedObj.LLT.calcMu{wakeNo}(iter,interpID(1)) = -1;
                        approxmatedObj.LLT.calcMu{wakeNo}(iter,interpID(2)) = 1;
                        iter = iter+1;
                    end
                    sd = (s(end)-s(1))*(cos(theta)./2+0.5)+s(1);
                    approxmatedObj.LLT.sinterp{wakeNo} = (sd(2:end)+sd(1:end-1))./2;
                    approxmatedObj.LLT.yinterp{wakeNo} = interp1(s,approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(:),2),approxmatedObj.LLT.sinterp{wakeNo},'linear','extrap');
                    approxmatedObj.LLT.zinterp{wakeNo} = interp1(s,approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(:),3),approxmatedObj.LLT.sinterp{wakeNo},'linear','extrap');
                    yd = interp1(s,approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(:),2),sd,'linear','extrap');
                    zd = interp1(s,approxmatedObj.tri.Points(approxmatedObj.wakeline{wakeNo}.edge(:),3),sd,'linear','extrap');
                    approxmatedObj.LLT.phiinterp{wakeNo} = atan((zd(2:end)-zd(1:end-1))./(yd(2:end)-yd(1:end-1)));
                    approxmatedObj.LLT.spanel{wakeNo} = (sd(2:end)-sd(1:end-1))./2;
    
                end
                if numel(approxmatedObj.wakeline)>0
                    approxmatedObj.LLT.Qij = approxmatedObj.Calc_Q(horzcat(obj.LLT.yinterp{:}),horzcat(obj.LLT.zinterp{:}),horzcat(obj.LLT.phiinterp{:}),horzcat(obj.LLT.spanel{:}),obj.halfmesh);
                end
            end

        end

        function approxmatedObj = makeApproximatedFemInstance(obj,modifiedVerts)
            approxmatedObj = obj.setFemVerts(modifiedVerts);

            %approxmatedObj.approxFemMat = [];
            approxmatedObj.femapproximated = 1;
            dv = approxmatedObj.femtri.Points-obj.femtri.Points;
            calcIndex = find(any(abs(dv)>sqrt(eps),2))';
            for i = calcIndex
                approxmatedObj.femLHS = approxmatedObj.femLHS+obj.approxFemMat.dLHSX{i}.*dv(i,1)+obj.approxFemMat.dLHSY{i}.*dv(i,2)+obj.approxFemMat.dLHSZ{i}.*dv(i,3);
            end
            
        end

        function obj = setProp(obj,propNo,ID,diameter,XZsliced)
            %IDのパネルタイプをプロペラ==4に変更
            for iter = 1:numel(propNo)
                obj = obj.setPanelType(ID(iter),'prop');
                obj.prop{propNo(iter)}.ID = ID(iter);
                nPanel = numel(obj.paneltype);
                %プロペラディスク⇒bodyパネルへの影響係数
                obj.prop{propNo(iter)}.diameter = diameter(iter);
                obj.prop{propNo(iter)}.area = pi*(obj.prop{propNo(iter)}.diameter/2)^2;
                obj.prop{propNo(iter)}.XZsliced  = XZsliced(iter);
                %ペラパネルの重心位置を計算
                mom = [0,0,0];
                propArea = 0;
                for i = 1:nPanel
                    if obj.surfID(i) == ID(iter)
                        darea = obj.vertex(obj.tri.Points(obj.tri(i,1),:),obj.tri.Points(obj.tri(i,2),:),obj.tri.Points(obj.tri(i,3),:));
                        propArea = propArea + darea;
                        mom = mom + darea.*obj.center(i,:);
                        if XZsliced(iter) == 1
                            mom = mom + darea.*[obj.center(i,1),-obj.center(i,2),obj.center(i,3)];
                        end
                    end
                end
                if XZsliced(iter) == 1
                    propArea = 2*propArea;
                end
                obj.prop{propNo(iter)}.normal = mean(obj.orgNormal(obj.surfID == ID(iter),:),1);
                obj.prop{propNo(iter)}.center = mom./propArea;
                obj = obj.setPropState(propNo(iter),sqrt(eps),sqrt(eps),sqrt(eps));
            end
            
            for iter = 1:numel(obj.prop)
                if not(isempty(obj.prop{iter}))
                    obj = makePropEquation(obj,iter);
                end
            end
        end


        function obj = makePropEquation(obj, propNo)
            % makePropEquationメソッドは、指定されたプロペラ番号に対してプロペラの影響行列を計算し、オブジェクトのプロパティに格納します。
            % 
            % 入力:
            %   - obj: UNLSIオブジェクト
            %   - propNo: プロペラ番号
            %
            % 出力:
            %   - obj: 更新されたUNLSIオブジェクト
            %
            % 例:
            %   obj = makePropEquation(obj, 1);
            %
            [VortexAi, VortexBi, VortexAo, VortexBo, propWakeA] = obj.propellerInfluenceMatrix(obj, propNo, obj.settingUNLSI.propWakeLength);
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

        function obj = solveFlow(obj,alpha,beta,Mach,Re,ID)
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
                if nargin > 4
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
                    if nargin > 4
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

                Vinf = zeros(nPanel,3);
                for i = 1:nPanel
                   rvec = obj.center(i,:)'-obj.rotCenter(i,:)';
                   Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.angularVelocity(i,:)./180.*pi,rvec(:)')'))';
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
                    u =  -(obj.LHS+obj.wakeLHS)\RHV;
                    %u = -lsqminnorm((obj.LHS+obj.wakeLHS),RHV,obj.settingUNLSI.lsqminnormTol);
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
                    obj.Cp{flowNo} = obj.valLimiter(obj.Cp{flowNo},obj.CpLimit(:,1),obj.CpLimit(:,2),obj.settingUNLSI.CpSlope);
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
                Cpc =  obj.Cp{flowNo}(:,iterflow);
                Cfec = obj.Cfe{flowNo};
                if nargin > 5
                    Cpc(not(any(obj.surfID == ID,2)),:) = 0;
                    Cfec(not(any(obj.surfID == ID,2)),:) = 0;
                end
                %Cp⇒力への変換
                dCA_p = (-Cpc.*obj.modNormal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-Cpc.*obj.modNormal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-Cpc.*obj.modNormal(:,3)).*obj.area./obj.SREF;
                dCA_f = (+Cfec.*s(:,1)).*obj.area./obj.SREF;
                dCY_f = (+Cfec.*s(:,2)).*obj.area./obj.SREF;
                dCN_f = (+Cfec.*s(:,3)).*obj.area./obj.SREF;
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
                Vinf = zeros(nPanel,3);
                for i = 1:nPanel
                   rvec = obj.center(i,:)'-obj.rotCenter(i,:)';
                   Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.angularVelocity(i,:)./180.*pi,rvec(:)')'))';
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
                    usolve =  -(obj.LHS+obj.wakeLHS)\RHV;
                    %usolve = -lsqminnorm((obj.LHS+obj.wakeLHS),RHV,obj.settingUNLSI.lsqminnormTol);
                    Rsolve = (obj.LHS+obj.wakeLHS)*usolve+RHV;
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

        function [obj,CT,Cp,Cq,efficiency] = solveSteadyProp(obj,propID,rpm,rotAxis,rotOrigin,alpha,beta,Mach,Re,fig,caxis)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %プロペラ計算を行う
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:size(propID,1)
                obj = obj.setRotation(propID(i,:),rotOrigin(i,:),[-rpm(i,:)*pi/30*180/pi,0,0]/obj.settingUNLSI.Vinf);
                obj = obj.setHelixWake(propID(i,:),0.02,rpm(i,:),rotAxis(i,:),rotOrigin(i,:));
            end
            obj = obj.solveFlow(alpha,beta,Mach,Re,propID(:)');%パネル法を解く
            data = obj.getAERODATA(alpha,beta,Mach,Re);
            if nargin > 9
                if nargin == 10
                    caxis = [-20,1];
                end
                obj.plotGeometry(fig,obj.getCp(alpha,beta,Mach,Re),caxis);%圧力係数のプロット
                obj.plotWakeShape(fig);
            end
            %TODO それぞれのペラについて求める
            CT = -data(15) * 0.5 * (obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)) ^2 * obj.SREF / ((rpm(1)/60)^2*obj.BREF^4);
            Cq = data(18) * 0.5 * (obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)) ^2 * obj.SREF * obj.BREF / ((rpm(1)/60)^2*obj.BREF^5);
            Cp = Cq * 2 * pi;
            J = obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)/(rpm/60)/obj.BREF;
            efficiency = CT*J/Cp;
        end
        
        function [obj,CT,Cp,Cq,efficiency] = solveUnsteadyProp(obj,propID,dt,rpm,rotAxis,rotOrigin,alpha,beta,Mach,Re,fig,caxis)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %プロペラ計算を行う
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj = obj.setRotation(propID,rotOrigin,[-rpm*pi/30*180/pi,0,0]/obj.settingUNLSI.Vinf);
            obj = obj.solveFlow(alpha,beta,Mach,Re,propID);%パネル法を解く
            data = obj.getAERODATA(alpha,beta,Mach,Re);
            obj = obj.marchWake(dt,alpha,beta,Mach,Re);
            if nargin > 10
                obj.plotGeometry(fig,obj.getCp(alpha,beta,Mach,Re),caxis);%圧力係数のプロット
                obj.plotWakeShape(fig);
            end
            
            rotAngle = -rpm*pi/30 * dt * 180 /pi;
            obj = obj.rotateVerts(propID,rotAngle,rotAxis,rotOrigin);

            CT = -data(15) * 0.5 * (obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)) ^2 * obj.SREF / ((rpm/60)^2*obj.BREF^4);
            Cq = data(18) * 0.5 * (obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)) ^2 * obj.SREF * obj.BREF / ((rpm/60)^2*obj.BREF^5);
            Cp = Cq * 2 * pi;
            J = obj.settingUNLSI.Vinf*cosd(alpha)*cosd(beta)/(rpm/60)/obj.BREF;
            efficiency = CT*J/Cp;
        end
        
        function obj = rotateVerts(obj,ID,rotAngle,rotAxis,rotOrigin)
            rotAxis = rotAxis ./ norm(rotAxis);
            K = [0, -rotAxis(3), rotAxis(2);
                 rotAxis(3), 0, -rotAxis(1);
                -rotAxis(2), rotAxis(1), 0];
            dcm = eye(3) + sind(rotAngle) * K + (1 - cosd(rotAngle)) * K^2;
            usedVerts = unique(obj.tri.ConnectivityList(any(obj.surfID == ID,2),:));
            p = obj.tri.Points';
            p(:,usedVerts) = dcm * (obj.tri.Points(usedVerts,:)'-rotOrigin(:))+rotOrigin(:);
            obj = obj.setVerts(p');
        end

        function [vel] = calcVelocity(obj,u,VmuWake,VmuBody)
            %%%%%%%%%%%%%LSIの求解%%%%%%%%%%%%%%%%%%%%%
            %ポテンシャルから力を求める
            %結果は配列に出力される
            % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt 7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14:E(翼効率) 15:CFx 16:CFy 17:CFz 18:CMx 19:CMy 20:CMz 21:CMl 22:CMm 23:CMn 24:FOpt 
            %上記で求めていないものは0が代入される
            %u: doubletの強さ
            %flowNo:解きたい流れのID
            %alpha:迎角[deg]
            %beta:横滑り角[deg]
            %omega:主流の回転角速度(deg/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vel(:,1) = -((VmuBody.X+VmuWake.X)*u)./4./pi;
            vel(:,2) = -((VmuBody.Y+VmuWake.Y)*u)./4./pi;
            vel(:,3) = -((VmuBody.Z+VmuWake.Z)*u)./4./pi;
        end
        
        function [AERODATA,Cp,Cfe,R,obj] = solveFlowForAdjoint(obj,u,flowNo,alpha,beta,ID)
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
                Vinf = zeros(nPanel,3);
                for i = 1:nPanel
                   rvec = obj.center(i,:)'-obj.rotCenter(i,:)';
                   Vinf(i,:) = (reshape(obj.deflAngle(find(obj.deflAngle(:,1)==obj.surfID(i,1)),6:14),[3,3])'*(T*[1;0;0])-(cross(obj.angularVelocity(i,:)./180.*pi,rvec(:)')'))';
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
                    obj.Cp{flowNo} = obj.valLimiter(obj.Cp{flowNo},obj.CpLimit(:,1),obj.CpLimit(:,2),obj.settingUNLSI.CpSlope);
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
                Cpc =  obj.Cp{flowNo}(:,iterflow);
                Cfec = obj.Cfe{flowNo};
                if nargin > 6
                    Cpc(not(any(obj.surfID == ID,2)),:) = 0;
                    Cfec(not(any(obj.surfID == ID,2)),:) = 0;
                end
                %Cp⇒力への変換
                dCA_p = (-Cpc.*obj.modNormal(:,1)).*obj.area./obj.SREF;
                dCY_p = (-Cpc.*obj.modNormal(:,2)).*obj.area./obj.SREF;
                dCN_p = (-Cpc.*obj.modNormal(:,3)).*obj.area./obj.SREF;
                dCA_f = (+Cfec.*s(:,1)).*obj.area./obj.SREF;
                dCY_f = (+Cfec.*s(:,2)).*obj.area./obj.SREF;
                dCN_f = (+Cfec.*s(:,3)).*obj.area./obj.SREF;
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
                    R(nbPanel*(iterflow-1)+1:nbPanel*iterflow,1) = (obj.LHS+obj.wakeLHS)*usolve+RHV;
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

        %動安定微係数推算
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
	            Yb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF/(2*mass)*dynCoefStruct{i}.Cyb;
	            Lb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF*obj.BREF/(2*Inatia(1,1))*dynCoefStruct{i}.Clb;
	            Nb = obj.settingUNLSI.rho*UREF*UREF*obj.SREF*obj.BREF/(2*Inatia(3,3))*dynCoefStruct{i}.Cnb;
	            Lp = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(1,1))*dynCoefStruct{i}.Clp;
	            Np = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(3,3))*dynCoefStruct{i}.Cnp;
	            Lr = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(1,1))*dynCoefStruct{i}.Clr;
	            Nr = obj.settingUNLSI.rho*UREF*obj.SREF*obj.BREF*obj.BREF/(4*Inatia(3,3))*dynCoefStruct{i}.Cnr;
                Za = obj.settingUNLSI.rho*UREF^2*obj.SREF/(2*mass)*dynCoefStruct{i}.Cza;
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
                    udw(1+nbPanel*(iterflow-1):nbPanel*iterflow,1) = -(obj.LHS+obj.wakeLHS)\rhsdw;
        
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
                    udw(1+nbPanel*(iterflow-1):nbPanel*iterflow,2) = -(obj.LHS+obj.wakeLHS)\rhsdw;
                end

                for jter = 1:3%p,q,rについて差分をとる
                    %if isempty(obj.settingUNLSI.angularVelocity)
                        omegadw = zeros(1,3);
                    %end
                    %else
                        %omegadw = obj.settingUNLSI.angularVelocity(iterflow,:);
                    %end
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
                    buff = -(obj.LHS+obj.wakeLHS)\rhsdw;
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
                    avorg = obj.angularVelocity;
                    for jter = 1:3%p,q,rについて差分をとる
                        omegadwf = zeros(1,3);
                        omegadwr = zeros(1,3);
                        omegadwf(jter) = omegadwf(jter)+dw/pi*180;
                        omegadwr(jter) = omegadwr(jter)-dw/pi*180;
                        obj.angularVelocity = avorg + repmat(omegadwf,[size(avorg,1),1]);
                        uf = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                        obj.angularVelocity = avorg + repmat(omegadwr,[size(avorg,1),1]);
                        ur = obj.solvePertPotential(flowNo,alpha(iterflow),beta(iterflow));
                        obj.settingUNLSI.angularVelocity = avorg;
                        udwf(1+nbPanel*(iterflow-1):nbPanel*iterflow,2+jter) = uf-u0;
                        udwr(1+nbPanel*(iterflow-1):nbPanel*iterflow,2+jter) = u0-ur;
                    end
                end
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
                avorg = obj.angularVelocity;
                for jter = 1:3%p,q,rについて差分をとる
                    omegadwf = zeros(1,3);
                    omegadwr = zeros(1,3);
                    omegadwf(jter) = omegadwf(jter)+dw.*180./pi;
                    omegadwr(jter) = omegadwr(jter)-dw.*180./pi;
                    obj.angularVelocity = avorg + repmat(omegadwf,[size(avorg,1),1]);
                    AERODATAf = obj.solveFlowForAdjoint(u0+udwf(:,2+jter),flowNo,alpha,beta);
                    obj.angularVelocity = avorg + repmat(omegadwr,[size(avorg,1),1]);
                    AERODATAr = obj.solveFlowForAdjoint(u0-udwr(:,2+jter),flowNo,alpha,beta);
                    obj.angularVelocity = avorg;
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
                for i = 1:6
                    for j = 1:size(dyndata(:,7:end),2)/6
                        obj.ppDyn{i,j+1} = griddedInterpolant(X1,X2,X3,dyndata2{6*(j-1)+i},"linear","linear");
                    end
                end
                obj.ppDyn{1,1} = betarange2;
                obj.ppDyn{2,1} = machrange2;
                obj.ppDyn{3,1} = alpharange2;
                obj.ppDyn{4,1} = Re;
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
                for i = 1:6
                    for j = 1:size(dyndata(:,7:end),2)/6
                        obj.ppDyn{i,j+1} = griddedInterpolant(X1,X2,dyndata2{6*(j-1)+i},"linear","linear");
                    end
                end
                obj.ppDyn{1,1} = betarange2;
                obj.ppDyn{2,1} = MachRange;
                obj.ppDyn{3,1} = alpharange2;
                obj.ppDyn{4,1} = Re;
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
                if strcmpi(dispCoeforDyn,"coef")
                    switch dispIndex
                        case 1
                            zlabel("CD");
                        case 2
                            zlabel("CY");
                        case 3
                            zlabel("CL");
                        case 4
                            zlabel("Cl");
                        case 5
                            zlabel("Cm");
                        case 6
                            zlabel("Cn");
                        otherwise
                            error("invalid INDEX");
                    end
                else
                    if dispIndex<7
                        error("invalid dispIndex");
                    end
                    switch floor(dispIndex/6)
                        case 1
                            derivstr = "b";
                        case 2
                            derivstr = "a";
                        case 3
                            derivstr = "p";
                        case 4
                            derivstr = "q";
                        case 5
                            derivstr = "r";
                        otherwise
                            derivstr = strcat("f",num2str(floor(dispIndex/6)-5));
                    end
                    switch mod(dispIndex,6)
                        case 1
                            coefstr = "CX";
                        case 2
                            coefstr = "CY";
                        case 3
                            coefstr = "CZ";
                        case 4
                            coefstr = "Cl";
                        case 5
                            coefstr = "Cm";
                        case 0
                            coefstr = "Cn";
                    end
                    zlabel(coefstr+derivstr);
                end
                
            end
        end
    
        function plotSurrogateModel(obj,ppCoef,ppDyn,alphaRange,betaRange,Mach,figureNo)
            for dispIndex = 1:numel(ppDyn)
                testData = [];
                if numel(ppCoef{1}.GridVectors) == 2
                    [x1,x2] = ndgrid(alphaRange,betaRange);
                    for a = 1:numel(alphaRange)
                        for b = 1:numel(betaRange)
                            if dispIndex<7
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
                             if dispIndex<7
                                testData(a,b) = ppCoef{dispIndex}(x1(a,b),x2(a,b),Mach);
                            else
                                testData(a,b) = ppDyn{dispIndex}(x1(a,b),x2(a,b),Mach);
                            end
                        end
                    end
                end
                figure(figureNo+dispIndex-1);clf;
                mesh(x1,x2,testData);
                xlabel("alpha(deg)");
                ylabel("beta(deg)");
                if dispIndex < 7
                    switch dispIndex
                        case 1
                            zlabel("CD");
                        case 2
                            zlabel("CY");
                        case 3
                            zlabel("CL");
                        case 4
                            zlabel("Cl");
                        case 5
                            zlabel("Cm");
                        case 6
                            zlabel("Cn");
                        otherwise
                            error("invalid INDEX");
                    end
                else
                    switch floor((dispIndex-1)/6)
                        case 1
                            derivstr = "b";
                        case 2
                            derivstr = "a";
                        case 3
                            derivstr = "p";
                        case 4
                            derivstr = "q";
                        case 5
                            derivstr = "r";
                        otherwise
                            derivstr = strcat("f",num2str(floor((dispIndex-1)/6)-5));
                    end
                    switch mod(dispIndex,6)
                        case 1
                            coefstr = "CX";
                        case 2
                            coefstr = "CY";
                        case 3
                            coefstr = "CZ";
                        case 4
                            coefstr = "Cl";
                        case 5
                            coefstr = "Cm";
                        case 0
                            coefstr = "Cn";
                    end
                    zlabel(coefstr+derivstr);
                end
            end
        end



        %FEM関連
        function obj = setFemMesh(obj,verts,connectivity,femID,aeroIDLink)
            obj.femutils = [];
            if nargin<5
                obj.femutils.aeroIDLink = unique(obj.surfID);
            else
                obj.femutils.aeroIDLink = aeroIDLink;
            end
            [verts,connectivity] = obj.mergeVerts(0.004,verts,connectivity);
            obj.femtri = triangulation(connectivity,verts);
            obj.femID = femID;

            obj.femcenter = (verts(obj.femtri.ConnectivityList(:,1),:)+verts(obj.femtri.ConnectivityList(:,2),:)+verts(obj.femtri.ConnectivityList(:,3),:))./3;
            obj.femarea = 0.5 * sqrt(sum(cross(obj.femtri.Points(obj.femtri.ConnectivityList(:, 2), :) - obj.femtri.Points(obj.femtri.ConnectivityList(:, 1), :), obj.femtri.Points(obj.femtri.ConnectivityList(:, 3), :) - obj.femtri.Points(obj.femtri.ConnectivityList(:, 1), :), 2).^2, 2));
            obj.femNormal = faceNormal(obj.femtri);

            obj = obj.checkFemMesh(obj.settingUNLSI.checkMeshTol,"delete");
                
            %取り急ぎ、解析用のメッシュと空力メッシュは同一。後でFEAメッシュ等からメッシュを入れられるようにする。
            usedVertsbuff = unique(obj.femtri.ConnectivityList);
            obj.femutils.usedVerts=usedVertsbuff(:)';
            obj.femutils.nbVerts = numel(obj.femutils.usedVerts);
            obj.femThn = ones(size(obj.femID)).*0.01;
            obj.femE = ones(size(obj.femID)).*100000000000;
            obj.femRho = ones(size(obj.femID)).*1000;
            obj.femVisc = ones(size(obj.femID)).*0;
            for iter = 1:size(obj.femtri.ConnectivityList,1)
                obj.femutils.IndexRow{iter} = zeros(18,18);
                obj.femutils.IndexCol{iter} = zeros(18,18);
                for i = 1:6
                    for j = 1:3
                        [r,c] = find(obj.femutils.usedVerts==obj.femtri.ConnectivityList(iter,j));
                        obj.femutils.IndexRow{iter}(3*(i-1)+j,:) = c+obj.femutils.nbVerts*(i-1);
                        obj.femutils.IndexCol{iter}(:,3*(i-1)+j) = c+obj.femutils.nbVerts*(i-1);
                    end
                end
            end
            %境界条件設定
            obj.femutils.InvMatIndex = [];
            obj.femutils.MatIndex = zeros(6*obj.femutils.nbVerts,1);
            for i = 1:numel(obj.femutils.usedVerts)
               if  abs(obj.femtri.Points(obj.femutils.usedVerts(i),2))<=0.01%0.01
                   obj.femutils.MatIndex(i,1) = 0;
                   obj.femutils.MatIndex(1*obj.femutils.nbVerts+i,1) = 0;
                   obj.femutils.MatIndex(2*obj.femutils.nbVerts+i,1) = 0;
                   obj.femutils.MatIndex(3*obj.femutils.nbVerts+i,1) = 0;
                   obj.femutils.MatIndex(4*obj.femutils.nbVerts+i,1) = 0;
                   obj.femutils.MatIndex(5*obj.femutils.nbVerts+i,1) = 0;
               else
                   obj.femutils.MatIndex(i,1) = 1;
                   obj.femutils.MatIndex(1*obj.femutils.nbVerts+i,1) = 1;
                   obj.femutils.MatIndex(2*obj.femutils.nbVerts+i,1) = 1;
                   obj.femutils.MatIndex(3*obj.femutils.nbVerts+i,1) = 1;%0930修正
                   obj.femutils.MatIndex(4*obj.femutils.nbVerts+i,1) = 1;%0930修正
                   obj.femutils.MatIndex(5*obj.femutils.nbVerts+i,1) = 0;
                   obj.femutils.InvMatIndex=[obj.femutils.InvMatIndex,i];
                %    disp(size([obj.femutils.InvMatIndex,i]))
               end
            end
            obj.femutils.InvMatIndex=[obj.femutils.InvMatIndex,1*obj.femutils.nbVerts+obj.femutils.InvMatIndex,2*obj.femutils.nbVerts+obj.femutils.InvMatIndex,3*obj.femutils.nbVerts+obj.femutils.InvMatIndex,4*obj.femutils.nbVerts+obj.femutils.InvMatIndex]';
            %obj.femutils.InvMatIndex=[obj.femutils.InvMatIndex,1*obj.femutils.nbVerts+obj.femutils.InvMatIndex,2*obj.femutils.nbVerts+obj.femutils.InvMatIndex]';
            %obj.femutils.InvMatIndex=[obj.femutils.InvMatIndex,1*obj.femutils.nbVerts+obj.femutils.InvMatIndex,2*obj.femutils.nbVerts+obj.femutils.InvMatIndex,3*obj.femutils.nbVerts+obj.femutils.InvMatIndex,4*obj.femutils.nbVerts+obj.femutils.InvMatIndex,5*obj.femutils.nbVerts+obj.femutils.InvMatIndex]';
            
            %構造メッシュと空力メッシュの変換行列の作成
            %Rstrを計算

            usedAerobuff = unique(obj.tri.ConnectivityList(any(obj.surfID == obj.femutils.aeroIDLink(:)',2),:));
            obj.femutils.usedAeroVerts =  usedAerobuff(:);
            femMeshScaling = abs(max(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1)-min(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1));
            %femMeshScaling = [1,1,1];
            femMeshScaling(femMeshScaling==0) = 1;
            Rss = obj.calcRMat(obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling,obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling);
            phi = @(r)obj.phi1(r,0.001);%変位計算の場合はphi1の方がうまくいく
            Mss = phi(Rss);
            Pss = [ones(size(obj.femtri.Points(obj.femutils.usedVerts,:),1),1),obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling ]';
            Ras = obj.calcRMat(obj.tri.Points(obj.femutils.usedAeroVerts,:)./femMeshScaling,obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling);
            Aas = [ones(size(obj.tri.Points(obj.femutils.usedAeroVerts,:),1),1),obj.tri.Points(obj.femutils.usedAeroVerts,:)./femMeshScaling,phi(Ras)];
            invM = pinv(Mss);
            Mp = pinv(Pss*invM*Pss');
            obj.fem2aeroMat = Aas * [Mp*Pss*invM;invM-invM*Pss'*Mp*Pss*invM];
            

            %femのパネル中心座標からメッシュノードへ投影する変換行列の作成
            phi = @(r)obj.phi3(r,1);
            Avv = phi(obj.calcRMat(obj.femtri.Points,obj.femtri.Points));
            Acv = phi(obj.calcRMat(obj.femcenter,obj.femtri.Points));
            obj.femverts2centerMat = Acv*pinv(Avv);

            %基準メッシュから構造メッシュへの投影UNGRADE用
            phi = @(r)obj.phi3(r,1);%変位計算の場合はphi1の方がうまくいく
            meshScaling = abs(max(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1)-min(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1));
            meshScaling(meshScaling==0) = 1;
            Agg = phi(obj.calcRMat(obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling,obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling));
            Amg = phi(obj.calcRMat(obj.femtri.Points(obj.femutils.usedVerts,:)./meshScaling,obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling));
            obj.aero2femMat = Amg*pinv(Agg);

        end

        function obj = setFemVerts(obj,modVerts,skipTransMat)
            if nargin < 3
                skipTransMat = 0;
            end
            obj.femtri = triangulation(obj.femtri.ConnectivityList,modVerts);

            obj.femcenter = (modVerts(obj.femtri.ConnectivityList(:,1),:)+modVerts(obj.femtri.ConnectivityList(:,2),:)+modVerts(obj.femtri.ConnectivityList(:,3),:))./3;
            obj.femarea = 0.5 * sqrt(sum(cross(obj.femtri.Points(obj.femtri.ConnectivityList(:, 2), :) - obj.femtri.Points(obj.femtri.ConnectivityList(:, 1), :), obj.femtri.Points(obj.femtri.ConnectivityList(:, 3), :) - obj.femtri.Points(obj.femtri.ConnectivityList(:, 1), :), 2).^2, 2));
            obj.femNormal = faceNormal(obj.femtri);

            if skipTransMat == 0
                %構造メッシュと空力メッシュの変換行列の作成
                %Rstrを計算
    
                usedAerobuff = unique(obj.tri.ConnectivityList(any(obj.surfID == obj.femutils.aeroIDLink(:)',2),:));
                obj.femutils.usedAeroVerts =  usedAerobuff(:);
                femMeshScaling = abs(max(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1)-min(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1));
                %femMeshScaling = [1,1,1];
                femMeshScaling(femMeshScaling==0) = 1;
                Rss = obj.calcRMat(obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling,obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling);
                phi = @(r)obj.phi1(r,0.001);%変位計算の場合はphi1の方がうまくいく
                Mss = phi(Rss);
                Pss = [ones(size(obj.femtri.Points(obj.femutils.usedVerts,:),1),1),obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling ]';
                Ras = obj.calcRMat(obj.tri.Points(obj.femutils.usedAeroVerts,:)./femMeshScaling,obj.femtri.Points(obj.femutils.usedVerts,:)./femMeshScaling);
                Aas = [ones(size(obj.tri.Points(obj.femutils.usedAeroVerts,:),1),1),obj.tri.Points(obj.femutils.usedAeroVerts,:)./femMeshScaling,phi(Ras)];
                invM = pinv(Mss);
                Mp = pinv(Pss*invM*Pss');
                obj.fem2aeroMat = Aas * [Mp*Pss*invM;invM-invM*Pss'*Mp*Pss*invM];
    
                %femのパネル中心座標からメッシュノードへ投影する変換行列の作成
                phi = @(r)obj.phi3(r,1);
                Avv = phi(obj.calcRMat(obj.femtri.Points,obj.femtri.Points));
                Acv = phi(obj.calcRMat(obj.femcenter,obj.femtri.Points));
                obj.femverts2centerMat = Acv*pinv(Avv);
    
                %基準メッシュから構造メッシュへの投影UNGRADE用
                meshScaling = abs(max(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1)-min(obj.tri.Points(obj.femutils.usedAeroVerts,:),[],1));
                meshScaling(meshScaling==0) = 1;
                Agg = phi(obj.calcRMat(obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling,obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling));
                Amg = phi(obj.calcRMat(obj.femtri.Points(obj.femutils.usedVerts,:)./meshScaling,obj.tri.Points(obj.femutils.usedAeroVerts,:)./meshScaling));
                obj.aero2femMat = Amg*pinv(Agg);
            end

        end
        function plotFemMesh(obj,figID,delta,plotID)
            if not(isa(figID, 'matlab.graphics.axis.Axes'))
                figure(figID);clf;
                ax = axes(gcf);
            else
                ax = figID;
            end
            if nargin<4
                plotID = unique(obj.femID);
            end
            if isempty(delta)
                trisurf(obj.femtri.ConnectivityList(any(obj.femID == plotID(:)',2),:),obj.femtri.Points(:,1),obj.femtri.Points(:,2),obj.femtri.Points(:,3),"FaceAlpha",0,'Parent',ax);
                if obj.halfmesh == 1
                     hold(ax,"on");
                     trisurf(obj.femtri.ConnectivityList(any(obj.femID == plotID(:)',2),:),obj.femtri.Points(:,1),-obj.femtri.Points(:,2),obj.femtri.Points(:,3),"FaceAlpha",0,'Parent',ax);
                end
                axis(ax,"equal");grid(ax,"on");hold(ax,"off");

            else
                modVerts = obj.femtri.Points + delta(:,1:3);
                caxis = vecnorm(delta(:,1:3),2,2);
                trisurf(obj.femtri.ConnectivityList(any(obj.femID == plotID(:)',2),:),modVerts(:,1),modVerts(:,2),modVerts(:,3),caxis,'FaceColor','interp','EdgeAlpha',0.15,'Parent',ax);
                if obj.halfmesh == 1
                     hold(ax,"on");
                     trisurf(obj.femtri.ConnectivityList(any(obj.femID == plotID(:)',2),:),modVerts(:,1),-modVerts(:,2),modVerts(:,3),caxis,'FaceColor','interp','EdgeAlpha',0.15,'Parent',ax);
                end
                axis(ax,"equal");grid(ax,"on");hold(ax,"off");colormap(ax,"jet");
            end
            axis(ax,"equal");xlabel(ax,"x");ylabel(ax,"y");zlabel(ax,"z");
            drawnow();pause(0.1);
        end

        function [obj,weight] = setFemMaterials(obj,femID,thickness,youngModulus,rho,viscous)
            for i = 1:numel(femID)
                obj.femThn(obj.femID==femID(i),1) = thickness(i);
                obj.femE(obj.femID==femID(i),1) = youngModulus(i);
                obj.femRho(obj.femID==femID(i),1) = rho(i);
            end
            if nargout>1
                weight = zeros(1,numel(femID));
                for i = 1:numel(femID)
                    weight(i) = sum(obj.femRho(obj.femID==femID(i),1).*obj.femarea(obj.femID==femID(i),1).*obj.femThn(obj.femID==femID(i),1));
                end
            end
            if nargin>5
                for i = 1:numel(femID)
                    obj.femVisc(obj.femID==femID(i),1) = viscous(i);
                end
            end
        end

        

        

        function obj = makeFemEquation(obj,calcID)
            if nargin<2
                calcID = 1:size(obj.femtri.ConnectivityList,1);
                skipFlag = 0;
            else
                skipFlag = 1;
            end
            nu = 0.34;
            qps = [1/6,1/6;2/3,1/6;1/6,2/3];
            obj.femLHS = sparse(6*numel(obj.femutils.usedVerts),6*numel(obj.femutils.usedVerts));
            obj.femMass = sparse(6*numel(obj.femutils.usedVerts),6*numel(obj.femutils.usedVerts));
            obj.femDamp = sparse(6*numel(obj.femutils.usedVerts),6*numel(obj.femutils.usedVerts));
            for iter = calcID
                Dp(1,1) = 1.0; Dp(1,2) = nu;
                Dp(2,1) = nu;  Dp(2,2) = 1.0;
                Dp(3,3) = (1.0-nu)/2.0;
                Dm = Dp;
                Dm = Dm.*obj.femE(iter,1)/(1.0-nu*nu);
                Dp = Dp.*(obj.femE(iter,1)*obj.femThn(iter,1)^3/(12.0*(1.0-nu*nu)));
                U = obj.femtri.Points(obj.femtri.ConnectivityList(iter,2),:)-obj.femtri.Points(obj.femtri.ConnectivityList(iter,1),:);
                V = obj.femtri.Points(obj.femtri.ConnectivityList(iter,3),:)-obj.femtri.Points(obj.femtri.ConnectivityList(iter,1),:);
                transUV(:,1) = U;
                transUV(:,2) = V;
                W = cross(U,V);
                trafo(1,:) = U./norm(U);
                trafo(3,:) = W./norm(W);
                trafo(2,:) = cross(trafo(3,:),trafo(1,:));
                transUV = trafo*transUV;
                dphi(1,1) = -transUV(1,1); % x12 = x1-x2 = 0-x2 = -x2
                dphi(2,1) =  transUV(1,2); % x31 = x3-x1 = x3-0 = x3
                dphi(3,1) =  transUV(1,1)-transUV(1,2); % x23 = x2-x3
                dphi(1,2) = -transUV(2,1); % y12 = y1-y2 = -y2 = 0 (stays zero, as node B and A lie on local x-axis and therefore)
                dphi(2,2) =  transUV(2,2); % y31 = y3-y1 = y3-0 = y3
                dphi(3,2) =  transUV(2,1)-transUV(2,2); % y23 = y2-y3 = 0-y3 = -y3

                B_m(1,1) =  dphi(3,2); %  y23
                B_m(1,3) =  dphi(2,2); %  y31
                B_m(1,5) =  dphi(1,2); %  y12
                B_m(2,2) = -dphi(3,1); % -x23
                B_m(2,4) = -dphi(2,1); % -x31
                B_m(2,6) = -dphi(1,1); % -x12
                B_m(3,1) = -dphi(3,1); % -x23
                B_m(3,2) =  dphi(3,2); %  y23
                B_m(3,3) = -dphi(2,1); % -x31
                B_m(3,4) =  dphi(2,2); %  y31
                B_m(3,5) = -dphi(1,1); % -x12
                B_m(3,6) =  dphi(1,2); %  y12
                B_m = B_m.*1.0/(2.0*(obj.femarea(iter,1)));
                Ke_m = B_m'*Dm*B_m.*obj.femThn(iter,1)*obj.femarea(iter,1);
                M_m= diag(ones(1,6).*0.5);
                M_m(1,3) = 0.25;
                M_m(1,5) = 0.25;
                M_m(2,4) = 0.25;
                M_m(2,6) = 0.25;
                M_m(3,1) = 0.25;
                M_m(3,5) = 0.25;
                M_m(4,2) = 0.25;
                M_m(4,6) = 0.25;
                M_m(5,1) = 0.25;
                M_m(5,3) = 0.25;
                M_m(6,2) = 0.25;
                M_m(6,4) = 0.25;
                M_m = M_m .* obj.femarea(iter,1).*obj.femRho(iter,1).*obj.femThn(iter,1)./3;

                
                sidelen(1) = dphi(1,1)^2+dphi(1,2)^2;
                sidelen(2) = dphi(2,1)^2+dphi(2,2)^2;
                sidelen(3) = dphi(3,1)^2+dphi(3,2)^2;
                Y(1,1) = dphi(3,2)^2.0;
                Y(1,2) = dphi(2,2)^2.0;
                Y(1,3) = dphi(3,2)*dphi(2,2);
                Y(2,1) = dphi(3,1)^2.0;
                Y(2,2) = dphi(2,1)^2.0;
                Y(2,3) = dphi(2,1)*dphi(3,1);
                Y(3,1) = -2.0*dphi(3,1)*dphi(3,2);
                Y(3,2) = -2.0*dphi(2,1)*dphi(2,1);
                Y(3,3) = -dphi(3,1)*dphi(2,2)-dphi(2,1)*dphi(3,2);
                Y = Y.* 1.0/(4.0*obj.femarea(iter,1)^2.0);
                Ainv = zeros(9,9);
                x12 = -transUV(1,1);
                y12 = -transUV(2,1);
                x31 = transUV(1,2);
                y31 = transUV(2,2);
                x23 = transUV(1,1)-transUV(1,2); 
                y23 =  transUV(2,1)-transUV(2,2);
                Ainv(1,1) = 1;
                Ainv(2,4) = 1;
                Ainv(3,7) = 1;
                Ainv(4,1) = -1;
                Ainv(4,4) = 1;
                Ainv(4,5) = y12;
                Ainv(4,6) = -x12;
                Ainv(5,4) = -1;
                Ainv(5,7) = 1;
                Ainv(5,8) = y23;
                Ainv(5,9) = -x23;
                Ainv(6,1) = 1;
                Ainv(6,2) = y31;
                Ainv(6,3) = -x31;
                Ainv(6,7) = -1;
                Ainv(7,1) = 2;
                Ainv(7,2) = -y12;
                Ainv(7,3) = x12;
                Ainv(7,4) = -2;
                Ainv(7,5) = -y12;
                Ainv(7,6) = x12;
                Ainv(8,4) = 2;
                Ainv(8,5) = -y23;
                Ainv(8,6) = x23;
                Ainv(8,7) = -2;
                Ainv(8,8) = -y23;
                Ainv(8,9) = x23;
                Ainv(9,1) = -2;
                Ainv(9,2) = -y31;
                Ainv(9,3) =  x31;
                Ainv(9,7) = 2;
                Ainv(9,8) = -y31;
                Ainv(9,9) = x31;
                Ke_p = zeros(9,9);
                M_p = zeros(9,9);
                M_m2 = zeros(6,6);
                for j = 1:3
                    B = obj.evalBTri(sidelen, qps(j,1), qps(j,2), dphi);
                    N = obj.evalNTri(sidelen, qps(j,1), qps(j,2), Ainv);
                    Nm = obj.evalNmTri(sidelen, qps(j,1), qps(j,2), Ainv);
                    temp = B'*Y'*Dp*Y*B./6;%temp1-3 101式
                    temp2 = N*N'./6;
                    temp3 = Nm'*Nm./6;
                    Ke_p = Ke_p+temp;
                    M_p = M_p+temp2;
                    M_m2 = M_m2+temp3;
                end
                Ke_p = Ke_p.*2.*obj.femarea(iter,1);
                M_p = M_p*obj.femRho(iter,1)*obj.femThn(iter,1)*2.*obj.femarea(iter,1);
                K_out = zeros(18,18);
                M_out = zeros(18,18);
                for i=0:2%%%%要素剛性行列と要素質量行列の作成
                    for j = 0:2 
                        K_out(  6*i+1,    6*j+1)   = Ke_m(2*i+1,  2*j+1);   % uu
                        K_out(  6*i+1,    6*j+1+1) = Ke_m(2*i+1,  2*j+1+1); % uv
                        K_out(  6*i+1+1,  6*j+1)   = Ke_m(2*i+1+1,2*j+1);   % vu
                        K_out(  6*i+1+1,  6*j+1+1) = Ke_m(2*i+1+1,2*j+1+1); % vv
                        K_out(2+6*i+1,  2+6*j+1)   = Ke_p(3*i+1,  3*j+1);   % ww
                        K_out(2+6*i+1,  2+6*j+1+1) = Ke_p(3*i+1,  3*j+1+1); % wx
                        K_out(2+6*i+1,  2+6*j+2+1) = Ke_p(3*i+1,  3*j+2+1); % wy
                        K_out(2+6*i+1+1,2+6*j+1)   = Ke_p(3*i+1+1,3*j+1);   % xw
                        K_out(2+6*i+1+1,2+6*j+1+1) = Ke_p(3*i+1+1,3*j+1+1); % xx
                        K_out(2+6*i+1+1,2+6*j+2+1) = Ke_p(3*i+1+1,3*j+2+1); % xy
                        K_out(2+6*i+2+1,2+6*j+1)   = Ke_p(3*i+2+1,3*j+1);   % yw
                        K_out(2+6*i+2+1,2+6*j+1+1) = Ke_p(3*i+2+1,3*j+1+1); % yx
                        K_out(2+6*i+2+1,2+6*j+2+1) = Ke_p(3*i+2+1,3*j+2+1); % yy

                        M_out(  6*i+1,    6*j+1)   = M_m(2*i+1,  2*j+1);   % uu
                        M_out(  6*i+1,    6*j+1+1) = M_m(2*i+1,  2*j+1+1); % uv
                        M_out(  6*i+1+1,  6*j+1)   = M_m(2*i+1+1,2*j+1);   % vu
                        M_out(  6*i+1+1,  6*j+1+1) = M_m(2*i+1+1,2*j+1+1); % vv
                        M_out(2+6*i+1,  2+6*j+1)   = M_p(3*i+1,  3*j+1);   % ww
                        M_out(2+6*i+1,  2+6*j+1+1) = M_p(3*i+1,  3*j+1+1); % wx
                        M_out(2+6*i+1,  2+6*j+2+1) = M_p(3*i+1,  3*j+2+1); % wy
                        M_out(2+6*i+1+1,2+6*j+1)   = M_p(3*i+1+1,3*j+1);   % xw
                        M_out(2+6*i+1+1,2+6*j+1+1) = M_p(3*i+1+1,3*j+1+1); % xx
                        M_out(2+6*i+1+1,2+6*j+2+1) = M_p(3*i+1+1,3*j+2+1); % xy
                        M_out(2+6*i+2+1,2+6*j+1)   = M_p(3*i+2+1,3*j+1);   % yw
                        M_out(2+6*i+2+1,2+6*j+1+1) = M_p(3*i+2+1,3*j+1+1); % yx
                        M_out(2+6*i+2+1,2+6*j+2+1) = M_p(3*i+2+1,3*j+2+1); % yy
                    end
                end
                TSub = [trafo,zeros(3,3);zeros(3,3),trafo];
                KeSub = zeros(6,6);
                KeNew = zeros(18,18);
                MSub = zeros(6,6);
                MNew = zeros(18,18);
                %Kl = K_out;
                Kg = zeros(18,18);
                Mg = zeros(18,18);
                for i=0:2
                    for j = 0:2
                        % copy values into temporary sub-matrix for correct format to transformation
                        for k=0:5
                            for l=0:5
                                KeSub(k+1,l+1) = K_out(i*6+k+1,j*6+l+1);
                                MSub(k+1,l+1) = M_out(i*6+k+1,j*6+l+1);
                            end
                        end
                        % the actual transformation step
                        KeSub=TSub'*KeSub*TSub;
                        MSub=TSub'*MSub*TSub;
                        % copy transformed values into new global stiffness matrix
                        for k=0:5
                            for l=0:5
                                KeNew(i*6+k+1,j*6+l+1) = KeSub(k+1,l+1);
                                MNew(i*6+k+1,j*6+l+1) = MSub(k+1,l+1);
                            end
                        end
                    end
                end
        
                for alpha = 0:5
                    for beta=0:5
                        for i=0:2
                            for j=0:2
                                Kg(3*alpha+i+1,3*beta+j+1) = KeNew(6*i+alpha+1,6*j+beta+1);
                                Mg(3*alpha+i+1,3*beta+j+1) = MNew(6*i+alpha+1,6*j+beta+1);
                            end  
                        end
                    end
                end
                linearidx = zeros(1,size(obj.femutils.IndexRow{iter}(:),1));
                for i = 1:size(obj.femutils.IndexRow{iter}(:),1)
                    linearidx(i) = sub2ind(size(obj.femLHS),obj.femutils.IndexRow{iter}(i),obj.femutils.IndexCol{iter}(i));
                end
                obj.femLHS(linearidx) = obj.femLHS(linearidx)+Kg(:)';
                obj.femMass(linearidx) = obj.femMass(linearidx)+Mg(:)';
                obj.femDamp(linearidx) = obj.femDamp(linearidx)+Mg(:)'./obj.femRho(iter,1).*obj.femVisc(iter,1);
                
            end
            obj.femLHS(obj.femutils.MatIndex==0,:)=[];
            obj.femLHS(:,obj.femutils.MatIndex==0)=[];
            obj.femMass(obj.femutils.MatIndex==0,:)=[];
            obj.femMass(:,obj.femutils.MatIndex==0)=[];
            obj.femDamp(obj.femutils.MatIndex==0,:)=[];
            obj.femDamp(:,obj.femutils.MatIndex==0)=[];

            if skipFlag == 0
                %すべてが0の項の処理
                delIndex = or(all(obj.femMass==0), all(obj.femLHS==0));
                obj.femLHS(delIndex,:) = [];
                obj.femLHS(:,delIndex) = [];
                obj.femMass(delIndex,:) = [];
                obj.femMass(:,delIndex) = [];
                obj.femDamp(delIndex,:) = [];
                obj.femDamp(:,delIndex) = [];
                obj.femutils.InvMatIndex(delIndex,:) = [];
                delIndexNo = find(delIndex);
                si = 1;
                delIndex2 = [];
                for i = 1:size(obj.femutils.MatIndex,1)
                    if obj.femutils.MatIndex(i) == 1
                        if any(si==delIndexNo)
                            delIndex2 = [delIndex2,i];
                        end
                         si = si+1;
                    end
                end
                obj.femutils.MatIndex(delIndex2) = 0;
                [U, S, V] = svds(obj.femLHS,size(obj.femLHS,1));
                obj.femLHS = U*S*V';
                [U, S, V] = svds(obj.femMass,size(obj.femMass,1));
                obj.femMass = U*S*V';
                [U, S, V] = svds(obj.femDamp,size(obj.femDamp,1));
                obj.femDamp = U*S*V';
            end

        end

        function obj = modifyMassDampingMatrix(obj)
            nu = 0.34;
            qps = [1/6,1/6;2/3,1/6;1/6,2/3];
            obj.femMass = sparse(6*numel(obj.femutils.usedVerts),6*numel(obj.femutils.usedVerts));
            obj.femDamp = sparse(6*numel(obj.femutils.usedVerts),6*numel(obj.femutils.usedVerts));
            for iter = 1:size(obj.femtri.ConnectivityList,1)
                Dp(1,1) = 1.0; Dp(1,2) = nu;
                Dp(2,1) = nu;  Dp(2,2) = 1.0;
                Dp(3,3) = (1.0-nu)/2.0;
                Dp = Dp.*(obj.femE(iter,1)*obj.femThn(iter,1)^3/(12.0*(1.0-nu*nu)));
                U = obj.femtri.Points(obj.femtri.ConnectivityList(iter,2),:)-obj.femtri.Points(obj.femtri.ConnectivityList(iter,1),:);
                V = obj.femtri.Points(obj.femtri.ConnectivityList(iter,3),:)-obj.femtri.Points(obj.femtri.ConnectivityList(iter,1),:);
                transUV(:,1) = U;
                transUV(:,2) = V;
                W = cross(U,V);
                trafo(1,:) = U./norm(U);
                trafo(3,:) = W./norm(W);
                trafo(2,:) = cross(trafo(3,:),trafo(1,:));
                transUV = trafo*transUV;
                dphi(1,1) = -transUV(1,1); % x12 = x1-x2 = 0-x2 = -x2
                dphi(2,1) =  transUV(1,2); % x31 = x3-x1 = x3-0 = x3
                dphi(3,1) =  transUV(1,1)-transUV(1,2); % x23 = x2-x3
                dphi(1,2) = -transUV(2,1); % y12 = y1-y2 = -y2 = 0 (stays zero, as node B and A lie on local x-axis and therefore)
                dphi(2,2) =  transUV(2,2); % y31 = y3-y1 = y3-0 = y3
                dphi(3,2) =  transUV(2,1)-transUV(2,2); % y23 = y2-y3 = 0-y3 = -y3

                M_m= diag(ones(1,6).*0.5);
                M_m(1,3) = 0.25;
                M_m(1,5) = 0.25;
                M_m(2,4) = 0.25;
                M_m(2,6) = 0.25;
                M_m(3,1) = 0.25;
                M_m(3,5) = 0.25;
                M_m(4,2) = 0.25;
                M_m(4,6) = 0.25;
                M_m(5,1) = 0.25;
                M_m(5,3) = 0.25;
                M_m(6,2) = 0.25;
                M_m(6,4) = 0.25;
                M_m = M_m .* obj.femarea(iter,1).*obj.femRho(iter,1).*obj.femThn(iter,1)./3;

                
                sidelen(1) = dphi(1,1)^2+dphi(1,2)^2;
                sidelen(2) = dphi(2,1)^2+dphi(2,2)^2;
                sidelen(3) = dphi(3,1)^2+dphi(3,2)^2;
                Y(1,1) = dphi(3,2)^2.0;
                Y(1,2) = dphi(2,2)^2.0;
                Y(1,3) = dphi(3,2)*dphi(2,2);
                Y(2,1) = dphi(3,1)^2.0;
                Y(2,2) = dphi(2,1)^2.0;
                Y(2,3) = dphi(2,1)*dphi(3,1);
                Y(3,1) = -2.0*dphi(3,1)*dphi(3,2);
                Y(3,2) = -2.0*dphi(2,1)*dphi(2,1);
                Y(3,3) = -dphi(3,1)*dphi(2,2)-dphi(2,1)*dphi(3,2);
                Y = Y.* 1.0/(4.0*obj.femarea(iter,1)^2.0);
                Ainv = zeros(9,9);
                x12 = -transUV(1,1);
                y12 = -transUV(2,1);
                x31 = transUV(1,2);
                y31 = transUV(2,2);
                x23 = transUV(1,1)-transUV(1,2); 
                y23 =  transUV(2,1)-transUV(2,2);
                Ainv(1,1) = 1;
                Ainv(2,4) = 1;
                Ainv(3,7) = 1;
                Ainv(4,1) = -1;
                Ainv(4,4) = 1;
                Ainv(4,5) = y12;
                Ainv(4,6) = -x12;
                Ainv(5,4) = -1;
                Ainv(5,7) = 1;
                Ainv(5,8) = y23;
                Ainv(5,9) = -x23;
                Ainv(6,1) = 1;
                Ainv(6,2) = y31;
                Ainv(6,3) = -x31;
                Ainv(6,7) = -1;
                Ainv(7,1) = 2;
                Ainv(7,2) = -y12;
                Ainv(7,3) = x12;
                Ainv(7,4) = -2;
                Ainv(7,5) = -y12;
                Ainv(7,6) = x12;
                Ainv(8,4) = 2;
                Ainv(8,5) = -y23;
                Ainv(8,6) = x23;
                Ainv(8,7) = -2;
                Ainv(8,8) = -y23;
                Ainv(8,9) = x23;
                Ainv(9,1) = -2;
                Ainv(9,2) = -y31;
                Ainv(9,3) =  x31;
                Ainv(9,7) = 2;
                Ainv(9,8) = -y31;
                Ainv(9,9) = x31;
                M_p = zeros(9,9);
                M_m2 = zeros(6,6);
                for j = 1:3
                    N = obj.evalNTri(sidelen, qps(j,1), qps(j,2), Ainv);
                    Nm = obj.evalNmTri(sidelen, qps(j,1), qps(j,2), Ainv);
                    temp2 = N*N'./6;
                    temp3 = Nm'*Nm./6;
                    M_p = M_p+temp2;
                    M_m2 = M_m2+temp3;
                end
                M_p = M_p*obj.femRho(iter,1)*obj.femThn(iter,1)*2.*obj.femarea(iter,1);
                M_out = zeros(18,18);
                for i=0:2
                    for j = 0:2 
                        M_out(  6*i+1,    6*j+1)   = M_m(2*i+1,  2*j+1);   % uu
                        M_out(  6*i+1,    6*j+1+1) = M_m(2*i+1,  2*j+1+1); % uv
                        M_out(  6*i+1+1,  6*j+1)   = M_m(2*i+1+1,2*j+1);   % vu
                        M_out(  6*i+1+1,  6*j+1+1) = M_m(2*i+1+1,2*j+1+1); % vv
                        M_out(2+6*i+1,  2+6*j+1)   = M_p(3*i+1,  3*j+1);   % ww
                        M_out(2+6*i+1,  2+6*j+1+1) = M_p(3*i+1,  3*j+1+1); % wx
                        M_out(2+6*i+1,  2+6*j+2+1) = M_p(3*i+1,  3*j+2+1); % wy
                        M_out(2+6*i+1+1,2+6*j+1)   = M_p(3*i+1+1,3*j+1);   % xw
                        M_out(2+6*i+1+1,2+6*j+1+1) = M_p(3*i+1+1,3*j+1+1); % xx
                        M_out(2+6*i+1+1,2+6*j+2+1) = M_p(3*i+1+1,3*j+2+1); % xy
                        M_out(2+6*i+2+1,2+6*j+1)   = M_p(3*i+2+1,3*j+1);   % yw
                        M_out(2+6*i+2+1,2+6*j+1+1) = M_p(3*i+2+1,3*j+1+1); % yx
                        M_out(2+6*i+2+1,2+6*j+2+1) = M_p(3*i+2+1,3*j+2+1); % yy
                    end
                end
                TSub = [trafo,zeros(3,3);zeros(3,3),trafo];
                MSub = zeros(6,6);
                MNew = zeros(18,18);
                %Kl = K_out;
                Mg = zeros(18,18);
                for i=0:2
                    for j = 0:2
                        % copy values into temporary sub-matrix for correct format to transformation
                        for k=0:5
                            for l=0:5
                                MSub(k+1,l+1) = M_out(i*6+k+1,j*6+l+1);
                            end
                        end
                        % the actual transformation step
                        MSub=TSub'*MSub*TSub;
                        % copy transformed values into new global stiffness matrix
                        for k=0:5
                            for l=0:5
                                MNew(i*6+k+1,j*6+l+1) = MSub(k+1,l+1);
                            end
                        end
                    end
                end
        
                for alpha = 0:5
                    for beta=0:5
                        for i=0:2
                            for j=0:2
                                Mg(3*alpha+i+1,3*beta+j+1) = MNew(6*i+alpha+1,6*j+beta+1);
                            end  
                        end
                    end
                end
                linearidx = zeros(1,size(obj.femutils.IndexRow{iter}(:),1));
                for i = 1:size(obj.femutils.IndexRow{iter}(:),1)
                    linearidx(i) = sub2ind(size(obj.femLHS),obj.femutils.IndexRow{iter}(i),obj.femutils.IndexCol{iter}(i));
                    %obj.femLHS(obj.femutils.IndexRow{iter}(i),obj.femutils.IndexCol{iter}(i)) = obj.femLHS(obj.femutils.IndexRow{iter}(i),obj.femutils.IndexCol{iter}(i))+Kg(i);
                end
                obj.femMass(linearidx) = obj.femMass(linearidx)+Mg(:)';
                obj.femDamp(linearidx) = obj.femDamp(linearidx)+Mg(:)'./obj.femRho(iter,1).*obj.femVisc(iter,1);
                if mod(iter,floor(size(obj.femtri.ConnectivityList,1)/10))==1
                    fprintf("%d / %d \n",iter,size(obj.femtri.ConnectivityList,1));
                end
                
            end
            obj.femMass(obj.femutils.MatIndex==0,:)=[];
            obj.femMass(:,obj.femutils.MatIndex==0)=[];
            obj.femDamp(obj.femutils.MatIndex==0,:)=[];
            obj.femDamp(:,obj.femutils.MatIndex==0)=[];

            %すべてが0の項の処理
            delIndex = or(all(obj.femMass==0), all(obj.femLHS==0));
            obj.femLHS(delIndex,:) = [];
            obj.femLHS(:,delIndex) = [];
            obj.femMass(delIndex,:) = [];
            obj.femMass(:,delIndex) = [];
            obj.femDamp(delIndex,:) = [];
            obj.femDamp(:,delIndex) = [];
            obj.femutils.InvMatIndex(delIndex,:) = [];
            delIndexNo = find(delIndex);
            si = 1;
            delIndex2 = [];
            for i = 1:size(obj.femutils.MatIndex,1)
                if obj.femutils.MatIndex(i) == 1
                    if any(si==delIndexNo)
                        delIndex2 = [delIndex2,i];
                    end
                     si = si+1;
                end
            end
            obj.femutils.MatIndex(delIndex2) = 0;

            [U, S, V] = svds(obj.femLHS,size(obj.femLHS,1));
            obj.femLHS = U*S*V';
            [U, S, V] = svds(obj.femMass,size(obj.femMass,1));
            obj.femMass = U*S*V';
            [U, S, V] = svds(obj.femDamp,size(obj.femDamp,1));
            obj.femDamp = U*S*V';            

        end

        function [delta,deltadot,S] = solveFem(obj,distLoad,selfLoadFlag)
            %distLoad:空力解析メッシュにかかる分布加重。圧力など
            if nargin < 3
                selfLoadFlag = 1; %defaultは自重による荷重を含める。
            end
            %空力解析メッシュから構造解析メッシュへ圧力分布を投影
            S = [];
            for iter = 1:size(distLoad,2)
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad(:,iter).*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
    
                femRHSp = Fp(obj.femutils.MatIndex==1,1);
                delta_p = lsqminnorm(obj.femLHS,femRHSp,obj.settingUNLSI.lsqminnormTol);
                %delta_p = obj.femMass\femRHSp;
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=delta_p;
                delta{iter} = zeros(size(obj.femtri.Points,1),6);
                delta{iter}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                deltadot{iter} = zeros(size(obj.femtri.Points,1),6);
                Ssolve = obj.femLHS*delta_p - femRHSp;
                S = [S;Ssolve];
            end
        end

        function [delta,deltadot,z,zdot,S] = solveModalFem(obj,distLoad,selfLoadFlag)
            modalNo = size(obj.femEigenVec,2);
            Fmodal = obj.femEigenVec(:,1:modalNo);
            %distLoad:空力解析メッシュにかかる分布加重。圧力など
            if nargin < 3
                selfLoadFlag = 1; %defaultは自重による荷重を含める。
            end
            %空力解析メッシュから構造解析メッシュへ圧力分布を投影
            S = [];
            for iter = 1:size(distLoad,2)
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad(:,iter).*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
    
                femRHSp = Fp(obj.femutils.MatIndex==1,1);
                fmode = Fmodal'*femRHSp;
                
                z{iter} = lsqminnorm(Fmodal'*obj.femLHS*Fmodal,fmode,obj.settingUNLSI.lsqminnormTol);
                zdot{iter} = zeros(size(z{iter}));
                delta_p = Fmodal*z{iter};
                %delta_p = obj.femMass\femRHSp;
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=delta_p;
                delta{iter} = zeros(size(obj.femtri.Points,1),6);
                delta{iter}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                deltadot{iter} = zeros(size(obj.femtri.Points,1),6);
                Ssolve = obj.femLHS*delta_p - femRHSp;
                S = [S;Ssolve];
            end
        end

        function [delta,deltadot,S,obj] = solveFemForAdjoint(obj,delta,distLoad,selfLoadFlag)
            %distLoad:空力解析メッシュにかかる分布加重。圧力など
            if nargin < 4
                selfLoadFlag = 1; %defaultは自重による荷重を含める。
            end
            %空力解析メッシュから構造解析メッシュへ圧力分布を投影
            S = [];
            for iter = 1:size(distLoad,2)
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad(:,iter).*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
    
                femRHSp = Fp(obj.femutils.MatIndex==1,1);
                yb = sparse(vertcat(delta{iter}(:,1),delta{iter}(:,2),delta{iter}(:,3),delta{iter}(:,4),delta{iter}(:,5),delta{iter}(:,6)));
                delta_p = yb(obj.femutils.MatIndex==1,1);
                Ssolve = obj.femLHS*delta_p - femRHSp;
                %delta_p = lsqminnorm(obj.femLHS,femRHSp,obj.settingUNLSI.lsqminnormTol);
                %delta_p = obj.femMass\femRHSp;
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=delta_p;
                delta{iter} = zeros(size(obj.femtri.Points,1),6);
                delta{iter}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                deltadot{iter} = zeros(size(obj.femtri.Points,1),6);
                S = [S;Ssolve];
            end
        end

        function obj = femModalAnalysis(obj,modalNo)
            [eigenVecBuff,eigenValBuff] = eigs(obj.femLHS,obj.femMass,modalNo,"smallestabs");

            [obj.femEigenVal,I] = sort(full(diag(eigenValBuff)));
            obj.femEigenVec = eigenVecBuff(:,I);
            nanIndex = isnan(obj.femEigenVal);
            obj.femEigenVal(nanIndex) = [];
            obj.femEigenVec(:,nanIndex) = [];
        end

        function [delta,deltadot]  = solveAeroelastic(obj,tspan,delta,deltadot,distLoad,selfLoadFlag)
            %distLoad:空力解析メッシュにかかる分布加重。圧力など
            if nargin < 6
                selfLoadFlag = 1; %defaultは自重による荷重を含める。
            end
            %空力解析メッシュから構造解析メッシュへ圧力分布を投影
            if isempty(delta)
                for iter = 1:size(distLoad,2)
                    delta{iter} = zeros(obj.femutils.nbVerts,6);
                end
            end
            if isempty(deltadot)
                for iter = 1:size(distLoad,2)
                    deltadot{iter} = zeros(obj.femutils.nbVerts,6);
                end
            end

            for iter = 1:size(distLoad,2)
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad(:,iter).*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
                femRHSp = Fp(obj.femutils.MatIndex==1,1);
                MassEx = blkdiag(eye(size(femRHSp,1)),obj.femMass);
                AEx = [zeros(size(femRHSp,1)),eye(size(femRHSp,1));-obj.femLHS,-obj.femDamp];
                fEx = [zeros(size(femRHSp,1),1);femRHSp];
                options = odeset('Mass',MassEx,"Jacobian",AEx,"Vectorized","on","RelTol",obj.settingUNLSI.ode23RelTol,"AbsTol",obj.settingUNLSI.ode23AbsTol,"NormControl","off");
                %deltaから初期値の復元
                yb = vertcat(delta{iter}(:,1),delta{iter}(:,2),delta{iter}(:,3),delta{iter}(:,4),delta{iter}(:,5),delta{iter}(:,6));
                y = yb(obj.femutils.MatIndex==1,1);
                ydb = vertcat(deltadot{iter}(:,1),deltadot{iter}(:,2),deltadot{iter}(:,3),deltadot{iter}(:,4),deltadot{iter}(:,5),deltadot{iter}(:,6));
                yd = ydb(obj.femutils.MatIndex==1,1);
                warning('off','MATLAB:singularMatrix')
                %try
                %    [~,ybuff] = ode15s(@(t,d)(AEx*d+fEx),tspan,[y;yd],options);
                %catch
                    [~,ybuff] = ode23(@(t,d)(AEx*d+fEx),tspan,[y;yd],options);
                %end
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=ybuff(end,1:numel(y))';
                disp_buff(isnan(disp_buff)) = 0;
                delta{iter} = zeros(size(obj.femtri.Points,1),6);    
                delta{iter}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                delta{iter}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=ybuff(end,numel(y)+1:end)';
                deltadot{iter} = zeros(size(obj.femtri.Points,1),6);
                deltadot{iter}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                deltadot{iter}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                deltadot{iter}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                deltadot{iter}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                deltadot{iter}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                deltadot{iter}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                warning('on')
            end
        end

        function [deltaLoad2] = removeFunc(obj,deltaLoad)
            if size(deltaLoad,2) == 1
                deltaLoad(obj.femutils.MatIndex==0,:)=[];
            else
                deltaLoad(obj.femutils.MatIndex==0,:)=[];
                deltaLoad(:,obj.femutils.MatIndex==0)=[];
            end
            % deltaLoad(obj.femutils.MatIndex==0,:)=[];
            % deltaLoad(:,obj.femutils.MatIndex==0)=[];
            deltaLoad2 = deltaLoad;
        end

        function [z,zdot,delta,deltadot]  = solveModalAeroelastic(obj,tspan,z,zdot,distLoad,selfLoadFlag)

            modalNo = size(obj.femEigenVec,2);
            %distLoad:空力解析メッシュにかかる分布加重。圧力など
            %空力解析メッシュから構造解析メッシュへ圧力分布を投影
            Fmodal = obj.femEigenVec(:,1:modalNo);
            if isempty(z)
                for iter = 1:size(distLoad,2)
                    z{iter} = zeros(modalNo,1);
                end
            end
            if isempty(zdot)
                for iter = 1:size(distLoad,2)
                    zdot{iter} = zeros(modalNo,1);
                end
            end

            for iter = 1:size(distLoad,2)
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad(:,iter).*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
                femRHSp = Fp(obj.femutils.MatIndex==1,1);
                fmode = Fmodal'*femRHSp;
                MassEx = blkdiag(eye(size(fmode,1)),Fmodal'*obj.femMass*Fmodal);
                AEx = [zeros(size(fmode,1)),eye(size(fmode,1));-Fmodal'*obj.femLHS*Fmodal,-Fmodal'*obj.femDamp*Fmodal];
                fEx = [zeros(size(fmode,1),1);fmode];
                options = odeset('Mass',MassEx,"Jacobian",AEx,"Vectorized","on","RelTol",obj.settingUNLSI.ode23RelTol,"AbsTol",obj.settingUNLSI.ode23AbsTol,"NormControl","off");
                %deltaから初期値の復元
                warning('off','MATLAB:singularMatrix')
                %try
                %    [~,ybuff] = ode15s(@(t,d)(AEx*d+fEx),tspan,[y;yd],options);
                %catch
                    [~,ybuff] = ode45(@(t,d)(AEx*d+fEx),tspan,[z{iter}(:,1);zdot{iter}(:,1)],options);
                %end
                z{iter}(:,1) = ybuff(end,1:modalNo)';
                zdot{iter}(:,1) = ybuff(end,modalNo+1:end);
                warning('on')
                if nargout>2
                    delta{iter} = obj.femSol2Delta(Fmodal*z{iter});
                    deltadot{iter} = obj.femSol2Delta(Fmodal*zdot{iter});
                end
            end
        end

        function [f,jac] = calcStrongCouplingJaccobian(obj,u,delta,alpha,beta,dynPress)
            %deltaから機体メッシュを再生成
            nflowVar = size(obj.LHS,1);
            dS_du = [];
            dR_ddelta = [];
            calcCount = 1;
            R0all = [];
            S0all = [];
            for iter = 1:numel(obj.flow)
                alphabuff = alpha;
                betabuff = beta;
                dynbuff = dynPress;
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
                    [~,~,S0] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNLSI.selfLoadFlag);
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
                            [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNLSI.selfLoadFlag);
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
                                [~,~,Sf] = obj.solveFemForAdjoint(deltaIn,Cp{iter}.*dynbuff(jter),obj.settingUNLSI.selfLoadFlag);
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

            for i = 1:1
                % if i == 1
                    % if obj.flowNoList(i,3) == 1
                        dR_du = (obj2.LHS+obj2.wakeLHS); %亜音速
                    % else
                        % dR_du = eye(nbPanel);
                    % end
                % else
                %     if obj.flowNoList(i,3) == 1
                %         dR_du = blkdiag(dR_du,(obj2.LHS+obj2.wakeLHS));
                %     else
                %         nbPanel = sum(obj.paneltype == 1);
                %         dR_du = blkdiag(dR_du,eye(nbPanel));
                %     end
                % end
            end
            if nargout>1
                jac= [dR_du,dR_ddelta;dS_du./10000,dS_ddelta./10000];
            end
            f = [R0all;S0all./10000];

        end


        %空力弾性のフラッター速度解析用のdelta偏微分を取得する
        function obj = calcAeroelasticDerivative(obj,obj2,delta,alpha,beta,Mach,Re)
            %deltaから機体メッシュを再生成
            calcCount = 1;
            for iter = 1:numel(alpha)
                %delta微分の計算
                %dR_ddelta
                deltap0_buff = delta{calcCount}(:);
                deltap0 = deltap0_buff(obj.femutils.MatIndex==1,1);
                %基準の解析
                deltapf = deltap0;
                disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
                disp_buff(obj.femutils.InvMatIndex,1)=deltapf;
                delta0 = delta;
                delta0{calcCount} = zeros(size(obj.femtri.Points,1),6);
                delta0{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
                delta0{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
                delta0{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
                delta0{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
                delta0{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
                delta0{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
                modVerts = obj.calcModifiedVerts(delta0{calcCount});
                obj3 = obj2.makeApproximatedInstance(modVerts,1);
                obj3 = obj3.solveFlow(alpha(iter),beta(iter),Mach(iter),Re(iter));
                distLoad = obj3.getCp(alpha(iter),beta(iter),Mach(iter),Re(iter));
                Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                interpLoad = obj.verts2centerMat'*(distLoad.*obj.area);
                vNormal = vertexNormal(obj.tri);
                %femPointLoad = obj.fem2aeroMat'*interpLoad;
                vertsLoad = vNormal.*interpLoad;
                selfLoadFlag = 1;
                if selfLoadFlag == 1
                    selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                end
                femPointLoad = zeros(size(obj.femtri.Points,1),3);
                femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                if selfLoadFlag == 1
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                else
                    femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                end
                for i = 1:size(obj.femutils.usedVerts,2) 
                    Fp(i,1) = -femPointLoad(i,1);
                    Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                    Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                end
                femRHSp0 = Fp(obj.femutils.MatIndex==1,1);
                pert = sqrt(eps);
                for k = 1:size(deltap0,1)
                    if k <= sum(obj.femutils.MatIndex(1:sub2ind(size(delta{calcCount}),1,4)))
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
                        obj3 = obj3.solveFlow(alpha(iter),beta(iter),Mach(iter),Re(iter));
                        distLoad = obj3.getCp(alpha(iter),beta(iter),Mach(iter),Re(iter));
                        Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
                        interpLoad = obj.verts2centerMat'*(distLoad.*obj.area);
                        vNormal = vertexNormal(obj.tri);
                        %femPointLoad = obj.fem2aeroMat'*interpLoad;
                        vertsLoad = vNormal.*interpLoad;
                        if selfLoadFlag == 1
                            selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
                        end
                        femPointLoad = zeros(size(obj.femtri.Points,1),3);
                        femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
                        femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
                        if selfLoadFlag == 1
                            femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
                        else
                            femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
                        end
                        for i = 1:size(obj.femutils.usedVerts,2) 
                            Fp(i,1) = -femPointLoad(i,1);
                            Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
                            Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
                        end
                        femRHSp = Fp(obj.femutils.MatIndex==1,1);
                        obj.df_ddelta(:,k) = (femRHSp-femRHSp0)./pert;
                        disp(k);
                    else
                        obj.df_ddelta(:,k) = 0;
                    end
                end
                calcCount = calcCount + 1;
            end

        end

        % function obj = calcfemModalDerivative(obj,obj2,z0,alpha,beta,Mach,Re)
        %     %deltaから機体メッシュを再生成
        %     calcCount = 1;
        %     for iter = 1:numel(alpha)
        %         %delta微分の計算
        %         %dR_ddelta
        %         deltap0_buff = z0;%基準z0
        %         deltap0 = obj.femEigenVec*z0;%V*z0
        %         %基準の解析
        %         deltapf = deltap0;
        %         disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
        %         disp_buff(obj.femutils.InvMatIndex,1)=deltapf;
        %         delta0 = delta;
        %         delta0{calcCount} = zeros(size(obj.femtri.Points,1),6);
        %         delta0{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
        %         delta0{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
        %         delta0{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
        %         delta0{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
        %         delta0{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
        %         delta0{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
        %         modVerts = obj.calcModifiedVerts(delta0{calcCount});
        %         obj3 = obj2.makeApproximatedInstance(modVerts,1);
        %         obj3 = obj3.solveFlow(alpha(iter),beta(iter),Mach(iter),Re(iter));
        %         distLoad = obj3.getCp(alpha(iter),beta(iter),Mach(iter),Re(iter));
        %         Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
        %         interpLoad = obj.verts2centerMat'*(distLoad.*obj.area);
        %         vNormal = vertexNormal(obj.tri);
        %         %femPointLoad = obj.fem2aeroMat'*interpLoad;
        %         vertsLoad = vNormal.*interpLoad;
        %         selfLoadFlag = 1;
        %         if selfLoadFlag == 1
        %             selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
        %         end
        %         femPointLoad = zeros(size(obj.femtri.Points,1),3);
        %         femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
        %         femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
        %         if selfLoadFlag == 1
        %             femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
        %         else
        %             femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
        %         end
        %         for i = 1:size(obj.femutils.usedVerts,2) 
        %             Fp(i,1) = -femPointLoad(i,1);
        %             Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
        %             Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
        %         end
        %         femRHSp0 = Fp(obj.femutils.MatIndex==1,1);
        %         pert = sqrt(eps);
        %         for k = 1:size(deltap0,1)
        %             if k <= sum(obj.femutils.MatIndex(1:sub2ind(size(delta{calcCount}),1,4)))
        %                 deltapf = deltap0;
        %                 deltapf(k,1) = deltap0(k,1)+pert;%pert→(z0+pert)*V
        %                 disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
        %                 disp_buff(obj.femutils.InvMatIndex,1)=deltapf;
        %                 deltaf = delta;
        %                 deltaf{calcCount} = zeros(size(obj.femtri.Points,1),6);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
        %                 deltaf{calcCount}(obj.femutils.usedVerts,6)=disp_buff(5*obj.femutils.nbVerts+1:6*obj.femutils.nbVerts,1);
        %                 modVerts = obj.calcModifiedVerts(deltaf{calcCount});
        %                 obj3 = obj2.makeApproximatedInstance(modVerts,1);
        %                 obj3 = obj3.solveFlow(alpha(iter),beta(iter),Mach(iter),Re(iter));
        %                 distLoad = obj3.getCp(alpha(iter),beta(iter),Mach(iter),Re(iter));
        %                 Fp = sparse(6*size(obj.femutils.usedVerts,2),1);
        %                 interpLoad = obj.verts2centerMat'*(distLoad.*obj.area);
        %                 vNormal = vertexNormal(obj.tri);
        %                 %femPointLoad = obj.fem2aeroMat'*interpLoad;
        %                 vertsLoad = vNormal.*interpLoad;
        %                 if selfLoadFlag == 1
        %                     selfLoad = obj.femverts2centerMat'*(9.8.*obj.femarea.*obj.femThn.*obj.femRho);
        %                 end
        %                 femPointLoad = zeros(size(obj.femtri.Points,1),3);
        %                 femPointLoad(obj.femutils.usedVerts,1) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,1);
        %                 femPointLoad(obj.femutils.usedVerts,2) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,2);
        %                 if selfLoadFlag == 1
        %                     femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3) + selfLoad(obj.femutils.usedVerts,1);
        %                 else
        %                     femPointLoad(obj.femutils.usedVerts,3) = obj.fem2aeroMat'*vertsLoad(obj.femutils.usedAeroVerts,3);
        %                 end
        %                 for i = 1:size(obj.femutils.usedVerts,2) 
        %                     Fp(i,1) = -femPointLoad(i,1);
        %                     Fp(i+size(obj.femtri.Points,1),1) = -femPointLoad(i,2);
        %                     Fp(i+2*size(obj.femtri.Points,1),1) = -femPointLoad(i,3);
        %                 end
        %                 femRHSp = Fp(obj.femutils.MatIndex==1,1);
        %                 obj.df_ddelta(:,k) = (femRHSp-femRHSp0)./pert;
        %                 disp(k);
        %             else
        %                 obj.df_ddelta(:,k) = 0;
        %             end
        %         end
        %         calcCount = calcCount + 1;
        %     end


        function delta = femSol2Delta(obj,femSol)
            disp_buff = zeros(size(obj.femtri.Points,1)*6,1);
            disp_buff(obj.femutils.InvMatIndex,1)=femSol;
            delta = zeros(size(obj.femtri.Points,1),6);
            delta(obj.femutils.usedVerts,1)=disp_buff(1:obj.femutils.nbVerts,1);
            delta(obj.femutils.usedVerts,2)=disp_buff(obj.femutils.nbVerts+1:2*obj.femutils.nbVerts,1);
            delta(obj.femutils.usedVerts,3)=disp_buff(2*obj.femutils.nbVerts+1:3*obj.femutils.nbVerts,1);
            delta(obj.femutils.usedVerts,4)=disp_buff(3*obj.femutils.nbVerts+1:4*obj.femutils.nbVerts,1);
            delta(obj.femutils.usedVerts,5)=disp_buff(4*obj.femutils.nbVerts+1:5*obj.femutils.nbVerts,1);
        end

        function deltaInterp = interpFemDelta(obj,delta,xyzinAeroMesh)
            deltaAero = zeros(size(obj.tri.Points));
            deltaAero(obj.femutils.usedAeroVerts,1) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,1);
            deltaAero(obj.femutils.usedAeroVerts,2) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,2);
            deltaAero(obj.femutils.usedAeroVerts,3) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,3);
            deltaAero(obj.femutils.usedAeroVerts,4) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,4);
            deltaAero(obj.femutils.usedAeroVerts,5) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,5);
            deltaAero(obj.femutils.usedAeroVerts,6) = obj.fem2aeroMat*delta(obj.femutils.usedVerts,6);

            %空力メッシュから補間点への補間
            phi = @(r)obj.phi3(r,1);
            Avv = phi(obj.calcRMat(obj.tri.Points,obj.tri.Points));
            Acv = phi(obj.calcRMat(obj.tri.Points,xyzinAeroMesh));
            
            interpMat = Acv'*pinv(Avv);
            deltaInterp(:,1)  = interpMat*deltaAero(:,1);
            deltaInterp(:,2)  = interpMat*deltaAero(:,2);
            deltaInterp(:,3)  = interpMat*deltaAero(:,3);
            deltaInterp(:,4)  = interpMat*deltaAero(:,4);
            deltaInterp(:,5)  = interpMat*deltaAero(:,5);
            deltaInterp(:,6)  = interpMat*deltaAero(:,6);

        end

        function modVerts = calcModifiedVerts(obj,delta)
            %メッシュ変形を使って変位後の解析メッシュを計算する
            modVerts = obj.tri.Points;
            modVerts(obj.femutils.usedAeroVerts,1) = obj.tri.Points(obj.femutils.usedAeroVerts,1) + obj.fem2aeroMat*delta(obj.femutils.usedVerts,1);
            modVerts(obj.femutils.usedAeroVerts,2) = obj.tri.Points(obj.femutils.usedAeroVerts,2) + obj.fem2aeroMat*delta(obj.femutils.usedVerts,2);
            modVerts(obj.femutils.usedAeroVerts,3) = obj.tri.Points(obj.femutils.usedAeroVerts,3) + obj.fem2aeroMat*delta(obj.femutils.usedVerts,3);
        end
        
        function  exportDrawing(obj,filename,exportID,intersectYpos,outputYpos)
           exportTri.faces = obj.tri.ConnectivityList(any(obj.surfID==exportID(:)',2),:);
           exportTri.vertices = obj.tri.Points;
           cutBox = [min(exportTri.vertices(unique(exportTri.faces(:)),[1,3]));max(exportTri.vertices(unique(exportTri.faces(:)),[1,3]))];
           cutSurf.vertices = [cutBox(1,1),intersectYpos,cutBox(1,2);cutBox(2,1),intersectYpos,cutBox(1,2);cutBox(1,1),intersectYpos,cutBox(2,2);cutBox(2,1),intersectYpos,cutBox(2,2)];
           cutSurf.faces = [1,2,4;1,4,3];
           [intMatrix, intSurface] = obj.SurfaceIntersection(obj,exportTri, cutSurf);
           orgChord = [intSurface.vertices(:,1),intSurface.vertices(:,3)];
           [~,LEindex] = min(orgChord(:,1));
           xyorg = orgChord(LEindex,:);
           [~,TEindex] = max(orgChord(:,1));
           xymax = orgChord(TEindex,:);
           alfa = -atan2(xymax(2)-xyorg(2),xymax(1)-xyorg(1));
           chord = norm(xymax-xyorg);
           modChord = (orgChord-xyorg)*[cos(alfa),sin(alfa);-sin(alfa),cos(alfa)]./chord;
           %
           A = [modChord(:,1).^4,modChord(:,1).^3,modChord(:,1).^2,modChord(:,1)];
           c = [1,1,1,1];
           g = [2.*A'*A,c';c,0];
           b = [2.*A'*modChord(:,2);0];
           w = g\b;
           y = A*w(1:4);
           isUpper = (modChord(:,2)-y)>0;
           isLower = (modChord(:,2)-y)<0;

           figure(1);clf;hold on;
           plot(modChord(:,1),modChord(:,2));
           plot(modChord(:,1),y);
           %前縁後縁の点を配して
           x0 = [-0.33,0.2,alfa,chord,0.126,0.3669,-0.13634,0.78991,-0.37624,0.71142,-0.02436,0.2844,0.3669,-0.13634,0.78991,-0.37624,0.71142,-0.02436,0.2844];
           test = fminunc(@(x)obj.fittingCST(obj,x,orgChord,isUpper,isLower,xyorg,alfa,chord),x0);
            x = linspace(0,1,100)';
            wU = test(5:12);
            wL = [-test(5),test(13:19)];
            yU = obj.CSTClassShape(wU,x,0.5,1,0);
            yL = obj.CSTClassShape(wL,x,0.5,1,0);
            figure(1);clf;hold on;
            plot(modChord(:,1),modChord(:,2));
            plot(x,yL);
            plot(x,yU);
        end


    end        

    methods(Static)

        function res = fittingCST(obj,x,orgChord,isUpper,isLower,xyorg,alfa,chord)
            xyorg = x(1:2);
            alfa = x(3);
            chord = x(4);
            wU = x(5:12);
            wL = [-x(5),x(13:19)];
            %変更上面下面の判定
            modChord = (orgChord-xyorg)*[cos(alfa),sin(alfa);-sin(alfa),cos(alfa)]./chord;
            yU = obj.CSTClassShape(wU,modChord(isUpper,1),0.5,1,0);
            yL = obj.CSTClassShape(wL,modChord(isLower,1),0.5,1,0);
            %canberlineを探す
            res = sum((yU-modChord(isUpper,2)).^2)+sum((yL-modChord(isLower,2)).^2);
            disp(res)
            figure(1);clf;hold on;
            plot(modChord(:,1),modChord(:,2));
            plot(modChord(isUpper,1),yU);
            plot(modChord(isLower,1),yL);
        end
        function [coord] = CSTAirfoil(wl,wu,dz,N)
            %% Function to calculate class and shape function

            % Description : Create a set of airfoil coordinates using CST parametrization method 
            % Input  : wl = CST weight of lower surface
            %          wu = CST weight of upper surface
            %          dz = trailing edge thickness
            % Output : coord = set of x-y coordinates of airfoil generated by CST
            
            
            % Create x coordinate
            x=ones(N+1,1);y=zeros(N+1,1);zeta=zeros(N+1,1);
            for i=1:N+1
                zeta(i)=2*pi/N*(i-1);
                x(i)=0.5*(cos(zeta(i))+1);
            end
            
            % N1 and N2 parameters (N1 = 0.5 and N2 = 1 for airfoil shape)
            N1 = 0.5;
            N2 = 1;
            
            zerind = find(x(:,1) == 0); % Used to separate upper and lower surfaces
            
            xl= x(1:zerind-1); % Lower surface x-coordinates
            xu = x(zerind:end); % Upper surface x-coordinates
            
            [yl] = ClassShape(wl,xl,N1,N2,-dz); % Call ClassShape function to determine lower surface y-coordinates
            [yu] = ClassShape(wu,xu,N1,N2,dz);  % Call ClassShape function to determine upper surface y-coordinates
            
            y = [yl;yu]; % Combine upper and lower y coordinates
            
            coord = [x y]; % Combine x and y into single output
        end
        
        function [y] = CSTClassShape(w,x,N1,N2,dz)
                % Class function; taking input of N1 and N2
                for i = 1:size(x,1)
                    C(i,1) = x(i)^N1*((1-x(i))^N2);
                end
                
                % Shape function; using Bernstein Polynomials
                n = size(w,2)-1; % Order of Bernstein polynomials
                
                for i = 1:n+1
                     K(i) = factorial(n)/(factorial(i-1)*(factorial((n)-(i-1))));
                end
                
                for i = 1:size(x,1)
                    S(i,1) = 0;
                    for j = 1:n+1
                        S(i,1) = S(i,1) + w(j)*K(j)*x(i)^(j-1)*((1-x(i))^(n-(j-1)));
                    end
                end
                
                % Calculate y output
                for i = 1:size(x,1)
                   y(i,1) = C(i,1)*S(i,1) + x(i)*dz;
                end
            end




        function res = softplus(x,beta)
            res = 1./beta * log(1+exp(beta.*x));
        end

        
        function val = valLimiter(val,minVal,maxVal,a)
            sigmoid = @(x)1./(1 + exp(-a.*(x)));
            maxCoef = 1-sigmoid(val - maxVal);
            val = val .* maxCoef + maxVal .* (1-maxCoef);
            minCoef = 1-sigmoid(-val + minVal);
            val = val .* minCoef + minVal .* (1-minCoef);
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
                if not(isempty(obj.prop{propNo}))
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
                dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                
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
                dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                
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
                dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    % fprintf("size of POI.X %d %d\n",size(POI.X));
                    % fprintf("size of c.X %d %d\n",size(c.X));
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                        
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
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0; 
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
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
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

        function [VortexA] = wakeInfluenceMatrix(obj,wakeNo,edgeNo,rowIndex,wakeShape)
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
             for i = 1:size(wakeShape,1)
                if obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) == 1
                    if i == 1
                        wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:);
                        wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:);
                        wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i,:);
                        wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i,:);
                    else
                        wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i-1,:);
                        wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i-1,:);
                        wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i,:);
                        wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i,:);
                    end
                    
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
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
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
                    snorm = obj.matrix_norm(s);
                    Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                    PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                    PB = PA-Al.*smdot;
                    dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                    phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                    phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                    
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
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                        
                        
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
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                        
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
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                        
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
                        snorm = obj.matrix_norm(s);
                        Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
                        PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
                        PB = PA-Al.*smdot;
                        dnom = (PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2);
                        phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),dnom);
                        phiV(abs(dnom) == 0 & abs(PN) < obj.settingUNLSI.pnThreshold * snorm) = 0;
                        
                        VortexA = VortexA+phiV;
                    end
                end
            end

        end
        
        function VelocityA = velocityInfluence(obj,controlPoint)
            %%%%%%%%%%%%%%%%%%%影響係数の計算　パネル⇒パネル%%%%%%%%%%%%%%%
            %rowIndex : 計算する行
            %colIndex : 計算する列
            %Velocity~r : numel(rowIndex)×nbPanelの影響係数
            %Velocity~c : nbPanel×numel(rowIndex)の影響係数
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            verts = obj.tri.Points;
            con  = obj.tri.ConnectivityList(obj.paneltype==1,:);
            center = obj.center(obj.paneltype==1,:);
            normal = obj.orgNormal(obj.paneltype==1,:);

            nbPanel = size(con,1);

            POI.X(:,1) = controlPoint(:,1);
            POI.Y(:,1) = controlPoint(:,2);
            POI.Z(:,1) = controlPoint(:,3);
            
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

            c.X = repmat(c.X,[size(controlPoint,1),1]);
            c.Y = repmat(c.Y,[size(controlPoint,1),1]);
            c.Z = repmat(c.Z,[size(controlPoint,1),1]);
            n.X = repmat(n.X,[size(controlPoint,1),1]);
            n.Y = repmat(n.Y,[size(controlPoint,1),1]);
            n.Z = repmat(n.Z,[size(controlPoint,1),1]);
            N1.X = repmat(N1.X,[size(controlPoint,1),1]);
            N1.Y = repmat(N1.Y,[size(controlPoint,1),1]);
            N1.Z = repmat(N1.Z,[size(controlPoint,1),1]);
            N2.X = repmat(N2.X,[size(controlPoint,1),1]);
            N2.Y = repmat(N2.Y,[size(controlPoint,1),1]);
            N2.Z = repmat(N2.Z,[size(controlPoint,1),1]);
            N3.X = repmat(N3.X,[size(controlPoint,1),1]);
            N3.Y = repmat(N3.Y,[size(controlPoint,1),1]);
            N3.Z = repmat(N3.Z,[size(controlPoint,1),1]);

%             n12.X = (N1.X+N2.X)./2;
%             n12.Y = (N1.Y+N2.Y)./2;
%             n12.Z = (N1.Z+N2.Z)./2;

            pjk.X = POI.X-c.X;
            pjk.Y = POI.Y-c.Y;
            pjk.Z = POI.Z-c.Z;
            %PN = obj.matrix_dot(pjk,n);

            %1回目
            a.X = POI.X-N1.X;
            a.Y = POI.Y-N1.Y;
            a.Z = POI.Z-N1.Z;
            b.X = POI.X-N2.X;
            b.Y = POI.Y-N2.Y;
            b.Z = POI.Z-N2.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            %m = obj.getUnitVector(c,n12);
            %l = obj.matrix_cross(m,n);
            s.X = N2.X-N1.X;
            s.Y = N2.Y-N1.Y;
            s.Z = N2.Z-N1.Z;
            %smdot = obj.matrix_dot(s,m);
            %sldot = obj.matrix_dot(s,l);
            snorm = obj.matrix_norm(s);
            %Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
            %PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
            %PB = PA-Al.*smdot;
            %phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            Vmu = obj.matrix_cross(a,b);
            VelocityA.X = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
            VelocityA.Y = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
            VelocityA.Z = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%             VelocityB.X =  (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%             VelocityB.Y =  (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%             VelocityB.Z =  (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;

            
            %2回目
            a.X = POI.X-N2.X;
            a.Y = POI.Y-N2.Y;
            a.Z = POI.Z-N2.Z;
            b.X = POI.X-N3.X;
            b.Y = POI.Y-N3.Y;
            b.Z = POI.Z-N3.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
            %m = obj.getUnitVector(c,n12);
            %l = obj.matrix_cross(m,n);
            s.X = N3.X-N2.X;
            s.Y = N3.Y-N2.Y;
            s.Z = N3.Z-N2.Z;
            %smdot = obj.matrix_dot(s,m);
%             sldot = obj.matrix_dot(s,l);
            snorm = obj.matrix_norm(s);
%             Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
%             PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
%             PB = PA-Al.*smdot;
            %phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            Vmu = obj.matrix_cross(a,b);
            VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
            VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
            VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%             VelocityB.X = VelocityB.X + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%             VelocityB.Y = VelocityB.Y + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%             VelocityB.Z = VelocityB.Z + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;
            %3回目
            a.X = POI.X-N3.X;
            a.Y = POI.Y-N3.Y;
            a.Z = POI.Z-N3.Z;
            b.X = POI.X-N1.X;
            b.Y = POI.Y-N1.Y;
            b.Z = POI.Z-N1.Z;
            anorm = obj.matrix_norm(a);
            bnorm = obj.matrix_norm(b);
%             m = obj.getUnitVector(c,n12);
%             l = obj.matrix_cross(m,n);
            s.X = N1.X-N3.X;
            s.Y = N1.Y-N3.Y;
            s.Z = N1.Z-N3.Z;
%             smdot = obj.matrix_dot(s,m);
%             sldot = obj.matrix_dot(s,l);
            snorm = obj.matrix_norm(s);
%             Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
%             PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
%             PB = PA-Al.*smdot;
%             phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
            Vmu = obj.matrix_cross(a,b);
            VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
            VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
            VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%             VelocityB.X = VelocityB.X + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%             VelocityB.Y = VelocityB.Y + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%             VelocityB.Z = VelocityB.Z + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;

            %半裁の考慮
            if obj.halfmesh == 1
                POI.Y = -POI.Y;
                pjk.X = POI.X-c.X;
                pjk.Y = POI.Y-c.Y;
                pjk.Z = POI.Z-c.Z;
                %PN = obj.matrix_dot(pjk,n);
                %1回目
                a.X = POI.X-N1.X;
                a.Y = POI.Y-N1.Y;
                a.Z = POI.Z-N1.Z;
                b.X = POI.X-N2.X;
                b.Y = POI.Y-N2.Y;
                b.Z = POI.Z-N2.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
%                 m = obj.getUnitVector(c,n12);
%                 l = obj.matrix_cross(m,n);
                s.X = N2.X-N1.X;
                s.Y = N2.Y-N1.Y;
                s.Z = N2.Z-N1.Z;
%                 smdot = obj.matrix_dot(s,m);
%                 sldot = obj.matrix_dot(s,l);
                snorm = obj.matrix_norm(s);
%                 Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
%                 PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
%                 PB = PA-Al.*smdot;
%                 phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                Vmu = obj.matrix_cross(a,b);
                VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%                 VelocityB.X = VelocityB.X + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%                 VelocityB.Y = VelocityB.Y - (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%                 VelocityB.Z = VelocityB.Z + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;
                %2回目
                a.X = POI.X-N2.X;
                a.Y = POI.Y-N2.Y;
                a.Z = POI.Z-N2.Z;
                b.X = POI.X-N3.X;
                b.Y = POI.Y-N3.Y;
                b.Z = POI.Z-N3.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
%                 m = obj.getUnitVector(c,n12);
%                 l = obj.matrix_cross(m,n);
                s.X = N3.X-N2.X;
                s.Y = N3.Y-N2.Y;
                s.Z = N3.Z-N2.Z;
%                 smdot = obj.matrix_dot(s,m);
%                 sldot = obj.matrix_dot(s,l);
                snorm = obj.matrix_norm(s);
%                 Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
%                 PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
%                 PB = PA-Al.*smdot;
%                 phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                Vmu = obj.matrix_cross(a,b);
                VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%                 VelocityB.X = VelocityB.X + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%                 VelocityB.Y = VelocityB.Y - (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%                 VelocityB.Z = VelocityB.Z + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;
                %3回目
                a.X = POI.X-N3.X;
                a.Y = POI.Y-N3.Y;
                a.Z = POI.Z-N3.Z;
                b.X = POI.X-N1.X;
                b.Y = POI.Y-N1.Y;
                b.Z = POI.Z-N1.Z;
                anorm = obj.matrix_norm(a);
                bnorm = obj.matrix_norm(b);
%                 m = obj.getUnitVector(c,n12);
%                 l = obj.matrix_cross(m,n);
                s.X = N1.X-N3.X;
                s.Y = N1.Y-N3.Y;
                s.Z = N1.Z-N3.Z;
%                 smdot = obj.matrix_dot(s,m);
%                 sldot = obj.matrix_dot(s,l);
                snorm = obj.matrix_norm(s);
%                 Al = obj.matrix_dot(n,obj.matrix_cross(s,a));
%                 PA = obj.matrix_dot(a,obj.matrix_cross(l,obj.matrix_cross(a,s)));
%                 PB = PA-Al.*smdot;
%                 phiV = atan2((smdot.*PN.*(bnorm.*PA-anorm.*PB)),(PA.*PB+PN.^2.*anorm.*bnorm.*(smdot).^2));
                Vmu = obj.matrix_cross(a,b);
                VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
%                 VelocityB.X = VelocityB.X + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.X-sldot.*m.X)+phiV.*n.X;
%                 VelocityB.Y = VelocityB.Y - (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Y-sldot.*m.Y)+phiV.*n.Y;
%                 VelocityB.Z = VelocityB.Z + (log(((anorm+bnorm+snorm)./(anorm+bnorm-snorm)))./snorm).*(smdot.*l.Z-sldot.*m.Z)+phiV.*n.Z;
            end
        end

        function [VelocityA] = wakeVelocityInfluence(obj,wakeNo,edgeNo,controlPoint,wakeShape)
            %%%%%%%%%%%%wakeからパネルへの影響関数%%%%%%%%%%%%%%%%%
            %wakeNoとedgeNoを指定して全パネルへの影響を計算
            %
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            POI.X = controlPoint(:,1);
            POI.Y = controlPoint(:,2);
            POI.Z = controlPoint(:,3);
            nrow = size(controlPoint,1);
            VelocityA.X = zeros(nrow,1);
            VelocityA.Y = zeros(nrow,1);
            VelocityA.Z = zeros(nrow,1);
            for i = 1:size(wakeShape,1)
                if obj.wakeline{wakeNo}.validPanel{edgeNo}(i,1) == 1
                    if i == 1
                        wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:);
                        wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:);
                        wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i,:);
                        wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i,:);
                    else
                        wakepos(1,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i-1,:);
                        wakepos(2,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i-1,:);
                        wakepos(3,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(2,edgeNo),:)+wakeShape(i,:);
                        wakepos(4,:) = obj.tri.Points(obj.wakeline{wakeNo}.validedge(1,edgeNo),:)+wakeShape(i,:);
                    end
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
                    %PN = obj.matrix_dot(pjk,n);
                    
                    %1回目
    %                 a.X = POI.X-Nw1.X;
    %                 a.Y = POI.Y-Nw1.Y;
    %                 a.Z = POI.Z-Nw1.Z;
    %                 b.X = POI.X-Nw2.X;
    %                 b.Y = POI.Y-Nw2.Y;
    %                 b.Z = POI.Z-Nw2.Z;
    %                 s.X = Nw2.X-Nw1.X;
    %                 s.Y = Nw2.Y-Nw1.Y;
    %                 s.Z = Nw2.Z-Nw1.Z;
    %                 anorm = obj.matrix_norm(a);
    %                 bnorm = obj.matrix_norm(b);
    %                 Vmu = obj.matrix_cross(a,b);
    %                 VelocityA.X = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.X;
    %                 VelocityA.Y = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Y;
    %                 VelocityA.Z = (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Z;
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
                    snorm = obj.matrix_norm(s);
                    Vmu = obj.matrix_cross(a,b);
    %                 VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.X;
    %                 VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Y;
    %                 VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Z;
                    VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                    VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                    VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
                    %3回目
    %                 a.X = POI.X-Nw3.X;
    %                 a.Y = POI.Y-Nw3.Y;
    %                 a.Z = POI.Z-Nw3.Z;
    %                 b.X = POI.X-Nw4.X;
    %                 b.Y = POI.Y-Nw4.Y;
    %                 b.Z = POI.Z-Nw4.Z;
    %                 s.X = Nw4.X-Nw3.X;
    %                 s.Y = Nw4.Y-Nw3.Y;
    %                 s.Z = Nw4.Z-Nw3.Z;
    %                 anorm = obj.matrix_norm(a);
    %                 bnorm = obj.matrix_norm(b);
    %                 Vmu = obj.matrix_cross(a,b);
    %                 VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.X;
    %                 VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Y;
    %                 VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Z;
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
                    snorm = obj.matrix_norm(s);
                    Vmu = obj.matrix_cross(a,b);
                    VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                    VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                    VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
    
                    %半裁
                    if obj.halfmesh == 1
                        POI.Y = -POI.Y;
                        pjk.X = POI.X-c.X;
                        pjk.Y = POI.Y-c.Y;
                        pjk.Z = POI.Z-c.Z;
                        %PN = obj.matrix_dot(pjk,n);
                        
                        %1回目
    %                     a.X = POI.X-Nw1.X;
    %                     a.Y = POI.Y-Nw1.Y;
    %                     a.Z = POI.Z-Nw1.Z;
    %                     b.X = POI.X-Nw2.X;
    %                     b.Y = POI.Y-Nw2.Y;
    %                     b.Z = POI.Z-Nw2.Z;
    %                     s.X = Nw2.X-Nw1.X;
    %                     s.Y = Nw2.Y-Nw1.Y;
    %                     s.Z = Nw2.Z-Nw1.Z;
    %                     anorm = obj.matrix_norm(a);
    %                     bnorm = obj.matrix_norm(b);
    %                     Vmu = obj.matrix_cross(a,b);
    %                     VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.X;
    %                     VelocityA.Y = VelocityA.Y - (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Y;
    %                     VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Z;
                        
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
                        snorm = obj.matrix_norm(s);
                        Vmu = obj.matrix_cross(a,b);
                        VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                        VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                        VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
                        %3回目
    %                     a.X = POI.X-Nw3.X;
    %                     a.Y = POI.Y-Nw3.Y;
    %                     a.Z = POI.Z-Nw3.Z;
    %                     b.X = POI.X-Nw4.X;
    %                     b.Y = POI.Y-Nw4.Y;
    %                     b.Z = POI.Z-Nw4.Z;
    %                     s.X = Nw4.X-Nw3.X;
    %                     s.Y = Nw4.Y-Nw3.Y;
    %                     s.Z = Nw4.Z-Nw3.Z;
    %                     anorm = obj.matrix_norm(a);
    %                     bnorm = obj.matrix_norm(b);
    %                     Vmu = obj.matrix_cross(a,b);
    %                     VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.X;
    %                     VelocityA.Y = VelocityA.Y - (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Y;
    %                     VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))).*Vmu.Z;
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
                        snorm = obj.matrix_norm(s);
                        Vmu = obj.matrix_cross(a,b);
                        VelocityA.X = VelocityA.X + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.X;
                        VelocityA.Y = VelocityA.Y + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Y;
                        VelocityA.Z = VelocityA.Z + (anorm+bnorm)./(anorm.*bnorm.*(anorm.*bnorm+obj.matrix_dot(a,b))+obj.settingUNLSI.deltaVortexCore^2.*snorm.*snorm).*Vmu.Z;
                    end
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
        
        function out = evalNTri(C,L1,L2,Ainv)

            mu1 = (C(1)-C(2))/C(3);
            mu2 = (C(3)-C(1))/C(2);
            mu3 = (C(2)-C(3))/C(1);
        
	        % some abbreviations to shorten the following terms
            L3 = 1-L1-L2;
            f13mu1 = 1+3*mu1;
            f13mu2 = 1+3*mu2;
            f13mu3 = 1+3*mu3;
            f1m3mu3 = 1-3*mu3;
            fm13mu2 = -1+3*mu2;
            fm1m3mu3 = -1-3*mu3;
            f1mmu1 = 1-mu1;
            f1mmu2 = 1-mu2;
            f1mmu3 = 1-mu3;
        
            a = 3*f1mmu3*L1-f13mu3*L2+f13mu3*L3;
            b = 3*f1mmu2*L3-f13mu2*L1+f13mu2*L2;
            c = 3*f1mmu1*L2-f13mu1*L3+f13mu1*L1;
            
            x = zeros(9,1);
            x(1,1) = L1;
            x(2,1) = L2;
            x(3,1) = L3;
            x(4,1) = L1*L2;
            x(5,1) = L2*L3;
            x(6,1) = L3*L1;
            x(7,1) = L2*L1^2+0.5*L1*L2*L3*a;
            x(8,1) = L3*L2^2+0.5*L1*L2*L3*c;
            x(9,1) = L1*L3^2+0.5*L1*L2*L3*b;

            out = Ainv'*x;


        end
        
        function out = evalNmTri(C,L1,L2,Ainv)
	        % some abbreviations to shorten the following terms
            L3 = 1-L1-L2;

            out = [L1*eye(2),L2*eye(2),L3*eye(2)];



        end

        function out = evalBTri(C,L1,L2,dphi)

            mu1 = (C(1)-C(2))/C(3);
            mu2 = (C(3)-C(1))/C(2);
            mu3 = (C(2)-C(3))/C(1);
        
	        % some abbreviations to shorten the following terms
            L3 = 1-L1-L2;
            f13mu1 = 1+3*mu1;
            f13mu2 = 1+3*mu2;
            f13mu3 = 1+3*mu3;
            f1m3mu3 = 1-3*mu3;
            fm13mu2 = -1+3*mu2;
            fm1m3mu3 = -1-3*mu3;
            f1mmu1 = 1-mu1;
            f1mmu2 = 1-mu2;
            f1mmu3 = 1-mu3;
        
            a = 3*f1mmu3*L1-f13mu3*L2+f13mu3*L3;
            b = 3*f1mmu2*L3-f13mu2*L1+f13mu2*L2;
            c = 3*f1mmu1*L2-f13mu1*L3+f13mu1*L1;
        
	        % see page 38f of the thesis:
	        % the following terms contains second order derivatives of the 9 shape functions
	        % wrt the triangle coordinates L1 and L2
            out(1,1) = 6 + L2*(-4-2*a) + 4*f1m3mu3*(L2*L3-L1*L2) - 12*L1 + 2*L2*b + 8*(L2*L3-L1*L2);
        
            out(1,2) = -dphi(2,2)*(-2+6*L1+4*L2-L2*b-4*L2*L3+4*L1*L2) ...
                       -dphi(1,2)*(2*L2-L2*a+L2*L3*2*f1m3mu3-L1*L2*2*f1m3mu3);
        
            out(1,3) =  dphi(2,1)*(-2+6*L1+4*L2-L2*b-4*L2*L3+4*L1*L2) ...
                       +dphi(1,1)*(2*L2-L2*a+L2*L3*2*f1m3mu3-L1*L2*2*f1m3mu3);
        
            out(1,4) = -2*L2*c + 4*f13mu1*(L2*L3-L1*L2) - 4*L2 + 2*L2*a + 4*f1m3mu3*(-L2*L3+L1*L2);
        
            out(1,5) = -dphi(1,2)*(2*L2-L2*a+L2*L3*2*f1m3mu3-L1*L2*2*f1m3mu3) ...
                       -dphi(3,2)*(-L2*c+L2*L3*2*f13mu1-L1*L2*2*f13mu1);
        
            out(1,6) = dphi(1,1)*(2*L2-L2*a+L2*L3*2*f1m3mu3-L1*L2*2*f1m3mu3) ...
                      +dphi(3,1)*(-L2*c+L2*L3*2*f13mu1-L1*L2*2*f13mu1);
        
            out(1,7) = -6 + 12*L1 + 8*L2 - 2*L2*b + 8*(L1*L2-L2*L3) + 2*L2*c + 4*f13mu1*(L1*L2-L2*L3);
        
            out(1,8) = -dphi(3,2)*(-L2*c+L2*L3*2*f13mu1-L1*L2*2*f13mu1) ...
                       -dphi(2,2)*(-4+6*L1+4*L2-L2*b-4*L2*L3+4*L1*L2);
        
            out(1,9) = dphi(3,1)*(-L2*c+L2*L3*2*f13mu1-L1*L2*2*f13mu1) ...
                      +dphi(2,1)*(-4+6*L1+4*L2-L2*b-4*L2*L3+4*L1*L2);
        
            out(2,1) = -2*L1*a + 2*L1*L3*2*fm1m3mu3 - 2*L1*L2*2*fm1m3mu3 - 4*L1+2*L1*b - 2*L1*L3*2*fm13mu2 + 2*L1*L2*2*fm13mu2;
        
            out(2,2) = -dphi(2,2)*(2*L1-1*L1*b+1*L1*L3*2*fm13mu2-1*L1*L2*2*fm13mu2) ...
                       -dphi(1,2)*(-1*L1*a+1*L1*L3*2*fm1m3mu3-1*L1*L2*2*fm1m3mu3);
        
            out(2,3) = dphi(2,1)*(2*L1-1*L1*b+1*L1*L3*2*fm13mu2-1*L1*L2*2*fm13mu2) ...
                      +dphi(1,1)*(-1*L1*a+1*L1*L3*2*fm1m3mu3-1*L1*L2*2*fm1m3mu3);
        
            out(2,4) = 6 - 12*L2 - 4*L1-2*L1*c + 8*L3*L1 - 8*L1*L2 + 2*L1*a - 2*L1*L3*2*fm1m3mu3 + 2*L1*L2*2*fm1m3mu3;
        
            out(2,5) = -dphi(1,2)*(-1*L1*a+1*L1*L3*2*fm1m3mu3-1*L1*L2*2*fm1m3mu3) ...
                       -dphi(3,2)*(-6*L2+2-2*L1-1*L1*c+4*L3*L1-4*L1*L2);
        
            out(2,6) = dphi(1,1)*(-1*L1*a+1*L1*L3*2*fm1m3mu3-1*L1*L2*2*fm1m3mu3) ...
                      +dphi(3,1)*(-6*L2+2-2*L1-1*L1*c+4*L3*L1-4*L1*L2);
        
            out(2,7) = -6 + 8*L1 - 2*L1*b + 2*L1*L3*2*fm13mu2 - 2*L1*L2*2*fm13mu2 + 12*L2 + 2*L1*c - 8*L3*L1 +  8*L1*L2;
        
            out(2,8) = -dphi(3,2)*(-6*L2+4-2*L1-1*L1*c+4*L3*L1-4*L1*L2) ...
                       -dphi(2,2)*(2*L1-1*L1*b+1*L1*L3*2*fm13mu2-1*L1*L2*2*fm13mu2);
        
            out(2,9) = dphi(3,1)*(-6*L2+4-2*L1-1*L1*c+4*L3*L1-4*L1*L2) ...
                      +dphi(2,1)*(2*L1-1*L1*b+1*L1*L3*2*fm13mu2-1*L1*L2*2*fm13mu2);
        
            out(3,1) = 2 - 4*L1 + L3*a - L2*a + L2*L3*2*fm1m3mu3 - L1*a - L1*L2*2*fm1m3mu3 + L1*L3*2*f1m3mu3 - L1*L2*2*f1m3mu3 ...
                         - 4*L2 - L3*b + L2*b - L2*L3*2*fm13mu2  + L1*b + L1*L2*2*fm13mu2  + 4*L3*L1         - 4*L1*L2;
        
            out(3,2) = -dphi(2,2)*(-1 + 4*L1 + 2*L2 + 0.5*L3*b - 0.5*L2*b + 0.5*L2*L3*2*fm13mu2 ...
                                     - 0.5*L1*b - 0.5*L1*L2*2*fm13mu2 - 2*L3*L1 + 2*L1*L2) ...
                       -dphi(1,2)*(2*L1 + 0.5*L3*a - 0.5*L2*a + 0.5*L2*L3*2*fm1m3mu3 - 0.5*L1*a ...
                                     - 0.5*L1*L2*2*fm1m3mu3 + 0.5*L1*L3*2*f1m3mu3 - 0.5*L1*L2*2*f1m3mu3);
        
            out(3,3) =  dphi(2,1)*(-1 + 4*L1 + 2*L2 + 0.5*L3*b - 0.5*L2*b + 0.5*L2*L3*2*fm13mu2 ...
                                     - 0.5*L1*b - 0.5*L1*L2*2*fm13mu2 - 2*L3*L1 + 2*L1*L2) ...
                       +dphi(1,1)*(2*L1 + 0.5*L3*a - 0.5*L2*a + 0.5*L2*L3*2*fm1m3mu3 - 0.5*L1*a ...
                                     - 0.5*L1*L2*2*fm1m3mu3 + 0.5*L1*L3*2*f1m3mu3 - 0.5*L1*L2*2*f1m3mu3);
        
            out(3,4) = 2 - 4*L2 + L3*c - L2*c + 4*L2*L3 - L1*c - 4*L1*L2 + L1*L3*2*f13mu1 - L1*L2*2*f13mu1 ...
                         - 4*L1 - L3*a + L2*a + L1*a - L2*L3*2*fm1m3mu3 + L1*L2*2*fm1m3mu3 - L1*L3*2*f1m3mu3 ...
                         + L1*L2*2*f1m3mu3;
        
            out(3,5) = -dphi(1,2)*(2*L1 ...
                           +0.5*L3*a ...
                           -0.5*L2*a ...
                           +0.5*L2*L3*2*fm1m3mu3 ...
                           -0.5*L1*a ...
                           -0.5*L1*L2*2*fm1m3mu3 ...
                           +0.5*L1*L3*2*f1m3mu3 ...
                           -0.5*L1*L2*2*f1m3mu3 ...
                           -1) ...
                     -dphi(3,2)*(-2*L2 ...
                           +0.5*L3*c ...
                           -0.5*L2*c ...
                           +2*L2*L3 ...
                           -0.5*L1*c ...
                           -2*L1*L2 ...
                           +0.5*L1*L3*2*f13mu1 ...
                           -0.5*L1*L2*2*f13mu1 ...
                           );
        
            out(3,6) = dphi(1,1)*(2*L1 ...
                          +0.5*L3*a ...
                          -0.5*L2*a ...
                          +0.5*L2*L3*2*fm1m3mu3 ...
                          -0.5*L1*a ...
                          -0.5*L1*L2*2*fm1m3mu3 ...
                          +0.5*L1*L3*2*f1m3mu3 ...
                          -0.5*L1*L2*2*f1m3mu3 ...
                          -1) ...
                     +dphi(3,1)*(-2*L2 ...
                           +0.5*L3*c ...
                           -0.5*L2*c ...
                           +2*L2*L3 ...
                           -0.5*L1*c ...
                           -2*L1*L2 ...
                           +0.5*L1*L3*2*f13mu1 ...
                           -0.5*L1*L2*2*f13mu1 ...
                           );
        
            out(3,7) = -4 ...
                     +8*L1 ...
                     +8*L2 ...
                     +L3*b ...
                     -L2*b ...
                     +L2*L3*2*fm13mu2 ...
                     -L1*b ...
                     -L1*L2*2*fm13mu2 ...
                     -4*L3*L1 ...
                     +8*L1*L2 ...
                     -L3*c ...
                     +L2*c ...
                     -4*L2*L3 ...
                     +L1*c ...
                     -L1*L3*2*f13mu1 ...
                     +L1*L2*2*f13mu1;
        
            out(3,8) = -dphi(3,2)*(-2*L2 ...
                           +0.5*L3*c ...
                           -0.5*L2*c ...
                           +2*L2*L3 ...
                           -0.5*L1*c ...
                           -2*L1*L2 ...
                           +0.5*L1*L3*2*f13mu1 ...
                           -0.5*L1*L2*2*f13mu1 ...
                           +1) ...
                     -dphi(2,2)*(-2  ...
                           +4*L1 ...
                           +2*L2 ...
                           +0.5*L3*b ...
                           -0.5*L2*b ...
                           +0.5*L2*L3*2*fm13mu2 ...
                           -0.5*L1*b ...
                           -0.5*L1*L2*2*fm13mu2 ...
                           -2*L3*L1 ...
                           +2*L1*L2 ...
                           );
        
            out(3,9) = dphi(3,1)*(-2*L2  ...
                          +0.5*L3*c ...
                          -0.5*L2*c ...
                          +2*L2*L3 ...
                          -0.5*L1*c ...
                          -2*L1*L2 ...
                          +0.5*L1*L3*2*f13mu1 ...
                          -0.5*L1*L2*2*f13mu1 ...
                          +1 ...
                          ) ...
                    +dphi(2,1)*(-2 ...
                          +4*L1 ...
                          +2*L2 ...
                          +0.5*L3*b ...
                          -0.5*L2*b ...
                          +0.5*L2*L3*2*fm13mu2 ...
                          -0.5*L1*b ...
                          -0.5*L1*L2*2*fm13mu2 ...
                          -2*L3*L1 ...
                          +2*L1*L2 ...
                         );
	        % the last row of the matrix must be multipled by 2 (this way, the upper terms gets a bit shorter ...)
            for i=1:9
                out(3,i) = out(3,i)*2.0;
            end

end

        function res = pmsolve(M1,M2,nu,kappa)
            M1nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M1^2-1)))-atan(sqrt(M1^2-1));
            M2nu = M1nu+nu;
            res = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M2^2-1)))-atan(sqrt(M2^2-1))-M2nu;
        end

        function R = calcRMat(X1,X2)
            R = sum((X1 - permute(X2, [3, 2, 1])).^2, 2);
            R = sqrt(squeeze(R));
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
            elseif rbfMode == 3
                phi = @(r,r0)obj.phi3(r,r0);
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
            %%%%%%%%%%%%%%%%%%ここ大丈夫？
            r = real(sqrt(H'-2.*X'*X+H));
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
        
        function [pp] = RbfinvRMake(obj,xd,rbfMode,r0)
            if rbfMode ==4
                phi = @(r,r0)obj.phi4(r,r0);
            elseif rbfMode == 1
                phi = @(r,r0)obj.phi1(r,r0);
            elseif rbfMode == 2
                phi = @(r,r0)obj.phi2(r,r0);
            elseif rbfMode == 3
                phi = @(r,r0)obj.phi3(r,r0);
            else
                error('未実装')
            end
        
            for i = 1:size(xd,2)
                scaleShift(i,1) = min(xd(:,i));
                colmax(i,1) = max(xd(:,i));
                scaleWeight(i,1) = colmax(i,1)-scaleShift(i,1);
                if scaleWeight(i,1) == 0
                    scaleWeight(i,1) = 1;
                end 
            end
            %xdのスケールの変更
            for i = 1:size(xd,2)
                xd(:,i) = (xd(:,i)-scaleShift(i,1))./scaleWeight(i,1);
            end
            
            %X = xd';
            %H = sum(xd.^2,2);
            %H = repmat(H,[1,size(X,2)]);
            r = obj.calcRMat(xd,xd);
            a = phi(r,r0);
            pp.invR = pinv(a);
            pp.rbfMode = rbfMode;
            pp.nSample = size(xd,1);
            pp.nDesign = size(xd,2);
            pp.val_samp = xd;
            pp.R0 = r0;
            pp.scaleShift = scaleShift;
            pp.scaleWeight = scaleWeight;
        end
        
        function [pp] = invRppMake(pp,fd)
            pp.w = pp.invR*fd;
            pp.res_samp = fd;
        end
    
        function phi = phi1(r,r0)
            phi = sqrt(r.*r+r0.*r0);
        end
        
        function phi = phi2(r,r0)
            phi = 1./sqrt(r.*r+r0.*r0);
        end
        
        function phi = phi3(r,r0)
            phi = r.^2.*log(r./r0);
            phi(r<=0) = 0;
        end
        
        function phi = phi4(r,r0)
            phi = exp(-0.5.*r.^2/r0.^2);
        end

        function phi = phi5(r,r0)
            phi = (1-r./r0).^4.*(4.*r./r0+1);
        end

        function fi = execRbfInterp(obj,pp,val_interp)
            if pp.rbfMode ==4
                phi = @(r,r0)obj.phi4(r,r0);
            elseif pp.rbfMode == 1
                phi = @(r,r0)obj.phi1(r,r0);
            elseif pp.rbfMode == 2
                phi = @(r,r0)obj.phi2(r,r0);
            elseif pp.rbfMode == 3
                phi = @(r,r0)obj.phi3(r,r0);
            else
                error('未実装')
            end
            nSamp = size(pp.val_samp,1);
            nInterp = size(val_interp,1);
            val_interp= (val_interp-repmat(pp.scaleShift(:)',[nInterp,1]))./repmat(pp.scaleWeight(:)',[nInterp,1]);
            %Xi = sum(val_interp.^2,2);
            %H1 = repmat(Xi,[1,nSamp]);
            %Xs = sum(pp.val_samp.^2,2);
            %H2 = repmat(Xs',[nInterp,1]);
            %M = val_interp*pp.val_samp';
            %r = sqrt(H1-2.*M+H2);
            r = obj.calcRMat(pp.val_samp,val_interp);
            fi = phi(r,pp.R0)'*pp.w;
        end

        function [intMatrix, intSurface] = SurfaceIntersection(obj,surface1, surface2, varargin)
        %SURFACEINTERSECTION intersection of 2 surfaces
        % [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2)
        % calculates the intersection of surfaces 1 and 2. Code can either return
        % just the matrix indicating which face of surface1 intersected with face
        % of surface2, which is calculated using Tomas Moller algorithm, or can
        % also return the actual line of intersection. In case when parts of the
        % surface 1 and 2 lay on the same plane the intersection is a 2D area
        % instead of 1D edge. In such a case the intersection area will be
        % triangulated and intSurface.edges will hold the edges of the
        % triangulation surface and intSurface.faces will hold the faces.
        %
        % INPUT:
        %  * surface1 & surface2 - two surfaces defined as structs or classes.
        %    Several inputs are possible:
        %    - struct with "faces" and "vertices" fields
        %    - 'triangulation' class (only the boundary surface will be used)
        %    - 'delaunayTriangulation' class
        %
        % OUTPUT:
        % * intMatrix - sparse Matrix with n1 x n2 dimension where n1 and n2 are
        %               number of faces in surfaces
        % * intSurface - a structure with following fields:
        %     intSurface.vertices - N x 3 array of unique points
        %     intSurface.edges    - N x 2 array of edge vertex ID's
        %     intSurface.faces    - N x 3 array of face vertex ID's
        %
        % ALGORITHM:
        % Based on Triangle/triangle intersection test routine by Tomas Mler, 1997.
        %  See article "A Fast Triangle-Triangle Intersection Test",
        %  Journal of Graphics Tools, 2(2), 1997
        %  http://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf
        %  http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt
        
        %% Get FACES and VERTICES inputs
        if isa(surface1, 'triangulation')
          [surface1.faces, surface1.vertices] = freeBoundary(surface1);
        elseif isa(surface1, 'delaunayTriangulation')
          S = surface1;
          surface1 = [];
          surface1.faces    = S.ConnectivityList;
          surface1.vertices = S.Points;
          clear S
        end
        if isa(surface2, 'triangulation')
          [surface2.faces, surface1.vertices] = freeBoundary(surface2);
        elseif isa(surface2, 'delaunayTriangulation')
          S = surface2;
          surface2 = [];
          surface2.faces    = S.ConnectivityList;
          surface2.vertices = S.Points;
          clear S
        end
        ok1 = isstruct(surface1) && isfield(surface1, 'vertices') && isfield(surface1, 'faces');
        ok2 = isstruct(surface2) && isfield(surface2, 'vertices') && isfield(surface2, 'faces');
        assert(ok1, 'Surface #1 must be a struct with "faces" and "vertices" fields' );
        assert(ok2, 'Surface #2 must be a struct with "faces" and "vertices" fields' );
        
        %% Flip dimentions if necessery
        if size(surface1.faces,1)==3 && size(surface1.faces,2)~=3
          surface1.faces = surface1.faces';
        end
        if size(surface1.vertices,1)==3 && size(surface1.vertices,2)~=3
          surface1.vertices = surface1.vertices';
        end
        if size(surface2.faces,1)==3 && size(surface2.faces,2)~=3
          surface2.faces = surface2.faces';
        end
        if size(surface2.vertices,1)==3 && size(surface2.vertices,2)~=3
          surface2.vertices = surface2.vertices';
        end
        
        %% Parse extra parameters
        getIntersection = (nargout>1);
        debug = false;
        PointRoundingTol = 1e6;
        algorithm = 'moller';
        k=1;
        nVarargs = length(varargin);
        while (k<=nVarargs)
          assert(ischar(varargin{k}), 'Incorrect input parameters')
          switch lower(varargin{k})
            case 'debug'
              debug = varargin{k+1}~=0;
              k = k+1;
            case 'algorithm'
              algorithm = lower(strtrim(varargin{k+1}));
              k = k+1;
            case 'pointroundingtol'
              PointRoundingTol = varargin{k+1};
              k = k+1;
          end
          k = k+1;
        end
        
        %% Initialize variables
        epsilon = eps;
        nFace1 = size(surface1.faces,1);
        nFace2 = size(surface2.faces,1);
        nVert1 = size(surface1.vertices,1);
        nVert2 = size(surface2.vertices,1);
        
        %% create strip down versions of MATLAB cross and dot function
        cross_prod = @(a,b) [...
          a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
          a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
          a(:,1).*b(:,2)-a(:,2).*b(:,1)];
        dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
        normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));
        
        %% Initialize output variables
        % intersect is a nFace1 x nFace2 matrix. Possible values: -2 (do not know),
        % -1 (coplanar with unknown overlap), 0 (no intersections), 1 (intersects).
        % Negative values are internal only.
        intMatrix  = zeros([nFace1,nFace2], 'int8')-2; % -2 indicates that there was no succesful test yet
        intSurface.vertices = [];
        intSurface.faces    = [];
        intSurface.edges    = [];
        
        % =======================================================================
        %% === Stage 1 ==========================================================
        % =======================================================================
        % Each triangle is a subset of the plane it lies in, so for two triangles
        % to intersect they must overlap along the line of intersection of their
        % planes. Hence, a necessary condition for intersection is that each
        % triangle must intersect the plane of the other.
        % Mler痴 method begins by checking the mutual intersection of each
        % triangle with the plane of the other. To do so, it determines for each
        % triangle on which side of the other triangle痴 supporting plane its
        % vertices lie. Now, if all vertices of one triangle lie on the same side
        % and no vertex is on the plane, the intersection is rejected.
        
        %% compute plane equations for each triangle of the surface #1
        % plane equation #1: N1.X-d1=0
        V1 = surface1.vertices(surface1.faces(:,1),:);
        V2 = surface1.vertices(surface1.faces(:,2),:);
        V3 = surface1.vertices(surface1.faces(:,3),:);
        N1 = cross_prod(V2-V1,V3-V1); % array size nFace1 x 3
        N1 = normalize(N1);
        d1 = dot_prod(N1,V1);         % array size nFace1 x 1
        
        %% Distance from surface #2 vertices to planes of surface #1
        % Calculate signed distance from all vertices of surface #2 to each plane
        % of of surface #1
        du = zeros(nFace1,nVert2);
        for iVert2 = 1:nVert2
          p = surface2.vertices(iVert2,:);
          du(:,iVert2) = N1(:,1)*p(1) + N1(:,2)*p(2) + N1(:,3)*p(3) - d1;
        end
        if debug
          assert(all(size(du)==[nFace1,nVert2]), 'Incorrect array dimensions: dv')
        end
        du(abs(du)<epsilon)=0; % robustness check
        % Distances from vertex 1, 2 & 3 of faces of surface #2 to planes of surface #1
        du1 = du(:,surface2.faces(:,1));
        du2 = du(:,surface2.faces(:,2));
        du3 = du(:,surface2.faces(:,3));
        if debug
          assert(all(size(du1)==size(intMatrix)), 'Incorrect array dimensions: du1')
        end
        clear du
        intMatrix(du1.*du2>0 & du1.*du3>0) = 0;   % same sign on all of them & not equal 0
        if(all(intMatrix==0)), return; end        % no intersections
        intMatrix(du1==0 & du2==0 & du3==0) = -1; % coplanar with unknown overlap
        
        %% compute plane of triangle (U0,U1,U2)
        % plane equation 2: N2.X-d2=0
        U1 = surface2.vertices(surface2.faces(:,1),:);
        U2 = surface2.vertices(surface2.faces(:,2),:);
        U3 = surface2.vertices(surface2.faces(:,3),:);
        N2 = cross_prod(U2-U1,U3-U1); % array size nFace1 x 3
        N2 = normalize(N2);
        d2 = dot_prod(N2,U1);        % array size nFace1 x 1
        
        %% Distance from surface #1 vertices to planes of surface #2
        % Calculate signed distance from all vertices of surface #1 to each plane
        % of of surface #2
        dv = zeros(nFace2,nVert1);
        for iVert1 = 1:nVert1
          p = surface1.vertices(iVert1,:);
          dv(:,iVert1) = N2(:,1)*p(1) + N2(:,2)*p(2) + N2(:,3)*p(3) - d2;
        end
        if debug
          assert(all(size(dv)==[nFace2,nVert1]), 'Incorrect array dimensions: dv')
        end
        dv(abs(dv)<epsilon)=0; % robustness check
        % Distances from vertex 1, 2 & 3 of faces of surface #1 to planes of surface #2
        dv1 = dv(:,surface1.faces(:,1))';
        dv2 = dv(:,surface1.faces(:,2))';
        dv3 = dv(:,surface1.faces(:,3))';
        if debug
          assert(all(size(dv1)==size(intMatrix)), 'Incorrect array dimensions: dv1')
        end
        clear dv
        intMatrix(dv1.*dv2>0 & dv1.*dv3>0) = 0;   % same sign on all of them & not equal 0
        if(all(intMatrix==0)), return; end        % no intersections
        intMatrix(dv1==0 & dv2==0 & dv3==0) = -1; % coplanar with unknown overlap
        
        % =======================================================================
        %% === Stage 2 ==========================================================
        % =======================================================================
        
        %% Process remaining (non-coplanar) triangle pairs
        tMsk = (intMatrix==-2);
        n = nnz(tMsk);
        if n>0
          [face1, face2] = find(tMsk);
          switch lower(algorithm)
            case 'moller'
              if size(dv1(tMsk),1)==1
                dv = [dv1(tMsk)', dv2(tMsk)', dv3(tMsk)'];
                du = [du1(tMsk)', du2(tMsk)', du3(tMsk)'];
              else
                dv = [dv1(tMsk), dv2(tMsk), dv3(tMsk)];
                du = [du1(tMsk), du2(tMsk), du3(tMsk)];
              end
              
              [intMatrix(tMsk), intSurface] = obj.TriangleIntersection3D_Moller(obj,...
                V1(face1,:), V2(face1,:), V3(face1,:), N1(face1,:), d1(face1,:), dv, ...
                U1(face2,:), U2(face2,:), U3(face2,:), N2(face2,:), d2(face2,:), du, ...
                getIntersection, debug);
            otherwise
              error('Unknown algorithm name');
          end
        end % if
        
        %% Process coplanar triangle pairs. Pass #1:
        % compare the overlap of the bounding boxes
        tMsk = (intMatrix==-1);
        if nnz(tMsk)>0
          [face1, face2] = find(tMsk);
          overlap = true;
          for idim = 1:3
            v = [V1(face1,idim), V2(face1,idim), V3(face1,idim)];
            u = [U1(face2,idim), U2(face2,idim), U3(face2,idim)];
            t1 = min(v,[],2);
            t2 = max(v,[],2);
            s1 = min(u,[],2);
            s2 = max(u,[],2);
            overlap = overlap & (s1<=t2 & t1<=s2);
          end
          % if overlap intMatrix will remain "-1" otherwise it will change to "0"
          intMatrix(tMsk) = -1*overlap;
          clear v u t1 t2 s1 s2 overlap
        end
        
        %% Process coplanar triangle pairs. Pass #2:
        % use edge-edge intersections
        tMsk = (intMatrix==-1);
        if nnz(tMsk)>0
          [face1, face2] = find(tMsk);
          
          % repack data prior to function call
          V(:,:,1)=V1(face1,:); V(:,:,2)=V2(face1,:); V(:,:,3)=V3(face1,:);
          U(:,:,1)=U1(face2,:); U(:,:,2)=U2(face2,:); U(:,:,3)=U3(face2,:);
          [intMatrix(tMsk), intSurface2] = obj.TriangleIntersection2D(obj,V, U, ...
            N1(face1,:), getIntersection, debug);
          
          %% Merge surfaces
          if getIntersection
            np = size(intSurface.vertices,1);
            intSurface.vertices = [intSurface.vertices; intSurface2.vertices];
            intSurface.faces    = [intSurface.faces;    intSurface2.faces+np];
            intSurface.edges    = [intSurface.edges;    intSurface2.edges+np];
            if debug
              np = size(intSurface.vertices,1);
              assert(max(intSurface.faces(:))<=np, 'Bad surface definition')
              assert(max(intSurface.edges(:))<=np, 'Bad surface definition')
            end
          end
        end
        
        %% Clean up the outputs
        intMatrix = sparse(double(intMatrix));
        if(getIntersection)
          % make point array unique
          P = round(intSurface.vertices*PointRoundingTol)/PointRoundingTol;
          [~,ia,ic] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
          intSurface.vertices = intSurface.vertices(ia,:);
          intSurface.faces = ic(intSurface.faces);
          intSurface.edges = ic(intSurface.edges);
        end
        end % function
     
        function [iMsk, intSurface] = TriangleIntersection3D_Moller(obj,...
          V1, V2, V3, N1, d1, dv, ...
          U1, U2, U3, N2, d2, du, ...
          getIntersection, debug)
        %TriangleIntersection3D tests if 2 triangles defined in 3D intersect.
        % This is a secondary test following Tomas Moller algorithm
        %
        % INPUTS:
        %   V1, V2, V3, - Nx3 array of surface 1 triangle vertex coordinates
        %   U1, U2, U3, - Nx3 array of surface 2 triangle vertex coordinates
        %   N1, d1      - Nx3 array of surface 1 triangle plane equations N1.X-d1=0
        %   N2, d2      - Nx3 array of surface 2 triangle plane equations N2.X-d2=0
        %   dv          - Nx3 array of distances of surface 1 triangle vertices to surface 2 planes
        %   du          - Nx3 array of distances of surface 2 triangle vertices to surface 1 planes
        %   getIntersection - do we need to output the intersecting surface?
        %      Algorithm is much simpler if we do not.
        %   debug       - In the debugging mode much more extra "sanity check" test
        %      are performed.
        %
        % OUTPUT:
        %   iMsk - N x 1 intersection boolean mask marking which triangles overlap
        %   intSurface - intersection surface
        %
        % ALGORITHM:
        % The input triangles are guaranteed to intersect the line of intersection
        % of the two planes. Furthermore, these intersections form intervals on
        % this line, and the triangles overlap iff these intervals overlap as well.
        % Hence, the last part of  the algorithm computes a parametric equation
        % L(t) of the line of intersection of the two planes, finds the intervals
        % (i.e. scalar intervals on L(t)) for which the line lies inside each
        % triangle and performs a one-dimensional interval overlap test.
        if debug
          ok = size(N1,2)==3 && size(N2,2)==3 && size(dv,2)==3 && size(du,2)==3 && ...
            size(V1,2)==3 && size(V2,2)==3 && size(V3,2)==3 && ...
            size(U1,2)==3 && size(U2,2)==3 && size(U3,2)==3;
          assert(ok, 'Incorrect array dimensions');
        end
        
        %% create strip down versions of MATLAB cross and dot function
        cross_prod = @(a,b) [...
          a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
          a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
          a(:,1).*b(:,2)-a(:,2).*b(:,1)];
        dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
        normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));
        
        %% Find intervals of surface 1 and 2 triangles
        % compute the scalar intervals on L(t) for which the line lies inside each
        % triangle
        
        % Plane creates two open half-spaces. Find the odd vertex, which:
        % 1) if no or two vertices are on the plane than pick the vertex which is
        %    by itself in its half-space
        % 2) if one vertex is on the plane and the other two occupy the same
        %    half-space than pick the vertex on the plane
        % 3) if one vertex is on the plane and the other two occupy different
        %    half-spaces than pick one of the vertices off the plane
        % Find vertex using a look-up table "lut" with key calculated based on
        % sign of dv and du arrays
        lut = [0;3;3;2;1;3;2;2;1;1;2;3;3;0;3;3;2;1;1;2;2;3;1;2;3;3;0];
        n = numel(d1);
        rows = (1:n)';
        
        %% order surface 1 triangle vertices
        a1 = lut(sign(dv)*[9; 3; 1] + 14); % calculate the key and call the look-up table
        [b1, c1] = obj.otherDim(a1);
        if debug
          assert(all(a1>0), 'Something Wrong: triangles are coplanar')
        end
        a1 = sub2ind([n,3],rows,a1); % convert row and column IDs to array indecies
        b1 = sub2ind([n,3],rows,b1);
        c1 = sub2ind([n,3],rows,c1);
        
        %% order surface 2 triangle vertices
        a2 = lut(sign(du)*[9; 3; 1] + 14); % calculate the key and call the look-up table
        [b2, c2] = obj.otherDim(a2);
        if debug
          assert(all(a2>0), 'Something Wrong: triangles are coplanar')
        end
        a2 = sub2ind([n,3],rows,a2);
        b2 = sub2ind([n,3],rows,b2);
        c2 = sub2ind([n,3],rows,c2);
        
        %% compute direction of L the line of intersection of 2 planes
        % containing 2 triangles. Line L parametric equation: t*D+O=0
        D = cross_prod(N1,N2);    % D must be perpendicular to both N1 and N2
        [~, maxDim] = max(abs(D),[],2); % compute and index to the largest component of D
        if(getIntersection)
          D = normalize(D);
          O = zeros(n,3);
          d = [d1, d2, zeros(n,1)];
          for r =1:n
            N = [N1(r,:); N2(r,:); 0, 0, 0];
            N(3,maxDim(r)) = 1;
            dd = d(r,:)';
            O(r,:) = (N\dd)'; %Solve systems of linear equations N*D3 = d for D3
          end
          clear N d dd
        end
        
        %% projection of triangle(V1,V2,V3) and triangle(U1,U2,U3) onto intersection line
        % Vp and Up are Nx3 arrays with columns indicating corners of triangles 1 and 2
        if(getIntersection)
          Vp=[dot_prod(V1-O,D), dot_prod(V2-O,D), dot_prod(V3-O,D)];
          Up=[dot_prod(U1-O,D), dot_prod(U2-O,D), dot_prod(U3-O,D)];
        else
          % Project on one of the axis (closest to the intersection line) instead.
          % Simplified projection is faster and sufficient if we do not need
          % intersection line
          idx = sub2ind([n,3],rows,maxDim);
          Vp = [V1(idx), V2(idx), V3(idx)];
          Up = [U1(idx), U2(idx), U3(idx)];
        end
        clear V1 V2 V3 U1 U2 U3
        
        %% Calculate surface 1 and 2 triangle intervals
        % t1 and t2 are intersection points of surface 1 with the intersection line
        % t*D+O=0, and s1 & s2 are intersection points of surface 2 with the same
        % line. Tomas Moller algorithm made this section much more complicated
        % trying to avoid divisions. However, I could not detect any speed-up.
        % Operations (ADD: 12; MUL:4 ; DIV:4 )
        t1 = Vp(a1) - (Vp(b1)-Vp(a1)).*dv(a1)./(dv(b1)-dv(a1));
        t2 = Vp(a1) - (Vp(c1)-Vp(a1)).*dv(a1)./(dv(c1)-dv(a1));
        s1 = Up(a2) - (Up(b2)-Up(a2)).*du(a2)./(du(b2)-du(a2));
        s2 = Up(a2) - (Up(c2)-Up(a2)).*du(a2)./(du(c2)-du(a2));
        
        %% Order the intervals as to t1<t2 and s1<s2
        msk = t2<t1; % order t1 and t2 so t1<t2
        t = t1(msk); t1(msk)=t2(msk); t2(msk)=t; % swap
        msk = s2<s1; % order s1 and s2 so s1<s2
        t = s1(msk); s1(msk)=s2(msk); s2(msk)=t; % swap
        
        %% Perform THE final test we were preparying for.
        % It test for the overlap of 2 1D intervals s1->s2 and t1->t2
        iMsk = (s1<t2 & t1<s2);
        
        %% calculate intersection segments
        n = nnz(iMsk);
        if(getIntersection && n>0)
          % p1 = D*max(t1,s1) + O;    p2 = D*min(t2,s2) + O
          p1 = bsxfun(@times,D(iMsk,:),max(t1(iMsk),s1(iMsk))) + O(iMsk,:);
          p2 = bsxfun(@times,D(iMsk,:),min(t2(iMsk),s2(iMsk))) + O(iMsk,:);
          intSurface.vertices = [p1; p2];
          intSurface.faces    = [1:n; n+1:2*n; n+1:2*n]';
          intSurface.edges    = intSurface.faces(:,1:2);
        else
          intSurface.vertices = [];
          intSurface.faces    = [];
          intSurface.edges    = [];
        end % if
        end % function
           
        function [overlap, intSurface] = TriangleIntersection2D(obj,V, U, N, ...
          getIntersection, debug)
        % Triangles V(V0,V1,V2) and U(U0,U1,U2) are are coplanar. Do they overlap?
        % INPUTS:
        % N - array(n,3) of surface normals where V(i,:,:) and U(i,:,:) are on the same plane
        % V - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #1 vertices
        % U - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #2 vertices
        %
        % OUTPUT:
        %   iMsk - N x 1 intersection boolean mask marking which triangles overlap
        %   intSurface - intersection surface
        
        
        %  * parameters: vertices of triangle 1: V0,V1,V2
        %  *             vertices of triangle 2: U0,U1,U2
        %  * result    : returns 1 if the triangles intersect, otherwise 0
        
        %% Constants needed for creating a mesh based on 3 to 6 points in a circle
        tri_mesh{6}  = [1 2 6; 2 4 6; 2 3 4; 4 5 6];
        tri_mesh{5}  = [1 2 3; 1 3 4; 4 5 1];
        tri_mesh{4}  = [1 2 3; 1 3 4];
        tri_mesh{3}  = 1:3;
        vertices = [];
        faces    = [];
        pairs    = [];  % each row corresponds to pair of faces. match row number with face number
        nVert    = 0;
        
        %% use edge-edge intersections
        overlap = false(size(N,1),1);
        i1Idx = [1 1 1 2 2 2 3 3 3];
        i2Idx = [3 3 3 1 1 1 2 2 2];
        j1Idx = [1 2 3 1 2 3 1 2 3];
        j2Idx = [3 1 2 3 1 2 3 1 2];
        for row = 1:size(N,1)
          % When it is necesary to project 3D plane on 2D, dIdx will be the optimal
          % dimensions to use.
          [~, a] = max(abs(N(row,:))); 
          [b, c] = obj.otherDim(a); 
          dIdx = [b, c]; 
          order = [];
        
          %% test all edges of triangle 1 against the edges of triangle 2
          % triangles overlap if edges cross
          [edgeMat, P] = obj.EdgesIntersect3D(...
            squeeze(V(row,:,i1Idx))',squeeze(V(row,:,i2Idx))', ...
            squeeze(U(row,:,j1Idx))',squeeze(U(row,:,j2Idx))');
          overlap(row) = any(edgeMat);
          if ~getIntersection && overlap(row), continue; end
          
          if ~overlap(row)
            %% project onto an axis-aligned plane, that maximizes the area
            % of the triangles, compute indices: dIdx which correspond to 2 smallest N1
            % components.
            V2d = [V(row,dIdx,1); V(row,dIdx,2); V(row,dIdx,3)]; % each row is a 2D vertex
            U2d = [U(row,dIdx,1); U(row,dIdx,2); U(row,dIdx,3)];
            
            %% test if tri1 is totally contained in tri2 or vice varsa
            if obj.PointInTriangle2D(V2d(1,:), U2d) % tri1 is totally contained in tri2
              overlap(row) = true;
              order = 1:3;
            elseif obj.PointInTriangle2D(U2d(1,:), V2d) % tri2 is totally contained in tri1
              overlap(row) = true;
              order = 4:6;
            end
            if overlap(row) && ~getIntersection, continue; end
            clear V2d U2d
          end
          
          %% Build the intersection surface
          if getIntersection && overlap(row)
            %Assemble all the points which might be needed for desining
            %intersection polygon: Intersection points and points from triangle 1
            %and 2
            points   = [P(edgeMat,:); squeeze(V(row,:,:))'; squeeze(U(row,:,:))'];
            if isempty(order) % when one tri is totally contained in the other tri then order is set
              order = obj.IntersectionPolygon(obj,edgeMat>0, points, dIdx, debug);
              if isempty(order), continue; end
            end
            nPoint   = length(order);    % how many points will be added?
            nFace    = nPoint-2;         % how many faces will be added?
            vertices = [vertices; points(order,:)]; %#ok<*AGROW>
            faces    = [faces; nVert+tri_mesh{nPoint} ];
            pairs    = [pairs; row+zeros(nFace,1)];  % each row corresponds to pair of faces. match row number with face number
            nVert    = nVert + nPoint;
            if debug
              assert(max(faces(:))<=size(vertices,1), 'Bad surface definition')
            end
          end
        end % for
        
        %% Prepare outputs
        intSurface.vertices = vertices;
        intSurface.faces    = faces;
        if isempty(faces)
          intSurface.edges = [];
        else
          intSurface.edges = [faces(:,1:2); faces(:,2:3); faces(:,[1,3])];
        end
        end % function
               
        function polygon = IntersectionPolygon(obj,edgeMat, points, dIdx, debug)
        % edgeMat is an edge intersection matrix with 3 rows for edges between
        % the points 1-3, 1-2, & 2-3 of the triangle 1 and 3 columns for the same
        % edges of the triangle 2. If 2 edges intersect a point of intersection
        % is calculated and stored in array "points" followed by points of the
        % triangles 1 & 2.  This function calculates the polygon of the intersection
        % between 2 triangles.
        
        persistent orderLUT verified
        if isempty(orderLUT) || isempty(orderLUT{3})
          % This pre-calculated look-up table is used to quickly look up the order of
          % the vertices in array "points" which make up polygon of the intersection
          % between 2 triangles. A unique key is calculated for each edgeMat using
          % dot product between edgeMat(:) and [256 128 64 32 16 8 4 2 1], which is
          % used to look up point order around the polygon. Negative numbers in the
          % LUT indicate values which were not observed yet so they were not
          % independently verified.
          % reshape(sprintf('%09s',dec2base(key, 2)),3,3) will convert from the key
          % to matrix.
          OrderLUT = zeros(432,1);  
          OrderLUT(003) = 127;
          OrderLUT(005) = 128;
          OrderLUT(006) = 126;
          OrderLUT(009) = 124;
          OrderLUT(010) = 1427;
          OrderLUT(012) = 1428;
          OrderLUT(017) = 1427;
          OrderLUT(018) = 124;
          OrderLUT(020) = 1426;
          OrderLUT(024) = 127;
          OrderLUT(027) = 1243;
          OrderLUT(029) = 12438;
          OrderLUT(030) = 12034;
          OrderLUT(033) = 1428;
          OrderLUT(034) = 1426;
          OrderLUT(036) = 124;
          OrderLUT(040) = 128;
          OrderLUT(043) = 21834;
          OrderLUT(045) = 1243;
          OrderLUT(046) = 21349;
          OrderLUT(048) = 126;
          OrderLUT(051) = 12340;
          OrderLUT(053) = 12943;
          OrderLUT(054) = 1243;
          OrderLUT(065) = 125;
          OrderLUT(066) = 1527;
          OrderLUT(068) = 1825;
          OrderLUT(072) = 123;
          OrderLUT(080) = 1327;
          OrderLUT(083) = 15234;
          OrderLUT(085) = -15234;
          OrderLUT(086) = -15243;
          OrderLUT(090) = 13247;
          OrderLUT(092) = -13247;
          OrderLUT(096) = 1328;
          OrderLUT(099) = 152834;
          OrderLUT(101) = 15234;
          OrderLUT(102) = 152349;
          OrderLUT(106) = 132847;
          OrderLUT(108) = 13247;
          OrderLUT(114) = 102347;
          OrderLUT(116) = -13247;
          OrderLUT(129) = 1527;
          OrderLUT(130) = 125;
          OrderLUT(132) = 1526;
          OrderLUT(136) = 1327;
          OrderLUT(139) = 15243;
          OrderLUT(141) = 152438;
          OrderLUT(142) = 152034;
          OrderLUT(144) = 123;
          OrderLUT(153) = 12347;
          OrderLUT(156) = 123047;
          OrderLUT(160) = 1326;
          OrderLUT(163) = -152043;
          OrderLUT(165) = 13247;
          OrderLUT(166) = 15234;
          OrderLUT(169) = -182347;
          OrderLUT(172) = 193247;
          OrderLUT(177) = -132047;
          OrderLUT(180) = 13247;
          OrderLUT(192) = 127;
          OrderLUT(195) = 1243;
          OrderLUT(197) = 12438;
          OrderLUT(198) = 12034;
          OrderLUT(202) = 12364;
          OrderLUT(204) = 123648;
          OrderLUT(209) = 21364;
          OrderLUT(212) = -21364;
          OrderLUT(216) = 1243;
          OrderLUT(225) = -124638;
          OrderLUT(226) = 120364;
          OrderLUT(232) = 12438;
          OrderLUT(238) = 124356;
          OrderLUT(240) = 12034;
          OrderLUT(245) = -214356;
          OrderLUT(257) = 1528;
          OrderLUT(258) = 1526;
          OrderLUT(260) = 125;
          OrderLUT(264) = 1328;
          OrderLUT(267) = -152438;
          OrderLUT(269) = 15243;
          OrderLUT(270) = -152943;
          OrderLUT(272) = 1326;
          OrderLUT(275) = 152340;
          OrderLUT(277) = 152943;
          OrderLUT(278) = 15243;
          OrderLUT(281) = 182347;
          OrderLUT(282) = -103247;
          OrderLUT(288) = 123;
          OrderLUT(297) = 12347;
          OrderLUT(298) = -123947;
          OrderLUT(305) = 123947;
          OrderLUT(306) = 12347;
          OrderLUT(320) = 128;
          OrderLUT(323) = 21834;
          OrderLUT(325) = 1243;
          OrderLUT(326) = 21349;
          OrderLUT(330) = -123648;
          OrderLUT(332) = 12364;
          OrderLUT(337) = 183642;
          OrderLUT(340) = -129364;
          OrderLUT(344) = 21834;
          OrderLUT(350) = -124365;
          OrderLUT(353) = 12463;
          OrderLUT(354) = 136492;
          OrderLUT(360) = 1243;
          OrderLUT(368) = 12943;
          OrderLUT(371) = 126543;
          OrderLUT(384) = 126;
          OrderLUT(387) = 12340;
          OrderLUT(389) = 12943;
          OrderLUT(390) = 1243;
          OrderLUT(394) = -103642;
          OrderLUT(396) = 129364;
          OrderLUT(401) = 123640;
          OrderLUT(404) = 12364;
          OrderLUT(408) = 12340;
          OrderLUT(413) = 215643;
          OrderLUT(417) = -136492;
          OrderLUT(418) = 12463;
          OrderLUT(424) = 13492;
          OrderLUT(427) = -213456;
          OrderLUT(432) = 1342;
          
          % Convert to more convinient format
          orderLUT = cell(size(OrderLUT));
          for i = 1:size(OrderLUT,1)
            polygon = abs(OrderLUT(i));
            if polygon>0
              polygon = num2str(polygon)-48; % Convert from a single number to array of digits
              polygon(polygon==0) = 10;      % 0 stands for 10
              orderLUT{i} = polygon;
            end
          end
          % Negative numbers in the LUT indicate values which were not observed yet
          % so they were not independently verified.
          verified = OrderLUT>0;
          clear OrderLUT
        end
        
        %% Calculate unique key for each edgeMat configuration
        key = dot(1*edgeMat(:)', [256 128 64 32 16 8 4 2 1]);
        assert(key<=432, 'Error: in IntersectionPolygon: key is out of bound');
        
        %% Look up the point order around the polygon
        polygon = orderLUT{key};
        if (isempty(polygon))
          return
        end
        
        %% in a rare case of 2 intersections there is ambiguity if one or two
        % vertices of the triangle lay inside the other triangle. OrderLUT stores
        % only the single vertex cases.
        nx = nnz(edgeMat(:));
        if nx==2
          pList = polygon;       % list of vertices to check
          pList(pList<=nx) = []; % keep only the triangle points of the polygon
          flip = false;    % was there a flip from single vertex to vertices case?
          for ip = 1:length(pList)
            p = pList(ip);                 % point to check
            t = floor((p-nx-1)/3);         % does it belong to triangle 0 or 1 (actually 1 or 2)
            tri = (1:3) + nx + 3*abs(1-t); % Points belonging to the other triangle
            if ~obj.PointInTriangle2D(points(p,dIdx), points(tri,dIdx))
              d = nx+t*3;    % offset
              % "p-d" is vertex number of point just tested: 1, 2, or 3. "b, c" are
              % the other 2 vertices
              [b, c] = obj.otherDim(p-d);
              polygon = [polygon(polygon~=p), b+d, c+d]; % remove i2 and add i0 and i1
              flip = true;
            end
          end
          if flip
            % if ther were any flips than use existing codes to figure out the
            % order of the points around the polygon
            DT = delaunayTriangulation(points(polygon,dIdx));
            idx = freeBoundary(DT)';
            idx(2,:) = [];
            polygon = polygon(idx);
          end
        end
        
        %% Check to duplicate points
        tol = 1e6;
        P = round(points(polygon,:)*tol)/tol;
        [~,ia] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
        polygon = polygon(sort(ia));
        
        %% Test the results using more expensive function
        doPlot = (~verified(key));
        if debug && length(polygon)>3
          DT = delaunayTriangulation(points(polygon,dIdx));
          idx = freeBoundary(DT)';
          idx(2,:) = [];
          k = max(abs(diff(idx)));
          %doPlot = (k>1 && k<(length(idx)-1)) || (~verified(key));
          assert(k==1 || k==(length(idx)-1), 'Two triangle intersection polygon is not convex')
        end
        if debug && doPlot % plot the interesting cases
          obj.PlotTwoTriangles(points, polygon, 'm')
          title(sprintf('key = %i', key));
        end 
        
        end % function
                
        function PlotTwoTriangles(points, polygon, color)
        % Plotting function used for debugging
        nx = size(points,1)-6;
        d = (max(points,[],1)-min(points,[],1))/200;
        figure(2)
        clf
        hold on
        line( points(nx+(1:2),1), points(nx+(1:2),2), points(nx+(1:2),3), 'Color', 'g');
        line( points(nx+(2:3),1), points(nx+(2:3),2), points(nx+(2:3),3), 'Color', 'g');
        line( points(nx+[1,3],1), points(nx+[1,3],2), points(nx+[1,3],3), 'Color', 'g');
        line( points(nx+(4:5),1), points(nx+(4:5),2), points(nx+(4:5),3), 'Color', 'b');
        line( points(nx+(5:6),1), points(nx+(5:6),2), points(nx+(5:6),3), 'Color', 'b');
        line( points(nx+[4,6],1), points(nx+[4,6],2), points(nx+[4,6],3), 'Color', 'b');
        plot3( points(:,1), points(:,2), points(:,3), 'm.');
        if (length(polygon)>2)
          idx = polygon([1:end, 1]);
          plot3( points(idx,1), points(idx,2),points(idx,3), 'Color', color, 'LineWidth', 1);
        end
        for i = 1:nx+6
          text(points(i,1)+d(1), points(i,2)+d(2), points(i,3), num2str(i))
        end
        
        end % function
        
        function [intersect, X] = EdgesIntersect3D(V1,V2, U1,U2)
        %EdgesIntersectPoint3D calculates point of intersection of 2 coplanar
        % segments in 3D
        %
        % INPUTS:
        %   V1,V2 - 1 x 3 coordinates of endpoints of edge 1
        %   U1,U2 - 1 x 3 coordinates of endpoints of edge 2
        % OUTPUT:
        %   X - 1 x 3 coordinates of the intersection point
        A = V2-V1;
        B = U1-U2;
        C = U1-V1;
        %% Solve system of equations [A,B,1] * [d;e;0] = C for d and e
        det3 = @(a,b) ... % determinant of a matrix with columns: [a, b, 1]
          a(:,1).*b(:,2)-a(:,3).*b(:,2) + ...
          a(:,2).*b(:,3)-a(:,2).*b(:,1) + ...
          a(:,3).*b(:,1)-a(:,1).*b(:,3);
        f=det3(A,B); % https://en.wikipedia.org/wiki/Cramer%27s_rule#Explicit_formulas_for_small_systems
        t=det3(C,B)./f; % use Cramer's rule
        s=det3(A,C)./f;
        intersect = (t>=0 & t<=1 & s>=0 & s<=1);
        X = V1 + bsxfun(@times,A,t);
        end % function
        
        function inside = PointInTriangle2D(V1, U)
        % check if V1 is inside triangle U (U1,U2,U3)
        % Algorithm is checking on which side of the half-plane created by the
        % edges the point is. It uses sign of determinant to calculate orientation
        % of point triplets.
        % INPUTS:
        %   V1 - 1 x 2 coordinates of a point
        %   U  - 3 x 2 coordinates of endpoints of 3 edges of a triangle
        % OUTPUT:
        %   inside - a boolean or boolean array
        det2 = @(A,B,C) (A(:,1)-C(:,1))*(B(:,2)-C(:,2)) - (B(:,1)-C(:,1))*(A(:,2)-C(:,2));
        b1 = (det2(U(1,:), U(2,:), V1) > 0);
        b2 = (det2(U(2,:), U(3,:), V1) > 0);
        b3 = (det2(U(3,:), U(1,:), V1) > 0);
        inside = ((b1 == b2) & (b2 == b3)); % inside if same orientation for all 3 edges
        end % function
        
        function [b, c] = otherDim(a)
        % return set [1 2 3] without k
        b = mod(a+1,3)+1;  % b and c are vertices which are on the same side of the plane
        c = 6-a-b;         % a+b+c = 6
        end
       

        function [pVError, nError] = distanceVertex2Mesh(mesh, vertex)
        % DISTANCEVERTEX2MESH - calculate the distance between vertices and a mesh
        %
        % Syntax: [pVError, nError] = distanceVertex2Mesh(mesh, vertice)
        %
        % Inputs:
        %   mesh -  the mesh which is used as a reference for the distance
        %           calculation. 'mesh' needs to be a structure with two fields
        %           called 'vertices' and 'faces', where 'vertices' is a n x 3
        %           matrix defining n vertices in 3D space, and faces is a m x 3
        %           matrix defining m faces with 3 vertice ids each.
        %   vertice(s) -    vertices is a q x 3 matrix defining q vertices in 3D
        %                   space
        %
        % Outputs:
        %   pVError -   q x 1 array containing the shortest distance for each of
        %               the q input vertices to the surface of the 'mesh'
        %   nError -    average normalized error: the error is normalized for a
        %               densely sampled mesh of a unit sphere of radius 1 and
        %               center 0,0,0
        %
        %
        % Example:
        % [x,y,z] = sphere(20); % create a unit sphere of 441 samples
        % ball = surf2patch(x,y,z,'triangles');
        %                       % create triangulation of vertices
        %                       % ball is a structure with faces and vertices
        % vertices = 1.1 .* [x(:),y(:),z(:)];
        %                       % create list of vertices, where all vertices
        %                       % are shifted by 0.1 outwards
        % [pVError, nError] = distanceVertex2Mesh(ball, vertices);
        %                       % pVError is a list of 441 x 1 with 0.1 per entry
        %                       % nError is a single value of 0.1
        %
        % Other m-files required: none
        % Subfunctions:
        %   distance3DP2P (calculate distance of point to vertex)
        %   distance3DP2E (calculate distance of point to edge)
        %   distance3DP2F (calculate distance of point to face)
        % MAT-files required: none
        %
        % See also: surf2patch
        %
        % Author: Christopher Haccius
        % Telecommunications Lab, Saarland University, Germany
        % email: haccius@nt.uni-saarland.de
        % March 2015; Last revision: 26-March-2015
                % Point-to-Point Distance in 3D Space
        function dist = distance3DP2V(v1,v2)
            dist = norm(v1-v2);% Euclidean distance
        end
        
        % Point-to-LineSegment Distance in 3D Space, Line defined by 2 Points (2,3)
        function dist = distance3DP2E(v1,v2,v3)
            d = norm(cross((v3-v2),(v2-v1)))/norm(v3 - v2);
            % check if intersection is on edge
            s = - (v2-v1)*(v3-v2)' / (norm(v3-v2))^2;
            if (s>=0 && s<=1)
                dist = d;
            else
                dist = inf;
            end
        end
        
        % Point-to-Face Distance in 3D Space, Face defined by 3 Points (2,3,4)
        function dist = distance3DP2F(v1,v2,v3,v4)
            n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
            d = abs(n * (v1 - v2)');
            % check if intersection is on face
            n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
            f1 = v1 + d * n;
            f2 = v1 - d * n;
            m = [v3-v2;v4-v2]';
            try 
                r1 = m\(f1-v2)';
            catch
                r1 = [inf;inf];
            end
            try
                r2 = m\(f2-v2)';
            catch
                r2 = [inf;inf];
            end
            if ((sum(r1)<=1 && sum(r1)>=0 && all(r1 >=0) && all(r1 <=1)) || ...
                    (sum(r2)<=1 && sum(r2)>=0 && all(r2 >=0) && all(r2 <=1)))
                dist = d;
            else
                dist = inf;
            end
        end
        warning ('off','MATLAB:rankDeficientMatrix'); % turn off warnings for 
                    % defficient ranks occuring in point-to-face distance
                    % calculation
        
        %vertices = mesh.vertices;
        %faces = mesh.faces;
        %morita Modified
        vertices = mesh.Points;
        faces = mesh.ConnectivityList;
        
        [numV,dim] = size(vertices);
        [numF,pts] = size(faces);
        
        [tV,dimT] = size(vertex);
        
        if(pts~=3)
            error('Only Triangulations allowed (Faces do not have 3 Vertices)!');
        elseif (dim~=dimT || dim~=3)
            error('Mesh and Vertices must be in 3D space!');
        end
        
        % initialie minimal distance to infinty
        d_min = Inf(tV,1);
        
        % first check: find closest vertex
        for c1 = 1:tV % iterate over all test vertices
            for c2 = 1:numV % iterate over all mesh vertices
                v1 = vertex(c1,:);
                v2 = vertices(c2,:);
                d = distance3DP2V(v2,v1); % Euclidean distance
                if d < d_min(c1)
                    d_min(c1) = d;
                end
            end
        end
        
        % second check: find closest edge
        for c1 = 1:tV % iterate over all test vertices
            for c2 = 1:numF % iterate over all faces
                for c3 = 1:2 % iterate over all edges of the face
                    for c4 = c3+1:3
                        v1 = vertex(c1,:);
                        v2 = vertices(faces(c2,c3),:);
                        v3 = vertices(faces(c2,c4),:);
                        % check if edge is possible
                        if ( all(min(v2,v3) < (v1 + d_min(c1))) && ...
                             all(max(v2,v3) > (v1 - d_min(c1))) )
                            d = distance3DP2E(v1,v2,v3);
                            if d < d_min(c1) % d is shorter than previous shortest
                                d_min(c1) = d;
                            end
                        end
                    end
                end
            end
        end
        
        % third check: find closest face
        for c1 = 1:tV % iterate over all test vertices
            for c2 = 1:numF % iterate over all faces
                v1 = vertex(c1,:);
                v2 = vertices(faces(c2,1),:);
                v3 = vertices(faces(c2,2),:);
                v4 = vertices(faces(c2,3),:);
                % check if face is possible
                if ( all(min([v2;v3;v4]) < (v1 + d_min(c1))) && ...
                     all(max([v2;v3;v4]) > (v1 - d_min(c1))) )
                    d = distance3DP2F(v1,v2,v3,v4);
                    if d < d_min(c1)
                        d_min(c1) = d;
                    end
        
                end
            end
        end
        
        pVError = d_min; % output error per vertex is d_min
        minS = min(vertices); % get size of mesh
        maxS = max(vertices);
        s = 1/sqrt(3) * norm(0.5 * (maxS - minS)); % calculate size of mesh
        nError = sum(pVError) / (s * tV); % average and normalize error
        end % end of function
        
       
    
    end
end