%
%%%%%%%%ここから
%
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata); %コンストラクタの実行
wing = wing.setREFS(80,20,4,[0,0,0]); %基準面積 基準長の設定
wing = wing.setUNLSISettings("nCalcDivide",2);%パネル法行列の作成における分割数を設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.setWakeShape([1,0,0]);
wing = wing.makeEquation(); %パネル法行列の作成
%}
%}
%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 10;
Vinf =250;
Re = Vinf * 4 / 1.512 * 1e5;
wing = wing.setUNLSISettings("Vinf",Vinf);
wing = wing.solveFlow(alpha,0,0.001,Re);%パネル法を解く
wing.plotGeometry(1,wing.getCp(alpha,0,0.001,Re),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha,0));
[con,verts,femID] = readFemMesh('wing_WingGeom_Struct0.msh');
wing = wing.setFemMesh(verts,con,femID);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
[wing,weight] = wing.setFemMaterials([1,2,3],[0.01,0.001,0.001],[73500000000,73500000000,73500000000],[2700,2700,2700],[1000,1000,1000]);
disp("weight");
disp(weight);
wing = wing.makeFemEquation();

wing = wing.femModalAnalysis(4);%モード解析用の固有値解析
%wing = wing.calcApproximatedEquation();%空力近似マトリックスの作成
%モード形状
for i = 1:4
    delta{1} = wing.femSol2Delta(wing.femEigenVec(:,i)./max(abs(wing.femEigenVec(:,i))));
    modVerts = wing.calcModifiedVerts(delta{1}.*10); %構造メッシュの変形に従って空力メッシュを変形
    wing2 = wing.setVerts(modVerts); %節点の移動のみ
    wing2.plotGeometry(i+2,wing.getCp(alpha,0,0.001,Re),[-3,1.5]);
end
%
%%%%%%%以下空力弾性計算
dt = 0.01;
wing2 = wing;
for i = 1:500
    disp(i)
    tic;
    if i == 1
        %[delta,deltadot]  = wing.solveAeroelastic([0,dt],[],[],wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%直接空力弾性計算。deltaとdeltadotに空行列を渡すと、初期値0,0からスタート
        [z,zdot,delta,deltadot]  = wing.solveModalAeroelastic([0,dt],[],[],wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%モード分解による空力弾性計算。zとzdotに空行列を渡すと、初期値0,0からスタート。低次3モードを積分
    else
        %[delta,deltadot]  = wing.solveAeroelastic([0,dt],delta,deltadot,wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%初期値deltaとdeltadotから、tspan間での空力弾性応答を計算
        [z,zdot,delta,deltadot]  = wing.solveModalAeroelastic([0,dt],z,zdot,wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%モード分解による空力弾性計算。
    end
    toc;
    modVerts = wing.calcModifiedVerts(delta{1});
    
    %wing2 = wing.makeApproximatedInstance(modVerts);%空力近似マトリックスを用いたインスタンスの作成
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.marchWake(dt,alpha,0,0.001,Re);
    wing.wakeline = wing2.wakeline;
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(2,wing2.getCp(alpha,0,0.001,Re),[-3,1.5]);
    %wing2.plotWakeShape(2);
    M(i) = getframe(gcf);
    videoMaker(M,"wingCalculating",1/dt)
end
movie(M,1,1/dt);
%}