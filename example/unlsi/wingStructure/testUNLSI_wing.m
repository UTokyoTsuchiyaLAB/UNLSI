%
%%%%%%%%ここから
%
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
wing = wing.setREFS(80,20,4,[0,0,0]); %基準面積 基準長の設定
wing = wing.setUNLSISettings("nCalcDivide",2);%パネル法行列の作成における分割数を設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%}
%}
%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 10;
Vinf = 100;
Re = Vinf * 4 / 1.512 * 1e5;
wing = wing.solveFlow(alpha,0,0.001,Re);%パネル法を解く
wing.plotGeometry(1,wing.getCp(alpha,0,0.001,Re),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha,0));
[con,verts,femID] = readFemMesh('wing_WingGeom_Struct0.msh');
wing = wing.setFemMesh(verts,con,femID);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
[wing,weight] = wing.setFemMaterials([1,2,3],[0.002,0.003,0.01],[73500000000,73500000000,73500000000],[2700,2700,2700]);
disp("weight");
disp(weight);
wing = wing.makeFemEquation();
wing = wing.femModalAnalysis(5);
%

wing2 = wing;
for iter = 1
    delta = wing.solveModalFem(wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);
    %delta = wing.solveFem(wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);
    deltaInterp = wing.interpFemDelta(delta{1},[0,10,0;4,10,0]);
    disp("Tip twist(deg)");
    disp(atan2d(deltaInterp(1,3)-deltaInterp(2,3),4));
    modVerts = wing.calcModifiedVerts(delta{1});
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(2,wing2.getCp(alpha,0,0.001,Re),[-3,1.5]);%圧力係数のプロット
end
