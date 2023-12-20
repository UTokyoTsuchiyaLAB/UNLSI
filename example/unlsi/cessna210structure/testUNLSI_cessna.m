%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "Cessna-210.vspgeom", 0);%形状の読み込み
[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh('Cessna-210_NormalWing_Struct0.msh');
wing = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
wing = wing.setREFS(175,36.75,4.91); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%}
%
wing = wing.setFemMesh(vtx,srf,fc);
wing = wing.setFemMaterials([1,2,3],[0.002,0.1,0.01],[5000000000,5000000000,5000000000]);
wing = wing.connectAeroIDandFemID(5,1);
wing = wing.makeFemEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}
alpha = 0;
wing = wing.solveFlow(alpha,0,0.001,2.3*10^6); %パネル法を解く
disp(wing.getAERODATA(alpha));
%wing.plotGeometry(1,wing.getCp(0,0,0.001,2.3*10^6),[-2,1],"exact"); %圧力係数のプロット

wing2 = wing;%変位後の空力解析用のインスタンスを準備
for iter = 1:3
    delta = wing.solveFem(wing2.getCp(alpha).*50.^2.*1.225.*0.5);
    modVerts = wing.calcModifiedVerts(delta{1});
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,2.3*10^6);%パネル法を解く
    disp(wing2.getAERODATA(alpha));
    wing2.plotGeometry(1,wing2.getCp(alpha),[-2,1]);%圧力係数のプロット
end