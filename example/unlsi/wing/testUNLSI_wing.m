%
%%%%%%%%ここから
%
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata); %コンストラクタの実行
wing = wing.setREFS(80,20,4); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.setUNLSISettings("nCalcDivide",2);
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%}
%}
%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 10;
wing = wing.solveFlow(alpha,0,0.001,500000);%パネル法を解く
disp(wing.getAERODATA(alpha,0));
[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh('wing_WingGeom_Struct0.msh');
wing = wing.setFemMesh(vtx,srf,fc);
wing = wing.setFemMaterials([1,2,3],[0.002,0.01,0.05],[5000000000,5000000000,5000000000]);
wing = wing.connectAeroIDandFemID([1,2],1);
wing = wing.makeFemEquation();
%}
wing2 = wing;
for iter = 1:4
    delta = wing.solveFem(wing2.getCp(alpha,0,0.001,500000).*30.^2.*1.225.*0.5);
    modVerts = wing.calcModifiedVerts(delta{1});
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,500000);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(1,wing2.getCp(alpha,0,0.001,500000),[-2,1]);%圧力係数のプロット
end
%}