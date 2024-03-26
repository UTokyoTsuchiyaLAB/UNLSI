%
%%%%%%%%ここから
%
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata,1,0.001); %コンストラクタの実行
wing = wing.setREFS(80,20,4); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.setUNLSISettings("nCalcDivide",2);
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%}
%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 3;
%wing = wing.solveFlow(alpha,0,0.001,500000);%パネル法を解く
wing = wing.flowCondition(1,0.001,500000);
u0 = wing.solvePertPotential(1,alpha,0);
[~,~,~,~,wing] = wing.solveFlowForAdjoint(u0,1,alpha,0);
disp(wing.getAERODATA(alpha,0));
wing.plotGeometry(1,wing.getCp(alpha,0,0.001,500000),[-2,1]);%圧力係数のプロット