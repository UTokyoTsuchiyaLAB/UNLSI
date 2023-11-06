%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata); %コンストラクタの実行
wing = wing.setREFS(80,20,4); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}


wing = wing.solveFlow([0,5,10],[0,0,0],0.001,500000);%パネル法を解く
AERODATA = wing.getAERODATA(5,0);
wing.plotGeometry(1,wing.getCp(0,0,0.001,400000),[-2,1]);%圧力係数のプロット
