%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata); %コンストラクタの実行
wing = wing.setREFS(80,20,4,[0,0,0]); %基準面積 基準長の設定
wing = wing.setUNLSISettings("propCalcFlag",1,"Vinf",50,"propWakeLength",0.3);
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.setProp(1,2,5,0);
wing = wing.makeEquation(); %パネル法行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}

[wing,thrust,power,Jref] = wing.setPropState(1,0.4,0.6,2000);
wing = wing.solveFlow(5,0,0.001,500000);%パネル法を解く
disp(wing.AERODATA{1}); %結果の表示
wing.plotGeometry(1,wing.Cp{1}(:,1),[-2,1]);%圧力係数のプロット