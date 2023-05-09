%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
wing = wing.checkMesh(sqrt(eps)); %メッシュ状態のチェック
wing = wing.makeCluster(50,50); %速度分布を求めるためのパネルクラスターを作成
wing = wing.setProp(1,2,5,0);
wing = wing.makePropEquation(100);
wing = wing.makeEquation(20,5,10); %パネル法行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}

wing = wing.flowCondition(1,0.0001); %flowNoのマッハ数を指定 
wing = wing.setREFS(80,20,4); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.setCf(1,500000,4,0.052*(10^-5),0); %摩擦係数パラメータの指定
wing = wing.setPropState(100,1.225,0.4,0.6,2000);
wing = wing.solveFlow(1,10,0,[],1);%パネル法を解く
disp(wing.AERODATA{1}); %結果の表示
wing.plotGeometry(1,wing.Cp{1}(:,1),[-2,1]);%圧力係数のプロット