%
%%%%%%%%ここから
%
clear;
dt = 0.02;
theta = 0;
rpm = 135;
Vinf = 7;

[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata,0,0.001); %コンストラクタの実行
wing = wing.setREFS(1,3,0.0001,[0,0,0]); %基準面積 基準長の設定
wing = wing.setUNLSISettings("nCalcDivide",2);
wing = wing.setUNLSISettings("CpLimit",[-1000,1000],"Vinf",Vinf,"nWakeMax",30);
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%{
[wing,CT,Cp,Cq,efficiency] = wing.solveSteadyProp(1:2,rpm,[1,0,0],[0,0,0],0,0,0.001,500000,2);
disp([CT*(rpm/60)^2*wing.BREF^4*1.225,Cp*(rpm/60)^3*wing.BREF^5*1.225,Cq*(rpm/60)^2*wing.BREF^5*1.225,efficiency]);
%}

%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 0;
wing = wing.setWakeShape([0.1,0,0]);
for i = 1:5000
    [wing,CT,Cp,Cq,efficiency] = wing.solveUnsteadyProp(1:2,dt,rpm,[1,0,0],[0,0,0],alpha,0,0.001,100000,3,[-20,1]);
    disp([CT*(rpm/60)^2*wing.BREF^4*1.225,Cp*(rpm/60)^3*wing.BREF^5*1.225,Cq*(rpm/60)^2*wing.BREF^5*1.225,efficiency]);
    wing = wing.makeEquation(); %パネル法行列の作成
end
%}