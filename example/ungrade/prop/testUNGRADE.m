
%
clear;
orgVal = [ones(1,6)*30, ones(1,5)*0.2]; %設計変数の初期値

Machrange = [0.001]; %解析するマッハ数
alpharange = [0]; %解析する迎角
betarange = [0]; %解析する横滑り角
geomErr = 0.005;

lb = [ones(1,6)* 0, ones(1,5)*0.05]; %設計変数の下限値
ub = [ones(1,6)*90, ones(1,5)*0.3]; %設計変数の上限値
cmin = [25/100]';%制約条件の下限値 //推力
cmax = [30/100]';%制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"prop","org.des"),@(x)vspGeomGen(x,"prop","org.des"),orgVal,lb,ub,0);%コンストラクタの実行
ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNLSISettings("Vinf",7,"nWakeMax",30);%主流速度の設定
ungradetest = ungradetest.setUNGRADESettings("gradientCalcMethod","chain");%設計変数の数が多いのでchainを選択
ungradetest = ungradetest.setRotation(1:2,[0,0,0],[-135*pi/30*180/pi,0,0]/ungradetest.settingUNLSI.Vinf);%回転角速度の設定
ungradetest = ungradetest.setHelixWake(1:2,0.02,135,[1,0,0],[0,0,0]);%螺旋wakeの設定
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setFlowCondition(alpharange,betarange,Machrange,500000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun); %評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-20,1]); %機体形状と圧力係数をプロットする
M = getframe(gcf);
ungradetest.plotWakeShape(1); %↑の結果をプロットする
%}

%
for i = 1:200
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    ungradetest = ungradetest.setHelixWake(1:2,0.02,135,[1,0,0],[0,0,0]);
    ungradetest = ungradetest.setRotation(1:2,[0,0,0],[-135*pi/30*180/pi,0,0]/ungradetest.settingUNLSI.Vinf);
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-20,1]); %↑の結果をプロットする
        M = [M,getframe(gcf)];
        videoMaker(M,"propGeom",10);
        ungradetest.plotWakeShape(1); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2,"propConv"); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
save prop.mat
%}


