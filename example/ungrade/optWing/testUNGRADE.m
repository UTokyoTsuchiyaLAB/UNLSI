clear;
orgVal = [0.47,31.9,0.2,-10,6,0,0.39,0.15]; %設計変数の初期値

lb = [0.4,25,0.1,-20,0,-0.04,0.2,0.1]; %設計変数の下限値
ub = [0.6,45,0.3,  0,10, 0.04,0.4,0.2]; %設計変数の上限値
cmin = [1.0, -100, 1.2/2]'; %制約条件の下限値
cmax = [1.1, -0.4, 1.3/2]'; %制約条件の上限値
geomErr = 0.001;

ungradetest = UNGRADE(@(x)vspMeshGen(x,"optWing2","org.des"),@(x)vspGeomGen(x,"optWing2","org.des"),orgVal,lb,ub,1); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNLSISettings("propCalcFlag",1,"Vinf",15);
ungradetest = ungradetest.setUNGRADESettings("dynCoefFlag",1);
ungradetest = ungradetest.setProp(1,3,0.15,1);
[ungradetest,thrust,power] = ungradetest.setPropState(1,0.4,0.6,8000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ungradetest = ungradetest.setFlowCondition(0,0,0.001,500000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


for i = 1:60
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.1);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
