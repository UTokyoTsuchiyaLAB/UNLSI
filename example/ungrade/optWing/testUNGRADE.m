clear;
orgVal = modifyDesFile("init.des");%設計変数の初期値

lb = [0.30,0.00 ,0.02,0.00 ,0.02,-0.4,-0.4, 0.0,-0.1,-0.1,-5,0.15,-10,-10,-0.4,-0.4,0.0,-0.1,-0.1,0.30]; %設計変数の下限値
ub = [0.60,0.01 ,0.10,0.01 ,0.10, 0.1, 0.1, 0.4, 0.4, 0.4,25,0.45,  0, 10, 0.1, 0.1,0.4, 0.4, 0.4,0.45]; %設計変数の上限値
cmin = [0.1,  -0.80]'; %制約条件の下限値
cmax = [0.11, -0.40]'; %制約条件の上限値
geomErr = 0.001;

ungradetest = UNGRADE(@(x)vspMeshGen(x,"optWing2","org.des"),@(x)vspGeomGen(x,"optWing2","org.des"),orgVal,lb,ub,1,0.001); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNGRADESettings('TrustRegion',0.15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ungradetest = ungradetest.setFlowCondition([0,3],[0,0],[0.001,0.001],200000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

M = getframe(gcf);
for i = 1:240
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax);%次の設計変数を計算する
        ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        M = [M,getframe(gcf)];
        videoMaker(M,"geometry",10);
        ungradetest.plotOptimizationState(2,"optWing2"); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
save optWing2.mat
