clear;
orgVal = [0.04,0.12,0,0.04,0.12,5.5,0,4,2]; %設計変数の初期値

lb = [-0.05,0.1,0,-0.05,0.10,4.0,-5,2,0]; %設計変数の下限値
ub = [ 0.05,0.14,6,0.05,0.14,8.0, 5,8,6]; %設計変数の上限値
cmin = [4.5]'; %制約条件の下限値
cmax = [5.0]'; %制約条件の上限値
geomErr = 0.15;

ungradetest = UNGRADE(@(x)vspMeshGen(x,"Cessna-210","org.des"),@(x)vspGeomGen(x,"Cessna-210","org.des"),orgVal,lb,ub,1); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNGRADESettings('updateMethod','Levenberg–Marquardt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ungradetest = ungradetest.setFlowCondition(0,0,0.262,1000000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

for i = 1:60
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
save cessna.mat
