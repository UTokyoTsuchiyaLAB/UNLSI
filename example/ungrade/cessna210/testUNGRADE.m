clear;
orgVal = [0.04,0.4,0.12,0,0.04,0.4,0.12,5.5,0,4,2]; %設計変数の初期値

lb = [-0.05,0.2,0.1,0,-0.05,0.2,0.1,4.0,0,2,0]; %設計変数の下限値
ub = [ 0.05,0.5,0.14,6,0.05,0.5,0.14,8,20,8,6]; %設計変数の上限値
cmin = [40/100, 0.2 -0.3]'; %制約条件の下限値
cmax = [45/100, 1.0  0.3]'; %制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"Cessna-210","org.des"),@(x)vspGeomGen(x,"Cessna-210","org.des"),orgVal,lb,ub,1,[0.001,0.001,0.001],[0,5,10],[0,0,0]); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1); %摩擦係数関連のパラメータのセット
ungradetest = ungradetest.solveAnalysis(1,0,0); %現在の機体形状で解析を実行する。（この中でパネル法行列が作成される。この関数を実行しなくても必要とあらば作成される）
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5); %次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar); %次の設計変数を用いて機体形状を更新する
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.5);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.0);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
%}
