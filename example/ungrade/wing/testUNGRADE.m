clear;
orgVal = [4,4,4,4,4]; %設計変数の初期値

Machrange = [0.001,0.001,0.001]; %解析するマッハ数
alpharange = [2,4,6]; %解析する迎角
betarange = [0,0,0]; %解析する横滑り角

ub = [  4,   4,   4,   4,   4]; %設計変数の下限値
lb = [  1,   1,   1,   1,   1]; %設計変数の上限値
cmin = [0,0,0,0]';%制約条件の下限値
cmax = [10,10,10,10]';%制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1,Machrange,alpharange,betarange);%コンストラクタの実行
ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1); %摩擦係数関連のパラメータのセット
ungradetest = ungradetest.setOptions('n_divide',3); %設定値の変更 n_divide:パネル法行列を計算するときの分割数
ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1)); %現在の機体形状で解析を実行する。（この中でパネル法行列が作成される。この関数を実行しなくても必要とあらば作成される）
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
    ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.5);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    %ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));%表示等を行わない場合
    %ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    %ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.0);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    %ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
    %ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    %ungradetest.plotOptimizationState(2);
end

ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
ungradetest.plotOptimizationState(2);
%}
