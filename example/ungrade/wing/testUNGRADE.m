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
[Iorg,~,ungradetest] = ungradetest.evaluateObjFun(@objFun); %評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


%%%%%%%%%%%%%%%%%初級者向け -- とりあえず実行%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:15
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        ungradetest = ungradetest.solveAnalysis(1,alpharange(1),0); %現在の形状に対し任意の迎角・横滑り角・各速度で解析を行う
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end

%%%%%%%%%%%%%%%%%上級者向け -- calcNextVariablesの中身を個別実行%%%%%%%%%%%%%%%%%%
for i = 1:30
    ungradetest = ungradetest.calcLagrangianGradient(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%次の設計変数を計算する
    ungradetest = ungradetest.updateHessian();%準ニュートン法のヘッシアン更新
    nextVar = ungradetest.descentLagrangian();%次の設計変数を計算する
    ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
    [Idx,~,ungradetestdx] = ungradetestdx.evaluateObjFun(@objFun);
    if Idx<Iorg
        ungradetest = ungradetestdx;
        Iorg = Idx;
    else
        nextVar = ungradetest.descentLagrangian("betaLM",2.0);%次の設計変数を計算する。ここでのsettingの変更は残らない
        ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
        [Iorg,~,ungradetest] = ungradetest.evaluateObjFun(@objFun);
    end
    
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end

