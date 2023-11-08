    clear;
    orgVal = [2,2,2,2,2]; %設計変数の初期値
    
    Machrange = [0.001]; %解析するマッハ数
    alpharange = [5]; %解析する迎角
    betarange = [0]; %解析する横滑り角
    
    ub = [  4,   4,   4,   4,   4]; %設計変数の下限値
    lb = [  1,   1,   1,   1,   1]; %設計変数の上限値
    cmin = []';%制約条件の下限値
    cmax = []';%制約条件の上限値
    
    ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1);%コンストラクタの実行
    ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
    %%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
    ungradetest = ungradetest.setUNLSISettings("propCalcFlag",1,"Vinf",50);
    ungradetest = ungradetest.setProp(1,2,2,0);
    [ungradetest,thrust,power] = ungradetest.setPropState(1,0.04,0.02,5000);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ungradetest = ungradetest.setFlowCondition(alpharange,betarange,Machrange,500000);
    ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1)); %現在の機体形状で解析を実行する。（この中でパネル法行列が作成される。この関数を実行しなくても必要とあらば作成される）
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする
    [Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する

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