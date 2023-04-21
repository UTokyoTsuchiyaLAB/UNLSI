clear;
orgVal = [0.04,0.12,0,0.04,0.12,5.5,0,4,2]; %設計変数の初期値

lb = [-0.05,0.1,0,-0.05,0.1,4.0,0,2,0]; %設計変数の下限値
ub = [ 0.05,0.14,6,0.05,0.14,8,20,8,6]; %設計変数の上限値
cmin = [40/100]'; %制約条件の下限値
cmax = [45/100]'; %制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"Cessna-210","org.des"),@(x)vspGeomGen(x,"Cessna-210","org.des"),orgVal,lb,ub,1,[0.001],[0],[0]); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1); %摩擦係数関連のパラメータのセット
[Iorg,~,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


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

