clear;
orgVal = [2,2,2,2,2,2]; %設計変数の初期値

Machrange = [0.001]; %解析するマッハ数
alpharange = [2]; %解析する迎角
betarange = [0]; %解析する横滑り角
geomErr = 0.005;

lb = [  1,   1,   1,   1,   1,   1]; %設計変数の下限値
ub = [  4,   4,   4,   4,   4,   4]; %設計変数の上限値
cmin = [0.15]';%制約条件の下限値
cmax = [0.17]';%制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1);%コンストラクタの実行
ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setFlowCondition(alpharange,betarange,Machrange,500000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun); %評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


%%%%%%%%%%%%%%%%%初級者向け -- とりあえず実行%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:5
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,1);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        ungradetest = ungradetest.solveAnalysis(1,alpharange(1),0); %現在の形状に対し任意の迎角・横滑り角・各速度で解析を行う
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end

%%%%%%%%%%%%%%%%%上級者向け -- calcNextVariablesの中身を個別実行%%%%%%%%%%%%%%%%%%
for i = 1:30
    ungradetest = ungradetest.calcLagrangianGradient(@objFun,cmin,cmax,"TrustRegion",0.2);%次の設計変数を計算する
    ungradetest = ungradetest.updateHessian();%準ニュートン法のヘッシアン更新
    [nextVar,dxNorm] = ungradetest.descentLagrangian();%次の設計変数を計算する
    ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%MeshDeformationにより次の設計変数を用いて機体形状を更新する。
    if max(distanceVertex2Mesh(ungradetest.orgGeom, ungradetest.orgMesh.Points))>geomErr
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%Mesh計算を行い、次の機体形状を更新する
    end
    [Idx,condx,ungradetestdx] = ungradetestdx.evaluateObjFun(@objFun);
    penartyorg = 100 * sum(max(max(0,cmin-conorg),max(conorg-cmax)));
    penartydx = 100 * sum(max(max(0,cmin-condx),max(condx-cmax)));
    if Idx + penartydx<Iorg + penartyorg
        ungradetest = ungradetestdx;
        Iorg = Idx;
        conorg = condx;
    else
        nextVar = ungradetest.descentLagrangian("TrustRegion",dxNorm*0.5);%次の設計変数を計算する。ここでのsettingの変更は残らない
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%%MeshDeformationにより次の設計変数を用いて機体形状を更新する。
        if max(distanceVertex2Mesh(ungradetest.orgGeom, ungradetest.orgMesh.Points))>geomErr
            ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%%Mesh計算を行い、次の機体形状を更新する
        end
        ungradetest = ungradetestdx;
        [Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);
    end
    
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end

