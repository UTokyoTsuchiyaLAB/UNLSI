clear;
orgVal = [0.04,0.12,0,0.04,0.12,5.5,0,4,2]; %設計変数の初期値

lb = [-0.05,0.1,0,-0.05,0.1,4.0,0,2,0]; %設計変数の下限値
ub = [ 0.05,0.14,6,0.05,0.14,8,20,8,6]; %設計変数の上限値
cmin = [40]'; %制約条件の下限値
cmax = [45]'; %制約条件の上限値
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
    ungradetest = ungradetest.calcLagrangianGradient(@objFun,cmin,cmax,"TrustRegion",0.05);%次の設計変数を計算する
    ungradetest = ungradetest.updateHessian();%準ニュートン法のヘッシアン更新
    [nextVar,dxNorm] = ungradetest.descentLagrangian();%次の設計変数を計算する
    ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%次の設計変数を用いて機体形状を更新する
    if max(distanceVertex2Mesh(ungradetest.orgGeom, ungradetest.orgMesh.Points))>geomErr
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
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
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%次の設計変数を用いて機体形状を更新する
        if max(distanceVertex2Mesh(ungradetest.orgGeom, ungradetest.orgMesh.Points))>geomErr
            ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
        end
        ungradetest = ungradetestdx;
        [Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);
    end
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
