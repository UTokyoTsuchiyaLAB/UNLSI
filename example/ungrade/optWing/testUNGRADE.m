clear;
orgVal = [0.47,31.9,0.2,-10,6,0,0.39,0.15]; %設計変数の初期値

lb = [0.4,25,0.1,-20,0,-0.04,0.2,0.1]; %設計変数の下限値
ub = [0.6,45,0.3,  0,10, 0.04,0.4,0.2]; %設計変数の上限値
cmin = [1.0, -100, 1.2/2]'; %制約条件の下限値
cmax = [1.1, -0.2, 1.3/2]'; %制約条件の上限値
geomErr = 0.001;

ungradetest = UNGRADE(@(x)vspMeshGen(x,"optWing2","org.des"),@(x)vspGeomGen(x,"optWing2","org.des"),orgVal,lb,ub,1,0.001,0,0); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
ungradetest = ungradetest.setProp(1,3,0.15,1);
ungradetest = ungradetest.setDeflAngle(4,[0,1,0],0);
ungradetest = ungradetest.setDeflGroup(1,"elev",4,1);
[ungradetest,thrust,power,Jref] = ungradetest.setPropState(15,1.225,0.4,0.6,8000);
ungradetest = ungradetest.setCfParameter(100000,4,0.052*(10^-5),0,1); %摩擦係数関連のパラメータのセット
ungradetest = ungradetest.setOptions('n_divide',1,'updateMethod','Levenberg–Marquardt');
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


for i = 1:200
    ungradetest = ungradetest.calcLagrangianGradient(@objFun,cmin,cmax,"TrustRegion",0.05);%次の設計変数を計算する
    ungradetest = ungradetest.updateHessian();%準ニュートン法のヘッシアン更新
    [nextVar,dxNorm] = ungradetest.descentLagrangian();%次の設計変数を計算する
    if mod(ungradetest.iteration,20)==0
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
    else
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%次の設計変数を用いて機体形状を更新する
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
        
        if mod(ungradetest.iteration,20)==0
            ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar);%次の設計変数を用いて機体形状を更新する
        else
            ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%次の設計変数を用いて機体形状を更新する
        end
        ungradetest = ungradetestdx;
        [Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);
    end
    
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end

