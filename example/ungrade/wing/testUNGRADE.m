clear;
orgVal = [3,3,3,3,3,3]; %設計変数の初期値

Machrange = [0.001]; %解析するマッハ数
alpharange = [5]; %解析する迎角
betarange = [0]; %解析する横滑り角
geomErr = 0.005;

lb = [  1,     1,   1,1,1,1]; %設計変数の下限値
ub = [  4,     4,   4,4 4,4]; %設計変数の上限値
cmin = [2.2, 0, 0, 0, 0, 0]';%制約条件の下限値
cmax = [2.4 10,10,10,10,10]';%制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1);%コンストラクタの実行
ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setFlowCondition(alpharange,betarange,Machrange,500000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun); %評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


%%%%%%%%%%%%%%%%%初級者向け -- とりあえず実行%%%%%%%%%%%%%%%%%%%%%%%%%
%
for i = 1:60
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.1);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        ungradetest.plotOptimizationState(2); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
%}

%%%%%%%%%%%%%%%%%上級者向け -- calcNextVariablesの中身を個別実行%%%%%%%%%%%%%%%%%%
%{
for iter = 1:60
    TR = 0.05;
    ungradetest = ungradetest.calcLagrangianGradient(@objFun,cmin,cmax,"TrustRegion",TR,"betaLM",0.1);%次の設計変数を計算する
    ungradetest = ungradetest.updateHessian();%準ニュートン法のヘッシアン更新
    [nextVar,dxNorm] = ungradetest.descentLagrangian();%次の設計変数を計算する
    for i = 1:30  
        ungradetestdx= ungradetest.updateMeshGeomfromVariables(nextVar,1);%Mesh作成により次の設計変数を用いて機体形状を更新する。
        [Ldx,Idx,condx,ungradetestdx] = ungradetestdx.calcLagrangian(@objFun,cmin,cmax);
        %penartyorg = 100 * sum(max(max(0,cmin-conorg),max(conorg-cmax)));
        %penartydx = 100 * sum(max(max(0,cmin-condx),max(condx-cmax)));
        disp("Ldx,Lorg")
        disp([Ldx,ungradetest.LagrangianInfo.Lorg])
        if Ldx <= ungradetest.LagrangianInfo.Lorg || TR < 1e-3
            ungradetest = ungradetestdx;
            Iorg = Idx;
            conorg = condx;
            break;
        else
            for j = 1:30
                TR = TR * 0.9; %信頼領域法の精度検証に相当する。
                nextVar2 = ungradetest.descentLagrangian("TrustRegion",TR);%次の設計変数を計算する。ここでのsettingの変更は残らない
                if any(abs(nextVar2-nextVar)>0.0011) 
                    nextVar = nextVar2;
                    break
                else
                    nextVar = nextVar2;
                end
            end
        end
    end
    %ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,1);%Mesh作成により次の設計変数を用いて機体形状を更新する。
    %[Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax);
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
%}

