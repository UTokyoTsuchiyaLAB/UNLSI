clear;
orgVal = modifyDesFile("initial.des");

lb = [2, 0,2,0,-ones(1,4).*0.25 ,-ones(1,4).*0.0,0,-ones(1,4).*0.25 ,-ones(1,4).*0.0]; %設計変数の下限値
ub = [8,10,8,5, ones(1,4).*0.0 , ones(1,4).*0.25,5, ones(1,4).*0.0 , ones(1,4).*0.25]; %設計変数の上限値
cmin = [30]'; %制約条件の下限値
cmax = [35]'; %制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"Cessna-210","org.des"),@(x)vspGeomGen(x,"Cessna-210","org.des"),orgVal,lb,ub,1,0.1); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNGRADESettings("femCouplingFlag",0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ungradetest = ungradetest.setFlowCondition(0,0,0.262,1000000);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

M = getframe(gcf);
for i = 1:400
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax);%次の設計変数を計算する
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    %%%%%%%%%ここから
        [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
        ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
        M = [M,getframe(gcf)];
        videoMaker(M,"cessnaGeom",10);
        ungradetest.plotOptimizationState(2,"cessnaConv"); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
save cessna.mat