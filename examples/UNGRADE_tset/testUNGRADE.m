%
%%%%%%%%ここから
clear;
checkVariables = 1;
orgVal = [4,4,4,4,4,0];
%orgVal = modifyDesFile("org.des","org.des");


ub = [  4,   4,   4,   4,   4, 3];
lb = [  1,   1,   1,   1,   1, -3];
cmin = [0.325,-4,-4,-4,-4]';
cmax = [0.375, 0, 0, 0, 0]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal);
ungradetest = UNGRADE(@(x)vspSurfGen(x,"wing","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
ungradetest.checkSurfGenWork(0.5);
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"wing","org.des"));%基準メッシュ生成関数の指定
ungradetest = ungradetest.setOptFlowCondition(0.0001,[5,10],[0,0],20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptimization(eye(numel(orgVal)).*1,0.3,0.001,0.6,"BFGS",6);

%ungradetest = ungradetest.updateVariables(@objFun,"nonlin",cmin,cmax);%設計変数微分の計算方法の指定
for i = 1:300
    [dx,ungradetest] = ungradetest.updateVariables(@objFun,"direct",cmin,cmax);%設計変数微分の計算方法の指定
end
%ungradetest = ungradetest.solve

%}