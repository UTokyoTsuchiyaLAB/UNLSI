%
%%%%%%%%ここから
clear;
checkVariables = 1;
orgVal = [4,4,4,4,4,3];
%orgVal = modifyDesFile("org.des","org.des");


ub = [  4,   4,   4,   4,   4, 6];
lb = [  1,   1,   1,   1,   1, 0];
cmin = [0.525,-4,-4,-4,-4]';
cmax = [0.575, 0, 0, 0, 0]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal);
ungradetest = UNGRADE(@(x)vspSurfGen(x,"wing","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
if checkVariables == 1
    ungradetest.checkSurfGenWork(0.5,1);
end
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"wing","org.des"));%基準メッシュ生成関数の指定
ungradetest = ungradetest.setOptFlowCondition(0.0001,5,0,20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptimization(eye(numel(orgVal)).*1,0.3,0.001,0.6,"BFGS",6);

%ungradetest = ungradetest.updateVariables(@objFun,"nonlin",cmin,cmax);%設計変数微分の計算方法の指定
for i = 1:300
    [dx,ungradetest] = ungradetest.updateVariables(@objFun,"direct",cmin,cmax);%設計変数微分の計算方法の指定
end
%ungradetest = ungradetest.solve

%}