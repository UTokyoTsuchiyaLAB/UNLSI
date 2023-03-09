%
%%%%%%%%ここから
clear;
orgVal = [4.0000,4.0000, 4.0000,4.0000,4.0000];
%orgVal = modifyDesFile("org.des","org.des");

ub = [  4,   4,   4,   4,   4];
lb = [0.5, 0.5, 0.5, 0.5, 0.5];
cmin = [0.45, -4,-4,-4,-4]';
cmax = [0.55, 0, 0, 0, 0]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal);
ungradetest = UNGRADE(@(x)vspSurfGen(x,"wing","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"wing","org.des"));%基準メッシュ生成関数の指定
ungradetest = ungradetest.setOptFlowCondition(0.0001,5,0,20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setHessianUpdate(eye(numel(orgVal)).*1,0.1,"SR1",12);


for i = 1:300
    ungradetest = ungradetest.updateVariables(@objFun,"direct",cmin,cmax);%設計変数微分の計算方法の指定
end
%ungradetest = ungradetest.solve

%}