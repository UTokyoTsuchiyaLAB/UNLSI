%
%%%%%%%%ここから
clear;
orgVal = [0.04,0.4,0.12,0,0.04,0.4,0.12,5.5,0,4,2];
%orgVal = modifyDesFile("org.des","org.des");

lb = [-0.05,0.2,0.1,0,-0.05,0.2,0.1,4.0,0,2,0];
ub = [ 0.05,0.5,0.14,6,0.05,0.5,0.14,8,20,8,6];
cmin = [40/100, 0.2 -0.3]';
cmax = [45/100, 1.0  0.3]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal,"Cessna-210","org.des");
ungradetest = UNGRADE(@(x)vspSurfGen(x,"Cessna-210","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
ungradetest.checkSurfGenWork(0.5);
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"Cessna-210","org.des"));%基準メッシュ生成関数の指定
ungradetest = ungradetest.setOptFlowCondition([0.001,0.001,0.001],[0,5,10],[0,0,0],20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptimization(eye(numel(orgVal)).*1,0.3,0.001,0.6);

for i = 1:300
    [dx,convergenceFlag,ungradetest] = ungradetest.updateVariables(@objFun,0.0001,"direct","BFGS",20,cmin,cmax);%設計変数微分の計算方法の指定
    ungradetest.plotOptimizationState(2);
    if convergenceFlag == 1
        break;
    end
end

%}