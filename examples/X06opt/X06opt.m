%
%%%%%%%%ここから
clear;
orgVal = [1.9,58,0,2.9,2.1,50,1.0,0];
%orgVal = modifyDesFile("org.des","org.des");

ub = [2.5,65, 5, 3.5,2.5,65,3.0,5];
lb = [1.5,50,-5, 2.0,1.5,50,0.5,-5];
cmin = [0.1]';
cmax = [0.2]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal,"X06","org.des");
ungradetest = UNGRADE(@(x)vspSurfGen(x,"X06","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"X06","org.des"));%基準メッシュ生成関数の指定
ungradetest = ungradetest.setOptFlowCondition(0.0001,5,0,20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptimization(eye(numel(orgVal)).*1,0.1,0.001,1,"SSR1",12);

%ungradetest = ungradetest.updateVariables(@objFun,"nonlin",cmin,cmax);%設計変数微分の計算方法の指定
for i = 1:300
    [dx,ungradetest] = ungradetest.updateVariables(@objFun,"direct",cmin,cmax);%設計変数微分の計算方法の指定
    if norm(dx./ungradetest.designScale) < 0.001
        ungradetest = ungradetest.updateVariables(@objFun,"nonlin",cmin,cmax);%設計変数微分の計算方法の指定
    end
end
%ungradetest = ungradetest.solve

%}