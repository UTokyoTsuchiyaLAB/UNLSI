%
%%%%%%%%ここから
clear;
orgVal = [4.0000    9.0000         0    1.0000         0         0    0.1200];
%orgVal = modifyDesFile("org.des","org.des");
for i  = 1:5
    orgVal = modifyDesFile("org.des","org.des",orgVal);
    disp(orgVal)
    %
    CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
    [~,~] = system(CommandVSP);
    [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
    [orgSurf,orgCon] = vspSurfGen(orgVal,"wing","org.des");
    
    lb = [3, 8, 0,0.5,-5,-0.02,0.10];
    ub = [5,10,30,  1, 5, 0.02,0.14];
    unmeshtest = UNMESH(p',con',orgSurf,orgCon,orgVal,lb,ub);
    pert = [0.001,0.001,0.001,0.001,0.001,0.001,0.001];
    unmeshtest = unmeshtest.makeMeshGradient(@(x)vspSurfGen(x,"wing","org.des"));
    
    wing = UNLSI(p',con',id',wedata,1);
    wing.checkMesh(sqrt(eps));
    wing = wing.makeCluster(50,50);
    wing = wing.makeEquation(20,5,3);
    wing = wing.calcApproximatedEquation();
    %}
    [obj0, dobj_dx, con0, dcons_dx] = unmeshtest.calcObjandConsGradients(@(x)objFun(x,unmeshtest,wing));
    
    orgVal = orgVal - dobj_dx.*0.01;
    wing = wing.flowCondition(1,0.0001);
    wing = wing.setREFS(72,18,4);
    wing = wing.setRotationCenter([0,0,0]);
    wing = wing.setCf(1,500000,0.2,0.052*(10^-5),0);
    wing = wing.solveFlow(1,10,0);
    wing.plotGeometry(1,wing.Cp,[-2,1]);
end
      

%}