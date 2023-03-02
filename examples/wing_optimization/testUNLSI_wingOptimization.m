%
%%%%%%%%ここから
clear;
%orgVal = [4.0000    9.0000         0    1.0000         0         0    0.1200];
orgVal = modifyDesFile("org.des","org.des");

nextVal = orgVal;
ub = [  4,   4,   4,   4,   4,   5];
lb = [0.5, 0.5, 0.5, 0.5, 0.5,  -5];
unmeshtest = UNMESH(lb,ub);


%}
for i  = 1:20
    fclose all;
    pause(1);
    [~,nextVal] = modifyDesFile("org.des","org.des",nextVal);
    pause(1);
    disp(nextVal);
    %
    CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript -des org.des"); %手順2のコマンド環境に合わせて各自書き換え
    [~,~] = system(CommandVSP);
    [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
    [orgSurf,orgCon] = vspSurfGen(nextVal,"wing","org.des");
    

    unmeshtest = unmeshtest.updateMesh(p',con',orgSurf,orgCon,nextVal);

    wing = UNLSI(p',con',id',wedata,1);
    wing.checkMesh(sqrt(eps));
    wing = wing.makeCluster(50,50);
    wing = wing.makeEquation(20,5,3);
    wing = wing.flowCondition(1,0.0001);
    wing = wing.setREFS(80,20,4);
    wing = wing.setRotationCenter([0,0,0]);
    wing = wing.setCf(1,500000,0.2,0.052*(10^-5),0);
    wing = wing.solveFlow(1,5,0);
    wing.plotGeometry(1,wing.Cp,[-2,1]);

    wing = wing.calcApproximatedEquation();
    %}
    unmeshtest = unmeshtest.makeMeshGradient(@(x)vspSurfGen(x,"wing","org.des"));
    unmeshtest = unmeshtest.calcObjandConsGradients(@(x)objFun(x,unmeshtest,wing),[0.5],[0.55]);
    [dx,unmeshtest] = unmeshtest.updateVariables();
    
    nextVal = nextVal + dx;
end
      

%}