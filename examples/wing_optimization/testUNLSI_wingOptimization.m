
%%%%%%%%ここから
clear;

orgVal = modifyDesFile("org.des","org.des");
%
CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
[~,~] = system(CommandVSP);
[~, orgSurf] = readvspgeom( "orgSurf.vspgeom", 0);

for iter = 1:3
    [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
    wing = UNLSI(p',con',id',wedata,1);
    wing.checkMesh(sqrt(eps));
    wing = wing.makeCluster(50,50);
    wing = wing.makeEquation(20,5,3);
    
    wing = wing.flowCondition(1,0.0001);
    wing = wing.flowCondition(2,4.0);
    wing = wing.setREFS(72,18,4);
    wing = wing.setRotationCenter([0,0,0]);
    wing = wing.setCf(1,500000,0.2,0.052*(10^-5),0);
    wing = wing.setCf(2,500000,0.2,0.052*(10^-5),0);
    wing = wing.solveFlow(1,10,3);
    wing.plotGeometry(1,wing.Cp,[-2,1]);
    wing = wing.calcApproximatedEquation();
    
    
    %設計変数勾配による表面近似
    ndim = numel(orgVal);
    designScale = [2,1,10,0.5,5,0.02,0.02];
    pert = 0.01;
    for i = 1:ndim
        sampleDes = orgVal;
        sampleDes(i) = orgVal(i) + pert*(designScale(i));
        modifyDesFile("org.des","mod.des",sampleDes);
        CommandVSP = strcat("vsp wing.vsp3 -script makeModSurf.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
        [~,~] = system(CommandVSP);
        [~, modSurf] = readvspgeom( "modSurf.vspgeom", 0);
        dmodSurf = modSurf'-orgSurf';
        sampleSurff = dmodSurf(:);
        sampleDes = orgVal;
        sampleDes(i) = orgVal(i) - pert*(designScale(i));
        modifyDesFile("org.des","mod.des",sampleDes);
        CommandVSP = strcat("vsp wing.vsp3 -script makeModSurf.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
        [~,~] = system(CommandVSP);
        [~, modSurf] = readvspgeom( "modSurf.vspgeom", 0);
        dmodSurf = modSurf'-orgSurf';
        sampleSurfr = dmodSurf(:);
        gradSurf(:,i) = (sampleSurff-sampleSurfr)./2./(pert*(designScale(i)));
    end
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5,'Display','iter');
    optVal = fmincon(@(x)testObj(x,wing,orgSurf',gradSurf),zeros(1,7),[],[],[],[],-designScale,designScale,[],options);
    orgVal = modifyDesFile("org.des","org.des",optVal+orgVal);
    CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
    [~,~] = system(CommandVSP);
    [~, orgSurf] = readvspgeom( "orgSurf.vspgeom", 0);
end
    
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
wing = UNLSI(p',con',id',wedata,1);
wing.checkMesh(sqrt(eps));
wing = wing.makeCluster(50,50);
wing = wing.makeEquation(20,5,3);

wing = wing.flowCondition(1,0.0001);
wing = wing.flowCondition(2,4.0);
wing = wing.setREFS(72,18,4);
wing = wing.setRotationCenter([0,0,0]);
wing = wing.setCf(1,500000,0.2,0.052*(10^-5),0);
wing = wing.setCf(2,500000,0.2,0.052*(10^-5),0);
wing = wing.solveFlow(1,10,3);
wing.plotGeometry(1,wing.Cp,[-2,1]);
      

%}