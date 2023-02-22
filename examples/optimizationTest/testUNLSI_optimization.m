
%
%%%%%%%%ここから
clear;

%opt.desの書き換え

CommandVSP = strcat("vsp opt.vsp3 -script opt.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
[~,~] = system(CommandVSP);
%}
%


[con3, p3, uv1, uv2, uv3, wedata, id] = readvspgeom( "opt.vspgeom", 0);
anemoi = UNLSI(p3',con3',id',wedata,1);
anemoi.checkMesh(sqrt(eps));
anemoi = anemoi.makeCluster(50,50);
anemoi.plotGeometry(2);

%メッシュ変形
[con1, p1, uv1, uv2, uv3, wedata, id] = readvspgeom( "des.vspgeom", 0);
[con2, p2, uv1, uv2, uv3, wedata, id] = readvspgeom( "org.vspgeom", 0);
Ftest.x = scatteredInterpolant(p2',p1(1,:)'-p2(1,:)','linear','linear');
Ftest.y = scatteredInterpolant(p2',p1(2,:)'-p2(2,:)','linear','linear');
Ftest.z = scatteredInterpolant(p2',p1(3,:)'-p2(3,:)','linear','linear');
p1dvertex(:,1) = Ftest.x(p3');
p1dvertex(:,2) = Ftest.y(p3');
p1dvertex(:,3) = Ftest.z(p3');

apxVerts = p3'+ p1dvertex;
anemoi = anemoi.setVerts(apxVerts);
anemoi.plotGeometry(3);
%}
%

anemoi = anemoi.makeEquation(20,5,10);
%%%%%%%%ここまでは一度インスタンスを作成したらスキップできる
%

%
anemoi = anemoi.flowCondition(1,0.0001);
anemoi = anemoi.setREFS(0.8,4,0.2);
anemoi = anemoi.setRotationCenter([-0.07,0,0]);
anemoi = anemoi.setCf(1,50000,0.2,0.052*(10^-5),0,1.5);
anemoi = anemoi.solveFlow(1,0,0);
anemoi.plotGeometry(1,anemoi.Cp,[-2,1]);
%}