clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);


wing = UNLSI(p',con',id',wedata,1);
wing.checkMesh(sqrt(eps));
wing = wing.makeCluster(50,50);
wing = wing.makeEquation(20,5,3);
wing = wing.calcApproximatedEquation();

orgVal = [4.0000,4.0000, 4.0000,4.0000,4.0000];

nextVal = orgVal;
ub = [  4,   4,   4,   4,   4];
lb = [0.5, 0.5, 0.5, 0.5, 0.5];
cmin = [ 0.4]';
cmax = [0.45]';
unmeshtest = UNMESH(lb,ub);
[orgSurf,orgCon] = vspSurfGen(orgVal,"wing","org.des");
unmeshtest = unmeshtest.updateMesh(p',con',orgSurf,orgCon,nextVal);
unmeshtest = unmeshtest.makeMeshGradient(@(x)vspSurfGen(x,"wing","org.des"));

pert = logspace(-9,-6,100);
for i = 1:numel(pert)
    x = orgVal;
    x(2) = orgVal(2)+pert(i);
    modSurf = unmeshtest.makeSurffromVariables(x);
    modifiedVerts = unmeshtest.meshDeformation(modSurf);
    wing2 = UNLSI(modifiedVerts,con',id',wedata,1);
    wing2.checkMesh(sqrt(eps));
    wing2 = wing2.makeCluster(50,50);
    wing2 = wing2.makeEquation(20,5,3);
    val(:,i) = wing2.LHS(:,1);
    approxunlsi = wing.makeAproximatedInstance(modifiedVerts);
    val2(:,i) = approxunlsi.LHS(:,1);
end
%{
wing = wing.calcApproximatedEquation();
%modifiedVerts = wing.tri.Points;
modifiedVerts(1,1) = modifiedVerts(1,1)+0.0001;
%approxunlsi = wing.makeAproximatedInstance(modifiedVerts);

wing2 = UNLSI(modifiedVerts,con',id',wedata,1);
wing2 = wing2.makeCluster(50,50);
wing2 = wing2.makeEquation(20,5,3);

errorLHS = approxunlsi.LHS-wing2.LHS;
errorRHS = approxunlsi.RHS-wing2.RHS;
%}