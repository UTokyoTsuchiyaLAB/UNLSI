clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
wing = UNLSI(p',con',id',wedata,1);
wing.checkMesh(sqrt(eps));
wing = wing.makeCluster(50,50);
wing = wing.makeEquation(20,5,3);

wing = wing.calcApproximatedEquation();
modifiedVerts = wing.tri.Points;
modifiedVerts(1,1) = modifiedVerts(1,1)+0.1;
approxunlsi = wing.makeAproximatedInstance(modifiedVerts);

wing2 = UNLSI(modifiedVerts,con',id',wedata,1);
wing2 = wing2.makeCluster(50,50);
wing2 = wing2.makeEquation(20,5,3);

%approxunlsi.LHS-wing2.LHS
%approxunlsi.RHS-wing2.RHS