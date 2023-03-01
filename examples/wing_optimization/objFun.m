function [res,con] = objFun(x,unmesh,unlsi)
    
    modSurf = unmesh.makeSurffromVariables(x);
    modifiedVerts = unmesh.meshDeformation(modSurf);
    approxunlsi = unlsi.makeAproximatedInstance(modifiedVerts);
    
    approxunlsi = approxunlsi.flowCondition(1,0.0001);
    approxunlsi = approxunlsi.setREFS(72,18,4);
    approxunlsi = approxunlsi.setRotationCenter([0,0,0]);
    approxunlsi = approxunlsi.setCf(1,500000,0.2,0.052*(10^-5),0);
    approxunlsi = approxunlsi.solveFlow(1,10,0);
    res = -approxunlsi.AERODATA(12);
    con = [];
end