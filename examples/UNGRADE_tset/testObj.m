function res = objFun(x,unlsi,orgVal,orgSurf,gradSurf)
    
    modSurf = orgSurf + reshape(gradSurf*(x(:)-orgVal(:)),size(orgSurf));
    modifiedVerts = unlsi.meshDeformation(orgSurf,modSurf);
    approxunlsi = unlsi.makeAproximatedInstance(modifiedVerts);
    
    approxunlsi = approxunlsi.flowCondition(1,0.0001);
    approxunlsi = approxunlsi.flowCondition(2,4.0);
    approxunlsi = approxunlsi.setREFS(72,18,4);
    approxunlsi = approxunlsi.setRotationCenter([0,0,0]);
    approxunlsi = approxunlsi.setCf(1,500000,0.2,0.052*(10^-5),0);
    approxunlsi = approxunlsi.setCf(2,500000,0.2,0.052*(10^-5),0);
    approxunlsi = approxunlsi.solveFlow(1,10,3);
    res = -approxunlsi.AERODATA(12);
end