function [res,con] = objFun(x,unmesh,unlsi)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CDo 7:CDi 8:CDtot 9:CDt 10:CDtot_t 11:CS 12:L/D E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    modSurf = unmesh.makeSurffromVariables(x);
    modifiedVerts = unmesh.meshDeformation(modSurf);
    approxunlsi = unlsi.makeAproximatedInstance(modifiedVerts);
    
    approxunlsi = approxunlsi.flowCondition(1,0.0001);
    approxunlsi = approxunlsi.setREFS(72,18,4);
    approxunlsi = approxunlsi.setRotationCenter([0,0,0]);
    approxunlsi = approxunlsi.setCf(1,500000,0.2,0.052*(10^-5),0);
    approxunlsi = approxunlsi.solveFlow(1,10,0);
    res = approxunlsi.AERODATA(8);
    con = approxunlsi.AERODATA(5);
end