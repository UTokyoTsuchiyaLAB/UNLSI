function [res,con] = objFun(x,valid,unmesh,unlsi)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CDo 7:CDi 8:CDtot 9:CDt 10:CDtot_t 11:CS 12:L/D E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    modSurf = unmesh.makeSurffromVariables(x);
    modifiedVerts = unmesh.meshDeformation(modSurf);
    approxunlsi = unlsi.makeAproximatedInstance(modifiedVerts);
    if valid == 1
        approxunlsi = approxunlsi.makeEquation(20,5,3);
    end
    approxunlsi = approxunlsi.solveFlow(1,5,0);
    res = approxunlsi.AERODATA(7)*10;
    con(1,1) = approxunlsi.AERODATA(5);
    %con(2,1) = x(3)-x(2);
    %con(3,1) = x(5)-x(3);
    %con(4,1) = x(6)-x(5);
    %con(5,1) = x(1)-x(6);
    %con(6,1) = x(4)-x(1);
end