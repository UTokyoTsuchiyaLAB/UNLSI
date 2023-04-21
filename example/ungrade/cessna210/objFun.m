function [res,con] = objFun(x,AERODATA,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    res = sum(AERODATA{1}(:,11))*SREF/100;
    con(1,1) = AERODATA{1}(1,6)*SREF/100;
end