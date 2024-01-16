function [res,con] = objFun(x,AERODATA,DYNCOEF,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    res = AERODATA{1}(1,10)*SREF;
    con(1,1) = AERODATA{1}(1,6)*SREF/10;
    con(2,1) = x(2)-x(3);
    con(3,1) = x(3)-x(5);
    con(4,1) = x(5)-x(6);
    con(5,1) = x(6)-x(1);
    con(6,1) = x(1)-x(4);
end