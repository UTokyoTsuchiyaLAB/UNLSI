function [res,con] = objFun(x,AERODATA,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CDo 7:CDi 8:CDtot 9:CDt 10:CDtot_t 11:CS 12:L/D E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    res = AERODATA(7)*100;
    %con = AERODATA(5);
    con(1,1) = AERODATA(5);
    con(2,1) = x(4)-x(2);
    con(3,1) = x(5)-x(4);
    con(4,1) = x(1)-x(5);
    con(5,1) = x(3)-x(1);
end