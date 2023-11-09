function [res,con] = objFun(x,AERODATA,DYNCOEF,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    res = SREF*AERODATA{1}(1,11)*0.5*1.225*15^2;
    con(1,1) = SREF*AERODATA{1}(1,6)*0.5*1.225*15^2/9.8;
    xcg = -AERODATA{1}(1,19)/AERODATA{1}(1,6);
    con(2,1) = (-DYNCOEF{1,1}(3,3)*xcg+DYNCOEF{1,1}(5,3))*100;
    con(3,1) = x(1)+x(8);
    %con(1,1) = x(2)-x(4);
    %con(2,1) = x(4)-x(5);
    %con(3,1) = x(5)-x(1);
    %con(4,1) = x(1)-x(3);
end