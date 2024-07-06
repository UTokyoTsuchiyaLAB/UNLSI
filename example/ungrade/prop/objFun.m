function [res,con] = objFun(x,AERODATA,DYNCOEF,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    rpm = 135;
    CT = -AERODATA{1}(15) * 0.5 * 7 ^2 * SREF / ((rpm/60)^2*BREF^4);
    Cq = AERODATA{1}(18) * 0.5 * 7 ^2 * SREF * BREF / ((rpm/60)^2*BREF^5);
    Cp = Cq * 2 * pi;
    J = 7/(rpm/60)/BREF;
    efficiency = CT*J/Cp;
    res = -efficiency;%効率の最大化
    con(1,1) = -AERODATA{1}(15) * 0.5 * 7 ^2 * SREF * 1.225/100;%推力の制約
end