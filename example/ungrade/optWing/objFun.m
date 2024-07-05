function [res,con] = objFun(x,AERODATA,DYNCOEF,Cp,Cfe,SREF,BREF,CREF,XYZREF,argin_x)
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    res = SREF*AERODATA{1}(1,11)*0.5*1.225*15^2; %抗揚比
    con(1,1) = SREF*AERODATA{1}(1,6)*0.5*1.225*15^2/9.8/10; %揚力条件
    xcg = -AERODATA{1}(1,19)/AERODATA{1}(1,6);%重心位置
    Cma = ((AERODATA{1}(2,19)-AERODATA{1}(1,19))+xcg*(AERODATA{1}(2,6)-AERODATA{1}(1,6)))/(3*pi/180);%静安定条件
    CLa = (AERODATA{1}(2,6)-AERODATA{1}(1,6))/(3*pi/180);%静安定条件
    con(2,1) = Cma/CLa*10;
    %con(2,1) = x(4)-x(5);
    %con(3,1) = x(5)-x(1);
    %con(4,1) = x(1)-x(3);
end