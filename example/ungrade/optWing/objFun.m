function [res,con] = objFun(geom,aero,fem)
        %geom:機体構造関連
    %   geom.x:設計変数
    %   geom.aeroMesh:空力解析の解析メッシュ
    %   geom.surfID:ID
    %   geom.SREF,BREF,CREF,XYZREF,argin_x
    %aero:空力解析関連
    %   aero.DATA;空力解析結果
    %       1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    %   aero.DYNCOEF:動安定微係数（フラグを1にした場合のみ）
    %   aero.Cp([x,y,z]):CpをRBFN補間する関数
    %   aero.Cfe([x,y,z]):CfeをRBFN補間する関数
    %fem:構造解析関連
    %   fem.mesh:構造解析メッシュ
    %   fem.delta([x,y,z]):変位をRBFN補間する関数。空力メッシュに投影したものなので注意
    %
    lepos = geom.x(end)-0.12+0.55*tand(geom.x(6));
    deformation = fem.delta{1,1}([lepos,0.7,0]);
    res = geom.SREF*aero.DATA{1}(1,11)*0.5*1.225*15^2+deformation(1,3)^2*10000; %抗揚比
    con(1,1) = geom.SREF*aero.DATA{1}(1,6)*0.5*1.225*15^2/9.8/10; %揚力条件
    xcg = -aero.DATA{1}(1,19)/aero.DATA{1}(1,6);
    con(2,1) = -aero.DYNCOEF{1}(3,3)*xcg+aero.DYNCOEF{1}(5,3)*10;
%     xcg = -AERODATA{1}(1,19)/AERODATA{1}(1,6);%重心位置
%     Cma = ((AERODATA{1}(2,19)-AERODATA{1}(1,19))+xcg*(AERODATA{1}(2,6)-AERODATA{1}(1,6)))/(3*pi/180);%静安定条件
%     CLa = (AERODATA{1}(2,6)-AERODATA{1}(1,6))/(3*pi/180);%静安定条件
%     con(2,1) = Cma/CLa*10;
    %con(2,1) = x(4)-x(5);
    %con(3,1) = x(5)-x(1);
    %con(4,1) = x(1)-x(3);
end