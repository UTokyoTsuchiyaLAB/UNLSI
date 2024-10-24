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
    % 1:Beta 2:Mach 3:AoA 4:Re/1e6 5:CL 6:CLt  7:CDo 8:CDi 9:CDtot 10:CDt 11:CDtot_t 12:CS 13:L/D 14E CFx CFy CFz CMx CMy       CMz       CMl       CMm       CMn      FOpt 
    rpm = 135;
    CT = -aero.DATA{1}(15) * 0.5 * 7 ^2 * geom.SREF / ((rpm/60)^2*geom.BREF^4);
    Cq = aero.DATA{1}(18) * 0.5 * 7 ^2 * geom.SREF * geom.BREF / ((rpm/60)^2*geom.BREF^5);
    Cp = Cq * 2 * pi;
    J = 7/(rpm/60)/geom.BREF;
    efficiency = CT*J/Cp;
    res = -efficiency;
    con(1,1) = -aero.DATA{1}(15) * 0.5 * 7 ^2 * geom.SREF * 1.225/100;
    con(2,1) = geom.x(8) * geom.x(7)*100;
    con(3,1) = geom.x(11) * 0.090*100;
end