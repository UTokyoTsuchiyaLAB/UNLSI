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
    con(1,1) = aero.DATA{1}(1,6)*geom.SREF/10;
%     chord = [0,geom.x(2);
%              2,geom.x(3);
%              4,geom.x(5);
%              6,geom.x(6);
%              8,geom.x(1);
%              10,geom.x(4)];
    %[L,R,k] = curvature(chord);
    res = -aero.DATA{1}(1,14);
    con(2,1) = geom.x(2)-geom.x(3);
    con(3,1) = geom.x(3)-geom.x(5);
    con(4,1) = geom.x(5)-geom.x(6);
    con(5,1) = geom.x(6)-geom.x(1);
    con(6,1) = geom.x(1)-geom.x(4);
    deformation = fem.delta{1,1}([0,10,0]);
    con(7,1) = deformation(1,3);
end

function [R,M,k] = circumcenter(A,B,C)
% Center and radius of the circumscribed circle for the triangle ABC
%  A,B,C  3D coordinate vectors for the triangle corners
%  R      Radius
%  M      3D coordinate vector for the center
%  k      Vector of length 1/R in the direction from A towards M
%         (Curvature vector)
  D = cross(B-A,C-A);
  b = norm(A-C);
  c = norm(A-B);
  if nargout == 1
    a = norm(B-C);     % slightly faster if only R is required
    R = a*b*c/2/norm(D);
    if norm(D) == 0
      R = Inf;
    end
    return
  end
  E = cross(D,B-A);
  F = cross(D,C-A); 
  G = (b^2*E-c^2*F)/norm(D)^2/2;
  M = A + G;
  R = norm(G);  % Radius of curvature
  if R == 0
    k = G;
  elseif norm(D) == 0
    R = Inf;
    k = D;
  else
    k = G'/R^2;   % Curvature vector
  end
end

function [L,R,k] = curvature(X)
% Radius of curvature and curvature vector for 2D or 3D curve
%  [L,R,k] = curvature(X)
%   X:   2 or 3 column array of x, y (and possibly z) coordiates
%   L:   Cumulative arc length
%   R:   Radius of curvature
%   k:   Curvature vector
% The scalar curvature value is 1./R
% Version 2.6: Calculates end point values for closed curve

  N = size(X,1);
  dims = size(X,2);
  if dims == 2
    X = [X,zeros(N,1)];  % Use 3D expressions for 2D as well
  end
  L = zeros(N,1);
  R = NaN(N,1);
  k = NaN(N,3);
  for i = 2:N-1
    [R(i),~,k(i,:)] = circumcenter(X(i,:)',X(i-1,:)',X(i+1,:)');
    L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  end
  if norm(X(1,:)-X(end,:)) < 1e-10 % Closed curve. 
    [R(1),~,k(1,:)] = circumcenter(X(end-1,:)',X(1,:)',X(2,:)');
    R(end) = R(1);
    k(end,:) = k(1,:);
    L(end) = L(end-1) + norm(X(end,:)-X(end-1,:));
  end
  i = N;
  L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  if dims == 2
    k = k(:,1:2);
  end
end
