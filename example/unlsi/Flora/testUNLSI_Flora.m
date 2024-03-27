%
%%%%%%%%ここから
clear;tic
%%%%%%%%%
%パイプを等価な剛性を持つspar2本に変換
[tspar,b,Rarea] = InatiaPipe2biRec(30,1.5,50); %mm単位

%%%%%%%%%
Vinf = 10;%7.2;
Re = Vinf * 0.375 /1.512 *1e5;
alpha = 0;%[-3,0,3];
beta = 0;

massPropName = "flora_MassProps.txt";
[mass,cgPoint,Inatia] = readMassPropResult(massPropName);
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "flora.vspgeom", 0);
% [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "FWpeller2.vspgeom", 0);
wing = UNLSI(p',con',id',wedata,1);
wing = wing.setREFS(3.675,9.6,0.375);
wing = wing.setRotationCenter([-0.1,0,0]);%[-0.1,0,0]
wing = wing.setUNLSISettings("propCalcFlag",0,"laminarRatio",0.4,"propWakeLength",1,"Vinf",Vinf,"nGriddedInterp",10);
% wing = wing.setProp(1,4,1.6,1);

% wing = wing.setDeflAngle(6,[0,1,0],0);%6番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
% wing = wing.setDeflAngle(7,[0,1,0],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
% wing = wing.setDeflAngle(8,[0,0,1],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
% wing = wing.setDeflGroup(1,"elev",[6,7],[1,1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID1番） 
% wing = wing.setDeflGroup(2,"ail",[6,7],[1,-1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID2番wing = wing.setDeflGroup(2,"ail",[6,7],[1,-1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID2番
% wing = wing.setDeflGroup(3,"rud",8,1); 
% wing = wing.solveFlow(0,0,0.001,200000);
wing = wing.makeCluster();
wing = wing.makeEquation();
%}
%
% disp('start')
[confem,vertsfem,idfem] = readFemMesh('flora_WingGeom_Struct0.msh'); %FEMメッシュを読み込み
% [confem,vertsfem,idfem] = readFemMesh('FWpeller2_WingGeom_Struct0.msh'); %FEMメッシュを読み込み

wing = wing.setFemMesh(vertsfem,confem,idfem,[1 3 5]); %空力メッシュID1番をFEMメッシュに対応付ける
% [wing,weight] = wing.setFemMaterials([1,2,3,4],[0.0007,0.01,tspar/1000,tspar/1000],[1313000000,3400000000/37,400000000000*0.33,400000000000*0.33],[10,100,1600*Rarea,1600*Rarea]);%物性値をセット 肉厚,ヤング率,密度,減衰パラメータ
[wing,weight] = wing.setFemMaterials([1,2,3,4],[0.0007,0.001,tspar/1000,tspar/1000],[3000000000,3000000000,40*1e9*0.33,40*1e9*0.33],[100,100,6400*Rarea,6400*Rarea],[1000,1000,1000,1000]);%skin,rib*4,spar,spar
disp("weight");
disp(weight);
wing = wing.makeFemEquation(); %FEM行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
% wing = wing.setPropState(1,0.4,0.8,160); %プロペラの回転状況を指定
wing = wing.solveFlow(alpha,beta,0.001,Re); %パネル法を解く
wing.plotGeometry(1,wing.getCp(0),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha));
wing2 = wing;%変位後の空力解析用のインスタンスを準備

%
%空力係数の補間曲面の作成
% wing = wing.makeSurrogateModel([-10,-5,0,5,10],[-10,-5,0,5,10],0.0001,300000);%
% [ppCoef,ppDyn] = wing.getSurrogateModel();
% wing.plotSurrogateModel(ppCoef,ppDyn,linspace(-10,10,100),linspace(-10,10,100),0.001,1);
% UREF = Vinf;
% Mach = 0.001;
% REFS = [wing.SREF,wing.BREF,wing.CREF];
% save aeroRBF.mat ppCoef ppDyn REFS mass Inatia UREF Mach -mat;
%

%
dt = 0.01;
% VideoWriter オブジェクトを作成
v = VideoWriter('AeroAnalysis4.mp4', 'MPEG-4');
% 時間区切りからフレームレートの計算と適用
framerate = 1.0/dt;
v.FrameRate = framerate;
% 保存する動画の画質。数字の大きいほうが高画質.[0~100]
v.Quality = 100;
% ビデオの書き込みを開始
open(v);


%     % 最後のflag0→機体に働く重力無視，後ろから二番目の荷重を100倍とかにして変形が普通じゃない→lsqminnormの精度不足
%     % solveaeroelasticに変えればdtでるのでcm
%     % 空力弾性やりたいときはUNLSIのクラスをコピーしてそこに変形後の空力メッシュを入れて積分していくその時剛性マトリクスは変形後にするとそこが中立点になってしまうので元のものを入れる
%     % でも今の変形は半裁で中央を抑えているけれど実際はちゅうおうもグネグネしているので運動方程式との練成を考える必要がある．




tic
for i = 1:5
    disp(i);
    % tic;
    if i == 1
        [delta,deltadot]  = wing2.solveAeroelastic([0,dt],[],[],wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1); %空力弾性計算　初期値0
    else
        [delta,deltadot]  = wing2.solveAeroelastic([0,dt],delta,deltadot,wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%空力弾性計算
    end
    % toc;
    modVerts = wing.calcModifiedVerts(delta{1}); %構造メッシュの変形に従って空力メッシュを変形
    wing2 = wing2.setVerts(modVerts); %節点の移動のみ
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(2,wing2.getCp(alpha,0,0.001,Re),[-3,1.5]);
    M(i) = getframe; %動画用のフレームを取得
     % フレームを作成
    frame = getframe(gcf);
    % フレームをビデオに書き込み
    writeVideo(v, frame);
end
disp('end')
movie(M,1,1/dt);
close(v);
toc
%}