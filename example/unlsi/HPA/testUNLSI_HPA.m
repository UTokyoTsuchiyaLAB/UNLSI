%
%%%%%%%%ここから
clear;massPropName = "FWpeller2_MassProps.txt";
[mass,cgPoint,Inatia] = readMassPropResult(massPropName);
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "FWpeller2.vspgeom", 0);

wing = UNLSI(p',con',id',wedata,0,0.001);
wing = wing.setREFS(25.44,33.75,0.766,cgPoint);
wing = wing.setUNLSISettings("laminarRatio",0.4,"nWakeMax",30,"Vinf",7.2);
wing = wing.setProp(1,6,1.6,0);
wing = wing.setCpLimit([-2.2,1]);
wing = wing.makeCluster();
wing = wing.makeEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
wing = wing.setPropState(1,0.4,0.8,160); %プロペラの回転状況を指定
wing = wing.setDeflAngle(3,[0,1,0],0);%6番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(4,[0,1,0],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(5,[0,0,1],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflGroup(1,"elev",[3,4],[1,1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID1番） 
wing = wing.setDeflGroup(2,"rud",[5],[1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID2番
alpha = [-3,0,3];
beta = [0,0,0];
wing = wing.solveFlow(alpha,beta,0.001,300000);
wing.plotGeometry(1,wing.getCp(0,0,0.001,300000),[-2,1]);%圧力係数のプロット
disp(wing.getAERODATA(alpha,beta,0.001,300000));
%
%空力係数の補間曲面の作成
wing = wing.makeSurrogateModel([-10,-5,0,5,10],[-10,-5,0,5,10],0.0001,300000);
[ppCoef,ppDyn] = wing.getSurrogateModel();
wing.plotSurrogateModel(ppCoef,ppDyn,linspace(-10,10,100),linspace(-10,10,100),0.001,1);
UREF = 7.2;
Mach = 0.001;
REFS = [wing.SREF,wing.BREF,wing.CREF];
save aeroRBF.mat ppCoef ppDyn REFS mass Inatia UREF Mach -mat;
%}
