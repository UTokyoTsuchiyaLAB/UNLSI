%
%%%%%%%%ここから
clear;massPropName = "FWpeller2_MassProps.txt";
[mass,cgPoint,Inatia] = readMassPropResult(massPropName);
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "FWpeller2.vspgeom", 0);

wing = UNLSI(p',con',id',wedata,0,0.001);
wing = wing.setREFS(25.44,33.75,0.766,cgPoint);
wing = wing.setUNLSISettings("laminarRatio",0.4,"nWakeMax",30,"Vinf",7.2);
wing = wing.setCpLimit([-2.2,1]);%全体のCpにリミットをかける
wing = wing.setCpLimit([-20,1],6:7); %ID6と7のCpにリミットをかける（上書き）
wing = wing.makeCluster();
wing = wing.makeEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
wing = wing.setDeflAngle(3,[0,1,0],0);%6番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(4,[0,1,0],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(5,[0,0,1],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflGroup(1,"elev",[3,4],[1,1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID1番） 
wing = wing.setDeflGroup(2,"rud",[5],[1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID2番
alpha = [-3,0,3];
beta = [0,0,0];
%}
%
rpm = 135;
dt = 1/20;
wing = wing.setWakeShape([0.1,0,0]);
for i = 1:100
    disp(i)
    [wing,CT,Cp,Cq,efficiency] = wing.solveUnsteadyProp(6:7,dt,rpm,[1,0,0],[0,0,0],0,0,0.001,300000,3,[-2,1]);
    disp([CT*(rpm/60)^2*wing.BREF^4*1.225,Cp*(rpm/60)^3*wing.BREF^5*1.225,Cq*(rpm/60)^2*wing.BREF^5*1.225,efficiency]);
    % フレームを作成
    M(i) = getframe(gcf);
    % フレームをビデオに書き込み
    videoMaker(M,"hpaCalculating",1/dt)
end
%}
%{
%%%定常プロペラ計算
%wakeと主翼の位置関係によって大きく空力係数が変動するので注意
wing = wing.rotateVerts(6:7,60,[1,0,0],[0,0,0]);
wing = wing.makeEquation();
[wing,CT,Cp,Cq,efficiency] = wing.solveSteadyProp(6:7,140,[1,0,0],[0,0,0],0,0,0.001,300000,2,[-2,1]);
disp([CT*(rpm/60)^2*wing.BREF^4*1.225,Cp*(rpm/60)^3*wing.BREF^5*1.225,Cq*(rpm/60)^2*wing.BREF^5*1.225,efficiency]);
wing.plotGeometry(1,wing.getCp(0,0,0.001,300000),[-2,1]);%圧力係数のプロット
disp(wing.getAERODATA(alpha,beta,0.001,300000));
%}

