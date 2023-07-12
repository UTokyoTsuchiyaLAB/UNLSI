%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "optWing2.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata,0); %コンストラクタの実行 %動安定微係数用なので半裁にしない
wing.checkMesh(sqrt(eps)); %メッシュ状態のチェック
wing = wing.setProp(1,5,0.15,0);
wing = wing.makeCluster(50,50); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makePropEquation(100);
wing = wing.makeEquation(20,5,10); %パネル法行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}

wing = wing.flowCondition(1,0.0001); %flowNoのマッハ数を指定 

wing = wing.setREFS(0.379,1.204,0.357); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.setCf(1,100000,4,0.052*(10^-5),0); %摩擦係数パラメータの指定
wing = wing.setPropState(15,1.225,0.08,0.04,8000); %プロペラの回転状況を指定

wing = wing.setDeflAngle(6,[0,1,0],0);%6番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(7,[0,1,0],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflGroup(1,"elev",[6,7],[1,1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID1番） 
wing = wing.setDeflGroup(2,"ail",[6,7],[1,-1]); % ID6番7番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID2番）
wing = wing.calcDynCoef(1,0,0,[],0,1); %動安定微係数を計算（deflDerivFlagを1にすると動翼を微小量動かして数値微分を取ってくれる）
disp(wing.DYNCOEF{1,1});%DYNCOEFを表示

%空力係数の補間曲面の作成
massPropName = "optWing2_MassProps.txt";
wing = UNLSI2FEE(wing,massPropName,0,0.001,0,0,15);