%
%%%%%%%%ここから
clear;
winddata.alpha = -4:2:16; %解析迎角
winddata.CL = [-0.15,0.03,0.2,0.38,0.57,0.74,0.9,1.08,1.24,1.37,1.49];%風洞試験結果 CL
winddata.CD = [0.03,0.02,0.03,0.04,0.05,0.07,0.09,0.11,0.13,0.15,0.19];%風洞試験結果 CD

[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "Cessna-210.vspgeom", 0);%形状の読み込み
wing = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
wing = wing.setREFS(175,36.75,4.91,[0,0,0]); %基準面積 基準長の設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}

ad = [];
for i = 1:numel(winddata.alpha)
    wing = wing.solveFlow(winddata.alpha(i),0,0.001,2.3*10^6); %パネル法を解く
    if winddata.alpha(i) == 0
        wing.plotGeometry(1,wing.getCp(0,0,0.001,2.3*10^6),[-2,1]); %圧力係数のプロット
    end
    ad = [ad;wing.AERODATA{1}];
end


%%%%%風洞試験結果との比較
figure(2);clf,grid on;hold on;
set(gca,"FontSize",14);
plot(winddata.alpha,winddata.CL,'-o');
plot(winddata.alpha,ad(:,5),'-o');
plot(winddata.alpha,ad(:,6),'-o');
legend("WT","Cp Integral","Trefftz")
xlabel("AoA deg");
ylabel("CL");
figure(3);clf,grid on;hold on;
set(gca,"FontSize",14);
plot(winddata.alpha,winddata.CD,'-o');
plot(winddata.alpha,ad(:,9),'-o');
plot(winddata.alpha,ad(:,11),'-o');
legend("WT","Cp Integral","Trefftz")
xlabel("AoA deg");
ylabel("CD");