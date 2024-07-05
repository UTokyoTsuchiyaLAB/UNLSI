%楕円翼を計算する
b = 20;
g0 = 3.7330;
y = linspace(0,10,100);
gamma = g0 * sqrt(1-(y./(b/2)).^2);
CL = 1;
chord = gamma/CL;
%翼弦長分布
figure(3);clf;hold on;ylim([0,4]);grid on;
plot([fliplr(-y),y],[fliplr(chord),chord],"r--","LineWidth",2);
%var = modifyDesFile("org.des");
yvsp = [0,2,4,6,8,10];
chordvsp = ungradetest.unscaledVar([2,3,5,6,1,4]);
chordvsp0 = ungradetest2.unscaledVar([2,3,5,6,1,4]);
plot([fliplr(-yvsp),yvsp],[fliplr(chordvsp0),chordvsp0],"b-.","LineWidth",2);
plot([fliplr(-yvsp),yvsp],[fliplr(chordvsp),chordvsp],"b-","LineWidth",2);
set(gca,"FontSize",14);ylabel("chord");xlabel("spanwise");
legend("Analytical","Initial","Optimized");

%solveFlow
ungradetest = ungradetest.solveFlow(linspace(0,10,20),zeros(1,20),0.001,500000);
ungradetest2 = ungradetest2.solveFlow(linspace(0,10,20),zeros(1,20),0.001,500000);
%解析解
AR = ungradetest.BREF^2/ungradetest.SREF;
CLa = 2*pi/(1+2/AR);
CLi = CLa * (linspace(0,10,20).*pi/180);
CDi = CLi .^ 2 /(pi*1*AR);
%CL
figure(4);clf;hold on;grid on;
plot(ungradetest.AERODATA{1}(:,3),CLi,"r--","LineWidth",2);
plot(ungradetest2.AERODATA{1}(:,3),ungradetest2.AERODATA{1}(:,6),"b-.","LineWidth",2);
plot(ungradetest.AERODATA{1}(:,3),ungradetest.AERODATA{1}(:,6),"b-","LineWidth",2);
set(gca,"FontSize",14);ylabel("CL");xlabel("AoA(deg)");
legend("Analytical","Initial","Optimized");
%CD
figure(5);clf;hold on;grid on;
plot(ungradetest.AERODATA{1}(:,3),CDi,"r--","LineWidth",2);
plot(ungradetest2.AERODATA{1}(:,3),ungradetest2.AERODATA{1}(:,10),"b-.","LineWidth",2);
plot(ungradetest.AERODATA{1}(:,3),ungradetest.AERODATA{1}(:,10),"b-","LineWidth",2);
set(gca,"FontSize",14);ylabel("CDt");xlabel("AoA(deg)");
legend("Analytical","Initial","Optimized");