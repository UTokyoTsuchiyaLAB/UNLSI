function vehicle = NAL0Design(p,wgsname,solver)
%NAL0形状の作成

para.wh = p(1); %胴体幅
para.hu = p(2);%胴体上部高さ
para.hl = p(3);%胴体下部高さ
para.ln = p(4);%ノーズ長さ
para.Cr = p(5);%翼根翼弦長
para.Ck = p(6);%キンク翼弦長
para.Ct = p(7);%翼端翼弦長
para.semispan = p(8);%セミスパン
para.Rkink = p(9);%キンク割合
para.hvert = p(10);%垂直尾翼高さ
para.Crv = p(11);%垂直尾翼翼根
para.Ctv = p(12);%垂直尾翼翼端
para.phivert = p(13);%垂直尾翼傾き
para.HR = p(14);%拡大比

%ノーズの作成
n_nose = 15;
s_body = linspace(pi/2,0,7);
[Haa] = Haack(100);
x_nose = linspace(0,para.ln,n_nose)'; 
%胴体上部楕円
%rel = sqrt((para.wh.*cos(s_body)).^2+(para.hu.*sin(s_body)).^2);
for i = 1:n_nose
    vehicle.surface.nose{1}(i,:) = x_nose(i).*ones(size(s_body));
    vehicle.surface.nose{2}(i,:) = (interp1(Haa(:,1).*para.ln,Haa(:,2),x_nose(i)).*para.wh).*cos(s_body);
    vehicle.surface.nose{3}(i,:) = (interp1(Haa(:,1).*para.ln,Haa(:,2),x_nose(i)).*para.wh).*sin(s_body).*para.hu./para.wh;
end
s_body = linspace(0,-pi/2,7);
for i = 1:n_nose
    vehicle.surface.nose{1}(i,7:13) = x_nose(i).*ones(size(s_body));
    vehicle.surface.nose{2}(i,7:13) = (interp1(Haa(:,1).*para.ln,Haa(:,2),x_nose(i)).*para.wh).*cos(s_body);
    vehicle.surface.nose{3}(i,7:13) = (interp1(Haa(:,1).*para.ln,Haa(:,2),x_nose(i)).*para.wh).*sin(s_body).*para.hl./para.wh;
end
%翼作成
%翼型
n_wing = 15;
n_kink = 5;
n_foil = 12;
thn = 0.025;
%翼型の指定はここ
for i = 1:n_wing
    chord{i}(1:n_foil,1) = (linspace(1,0,n_foil)').^1.4;
    chord{i}(1:n_foil,2)= -thn./0.2.*(0.29690.*sqrt(chord{i}(:,1))-0.12600.*chord{i}(:,1)-0.35160.*chord{i}(:,1).^2+0.28430.*chord{i}(:,1).^3-0.1036.*chord{i}(:,1).^4);
    chord{i}(n_foil:2*n_foil-1,1) = (linspace(0,1,n_foil)').^1.4;
    chord{i}(n_foil:2*n_foil-1,2)= thn./0.2.*(0.29690.*sqrt(chord{i}(n_foil:2*n_foil-1,1))-0.12600.*chord{i}(n_foil:2*n_foil-1,1)-0.35160.*chord{i}(n_foil:2*n_foil-1,1).^2+0.28430.*chord{i}(n_foil:2*n_foil-1,1).^3-0.1036.*chord{i}(n_foil:2*n_foil-1,1).^4);
end
y_span(1:n_kink) = linspace(para.wh,(para.semispan-para.wh).*para.Rkink+para.wh,n_kink);
y_span(n_kink:n_wing) = linspace((para.semispan-para.wh).*para.Rkink+para.wh,para.semispan,n_wing-n_kink+1);
Csec(1:n_kink) = (para.Cr-para.Ck)./(-(para.semispan-para.wh).*para.Rkink).*(y_span(1:n_kink)-para.wh)+para.Cr;
Csec(n_kink:n_wing) = (para.Ct-para.Ck)./(para.semispan-((para.semispan-para.wh).*para.Rkink+para.wh)).*(y_span(n_kink:end)-para.semispan)+para.Ct;
for i = 1:n_wing
    vehicle.surface.wing{1}(:,i) = (chord{i}(:,1)-1).*Csec(i)+(para.ln+para.Cr);
    vehicle.surface.wing{2}(:,i) = ones(size(chord{i}(:,1))).*y_span(i);
    vehicle.surface.wing{3}(:,i) = (chord{i}(:,2)).*Csec(i);
end
%エレボンをつけることを考慮してwingを分けておく
for i = 1:3
    vehicle.surface.wingkink{i} = vehicle.surface.wing{i}(:,1:n_kink);
    vehicle.surface.wingout{i} = vehicle.surface.wing{i}(:,n_kink:end);
end

%横のラインに対して胴体の楕円との交点を求め，ポイント数を維持したまま補間をする
R_body = para.hu;
b_coef = para.wh/para.hu;
for i = n_foil:size(vehicle.surface.wingkink{1},1)
        y_WB(i) = bisection(@(y)WBinteraction(y,R_body,b_coef,[vehicle.surface.wingkink{1}(i,:);vehicle.surface.wingkink{2}(i,:);vehicle.surface.wingkink{3}(i,:)]),0,para.wh);
        z_WB(i) = interp1(vehicle.surface.wingkink{2}(i,:),vehicle.surface.wingkink{3}(i,:),y_WB(i),'linear');
        x_WB(i) = interp1(vehicle.surface.wingkink{2}(i,:),vehicle.surface.wingkink{1}(i,:),y_WB(i),'linear');
        vehicle.surface.wingkink{2}(i,:) = linspace(y_WB(i),vehicle.surface.wing{2}(i,n_kink),n_kink);
        vehicle.surface.wingkink{1}(i,:) = interp1(vehicle.surface.wing{2}(i,:),vehicle.surface.wing{1}(i,:),vehicle.surface.wingkink{2}(i,:),'linear','extrap');
        vehicle.surface.wingkink{3}(i,:) = interp1(vehicle.surface.wing{2}(i,:),vehicle.surface.wing{3}(i,:),vehicle.surface.wingkink{2}(i,:),'linear','extrap');
end
%胴体設計
for i=1:n_foil
    s_body = linspace(pi/2,acos(y_WB(i+n_foil-1)/para.wh),size(s_body,2));
    vehicle.surface.bodyup{1}(i,1:size(s_body,2)) = ones(size(s_body)).*vehicle.surface.wingkink{1}(i+n_foil-1,1);
    vehicle.surface.bodyup{2}(i,1:size(s_body,2)) = para.wh.*cos(s_body);
    vehicle.surface.bodyup{3}(i,1:size(s_body,2)) = para.hu.*sin(s_body);
end

R_body = para.hl;
b_coef = para.wh/para.hl;
for i = 1:n_foil
        y_WB(i) = bisection(@(y)WBinteraction(y,R_body,b_coef,[vehicle.surface.wingkink{1}(i,:);vehicle.surface.wingkink{2}(i,:);vehicle.surface.wingkink{3}(i,:)]),0,para.wh);
        z_WB(i) = interp1(vehicle.surface.wingkink{2}(i,:),vehicle.surface.wingkink{3}(i,:),y_WB(i),'linear');
        x_WB(i) = interp1(vehicle.surface.wingkink{2}(i,:),vehicle.surface.wingkink{1}(i,:),y_WB(i),'linear');
        vehicle.surface.wingkink{2}(i,:) = linspace(y_WB(i),vehicle.surface.wing{2}(i,n_kink),n_kink);
        vehicle.surface.wingkink{1}(i,:) = interp1(vehicle.surface.wing{2}(i,:),vehicle.surface.wing{1}(i,:),vehicle.surface.wingkink{2}(i,:),'linear','extrap');
        vehicle.surface.wingkink{3}(i,:) = interp1(vehicle.surface.wing{2}(i,:),vehicle.surface.wing{3}(i,:),vehicle.surface.wingkink{2}(i,:),'linear','extrap');
end
for i=1:n_foil
    s_body = linspace(pi/2,acos(y_WB(n_foil-i+1)/para.wh),size(s_body,2));
    vehicle.surface.body3{1}(i,1:size(s_body,2)) = ones(size(s_body)).*vehicle.surface.wingkink{1}(n_foil-i+1,1);
    vehicle.surface.body3{2}(i,1:size(s_body,2)) = para.wh.*cos(s_body);
    vehicle.surface.body3{3}(i,1:size(s_body,2)) = -para.hl.*sin(s_body);
end
%垂直尾翼設計
n_vert = 5;
n_vfoil = 4;
vthn = 0.07;
%翼型の指定はここ
for i = 1:n_vert
    vchord{i}(1:n_vfoil,1) = linspace(1,0,n_vfoil)';
    vchord{i}(1:n_vfoil,2)= -vthn./0.2.*(0.29690.*sqrt(vchord{i}(:,1))-0.12600.*vchord{i}(:,1)-0.35160.*vchord{i}(:,1).^2+0.28430.*vchord{i}(:,1).^3-0.1036.*vchord{i}(:,1).^4);
    vchord{i}(n_vfoil:2*n_vfoil-1,1) = linspace(0,1,n_vfoil)';
    vchord{i}(n_vfoil:2*n_vfoil-1,2)= vthn./0.2.*(0.29690.*sqrt(vchord{i}(n_vfoil:2*n_vfoil-1,1))-0.12600.*vchord{i}(n_vfoil:2*n_vfoil-1,1)-0.35160.*vchord{i}(n_vfoil:2*n_vfoil-1,1).^2+0.28430.*vchord{i}(n_vfoil:2*n_vfoil-1,1).^3-0.1036.*vchord{i}(n_vfoil:2*n_vfoil-1,1).^4);
end
v_span = linspace(0,1,n_vert);
Cvert = -(para.Crv-para.Ctv).*(v_span)+para.Crv;
%楕円極座標表示
r = sqrt(1/((cos(para.phivert^2)/para.wh^2)+(sin(para.phivert^2)/para.hu^2)));
for i = 1:n_vert
    vehicle.surface.vert{1}(:,i) = (vchord{i}(:,1)-1).*Cvert(i)+(para.ln+para.Cr);
    for j = 1:2*n_vfoil-1
        vehicle.surface.vert{2}(j,i) = v_span(i).*para.hvert*cos(para.phivert)-(vchord{i}(j,2)).*Cvert(i)*sin(para.phivert)+r*cos(para.phivert);
        vehicle.surface.vert{3}(j,i) = v_span(i).*para.hvert*sin(para.phivert)+(vchord{i}(j,2)).*Cvert(i)*cos(para.phivert)+r*sin(para.phivert);
    end
end

R_body = para.hu;
b_coef = para.wh/para.hu;
for i = 1:2*n_vfoil-1
        y_VB(i) = bisection(@(y)WBinteraction(y,R_body,b_coef,[vehicle.surface.vert{1}(i,:);vehicle.surface.vert{2}(i,:);vehicle.surface.vert{3}(i,:)]),0,para.wh);
        z_VB(i) = interp1(vehicle.surface.vert{2}(i,:),vehicle.surface.vert{3}(i,:),y_VB(i),'linear','extrap');
        x_VB(i) = interp1(vehicle.surface.vert{2}(i,:),vehicle.surface.vert{1}(i,:),y_VB(i),'linear','extrap');
end
vehicle.surface.vert{2}(:,1) = y_VB';
vehicle.surface.vert{1}(:,1) = x_VB';
vehicle.surface.vert{3}(:,1) = z_VB';
%本当は垂直尾翼パネルを補完したほうがよい。形を戻す→再生成→戻す

%胴体上部を変更
for i = 1:3
    vehicle.surface.body1{i} = vehicle.surface.bodyup{i}(:,1:2);
    vehicle.surface.body2{i} = vehicle.surface.bodyup{i}(:,2:end);
    vehicle.surface.body1{i}(end-n_vfoil+1:end,end) = vehicle.surface.vert{i}(n_vfoil:end,1);
    vehicle.surface.body2{i}(end-n_vfoil+1:end,1) = flipud(vehicle.surface.vert{i}(1:n_vfoil,1));
end

xwake = 1000;
n_wake = 2;
for i = 1:size(vehicle.surface.wingkink{1}(end,1:end),2)
    vehicle.wake.wingkink{1}(:,i) = linspace(vehicle.surface.wingkink{1}(end,i),vehicle.surface.wingkink{1}(end,i)+xwake,n_wake)';
    vehicle.wake.wingkink{2}(:,i) = vehicle.surface.wingkink{2}(end,i).*ones(n_wake,1);
    vehicle.wake.wingkink{3}(:,i) = vehicle.surface.wingkink{3}(end,i).*ones(n_wake,1);
end
for i = 1:size(vehicle.surface.wingout{1}(end,1:end),2)
    vehicle.wake.wingout{1}(:,i) = linspace(vehicle.surface.wingout{1}(end,i),vehicle.surface.wingout{1}(end,i)+xwake,n_wake)';
    vehicle.wake.wingout{2}(:,i)  = vehicle.surface.wingout{2}(end,i).*ones(n_wake,1);
    vehicle.wake.wingout{3}(:,i)  = vehicle.surface.wingout{3}(end,i).*ones(n_wake,1);
end
for i = 1:size(vehicle.surface.vert{1}(end,1:end),2)
    vehicle.wake.vert{1}(:,i) = linspace(vehicle.surface.vert{1}(end,i),vehicle.surface.vert{1}(end,i)+xwake,n_wake)';
    vehicle.wake.vert{2}(:,i)  = vehicle.surface.vert{2}(end,i).*ones(n_wake,1);
    vehicle.wake.vert{3}(:,i)  = vehicle.surface.vert{3}(end,i).*ones(n_wake,1);
end
%接続の確認
for i = 1:3
    vehicle.surface.body2{i}(:,end) = vehicle.surface.wingkink{i}(n_foil:end,1);
    vehicle.surface.body3{i}(:,end) = vehicle.surface.wingkink{i}(n_foil:-1:1,1);
    vehicle.surface.nose{i}(end,:) = [vehicle.surface.body1{i}(1,1:end-1),vehicle.surface.body2{i}(1,1:end-1),fliplr(vehicle.surface.body3{i}(1,1:end))];
    vehicle.surface.wingout{i}(end,:)=vehicle.surface.wingout{i}(1,:);
    vehicle.surface.wingkink{i}(end,:)=vehicle.surface.wingkink{i}(1,:);
    vehicle.surface.vert{i}(1,:) = vehicle.surface.vert{i}(end,:);
    vehicle.wake.wingout{i}(1,:) = vehicle.surface.wingout{i}(end,:);
    vehicle.wake.wingkink{i}(1,:) = vehicle.surface.wingkink{i}(end,:);
    vehicle.wake.wing{i} = [vehicle.wake.wingout{i}(:,1:end-1),vehicle.wake.wingkink{i}];
end

%翼型の端っことwake
for i = 1:3
    vehicle.surface.wingtip{i} = [vehicle.surface.wingout{i}(1:n_foil,end)';vehicle.surface.wingout{i}(end:-1:n_foil,end)'];
    vehicle.surface.verttip{i} = [vehicle.surface.vert{i}(1:n_vfoil,end)';vehicle.surface.vert{i}(end:-1:n_vfoil,end)'];
end

%ベース面
n_base = 4;
Rbase = linspace(1,0,n_base);
for i = 1:n_base
    vehicle.surface.base{1}(i,:) = [vehicle.surface.body1{1}(end,1:end-1),vehicle.surface.body2{1}(end,1:end-1),fliplr(vehicle.surface.body3{1}(end,1:end))];
    vehicle.surface.base{2}(i,:) = [vehicle.surface.body1{2}(end,1:end-1),vehicle.surface.body2{2}(end,1:end-1),fliplr(vehicle.surface.body3{2}(end,1:end))].*Rbase(i);
    vehicle.surface.base{3}(i,:) = [vehicle.surface.body1{3}(end,1:end-1),vehicle.surface.body2{3}(end,1:end-1),fliplr(vehicle.surface.body3{3}(end,1:end))].*Rbase(i);
end
%ベースwake
%vehicle.wake.base{1}  = [vehicle.surface.base{1}(1,1:end);vehicle.surface.base{1}(1,1:end)+xwake];
%vehicle.wake.base{2}  = [vehicle.surface.base{2}(1,1:end);vehicle.surface.base{2}(1,1:end)];
%vehicle.wake.base{3}  = [vehicle.surface.base{3}(1,1:end);vehicle.surface.base{3}(1,1:end)];

%REFS
SREF = (para.Cr+para.Ct)*para.semispan;
BREF = para.semispan*2;
CREF = para.ln+para.Cr;
vehicle.REFS = [0 0 0 0;SREF,BREF,CREF,1];

%{
figure(7);clf;hold on
mesh(vehicle.surface.nose{1},vehicle.surface.nose{2},vehicle.surface.nose{3})
mesh(vehicle.surface.wingkink{1},vehicle.surface.wingkink{2},vehicle.surface.wingkink{3})
mesh(vehicle.surface.wingout{1},vehicle.surface.wingout{2},vehicle.surface.wingout{3})
mesh(vehicle.surface.body1{1},vehicle.surface.body1{2},vehicle.surface.body1{3})
mesh(vehicle.surface.body2{1},vehicle.surface.body2{2},vehicle.surface.body2{3})
mesh(vehicle.surface.body3{1},vehicle.surface.body3{2},vehicle.surface.body3{3})
mesh(vehicle.surface.vert{1},vehicle.surface.vert{2},vehicle.surface.vert{3})
mesh(vehicle.surface.base{1},vehicle.surface.base{2},vehicle.surface.base{3})
mesh(vehicle.surface.wingtip{1},vehicle.surface.wingtip{2},vehicle.surface.wingtip{3})
mesh(vehicle.surface.verttip{1},vehicle.surface.verttip{2},vehicle.surface.verttip{3})
mesh(vehicle.wake.wing{1},vehicle.wake.wing{2},vehicle.wake.wing{3})
axis equal
%}
end
function [chord] = Haack(n)
    C=0;
    x = linspace(0,1,n)';
    h=acos(1-2.*x);
    r= sqrt((h-(1./2).*sin(2.*h)+(C.*((sin(h)).^3).*h))./pi);
    chord = [x r];
end

function res = WBinteraction(y,R_body,b_coef,line)
%x = interp1(line(2,:),line(1,:),y);
z = interp1(line(2,:),line(3,:),y,'linear','extrap');
res = y^2/(R_body*b_coef)^2+z^2/R_body^2-1;
end


function p = bisection(fun,a,b)

% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code. But may come handy
i = 1;
if fun(a)*fun(b)>0 
    p = (a + b)/2;
else
    p = (a + b)/2;
    err = abs(fun(p));
    while(err > 1e-7 && i<=10000) 
       if fun(a)*fun(p)<0 
           b = p;
       else
           a = p;          
       end
        p = (a + b)/2; 
       err = abs(fun(p));
       i = i+1;
    end
end
end
