function frameCpView(verts,frame,val,HR,figNo,subNo,viewVal,scale)
%surfIDÇÃêîÇéÊìæ
figure(figNo);
set(gcf,'Position',[100 100 840 630])
subplot(subNo(1),subNo(2),subNo(3));

hold on;set(gca, 'LineWidth',1.0, 'FontSize',10, 'Fontname','Arial')
plotVal(verts,frame,val,[1,2,3],HR);
view(viewVal);drawnow;axis equal;hold off;colorbar;caxis(scale);colormap jet;axis off

end

function plotVal(verts,frame,val,axisNo,HR)

for i = 1:size(frame.tri,1)
   if frame.isBody(i,1)== 1 || frame.isBase(i,1)== 1
       %trisurf([1,2,3],verts(frame.tri(i,:),axisNo(1)).*HR,verts(frame.tri(i,:),axisNo(2)).*HR,verts(frame.tri(i,:),axisNo(3)).*HR,val(i),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
       trisurf([1,2,3],verts(frame.tri(i,:),axisNo(1)).*HR,verts(frame.tri(i,:),axisNo(2)).*HR,verts(frame.tri(i,:),axisNo(3)).*HR,val(i))
   end
   %
   if frame.isStr(i,1)== 1
       symmmat =[1 1 1];
       symmmat(axisNo==2) =  -symmmat(axisNo==2);
       trisurf([1,2,3],symmmat(1).*verts(frame.tri(i,:),axisNo(1)).*HR,symmmat(2).*verts(frame.tri(i,:),axisNo(2)).*HR,symmmat(3).*verts(frame.tri(i,:),axisNo(3)).*HR,val(i)) 
       %trisurf([1,2,3],symmmat(1).*verts(frame.tri(i,:),axisNo(1)).*HR,symmmat(2).*verts(frame.tri(i,:),axisNo(2)).*HR,symmmat(3).*verts(frame.tri(i,:),axisNo(3)).*HR,val(i),'EdgeColor','none','LineStyle','none','FaceLighting','phong') 
   end
   %}
end
end

    