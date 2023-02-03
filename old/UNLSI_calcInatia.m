function [ina] = UNLSI_calcInatia(frame,cgpos,surfThn,surfRho)
ina.Ixx = 0;ina.Iyy = 0;ina.Izz = 0;ina.mass = 0;


for i = 1:numel(frame.surfID)
   if frame.isWake(i,1) == 0
       ina.mass = ina.mass + frame.area(i)*surfThn(i)*surfRho(i);
       ina.Ixx = ina.Ixx + frame.area(i)*surfThn(i)*surfRho(i)*((frame.center(i,2)-cgpos(2))^2+(frame.center(i,3)-cgpos(3))^2);
       ina.Iyy = ina.Iyy + frame.area(i)*surfThn(i)*surfRho(i)*((frame.center(i,1)-cgpos(1))^2+(frame.center(i,3)-cgpos(3))^2);
       ina.Izz = ina.Izz + frame.area(i)*surfThn(i)*surfRho(i)*((frame.center(i,1)-cgpos(1))^2+(frame.center(i,2)-cgpos(2))^2);
   end
end


end