function Hardness = airfoil_hardness(E,thn,C,chord,thn_pos)
chord = chord.*C;
n_foil = round(size(chord,1)/2);
Ycenter = 0.5*(interp1(chord(1:n_foil,1),chord(1:n_foil,2),thn_pos,'linear','extrap')+interp1(chord(n_foil:2*n_foil-1,1),chord(n_foil:2*n_foil-1,2),thn_pos,'linear','extrap'));
Hardness(1:2) = 0;
for i = 1:size(chord,1)-1
   Hardness(1) = Hardness(1)+sqrt((chord(i+1,1)-chord(i,1))^2+(chord(i+1,2)-chord(i,2))^2)*thn*(((chord(i+1,1)+chord(i,1))/2-thn_pos)^2);
   Hardness(2) = Hardness(2)+sqrt((chord(i+1,1)-chord(i,1))^2+(chord(i+1,2)-chord(i,2))^2)*thn*(((chord(i+1,2)+chord(i,2))/2-Ycenter)^2);
end
Hardness = Hardness.*E;
end