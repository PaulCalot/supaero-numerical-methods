function [x,y] = eval_Fa(r,xa,ya)
%% r points du segment de référence
%% (xa(i),ya(i)) = coordonnées de la i ème extrémité de l arete a
x = xa(1)*(1-r)+xa(2)*r;
y = ya(1)*(1-r)+ya(2)*r;

end 