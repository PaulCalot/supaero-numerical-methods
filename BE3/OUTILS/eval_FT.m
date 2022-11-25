function [x,y] = eval_FT(r,s,xT,yT)
%%(r,s) points du triangle de référence
%% (xT(i),yT(i)) = coordonnées du i ème sommet du triangle T

x = xT(1)*(1-r-s)+xT(2)*r+xT(3)*s;
y = yT(1)*(1-r-s)+yT(2)*r+yT(3)*s;

end 