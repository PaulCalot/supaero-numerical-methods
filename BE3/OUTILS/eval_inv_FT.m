function [r,s] = eval_inv_FT(x,y,xT,yT)
%%(x,y) points du triangle T
%% (xT(i),yT(i)) = coordonnées du i ème sommet du triangle T

A = [xT(2)-xT(1) xT(3)-xT(1);
      yT(2)-yT(1) yT(3)-yT(1)];
      
V = A\[x-xT(1);y-yT(1)];      
r = V(1);
s = V(2);
end 