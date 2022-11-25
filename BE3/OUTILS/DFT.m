function DFT = DFT(xT,yT)
%%Calcul de la matrice jacobienne associée au triangle T
%% (xT(i),yT(i)) = coordonnées du i ème sommet du triangle T

DFT= [xT(2)-xT(1) xT(3)-xT(1);
      yT(2)-yT(1) yT(3)-yT(1)];
end 