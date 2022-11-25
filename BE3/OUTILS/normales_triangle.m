function nT = normales_triangle(INV_DFT,nref)
 
 nT = INV_DFT*nref;%Calcul des normales unitaires sortantes à T
 for l =1:3 
  nT(:,l) = nT(:,l)/norm(nT(:,l));%normalisation des normales
 end 
 end