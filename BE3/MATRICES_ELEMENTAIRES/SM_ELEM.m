function ZELEM_SM = SM_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xa,ya,la)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% (xa,ya)(1:2 = coordonnées des extrémités de l'arete 
% la = longueur de l'arete 


N = (ordre +1)*(ordre+2)/2; 
ZELEM_SM = zeros(N,N);
 
for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:),ya(:)); 
  [r,s] = eval_inv_FT(x,y,xT,yT);
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  
  for l=1:N 
   fct_basel = fct_base(l,A,r,s,ordre);
   ZELEM_SM(k,l) = ZELEM_SM(k,l) +wg1D(n)*la*fct_basel*fct_basek;
   
  end 
 end 
end 


end 